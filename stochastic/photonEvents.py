"""Event-level L1 self-veto of photons, from the Geant4 photon sample that survives on this disk.

The only simulated events left locally are two AeffCalc analysis files of photons (oneshot_fidgammas_compact: 2 x 100 000 photons from the upper
hemisphere, log-uniform in 10 MeV - 100 GeV, compact-v2 layout; the rest of the samples were on a cluster that is no longer accessible). They hold, hit
by hit, the FIT mat hits (fitHitsData: side, plane, mat, energy deposit) and the PSD hits (psdHitsData: side, layer, tile, energy deposit), each tagged
with the photon's true energy and direction, which identify the event (the files have no event number); SaveEAnglesFiducial lists the photons of the
fiducial sample (thesis section 4.4).

For every photon this module redoes the trigger with the RoIs of the regenerated dump (data/patterns_allPlanes.txt, thesis-era code):

- L0: the 3IR patterns of the chosen selection whose three mats are hit (FIT hits above 10^-4 GeV, thesis section 4.3.2); as in TriggerArrays.cpp a
  side triggers if some double layer has a pattern in both views, and the patterns of its first triggered double layer are kept;
- L1: as in MultiRoiVeto.cpp, a side is clean if no PSD tile above 0.333 MeV (thesis section 2.1.3) fires in (U_i V_i) ∩ (U_j V_j), the masks of its
  even- and odd-view patterns, and the event passes if any triggered side is clean.

For photons every veto is a self-veto (backsplash, or the pair or a Compton electron reaching the PSD), so the fraction of L0-triggered photons passing
L1 is the efficiency the channel map keeps, with the real energy dependence and the real correlations between tiles. Simplifications: the calorimeter
requirement is taken per pattern (requireCalo), not per pair; the energy deposit of the calorimeter (L2) is not applied.
Naming: camelCase everywhere. Units: GeV.
"""
import glob
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analytic"))
import patternDump as pd

defaultFiles = sorted(glob.glob(str(pathlib.Path.home() / "AeffCalc/oneshot_fidgammas_compact/output_analysis/compact/*.analyzed.root")))
sideNames = ("Xpos", "Xneg", "Ypos", "Yneg", "TOP")   # Herd::RefFrame::Side, HerdSoftware include/analysis/common/RefFrame.h.
fitThreshold = 1e-4                                  # GeV. FIT mat activation threshold (thesis section 4.3.2).
psdThreshold = 3.33e-4                               # GeV. PSD tile activation threshold (thesis section 2.1.3).
generatedRange = (0.01, 100.0)                       # GeV. Log-uniform generation (thesis section 2.2.1).


def loadPhotons(files=defaultFiles):
    """[dict(energy, cosPolar, azimuth, fiducial, fitHits {(side, plane, mat)}, psdHits {(side, layer, tile)})] for every photon with a FIT hit."""
    import uproot
    events = []
    for fileIndex, path in enumerate(files):
        with uproot.open(path) as rootFile:
            fit = rootFile["fitHitsData"].arrays(library="np")
            psd = rootFile["psdHitsData"].arrays(library="np")
            fiducial = rootFile["SaveEAnglesFiducial"].arrays(library="np")
        fiducialKeys = set(zip(fiducial["eTrue"], fiducial["cosPolar"], fiducial["azimuth"]))
        byEvent = {}
        mask = fit["eDep"] > fitThreshold
        for e, c, a, s, l, m in zip(fit["eTrue"][mask], fit["cosPolar"][mask], fit["azimuth"][mask], fit["iSide"][mask], fit["iLayer"][mask], fit["iMat"][mask]):
            byEvent.setdefault((e, c, a), (set(), set()))[0].add((sideNames[int(s)], int(l), int(m)))
        mask = psd["eDep"] > psdThreshold
        for e, c, a, s, l, t in zip(psd["eTrue"][mask], psd["cosPolar"][mask], psd["azimuth"][mask], psd["iSide"][mask], psd["iLayer"][mask], psd["iTile"][mask]):
            if (e, c, a) in byEvent:
                byEvent[(e, c, a)][1].add((sideNames[int(s)], int(l), int(t)))
        for (e, c, a), (fitHits, psdHits) in byEvent.items():
            events.append(dict(energy=float(e), cosPolar=float(c), azimuth=float(a), fiducial=(e, c, a) in fiducialKeys, file=fileIndex,
                               fitHits=fitHits, psdHits=psdHits))
    return events


def triggerPhotons(events, patterns, tileKeys):
    """L0 for every event against a pattern set (rows of an incidence matrix). Returns the triggered events, each with
    sides = [(evenRows, oddRows)] (first triggered double layer of each triggered side) and fired = column indices of the PSD tiles above threshold."""
    row = {p.name: i for i, p in enumerate(patterns)}
    column = {key: j for j, key in enumerate(tileKeys)}
    triggered = []
    for event in events:
        hits = event["fitHits"]
        sides = []
        for side in sideNames:
            found = {}
            for (s, plane, mat) in hits:
                if s != side or plane > 9:
                    continue
                for m1 in (mat - 1, mat, mat + 1):
                    if (side, plane + 2, m1) not in hits:
                        continue
                    for m2 in (m1 - 1, m1, m1 + 1):
                        if (side, plane + 4, m2) in hits:
                            name = f"{side}{plane}{mat:02d}5{5 + m1 - mat}{5 + m2 - m1}"
                            if name in row:
                                found.setdefault(plane, []).append(row[name])
            doubleLayers = [plane for plane in sorted(found) if plane % 2 == 0 and plane + 1 in found]
            if doubleLayers:
                first = doubleLayers[0]
                sides.append((np.array(sorted(set(found[first]))), np.array(sorted(set(found[first + 1])))))
        if sides:
            fired = np.array(sorted(column[key] for key in event["psdHits"] if key in column), dtype=int)
            triggered.append(dict(event, sides=sides, fired=fired))
    return triggered


def passesL1(event, mask):
    """True if some triggered side has no fired tile inside (U_i V_i) ∩ (U_j V_j)."""
    if event["fired"].size == 0:
        return True
    for evenRows, oddRows in event["sides"]:
        region = mask[evenRows].any(axis=0) & mask[oddRows].any(axis=0)
        if not region[event["fired"]].any():
            return True
    return False


def _passProbability(event, mask, logSurvival):
    """P(L1 pass) of one event: a side passes if no fired tile is in its region and no accidental hit (independent per tile, log(1 - q_t) given)
    lands there; the event passes if any side passes (sides treated as independent for the accidentals)."""
    failAll = 1.0
    for evenRows, oddRows in event["sides"]:
        region = mask[evenRows].any(axis=0) & mask[oddRows].any(axis=0)
        clean = 0.0 if (event["fired"].size and region[event["fired"]].any()) else 1.0
        accidentalSurvival = 1.0 if logSurvival is None else float(np.exp(logSurvival[region].sum()))
        failAll *= 1.0 - clean * accidentalSurvival
    return 1.0 - failAll


def spectralWeights(sample, index=None):
    """Event weights: uniform (None: the log-uniform generated spectrum, as the thesis' tables) or reweighted to dN/dE ∝ E^-index."""
    energy = np.array([e["energy"] for e in sample])
    return np.ones_like(energy) if index is None else energy ** (1.0 - index)


def l1Efficiency(triggered, mask, energyBins=None, fiducialOnly=True, accidentalProbability=None, spectralIndex=None):
    """L1 efficiency of L0-triggered photons (by default the fiducial ones) with the given mask (patterns x tiles): the weighted fraction passing,
    with optional accidental hits (per-tile probability within the coincidence window, rateUtils) and spectral reweighting.
    Returns (efficiency, efficiency per energy bin or None, events per bin or None)."""
    sample = [e for e in triggered if e["fiducial"] or not fiducialOnly]
    logSurvival = None if accidentalProbability is None else np.log1p(-np.asarray(accidentalProbability, dtype=float))
    passed = np.array([_passProbability(e, mask, logSurvival) for e in sample])
    weights = spectralWeights(sample, spectralIndex)
    total = float((weights * passed).sum() / weights.sum())
    if energyBins is None:
        return total, None, None
    energy = np.array([e["energy"] for e in sample])
    index = np.digitize(energy, energyBins) - 1
    size = len(energyBins) - 1
    counts = np.bincount(index, minlength=size)[:size]
    good = np.bincount(index, weights=passed, minlength=size)[:size]
    return total, np.where(counts > 0, good / np.maximum(counts, 1), np.nan), counts


def tileFiringProbability(triggered, numberOfTiles, fiducialOnly=True, energyRange=None, priorStrength=20.0, spectralIndex=None):
    """The measured backsplash map: per tile, the (weighted) fraction of L0-triggered fiducial photons in which the tile fires, shrunk towards the
    mean over tiles with a pseudo-count of priorStrength events (the sample is a few hundred photons for 1328 tiles). Returns (probability, events)."""
    sample = [e for e in triggered if (e["fiducial"] or not fiducialOnly) and (energyRange is None or energyRange[0] <= e["energy"] < energyRange[1])]
    weights = spectralWeights(sample, spectralIndex)
    counts = np.zeros(numberOfTiles)
    for event, weight in zip(sample, weights):
        counts[event["fired"]] += weight
    total = weights.sum()
    mean = counts.sum() / max(total, 1e-300) / numberOfTiles
    return (counts + priorStrength * mean) / (total + priorStrength), len(sample)
