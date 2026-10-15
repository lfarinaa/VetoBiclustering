"""Measures the per-tile backsplash map from the surviving Geant4 photon sample and writes it to data/psdTileBacksplash.csv.

    python tools/measureBacksplash.py [--selection thesis]

For each PSD tile: the fraction of L0-triggered fiducial photons (stochastic/photonEvents.py) in which the tile fires above 0.333 MeV, with the
photons weighted as generated (log-uniform in 10 MeV - 100 GeV, as the thesis' energy-binned figures) and reweighted to dN/dE ∝ E^-2 (as thesis
Fig. 4.4; with this weighting the full-PSD veto keeps 0.90 of the photons, Table 4.2 quotes 0.89). Each map is shrunk towards its mean over tiles
with a pseudo-count of 20 photons: the sample is 340 fiducial photons for 1328 tiles, so single-tile values are noisy and the map is an estimate
of where backsplash lands, not a precision measurement. Columns: face, layer, tileId, logUniform, spectrumE2 (probabilities per photon).
"""
import argparse
import pathlib
import sys

root = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root / "analytic"))
sys.path.insert(0, str(root / "stochastic"))
import patternDump as pd
import photonEvents as pe


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--selection", default="thesis", choices=list(pd.triggerSelections))
    parser.add_argument("--output", type=pathlib.Path, default=root / "data" / "psdTileBacksplash.csv")
    args = parser.parse_args()
    patterns = pd.selectPatterns(pd.parsePatternDump(), **pd.triggerSelections[args.selection])
    tileKeys = pd.allTileKeys()
    triggered = pe.triggerPhotons(pe.loadPhotons(), patterns, tileKeys)
    logUniform, events = pe.tileFiringProbability(triggered, len(tileKeys))
    spectrum, _ = pe.tileFiringProbability(triggered, len(tileKeys), spectralIndex=2.0)
    lines = ["face,layer,tileId,logUniform,spectrumE2"] + [f"{f},{l},{t},{a:.4e},{b:.4e}" for (f, l, t), a, b in zip(tileKeys, logUniform, spectrum)]
    args.output.write_text("\n".join(lines) + "\n")
    print(f"{events} fiducial L0-triggered photons ({args.selection}); mean per tile: log-uniform {logUniform.mean():.2e}, E^-2 {spectrum.mean():.2e} -> {args.output}")


if __name__ == "__main__":
    main()
