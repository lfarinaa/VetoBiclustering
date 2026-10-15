"""Finds (or evaluates) a channel map for the ULEG L1 veto and lists the combination: which 3IR patterns share a veto decision (pattern groups), which
PSD tiles share a channel (OR of tiles), and which channels veto which group (routing). Notation: examples/channelBudgetShowcase.ipynb, section 1.

    python tools/findChannelMap.py --groups 16 --channels 64                       # optimise under a budget, all five FIT faces
    python tools/findChannelMap.py --groups 5 --channels 17 --groups-per-fit-face   # hardware segmented per FIT face
    python tools/findChannelMap.py --strategy face-rows --rows 2                    # evaluate the thesis' hybrid (section 4.2.2)
    python tools/findChannelMap.py --strategy lossless --fit-faces TOP --json top.json
    python tools/findChannelMap.py --groups 16 --channels 64 --window 1000 --selection calo8   # 1 us window, the 2023 trigger configuration

Strategies: optimise (default: greedy merges in three orders plus a one-group-per-FIT-face seed, each refined by local moves; the best is kept),
lossless (N_G*, N_C*), full-psd (one group, one channel), face-union (one group per FIT face, the union of its RoIs), face-rows (one group per FIT face,
its PSD face plus --rows rows of each adjacent face). Patterns: the regenerated dump (data/patterns_allPlanes.txt), restricted to a trigger
configuration (--selection; default 'thesis': CALO-intersecting, third hit adjacent to the first and no zig-zags, no upward-pointing patterns (upgoing: the particle moving towards the top), all planes). The cost of adding tile t to a
mask is -ln(1 - pi_t), with pi_t the probability that the tile vetoes a good photon: backsplash (measured per tile on Geant4 photons) or an accidental
hit from the measured tile rates within the coincidence window (--window, a placeholder of 100 ns); --cost-model excess counts added pairs instead,
and --cost takes any (patterns x 1328) matrix. The optimiser is a heuristic: it returns a good local optimum, not a proof.
"""
import argparse
import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analytic"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "stochastic"))
import patternDump as pd
import losslessUtils as lu
import mergeUtils as mu
import geometryUtils as gu
import rateUtils as ru


def buildStrategy(name, incidence, rows):
    matrix = incidence.matrix
    fitFaces = [p.fitSide for p in incidence.patterns]
    faces = sorted(set(fitFaces), key=pd.psdSides.index)
    groupOf = np.array([faces.index(face) for face in fitFaces])
    if name == "lossless":
        return lu.findLosslessChannelMap(matrix)
    if name == "full-psd":
        return lu.ChannelMap(np.zeros(matrix.shape[0], int), np.zeros(matrix.shape[1], int), np.ones((1, 1), bool))
    if name == "face-union":
        return mu.channelMapFromMasks(matrix, groupOf, np.array([matrix[groupOf == k].any(axis=0) for k in range(len(faces))]))
    if name == "face-rows":
        return mu.channelMapFromMasks(matrix, groupOf, np.array([gu.faceExtendedMask(incidence.tileKeys, face, rows) for face in faces]))
    raise ValueError(name)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--groups", type=int, help="N_G, the number of pattern groups (veto decisions)")
    parser.add_argument("--channels", type=int, help="N_C, the number of PSD channels (ORs of tiles)")
    parser.add_argument("--strategy", default="optimise", choices=["optimise", "lossless", "full-psd", "face-union", "face-rows"])
    parser.add_argument("--rows", type=int, default=2, help="rows of the adjacent faces in face-rows (thesis: up to two)")
    parser.add_argument("--selection", default="thesis", choices=list(pd.triggerSelections),
                        help="3IR pattern set (patternDump.triggerSelections; default: CALO-intersecting, no two-step codes, no upward-pointing (upgoing) patterns, all planes)")
    parser.add_argument("--fit-faces", nargs="+", default=list(pd.psdSides), choices=pd.psdSides, help="restrict to the patterns of these FIT faces")
    parser.add_argument("--groups-per-fit-face", action="store_true", help="a group may not mix patterns of different FIT faces")
    parser.add_argument("--channels-per-psd-face", action="store_true", help="a channel may not mix tiles of different PSD faces")
    parser.add_argument("--channels-per-layer", action="store_true", help="a channel may not mix tiles of different PSD layers")
    parser.add_argument("--cost-model", default="rates", choices=["rates", "excess"],
                        help="rates (default): C_it = -ln(1 - pi_t), backsplash plus accidental hits in the window; excess: 1 per added (pattern, tile)")
    parser.add_argument("--window", type=float, default=ru.defaultWindow * 1e9, help="coincidence window of the L1 veto in ns (default 100, a placeholder)")
    parser.add_argument("--cost", type=pathlib.Path, help="(patterns x 1328) .npy cost matrix in the order of the selected patterns; overrides --cost-model")
    parser.add_argument("--no-refine", action="store_true", help="skip the local refinement")
    parser.add_argument("--backsplash", default="spectrumE2",
                        help="backsplash per tile and photon: 'spectrumE2' (default) or 'logUniform' weighting of the measured map (data/psdTileBacksplash.csv), or a number")
    parser.add_argument("--events", action="store_true", help="also measure the L1 efficiency on the Geant4 fiducial photons (needs uproot and the AeffCalc files)")
    parser.add_argument("--dump", type=pathlib.Path, default=pd.defaultDumpPath)
    parser.add_argument("--output", type=pathlib.Path, help="write the listing here as well as to the terminal")
    parser.add_argument("--json", type=pathlib.Path, help="write the map as JSON (groups, channels, routing, metrics)")
    args = parser.parse_args()

    selected = [p for p in pd.selectPatterns(pd.parsePatternDump(args.dump), **pd.triggerSelections[args.selection]) if p.fitSide in args.fit_faces]
    full = pd.buildIncidence(selected)
    incidence = full.withoutEmptyColumns()
    keepColumns = [full.tileKeys.index(key) for key in incidence.tileKeys]
    backsplash = args.backsplash if args.backsplash in ("spectrumE2", "logUniform") else float(args.backsplash)
    vetoProbability = ru.tileVetoProbability(incidence.tileKeys, window=args.window * 1e-9, backsplash=backsplash)
    if args.cost is not None:
        cost = np.load(args.cost)[:, keepColumns]
    else:
        cost = ru.costMatrix(len(selected), vetoProbability) if args.cost_model == "rates" else None
    matrix = incidence.matrix

    rowLabels = np.array([p.fitSide for p in incidence.patterns]) if args.groups_per_fit_face else None
    colParts = []
    if args.channels_per_psd_face:
        colParts.append([key[0] for key in incidence.tileKeys])
    if args.channels_per_layer:
        colParts.append([str(key[1]) for key in incidence.tileKeys])
    colLabels = np.array(["/".join(parts) for parts in zip(*colParts)]) if colParts else None

    start = time.time()
    candidates = []
    if args.strategy == "optimise":
        if args.groups is None or args.channels is None:
            parser.error("--strategy optimise needs --groups and --channels")
        seeds = [("one group per FIT face", buildStrategy("face-union", incidence, args.rows))] if len(args.fit_faces) > 1 else []
        channelMap, candidates = mu.findChannelMap(matrix, args.groups, args.channels, cost, rowLabels, colLabels, seeds=seeds, refine=not args.no_refine)
    else:
        channelMap = buildStrategy(args.strategy, incidence, args.rows)

    reference = matrix.sum() if cost is None else (cost * matrix).sum()
    added, leaked = mu.computeMaskCost(matrix, channelMap, cost)
    if cost is not None:
        leaked = mu.computeMaskCost(matrix, channelMap)[1]
    pairs = mu.computePairMetrics(matrix, channelMap.vetoMask(), pd.evenOddBlocks(incidence.patterns), vetoProbability)
    metrics = dict(numberOfGroups=int(channelMap.routing.shape[0]), numberOfChannels=int(channelMap.routing.shape[1]),
                   routingLinks=int(channelMap.routing.sum()), maskExcess=added / reference, maskLeak=leaked / matrix.sum(),
                   eventExcess=pairs["excess"], eventLeak=pairs["leak"], relativeL1Efficiency=pairs.get("selfVeto"),
                   backsplashProbability=args.backsplash, windowNs=args.window, selection=args.selection, costModel=args.cost_model, seconds=time.time() - start)

    header = [f"Strategy {args.strategy}, selection {args.selection}, FIT faces {' '.join(args.fit_faces)}: {matrix.shape[0]} patterns x {matrix.shape[1]} tiles; "
              f"cost {args.cost_model if args.cost is None else args.cost}",
              f"per pattern:  excess xi = {metrics['maskExcess']:.3f}, leak kappa = {metrics['maskLeak']:.4f}",
              f"per event (Eq. 4.2 pairs):  excess = {metrics['eventExcess']:.3f}, leak = {metrics['eventLeak']:.4f}, "
              f"relative L1 efficiency = {metrics['relativeL1Efficiency']:.3f} (toy: independent tiles, backsplash {args.backsplash}, accidentals in a {args.window:g} ns window)"]
    if args.events:
        import photonEvents as pe
        triggered = pe.triggerPhotons(pe.loadPhotons(), incidence.patterns, incidence.tileKeys)
        accidental = ru.tileVetoProbability(incidence.tileKeys, window=args.window * 1e-9, backsplash=0.0)
        mask, exact = channelMap.vetoMask(), matrix
        for label, index in (("log-uniform", None), ("E^-2", 2.0)):
            mapped = pe.l1Efficiency(triggered, mask, accidentalProbability=accidental, spectralIndex=index)[0]
            reference = pe.l1Efficiency(triggered, exact, accidentalProbability=accidental, spectralIndex=index)[0]
            metrics[f"eventL1Efficiency_{label}"], metrics[f"eventL1EfficiencyExactRoI_{label}"] = mapped, reference
            header.append(f"Geant4 fiducial photons ({sum(e['fiducial'] for e in triggered)}), {label}: L1 efficiency {mapped:.3f} (exact RoI {reference:.3f})")
    if candidates:
        header.append("candidates (added cost / RoI cost): " + ", ".join(f"{name} {value / reference:.3f}" for name, value in candidates))
    listing = "\n".join(header) + "\n\n" + mu.describeChannelMap(channelMap, incidence.patterns, incidence.tileKeys)
    print(listing)
    if args.output:
        args.output.write_text(listing + "\n")
    if args.json:
        args.json.write_text(json.dumps(dict(
            metrics=metrics,
            groups={f"G{k}": [incidence.patterns[i].name for i in np.flatnonzero(channelMap.rowGroup == k)] for k in range(channelMap.routing.shape[0])},
            channels={f"C{l}": [list(map(lambda x: x if isinstance(x, str) else int(x), incidence.tileKeys[t])) for t in np.flatnonzero(channelMap.colGroup == l)]
                      for l in range(channelMap.routing.shape[1])},
            routing={f"G{k}": [f"C{l}" for l in np.flatnonzero(channelMap.routing[k])] for k in range(channelMap.routing.shape[0])}), indent=1))


if __name__ == "__main__":
    main()
