"""Draws the figures of the README and docs/HERD_SHOWCASE.md that are not notebook outputs, from the saved results in examples/results/ (no search is rerun):

    python tools/makeDocFigures.py

- docs/figures/psdRates.png: the PSD activation rate per tile, recovered from thesis Fig. 4.9, on the unfolded PSD;
- docs/figures/excessVsChannels.png: tiles added to the event RoIs by the optimised maps against the channel budget (optimisedGrid.pkl, from the
  channelBudgetShowcase notebook);
- docs/figures/constrainedMap3d.png: the (5, 8) map with one group per FIT face and every channel inside one PSD face, on the PSD in 3D from two sides
  (constrained.pkl, same notebook);
- docs/figures/jaccardMaps3d.png: the refined (5, 8) maps of the Jaccard merge and of the exact-cost greedy, in 3D from two sides
  (jaccardComparison.pkl, from tools/compareJaccardMerge.py);
- docs/figures/chipLayout.png: the best (5, 8) layout with front-end chips (45 chips of 32 inputs, boards contiguous on the cube), after the
  polish, on the unfolded PSD (chipLayout.pkl, from examples/chipLayoutShowcase.ipynb).
"""
import pathlib
import pickle
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch

root = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root / "analytic")); sys.path.insert(0, str(root / "stochastic"))
import patternDump as pd
import losslessUtils as lu
import mergeUtils as mu
import geometryUtils as gu
import plotUtils as pu
import rateUtils as ru
import chipLayout as cl

resultsDirectory, figureDirectory = root / "examples/results", root / "docs/figures"
views = (-60, 120)                                                         # Azimuths in degrees: together they show all five PSD faces.


def faceSummary(faces):
    """The one or two most common faces in `faces`, with "+" when more are present."""
    names, counts = np.unique(faces, return_counts=True)
    order = np.argsort(-counts)
    return " ".join(names[order[:2]]) + (" +" if len(names) > 2 else "")


def drawMapsIn3d(path, namedMaps, tileKeys, fitFace):
    """One row per map: the PSD cube from both views (layer 0), channels coloured, and a legend that names the pattern groups by their FIT faces
    and gives the groups each channel vetoes."""
    figure = plt.figure(figsize=(17, 5.5 * len(namedMaps)), constrained_layout=True)
    for row, (name, channelMap) in enumerate(namedMaps.items()):
        colourOf = {(face, tileId): pu.categoricalColours[channelMap.colGroup[column]]
                    for column, (face, layer, tileId) in enumerate(tileKeys) if layer == 0}
        for column, azimuth in enumerate(views):
            pu.plotPsdCube(figure.add_subplot(len(namedMaps), 3, 3 * row + column + 1, projection="3d"), colourOf, azimuth=azimuth,
                           title=f"seen from azimuth {azimuth}°")
        groupLines = [Patch(facecolor="none", edgecolor="none",
                            label=f"G{k}: {faceSummary(fitFace[channelMap.rowGroup == k])} ({int((channelMap.rowGroup == k).sum())} patterns)")
                      for k in range(channelMap.numberOfGroups)]
        channelLines = [Patch(facecolor=pu.categoricalColours[l], label=f"C{l} ({int((channelMap.colGroup == l).sum())} tiles) vetoes "
                              + " ".join(f"G{k}" for k in np.flatnonzero(channelMap.routing[:, l])))
                        for l in range(channelMap.numberOfChannels)]
        legendAxis = figure.add_subplot(len(namedMaps), 3, 3 * row + 3)
        legendAxis.axis("off")
        legendAxis.legend(handles=groupLines + channelLines, loc="center left", fontsize=10, frameon=False, title=name, handlelength=1.2)
    figure.savefig(path, dpi=100)
    plt.close(figure)


def main():
    pu.setStyle()
    blues = LinearSegmentedColormap.from_list("blues", pu.sequentialRamp)
    patterns = pd.selectPatterns(pd.parsePatternDump(pd.defaultDumpPath), **pd.triggerSelections["thesis"])
    incidence = pd.buildIncidence(patterns)
    matrix, tileKeys = incidence.matrix, incidence.tileKeys
    fitFace = np.array([p.fitSide for p in patterns])
    faceIndex = np.array([pd.psdSides.index(face) for face in fitFace])

    rates = dict(zip(tileKeys, ru.loadTileRates(tileKeys)))
    figure, axes = plt.subplots(1, 2, figsize=(11, 5.2), constrained_layout=True)
    norm = Normalize(min(rates.values()), max(rates.values()))
    for layer in (0, 1):
        collection = pu.plotUnfoldedPsd(axes[layer], rates, layer, norm, blues, f"activation rate, layer {layer}")
    figure.colorbar(collection, ax=axes, shrink=0.8, label=r"$r_t$ [s$^{-1}$ cm$^{-2}$]")
    figure.savefig(figureDirectory / "psdRates.png", dpi=110)
    plt.close(figure)

    optimised = pickle.loads((resultsDirectory / "optimisedGrid.pkl").read_bytes())
    channelBudgets = sorted({c for _, c in optimised}, reverse=True)
    groupBudgets = sorted({g for g, _ in optimised}, reverse=True)
    pairBlocks = pd.evenOddBlocks(patterns)
    excess = lambda channelMap: mu.computePairMetrics(matrix, channelMap.vetoMask(), pairBlocks)["excess"]
    strategies = {"full PSD veto": lu.ChannelMap(np.zeros(matrix.shape[0], int), np.zeros(matrix.shape[1], int), np.ones((1, 1), bool)),
                  "face + 2 rows": mu.channelMapFromMasks(matrix, faceIndex, np.array([gu.faceExtendedMask(tileKeys, face, 2) for face in pd.psdSides])),
                  "face, union of RoIs": mu.channelMapFromMasks(matrix, faceIndex, np.array([matrix[faceIndex == k].any(axis=0) for k in range(5)]))}
    markers = {"full PSD veto": ("s", (7, 5)), "face + 2 rows": ("D", (7, 6)), "face, union of RoIs": ("o", (7, -14))}
    figure, drawingAxis = plt.subplots(figsize=(8, 5.4), constrained_layout=True)
    for colour, groups in zip(pu.ordinalRamp, groupBudgets):
        drawingAxis.plot(channelBudgets, [excess(optimised[(groups, channels)]) for channels in channelBudgets], "-", color=colour,
                         label=f"$N_G$ = {groups}" + (" ($N_G^*$)" if groups == groupBudgets[0] else ""))
    for name, (marker, offset) in markers.items():
        channelMap = strategies[name]
        position = (channelMap.numberOfChannels, excess(channelMap))
        drawingAxis.plot(*position, marker, color=pu.inkPrimary, markersize=8, markerfacecolor="none" if marker == "o" else pu.inkPrimary)
        drawingAxis.annotate(f"{name} ($N_G$={channelMap.numberOfGroups})", position, xytext=offset, textcoords="offset points", fontsize=9, color=pu.inkPrimary)
    drawingAxis.set(xscale="log", yscale="log", xlabel="PSD channels $N_C$", ylabel=r"event-RoI excess $\xi_{\rm ev}$", title="Extra tiles in the event veto regions")
    drawingAxis.legend(title="pattern groups", loc="lower left")
    figure.savefig(figureDirectory / "excessVsChannels.png", dpi=110)
    plt.close(figure)

    constrained = pickle.loads((resultsDirectory / "constrained.pkl").read_bytes())
    drawMapsIn3d(figureDirectory / "constrainedMap3d.png",
                 {"one group per FIT face, each channel inside one PSD face": constrained[((5, 8), "groups per FIT face, channels per PSD face")][0]},
                 tileKeys, fitFace)
    chipRuns = pickle.loads((resultsDirectory / "chipLayout.pkl").read_bytes())
    settings = chipRuns["settings"]
    cost = ru.costMatrix(matrix.shape[0], ru.tileVetoProbability(tileKeys, window=100e-9))     # 100 ns: Placeholder window, as in the notebook.
    reference = float((cost * matrix).sum())
    neighbours = cl.tileNeighbours(tileKeys)
    budget = sum(cl.minimumChips(tileKeys, settings["capacity"]).values())
    candidates = {key: run for key, run in chipRuns["runs"].items() if run["chips"] <= budget and key[0] != "temperature from all moves"}
    bestKey = min(candidates, key=lambda key: candidates[key]["cost"])
    best, _, _ = cl.polishLayout(matrix, cost, tileKeys, candidates[bestKey]["layout"], settings["groups"], settings["boards"], settings["capacity"],
                                 settings["borderPenalty"] * reference, True, budget, neighbours=neighbours)
    veto = mu.computeMaskCost(matrix, best.channelMap(matrix), cost)[0] / reference
    figure, axes = plt.subplots(1, 2, figsize=(14, 6.6), constrained_layout=True)
    pu.plotChipLayout(axes, best, tileKeys, neighbours)
    figure.suptitle(f"(N_G, N_C) = ({settings['groups']}, {settings['boards']}): {budget} chips of {settings['capacity']} inputs; colour = board "
                    f"(veto channel), shade = chip; cost {veto:.3f} of sum C M", color=pu.inkPrimary)
    figure.savefig(figureDirectory / "chipLayout.png", dpi=100)
    plt.close(figure)
    comparison = pickle.loads((resultsDirectory / "jaccardComparison.pkl").read_bytes())
    drawMapsIn3d(figureDirectory / "jaccardMaps3d.png",
                 {"Jaccard merge, refined": comparison["refined"][("jaccard", (5, 8))]["map"],
                  "exact-cost greedy, refined": comparison["refined"][("greedy", (5, 8))]["map"]}, tileKeys, fitFace)


if __name__ == "__main__":
    main()
