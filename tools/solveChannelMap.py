"""The detector-agnostic entry point: a channel map for any trigger-pattern x veto-element incidence matrix, with no HERD code involved.

    python tools/solveChannelMap.py --matrix roi.npy --groups 8 --channels 16
    python tools/solveChannelMap.py --matrix roi.csv --cost cost.npy --row-labels rowFaces.txt --col-labels colBoards.txt \
                                    --groups 5 --channels 8 --json map.json

Inputs: --matrix, the incidence matrix M (rows: trigger patterns, columns: veto elements such as tiles or strips; 1 when the element is in the pattern's
veto region), as .npy or as a dense 0/1 text file (comma- or space-separated); --cost, optional, the cost C+ of adding element t to the veto of pattern i
(same shape, .npy; default 1: the cost counts added entries); --row-labels / --col-labels, optional text files with one label per line: a group or a
channel may not mix rows or columns with different labels (one board per detector face, for example); --row-names / --col-names, optional names
for the listing. The search is the no-leak limit (union routing: no veto-region entry is ever dropped), greedy merges in three orders refined by
local moves: a good local optimum, not a proof of optimality. HERD, the showcase, has its own front end with geometry and rates: findChannelMap.py.
"""
import argparse
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analytic"))
import losslessUtils as lu
import mergeUtils as mu


def loadMatrix(path):
    if path.suffix == ".npy":
        return np.load(path).astype(bool)
    text = path.read_text()
    return np.loadtxt(path, delimiter="," if "," in text else None, ndmin=2).astype(bool)


def loadLines(path, expected, what):
    if path is None:
        return None
    lines = np.array([line.strip() for line in path.read_text().splitlines() if line.strip()])
    if lines.size != expected:
        raise SystemExit(f"{path}: {lines.size} {what}, the matrix has {expected}")
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--matrix", type=pathlib.Path, required=True, help="incidence matrix M (.npy, or dense 0/1 text)")
    parser.add_argument("--groups", type=int, help="N_G, the number of pattern groups (default: the lossless N_G*)")
    parser.add_argument("--channels", type=int, help="N_C, the number of veto channels (default: the lossless N_C*)")
    parser.add_argument("--cost", type=pathlib.Path, help="cost of each added entry, same shape as M (.npy); default 1")
    parser.add_argument("--row-labels", type=pathlib.Path, help="one label per row; a group may not mix labels")
    parser.add_argument("--col-labels", type=pathlib.Path, help="one label per column; a channel may not mix labels")
    parser.add_argument("--row-names", type=pathlib.Path, help="one name per row, for the listing")
    parser.add_argument("--col-names", type=pathlib.Path, help="one name per column, for the listing")
    parser.add_argument("--no-refine", action="store_true", help="skip the local refinement")
    parser.add_argument("--json", type=pathlib.Path, help="write the map as JSON (groups, channels, routing, metrics)")
    args = parser.parse_args(argv)

    matrix = loadMatrix(args.matrix)
    cost = None if args.cost is None else np.load(args.cost).astype(float)
    if cost is not None and cost.shape != matrix.shape:
        raise SystemExit(f"{args.cost}: shape {cost.shape}, the matrix is {matrix.shape}")
    rowLabels = loadLines(args.row_labels, matrix.shape[0], "row labels")
    colLabels = loadLines(args.col_labels, matrix.shape[1], "column labels")
    rowNames = loadLines(args.row_names, matrix.shape[0], "row names")
    colNames = loadLines(args.col_names, matrix.shape[1], "column names")
    rowNames = np.arange(matrix.shape[0]).astype(str) if rowNames is None else rowNames
    colNames = np.arange(matrix.shape[1]).astype(str) if colNames is None else colNames

    exact = lu.findLosslessChannelMap(matrix)
    lossless = (exact.numberOfGroups, exact.numberOfChannels)
    groups = lossless[0] if args.groups is None else args.groups
    channels = lossless[1] if args.channels is None else args.channels
    best, candidates = mu.findChannelMap(matrix, groups, channels, cost=cost, rowLabels=rowLabels, colLabels=colLabels, refine=not args.no_refine)
    added, leaked = mu.computeMaskCost(matrix, best, cost)
    addedEntries = mu.computeMaskCost(matrix, best)[0]

    print(f"M: {matrix.shape[0]} patterns x {matrix.shape[1]} veto elements, {int(matrix.sum())} entries; lossless N_G* = {lossless[0]}, N_C* = {lossless[1]}")
    print(f"budget ({groups}, {channels}) -> map ({best.numberOfGroups}, {best.numberOfChannels}); added cost {added:.4g}, "
          f"added entries {int(addedEntries)} ({addedEntries / max(matrix.sum(), 1):.2%} of M), dropped entries {leaked}")
    print("candidates: " + ", ".join(f"{name} {cost_:.4g}" for name, cost_ in candidates))
    groupRows = [rowNames[best.rowGroup == k].tolist() for k in range(best.numberOfGroups)]
    channelColumns = [colNames[best.colGroup == l].tolist() for l in range(best.numberOfChannels)]
    for l, columns in enumerate(channelColumns):
        print(f"C{l}: {len(columns)} elements: {' '.join(columns[:20])}{' ...' if len(columns) > 20 else ''}")
    for k, rows in enumerate(groupRows):
        print(f"G{k}: {len(rows)} patterns, vetoed by " + " ".join(f"C{l}" for l in np.flatnonzero(best.routing[k])))

    if args.json:
        args.json.write_text(json.dumps(dict(
            shape=list(matrix.shape), lossless=dict(groups=lossless[0], channels=lossless[1]), budget=dict(groups=groups, channels=channels),
            groups=groupRows, channels=channelColumns, routing=best.routing.astype(int).tolist(),
            metrics=dict(addedCost=added, addedEntries=int(addedEntries), droppedEntries=leaked)), indent=1))
    return best


if __name__ == "__main__":
    main()
