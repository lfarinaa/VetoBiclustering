"""Sanity checks of the dump reader, the lossless channel map, the Boolean-rank bounds and the greedy merge.

Run all of them (about half a minute):   python tests/test_sanityChecks.py
With pytest, if installed:                pytest tests

Each check states what it compares. The numbers checked against the dump are those of data/patterns.txt; regenerate them in FINDINGS.md if the dump changes.
"""
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analytic"))
import patternDump as pd
import losslessUtils as lu
import mergeUtils as mu
import geometryUtils as gu
import rateUtils as ru

incidence = pd.buildIncidence(pd.parsePatternDump(pd.legacyDumpPath))   # The checks with fixed counts use the 2023 dump.


def test_dump_parses_completely():
    """426 patterns (78 + 82 + 80 + 82 + 104 by FIT side, the counts written in the dump), every tile ID inside its face, the per-pattern totals matching (checked by the parser)."""
    counts = {side: len(incidence.forFitSides(side).patterns) for side in pd.psdSides}
    assert counts == {"TOP": 104, "Xpos": 78, "Xneg": 82, "Ypos": 80, "Yneg": 82}, counts
    for pattern in incidence.patterns:
        assert all(0 <= tileId < pd.tilesPerLayer[side] and 0 <= layer < pd.psdLayers for side, layer, tileId in pattern.tiles)
        assert pattern.code[0] == "5" and all(abs(step) <= 1 for step in pattern.matSteps)


def test_dump_trigger_configuration():
    """The dump holds the 'Calo' configuration: no two-step codes (the zig-zags 546, 564 and the steep 544, 566), first planes 8 and 9 only."""
    assert {p.code for p in incidence.patterns}.isdisjoint({"544", "546", "564", "566"})
    assert {p.firstPlane for p in incidence.patterns} == {8, 9}


def test_lossless_map_is_exact_and_minimal():
    """The lossless map reproduces M exactly, its K and L are the numbers of distinct rows and columns, and no two groups (channels) share a mask (a column of S)."""
    for side in pd.psdSides:
        matrix = incidence.forFitSides(side).withoutEmptyColumns().matrix
        channelMap = lu.findLosslessChannelMap(matrix)
        assert np.array_equal(channelMap.vetoMask(), matrix)
        assert channelMap.numberOfGroups == np.unique(matrix, axis=0).shape[0]
        assert channelMap.numberOfChannels == np.unique(matrix, axis=1).shape[1]
        assert np.unique(channelMap.routing, axis=0).shape[0] == channelMap.numberOfGroups
        assert np.unique(channelMap.routing, axis=1).shape[1] == channelMap.numberOfChannels


def test_fooling_set_and_cover_are_valid():
    """The fooling set has no two entries in one all-ones rectangle, the greedy bicliques are all-ones and cover M exactly, so lower <= Boolean rank <= upper."""
    matrix = np.unique(np.unique(incidence.forFitSides("Xpos").withoutEmptyColumns().matrix, axis=0), axis=1)
    fooling = lu.computeFoolingSet(matrix)
    for a, (i, j) in enumerate(fooling):
        assert matrix[i, j]
        for k, l in fooling[a + 1:]:
            assert not (matrix[i, l] and matrix[k, j])
    covered = np.zeros_like(matrix)
    for rows, columns in lu.computeGreedyBicliqueCover(matrix):
        assert matrix[np.ix_(rows, columns)].all()
        covered[np.ix_(rows, columns)] = True
    assert np.array_equal(covered, matrix)
    lower, upper = lu.computeBooleanRankBounds(matrix)
    assert lower <= upper <= min(matrix.shape)


def test_small_boolean_rank_example():
    """A matrix of known Boolean rank 3 (the 3 x 3 complement of the identity has Boolean rank 3): the bounds bracket it."""
    matrix = ~np.eye(3, dtype=bool)
    lower, upper = lu.computeBooleanRankBounds(matrix)
    assert lower <= 3 <= upper


def test_union_routing_never_leaks_and_costs_are_consistent():
    """Greedy merges with union routing contain every RoI (zero leak), meet the budget, and the incrementally tracked cost equals the recomputed one."""
    matrix = incidence.forFitSides("Xneg").withoutEmptyColumns().matrix
    generator = np.random.default_rng(1)
    cost = generator.random(matrix.shape)
    channelMap, history = mu.mergeGreedily(matrix, 12, 40, cost=cost, record=True)
    added, leaked = mu.computeMaskCost(matrix, channelMap, cost)
    assert leaked == 0
    assert channelMap.numberOfGroups <= 12 and channelMap.numberOfChannels <= 40
    assert np.isclose(added, history[12])
    union = lu.computeUnionRouting(matrix, channelMap.rowGroup, channelMap.colGroup, *channelMap.routing.shape)
    assert np.array_equal(union.routing, channelMap.routing)


def test_inferred_geometry_matches_the_dump():
    """The RoIs of each FIT face reach only adjacent PSD faces, and at most 3 rows into them (rows 0-2), under the layout inferred in geometryUtils."""
    for side in pd.psdSides:
        used = incidence.forFitSides(side).matrix.any(axis=0)
        distances = [gu.rowsFromEdge(face, tileId, side) for (face, layer, tileId), u in zip(incidence.tileKeys, used) if u and face != side]
        assert None not in distances and max(distances) <= 2, (side, max(d for d in distances if d is not None))


def test_labels_constrain_every_search():
    """With groups confined to FIT faces and channels to PSD faces, no group and no channel crosses a face, whatever the seed."""
    matrix = incidence.matrix
    fit = np.array([p.fitSide for p in incidence.patterns]); psd = np.array([key[0] for key in incidence.tileKeys])
    faceIndex = np.array([pd.psdSides.index(f) for f in fit])
    seed = mu.channelMapFromMasks(matrix, faceIndex, np.array([matrix[faceIndex == k].any(axis=0) for k in range(5)]))   # its channels cross faces
    best, _ = mu.findChannelMap(matrix, 5, 8, rowLabels=fit, colLabels=psd, seeds=[("face", seed)], orders=("interleaved",))
    assert all(len(set(fit[best.rowGroup == k])) == 1 for k in range(best.routing.shape[0]))
    assert all(len(set(psd[best.colGroup == l])) == 1 for l in range(best.routing.shape[1]))
    assert mu.computeMaskCost(matrix, best)[1] == 0 and best.routing.shape[0] <= 5 and best.routing.shape[1] <= 8


def test_event_metrics_agree_with_explicit_sets():
    """computePairMetrics against explicit set intersections for a few pairs of the thesis hybrid (face + 2 rows), and the exact RoI has zero excess and leak."""
    matrix, blocks = incidence.matrix, pd.evenOddBlocks(incidence.patterns)
    faceIndex = np.array([pd.psdSides.index(p.fitSide) for p in incidence.patterns])
    hybrid = mu.channelMapFromMasks(matrix, faceIndex, np.array([gu.faceExtendedMask(incidence.tileKeys, face, 2) for face in pd.psdSides])).vetoMask()
    exactSum = keptSum = maskedSum = 0
    for even, odd in blocks:
        for i in even:
            for j in odd:
                t, v = matrix[i] & matrix[j], hybrid[i] & hybrid[j]
                exactSum += t.sum(); keptSum += (t & v).sum(); maskedSum += v.sum()
    metrics = mu.computePairMetrics(matrix, hybrid, blocks)
    assert np.isclose(metrics["excess"], (maskedSum - keptSum) / exactSum) and np.isclose(metrics["leak"], (exactSum - keptSum) / exactSum)
    exact = mu.computePairMetrics(matrix, matrix, blocks, 1e-4)
    assert exact["excess"] == 0 and exact["leak"] == 0 and np.isclose(exact["selfVeto"], 1.0)


def test_refinement_never_increases_the_cost():
    """Local refinement of a greedy map keeps union routing (no leak), the budgets, and does not raise the cost."""
    matrix = incidence.forFitSides("TOP").withoutEmptyColumns().matrix
    greedy = mu.mergeGreedily(matrix, 16, 64)
    refined = mu.refineLocally(matrix, greedy)
    assert mu.computeMaskCost(matrix, refined)[0] <= mu.computeMaskCost(matrix, greedy)[0] + 1e-9
    assert mu.computeMaskCost(matrix, refined)[1] == 0 and refined.routing.shape[0] <= 16 and refined.routing.shape[1] <= 64



def test_regenerated_dump_is_consistent():
    """The regenerated dump: 3820 patterns on planes 0-9 and all nine codes; its 'calo8' selection keeps the structure of the 2023 dump's configuration
    (planes 8-9, no two-step code, a calorimeter cube in every RoI); the only empty RoIs are six steep (544) patterns; the pairing gives one block per FIT face and double layer."""
    patterns = pd.parsePatternDump()
    assert len(patterns) == 3820 and {p.firstPlane for p in patterns} == set(range(10)) and len({p.code for p in patterns}) == 9
    calo8 = pd.selectPatterns(patterns, **pd.triggerSelections["calo8"])
    assert len(calo8) == 412 and all(p.caloCubes > 0 and p.firstPlane in (8, 9) and p.code not in pd.twoStepCodes for p in calo8)
    empty = [p for p in patterns if not p.tiles]                       # Patterns no tile can veto (FINDINGS 7): six steep 544s, removed by every selection that drops two-step codes.
    assert len(empty) == 6 and all(p.code == "544" for p in empty)
    assert all(p.tiles for p in pd.selectPatterns(patterns, **pd.triggerSelections["thesis"]))
    assert len(pd.evenOddBlocks(patterns)) == 5 * 5



def test_recovered_rates_match_the_event_summary():
    """The per-tile rates recovered from thesis Fig. 4.9 agree with the summary computed from the events (AeffCalc max_rates.txt: mean 1.889, std
    0.126, max 2.150 s^-1 cm^-2) within 1%, 5% and 0.5%; every tile of the layout has a rate."""
    rates = ru.loadTileRates(pd.allTileKeys())
    assert rates.size == 1328 and np.all(rates > 0)
    assert abs(rates.mean() / 1.889 - 1) < 0.01 and abs(rates.std() / 0.126 - 1) < 0.05 and abs(rates.max() / 2.150 - 1) < 0.005


def test_tile_probabilities_reduce_to_the_scalar_case():
    """With no window the per-tile veto probability is the uniform backsplash p, and the event metrics with a vector equal those with the scalar."""
    keys = incidence.tileKeys
    probability = ru.tileVetoProbability(keys, window=0.0, backsplash=7e-5)
    assert np.allclose(ru.tileVetoProbability(keys, window=0.0), ru.loadBacksplash(keys))
    assert np.allclose(probability, 7e-5)
    matrix, blocks = incidence.matrix, pd.evenOddBlocks(incidence.patterns)
    full = np.ones_like(matrix)
    assert np.isclose(mu.computePairMetrics(matrix, full, blocks, probability)["selfVeto"], mu.computePairMetrics(matrix, full, blocks, 7e-5)["selfVeto"])
    longer = ru.tileVetoProbability(keys, window=1e-6, backsplash=0.0)
    assert np.all(longer > ru.tileVetoProbability(keys, window=1e-7, backsplash=0.0))


def test_generic_front_end_on_a_toy_detector():
    """The detector-agnostic tool on a matrix with four planted blocks: N_G* = N_C* = 4, and a (4, 4) budget recovers the blocks at zero cost;
    with a column label that splits a block, no channel mixes labels and no veto-region entry is dropped."""
    import tempfile
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "tools"))
    import solveChannelMap
    matrix = np.zeros((12, 20), dtype=bool)
    for k in range(4):
        matrix[3 * k:3 * k + 3, 5 * k:5 * k + 5] = True
    with tempfile.TemporaryDirectory() as directory:
        path = pathlib.Path(directory) / "toy.npy"
        np.save(path, matrix)
        best = solveChannelMap.main(["--matrix", str(path), "--groups", "4", "--channels", "4"])
        assert (best.numberOfGroups, best.numberOfChannels) == (4, 4) and np.array_equal(best.vetoMask(), matrix)
        labels = pathlib.Path(directory) / "labels.txt"
        labels.write_text("\n".join(["a"] * 7 + ["b"] * 13))                # Block 1 (columns 5-9) now straddles two boards.
        constrained = solveChannelMap.main(["--matrix", str(path), "--groups", "4", "--channels", "4", "--col-labels", str(labels)])
        assert np.all(constrained.vetoMask() >= matrix)
        for l in range(constrained.numberOfChannels):
            assert len(set(np.array(["a"] * 7 + ["b"] * 13)[constrained.colGroup == l])) == 1


if __name__ == "__main__":
    failures = 0
    for name, check in [(name, value) for name, value in globals().items() if name.startswith("test_")]:
        start = time.time()
        try:
            check()
            print(f"PASS {name} ({time.time() - start:.1f} s)")
        except AssertionError as error:
            failures += 1
            print(f"FAIL {name}: {error}")
    sys.exit(1 if failures else 0)
