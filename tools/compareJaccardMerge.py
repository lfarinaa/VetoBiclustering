"""Compares the Jaccard agglomeration (analytic/jaccardMerge.py) with the exact-cost greedy search (analytic/mergeUtils.py) on the HERD showcase, for
few channels: N_G in {64, 16, 5} pattern groups and N_C in {32, 17, 8, 4} PSD channels, cost from the measured backsplash and accidental hits in a
100 ns window (rateUtils, the default of tools/findChannelMap.py).

    python tools/compareJaccardMerge.py                 # about 20 minutes; results in examples/results/jaccardComparison.pkl, figures in docs/figures

For every budget it reports the exact false-veto cost, relative to the cost-weighted size of the RoIs (sum C M), of
- the exact-cost greedy maps, channels first and groups first, unrefined;
- the Jaccard maps, three similarity variants x the same two orders, unrefined;
- the best map of each family after local refinement (mergeUtils.refineLocally);
and the L1 efficiency on the Geant4 photons (E^-2 weighting, accidentals in 100 ns) of the refined ones.

Progress along the optimisation is checked at the (5, 8) budget in two ways: the exact cost after every merge (the trajectory of each method), and,
every few merges of the Jaccard path, the cost the Jaccard choice adds against the cheapest merge available from the same state (the regret).
With --plot-only the report and figure are redrawn from the saved results.
"""
import pathlib
import pickle
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root / "analytic")); sys.path.insert(0, str(root / "stochastic"))
import patternDump as pd
import mergeUtils as mu
import jaccardMerge as jm
import rateUtils as ru
import photonEvents as pe
import plotUtils as pu

groupBudgets, channelBudgets = [64, 16, 5], [32, 17, 8, 4]
orders = ("channelsFirst", "groupsFirst")
trajectoryBudget = (5, 8)
regretEvery = 25                                                           # Check the Jaccard choice every this many merges (each check rebuilds the exact tables).
resultsPath = root / "examples/results/jaccardComparison.pkl"
figureDirectory = root / "docs/figures"


def greedyTrajectory(matrix, cost, budget, order):
    """Exact cost after every merge of the exact-cost greedy path to `budget`, in merge order."""
    state = mu._GreedyState(matrix, cost, None, None)
    channelHistory, groupHistory = {}, {}
    state.mergeTo(*budget, order, groupHistory=groupHistory, channelHistory=channelHistory)
    channelSteps = [channelHistory[n] for n in sorted(channelHistory, reverse=True)]
    groupSteps = [groupHistory[n] for n in sorted(groupHistory, reverse=True)]
    return np.array(channelSteps + groupSteps if order == "channelsFirst" else groupSteps + channelSteps)


def jaccardRegret(matrix, cost, variant, order, budget):
    """Every `regretEvery` merges of the Jaccard path: (step, side, cost added by the Jaccard choice, cheapest cost available from that state)."""
    records, step = [], [0]

    def inspect(state, side, keep, drop):
        step[0] += 1
        if step[0] % regretEvery:
            return
        current = state.channelMap()
        exact = mu._GreedyState(matrix, cost, None, None, start=current)
        cheapest = (exact.bestGroupMerge() if side == "group" else exact.bestChannelMerge())[0]
        trial = state.copy()
        (trial.applyGroupMerge if side == "group" else trial.applyChannelMerge)(keep, drop)
        records.append((step[0], side, trial.exactCost - state.exactCost, float(cheapest)))

    state = jm._JaccardState(matrix, np.asarray(cost, dtype=float), variant)
    state.mergeTo(*budget, order, inspect=inspect)
    return records


def main():
    patterns = pd.selectPatterns(pd.parsePatternDump(pd.defaultDumpPath), **pd.triggerSelections["thesis"])
    incidence = pd.buildIncidence(patterns)
    matrix, tileKeys = incidence.matrix, incidence.tileKeys
    cost = ru.costMatrix(matrix.shape[0], ru.tileVetoProbability(tileKeys, window=100e-9))   # 100 ns: Placeholder window (rateUtils.defaultWindow).
    reference = float((cost * matrix).sum())
    relativeCost = lambda channelMap: mu.computeMaskCost(matrix, channelMap, cost)[0] / reference

    results = dict(groupBudgets=groupBudgets, channelBudgets=channelBudgets, maps={}, seconds={}, trajectories={}, regret={})
    for order in orders:
        start = time.time()
        results["maps"][("greedy", order)] = mu.scanBudgets(matrix, channelBudgets, groupBudgets, cost=cost, order=order)[0]
        results["seconds"][("greedy", order)] = time.time() - start
        print(f"greedy {order}: {results['seconds'][('greedy', order)]:.0f} s", flush=True)
        for variant in jm.variants:
            start = time.time()
            results["maps"][(variant, order)] = jm.scanJaccardBudgets(matrix, groupBudgets, channelBudgets, cost=cost, variant=variant, order=order)[0]
            results["seconds"][(variant, order)] = time.time() - start
            print(f"{variant} {order}: {results['seconds'][(variant, order)]:.0f} s", flush=True)

    results["cost"] = {method: {budget: relativeCost(channelMap) for budget, channelMap in maps.items()} for method, maps in results["maps"].items()}
    results["refined"] = {}
    for budget in [(g, c) for g in groupBudgets for c in channelBudgets]:
        for family, methods in (("greedy", [("greedy", order) for order in orders]),
                                ("jaccard", [(variant, order) for variant in jm.variants for order in orders])):
            bestMethod = min(methods, key=lambda method: results["cost"][method][budget])
            refined = mu.refineLocally(matrix, results["maps"][bestMethod][budget], cost)
            results["refined"][(family, budget)] = dict(seedMethod=bestMethod, map=refined, cost=relativeCost(refined))
        print(f"refined {budget}: greedy {results['refined'][('greedy', budget)]['cost']:.3f}, "
              f"jaccard {results['refined'][('jaccard', budget)]['cost']:.3f}", flush=True)

    triggered = pe.triggerPhotons(pe.loadPhotons(), patterns, tileKeys)
    accidental = ru.tileVetoProbability(tileKeys, window=100e-9, backsplash=0.0)
    for entry in results["refined"].values():
        entry["efficiency"] = pe.l1Efficiency(triggered, entry["map"].vetoMask(), accidentalProbability=accidental, spectralIndex=2.0)[0]

    for order in orders:
        results["trajectories"][("greedy", order)] = greedyTrajectory(matrix, cost, trajectoryBudget, order) / reference
        for variant in jm.variants:
            trace = jm.mergeByJaccard(matrix, *trajectoryBudget, cost=cost, variant=variant, order=order)[1]
            results["trajectories"][(variant, order)] = np.array([record["exactCost"] for record in trace]) / reference
    bestJaccard = min([(variant, order) for variant in jm.variants for order in orders], key=lambda method: results["cost"][method][trajectoryBudget])
    results["regretMethod"] = bestJaccard
    results["regret"] = jaccardRegret(matrix, cost, *bestJaccard, trajectoryBudget)

    resultsPath.write_bytes(pickle.dumps({key: value for key, value in results.items() if key != "maps"}))
    report(results)
    plot(results)


def report(results):
    methods = [("greedy", order) for order in orders] + [(variant, order) for variant in jm.variants for order in orders]
    budgets = [(g, c) for g in groupBudgets for c in channelBudgets]
    print("\nExact false-veto cost relative to sum C M (lower is better), unrefined:")
    print("| method | " + " | ".join(f"({g}, {c})" for g, c in budgets) + " |")
    print("|---|" + "---|" * len(budgets))
    for method in methods:
        print(f"| {method[0]}, {method[1]} | " + " | ".join(f"{results['cost'][method][budget]:.3f}" for budget in budgets) + " |")
    for family in ("greedy", "jaccard"):
        print(f"| best {family}, refined | " + " | ".join(f"{results['refined'][(family, budget)]['cost']:.3f}" for budget in budgets) + " |")
    for family in ("greedy", "jaccard"):
        print(f"| eps_L1 E^-2, best {family} refined | " + " | ".join(f"{results['refined'][(family, budget)]['efficiency']:.3f}" for budget in budgets) + " |")
    print("seconds per scan:", {f"{m[0]} {m[1]}": round(s) for m, s in results["seconds"].items()})
    for side in ("group", "channel"):
        checks = np.array([(added, cheapest) for _, recordSide, added, cheapest in results["regret"] if recordSide == side]).reshape(-1, 2)
        costly = checks[checks[:, 1] > 1e-12]
        if not costly.size:
            print(f"Regret, {side} merges: {len(checks)} checks, every one free for both the Jaccard choice and the cheapest merge")
            continue
        ratio = costly[:, 0] / costly[:, 1]
        print(f"Regret, {side} merges along {results['regretMethod']} to {trajectoryBudget}: {len(costly)} checks with a cost; the Jaccard choice was "
              f"the cheapest in {np.mean(ratio <= 1 + 1e-9):.0%}; added/cheapest quartiles {np.percentile(ratio, [25, 50, 75]).round(2)}; summed "
              f"{costly[:, 0].sum():.3g} against {costly[:, 1].sum():.3g}")


def plot(results):
    pu.setStyle()
    figure, (trajectoryAxis, regretAxis) = plt.subplots(1, 2, figsize=(15, 5.2), constrained_layout=True)
    colours = dict(zip(["greedy"] + list(jm.variants), pu.categoricalColours))
    for (method, order), values in results["trajectories"].items():
        trajectoryAxis.plot(np.arange(1, values.size + 1), values, color=colours[method], linestyle="-" if order == "groupsFirst" else "--",
                            linewidth=2, label=f"{'exact-cost greedy' if method == 'greedy' else method}, {order}")
    trajectoryAxis.set_yscale("symlog", linthresh=0.01)
    trajectoryAxis.set(xlabel="merges from the lossless map", ylabel=r"exact false-veto cost / $\sum C M$",
                       title=f"Exact cost along each merge path to (N_G, N_C) = {trajectoryBudget}")
    trajectoryAxis.legend(fontsize=8, loc="upper left")
    costly = [(step, added / cheapest) for step, side, added, cheapest in results["regret"] if cheapest > 1e-12]
    steps, ratios = np.array(costly).T
    regretAxis.plot(steps, ratios, "o", color=colours[results["regretMethod"][0]], markersize=6)
    regretAxis.axhline(1.0, color=pu.inkPrimary, linewidth=1)
    regretAxis.annotate("the cheapest merge available", (steps.min(), 1.0), xytext=(0, 6), textcoords="offset points", fontsize=9, color=pu.inkPrimary)
    regretAxis.set(yscale="log", xlabel="merges from the lossless map (group phase)", ylabel="cost added by the Jaccard choice / cheapest available",
                   title=f"Jaccard choice / cheapest merge, {results['regretMethod'][0]}, {results['regretMethod'][1]}")
    figure.savefig(figureDirectory / "jaccardTrajectories.png", dpi=110)
    plt.close(figure)


if __name__ == "__main__":
    if "--plot-only" in sys.argv:                                          # Redraw from the saved results without rerunning.
        saved = pickle.loads(resultsPath.read_bytes())
        report(saved); plot(saved)
    else:
        main()
