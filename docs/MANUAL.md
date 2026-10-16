# User manual

How to run VetoBiclustering and where everything is. The core (`losslessUtils`, `mergeUtils`, `tools/solveChannelMap.py`) is detector-agnostic; the rest is the HERD showcase, in the design of the thesis and paper. For the problem and the headline numbers see the [README](../README.md); for results, [FINDINGS.md](../FINDINGS.md).

## 1. Install

```
pip install -r requirements.txt        # numpy, scipy, matplotlib, jupyterlab, nbformat, nbconvert, nbstripout (jax and optax for the planned relaxed optimisation)
nbstripout --install                   # once, so that notebook outputs are not committed
```

## 2. Layout

| Path | What it is |
|---|---|
| [`analytic/`](../analytic) | Reading the RoI map, the lossless channel maps, Boolean-rank bounds, the search for channel maps under a budget, the PSD layout, plots |
| [`stochastic/`](../stochastic) | Event level: `photonEvents.py` redoes L0 and L1 on the surviving Geant4 photons and measures the L1 efficiency of any channel map; the iterative and relaxed optimisations are planned (its README, ROADMAP) |
| [`examples/`](../examples) | Notebooks; each runs top to bottom and puts `../analytic` on its path |
| [`data/`](../data) | The RoI maps (regenerated: `patterns_allPlanes.txt`; February 2023: `patterns.txt`), the PSD tile rates, the measured backsplash, the regeneration recipe; provenance in its README |
| [`legacy/aeffcalc/`](../legacy/aeffcalc) | The C++ algorithms of the AeffCalc HERD analysis that define patterns and RoIs, with the datacards that configured them; reference only, not built here |
| [`tests/`](../tests) | Sanity checks (`python tests/test_sanityChecks.py`) |
| [`tools/solveChannelMap.py`](../tools/solveChannelMap.py) | The detector-agnostic program: any incidence matrix (`.npy` or 0/1 text), optional cost, row/column labels and names; lists the map as text and JSON |
| [`tools/findChannelMap.py`](../tools/findChannelMap.py) | The HERD showcase front end: find (or evaluate) a channel map for a budget, trigger selection, coincidence window and constraints, and list it (text, JSON); `--events` adds the measured L1 efficiency |
| [`tools/extractPsdRates.py`](../tools/extractPsdRates.py), [`tools/measureBacksplash.py`](../tools/measureBacksplash.py) | Rebuild `data/psdTileRates.csv` from thesis Fig. 4.9 and `data/psdTileBacksplash.csv` from the Geant4 photons |
| [`tools/checkCommitTrailers.py`](../tools/checkCommitTrailers.py) | The commit-message check for AI credit (see CONTRIBUTING.md) |

### The modules of `analytic/`

| Module | Contents |
|---|---|
| `patternDump.py` | `parsePatternDump` (a ProduceRoIMap dump → `TriggerPattern` records), `selectPatterns` and `triggerSelections` (the trigger configurations: CALO intersect, two-step codes, `pointsUpwards`, planes), `buildIncidence` (the matrix $M$), `evenOddBlocks` (the even/odd pairs of each double layer, Eq. 4.2) |
| `rateUtils.py` | Per-tile veto probability of a good photon: measured backsplash $b_t$ and accidentals $q_t = 1-e^{-r_t A\tau}$; the cost matrix $C_{it}=-\ln(1-\pi_t)$ |
| `losslessUtils.py` | `ChannelMap` (groups $g$, channels $c$, routing $S$, the mask $V$), `findLosslessChannelMap` ($N_G^*$, $N_C^*$), `computeUnionRouting`, fooling-set and greedy-cover bounds on the Boolean rank |
| `mergeUtils.py` | `findChannelMap` (the search: greedy in three orders + seeds + local refinement, optional labels a merge may not cross), `mergeGreedily`, `scanBudgets` (a grid of budgets), `refineLocally`, `computeMaskCost` ($\xi$, $\kappa$), `computePairMetrics` (event-level $\xi_{\rm ev}$, $\kappa_{\rm ev}$, relative $\varepsilon_{\rm L1}$), `channelMapFromMasks`, `describeChannelMap` |
| `geometryUtils.py` | The schematic PSD layout (adjacency inferred from the dump), `faceExtendedMask` (the thesis' face + n rows veto), tile quads for 3D sketches |
| `plotUtils.py` | Style; RoI per face; ordered incidence matrix and RoI overlap; the unfolded PSD of thesis Fig. 4.9; 3D cube sketches |

## 3. Quick start

```python
import sys; sys.path.insert(0, "analytic"); sys.path.insert(0, "stochastic")
import numpy as np, patternDump as pd, losslessUtils as lu, mergeUtils as mu, rateUtils as ru, photonEvents as pe

patterns = pd.selectPatterns(pd.parsePatternDump(), **pd.triggerSelections["thesis"])   # 1780 patterns
incidence = pd.buildIncidence(patterns); M, keys = incidence.matrix, incidence.tileKeys
exact = lu.findLosslessChannelMap(M)                                # N_G* = 1126, N_C* = 1326: reproduces every RoI
cost = ru.costMatrix(len(patterns), ru.tileVetoProbability(keys, window=100e-9))
best, tried = mu.findChannelMap(M, 16, 64, cost=cost)               # 16 groups, 64 channels
triggered = pe.triggerPhotons(pe.loadPhotons(), patterns, keys)     # Geant4 photons (needs uproot and the AeffCalc files)
print(pe.l1Efficiency(triggered, best.vetoMask(), spectralIndex=2.0)[0])
print(mu.describeChannelMap(best, patterns, keys, M))
```

Or from the shell: `python tools/findChannelMap.py --groups 16 --channels 64 --json map.json` (see `--help` for the constraints and strategies).
To use real costs, pass `cost=` (an array of the matrix's shape: expected false vetoes of good photons if tile $t$ joins the mask of pattern $i$) or `--cost file.npy`.

## 4. Notebooks

| Notebook | What it shows | Run time |
|---|---|---|
| [`vetoPatternShowcase.ipynb`](../examples/vetoPatternShowcase.ipynb) | A first look on the February 2023 dump: RoIs on the PSD, the ordered incidence matrix and RoI overlaps, RoIs per tile on the unfolded PSD, the lossless counts, a first look at budgets for the top FIT | about 1.5 min |
| [`channelBudgetShowcase.ipynb`](../examples/channelBudgetShowcase.ipynb) | **The main one**: notation and mathematics, the regenerated data, rates and backsplash, 3D sketches of the event RoI, the thesis' strategies on Geant4 photons, L1 efficiency against the number of channels, solutions constrained to a per-face layout | about 30 min; a minute with `BICLUSTERING_USE_CACHE=1` |

Run one with `jupyter nbconvert --to notebook --execute --inplace examples/<name>.ipynb`; the exit code is 0 when every cell ran.

## 5. The legacy code

`legacy/aeffcalc/algos/` holds the HERD EventAnalysis algorithms (Valerio Formato, Nicola Mori, Luis Fariña, 2019-2024) as copied from `~/AeffCalc/algos`:

| File | Role |
|---|---|
| `TriggerPattern.h` | The pattern: side, first plane, mats, holes; the recursive pattern finder; the ID scheme (`5 + mat step` per plane) |
| `UlegGeomPsdRoi.h` | The RoI geometry (`UlegGeomRoi`: double cone through the first and last mat), the tile and cube lists (`PsdVetoList`, `CaloVetoList`) and the pattern → list cache (`PatternVetoListMap`) |
| `ProduceRoIMap.{h,cpp}` | Enumerates every pattern with all mats hit and dumps its RoI (`dumpPatterns`): the source of `data/patterns.txt` |
| `TriggerArrays.{h,cpp}` | The trigger on simulated events: hit arrays, patterns, filters (ignored IDs, calorimeter in the RoI, upward patterns), first triggered planes |
| `PsdRoiVeto.{h,cpp}`, `MultiRoiVeto.{h,cpp}` | The veto on events: the intersection of the even-view and odd-view RoIs of the first triggers, and the energy threshold |
| `NewPsdTileAggregator.cpp` | The PSD geometric digitiser with tile aggregation factors: tile merging as it exists in HerdSoftware |

`legacy/aeffcalc/datacards/` holds the geometry (`geometry_compact.mac`), the digitisation card that runs ProduceRoIMap (commented out there) and the analysis card that chains the
trigger and veto configurations ("Base", "NoPlusMinus", "Calo", ...). `PatternTools.h` of AeffCalc is an older copy of `TriggerPattern.h` and was not copied.

## 6. Conventions

- camelCase everywhere, no unit suffixes; every number in the code carries its source or the word "Placeholder".
- Tolerances and counts checked by the tests are statements about the dump: change one only with a reason in FINDINGS.md.
- AI assistance is credited with an `Assisted-by: <harness>:<model>` trailer, never with a co-author or sign-off line; see [CONTRIBUTING.md](../CONTRIBUTING.md).
