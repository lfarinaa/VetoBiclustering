# legacy/aeffcalc

C++ algorithms and datacards of the AeffCalc HERD analysis (HerdSoftware / EventAnalysis plugin "Gamma"), copied on 2026-10-14 from `~/AeffCalc`, unchanged. They define the trigger
patterns and their PSD veto regions; they are kept as the reference for the Python code of this project and are not built here (they need the HerdSoftware data objects).

| Copied from | Files |
|---|---|
| `~/AeffCalc/algos/` | `TriggerPattern.h`, `UlegGeomPsdRoi.h`, `ProduceRoIMap.{h,cpp}`, `TriggerArrays.{h,cpp}`, `PsdRoiVeto.{h,cpp}`, `MultiRoiVeto.{h,cpp}`, `NewPsdTileAggregator.cpp` |
| `~/AeffCalc/oneshot_protons_compact/conffile_geometry/model/` | `datacards/geometry_compact.mac` |
| `~/AeffCalc/oneshot_gammas_compact/conffile_digitize/model/` | `datacards/digitizeMC_subset.eaconf` (ProduceRoIMap is commented out in it) |
| `~/AeffCalc/oneshot_gammas_compact/conffile_analysis/model/` | `datacards/analyzeMC_top.eaconf` (the trigger and veto configurations) |

What each file does: [docs/MANUAL.md](../../docs/MANUAL.md), section 5. Issues found in them: [FINDINGS.md](../../FINDINGS.md), section 4.
