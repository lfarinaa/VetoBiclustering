# data

All files describe the **compact-v2 HERD layout with the earlier PSD of 10 x 10 cm tiles** (14 x 14 on top, 13 x 9 per side, two layers). The later
PSD design with 5 x 30 cm tiles is a different detector: none of these files apply to it without new dumps and rates.

| File | What it is | Where it comes from |
|---|---|---|
| `patterns_allPlanes.txt` | **The RoI map used by default**: every 3IR pattern (start planes 0-9, all nine codes, 3820 patterns) with its PSD tiles and calorimeter cubes | Regenerated on 2026-10-15, see below |
| `patterns.txt` | The February 2023 RoI map (planes 8-9, the "Calo" configuration, 426 patterns) | `~/AeffCalc/protons_compact/FIT/output_digitize/input_job_file/patterns.txt`, unchanged; kept for reference |
| `psdTileRates.csv` | Activation rate per unit area of every PSD tile (s^-1 cm^-2), all species, 'NewData' flux model | Recovered from thesis Fig. 4.9 by `tools/extractPsdRates.py` |
| `psdTileBacksplash.csv` | Probability that each tile fires per L0-triggered fiducial photon (log-uniform and E^-2 weighting) | Measured on the surviving Geant4 photons by `tools/measureBacksplash.py` |
| `regeneration/` | Everything needed to regenerate `patterns_allPlanes.txt` | See below |

## The regenerated RoI map

- **Geometry**: `regeneration/geometry_compact.mac`, the card of the AeffCalc proton study (`protons_compact/FIT`, unchanged from March 2022 to the 2023 dump).
- **Simulation**: GGSPenny with HerdSoftware **0.3.2** (`regeneration/env_herdsoftware032.sh`: ROOT 6.22, Geant4 10.6.2, GGS 2.8.1, EventAnalysis 1.2.1), two photons
  (`regeneration/sim.mac`); ProduceRoIMap only needs the geometry stored in the file. Later versions (0.4.0, master) no longer accept the card's
  `psd/type tiles` and have other geometry defaults.
- **RoI code**: `legacy/aeffcalc/algos` (July 2024, the code behind the thesis results), built as a minimal plugin against 0.3.2 with two mechanical API
  changes (`regeneration/UlegGeomPsdRoi_hs032.patch`: `Line::Intersection` returns a point, `PsdGeoParams::Normal()`) and the include paths of 0.3.x;
  `regeneration/CMakeLists.txt`.
- **Job**: `regeneration/produceRoIMap_allPlanes.eaconf` (FIT digitiser with 12/10/8 mats, PSD geometric digitiser without aggregation, ProduceRoIMap with
  `startLayer 0`, no ignored codes, no calorimeter requirement). The trigger configurations are applied afterwards (`analytic/patternDump.py`,
  `triggerSelections`): the calorimeter requirement is `caloCubes > 0`, as ProduceRoIMap's `requireCalo`.
- **Checks**: the "calo8" selection of this dump equals a separate ProduceRoIMap run with `startLayer 8`, `ignoreIDs {544,546,564,566}`, `requireCalo true`;
  the anti-meteorite shield line of the geometry card changes no RoI.

**The 2023 dump is not reproduced**: with the same configuration the regenerated map has 412 patterns against 426, and 141 of the 398 shared patterns have
the same RoI (the regenerated RoIs are smaller; mean 260 -> 179 tiles for the differing ones). The February 2023 RoI code (`~/testSFTP/algos`) with
HerdSoftware 0.3.2 does not reproduce it either (238 of 412 identical). The HerdSoftware branch the cluster ran in February 2023 (0.3_bugfix) has no
geometry change after 0.3.2, and the geometry card is the same, so the difference is in the code or the simulation behind that dump, which could not
be recovered (FINDINGS 7).

## Format of the dumps

Per pattern: `Trigger pattern: <FullID>`, the total number of tiles, then for each PSD face `PSD side <face>: <n> tiles:` followed by
`Layer <i>: (<tile IDs>)` lines, then `Cubes in RoI: <n>`. FullID = FIT face, first plane, first mat (2 digits), code (one digit per plane: 5 + mat step).
Read with `analytic/patternDump.py`.
