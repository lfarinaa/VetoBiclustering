# Findings

The running log of results, and of the artifacts found along the way. Newest at the bottom.

## 1. The input: the ProduceRoIMap dump (2026-10-14)

- `data/patterns.txt` (from `~/AeffCalc/protons_compact/FIT/output_digitize/input_job_file/`, written 2023-02-13) holds 426 patterns: TOP 104, Xpos 78, Xneg 82, Ypos 80, Yneg 82, each with its PSD tiles by face and layer and the number of calorimeter cubes in the same RoI. Every tile of the PSD (1328) is in at least one RoI.
- It is the "Calo" configuration: codes 545, 554, 555, 556, 565 only (the two-step codes removed: zig-zags 546, 564 and steep diagonals 544, 566), a calorimeter cube required in the RoI, and **first planes 8 and 9 only**. The datacard next to it (`digitizeMC_subset.eaconf`, edited 2023-02-17, after the dump) says `startLayer 0`, so the dump was presumably made with `startLayer 8`, as the "Calo" trigger of the analysis card uses `layer_start 8`. To confirm before relying on the plane coverage.
- An RoI is large: 206 to 250 tiles on average depending on the FIT side, a sixth to a fifth of the PSD. On the top PSD a top-FIT pattern's RoI is a band of whole tile columns, since one view constrains one coordinate.

## 2. The 0%-loss case (2026-10-14)

*On the February 2023 dump (planes 8-9) and the earlier cost models; superseded for the numbers by sections 7-8.*

- With partitioning merges the lossless minimum is K* = distinct RoIs, L* = distinct tile signatures, reached together (proof in `analytic/losslessUtils.py`; checked exactly on every side by the tests). Over all sides K* = 259 of 426 patterns, L* = 1277 of 1328 tiles; per side see the README.
- **Tile merging is almost never lossless**: 4 of 548 tiles merge for the top FIT, 4 to 56 for the side FITs. Patterns merge much better (64 of 104 on top).
- The two views never share an RoI: the distinct RoIs of the even and odd planes add up to K* on every side (TOP 32 + 32 = 64).
- With fan-out (a tile in several OR channels) the minimum is the Boolean rank; the bounds are 40 to 64 on top and 105 to 259 overall (fooling set; the greedy cover did not beat min(K*, L*)). The gap is open.

## 3. A channel budget, greedy and with a placeholder cost (2026-10-14)

*On the February 2023 dump (planes 8-9) and the earlier cost models; superseded for the numbers by sections 7-8.*

- Union routing (a channel vetoes a group when any of its tiles is in any member's RoI) never leaks background; the cost is the (pattern, tile) pairs added to masks, relative to the pairs in the RoIs ("mask inflation"). This cost is geometric and treats every tile alike; it is not a photon efficiency.
- Top FIT, greedy (tile channels first, then groups): full tiles and K = 32 gives 6.5%; K = K* = 64 with L = 128 gives about 10%, L = 64 about 20%, L = 16 about 48%; K = 16, L = 64 gives 31.5%. Cutting the tiles costs more than cutting the patterns, as the lossless counts predict.

## 4. Issues in the legacy code (2026-10-14)

- **The veto is on the intersection of the two views.** `PsdRoiVeto` and `MultiRoiVeto` veto an event on the tiles in RoI(even-view pattern) ∩ RoI(odd-view pattern), not on one pattern's RoI. The lossless counts above still apply if the hardware takes the AND per tile channel, but an event-weighted cost has to be built on pattern pairs. The showcase uses single-view RoIs.
- `TriggerPattern::PointsAlmostUpwards` (`legacy/aeffcalc/algos/TriggerPattern.h:350`) assigns instead of comparing (`lastMatHeight = firstMatHeight`), so it returns true whenever the height is non-zero: `removeAlmostUpwardsPatterns` removes nearly every pattern. It is used by the "NoUpwardsPlus" configuration of the analysis card.
- `ProduceRoIMap`'s `patternsAndTiles` ntuple declares 9 variables but is filled with 10 (`ProduceRoIMap.cpp:80` and `:254`): the number of calorimeter cubes is silently dropped. The text dump is not affected.

## 5. Notation, geometry, and the thesis' strategies (2026-10-15)

*On the February 2023 dump (planes 8-9) and the earlier cost models; superseded for the numbers by sections 7-8.*

- **Notation** follows the thesis ($\{T\}_i$, even/odd patterns $i$, $j$, Eq. (4.2), $\varepsilon_{\rm L1}$). The budgets $K$, $L$ of sections 2-3 are renamed $N_G$, $N_C$ (L collides with the trigger levels L0-L2); $N_G^*$, $N_C^*$ are the old $K^*$, $L^*$.
- **PSD layout inferred from the dump**: every side face meets TOP along v = 0; on TOP, r grows towards +X and c towards +Y; on a side face of normal n, h grows along (-n_y, n_x). Under this layout the RoIs of each FIT face reach only the adjacent PSD faces, and at most rows 0-2 into them (test `test_inferred_geometry_matches_the_dump`). A first guess with h reversed made the side RoIs reach 12 rows into the neighbours, which is how the error showed.
- **Event RoI**: a single-view RoI holds 232 tiles on average, an even-odd intersection 140 (40 to 289); no pair of the dump has an empty intersection.
- **The thesis' hybrid** (face + 2 rows of the adjacent faces) is 5 groups and **17** channels (the overlapping face masks split the tiles into 17 signatures). Per view it misses 0.9% of the RoI entries (side RoIs reach a third row), at event level none. Three rows, or the union of each face's RoIs, are leak-free per view too. Event excess 2.09 (union: 2.06), relative $\varepsilon_{\rm L1}$ 0.982.
- **Calibration of the toy efficiency**: a per-tile backsplash probability p = 6.4e-5 (uniform, independent tiles) makes the full-PSD veto 0.89/0.96 of the RoI veto, as in Table 4.2 (Full PSD veto against All Sides Clean). An energy average; a placeholder until Geant4 gives per-tile rates.

## 6. Optimised maps against the budget (2026-10-15)

*On the February 2023 dump (planes 8-9) and the earlier cost models; superseded for the numbers by sections 7-8.*

- **Search**: greedy merges from $(N_G^*, N_C^*)$ in three orders (channels first, groups first, interleaved), seeds, and local single-move refinement. Merging all channels before any group is a poor order at small $N_G$: at (5, 17) it gives a per-pattern excess of 1.55 where groups-first reaches 0.81. Every order is now tried and the best kept.
- **At the thesis' budget (5, 17) nothing beats the union of each face's RoIs**, with or without the constraint that groups stay within a FIT face: the thesis' layout is a local optimum of its budget. Beyond it, groups are the bottleneck (N_G = 5 is flat from N_C = 17 up); at N_C ≤ 8 channels are (more than 16 groups change nothing).
- Relative $\varepsilon_{\rm L1}$ of optimised maps (calibrated p): 0.982 at (5, 17), 0.990 at (16, 64), 0.995 at (64, 64), 0.948 at N_C = 2 whatever N_G.
- **Channels confined to one PSD face** (one board per face) raise the per-pattern excess at (5, 17) from 0.81 to 1.37 and at (5, 8) from 1.52 to 2.71.
- A first version of the constrained search let a seed carry channels across faces; seeds are now cut along the labels before the search (test `test_labels_constrain_every_search`).

## 7. Regenerated RoIs, measured rates and backsplash, events as the figure of merit (2026-10-15)

- **RoI map regenerated for every start plane** (data/README.md): 3820 patterns. HerdSoftware 0.3.2 is the latest local version that accepts the compact-v2 card (`psd/type tiles` is gone in 0.4.0, which also changed geometry defaults); the 2024 RoI code builds against it with two mechanical API changes.
- **The February 2023 dump is not reproduced** by either the 2024 or the February 2023 RoI code on 0.3.2 with the same geometry card and configuration (141 and 238 of about 400 RoIs identical; the regenerated RoIs are strict subsets, about 80 tiles smaller where they differ). The cluster's 0.3_bugfix branch had no geometry change after 0.3.2, so the 2023 dump most likely came from different code or a different simulation; it is kept but no longer used.
- **Trigger selection**: the two-step codes are of two kinds: zig-zags 546, 564 (back to the first mat) and steep diagonals 544, 566 (third hit two mats away). Upward-pointing patterns (TriggerPattern::PointsUpwards; upgoing: the particle moves towards the top) are the even-plane side patterns with a net negative mat step: their RoIs sit at the bottom of the face (mean row 7.5 of 0-8 for 544) and never reach TOP, while net positive steps reach TOP in 98% of cases. Six even-plane side 544 patterns have empty RoIs (they extrapolate below the PSD): no tile could veto them; every selection without two-step codes removes them. The thesis selection (CALO intersect, third hit adjacent, no upward-pointing) keeps 1780 patterns: N_G* = 1126, N_C* = 1326 of 1328.
- **PSD rates recovered** from the published thesis Fig. 4.9 (the events are gone): 1328 tiles, mean 1.880, std 0.130, max 2.146 s^-1 cm^-2 against 1.889, 0.126, 2.150 from the events (AeffCalc max_rates.txt). The figure was drawn with an earlier unfolding (Xneg and Yneg one cell lower than the current psd_rates.py). Its tiles are the 10 x 10 cm compact-v2 ones (tileArea 100 cm^2); the later 5 x 30 cm design is not modelled.
- **The uniform backsplash calibration on Table 4.2 is withdrawn.** The surviving Geant4 photons (2 x 100 000, 10 MeV-100 GeV log-uniform, upper hemisphere; 340 L0-triggered fiducial photons with the thesis selection) give the event-level L1 directly. Weighted as generated, the full-PSD veto keeps 0.50 and the exact RoI veto 0.83 of the fiducial photons; weighted to E^-2 they keep 0.90 and 0.99, and the 0.90 matches Table 4.2's 0.89, so the table is most likely spectrum-weighted. Backsplash grows from about 1 fired tile per photon at 30 MeV to 21 at 30-100 GeV; the measured per-tile probability averages 1.9e-4 (E^-2), three times the old uniform value.
- **Accidentals**: q_t = 1 - exp(-r_t A tau) is 1.9e-5 per tile in 100 ns and 1.9e-4 in 1 us, as large as the E^-2 backsplash. Mask excess therefore costs efficiency at all energies once the window approaches a microsecond.
- **Speed**: the greedy now keeps the minimum of each row of the merge-cost matrix, so a merge costs O(N) instead of O(N^2); identical results, 3x faster on the 2023 dump, and the 1780-pattern scans run in a few minutes.

## 8. Channel budgets on Geant4 photons (2026-10-15)

Thesis selection (1780 patterns), cost from the measured backsplash and accidentals in 100 ns, L1 efficiency measured on the 340 L0-triggered fiducial photons
weighted to E^-2 (examples/channelBudgetShowcase.ipynb; README for the tables).

- **Thesis strategies**: exact RoI 0.988, full-PSD veto 0.880, face + 2 rows 0.947 (5 groups, 17 channels), face union 0.949; face + 0 rows 0.970 but leaks 11.5% of the event-RoI tiles. With a 1 us window: 0.970, 0.710, 0.881, 0.884.
- **The two-row hybrid misses 0.4% of the single-view RoI tiles and none of the event-RoI tiles.**
- **Optimised maps**: at (5, 17) nothing beats the face union, with or without face constraints; with 5 groups nothing improves above 17 channels; (16, 64) 0.969, (64, 64) 0.973, (256, 64) 0.984; at N_C = 8 about 0.945 whatever N_G >= 5. The greedy is not monotonic in the budget in a few cells (N_G = 16: 0.966 at 1326 channels, 0.969 at 128), a sign of its local optima.
- **Channels confined to one PSD face**: 0.949 -> 0.945 at (5, 17), 0.943 -> 0.920 at (5, 8).
- **Energy dependence**: exact RoI 1.00 up to 0.3 GeV, 0.86 at 3-10 GeV, 0.45 at 30-100 GeV; full PSD 0.83 already at 0.1-0.3 GeV and 0.06 at 30-100 GeV.
- **The E^-2 backsplash map is dominated by a few low-energy photons** (a handful of hot tiles): the per-tile values are a rough estimate.
- **Run time**: the budget grid takes about 25 minutes, the constrained solutions a few minutes with the groups-first order and the per-face seed; the interleaved order is too slow at this size (it rebuilds a full table per merge) and is left to the CLI.

## 9. Correction: the thesis' layout is 5 channels with fan-out (2026-10-15)

- The thesis allows overlap: each of the 5 face groups ORs its own tile set (its face plus adjacent rows), and edge tiles are wired to several groups. "17 channels, not 5" (sections 5, 8 and the first README) counted channels in the partition model, where each tile's signal goes to one channel only; it is not a flaw of the thesis' layout. With fan-out the two-row hybrid is 5 OR channels, 2048 tile-to-group links, 624 tiles fanned out to 2-3 groups.
- With fan-out a map is a grouping of the patterns with union masks, i.e. the N_C = N_C* column of the budget grid: 5 groups 0.949 (the face union; two rows 0.947), 16 groups 0.966 (5120 links, fan-out up to 9), 64 groups 0.977 (14 953, up to 23), 256 groups 0.987 (50 663, up to 69). The hardware cost then sits in the wiring, which the search does not yet constrain (ROADMAP).
