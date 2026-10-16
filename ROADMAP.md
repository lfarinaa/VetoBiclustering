# Roadmap

## Analytic (0% loss)

1. **Exact Boolean rank per FIT side**: an ILP (scipy's HiGHS) or SAT model of the minimum biclique cover on the deduplicated matrix (tens of rows, a few hundred columns), to close the gap between the fooling-set bound and N_G*.
2. ~~**The two-view veto**~~ done: the AND per channel realises the tile-level intersection exactly, and the figures of merit are evaluated on the even-odd pairs (FINDINGS 5).
3. **Hardware constraints of the trigger board**: groups per FIT face and channels per PSD face or layer are supported as labels (FINDINGS 6); still to add: channel counts per board, neighbour-only wiring, a fan-out limit and a link budget for the fan-out model (where a map is a grouping of the patterns and the wiring is the cost; FINDINGS 9), and an objective on the event RoI directly (the search minimises the per-pattern excess).
4. **All start planes**: a dump with `startLayer 0` and with the two-step codes (zig-zags 546, 564; steep 544, 566), to see how N_G* and N_C* grow with the trigger acceptance.

## Stochastic (event-weighted, iterative)

5. **Costs from events, and leaks allowed**: the false-veto cost C+_it = -w_i ln[(1 - b_it)(1 - q_t)] from Geant4 photons, with the backsplash as a pattern x tile map b_it (first-trigger patterns x PSD hits outside the RoI, above threshold; smoothed or parametrised by the distance from the pattern's track, since 1780 x 1328 entries need many photons) and the photon rate w_i of each pattern and the leak cost C-_it from protons and electrons; the AeffCalc samples first, a fresh production after. Then drop the no-leak limit: block routing by the weighted majority S_kl = [sum C- M >= sum C+ (1 - M)] in the merge costs and refinement, and sweep C-/C+ for the efficiency-rejection trade-off.
6. **Optimisers on the assignments**: simulated annealing / basin hopping over the group and channel labels, started from the greedy map, on training events, evaluated on test events.
7. **Relaxed optimisation**: soft assignments (softmax groups, sigmoid routing) and a smooth veto probability, trained by SGD in JAX; then rounding and the exact evaluation.
8. **Closing the loop**: generate, optimise, re-weight with the new map, iterate; the trade-off curve of photon efficiency against background acceptance for each channel budget.

## Differentiable trigger and veto (towards AutoDetOpt)

9. **Relaxed channel map**: soft assignments G (patterns x N_G, row softmax) and C (tiles x N_C, row softmax), sigmoid routing S, mask as a probabilistic OR, V_it = 1 - prod_{k,l} (1 - G_ik S_kl C_tl); event pass probability 1 - prod_sides [1 - prod_{t fired} (1 - V_it V_jt)] on the Geant4 photons; penalty on M (1 - V) for dropped RoI tiles. Adam in JAX, annealed temperatures or Gumbel-softmax, rounding and the exact event evaluation at the end. N_G, N_C stay an outer, discrete budget loop.
10. **Soft thresholds**: sigmoid of the energy deposit for the PSD (0.333 MeV) and FIT (1e-4 GeV) activations, so that the thresholds get gradients.
11. **Differentiable RoIs**: the RoI test (tile corners inside the cone through the first and last mats) as a sigmoid of the angular margin, giving gradients with respect to tile size, PSD-FIT distance, mat width and plane pitch; the L1 efficiency multiplies the effective area and the residual charged rate enters the background of AutoDetOpt's figures of merit. Geant4 events enter by reweighting, as in AutoDetOpt's stochastic part.
