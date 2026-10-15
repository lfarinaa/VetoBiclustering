# stochastic: event-weighted channel maps, optimised on simulated events

**Not started.** The analytic part (`../analytic`) works on the RoI map alone. Below the lossless channel counts what a merge costs depends on the events: how often each pattern
fires for a photon, how often each tile fires from backsplash or noise in those photons, and how often background crosses a tile of the RoI that a merge leaves out. This part is to
estimate those rates from Geant4 and to optimise the channel map on them.

## The objective

For a channel map (pattern groups g, tile channels c, routing S) and an event sample,

    loss = photons vetoed / photons triggered  +  lambda x background not vetoed / background triggered and crossing a tile

each term computed with the real veto of the event (the first triggered planes, the AND of the two views, the energy threshold: `../legacy/aeffcalc/algos/PsdRoiVeto.cpp`).
A first-order version is linear in the mask: C[p, t] = photon events with pattern p and tile t fired, which `../analytic/mergeUtils.py` already takes as its `cost`.

## Optimisers, roughly in the order to try them

1. **Greedy agglomeration** with the event cost (exists, `mergeGreedily(cost=...)`): the baseline.
2. **Local search on the labels**: move one pattern to another group or one tile to another channel, flip one routing bit; simulated annealing or basin hopping (scipy's `basinhopping` with a custom step on the integer labels) on the exact event loss.
3. **Relaxation and SGD**: soft group and channel assignments (softmax), soft routing (sigmoid), a smooth veto probability 1 - prod(1 - mask x hit); train with Adam (optax) on batches of events, anneal the temperature, round, evaluate exactly. Gumbel-softmax if the rounding gap is large.
4. **Iterate with generation**: optimise on one sample, regenerate or re-weight (the map changes which events are triggered and vetoed, not the physics), check on an independent sample.

## Data

Geant4 events of the compact-v2 layout already exist from the AeffCalc studies (`~/AeffCalc/oneshot_*_compact/`: photons, protons, electrons, positrons, alphas, with the FIT mat hits
and PSD hits digitised). A reader for those ROOT files (uproot) is the first step.
