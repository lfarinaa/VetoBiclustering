# The channel-map problem in detail

This page states the problem that VetoBiclustering solves, its mathematical class, how hard it is, the variants it leaves aside, and the notation of
the High Energy cosmic-Radiation Detection facility (HERD) showcase. The [README](../README.md) has the short version; results are in
[HERD_SHOWCASE.md](HERD_SHOWCASE.md).

## Statement

A region-of-interest (RoI) veto vetoes each trigger pattern with the veto tiles that a charged particle producing it must have crossed. The RoIs form
the incidence matrix $`M\in\{0,1\}^{P\times T}`$, with one row per trigger pattern and one column per tile. A channel map merges the patterns into
$`N_G`$ groups and the tiles into $`N_C`$ channels, and its cost is written as

```math
\mathcal{C}(g, c, S) = \sum_{i=1}^{P} \sum_{t=1}^{T} \Big[\, C^+_{it}\,(1 - M_{it})\,V_{it} \;+\; C^-_{it}\,M_{it}\,(1 - V_{it}) \Big],
\qquad V_{it} = S_{g(i)\,c(t)},
```

where
- $`M_{it}`$ is 1 when tile $`t`$ is in the RoI of pattern $`i`$ (written $`t\in\{T\}_i`$ in the thesis) and 0 otherwise;
- $`g(i)\in\{1,\dots,N_G\}`$ is the group of pattern $`i`$, and $`c(t)\in\{1,\dots,N_C\}`$ the channel of tile $`t`$;
- $`S_{kl}`$ is the routing, 1 when channel $`l`$ vetoes group $`k`$ and 0 otherwise;
- $`V_{it}`$ is the veto mask of the merged map, 1 when tile $`t`$ vetoes pattern $`i`$;
- $`C^+_{it}\ge 0`$ is the cost of a false veto (tile $`t`$ vetoes pattern $`i`$ outside its RoI, so photons are lost), dimensionless;
- $`C^-_{it}\ge 0`$ is the cost of a leak (an RoI tile no longer vetoes its pattern, so charged particles get through), dimensionless.

The search minimises $`\mathcal{C}`$ over $`g`$, $`c`$ and $`S`$ for given budgets $`N_G`$ and $`N_C`$. The groups and channels cut $`M`$ into an
$`N_G\times N_C`$ checkerboard of blocks on which the mask is constant, and for fixed partitions the routing decouples block by block. A block is
switched on when the leaks that leaving it off would let through cost at least as much as the false vetoes that switching it on adds,

```math
S_{kl} = 1 \iff \sum_{i\in g^{-1}(k)}\;\sum_{t\in c^{-1}(l)} C^-_{it}\,M_{it} \;\ge\; \sum_{i\in g^{-1}(k)}\;\sum_{t\in c^{-1}(l)} C^+_{it}\,(1 - M_{it}),
```

where $`g^{-1}(k)`$ is the set of patterns in group $`k`$ and $`c^{-1}(l)`$ the set of tiles in channel $`l`$. The search is therefore over the two
partitions only.

The version solved so far is the **no-leak limit** $`C^-\to\infty`$. A block is then switched on as soon as it holds one RoI tile, which is called
union routing,

```math
S_{kl} = \max_{i\in g^{-1}(k),\; t\in c^{-1}(l)} M_{it},
```

so that $`V\ge M`$ entry by entry and only the false-veto cost remains. Pricing leaks needs charged-particle events, and the HERD showcase has none left.

## Problem class

The problem is a **constrained co-clustering (biclustering)** of a binary matrix, which the table spells out.

| Aspect | This problem |
|---|---|
| Matrix | the binary pattern x tile incidence matrix $`M`$ (1780 x 1328 in the HERD showcase), with two non-negative weight matrices $`C^\pm`$ of the same shape |
| Variables | two hard partitions, of the rows (patterns into $`N_G`$ groups) and of the columns (tiles into $`N_C`$ channels), and the block values $`S`$ |
| Model | the Boolean matrix tri-factorisation $`V = G\,S\,H^\top`$, a matrix product in which sums are replaced by logical ORs; $`G\in\{0,1\}^{P\times N_G}`$ and $`H\in\{0,1\}^{T\times N_C}`$ are cluster-indicator matrices with one 1 per row, so $`V`$ is a block-constant binary approximation of $`M`$ |
| Constraints | at most $`N_G`$ row clusters and $`N_C`$ column clusters; optionally, labels that no cluster may mix (cannot-link constraints); in the no-leak limit, $`V\ge M`$ |
| Objective | a weighted, asymmetric Hamming distance between $`V`$ and $`M`$, in which false positives (lost signal) and false negatives (background) are priced separately |
| Exact solution | only at zero cost: $`V=M`$ is reachable when $`N_G\ge N_G^*`$, the number of distinct rows of $`M`$, and $`N_C\ge N_C^*`$, the number of distinct columns |

With one side fixed the problem reduces to an ordinary clustering. For fixed channels, grouping the patterns is a $`k`$-clustering of the rows of the
channel-level matrix into $`k=N_G`$ clusters, each represented by a binary centre (the OR of its members in the no-leak limit, a weighted majority in
general), and each pattern costs its weighted distance to that centre.

## How hard it is

Below the zero-cost point there is no exact solution, and the maps found are not guaranteed to be optimal, for four reasons.

- The general problem is NP-hard. With equal costs ($`C^+=C^-=1`$) it is the monochromatic bi-clustering of Wulff, Urner and Ben-David (2013), which
  they prove NP-hard; the weighted problem contains it. They also give an approximation scheme for the equal-cost case, which has not been checked
  for these weights.
- The no-leak limit, the version solved so far, has no proof of hardness here, and no polynomial algorithm is known to us; it is treated as hard.
- Exhaustive search is impossible. In the HERD showcase at $`(N_G, N_C) = (5, 8)`$ there are about $`10^{1242}`$ ways to group the 1780 patterns and
  $`10^{1195}`$ ways to split the 1328 tiles; with one group per tracker face and every channel inside one PSD face, the tile splits alone still exceed
  $`10^{200}`$.
- No lower bound exists below the zero-cost point, so the gap to the optimum is unknown. Different merge orders and seeds land on the same cost for
  every constrained HERD budget, and at 5 groups both searches recover the tracker faces on their own, but neither fact is a proof. The second
  algorithm (the Jaccard merge, [below](#the-jaccard-merge)) beats the first by about 1% at some budgets, which shows directly that the first is not
  optimal there.

The combinatorial problem is also only a proxy. The real figure of merit is the event-level efficiency, through the two-view AND described below, and
the real costs depend on the energy spectrum, the coincidence window and the background rates, all of which carry uncertainty.

## From analytic seeds to a stochastic optimisation

The analytic maps are the starting point, not the end product. The next stage uses them to seed a stochastic optimisation on simulated events
(simulated annealing or basin hopping over the group and channel labels, and a relaxed, differentiable map trained by stochastic gradient descent),
with a signal and noise model adjusted to the measured rates: photon losses from backsplash and accidental hits at the measured tile rates, and the
leak of charged background at the real fluxes. That stage either confirms the analytic maps (no move improves them) or adjusts them to the more
realistic objective, and sweeping the leak cost against the false-veto cost traces photon efficiency against background rejection for each budget
([ROADMAP](../ROADMAP.md)).

## Variants

Three variants lead to different problems and are only mentioned here.

- **Fan-out wiring**, the layout of the thesis, lets a tile feed several groups (each face's veto ORs its own tiles and the adjacent rows of its
  neighbours). Each group then has its own OR channel over any tile set (the union of its RoIs in the no-leak limit), so only the patterns are
  clustered: this is the special case $`N_C=T`$ above, with the budget moved into the wiring (the tile-to-group links). Covering $`M`$ exactly with the
  fewest all-ones blocks is a set cover, the minimum biclique cover or Boolean rank of $`M`$, and deciding whether the Boolean rank is at most a given number is NP-complete (Orlin 1977).
- **Hardware constraints** (a channel confined to one PSD face, a group to one tracker face) give the same problem with cannot-link constraints
  between rows or columns with different labels.
- **The event-level objective** replaces the per-entry costs $`C^\pm`$, which are its proxy, by the efficiency measured on simulated events.

## Notation of the HERD showcase

The notation follows the thesis. The ultra-low-energy gamma-ray (ULEG) trigger fires at its first level (L0) on 3IR patterns, three hits in a row in
consecutive same-view planes of the scintillating-fibre tracker (FIT), and vetoes at its second level (L1) when a tile of the plastic scintillator
detector (PSD) fires inside the RoI. A pattern in an even view is labelled $`i`$, one in the odd view of the same double layer $`j`$, and the veto
region of an event is the event RoI (thesis Eq. 4.2),

```math
\{T\}_{\rm event} = \Big(\bigcup_i \{T\}^{\rm even}_i\Big) \cap \Big(\bigcup_j \{T\}^{\rm odd}_j\Big),
```

where $`\{T\}^{\rm even}_i`$ and $`\{T\}^{\rm odd}_j`$ are the RoIs of the triggered even-view and odd-view patterns. Merging into channels preserves
this intersection exactly, because $`V_{it}`$ depends on the tile only through its channel: an AND per channel gives the same veto as the AND per
tile, $`\bigvee_t h_t V_{it} V_{jt} = \bigvee_l h_l S_{g(i)l} S_{g(j)l}`$, where $`\bigvee`$ is the logical OR over the index, $`h_t`$ is 1 when
tile $`t`$ fires, and $`h_l`$ is 1 when any tile of channel $`l`$ fires.

The false-veto cost of adding tile $`t`$ to the veto of pattern $`i`$ is modelled as

```math
C^+_{it} = -\,w_i\,\ln\big[(1 - b_{it})(1 - q_t)\big], \qquad q_t = 1 - e^{-r_t A \tau},
```

where
- $`w_i`$ is the relative photon rate of pattern $`i`$, the fraction of triggering photons that fire it (dimensionless);
- $`b_{it}`$ is the backsplash probability, the probability that tile $`t`$ fires on the shower of a photon that fired pattern $`i`$ (dimensionless,
  steeply rising with energy);
- $`q_t`$ is the probability of an accidental hit on tile $`t`$ within the coincidence window (dimensionless);
- $`r_t`$ is the activation rate of tile $`t`$ from the background, in s⁻¹ cm⁻² (about 1.9 s⁻¹ cm⁻², so about 190 Hz per tile; thesis Fig. 4.9);
- $`A`$ is the tile area, 100 cm² for the 10 x 10 cm tiles of the thesis design;
- $`\tau`$ is the coincidence window of the L1 veto, in s (100 ns, a placeholder).

So far the search uses $`w_i=1`$ and $`b_{it}\approx b_t`$, the backsplash averaged over all triggering photons, because 340 photons cannot fill a
1780 x 1328 map. The efficiencies do not share this approximation, since they replay the photons event by event, and the stochastic stage handles the
pattern dependence automatically. The figure of merit is the L1 efficiency $`\varepsilon_{\rm L1}`$, the fraction of the L0-triggered fiducial photons
of a Geant4 simulation that also pass L1 with the merged map, accidental hits included. Two diagnostics of a mask complement it: the excess $`\xi`$,
the tiles added to the RoIs as a fraction of the RoI tiles, and the leak $`\kappa`$, the RoI tiles missing from the mask as the same fraction.

## The Jaccard merge

The second algorithm is an agglomerative hierarchical clustering: it starts from the lossless map, merges the two most similar groups (or channels),
recomputes the similarities of the merged group with the others, and repeats. The reasoning is geometric. Patterns produced by nearby tracks cross
nearly the same tiles, so their RoIs overlap almost entirely, and a merged group is vetoed by the union of its members' RoIs. Two groups with
nearly identical masks therefore lose little by sharing one, and merging the most similar pair first postpones the costly merges to the end.
Repeating the step builds a tree of nested groupings, the dendrogram, in which each budget is a cut at a given number of groups.

Similarity is the weighted Jaccard index of the two veto masks, defined as

```math
J(a, b) = \frac{\sum_x w_x\,[x\in A]\,[x\in B]}{\sum_x w_x\,[x\in A\cup B]},
```

where
- $`a`$ and $`b`$ are two groups (or two channels), and $`A`$ and $`B`$ their masks: for a group, the channels its patterns are vetoed by; for a
  channel, the groups it vetoes;
- $`x`$ runs over the items of the other side (channels when groups are compared, groups when channels are compared);
- $`w_x`$ is the weight of item $`x`$: its number of tiles or patterns in the plain variant, its false-veto cost in the cost-weighted one;
- $`[\,\cdot\,]`$ is the indicator, 1 when the condition holds and 0 otherwise.

$`J`$ is 1 for identical masks and 0 for disjoint ones, and the pair with the smallest distance $`1-J`$ merges first. The distance is close to the
merge cost itself: with unit costs, merging two single patterns adds to their masks the tiles in one RoI but not the other, $`|A\cup B|\,(1-J)`$
tiles. It differs from the cost in two ways, which is why the Jaccard merge is a separate algorithm rather than an approximation of the greedy: it
divides by the size of the union, and it ignores how many patterns each group holds, although every pattern in a group pays for the tiles added to
the group's mask. The best variant multiplies the
cost-weighted distance by Ward's size penalty $`n_a n_b/(n_a + n_b)`$ (Ward 1963), where $`n_a`$ and $`n_b`$ are the cost weights of $`a`$ and $`b`$,
so that two large clusters merge last; this restores part of the dependence on group size. Merging all groups before any channel works best. A merge path is shared by every budget below where it stops,
so both sides can be merged once to a common trunk, such as $`(N_G, N_C) = (100, 100)`$, and every smaller budget branched from there.
