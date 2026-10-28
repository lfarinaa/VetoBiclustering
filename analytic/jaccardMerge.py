"""An alternative route to channel maps: agglomeration by Jaccard similarity instead of by exact merge cost.

The exact-cost greedy search (mergeUtils.mergeGreedily) merges, at every step, the two groups (or channels) whose merge adds the least false-veto cost.
This route instead merges the two most similar ones, similarity being measured on the veto masks, and leaves the cost out of the choice. It is the
classic agglomerative clustering with a Jaccard similarity: cheaper to explain, blind to how many patterns each merge actually penalises.

Notation as in mergeUtils: M[i, t] is 1 when tile t is in the region of interest (RoI) of pattern i; a map has the groups g(i), the channels c(t) and
the union routing S[k, l] (1 when any tile of channel l is in the RoI of any pattern of group k), so that no RoI tile is ever dropped.

The Jaccard similarity of two sets A and B is |A intersect B| / |A union B|, which is 1 for identical sets and 0 for disjoint ones. Here it is computed with
weights, as

    J(a, b) = sum_x w_x [x in A] [x in B] / sum_x w_x [x in A or B],

where
- a and b are two groups (or two channels) and A, B their masks: for a group, the channels its patterns are vetoed by (row a of S); for a channel,
  the groups it vetoes (column a of S);
- x runs over the elements of the other side (channels when comparing groups, groups when comparing channels);
- w_x is the weight of element x, the number of tiles (patterns) it holds, so that J is the Jaccard similarity of the tile-level (pattern-level) sets;
- [condition] is the indicator, 1 when the condition holds and 0 otherwise.

Variants, chosen with `variant`:
- "jaccard": the weights above (tile and pattern counts).
- "costJaccard": each tile weighted by its mean false-veto cost over patterns, each pattern by its mean over tiles, so that overlap is measured in
  units of cost rather than of tiles.
- "wardJaccard": the costJaccard distance 1 - J multiplied by n_a n_b / (n_a + n_b), where n_a and n_b are the cost weights of a and b. This is the
  size penalty of Ward's linkage (Ward 1963): it postpones merges of two large clusters, which add the most cost under union routing.

Orders: "channelsFirst" (merge channels down to the budget, then groups), "groupsFirst", and "alternating" (always the side furthest above its budget,
relative to the budget). The exact cost of the current map is carried along (sum_{k,l} S[k, l] Z[k, l], where Z[k, l] is the cost summed over the
entries of block (k, l) outside the RoIs), so every run records how the true objective evolves while the similarity drives the merges.
"""
import numpy as np

import losslessUtils as lu
import mergeUtils as mu

variants = ("jaccard", "costJaccard", "wardJaccard")


class _JaccardState:
    """Live groups and channels of a Jaccard agglomeration. Each side keeps a table of pairwise distances with cached row minima (mergeUtils._MinCache):
    a merge on one side changes only the merged item's row of that side's table, and leaves the other side's table stale until it is next used."""

    def __init__(self, matrix, cost, variant, rowLabels=None, colLabels=None):
        if variant not in variants:
            raise ValueError(f"variant {variant!r}: one of {variants}")
        self.variant = variant
        begin = mu._losslessWithLabels(matrix, rowLabels, colLabels)
        self.rowGroup, self.colGroup = begin.rowGroup.copy(), begin.colGroup.copy()
        self.routing = begin.routing.copy()
        numberOfGroups, numberOfChannels = self.routing.shape
        rowOneHot, colOneHot = mu._oneHot(self.rowGroup, numberOfGroups), mu._oneHot(self.colGroup, numberOfChannels)
        self.zeroCost = rowOneHot.T @ (cost * ~matrix) @ colOneHot          # Z[k, l]: cost of the block's entries outside the RoIs.
        countWeights = variant == "jaccard"
        patternWeight = np.ones(matrix.shape[0]) if countWeights else cost.mean(axis=1)
        tileWeight = np.ones(matrix.shape[1]) if countWeights else cost.mean(axis=0)
        self.groupWeight, self.channelWeight = rowOneHot.T @ patternWeight, colOneHot.T @ tileWeight
        self.groupLabel = np.zeros(numberOfGroups, int) if rowLabels is None else np.array([rowLabels[np.flatnonzero(self.rowGroup == k)[0]] for k in range(numberOfGroups)])
        self.channelLabel = np.zeros(numberOfChannels, int) if colLabels is None else np.array([colLabels[np.flatnonzero(self.colGroup == l)[0]] for l in range(numberOfChannels)])
        self.aliveGroups, self.aliveChannels = np.ones(numberOfGroups, bool), np.ones(numberOfChannels, bool)
        self.groupCache = self.channelCache = None                          # Built on first use, dropped when the other side merges.
        self.groupSize = self.channelSize = None
        self.exactCost = float((self.routing * self.zeroCost).sum())

    def _sideArrays(self, side):
        """(masks of the side's items over the other side, the other side's weights, the side's own weights, labels, alive) for 'group' or 'channel'."""
        if side == "group":
            return self.routing, self.channelWeight, self.groupWeight, self.groupLabel, self.aliveGroups
        return self.routing.T, self.groupWeight, self.channelWeight, self.channelLabel, self.aliveChannels

    def _distanceRow(self, side, item, intersectionRow, size):
        """Distances of `item` to every item of its side (inf where forbidden), from its weighted intersections and the mask sizes."""
        _, _, ownWeight, labels, alive = self._sideArrays(side)
        union = size + size[item] - intersectionRow
        with np.errstate(invalid="ignore", divide="ignore"):
            distance = np.where(union > 0, 1.0 - intersectionRow / union, 0.0)
        if self.variant == "wardJaccard":
            distance = distance * ownWeight * ownWeight[item] / np.maximum(ownWeight + ownWeight[item], 1e-300)
        distance[~alive | (labels != labels[item])] = np.inf
        distance[item] = np.inf
        return distance

    def _cache(self, side):
        cache = self.groupCache if side == "group" else self.channelCache
        if cache is not None:
            return cache
        masks, otherWeight, ownWeight, labels, alive = self._sideArrays(side)
        live = np.flatnonzero(alive)
        liveMasks = masks[live].astype(float)
        intersection = (liveMasks * otherWeight) @ liveMasks.T
        size = np.zeros(alive.size); size[live] = np.diag(intersection)
        union = size[live][:, None] + size[live][None, :] - intersection
        with np.errstate(invalid="ignore", divide="ignore"):
            liveDistance = np.where(union > 0, 1.0 - intersection / union, 0.0)
        if self.variant == "wardJaccard":
            liveDistance = liveDistance * np.outer(ownWeight[live], ownWeight[live]) / np.maximum(ownWeight[live][:, None] + ownWeight[live][None, :], 1e-300)
        distance = np.full((alive.size, alive.size), np.inf)
        distance[np.ix_(live, live)] = liveDistance
        blocked = (labels[:, None] != labels[None, :]) | ~alive[:, None] | ~alive[None, :]
        np.fill_diagonal(blocked, True)
        cache = mu._MinCache(distance / 2, blocked)                        # _MinCache symmetrises table + table.T.
        if side == "group":
            self.groupCache, self.groupSize = cache, size
        else:
            self.channelCache, self.channelSize = cache, size
        return cache

    def bestGroupMerge(self):
        return self._cache("group").best()

    def bestChannelMerge(self):
        return self._cache("channel").best()

    def _applyMerge(self, side, keep, drop):
        routingView = self.routing if side == "group" else self.routing.T
        zeroView = self.zeroCost if side == "group" else self.zeroCost.T
        ownWeight = self.groupWeight if side == "group" else self.channelWeight
        otherWeight = self.channelWeight if side == "group" else self.groupWeight
        before = (routingView[keep] * zeroView[keep]).sum() + (routingView[drop] * zeroView[drop]).sum()
        routingView[keep] |= routingView[drop]; routingView[drop] = False
        zeroView[keep] += zeroView[drop]; zeroView[drop] = 0.0
        ownWeight[keep] += ownWeight[drop]; ownWeight[drop] = 0.0
        self.exactCost += float((routingView[keep] * zeroView[keep]).sum() - before)
        alive = self.aliveGroups if side == "group" else self.aliveChannels
        alive[drop] = False
        assignment = self.rowGroup if side == "group" else self.colGroup
        assignment[assignment == drop] = keep
        cache, size = (self.groupCache, self.groupSize) if side == "group" else (self.channelCache, self.channelSize)
        if cache is not None:
            mergedMask = routingView[keep].astype(float)
            intersectionRow = routingView.astype(float) @ (otherWeight * mergedMask)
            size[keep] = intersectionRow[keep]; size[drop] = 0.0
            cache.merge(keep, drop, self._distanceRow(side, keep, intersectionRow, size))
        if side == "group":
            self.channelCache = None                                        # Channel masks over the groups changed.
        else:
            self.groupCache = None

    def applyGroupMerge(self, keep, drop):
        self._applyMerge("group", keep, drop)

    def applyChannelMerge(self, keep, drop):
        self._applyMerge("channel", keep, drop)

    def mergeTo(self, numberOfGroups, numberOfChannels, order="channelsFirst", trace=None, inspect=None):
        """Merge down to the budget in the given order; append one record per merge to `trace` (a list) if given, and call
        inspect(state, side, keep, drop) before every merge if given (to compare the choice with alternatives). A side whose labels forbid every
        remaining merge stops above its budget."""
        stuck = {"group": False, "channel": False}
        while True:
            groupExcess = self.aliveGroups.sum() / numberOfGroups if not stuck["group"] else 0.0
            channelExcess = self.aliveChannels.sum() / numberOfChannels if not stuck["channel"] else 0.0
            wantGroups, wantChannels = groupExcess > 1, channelExcess > 1
            if not (wantGroups or wantChannels):
                return
            if order == "channelsFirst":
                side = "channel" if wantChannels else "group"
            elif order == "groupsFirst":
                side = "group" if wantGroups else "channel"
            else:
                side = "group" if wantGroups and (not wantChannels or groupExcess >= channelExcess) else "channel"
            distance, keep, drop = self.bestGroupMerge() if side == "group" else self.bestChannelMerge()
            keep, drop = min(keep, drop), max(keep, drop)
            if not np.isfinite(distance):
                stuck[side] = True
                continue
            if inspect is not None:
                inspect(self, side, keep, drop)
            costBefore = self.exactCost
            (self.applyGroupMerge if side == "group" else self.applyChannelMerge)(keep, drop)
            if trace is not None:
                trace.append(dict(side=side, keep=int(keep), drop=int(drop), distance=float(distance), groups=int(self.aliveGroups.sum()),
                                  channels=int(self.aliveChannels.sum()), exactCost=self.exactCost, addedCost=self.exactCost - costBefore))

    def copy(self):
        twin = object.__new__(_JaccardState)
        twin.__dict__ = {key: (value.copy() if isinstance(value, np.ndarray) else value) for key, value in self.__dict__.items()}
        twin.groupCache = twin.channelCache = None                         # Cheaper to rebuild in the branch than to copy.
        return twin

    def channelMap(self):
        rowLive, colLive = np.flatnonzero(self.aliveGroups), np.flatnonzero(self.aliveChannels)
        rowIndex = np.full(self.aliveGroups.size, -1); rowIndex[rowLive] = np.arange(rowLive.size)
        colIndex = np.full(self.aliveChannels.size, -1); colIndex[colLive] = np.arange(colLive.size)
        return lu.ChannelMap(rowIndex[self.rowGroup], colIndex[self.colGroup], self.routing[np.ix_(rowLive, colLive)])


def mergeByJaccard(matrix, numberOfGroups, numberOfChannels, cost=None, variant="jaccard", order="channelsFirst", rowLabels=None, colLabels=None):
    """A union-routed map with at most N_G groups and N_C channels, merging the most similar pair at every step. Returns (map, trace), the trace
    holding one record per merge: the side merged, the pair (keep, drop, as indices of the lossless groups or channels), the Jaccard distance of the
    pair, the live counts and the exact cost after the merge."""
    cost = np.ones(matrix.shape) if cost is None else np.asarray(cost, dtype=float)
    state = _JaccardState(matrix, cost, variant, rowLabels, colLabels)
    trace = []
    state.mergeTo(numberOfGroups, numberOfChannels, order, trace)
    return state.channelMap(), trace


def scanJaccardBudgets(matrix, groupBudgets, channelBudgets, cost=None, variant="jaccard", order="channelsFirst", rowLabels=None, colLabels=None,
                       trunk=None, trunkOrder="groupsFirst"):
    """Jaccard maps on a grid of budgets in one pass for the sequential orders: the outer side (channels for channelsFirst) is merged through its
    budgets once, and each stop is copied and merged through the inner budgets. Returns ({(N_G, N_C): map}, {(N_G, N_C): trace}), each trace
    covering the whole path from the lossless map to that budget.

    trunk, optional (N_G, N_C) at or above every budget: both sides are first merged down to it once (in trunkOrder), and every budget branches
    from there. A merge path is shared by every budget below where it stops, so the trunk is computed once; it also changes the path (some channels
    merge before the groups reach their budget), which on the HERD showcase moved the cost by -16% to +11% per budget (FINDINGS 10)."""
    if order not in ("channelsFirst", "groupsFirst"):
        raise ValueError("scanJaccardBudgets takes the sequential orders; use mergeByJaccard for 'alternating'")
    cost = np.ones(matrix.shape) if cost is None else np.asarray(cost, dtype=float)
    state = _JaccardState(matrix, cost, variant, rowLabels, colLabels)
    maps, traces, outerTrace = {}, {}, []
    if trunk is not None:
        if trunk[0] < max(groupBudgets) or trunk[1] < max(channelBudgets):
            raise ValueError(f"trunk {trunk} is below a budget")
        state.mergeTo(*trunk, trunkOrder, outerTrace)
    channelsOuter = order == "channelsFirst"
    everything = (int(state.aliveGroups.sum()), int(state.aliveChannels.sum()))
    for outer in sorted(channelBudgets if channelsOuter else groupBudgets, reverse=True):
        state.mergeTo(everything[0] if channelsOuter else outer, outer if channelsOuter else everything[1], order, outerTrace)
        branch, branchTrace = state.copy(), list(outerTrace)
        for inner in sorted(groupBudgets if channelsOuter else channelBudgets, reverse=True):
            budget = (inner, outer) if channelsOuter else (outer, inner)
            branch.mergeTo(*budget, order, branchTrace)
            maps[budget], traces[budget] = branch.channelMap(), list(branchTrace)
    return maps, traces


def linkageFromTrace(trace, side, numberOfItems):
    """The merges of one side of a trace as a SciPy linkage matrix, for scipy.cluster.hierarchy.dendrogram. Leaves are the lossless groups (or
    channels), 0..numberOfItems-1; row r joins two clusters into cluster numberOfItems + r, at a height equal to the exact cost after that merge
    (non-decreasing under union routing). Returns (linkage, members), members[cluster] being the leaves of every cluster. A path that stops above
    one cluster gives fewer than numberOfItems - 1 rows (a forest): dendrogram needs the full tree, so merge to one item when the tree is wanted."""
    clusterOf = list(range(numberOfItems))                                  # The cluster that currently holds each lossless item (as keep).
    members = [[item] for item in range(numberOfItems)]
    rows = []
    for record in trace:
        if record["side"] != side:
            continue
        first, second = clusterOf[record["keep"]], clusterOf[record["drop"]]
        members.append(members[first] + members[second])
        rows.append([first, second, record["exactCost"], len(members[-1])])
        clusterOf[record["keep"]] = len(members) - 1
    return np.array(rows, dtype=float), members
