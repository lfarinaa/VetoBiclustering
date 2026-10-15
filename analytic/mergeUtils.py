"""Channel maps below the lossless counts: the merges of 3IR patterns into N_G groups and of PSD tiles into N_C channels that cost the least.

Notation (examples/channelBudgetShowcase.ipynb, section 1): M[i, t] = 1 when tile t is in the RoI {T}_i of pattern i; a channel map is the group g(i)
of every pattern, the channel c(t) of every tile and the routing S (N_G x N_C), and its veto mask is V[i, t] = S[g(i), c(t)].

The cost of a map is sum_{i,t} (V - M)[i, t] C[i, t], with C[i, t] the expected false vetoes of good photons when tile t joins the mask of pattern i
(pattern rate x backsplash probability of the tile, from Geant4 events eventually). By default C = 1: the cost counts the (pattern, tile) pairs added
to masks, the mask excess. With union routing (losslessUtils.computeUnionRouting) V contains M, so no RoI tile is ever dropped (no extra leak)
and, writing Cg[k, l] for the cost summed over group k and channel l, the cost is sum_{k,l} S[k, l] Cg[k, l] - sum M C.

Two searches, both heuristics (they reach local optima, not certified ones):

- mergeGreedily: start from the lossless map and merge, one pair at a time, the two channels (then the two groups) whose merge adds the least cost.
  The merge costs are exact and kept up to date incrementally (O(N_G N_C) per merge):
      channels l1, l2:  sum_k Cg[k, l2] S[k, l1] (1 - S[k, l2]) + (l1 <-> l2)
      groups k1, k2:    sum_l Cg[k2, l] S[k1, l] (1 - S[k2, l]) + (k1 <-> k2)
- refineLocally: from any map, move single tiles between channels and single patterns between groups while the exact cost drops (coordinate
  descent on the labels, union routing throughout), until no single move helps.

Both accept labels that a merge may not cross (for example the PSD face of each tile, the FIT face of each pattern), for hardware that
segments the veto per face. Naming: camelCase everywhere.
"""
import numpy as np

import losslessUtils as lu


def _oneHot(labels, size):
    matrix = np.zeros((labels.size, size))
    matrix[np.arange(labels.size), labels] = 1.0
    return matrix


def _losslessWithLabels(matrix, rowLabels, colLabels):
    """The lossless map, keeping apart rows (columns) with different labels so that later merges can respect them."""
    rowKey = matrix if rowLabels is None else np.hstack([matrix, _oneHot(np.unique(rowLabels, return_inverse=True)[1].ravel(), np.unique(rowLabels).size).astype(bool)])
    _, rowGroup = np.unique(rowKey, axis=0, return_inverse=True)
    rowGroup = rowGroup.ravel()
    groupRows = np.zeros((rowGroup.max() + 1, matrix.shape[1]), dtype=bool)
    groupRows[rowGroup] = matrix
    colKey = groupRows.T if colLabels is None else np.hstack([groupRows.T, _oneHot(np.unique(colLabels, return_inverse=True)[1].ravel(), np.unique(colLabels).size).astype(bool)])
    _, colGroup = np.unique(colKey, axis=0, return_inverse=True)
    return lu.computeUnionRouting(matrix, rowGroup, colGroup.ravel())


def computeMaskCost(matrix, channelMap, cost=None):
    """(added cost, leaked RoI entries) of a map: the cost of the mask entries outside the RoIs, and the RoI entries missing from the masks."""
    mask = channelMap.vetoMask()
    cost = np.ones(matrix.shape) if cost is None else cost
    return float((cost * (mask & ~matrix)).sum()), int((matrix & ~mask).sum())


class _MinCache:
    """A symmetric merge-cost matrix (inf where a merge is forbidden) with the minimum of every row cached, so that the cheapest merge is found in
    O(N) and a merge of keep and drop costs O(N) plus a full pass over the few rows whose minimum sat in those two columns."""

    def __init__(self, table, blocked):
        self.cost = table + table.T
        self.cost[blocked] = np.inf
        self.rowMin, self.rowArg = self.cost.min(axis=1), self.cost.argmin(axis=1)

    def best(self):
        a = int(np.argmin(self.rowMin))
        return self.rowMin[a], a, int(self.rowArg[a])

    def merge(self, keep, drop, keepRow):
        """keepRow: the new merge costs of keep with every other entry (inf where forbidden)."""
        self.cost[drop, :] = self.cost[:, drop] = np.inf
        self.cost[keep, :] = self.cost[:, keep] = keepRow
        self.rowMin[drop] = np.inf
        stale = np.flatnonzero((self.rowArg == keep) | (self.rowArg == drop))
        improved = keepRow < self.rowMin
        self.rowMin[improved], self.rowArg[improved] = keepRow[improved], keep
        stale = np.union1d(stale, [keep])
        stale = stale[stale != drop]
        self.rowMin[stale], self.rowArg[stale] = self.cost[stale].min(axis=1), self.cost[stale].argmin(axis=1)
        self.rowMin[drop] = np.inf


class _GreedyState:
    """Live groups and channels of an agglomeration, with both pairwise merge-cost tables kept current through every merge.
    Table entry [a, b] is the cost added to a's members when b's mask joins them; the merge cost of a and b is table[a, b] + table[b, a]."""

    def __init__(self, matrix, cost, rowLabels, colLabels, start=None):
        begin = _losslessWithLabels(matrix, rowLabels, colLabels) if start is None else lu.computeUnionRouting(matrix, start.rowGroup, start.colGroup)
        self.rowGroup, self.colGroup = begin.rowGroup.copy(), begin.colGroup.copy()
        self.routing = begin.routing.copy()
        numberOfGroups, numberOfChannels = self.routing.shape
        self.groupCost = _oneHot(self.rowGroup, numberOfGroups).T @ cost @ _oneHot(self.colGroup, numberOfChannels)
        self.groupLabel = np.zeros(numberOfGroups, int) if rowLabels is None else np.array([rowLabels[np.flatnonzero(self.rowGroup == k)[0]] for k in range(numberOfGroups)])
        self.channelLabel = np.zeros(numberOfChannels, int) if colLabels is None else np.array([colLabels[np.flatnonzero(self.colGroup == l)[0]] for l in range(numberOfChannels)])
        self.aliveGroups, self.aliveChannels = np.ones(numberOfGroups, bool), np.ones(numberOfChannels, bool)
        self.total = 0.0
        self.channelTable, self.groupTable = None, None                  # Built on first use; kept current incrementally only when interleaving.
        self.channelCache, self.groupCache = None, None                  # Merge costs with cached row minima, rebuilt from a table when it is rebuilt.
        self.incremental = False
        self.channelBlocked = self._blocked(self.channelLabel, self.aliveChannels)
        self.groupBlocked = self._blocked(self.groupLabel, self.aliveGroups)

    @staticmethod
    def _blocked(labels, alive):
        blocked = (labels[:, None] != labels[None, :]) | ~alive[:, None] | ~alive[None, :]
        np.fill_diagonal(blocked, True)
        return blocked

    def _groupContribution(self, l):      # Channel l's share of the group table.
        s = self.routing[:, l].astype(float)
        return np.outer((1.0 - s) * self.groupCost[:, l], s)

    def _channelContribution(self, k):    # Group k's share of the channel table.
        s = self.routing[k].astype(float)
        return np.outer(s, (1.0 - s) * self.groupCost[k])

    def _tables(self, channel=True, group=True):
        s, c = self.routing.astype(float), self.groupCost
        if channel and self.channelTable is None:
            self.channelTable = s.T @ ((1.0 - s) * c)
        if group and self.groupTable is None:
            self.groupTable = ((1.0 - s) * c) @ s.T

    def bestChannelMerge(self):
        self._tables(channel=True, group=self.incremental)
        if self.channelCache is None:
            self.channelCache = _MinCache(self.channelTable, self.channelBlocked)
        return self.channelCache.best()

    def bestGroupMerge(self):
        self._tables(channel=self.incremental, group=True)
        if self.groupCache is None:
            self.groupCache = _MinCache(self.groupTable, self.groupBlocked)
        return self.groupCache.best()

    def applyChannelMerge(self, added, keep, drop):
        self.total += added
        if self.incremental:
            self.groupTable -= self._groupContribution(keep) + self._groupContribution(drop)
        self.routing[:, keep] |= self.routing[:, drop]; self.routing[:, drop] = False
        self.groupCost[:, keep] += self.groupCost[:, drop]; self.groupCost[:, drop] = 0.0
        if self.incremental:
            self.groupTable += self._groupContribution(keep)
        else:
            self.groupTable = None
        self.groupCache = None
        self.colGroup[self.colGroup == drop] = keep
        self.aliveChannels[drop] = False
        self.channelBlocked[drop, :] = self.channelBlocked[:, drop] = True
        if self.channelTable is not None:
            s, c = self.routing.astype(float), self.groupCost
            self.channelTable[keep, :] = s[:, keep] @ ((1.0 - s) * c)
            self.channelTable[:, keep] = s.T @ ((1.0 - s[:, keep]) * c[:, keep])
            if self.channelCache is not None:
                keepRow = self.channelTable[keep, :] + self.channelTable[:, keep]
                keepRow[self.channelBlocked[keep]] = np.inf
                self.channelCache.merge(keep, drop, keepRow)
        else:
            self.channelCache = None

    def applyGroupMerge(self, added, keep, drop):
        self.total += added
        if self.incremental:
            self.channelTable -= self._channelContribution(keep) + self._channelContribution(drop)
        self.routing[keep] |= self.routing[drop]; self.routing[drop] = False
        self.groupCost[keep] += self.groupCost[drop]; self.groupCost[drop] = 0.0
        if self.incremental:
            self.channelTable += self._channelContribution(keep)
        else:
            self.channelTable = None
        self.channelCache = None
        self.rowGroup[self.rowGroup == drop] = keep
        self.aliveGroups[drop] = False
        self.groupBlocked[drop, :] = self.groupBlocked[:, drop] = True
        if self.groupTable is not None:
            s, c = self.routing.astype(float), self.groupCost
            self.groupTable[keep, :] = ((1.0 - s[keep]) * c[keep]) @ s.T
            self.groupTable[:, keep] = ((1.0 - s) * c) @ s[keep]
            if self.groupCache is not None:
                keepRow = self.groupTable[keep, :] + self.groupTable[:, keep]
                keepRow[self.groupBlocked[keep]] = np.inf
                self.groupCache.merge(keep, drop, keepRow)
        else:
            self.groupCache = None

    def mergeTo(self, numberOfGroups, numberOfChannels, order="interleaved", groupHistory=None, channelHistory=None):
        """Merge down to the budget: 'channelsFirst', 'groupsFirst', or 'interleaved' (always the cheaper of the two best merges)."""
        self.incremental = order == "interleaved"
        while True:
            wantChannels = self.aliveChannels.sum() > numberOfChannels
            wantGroups = self.aliveGroups.sum() > numberOfGroups
            channelMove = self.bestChannelMerge() if wantChannels else (np.inf, 0, 0)
            groupMove = self.bestGroupMerge() if wantGroups else (np.inf, 0, 0)
            if order == "channelsFirst" and np.isfinite(channelMove[0]):
                groupMove = (np.inf, 0, 0)
            if order == "groupsFirst" and np.isfinite(groupMove[0]):
                channelMove = (np.inf, 0, 0)
            if not (np.isfinite(channelMove[0]) or np.isfinite(groupMove[0])):
                return                                                    # Budget reached, or the labels forbid any further merge.
            if channelMove[0] <= groupMove[0]:
                self.applyChannelMerge(*channelMove)
                if channelHistory is not None:
                    channelHistory[int(self.aliveChannels.sum())] = self.total
            else:
                self.applyGroupMerge(*groupMove)
                if groupHistory is not None:
                    groupHistory[int(self.aliveGroups.sum())] = self.total

    def mergeChannelsTo(self, numberOfChannels, history=None):
        self.mergeTo(self.aliveGroups.size, numberOfChannels, "channelsFirst", channelHistory=history)

    def mergeGroupsTo(self, numberOfGroups, history=None):
        self.mergeTo(numberOfGroups, self.aliveChannels.size, "groupsFirst", groupHistory=history)

    def copy(self):
        twin = object.__new__(_GreedyState)
        twin.__dict__ = {key: (value.copy() if isinstance(value, np.ndarray) else value) for key, value in self.__dict__.items()}
        twin.channelTable = twin.groupTable = None                       # Cheaper to rebuild in the branch than to copy and keep current.
        twin.channelCache = twin.groupCache = None
        return twin

    def channelMap(self):
        rowLive, colLive = np.flatnonzero(self.aliveGroups), np.flatnonzero(self.aliveChannels)
        rowIndex = np.full(self.aliveGroups.size, -1); rowIndex[rowLive] = np.arange(rowLive.size)
        colIndex = np.full(self.aliveChannels.size, -1); colIndex[colLive] = np.arange(colLive.size)
        return lu.ChannelMap(rowIndex[self.rowGroup], colIndex[self.colGroup], self.routing[np.ix_(rowLive, colLive)])


def mergeGreedily(matrix, numberOfGroups, numberOfChannels, cost=None, rowLabels=None, colLabels=None, record=False, order="channelsFirst", start=None):
    """A union-routed map with at most N_G groups and N_C channels (fewer if the labels forbid the merges), merging from the lossless map (or from a
    `start` map) in the given order. With record=True also returns the added cost after each group merge, {N_G: cost}."""
    cost = np.ones(matrix.shape) if cost is None else np.asarray(cost, dtype=float)
    state = _GreedyState(matrix, cost, rowLabels, colLabels, start)
    history = {}
    if order == "channelsFirst":
        state.mergeChannelsTo(numberOfChannels)
        history[int(state.aliveGroups.sum())] = state.total
    state.mergeTo(numberOfGroups, numberOfChannels, order, groupHistory=history)
    return (state.channelMap(), history) if record else state.channelMap()


def _splitByLabels(channelMap, rowLabels, colLabels):
    """The map whose groups (channels) are those of channelMap cut along the row (column) labels."""
    rowGroup, colGroup = channelMap.rowGroup, channelMap.colGroup
    if rowLabels is not None:
        rowGroup = np.unique(np.stack([rowGroup, np.unique(rowLabels, return_inverse=True)[1].ravel()], axis=1), axis=0, return_inverse=True)[1].ravel()
    if colLabels is not None:
        colGroup = np.unique(np.stack([colGroup, np.unique(colLabels, return_inverse=True)[1].ravel()], axis=1), axis=0, return_inverse=True)[1].ravel()
    return lu.ChannelMap(rowGroup, colGroup, np.zeros((rowGroup.max() + 1, colGroup.max() + 1), bool))


def findChannelMap(matrix, numberOfGroups, numberOfChannels, cost=None, rowLabels=None, colLabels=None, seeds=(), refine=True,
                   orders=("channelsFirst", "groupsFirst", "interleaved")):
    """The generic search: greedy merges in every order, plus every seed map (merged down to the budget if it is larger), each refined locally;
    returns (best map, [(name, added cost)] for every candidate). Seeds encode hardware-motivated layouts, such as one group per FIT face."""
    cost = np.ones(matrix.shape) if cost is None else np.asarray(cost, dtype=float)
    candidates = []
    starts = [(order, None, order) for order in orders] + [(f"seed {name}", seed, "interleaved") for name, seed in seeds]
    for name, seed, order in starts:
        if seed is not None:
            seed = _splitByLabels(seed, rowLabels, colLabels)            # A seed must not carry merges the labels forbid.
        channelMap = mergeGreedily(matrix, numberOfGroups, numberOfChannels, cost, rowLabels, colLabels, order=order, start=seed)
        if refine:
            channelMap = refineLocally(matrix, channelMap, cost, rowLabels, colLabels)
        candidates.append((name, computeMaskCost(matrix, channelMap, cost)[0], channelMap))
    best = min(candidates, key=lambda candidate: candidate[1])
    return best[2], [(name, added) for name, added, _ in candidates]


def scanBudgets(matrix, channelBudgets, groupBudgets, cost=None, rowLabels=None, colLabels=None, order="channelsFirst"):
    """Greedy maps on a grid of budgets, sharing the merges of the outer dimension ('channelsFirst': channels outer; 'groupsFirst': groups outer).
    Returns ({(N_G, N_C): ChannelMap}, {N: added cost} along the outer dimension alone). A budget the labels cannot reach gives the smallest map they allow."""
    cost = np.ones(matrix.shape) if cost is None else np.asarray(cost, dtype=float)
    state = _GreedyState(matrix, cost, rowLabels, colLabels)
    outerCurve, maps = {}, {}
    channelsOuter = order == "channelsFirst"
    for outer in sorted(channelBudgets if channelsOuter else groupBudgets, reverse=True):
        (state.mergeChannelsTo if channelsOuter else state.mergeGroupsTo)(outer, outerCurve)
        branch = state.copy()
        for inner in sorted(groupBudgets if channelsOuter else channelBudgets, reverse=True):
            (branch.mergeGroupsTo if channelsOuter else branch.mergeChannelsTo)(inner)
            maps[(inner, outer) if channelsOuter else (outer, inner)] = branch.channelMap()
    return maps, outerCurve


def _sweepMoves(assignment, features, weights, numberOfBins, itemLabels=None, binLabels=None):
    """One pass of single-item moves between bins, each accepted if it lowers sum_bins sum_p (count[p] > 0) sum[p].
    features (items x P, bool) and weights (items x P) describe each item against the other side's partition. Returns (assignment, improvement)."""
    counts = np.zeros((numberOfBins, features.shape[1])); sums = np.zeros((numberOfBins, features.shape[1]))
    np.add.at(counts, assignment, features); np.add.at(sums, assignment, weights)
    binCost = lambda c, s: ((c > 0) * s).sum(axis=-1)
    improvement = 0.0
    for item in np.random.default_rng(0).permutation(assignment.size):
        source = assignment[item]
        removedCost = binCost(counts[source] - features[item], sums[source] - weights[item]) - binCost(counts[source], sums[source])
        addedCost = binCost(counts + features[item], sums + weights[item]) - binCost(counts, sums)
        delta = removedCost + addedCost
        delta[source] = 0.0
        if binLabels is not None:
            delta[binLabels != itemLabels[item]] = np.inf
        target = int(np.argmin(delta))
        if delta[target] < -1e-12:
            counts[source] -= features[item]; sums[source] -= weights[item]
            counts[target] += features[item]; sums[target] += weights[item]
            assignment[item] = target
            improvement -= delta[target]
    return assignment, improvement


def refineLocally(matrix, channelMap, cost=None, rowLabels=None, colLabels=None, maxSweeps=20):
    """Coordinate descent from a map: alternately move tiles between channels and patterns between groups (union routing) while the cost drops.
    The budgets never grow (a bin may empty). Returns the refined map."""
    cost = np.ones(matrix.shape) if cost is None else np.asarray(cost, dtype=float)
    rowGroup, colGroup = channelMap.rowGroup.copy(), channelMap.colGroup.copy()
    numberOfGroups, numberOfChannels = channelMap.routing.shape
    groupLabels = None if rowLabels is None else np.array([rowLabels[np.flatnonzero(rowGroup == k)[0]] if (rowGroup == k).any() else None for k in range(numberOfGroups)], dtype=object)
    channelLabels = None if colLabels is None else np.array([colLabels[np.flatnonzero(colGroup == l)[0]] if (colGroup == l).any() else None for l in range(numberOfChannels)], dtype=object)
    for _ in range(maxSweeps):
        rowOneHot = _oneHot(rowGroup, numberOfGroups)
        tileFeatures = (rowOneHot.T @ matrix) > 0                           # [k, t]: is t in the RoI of any pattern of group k.
        colGroup, gainTiles = _sweepMoves(colGroup, tileFeatures.T.astype(float), (rowOneHot.T @ cost).T, numberOfChannels,
                                          None if colLabels is None else np.asarray(colLabels, dtype=object), channelLabels)
        colOneHot = _oneHot(colGroup, numberOfChannels)
        patternFeatures = (matrix @ colOneHot) > 0                          # [i, l]: does channel l hold a tile of pattern i's RoI.
        rowGroup, gainPatterns = _sweepMoves(rowGroup, patternFeatures.astype(float), cost @ colOneHot, numberOfGroups,
                                             None if rowLabels is None else np.asarray(rowLabels, dtype=object), groupLabels)
        if gainTiles + gainPatterns < 1e-9:
            break
    _, rowGroup = np.unique(rowGroup, return_inverse=True); _, colGroup = np.unique(colGroup, return_inverse=True)
    return lu.computeUnionRouting(matrix, rowGroup.ravel(), colGroup.ravel())


def computePairMetrics(matrix, mask, pairBlocks, backsplashProbability=None):
    """Figures of merit on the event RoI of the thesis, Eq. (4.2): for an even-view pattern i and an odd-view pattern j of the same double layer
    and FIT face the veto region is {T}_i ∩ {T}_j, and V_i ∩ V_j for the mask (with channel masks, the AND per channel equals the AND per tile).
    pairBlocks: [(evenRows, oddRows)] per FIT face (patternDump.evenOddBlocks). Averages are uniform over all pairs:
      excess    xi    = sum |V_ij \\ T_ij| / sum |T_ij|   (tiles added to the event veto regions)
      leak      kappa = sum |T_ij \\ V_ij| / sum |T_ij|   (event-RoI tiles that no longer veto)
      selfVeto  given the probability that a tile vetoes a good photon (a scalar p, or one value per tile, rateUtils.tileVetoProbability), the L1
                efficiency relative to the exact RoI veto, <prod_{t in V_ij} (1 - pi_t)> / <prod_{t in T_ij} (1 - pi_t)>."""
    m, v = matrix.astype(float), mask.astype(float)
    logSurvival = None
    if backsplashProbability is not None:
        logSurvival = np.log1p(-np.broadcast_to(np.asarray(backsplashProbability, dtype=float), (matrix.shape[1],)))
    exact, masked, kept, survivalMask, survivalExact = [], [], [], [], []
    for evenRows, oddRows in pairBlocks:
        exact.append((m[evenRows] @ m[oddRows].T).ravel())
        masked.append((v[evenRows] @ v[oddRows].T).ravel())
        kept.append(((m * v)[evenRows] @ (m * v)[oddRows].T).ravel())
        if logSurvival is not None:
            survivalMask.append(np.exp((v[evenRows] * logSurvival) @ v[oddRows].T).ravel())
            survivalExact.append(np.exp((m[evenRows] * logSurvival) @ m[oddRows].T).ravel())
    exact, masked, kept = map(np.concatenate, (exact, masked, kept))
    result = dict(excess=float((masked - kept).sum() / exact.sum()), leak=float((exact - kept).sum() / exact.sum()), pairs=int(exact.size),
                  meanExactTiles=float(exact.mean()), meanMaskTiles=float(masked.mean()))
    if logSurvival is not None:
        result["selfVeto"] = float(np.concatenate(survivalMask).mean() / np.concatenate(survivalExact).mean())
        result["absoluteSurvivalExact"] = float(np.concatenate(survivalExact).mean())
    return result


def channelMapFromMasks(matrix, groupOf, groupMasks):
    """A map that realises given per-group tile masks exactly (groupMasks: N_G x tiles): tiles with the same column in groupMasks share a channel."""
    _, colGroup = np.unique(groupMasks.T, axis=0, return_inverse=True)
    colGroup = colGroup.ravel()
    routing = (groupMasks.astype(float) @ _oneHot(colGroup, colGroup.max() + 1)) > 0
    return lu.ChannelMap(np.asarray(groupOf), colGroup, routing)


def describeChannelMap(channelMap, patterns, tileKeys, matrix=None, cost=None):
    """A plain-text listing of a map: every group with its patterns, every channel with its tiles (compressed per face and layer), the routing."""
    lines = []
    numberOfGroups, numberOfChannels = channelMap.routing.shape
    lines.append(f"N_G = {numberOfGroups} pattern groups, N_C = {numberOfChannels} PSD channels, {int(channelMap.routing.sum())} routing links")
    if matrix is not None:
        added, leaked = computeMaskCost(matrix, channelMap, cost)
        lines.append(f"mask excess {added / (matrix.sum() if cost is None else (np.asarray(cost) * matrix).sum()):.2%}, leaked RoI entries {leaked}")
    lines.append("\n# Pattern groups (3IR pattern FullIDs)")
    for k in range(numberOfGroups):
        names = [patterns[i].name for i in np.flatnonzero(channelMap.rowGroup == k)]
        lines.append(f"G{k}: {len(names)} patterns: {' '.join(names)}")
    lines.append("\n# PSD channels (face/layer: tile IDs)")
    for l in range(numberOfChannels):
        byFaceLayer = {}
        for column in np.flatnonzero(channelMap.colGroup == l):
            face, layer, tileId = tileKeys[column]
            byFaceLayer.setdefault((face, layer), []).append(tileId)
        text = "; ".join(f"{face}/L{layer}: {_ranges(sorted(ids))}" for (face, layer), ids in sorted(byFaceLayer.items()))
        lines.append(f"C{l}: {sum(map(len, byFaceLayer.values()))} tiles: {text}")
    lines.append("\n# Routing: channels that veto each group")
    for k in range(numberOfGroups):
        lines.append(f"G{k} <- {' '.join(f'C{l}' for l in np.flatnonzero(channelMap.routing[k]))}")
    return "\n".join(lines)


def _ranges(ids):
    parts, start = [], None
    for index, value in enumerate(ids):
        if start is None:
            start = value
        if index + 1 == len(ids) or ids[index + 1] != value + 1:
            parts.append(str(start) if start == value else f"{start}-{value}")
            start = None
    return ",".join(parts)
