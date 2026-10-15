"""The 0%-loss case: channel maps that reproduce every pattern's veto region exactly, and how few channels they can use.

The incidence matrix M (patterns x tiles, True when the tile is in the pattern's RoI) is to be realised by two merges and a routing:

    pattern groups   g(p) in 1..K          (trigger signals merged: one veto decision per group)
    tile channels    c(t) in 1..L          (veto signals merged: tiles OR-ed into one channel)
    routing          S, K x L boolean      (channel l vetoes group k)

and the veto mask of pattern p is V[p, t] = S[g(p), c(t)]. Lossless means V = M: no tile outside a pattern's RoI can veto it (no extra gamma
loss) and every tile inside does (no extra background leak). With the merges as partitions, as in hardware where a signal goes to one channel:

- two patterns can share a group only if their RoIs are identical, so K >= the number of distinct rows of M;
- two tiles can share a channel only if they are in exactly the same RoIs, so L >= the number of distinct columns of M;
- grouping the identical rows and the identical columns reaches both bounds at once. The minimum is exact and closed-form: K* and L*.

If a tile may feed several OR channels (fan-out in the trigger logic) the channels need not partition the tiles. The fewest OR channels
that reproduce M is then its Boolean rank, the size of a minimum biclique cover of the pattern-tile graph (NP-hard in general). It is at most
min(K*, L*); computeBooleanRankBounds gives a fooling-set lower bound and a greedy cover as the upper bound.

Naming: camelCase everywhere.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class ChannelMap:
    rowGroup: np.ndarray   # (patterns,) group index of each pattern, 0..K-1
    colGroup: np.ndarray   # (tiles,) channel index of each tile, 0..L-1
    routing: np.ndarray    # (K, L) bool

    @property
    def numberOfGroups(self):
        return self.routing.shape[0]

    @property
    def numberOfChannels(self):
        return self.routing.shape[1]

    def vetoMask(self):
        """V[p, t]: does tile t veto pattern p under this map."""
        return self.routing[np.ix_(self.rowGroup, self.colGroup)]


def findLosslessChannelMap(matrix):
    """The map with the fewest groups and channels that reproduces the matrix exactly: identical rows grouped, identical columns grouped."""
    uniqueRows, rowGroup = np.unique(matrix, axis=0, return_inverse=True)
    uniqueBoth, colGroup = np.unique(uniqueRows, axis=1, return_inverse=True)
    return ChannelMap(rowGroup.ravel(), colGroup.ravel(), uniqueBoth.astype(bool))


def computeUnionRouting(matrix, rowGroup, colGroup, numberOfGroups=None, numberOfChannels=None):
    """The 'safe' routing for given merges: channel l vetoes group k if any tile of l is in the RoI of any pattern of k.
    The mask then contains every RoI (V >= M): merging can only add tiles, so it costs gamma efficiency but never lets background through."""
    numberOfGroups = rowGroup.max() + 1 if numberOfGroups is None else numberOfGroups
    numberOfChannels = colGroup.max() + 1 if numberOfChannels is None else numberOfChannels
    rowOneHot = np.zeros((matrix.shape[0], numberOfGroups)); rowOneHot[np.arange(matrix.shape[0]), rowGroup] = 1.0
    colOneHot = np.zeros((matrix.shape[1], numberOfChannels)); colOneHot[np.arange(matrix.shape[1]), colGroup] = 1.0
    return ChannelMap(rowGroup, colGroup, (rowOneHot.T @ matrix.astype(float) @ colOneHot) > 0)


def computeFoolingSet(matrix, seed=0, restarts=4):
    """A fooling set of 1-entries (greedy, best of a few random orders): no two of them fit in one all-ones rectangle, so each needs its own
    biclique and the size is a lower bound on the Boolean rank. Returns the list of (row, column)."""
    rows, cols = np.nonzero(matrix)
    generator = np.random.default_rng(seed)
    best = []
    for restart in range(restarts):
        order = np.argsort(matrix.sum(axis=1)[rows] * matrix.sum(axis=0)[cols] + generator.random(rows.size)) if restart == 0 else generator.permutation(rows.size)
        chosenRows, chosenCols = [], []
        for index in order:
            i, j = rows[index], cols[index]
            # (i, j) and (i', j') fool each other unless M[i, j'] and M[i', j] are both 1; keep (i, j) only if it fools every chosen one.
            if chosenRows and np.any(matrix[i, chosenCols] & matrix[chosenRows, j]):
                continue
            chosenRows.append(i); chosenCols.append(j)
        if len(chosenRows) > len(best):
            best = list(zip(chosenRows, chosenCols))
    return best


def computeGreedyBicliqueCover(matrix):
    """Bicliques (row set, column set) whose union is exactly the 1-entries of the matrix, chosen greedily: each one starts from the row whose
    RoI has the most uncovered tiles and adds every row whose RoI contains it. Its length is an upper bound on the Boolean rank."""
    uncovered = matrix.copy()
    bicliques = []
    while uncovered.any():
        seed = np.argmax(uncovered.sum(axis=1))
        columns = matrix[seed].copy()
        bestRows, bestColumns, bestGain = None, None, -1
        # Try the seed's full RoI, then shrink it towards the uncovered part, keeping the biclique that covers the most new entries.
        for candidateColumns in (columns, columns & uncovered[seed]):
            rows = np.all(matrix[:, candidateColumns], axis=1)
            gain = uncovered[np.ix_(rows, candidateColumns)].sum()
            if gain > bestGain:
                bestRows, bestColumns, bestGain = rows, candidateColumns, gain
        uncovered[np.ix_(bestRows, bestColumns)] = False
        bicliques.append((np.flatnonzero(bestRows), np.flatnonzero(bestColumns)))
    return bicliques


def computeBooleanRankBounds(matrix):
    """(lower, upper) bounds on the Boolean rank: the fooling set and the smaller of the greedy cover and min(distinct rows, distinct columns)."""
    deduplicated = np.unique(np.unique(matrix, axis=0), axis=1)
    lower = len(computeFoolingSet(deduplicated))
    upper = min(len(computeGreedyBicliqueCover(deduplicated)), *deduplicated.shape)
    return lower, upper
