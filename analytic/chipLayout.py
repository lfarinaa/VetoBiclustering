"""Channel maps with front-end chips: veto tiles are read by chips of limited input count, each giving one fast-OR output, and the chip outputs are
grouped into the boards that form the veto channels.

The readout has three levels:

    tiles t --(at most `capacity` neighbouring tiles of one face)--> chips --(any chips)--> boards (the N_C channels) --(routing S)--> pattern groups

A chip ORs the tiles wired to it, so a board, the OR of its chips, is a union of whole chips. The chips restrict which tile sets a board can hold:
without a limit on their number the best layout would be one tile per chip, so the number of chips is a budget (by default the minimum, the tiles of
each face divided by the capacity, rounded up). The cost of a layout is the false-veto cost of losslessUtils/mergeUtils under union routing, with
the channel of a tile being the board of its chip.

The search merges up to a point and then anneals:
1. pattern groups: the groups-first exact-cost greedy or the Ward-scaled Jaccard merge, with the tiles unmerged;
2. chips: merging neighbouring tiles of one face by exact cost (with the groups fixed) under the capacity, which stalls above the budget with chips
   too full to join, then a repair that dissolves the smallest chips through chains of neighbouring chips (repairChips); alternatives are a
   similarity-weighted METIS partition of each face (metisChips) and a serpentine tiling that ignores the cost;
3. boards: the exact-cost greedy merging whole chips;
4. simulated annealing over three moves: a tile to a neighbouring chip (capacity respected, the chip it leaves stays connected), a chip to another
   board, a pattern to another group. A move that lowers the cost is kept; one that raises it by Delta is kept with probability exp(-Delta / T),
   with the temperature T lowered geometrically.

Defaults are placeholders until the hardware is specified: 32 inputs per chip, one fast-OR output, neighbours in the 4-connected grid of a layer and
across the two layers (same or adjacent cell). Naming: camelCase everywhere.
"""
from dataclasses import dataclass
import heapq

import numpy as np
from numba import njit

import losslessUtils as lu
import mergeUtils as mu
import patternDump as pd
import geometryUtils as gu

defaultCapacity = 32                                                       # Inputs per front-end chip. Placeholder: settable, no datasheet yet.


# ----------------------------------------------------------------------------------------------------------------------------------------------
# Geometry: the neighbour graph and the chip budget
# ----------------------------------------------------------------------------------------------------------------------------------------------

def tileNeighbours(tileKeys, crossLayer=True):
    """The neighbour graph of the tiles, in compressed sparse row form (indexPointer, neighbourIndices): two tiles of the same face are neighbours
    when their grid cells share an edge in the same layer, or (crossLayer) when they sit in the same or edge-sharing cells of the two layers.
    Faces are never joined: a chip reads one face."""
    position = {key: index for index, key in enumerate(tileKeys)}
    neighbours = [[] for _ in tileKeys]
    steps = ((1, 0), (-1, 0), (0, 1), (0, -1))
    for index, (face, layer, tileId) in enumerate(tileKeys):
        row, column = gu.gridIndex(face, tileId)
        rows, columns = gu.gridShape[face]
        cells = [(row + dr, column + dc, layer) for dr, dc in steps]
        if crossLayer:
            otherLayer = 1 - layer
            cells += [(row, column, otherLayer)] + [(row + dr, column + dc, otherLayer) for dr, dc in steps]
        for cellRow, cellColumn, cellLayer in cells:
            if 0 <= cellRow < rows and 0 <= cellColumn < columns:
                neighbourId = cellRow * columns + cellColumn
                other = position.get((face, cellLayer, neighbourId))
                if other is not None:
                    neighbours[index].append(other)
    indexPointer = np.zeros(len(tileKeys) + 1, dtype=np.int64)
    indexPointer[1:] = np.cumsum([len(entry) for entry in neighbours])
    return indexPointer, np.array([other for entry in neighbours for other in entry], dtype=np.int64)


def minimumChips(tileKeys, capacity=defaultCapacity):
    """{face: the fewest chips that can read its tiles}, the tiles of the face divided by the capacity and rounded up."""
    faces = np.array([key[0] for key in tileKeys])
    return {face: int(np.ceil((faces == face).sum() / capacity)) for face in pd.psdSides if (faces == face).any()}


def serpentineChips(tileKeys, chipsPerFace):
    """A contiguous chip tiling: the cells of each face are walked row by row in alternating directions (both layers of a cell together), and the walk
    is cut into chipsPerFace[face] chunks of nearly equal size. Consecutive cells of the walk share an edge, so every chunk is connected under the
    cross-layer neighbour graph. Returns the chip of every tile."""
    chipOf = np.full(len(tileKeys), -1)
    nextChip = 0
    for face, numberOfChips in chipsPerFace.items():
        rows, columns = gu.gridShape[face]
        walk = []
        for row in range(rows):
            for column in (range(columns) if row % 2 == 0 else range(columns - 1, -1, -1)):
                walk.append(row * columns + column)
        tilesInWalk = [index for tileId in walk for index, key in enumerate(tileKeys) if key[0] == face and key[2] == tileId]
        for chunk in np.array_split(np.array(tilesInWalk), numberOfChips):
            chipOf[chunk] = nextChip
            nextChip += 1
    assert (chipOf >= 0).all(), "a tile outside the face grids"
    return chipOf


# ----------------------------------------------------------------------------------------------------------------------------------------------
# Merge stage: neighbouring tiles into chips, chips into boards
# ----------------------------------------------------------------------------------------------------------------------------------------------

def mergeTilesIntoChips(matrix, cost, tileKeys, rowGroup, chipsPerFace, capacity=defaultCapacity, neighbours=None):
    """Merge neighbouring clusters of one face, the pair that adds the least exact cost first (a chip treated as a channel, with the pattern groups
    `rowGroup` fixed), ties broken towards the smaller merged chip, under the capacity, until each face is down to chipsPerFace[face] chips or no
    allowed merge is left. Returns (chipOf, {face: chips reached}); a face above its budget has fragmented."""
    indexPointer, neighbourIndices = tileNeighbours(tileKeys) if neighbours is None else neighbours
    numberOfGroups = int(rowGroup.max()) + 1
    groupOneHot = mu._oneHot(rowGroup, numberOfGroups)
    touches = (groupOneHot.T @ matrix) > 0                                  # [k, t]: is tile t in the RoI of any pattern of group k.
    zeroCost = groupOneHot.T @ (cost * ~matrix)                             # [k, t]: false-veto cost of tile t for group k.
    faces = np.array([key[0] for key in tileKeys])
    numberOfTiles = len(tileKeys)
    clusterTouches = {t: touches[:, t].copy() for t in range(numberOfTiles)}
    clusterZero = {t: zeroCost[:, t].copy() for t in range(numberOfTiles)}
    clusterSize = {t: 1 for t in range(numberOfTiles)}
    clusterNeighbours = {t: set(neighbourIndices[indexPointer[t]:indexPointer[t + 1]].tolist()) for t in range(numberOfTiles)}
    clusterOf = np.arange(numberOfTiles)
    remaining = {face: int((faces == face).sum()) for face in chipsPerFace}

    def mergeCost(a, b):
        return float((clusterTouches[a] & ~clusterTouches[b]) @ clusterZero[b] + (clusterTouches[b] & ~clusterTouches[a]) @ clusterZero[a])

    heap = []
    for a in range(numberOfTiles):
        for b in clusterNeighbours[a]:
            if a < b:
                heapq.heappush(heap, (mergeCost(a, b), clusterSize[a] + clusterSize[b], a, b))
    while heap:
        added, size, a, b = heapq.heappop(heap)
        if a not in clusterSize or b not in clusterSize or clusterSize[a] + clusterSize[b] != size or b not in clusterNeighbours[a]:
            continue                                                        # Stale entry: a cluster merged or grew since it was pushed.
        face = faces[a]
        if size > capacity or remaining[face] <= chipsPerFace[face]:
            continue
        clusterTouches[a] |= clusterTouches.pop(b)
        clusterZero[a] += clusterZero.pop(b)
        clusterSize[a] += clusterSize.pop(b)
        for other in clusterNeighbours.pop(b):
            clusterNeighbours[other].discard(b)
            if other != a:
                clusterNeighbours[other].add(a)
                clusterNeighbours[a].add(other)
        clusterNeighbours[a].discard(a)
        clusterOf[clusterOf == b] = a
        remaining[face] -= 1
        for other in clusterNeighbours[a]:
            if clusterSize[a] + clusterSize[other] <= capacity:
                heapq.heappush(heap, (mergeCost(a, other), clusterSize[a] + clusterSize[other], min(a, other), max(a, other)))
    _, chipOf = np.unique(clusterOf, return_inverse=True)
    return chipOf.ravel(), dict(remaining)


def _chipProfiles(matrix, cost, rowGroup, chipOf):
    """Per chip, the groups whose RoIs it touches and its false-veto cost per group (the chip treated as a channel)."""
    numberOfGroups = int(rowGroup.max()) + 1
    groupOneHot = mu._oneHot(rowGroup, numberOfGroups)
    touches = (groupOneHot.T @ matrix).T                                    # [t, k]: patterns of group k with tile t in their RoI.
    zero = (groupOneHot.T @ (cost * ~matrix)).T                             # [t, k]: false-veto cost of tile t for group k.
    return touches, zero


def _isConnected(members, indexPointer, neighbourIndices):
    memberSet = set(int(member) for member in members)
    if not memberSet:
        return True
    first = next(iter(memberSet))
    seen, frontier = {first}, [first]
    while frontier:
        current = frontier.pop()
        for other in neighbourIndices[indexPointer[current]:indexPointer[current + 1]]:
            other = int(other)
            if other in memberSet and other not in seen:
                seen.add(other); frontier.append(other)
    return len(seen) == len(memberSet)


def repairChips(matrix, cost, tileKeys, rowGroup, chipOf, chipsPerFace, capacity=defaultCapacity, neighbours=None, maximumPathLength=6):
    """Bring a chip layout within the capacity and the budget chipsPerFace, keeping every chip connected. Oversized chips first shed boundary tiles,
    then the smallest chips of a face above its budget are dissolved. A tile goes to a neighbouring chip with room; when every neighbouring chip is
    full, room is made along a chain of neighbouring chips (each passing one boundary tile to the next, the last one having room), the shortest chain
    first, as an augmenting path in matching. Among the possible destinations the one whose groups the tile adds least false-veto cost to is
    chosen (the chip treated as a channel). A chip in several disconnected pieces is first split into one chip per piece. Returns (chipOf,
    {face: chips}); a face is left above its budget if no chain is found."""
    indexPointer, neighbourIndices = tileNeighbours(tileKeys) if neighbours is None else neighbours
    chipOf = np.unique(np.asarray(chipOf), return_inverse=True)[1].ravel().copy()
    faces = np.array([key[0] for key in tileKeys])
    nextChip = int(chipOf.max()) + 1
    for chip in np.unique(chipOf):                                          # A chip in several pieces becomes one chip per piece.
        remaining = set(np.flatnonzero(chipOf == chip).tolist())
        first = True
        while remaining:
            seed = remaining.pop()
            piece, frontier = {seed}, [seed]
            while frontier:
                current = frontier.pop()
                for other in neighbourIndices[indexPointer[current]:indexPointer[current + 1]]:
                    other = int(other)
                    if other in remaining:
                        remaining.discard(other); piece.add(other); frontier.append(other)
            if not first:
                chipOf[list(piece)] = nextChip
                nextChip += 1
            first = False
    tileTouches, tileZero = _chipProfiles(matrix, cost, rowGroup, chipOf)
    members = {chip: set(np.flatnonzero(chipOf == chip).tolist()) for chip in np.unique(chipOf)}

    def chipTouches(chip):
        return tileTouches[list(members[chip])].sum(axis=0) > 0

    def addedCost(tile, chip):                                              # Cost of tile t joining chip c: new groups pay the chip, the chip's groups pay t.
        touching = chipTouches(chip)
        tileTouching = tileTouches[tile] > 0
        return float((tileTouching & ~touching) @ tileZero[list(members[chip])].sum(axis=0) + (touching & ~tileTouching) @ tileZero[tile])

    def neighbourChips(tile, exclude):
        return {int(chipOf[other]) for other in neighbourIndices[indexPointer[tile]:indexPointer[tile + 1]] if chipOf[other] not in exclude}

    def move(tile, chip):
        members[int(chipOf[tile])].discard(tile)
        members[chip].add(tile)
        chipOf[tile] = chip

    def canGive(chip, tile):                                                # The chip stays connected without the tile.
        return _isConnected(members[chip] - {tile}, indexPointer, neighbourIndices)

    def placeTile(tile, sourceChip, banned):
        """Move `tile` out of sourceChip, making room along the shortest chain if needed. True on success."""
        direct = [chip for chip in neighbourChips(tile, banned | {sourceChip}) if len(members[chip]) < capacity]
        if direct:
            move(tile, min(direct, key=lambda chip: addedCost(tile, chip)))
            return True
        frontier = [(chip, [chip]) for chip in neighbourChips(tile, banned | {sourceChip})]
        visited = {chip for chip, _ in frontier} | banned | {sourceChip}
        while frontier:
            nextFrontier = []
            for chip, path in frontier:
                if len(path) > maximumPathLength:
                    continue
                for boundary in sorted(members[chip], key=lambda t: -len(neighbourChips(t, {chip}))):
                    for target in neighbourChips(boundary, visited - {path[-1]} | {chip} | banned | {sourceChip}):
                        if len(members[target]) < capacity and canGive(chip, boundary):
                            for position in range(len(path) - 1, -1, -1):  # Shift one tile along the chain, from its end back to its start.
                                giver = path[position]
                                receiver = target if position == len(path) - 1 else path[position + 1]
                                candidates = [t for t in members[giver] if receiver in neighbourChips(t, {giver}) and canGive(giver, t)]
                                if not candidates:
                                    return False
                                move(min(candidates, key=lambda t: addedCost(t, receiver)), receiver)
                            move(tile, path[0])
                            return True
                    for target in neighbourChips(boundary, visited | banned | {sourceChip, chip}):
                        if target not in visited:
                            visited.add(target)
                            nextFrontier.append((target, path + [target]))
            frontier = nextFrontier
        return False

    for chip in sorted(members, key=lambda chip: -len(members[chip])):     # Oversized chips shed boundary tiles.
        while len(members[chip]) > capacity:
            boundary = [t for t in members[chip] if neighbourChips(t, {chip}) and canGive(chip, t)]
            if not boundary or not placeTile(min(boundary, key=lambda t: -len(neighbourChips(t, {chip}))), chip, set()):
                break
    for face, budget in chipsPerFace.items():
        while True:
            faceChips = [chip for chip in members if members[chip] and faces[next(iter(members[chip]))] == face]
            if len(faceChips) <= budget:
                break
            dissolving = min(faceChips, key=lambda chip: len(members[chip]))
            snapshot = (chipOf.copy(), {chip: set(tiles) for chip, tiles in members.items()})
            success = True
            while members[dissolving]:
                border = [t for t in members[dissolving] if neighbourChips(t, {dissolving})]
                if not border or not placeTile(border[0], dissolving, {dissolving}):
                    success = False
                    break
            if not success:
                chipOf[:], members = snapshot[0], snapshot[1]
                break
            del members[dissolving]
    chipOf = np.unique(chipOf, return_inverse=True)[1].ravel()
    reached = {face: int(np.unique(chipOf[faces == face]).size) for face in chipsPerFace}
    return chipOf, reached


def metisChips(matrix, cost, tileKeys, chipsPerFace, capacity=defaultCapacity, neighbours=None, seed=0):
    """Partition each face into chipsPerFace[face] connected chips with METIS (pymetis, contiguous parts, imbalance allowed up to the capacity),
    cutting between dissimilar tiles: the weight of an edge is the cost-weighted Jaccard similarity of the two tiles' pattern sets (columns of M),
    so that tiles vetoing the same patterns stay together. Oversized or disconnected parts, if any, go through repairChips. Returns chipOf."""
    import pymetis
    indexPointer, neighbourIndices = tileNeighbours(tileKeys) if neighbours is None else neighbours
    faces = np.array([key[0] for key in tileKeys])
    patternWeight = np.asarray(cost, dtype=float).mean(axis=1)
    weightedColumns = matrix * patternWeight[:, None]
    chipOf = np.full(len(tileKeys), -1)
    nextChip = 0
    for face, numberOfChips in chipsPerFace.items():
        tiles = np.flatnonzero(faces == face)
        local = {int(tile): position for position, tile in enumerate(tiles)}
        adjacency, edgeWeights = [], []
        for tile in tiles:
            row, weights = [], []
            for other in neighbourIndices[indexPointer[tile]:indexPointer[tile + 1]]:
                if int(other) in local:
                    intersection = float(weightedColumns[:, tile] @ matrix[:, other])
                    union = float(weightedColumns[:, tile].sum() + weightedColumns[:, other].sum() - intersection)
                    row.append(local[int(other)])
                    weights.append(1 + int(round(1000 * (intersection / union if union > 0 else 1.0))))
            adjacency.append(row); edgeWeights.extend(weights)
        imbalance = max(int(np.floor(1000 * (capacity * numberOfChips / tiles.size - 1))), 1)   # METIS ufactor: allowed excess over the mean, per mille.
        xadj = np.concatenate([[0], np.cumsum([len(row) for row in adjacency])])
        options = pymetis.Options(ufactor=imbalance, seed=seed, contig=1)
        _, parts = pymetis.part_graph(numberOfChips, xadj=xadj, adjncy=np.concatenate(adjacency), eweights=edgeWeights, options=options)
        chipOf[tiles] = nextChip + np.asarray(parts)
        nextChip += numberOfChips
    sizes = np.bincount(chipOf)
    needsRepair = sizes.max() > capacity or not all(_isConnected(np.flatnonzero(chipOf == chip), indexPointer, neighbourIndices)
                                                    for chip in np.unique(chipOf))
    if needsRepair:
        rowGroup = np.zeros(matrix.shape[0], dtype=int)
        chipOf = repairChips(matrix, cost, tileKeys, rowGroup, chipOf, chipsPerFace, capacity, (indexPointer, neighbourIndices))[0]
    return chipOf


def mergeChipsIntoBoards(matrix, cost, rowGroup, chipOf, numberOfBoards, chipLabels=None):
    """Merge whole chips into numberOfBoards boards by exact cost (mergeUtils), the pattern groups kept. Returns the board of every chip."""
    numberOfChips = int(chipOf.max()) + 1
    start = lu.computeUnionRouting(matrix, rowGroup, chipOf)
    state = mu._GreedyState(matrix, cost, None, None if chipLabels is None else np.asarray(chipLabels)[chipOf], start)
    state.mergeTo(int(rowGroup.max()) + 1, numberOfBoards, "channelsFirst")
    tileBoard = state.channelMap().colGroup
    boardOfChip = np.zeros(numberOfChips, dtype=np.int64)
    boardOfChip[chipOf] = tileBoard
    return boardOfChip


# ----------------------------------------------------------------------------------------------------------------------------------------------
# Annealing stage
# ----------------------------------------------------------------------------------------------------------------------------------------------

@njit(cache=True)
def _blockCost(count, zero):
    return zero if count > 0 else 0.0


@njit(cache=True)
def _leavesChipConnected(tile, chip, chipOf, chipSize, indexPointer, neighbourIndices, visited, stack):
    """Whether chip `chip` stays connected once `tile` leaves it (breadth-first search over its other tiles)."""
    if chipSize[chip] <= 2:
        return True
    start = -1
    for position in range(indexPointer[tile], indexPointer[tile + 1]):
        other = neighbourIndices[position]
        if chipOf[other] == chip:
            start = other
            break
    if start < 0:
        return False
    visited[:] = False
    visited[tile] = True
    visited[start] = True
    stack[0] = start
    top, reached = 1, 1
    while top > 0:
        top -= 1
        current = stack[top]
        for position in range(indexPointer[current], indexPointer[current + 1]):
            other = neighbourIndices[position]
            if not visited[other] and chipOf[other] == chip:
                visited[other] = True
                stack[top] = other
                top += 1
                reached += 1
    return reached == chipSize[chip] - 1


@njit(cache=True)
def _anneal(roi, zeroCost, groupOf, chipOf, boardOf, chipSize, groupLabel, patternLabel, boardFace, chipFace,
            indexPointer, neighbourIndices, numberOfGroups, numberOfBoards, capacity, steps, startTemperature, endTemperature,
            moveShares, seed, traceEvery, proposalDeltas, chipBudget, chipPenalty, borderPenalty):
    """The annealing loop. With proposalDeltas of non-zero length, the cost change of every proposal is recorded there and nothing is accepted
    (a probe of the move sizes at the start state)."""
    np.random.seed(seed)
    probing = proposalDeltas.size > 0
    probed = 0
    numberOfPatterns, numberOfTiles = roi.shape
    numberOfChips = chipSize.size
    groupTile = np.zeros((numberOfGroups, numberOfTiles), dtype=np.int64)            # Patterns of group k with tile t in their RoI.
    groupTileZero = np.zeros((numberOfGroups, numberOfTiles))
    for i in range(numberOfPatterns):
        for t in range(numberOfTiles):
            groupTile[groupOf[i], t] += roi[i, t]
            groupTileZero[groupOf[i], t] += zeroCost[i, t]
    patternBoard = np.zeros((numberOfPatterns, numberOfBoards), dtype=np.int64)       # RoI tiles of pattern i on board b.
    patternBoardZero = np.zeros((numberOfPatterns, numberOfBoards))
    for t in range(numberOfTiles):
        board = boardOf[chipOf[t]]
        for i in range(numberOfPatterns):
            patternBoard[i, board] += roi[i, t]
            patternBoardZero[i, board] += zeroCost[i, t]
    blockCount = np.zeros((numberOfGroups, numberOfBoards), dtype=np.int64)
    blockZero = np.zeros((numberOfGroups, numberOfBoards))
    for i in range(numberOfPatterns):
        for b in range(numberOfBoards):
            blockCount[groupOf[i], b] += patternBoard[i, b]
            blockZero[groupOf[i], b] += patternBoardZero[i, b]
    liveChips = 0
    for c in range(numberOfChips):
        if chipSize[c] > 0:
            liveChips += 1
    total = chipPenalty * max(liveChips - chipBudget, 0)                   # Chips above the budget, if the start has any.
    for t in range(numberOfTiles):                                          # Chip borders: neighbouring tiles of different chips (each pair once).
        for position in range(indexPointer[t], indexPointer[t + 1]):
            other = neighbourIndices[position]
            if other > t and chipOf[other] != chipOf[t]:
                total += borderPenalty
    for k in range(numberOfGroups):
        for b in range(numberOfBoards):
            total += _blockCost(blockCount[k, b], blockZero[k, b])
    best = total
    bestGroupOf, bestChipOf, bestBoardOf = groupOf.copy(), chipOf.copy(), boardOf.copy()
    visited = np.zeros(numberOfTiles, dtype=np.bool_)
    stack = np.zeros(numberOfTiles, dtype=np.int64)
    traceLength = steps // traceEvery + 1
    trace = np.zeros((traceLength, 4))
    accepted = np.zeros(3, dtype=np.int64)
    proposed = np.zeros(3, dtype=np.int64)
    cooling = (endTemperature / startTemperature) ** (1.0 / max(steps - 1, 1))
    temperature = startTemperature
    chipTile = np.zeros(numberOfGroups, dtype=np.int64)
    chipZero = np.zeros(numberOfGroups)
    for step in range(steps):
        draw = np.random.random()
        if draw < moveShares[0]:                                            # A tile to a neighbouring chip.
            t = np.random.randint(numberOfTiles)
            degree = indexPointer[t + 1] - indexPointer[t]
            u = neighbourIndices[indexPointer[t] + np.random.randint(degree)] if degree > 0 else t
            fromChip, toChip = chipOf[t], chipOf[u]
            if fromChip != toChip and chipSize[toChip] < capacity and _leavesChipConnected(t, fromChip, chipOf, chipSize, indexPointer,
                                                                                            neighbourIndices, visited, stack):
                proposed[0] += 1
                fromBoard, toBoard = boardOf[fromChip], boardOf[toChip]
                delta = -chipPenalty if (chipSize[fromChip] == 1 and liveChips > chipBudget) else 0.0
                if borderPenalty > 0.0:                                     # Edges to the old chip become borders, edges to the new one stop being.
                    for position in range(indexPointer[t], indexPointer[t + 1]):
                        other = neighbourIndices[position]
                        if chipOf[other] == fromChip:
                            delta += borderPenalty
                        elif chipOf[other] == toChip:
                            delta -= borderPenalty
                if fromBoard != toBoard:
                    for k in range(numberOfGroups):
                        a, z = groupTile[k, t], groupTileZero[k, t]
                        delta += (_blockCost(blockCount[k, fromBoard] - a, blockZero[k, fromBoard] - z) - _blockCost(blockCount[k, fromBoard], blockZero[k, fromBoard])
                                  + _blockCost(blockCount[k, toBoard] + a, blockZero[k, toBoard] + z) - _blockCost(blockCount[k, toBoard], blockZero[k, toBoard]))
                if probing:
                    if probed < proposalDeltas.size:
                        proposalDeltas[probed] = delta
                        probed += 1
                elif delta <= 0.0 or np.random.random() < np.exp(-delta / temperature):
                    accepted[0] += 1
                    chipOf[t] = toChip
                    if chipSize[fromChip] == 1:
                        liveChips -= 1
                    chipSize[fromChip] -= 1
                    chipSize[toChip] += 1
                    if fromBoard != toBoard:
                        for k in range(numberOfGroups):
                            blockCount[k, fromBoard] -= groupTile[k, t]; blockZero[k, fromBoard] -= groupTileZero[k, t]
                            blockCount[k, toBoard] += groupTile[k, t]; blockZero[k, toBoard] += groupTileZero[k, t]
                        for i in range(numberOfPatterns):
                            patternBoard[i, fromBoard] -= roi[i, t]; patternBoardZero[i, fromBoard] -= zeroCost[i, t]
                            patternBoard[i, toBoard] += roi[i, t]; patternBoardZero[i, toBoard] += zeroCost[i, t]
                    total += delta
        elif draw < moveShares[0] + moveShares[1]:                          # A chip to another board.
            chip = np.random.randint(numberOfChips)
            fromBoard = boardOf[chip]
            toBoard = np.random.randint(numberOfBoards - 1) if numberOfBoards > 1 else 0
            if toBoard >= fromBoard:
                toBoard += 1
            if (chipSize[chip] > 0 and numberOfBoards > 1 and (boardFace[toBoard] < 0 or boardFace[toBoard] == chipFace[chip])):
                proposed[1] += 1
                chipTile[:] = 0
                chipZero[:] = 0.0
                for t in range(numberOfTiles):
                    if chipOf[t] == chip:
                        for k in range(numberOfGroups):
                            chipTile[k] += groupTile[k, t]
                            chipZero[k] += groupTileZero[k, t]
                delta = 0.0
                for k in range(numberOfGroups):
                    delta += (_blockCost(blockCount[k, fromBoard] - chipTile[k], blockZero[k, fromBoard] - chipZero[k]) - _blockCost(blockCount[k, fromBoard], blockZero[k, fromBoard])
                              + _blockCost(blockCount[k, toBoard] + chipTile[k], blockZero[k, toBoard] + chipZero[k]) - _blockCost(blockCount[k, toBoard], blockZero[k, toBoard]))
                if probing:
                    if probed < proposalDeltas.size:
                        proposalDeltas[probed] = delta
                        probed += 1
                elif delta <= 0.0 or np.random.random() < np.exp(-delta / temperature):
                    accepted[1] += 1
                    boardOf[chip] = toBoard
                    for k in range(numberOfGroups):
                        blockCount[k, fromBoard] -= chipTile[k]; blockZero[k, fromBoard] -= chipZero[k]
                        blockCount[k, toBoard] += chipTile[k]; blockZero[k, toBoard] += chipZero[k]
                    for t in range(numberOfTiles):
                        if chipOf[t] == chip:
                            for i in range(numberOfPatterns):
                                patternBoard[i, fromBoard] -= roi[i, t]; patternBoardZero[i, fromBoard] -= zeroCost[i, t]
                                patternBoard[i, toBoard] += roi[i, t]; patternBoardZero[i, toBoard] += zeroCost[i, t]
                    total += delta
        else:                                                               # A pattern to another group.
            i = np.random.randint(numberOfPatterns)
            fromGroup = groupOf[i]
            toGroup = np.random.randint(numberOfGroups - 1) if numberOfGroups > 1 else 0
            if toGroup >= fromGroup:
                toGroup += 1
            if numberOfGroups > 1 and (groupLabel[toGroup] < 0 or groupLabel[toGroup] == patternLabel[i]):
                proposed[2] += 1
                delta = 0.0
                for b in range(numberOfBoards):
                    a, z = patternBoard[i, b], patternBoardZero[i, b]
                    delta += (_blockCost(blockCount[fromGroup, b] - a, blockZero[fromGroup, b] - z) - _blockCost(blockCount[fromGroup, b], blockZero[fromGroup, b])
                              + _blockCost(blockCount[toGroup, b] + a, blockZero[toGroup, b] + z) - _blockCost(blockCount[toGroup, b], blockZero[toGroup, b]))
                if probing:
                    if probed < proposalDeltas.size:
                        proposalDeltas[probed] = delta
                        probed += 1
                elif delta <= 0.0 or np.random.random() < np.exp(-delta / temperature):
                    accepted[2] += 1
                    groupOf[i] = toGroup
                    for b in range(numberOfBoards):
                        blockCount[fromGroup, b] -= patternBoard[i, b]; blockZero[fromGroup, b] -= patternBoardZero[i, b]
                        blockCount[toGroup, b] += patternBoard[i, b]; blockZero[toGroup, b] += patternBoardZero[i, b]
                    for t in range(numberOfTiles):
                        groupTile[fromGroup, t] -= roi[i, t]; groupTileZero[fromGroup, t] -= zeroCost[i, t]
                        groupTile[toGroup, t] += roi[i, t]; groupTileZero[toGroup, t] += zeroCost[i, t]
                    total += delta
        if total < best - 1e-12:
            best = total
            bestGroupOf[:] = groupOf
            bestChipOf[:] = chipOf
            bestBoardOf[:] = boardOf
        if step % traceEvery == 0:
            trace[step // traceEvery, 0] = step
            trace[step // traceEvery, 1] = temperature
            trace[step // traceEvery, 2] = total
            trace[step // traceEvery, 3] = best
        temperature *= cooling
    return bestGroupOf, bestChipOf, bestBoardOf, best, trace, accepted, proposed


@dataclass
class ChipLayout:
    rowGroup: np.ndarray      # (patterns,) pattern group of each pattern
    chipOf: np.ndarray        # (tiles,) chip of each tile
    boardOfChip: np.ndarray   # (chips,) board (channel) of each chip

    def channelMap(self, matrix):
        """The ChannelMap of the layout under union routing: the channel of a tile is the board of its chip."""
        return lu.computeUnionRouting(matrix, self.rowGroup, self.boardOfChip[self.chipOf])

    @property
    def numberOfChips(self):
        return int(np.unique(self.chipOf).size)


def chipBorderLength(layout, tileKeys, neighbours=None):
    """The number of chip borders: pairs of neighbouring tiles (in the neighbour graph) that belong to different chips."""
    indexPointer, neighbourIndices = tileNeighbours(tileKeys) if neighbours is None else neighbours
    owner = layout.chipOf
    return int(sum(1 for t in range(len(tileKeys)) for other in neighbourIndices[indexPointer[t]:indexPointer[t + 1]] if other > t and owner[other] != owner[t]))


def checkLayout(layout, tileKeys, capacity=defaultCapacity, neighbours=None):
    """Raises AssertionError unless every chip is on one face, holds at most `capacity` tiles and is connected in the neighbour graph."""
    indexPointer, neighbourIndices = tileNeighbours(tileKeys) if neighbours is None else neighbours
    faces = np.array([key[0] for key in tileKeys])
    for chip in np.unique(layout.chipOf):
        members = np.flatnonzero(layout.chipOf == chip)
        assert members.size <= capacity, f"chip {chip}: {members.size} tiles"
        assert np.unique(faces[members]).size == 1, f"chip {chip}: several faces"
        seen, frontier = {members[0]}, [members[0]]
        memberSet = set(members.tolist())
        while frontier:
            current = frontier.pop()
            for other in neighbourIndices[indexPointer[current]:indexPointer[current + 1]]:
                if other in memberSet and other not in seen:
                    seen.add(other); frontier.append(other)
        assert len(seen) == members.size, f"chip {chip}: not connected"


def annealLayout(matrix, cost, tileKeys, start, numberOfGroups, numberOfBoards, capacity=defaultCapacity, steps=2_000_000, seed=0,
                 moveShares=(0.5, 0.2, 0.3), temperatures=None, startAcceptance=0.1, endTemperatureRatio=1e-3, chipBudget=None,
                 chipPenalty=None, borderPenalty=0.0, patternLabels=None, groupLabels=None, chipFaceBoards=None, neighbours=None):
    """Simulated annealing from `start` (a ChipLayout). Unless `temperatures` = (start, end) is given, the start temperature is set from the cost
    changes of moves proposed at the start state (none accepted), for each move type separately: a median uphill move of the finest type (the
    tile moves, by orders of magnitude) is kept with probability startAcceptance, so uphill chip and pattern moves are almost never kept. A random walk is
    not used for this: under union routing it drifts onto the plateau where every group is vetoed by every board, where most moves cost nothing.
    The end temperature is endTemperatureRatio times the start one. A start with more chips than chipBudget pays chipPenalty per extra chip (by
    default the cost of the whole start map, so that emptying a chip always pays); moving the last tile out of a chip removes it. borderPenalty is
    added per chip border (a pair of neighbouring tiles on different chips): the veto cost does not depend on the shape of a chip within its board,
    so a penalty much smaller than a tile move only decides between otherwise equal moves and, once the temperature is below it, makes the chips
    compact. The returned cost includes both penalties; mergeUtils.computeMaskCost gives the veto cost alone.
    patternLabels/groupLabels keep patterns in groups of their own label (for example one group per FIT face); chipFaceBoards, a face per board,
    keeps chips on boards of their own face. Returns (best ChipLayout, its cost, a dict with the trace and the acceptance per move type)."""
    indexPointer, neighbourIndices = tileNeighbours(tileKeys) if neighbours is None else neighbours
    roi = matrix.astype(np.int64)
    zeroCost = np.asarray(cost, dtype=float) * ~matrix
    faces = np.array([key[0] for key in tileKeys])
    faceIndex = {face: index for index, face in enumerate(pd.psdSides)}
    numberOfChips = int(start.chipOf.max()) + 1
    chipFace = np.full(numberOfChips, -1, dtype=np.int64)
    chipFace[start.chipOf] = [faceIndex[face] for face in faces]
    boardFace = np.full(numberOfBoards, -1, dtype=np.int64) if chipFaceBoards is None else np.array([faceIndex[face] for face in chipFaceBoards])
    patternLabel = np.zeros(matrix.shape[0], dtype=np.int64) if patternLabels is None else np.unique(patternLabels, return_inverse=True)[1].ravel()
    groupLabel = np.full(numberOfGroups, -1, dtype=np.int64)
    if groupLabels is not None:
        names = np.unique(patternLabels)
        groupLabel[:] = [int(np.flatnonzero(names == label)[0]) for label in groupLabels]
    chipSize = np.bincount(start.chipOf, minlength=numberOfChips).astype(np.int64)
    if chipPenalty is None:
        chipPenalty = float((zeroCost * (start.channelMap(matrix).vetoMask() & ~matrix)).sum()) if chipBudget is not None else 0.0
    shares = np.asarray(moveShares, dtype=float) / np.sum(moveShares)
    noDeltas = np.zeros(0)
    if temperatures is None:
        medianUphill = []
        for moveType in range(3):                                           # The three moves differ in scale by orders of magnitude.
            if shares[moveType] == 0:
                continue
            onlyThisMove = np.zeros(3); onlyThisMove[moveType] = 1.0
            deltas = np.full(20000, np.nan)
            _anneal(roi, zeroCost, start.rowGroup.copy(), start.chipOf.copy(), start.boardOfChip.copy(), chipSize.copy(), groupLabel,
                    patternLabel, boardFace, chipFace, indexPointer, neighbourIndices, numberOfGroups, numberOfBoards, capacity, 60000, 1.0, 1.0,
                    onlyThisMove, seed + 1_000_003 + moveType, 60000, deltas, numberOfChips, 0.0, 0.0)
            uphill = deltas[np.isfinite(deltas) & (deltas > 1e-12)]
            if uphill.size:
                medianUphill.append(float(np.median(uphill)))
        startTemperature = min(medianUphill) / np.log(1.0 / startAcceptance) if medianUphill else 1e-3
        temperatures = (startTemperature, startTemperature * endTemperatureRatio)
    bestGroupOf, bestChipOf, bestBoardOf, best, trace, accepted, proposed = _anneal(
        roi, zeroCost, start.rowGroup.copy(), start.chipOf.copy(), start.boardOfChip.copy(), chipSize, groupLabel, patternLabel, boardFace,
        chipFace, indexPointer, neighbourIndices, numberOfGroups, numberOfBoards, capacity, steps, temperatures[0], temperatures[1], shares, seed,
        max(steps // 400, 1), noDeltas, numberOfChips if chipBudget is None else chipBudget, chipPenalty, borderPenalty)
    _, compactChips = np.unique(bestChipOf, return_inverse=True)
    usedChips = np.unique(bestChipOf)
    layout = ChipLayout(bestGroupOf, compactChips.ravel(), bestBoardOf[usedChips])
    return layout, best, dict(trace=trace, accepted=accepted, proposed=proposed, temperatures=temperatures)


def findChipLayout(matrix, cost, tileKeys, numberOfGroups, numberOfBoards, capacity=defaultCapacity, chipsPerFace=None, groupSeeds=("greedy", "wardJaccard"),
                   chipSeeds=("repairedMerge", "metis", "serpentine"), annealingSeeds=(0, 1, 2), steps=20_000_000, relativeBorderPenalty=1e-4,
                   groupLabels=None, chipFaceBoards=None):
    """The whole search: for every pattern-group seed (groups-first exact-cost greedy, or the Ward-scaled Jaccard merge) and chip seed
    ("repairedMerge": neighbouring tiles merged by exact cost, then repaired to the budget; "merge": the same without repair, which usually ends
    above the budget; "metis": a similarity-weighted graph partition of each face; "serpentine": a contiguous tiling that ignores the cost), boards
    by the exact-cost greedy, then simulated annealing from several random seeds. The chip border penalty is relativeBorderPenalty times sum C M
    per border: at 1e-4 it halves the chip borders on the HERD showcase for about 0.5% of veto cost (compact chips instead of winding ones). chipsPerFace is the
    chip budget per face (default: the minimum). groupLabels, one label per pattern (for example its FIT face), keeps each group to one label, the
    groups being named by the labels in order (numberOfGroups must equal the number of labels). Returns (best ChipLayout within the budget,
    [dict per candidate: seeds, chips, cost before and after annealing, run time])."""
    import time
    import jaccardMerge as jm
    neighbours = tileNeighbours(tileKeys)
    chipsPerFace = minimumChips(tileKeys, capacity) if chipsPerFace is None else chipsPerFace
    chipBudget = sum(chipsPerFace.values())
    reference = float((cost * matrix).sum())
    relative = lambda layout: mu.computeMaskCost(matrix, layout.channelMap(matrix), cost)[0] / reference
    numberOfLosslessChannels = lu.findLosslessChannelMap(matrix).numberOfChannels
    groupStarts = {}
    for groupSeed in groupSeeds:
        if groupSeed == "greedy":
            groupStarts[groupSeed] = mu.mergeGreedily(matrix, numberOfGroups, numberOfLosslessChannels, cost, order="groupsFirst",
                                                      rowLabels=groupLabels).rowGroup
        else:
            groupStarts[groupSeed] = jm.mergeByJaccard(matrix, numberOfGroups, numberOfLosslessChannels, cost=cost, variant=groupSeed,
                                                       order="groupsFirst", rowLabels=groupLabels)[0].rowGroup
    labelOfGroup = None
    if groupLabels is not None:
        names = list(np.unique(groupLabels))
        for groupSeed, rowGroup in groupStarts.items():                    # Name the groups by their label, in the order of the labels.
            order = [names.index(groupLabels[np.flatnonzero(rowGroup == k)[0]]) for k in range(rowGroup.max() + 1)]
            groupStarts[groupSeed] = np.array(order)[rowGroup]
        labelOfGroup = names
    candidates, best, bestCost = [], None, np.inf
    for groupSeed, rowGroup in groupStarts.items():
        for chipSeed in chipSeeds:
            if chipSeed == "serpentine":
                chipOf = serpentineChips(tileKeys, chipsPerFace)
            elif chipSeed == "metis":
                chipOf = metisChips(matrix, cost, tileKeys, chipsPerFace, capacity, neighbours)
            else:
                chipOf = mergeTilesIntoChips(matrix, cost, tileKeys, rowGroup, chipsPerFace, capacity, neighbours)[0]
                if chipSeed == "repairedMerge":
                    chipOf = repairChips(matrix, cost, tileKeys, rowGroup, chipOf, chipsPerFace, capacity, neighbours)[0]
            boardOfChip = mergeChipsIntoBoards(matrix, cost, rowGroup, chipOf, numberOfBoards)
            start = ChipLayout(rowGroup, chipOf, boardOfChip)
            for annealingSeed in annealingSeeds:
                began = time.time()
                layout, _, info = annealLayout(matrix, cost, tileKeys, start, numberOfGroups, numberOfBoards, capacity, steps, annealingSeed,
                                               chipBudget=chipBudget, borderPenalty=relativeBorderPenalty * reference,
                                               patternLabels=groupLabels, groupLabels=labelOfGroup,
                                               chipFaceBoards=chipFaceBoards, neighbours=neighbours)
                checkLayout(layout, tileKeys, capacity, neighbours)
                result = dict(groupSeed=groupSeed, chipSeed=chipSeed, annealingSeed=annealingSeed, startChips=start.numberOfChips,
                              chips=layout.numberOfChips, startCost=relative(start), cost=relative(layout),
                              chipBorders=chipBorderLength(layout, tileKeys, neighbours), seconds=time.time() - began)
                candidates.append(result)
                if layout.numberOfChips <= chipBudget and result["cost"] < bestCost:
                    best, bestCost = layout, result["cost"]
    return best, candidates
