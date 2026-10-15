"""Reads the trigger-pattern -> veto-tile map dumped by ProduceRoIMap (legacy/aeffcalc/algos/ProduceRoIMap.cpp, "dumpPatterns true") into a
boolean incidence matrix: one row per trigger pattern, one column per PSD tile, True when the tile is in the pattern's veto region (RoI).

The dump in data/patterns.txt is the compact-v2 HERD layout (legacy/aeffcalc/datacards/geometry_compact.mac) with the "Calo" trigger
configuration: 3 hits in 3 same-view planes, at most one mat between consecutive planes, the codes 544, 546, 564, 566 removed and a
calorimeter cube required in the RoI (requireCalo). It only contains patterns starting in planes 8 and 9 (startLayer 8: an inference from the dump,
the datacard next to it says 0). Each RoI is the "direct" one of UlegGeomRoi (legacy/aeffcalc/algos/UlegGeomPsdRoi.h): the PSD tiles any corner of
which lies inside the double cone of straight lines through the first and the last mat of the pattern, above the first mat.

Pattern names follow TriggerPattern::FullID: <FIT side><first plane><first mat, 2 digits><pattern code>, the code being one digit per plane,
5 + (mat - previous mat), the first digit always 5. Example: TOP800555 is a vertical pattern starting at plane 8, mat 0, of the top FIT.
Tiles are keyed (PSD side, PSD layer, tile ID), the tile ID being the volume ID of PsdGeoParams (0..195 on top, 0..116 on each side).

Naming: camelCase everywhere, no unit suffixes.
"""
import pathlib
import re
from dataclasses import dataclass

import numpy as np

dataDirectory = pathlib.Path(__file__).resolve().parent.parent / "data"
defaultDumpPath = dataDirectory / "patterns_allPlanes.txt"     # Every start plane and code, regenerated 2026-10-15 (data/README.md).
legacyDumpPath = dataDirectory / "patterns.txt"                # The February 2023 dump: planes 8-9, "Calo" configuration, earlier RoI code.
zigzagCodes = ("546", "564")                                   # Patterns that step out and back to the first mat ("NoPlusMinus" in analyzeMC_top.eaconf).
steepCodes = ("544", "566")                                    # Straight but steep: the third hit two mats from the first ("plus-plus"). Allowed only by the
                                                               # thesis' "Allow 2x adjacent" variant (Table 4.1); the 3IR definition keeps the third hit in
                                                               # the module adjacent to the first (section 4.2.1).
twoStepCodes = zigzagCodes + steepCodes                        # ignoreIDs {544,546,564,566} of the "Calo" and "NoUpwards" triggers ("NoPlusPlus").

def pointsUpwards(pattern):
    """An upgoing pattern: the particle moves towards the top, i.e. arrives from below (the Earth). The thesis calls these "upward-pointing"
    (section 4.2.1) and TriggerArrays.cpp removes them with removeUpwardsPatterns (TriggerPattern::PointsUpwards): on a lateral face the inner mat
    lies above the outer one. Only the even planes of the side faces are segmented in height
    (8 mats); there the mat index grows downwards, as the RoIs show: patterns with a net negative mat step have RoIs low on the face that never reach
    TOP (tracks from below), those with a net positive step reach the top PSD (tracks from above). FINDINGS 7. Steep 544 patterns on the even planes of the side faces are upgoing (their RoIs never reach TOP; six of them have empty RoIs)."""
    return pattern.fitSide != "TOP" and pattern.firstPlane % 2 == 0 and sum(pattern.matSteps) < 0


triggerSelections = {
    "all": dict(),
    "calo": dict(requireCalo=True),                                       # A calorimeter cube in the RoI (requireCalo).
    "caloAdjacent": dict(requireCalo=True, excludeCodes=twoStepCodes),    # Third hit adjacent to the first, no zig-zags: the thesis' CALO-intersecting set.
    "calo8": dict(requireCalo=True, excludeCodes=twoStepCodes, planes=(8, 9)),   # The analysis card's "Calo" trigger (layer_start 8), as the 2023 dump.
    "thesis": dict(requireCalo=True, excludeCodes=twoStepCodes, removeUpwards=True),   # CALO intersect + no upwards (Table 4.1; upgoing patterns removed), all planes: the default.
    "thesis8": dict(requireCalo=True, excludeCodes=twoStepCodes, removeUpwards=True, planes=(8, 9)),   # The card's "NoUpwards" trigger (layer_start 8).
}


def selectPatterns(patterns, requireCalo=False, excludeCodes=(), planes=None, removeUpwards=False):
    """The patterns of a trigger configuration: requireCalo keeps RoIs holding a calorimeter cube, excludeCodes drops pattern codes, planes keeps
    only the given first planes, removeUpwards drops the upward-pointing patterns (pointsUpwards: upgoing, the particle moving towards the top). triggerSelections names the configurations used in the notebooks."""
    return [p for p in patterns if (not requireCalo or p.caloCubes > 0) and p.code not in excludeCodes and (planes is None or p.firstPlane in planes)
            and not (removeUpwards and pointsUpwards(p))]
psdSides = ("TOP", "Xpos", "Xneg", "Ypos", "Yneg")   # Order used for the columns. Names as in Herd::RefFrame::SideName.
psdLayers = 2                                         # nTopTileLayers 2, nSideTileLayers 2 (geometry_compact.mac).
tilesPerLayer = {"TOP": 14 * 14, "Xpos": 13 * 9, "Xneg": 13 * 9, "Ypos": 13 * 9, "Yneg": 13 * 9}   # nTopTilesAlongX x Y = 14 x 14; nSideTilesAlongH x V = 13 x 9 (geometry_compact.mac).
codeLength = 3                                        # Digits of the pattern code: patternSize 3 (ProduceRoIMap default).

_patternLine = re.compile(r"^Trigger pattern: (TOP|Xpos|Xneg|Ypos|Yneg)(\d+)$")
_sideLine = re.compile(r"^PSD side (\w+): (\d+) tiles:$")
_layerLine = re.compile(r"^Layer (\d+): \(([\d,]*)\)$")
_cubesLine = re.compile(r"^Cubes in RoI: (\d+)$")
_totalLine = re.compile(r"^Tiles in RoI: (\d+)$")


@dataclass
class TriggerPattern:
    name: str            # FullID, e.g. "TOP800555"
    fitSide: str
    firstPlane: int
    firstMat: int
    code: str            # e.g. "555"
    tiles: frozenset     # {(psdSide, psdLayer, tileID)}
    caloCubes: int       # Calorimeter cubes in the same RoI (from the dump; not used for the veto)

    @property
    def matSteps(self):
        """Mat change between consecutive planes, e.g. (0, +1) for code 556."""
        return tuple(int(digit) - 5 for digit in self.code[1:])


def parsePatternDump(path=defaultDumpPath):
    """Every pattern of the dump, in file order. Checks the tile counts written in the dump against the parsed lists."""
    patterns, current = [], None

    def close():
        if current is not None:
            assert len(current["tiles"]) == current["declaredTotal"], f"{current['name']}: {len(current['tiles'])} tiles parsed, {current['declaredTotal']} declared"
            patterns.append(TriggerPattern(current["name"], current["fitSide"], current["firstPlane"], current["firstMat"], current["code"], frozenset(current["tiles"]), current["caloCubes"]))

    psdSide = None
    for rawLine in pathlib.Path(path).read_text().splitlines():
        line = rawLine.strip()
        if match := _patternLine.match(line):
            close()
            fitSide, digits = match.groups()
            current = dict(name=fitSide + digits, fitSide=fitSide, firstPlane=int(digits[:-(2 + codeLength)]), firstMat=int(digits[-(2 + codeLength):-codeLength]),
                           code=digits[-codeLength:], tiles=set(), caloCubes=-1, declaredTotal=-1)
        elif current is None:
            continue
        elif match := _totalLine.match(line):
            current["declaredTotal"] = int(match.group(1))
        elif match := _sideLine.match(line):
            psdSide = match.group(1)
        elif match := _layerLine.match(line):
            layer, ids = int(match.group(1)), match.group(2)
            current["tiles"].update((psdSide, layer, int(tileId)) for tileId in ids.split(",") if tileId)
        elif match := _cubesLine.match(line):
            current["caloCubes"] = int(match.group(1))
    close()
    return patterns


def allTileKeys():
    """Every PSD tile of the layout, (side, layer, ID), in column order."""
    return [(side, layer, tileId) for side in psdSides for layer in range(psdLayers) for tileId in range(tilesPerLayer[side])]


@dataclass
class Incidence:
    matrix: np.ndarray      # (patterns, tiles) bool
    patterns: list          # TriggerPattern per row
    tileKeys: list          # (psdSide, psdLayer, tileID) per column

    @property
    def shape(self):
        return self.matrix.shape

    def restrict(self, rows=None, columns=None):
        rows = np.arange(self.shape[0]) if rows is None else np.asarray(rows)
        columns = np.arange(self.shape[1]) if columns is None else np.asarray(columns)
        return Incidence(self.matrix[np.ix_(rows, columns)], [self.patterns[i] for i in rows], [self.tileKeys[j] for j in columns])

    def forFitSides(self, *fitSides):
        return self.restrict(rows=[i for i, pattern in enumerate(self.patterns) if pattern.fitSide in fitSides])

    def withoutEmptyColumns(self):
        """Drops the tiles that are in no RoI: they never veto, whatever the channel map."""
        return self.restrict(columns=np.flatnonzero(self.matrix.any(axis=0)))


def buildIncidence(patterns=None, tileKeys=None):
    """The pattern x tile incidence matrix. All tiles of the layout by default (columns with no RoI included)."""
    patterns = parsePatternDump() if patterns is None else patterns
    tileKeys = allTileKeys() if tileKeys is None else tileKeys
    column = {key: j for j, key in enumerate(tileKeys)}
    matrix = np.zeros((len(patterns), len(tileKeys)), dtype=bool)
    for i, pattern in enumerate(patterns):
        matrix[i, [column[key] for key in pattern.tiles]] = True
    return Incidence(matrix, list(patterns), list(tileKeys))


def evenOddBlocks(patterns):
    """[(evenRows, oddRows)] per FIT face and double layer: the patterns starting on plane 2n and on plane 2n + 1, the two views the trigger pairs
    (TriggerArrays.cpp pairs plane k with k + 1 of the same double layer)."""
    blocks = []
    for side in psdSides:
        for doubleLayer in sorted({p.firstPlane // 2 for p in patterns if p.fitSide == side}):
            even = [i for i, p in enumerate(patterns) if p.fitSide == side and p.firstPlane == 2 * doubleLayer]
            odd = [i for i, p in enumerate(patterns) if p.fitSide == side and p.firstPlane == 2 * doubleLayer + 1]
            if even and odd:
                blocks.append((np.array(even), np.array(odd)))
    return blocks
