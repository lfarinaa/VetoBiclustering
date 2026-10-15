"""Schematic 3D layout of the PSD tiles, for sketches and for the face-based veto strategies of the thesis (section 4.2.2).

The tile grids are those of the volume IDs in the dump (patternDump.py): TOP ID = r x 14 + c, side ID = h x 9 + v. How each grid sits on the cube is
inferred from the dump itself, from which edge of each face the RoIs of the neighbouring FIT faces reach (FINDINGS 5):

- every side face touches TOP along v = 0 (v grows downwards);
- on TOP, r grows towards +X (Xpos borders r = 13, Xneg r = 0) and c towards +Y (Ypos borders c = 13, Yneg c = 0);
- on a side face of outward normal n = (nx, ny, 0), h grows along u = (-ny, nx, 0): Xpos h = 0 borders Yneg and h = 12 Ypos, Ypos h = 0 borders Xpos, and so on.

Sizes: 10 x 10 cm tiles (thesis section 1.8.1); top face 14 x 14 tiles, sides 13 x 9 (geometry_compact.mac). The two staggered layers, the gaps and the
real offsets between faces are not drawn: the layout is a sketch with the right adjacency, not the HerdSoftware geometry.
Naming: camelCase everywhere. Units: cm.
"""
import numpy as np

import patternDump as pd

tileSize = 10.0                                   # cm. "10 x 10 cm2 scintillator tiles" (thesis, section 1.8.1).
faceNormals = {"TOP": np.array([0.0, 0.0, 1.0]), "Xpos": np.array([1.0, 0.0, 0.0]), "Xneg": np.array([-1.0, 0.0, 0.0]),
               "Ypos": np.array([0.0, 1.0, 0.0]), "Yneg": np.array([0.0, -1.0, 0.0])}
gridShape = {"TOP": (14, 14), "Xpos": (13, 9), "Xneg": (13, 9), "Ypos": (13, 9), "Yneg": (13, 9)}   # (r, c) on TOP, (h, v) on the sides.
halfWidth = 14 * tileSize / 2.0                   # The side faces stand at the edge of the top face (schematic).
layerOffset = 1.5                                 # cm. Layer 1 drawn this far inside layer 0 (placeholder: the drawing only needs them apart).


def gridIndex(face, tileId):
    """(r, c) on TOP, (h, v) on a side face."""
    return (tileId // 14, tileId % 14) if face == "TOP" else (tileId // 9, tileId % 9)


def alongAxis(face):
    """Direction of growing h on a side face."""
    normal = faceNormals[face]
    return np.array([-normal[1], normal[0], 0.0])


def tileFrame(face, layer, tileId):
    """(centre, first axis, second axis) of a tile, the axes as half-size vectors spanning the tile."""
    i, j = gridIndex(face, tileId)
    inward = -faceNormals[face] * layer * layerOffset
    if face == "TOP":
        centre = np.array([(i - 6.5) * tileSize, (j - 6.5) * tileSize, 0.0]) + inward
        return centre, np.array([tileSize / 2, 0, 0]), np.array([0, tileSize / 2, 0])
    u = alongAxis(face)
    centre = faceNormals[face] * halfWidth + u * (i - 6) * tileSize + np.array([0, 0, -(j + 0.5) * tileSize]) + inward
    return centre, u * tileSize / 2, np.array([0, 0, tileSize / 2])


def tileQuad(face, layer, tileId):
    """The four corners of a tile, for Poly3DCollection."""
    centre, a, b = tileFrame(face, layer, tileId)
    return [centre - a - b, centre + a - b, centre + a + b, centre - a + b]


def rowsFromEdge(face, tileId, neighbour):
    """How many tile rows separate a tile of `face` from the edge it shares with face `neighbour` (0 = the row along the edge; None if not adjacent)."""
    if face == neighbour or np.dot(faceNormals[face], faceNormals[neighbour]) < -0.5:
        return None
    i, j = gridIndex(face, tileId)
    if face == "TOP":
        normal = faceNormals[neighbour]
        q = (i - 6.5) * normal[0] + (j - 6.5) * normal[1]
        return int(round(6.5 - q))
    if neighbour == "TOP":
        return j
    q = (i - 6) * np.dot(alongAxis(face), faceNormals[neighbour])
    return int(round(6 - q))


def faceExtendedMask(tileKeys, fitSide, rows):
    """The thesis' hybrid veto region for the patterns of one FIT face: every tile of the PSD face in front of it, plus `rows` rows of each adjacent
    face along the shared edge (section 4.2.2, "up to two adjacent rows or columns of tiles from neighboring faces"). One bool per tile key."""
    mask = np.zeros(len(tileKeys), dtype=bool)
    for column, (face, layer, tileId) in enumerate(tileKeys):
        if face == fitSide:
            mask[column] = True
        else:
            distance = rowsFromEdge(face, tileId, fitSide)
            mask[column] = distance is not None and distance < rows
    return mask
