"""Plot helpers: a pattern's veto region on the PSD faces, the incidence matrix, and the merging curves.

Colours: one sequential blue ramp for magnitudes (layers in the RoI, matrix entries) and an ordinal blue ramp for ordered series (the
number of tile channels L); steps of the dataviz reference palette, the ordinal ones no lighter than step 250 so that every curve stays visible.

Tile layout on each face, as read from the PsdGeoParams volume IDs in the dump: top ID = row x 14 + column (14 x 14); side ID = H x 9 + V
(13 x 9). The orientation of rows and columns on the real detector is not checked: the panels show the index grids, not a drawing of HERD.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

import patternDump as pd

inkPrimary, inkMuted, gridColour = "#1f1f1e", "#6b6a63", "#e4e3dc"
sequentialRamp = ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # dataviz palette, blue 100..700 (first step: near-surface tint).
ordinalRamp = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]      # dataviz palette, blue 250..700.
layerColours = ListedColormap(["#f4f8fd", "#86b6ef", "#1c5cab"])                      # 0, 1 or 2 PSD layers in the RoI.


def setStyle():
    plt.rcParams.update({"figure.dpi": 110, "font.size": 11, "axes.edgecolor": inkMuted, "axes.labelcolor": inkPrimary, "xtick.color": inkMuted,
                         "ytick.color": inkMuted, "axes.grid": True, "grid.color": gridColour, "grid.linewidth": 0.8, "axes.spines.top": False,
                         "axes.spines.right": False, "lines.linewidth": 2.0, "legend.frameon": False})


def faceGrid(side):
    """(rows, columns, index function) of the tile grid of one PSD face."""
    if side == "TOP":
        return 14, 14, lambda tileId: (tileId // 14, tileId % 14)
    return 9, 13, lambda tileId: (tileId % 9, tileId // 9)


def roiFaceImages(tiles):
    """Per PSD face, an image counting the PSD layers (0, 1, 2) whose tile at that grid position is in the set of tiles."""
    images = {}
    for side in pd.psdSides:
        rows, columns, index = faceGrid(side)
        image = np.zeros((rows, columns), dtype=int)
        for tileSide, layer, tileId in tiles:
            if tileSide == side:
                image[index(tileId)] += 1
        images[side] = image
    return images


def plotRoiOnFaces(pattern, axes=None, title=None):
    """The veto region of one pattern on the five PSD faces (shade: number of PSD layers in the RoI)."""
    if axes is None:
        figure, axes = plt.subplots(1, 5, figsize=(16, 3.4), constrained_layout=True)
    images = roiFaceImages(pattern.tiles)
    for drawingAxis, side in zip(axes, pd.psdSides):
        drawingAxis.imshow(images[side], cmap=layerColours, vmin=0, vmax=2, origin="lower", aspect="equal")
        drawingAxis.set(title=f"PSD {side}: {int((images[side] > 0).sum())} positions", xticks=[], yticks=[])
        drawingAxis.grid(False)
        drawingAxis.title.set_fontsize(10)
    axes[0].figure.suptitle(title or f"RoI of {pattern.name} (FIT {pattern.fitSide}, plane {pattern.firstPlane}, mat {pattern.firstMat}, code {pattern.code})", color=inkPrimary)
    return axes


def plotIncidence(matrix, drawingAxis, title):
    drawingAxis.imshow(matrix, cmap=ListedColormap([sequentialRamp[0], sequentialRamp[4]]), aspect="auto", interpolation="nearest")
    drawingAxis.set(title=title, xlabel="PSD tile (column)", ylabel="trigger pattern (row)")
    drawingAxis.grid(False)


# ---------------------------------------------------------------------------------------------------------------------------------------------
# 3D sketches of the PSD (geometryUtils: schematic layout, adjacency inferred from the dump). Only the outer layer is drawn; a tile position is
# coloured by what holds for either of its two layers.

categoricalColours = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]   # dataviz palette slots 1..8, fixed order.
neutralTile, edgeColour = "#eeede6", "#ffffff"


def plotPsdCube(drawingAxis, colourOf, azimuth=-60, elevation=24, title=None, legend=None):
    """Draw the PSD faces turned towards the camera. colourOf: {(face, tileId): colour}; other positions are neutral. legend: [(colour, label)]."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import geometryUtils as gu
    view = np.array([np.cos(np.radians(elevation)) * np.cos(np.radians(azimuth)), np.cos(np.radians(elevation)) * np.sin(np.radians(azimuth)), np.sin(np.radians(elevation))])
    quads, colours = [], []
    for face in pd.psdSides:
        if np.dot(gu.faceNormals[face], view) <= 0.05:
            continue
        for tileId in range(pd.tilesPerLayer[face]):
            quads.append(gu.tileQuad(face, 0, tileId))
            colours.append(colourOf.get((face, tileId), neutralTile))
    drawingAxis.add_collection3d(Poly3DCollection(quads, facecolors=colours, edgecolors=edgeColour, linewidths=0.4))
    for face in pd.psdSides:
        if np.dot(gu.faceNormals[face], view) > 0.05:
            anchor = np.array([0, 0, 22]) if face == "TOP" else gu.faceNormals[face] * (gu.halfWidth + 30) + np.array([0, 0, -45])
            drawingAxis.text(*anchor, face, color=inkMuted, fontsize=9, ha="center")
    limit = gu.halfWidth + 5
    drawingAxis.set(xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-100, 40))
    drawingAxis.set_box_aspect((2 * limit, 2 * limit, 140))
    drawingAxis.view_init(elev=elevation, azim=azimuth)
    drawingAxis.set_axis_off()
    if title:
        drawingAxis.set_title(title, fontsize=10, color=inkPrimary)
    if legend:
        from matplotlib.patches import Patch
        drawingAxis.legend(handles=[Patch(facecolor=colour, label=label) for colour, label in legend], loc="lower left", fontsize=8, frameon=False)


def colourBySets(tileKeys, columnSets):
    """{(face, tileId): colour} from [(bool per column, colour)], later sets painting over earlier ones; a position counts if either layer is set."""
    colourOf = {}
    for mask, colour in columnSets:
        for column in np.flatnonzero(mask):
            face, layer, tileId = tileKeys[column]
            colourOf[(face, tileId)] = colour
    return colourOf


# ---------------------------------------------------------------------------------------------------------------------------------------------
# Readable views of the incidence matrix: a geometric order for the patterns, a seriated order for the tiles, RoI similarity, and the unfolded PSD
# of the thesis (Fig. 4.9: TOP in the centre, each lateral face folded out along the edge it shares with TOP).

def patternOrder(patterns):
    """Rows sorted by FIT face, view (first plane), first mat, then slope (the pattern code)."""
    return np.array(sorted(range(len(patterns)), key=lambda i: (pd.psdSides.index(patterns[i].fitSide), patterns[i].firstPlane, patterns[i].firstMat, patterns[i].code)))


def tileOrder(matrix, tileKeys, rowOrder):
    """Columns grouped by PSD face (then layer), and within a face sorted by the mean position, in rowOrder, of the patterns whose RoI holds the tile
    (a one-step seriation: tiles used by neighbouring patterns sit together). Returns (order, [(face, first column, last column)])."""
    position = np.empty(len(rowOrder)); position[rowOrder] = np.arange(len(rowOrder))
    used = matrix.sum(axis=0)
    centroid = np.where(used > 0, (matrix * position[:, None]).sum(axis=0) / np.maximum(used, 1), np.inf)
    order, blocks = [], []
    for face in pd.psdSides:
        columns = [j for j, key in enumerate(tileKeys) if key[0] == face]
        columns.sort(key=lambda j: (centroid[j], tileKeys[j][1], tileKeys[j][2]))
        blocks.append((face, len(order), len(order) + len(columns) - 1))
        order.extend(columns)
    return np.array(order), blocks


def plotOrderedIncidence(drawingAxis, matrix, patterns, tileKeys, title):
    rows = patternOrder(patterns)
    columns, blocks = tileOrder(matrix, tileKeys, rows)
    drawingAxis.imshow(matrix[np.ix_(rows, columns)], cmap=ListedColormap([sequentialRamp[0], sequentialRamp[4]]), aspect="auto", interpolation="nearest")
    for face, first, last in blocks[1:]:
        drawingAxis.axvline(first - 0.5, color=inkMuted, linewidth=0.8)
    drawingAxis.set_xticks([(first + last) / 2 for _, first, last in blocks], [face for face, _, _ in blocks])
    rowFaces = [patterns[i].fitSide for i in rows]
    starts = [n for n in range(len(rows)) if n == 0 or rowFaces[n] != rowFaces[n - 1]]
    for n in starts[1:]:
        drawingAxis.axhline(n - 0.5, color=inkMuted, linewidth=0.8)
    drawingAxis.set_yticks([(a + b) / 2 for a, b in zip(starts, starts[1:] + [len(rows)])], [rowFaces[n] for n in starts])
    drawingAxis.set(title=title, xlabel="PSD tiles, by face (within a face: ordered by the patterns that use them)", ylabel="3IR patterns, by FIT face")
    drawingAxis.grid(False)


def roiJaccard(matrix):
    """|T_a ∩ T_b| / |T_a ∪ T_b| for every pair of RoIs."""
    m = matrix.astype(float)
    common = m @ m.T
    sizes = m.sum(axis=1)
    return common / (sizes[:, None] + sizes[None, :] - common)


def unfoldedPosition(face, tileId):
    """Centre of a tile, in tile units, on the unfolded cross."""
    import geometryUtils as gu
    i, j = gu.gridIndex(face, tileId)
    if face == "TOP":
        return np.array([i - 6.5, j - 6.5])
    normal = gu.faceNormals[face][:2]
    along = np.array([-normal[1], normal[0]])
    return normal * (7.0 + j + 0.5) + along * (i - 6)


def plotUnfoldedPsd(drawingAxis, valueOf, layer, norm, colourMap, title=None):
    """Thesis Fig. 4.9 layout: one square per tile of the given layer, coloured by valueOf[(face, layer, tileId)]."""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    squares, values = [], []
    for face in pd.psdSides:
        for tileId in range(pd.tilesPerLayer[face]):
            x, y = unfoldedPosition(face, tileId)
            squares.append(Rectangle((x - 0.5, y - 0.5), 1.0, 1.0))
            values.append(valueOf.get((face, layer, tileId), 0.0))
    collection = PatchCollection(squares, cmap=colourMap, norm=norm, edgecolor="white", linewidth=0.4)
    collection.set_array(np.array(values))
    drawingAxis.add_collection(collection)
    drawingAxis.set(xlim=(-17, 17), ylim=(-17, 17), aspect="equal", xticks=[], yticks=[])
    for face, (x, y) in {"Xpos": (12, -7.9), "Xneg": (-12, 7.9), "Ypos": (7.9, 12), "Yneg": (-7.9, -12)}.items():
        drawingAxis.text(x, y, face, ha="center", va="center", fontsize=8, color=inkMuted)
    drawingAxis.grid(False)
    for spine in drawingAxis.spines.values():
        spine.set_visible(False)
    if title:
        drawingAxis.set_title(title, fontsize=10)
    return collection
