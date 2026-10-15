"""Recovers the per-tile PSD activation rates of the thesis (Fig. 4.9, all species, 'NewData' flux model) from the published raster image, the only
copy left: the analysed events were on a cluster that is no longer accessible.

    python tools/extractPsdRates.py [image.png] [--output data/psdTileRates.csv]

Method: each tile's cell is located with the unfolding of AeffCalc/global_analyses/psd_rates.py (nbVPsdRepFromID, in the earlier form that drew the
thesis figure: Xneg at X = 8 - r and Yneg at Y = 8 - r; the later form shifts both by one cell and leaves those cells white). The dominant colour of the
cell interior is matched to the colorbar's own pixels, whose ticks (1.5 ... 2.1) give the scale; extrapolated, the colorbar spans 1.469 to 2.149
s^-1 cm^-2, and the top matches the maximum in AeffCalc/global_analyses/max_rates.txt (2.1496). Cells under the 'Preliminary' watermark (grey at alpha
0.3) are un-blended before matching. The output lists the rate per unit area, the matching residual (RGB distance) and how the cell was read.
Precision: one colorbar pixel is 0.0018 s^-1 cm^-2.
"""
import argparse
import pathlib
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "analytic"))
import patternDump as pd

defaultImage = pathlib.Path.home() / "MDPIpaper/instruments-2886224/Raw_PSD_rates_map_all_particles_NewData.png"
panelLeft = {0: 48, 1: 556}                  # px. Left edges of the two panels (layer 0, layer 1), found from the axes frames.
panelBottom, panelSize, gridCells = 471, 373, 34   # px, px, cells. hSelPsd2D is Int(0, 34) on both axes.
colorbarColumn, colorbarRows = 963, (99, 471)       # px. A column inside the colorbar, and its vertical extent.
tickRows, tickValues = (125, 454), (2.1, 1.5)       # px and s^-1 cm^-2: the outermost labelled ticks.
watermarkGrey, watermarkAlpha = np.array([128.0, 128.0, 128.0]), 0.3   # matplotlib "gray" at alpha 0.3 (psd_rates.py, the 'Preliminary' text).


def cellOf(face, tileId):
    """(X, Y) of a tile in the unfolded histogram, as the thesis figure was drawn."""
    sideV, topSize = 9, 14
    if face == "TOP":
        return sideV + 1 + tileId // topSize, sideV + 1 + tileId % topSize
    r, f = tileId % sideV, tileId // sideV
    return {"Xpos": (sideV + 2 + topSize + r, sideV + 1 + f), "Xneg": (sideV - r - 1, sideV + topSize - f),
            "Ypos": (sideV + topSize - f, sideV + 2 + topSize + r), "Yneg": (sideV + 1 + f, sideV - r - 1)}[face]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("image", nargs="?", type=pathlib.Path, default=defaultImage)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parent.parent / "data" / "psdTileRates.csv")
    args = parser.parse_args()
    image = np.asarray(Image.open(args.image).convert("RGB")).astype(float)
    rows = np.arange(*colorbarRows)
    lookup = image[rows, colorbarColumn]
    lookupValue = tickValues[0] + (rows - tickRows[0]) * (tickValues[1] - tickValues[0]) / (tickRows[1] - tickRows[0])

    def decode(colour):
        distance = np.linalg.norm(lookup - colour, axis=1)
        k = np.argmin(distance)
        return lookupValue[k], distance[k]

    cell = panelSize / gridCells
    lines = ["face,layer,tileId,ratePerArea,residual,reading"]
    for layer, left in panelLeft.items():
        for face in pd.psdSides:
            for tileId in range(pd.tilesPerLayer[face]):
                x, y = cellOf(face, tileId)
                cx, cy = left + (x + 0.5) * cell, panelBottom - (y + 0.5) * cell
                pixels = image[int(cy - 3):int(cy + 4), int(cx - 3):int(cx + 4)].reshape(-1, 3)
                colours, counts = np.unique(pixels, axis=0, return_counts=True)
                colour = colours[np.argmax(counts)]
                value, residual = decode(colour)
                reading = "direct"
                if residual > 12:
                    unblendedValue, unblendedResidual = decode((colour - watermarkAlpha * watermarkGrey) / (1 - watermarkAlpha))
                    if unblendedResidual < residual:
                        value, residual, reading = unblendedValue, unblendedResidual, "unblended"
                if colour.min() > 245:
                    reading = "white"
                lines.append(f"{face},{layer},{tileId},{value:.4f},{residual:.1f},{reading}")
    args.output.write_text("\n".join(lines) + "\n")
    values = np.array([float(line.split(",")[3]) for line in lines[1:]])
    readings = [line.split(",")[5] for line in lines[1:]]
    residuals = np.array([float(line.split(",")[4]) for line in lines[1:]])
    print(f"{len(values)} tiles -> {args.output}: mean {values.mean():.3f}, std {values.std():.3f}, min {values.min():.3f}, max {values.max():.3f} s^-1 cm^-2; "
          f"readings {dict(zip(*np.unique(readings, return_counts=True)))}; residual > 12 in {(residuals > 12).sum()} cells")


if __name__ == "__main__":
    main()
