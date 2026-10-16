"""Per-tile probabilities that a PSD tile vetoes a good photon: backsplash from the photon's own shower, and accidental (time-coincident) hits from the
charged and neutral background that fires every tile at about 190 Hz.

    pi_t = 1 - (1 - b_t) (1 - q_t),      q_t = 1 - exp(-R_t tau),      R_t = r_t A

- b_t: backsplash probability of tile t per triggered photon, measured on the surviving Geant4 photons (tools/measureBacksplash.py ->
  data/psdTileBacksplash.csv), weighted to an E^-2 spectrum by default (with it the full-PSD veto keeps 0.90 of fiducial photons; thesis Table 4.2:
  0.89) or as generated (log-uniform). A uniform value can be given instead.
- r_t: activation rate per unit area of each tile (thesis Fig. 4.9, all species, 'NewData' flux model), recovered from the published figure by
  tools/extractPsdRates.py into data/psdTileRates.csv; mean 1.88 s^-1 cm^-2 (AeffCalc max_rates.txt: 1.89 from the events).
- A: tile area, 100 cm^2 (10 x 10 cm tiles, thesis section 1.8.1; tileArea = 100 in AeffCalc/analysis_and_plots/psd_rates.py, the version that made
  the figure; the later 5 x 30 cm tile design, 150 cm^2, is a different PSD and is not modelled here).
- tau: the coincidence window of the L1 veto. Not fixed by the thesis: a parameter (a placeholder of 100 ns by default).

Accidental vetoes scale with tau and with the number of tiles in the mask, so mask excess costs photons twice: backsplash at high energy, accidentals
at all energies. Naming: camelCase everywhere. Units: s, cm.
"""
import pathlib

import numpy as np

defaultRatesPath = pathlib.Path(__file__).resolve().parent.parent / "data" / "psdTileRates.csv"
tileArea = 100.0                 # cm^2. 10 x 10 cm tiles (thesis, section 1.8.1).
defaultBacksplashPath = pathlib.Path(__file__).resolve().parent.parent / "data" / "psdTileBacksplash.csv"
defaultWindow = 100e-9           # s. Placeholder coincidence window of the L1 veto.


def loadTileRates(tileKeys, path=defaultRatesPath):
    """Activation rate per unit area r_t (s^-1 cm^-2) for each tile key, in column order."""
    table = {}
    for line in pathlib.Path(path).read_text().splitlines()[1:]:
        if line:
            face, layer, tileId, rate = line.split(",")[:4]
            table[(face, int(layer), int(tileId))] = float(rate)
    return np.array([table[key] for key in tileKeys])


def loadBacksplash(tileKeys, weighting="spectrumE2", path=defaultBacksplashPath):
    """Measured backsplash probability b_t per tile key ('spectrumE2' or 'logUniform' weighting of the photons)."""
    lines = pathlib.Path(path).read_text().splitlines()
    column = lines[0].split(",").index(weighting)
    table = {}
    for line in lines[1:]:
        if line:
            parts = line.split(",")
            table[(parts[0], int(parts[1]), int(parts[2]))] = float(parts[column])
    return np.array([table[key] for key in tileKeys])


def tileVetoProbability(tileKeys, window=defaultWindow, backsplash="spectrumE2", ratesPath=defaultRatesPath):
    """pi_t for each tile: the probability that the tile vetoes a good photon (backsplash or an accidental hit inside the window). backsplash is a
    weighting of the measured map ('spectrumE2', 'logUniform'), or a number (uniform per tile), or an array (one value per tile)."""
    accidental = -np.expm1(-loadTileRates(tileKeys, ratesPath) * tileArea * window)
    if isinstance(backsplash, str):
        backsplash = loadBacksplash(tileKeys, backsplash)
    return 1.0 - (1.0 - np.asarray(backsplash, dtype=float)) * (1.0 - accidental)


def costMatrix(numberOfPatterns, vetoProbability):
    """C_it = -ln(1 - pi_t), the same for every pattern: the expected log-loss of photons when tile t joins a mask, so that the cost of a mask is -ln
    of its survival probability. An approximation of C_it = -w_i ln[(1 - b_it)(1 - q_t)]: the backsplash really depends on the trigger pattern (b_it,
    a pattern x tile map) and patterns fire at different photon rates w_i; here b_it is averaged over all triggers and w_i = 1 (too few photons)."""
    return np.tile(-np.log1p(-np.asarray(vetoProbability)), (numberOfPatterns, 1))
