"""Optional intensity normalization and baseline estimation for visualization.

All functions return new NumPy arrays for plotting and NEVER modify the
original experimental intensity arrays.
"""

from __future__ import annotations

from typing import Tuple
import numpy as np
from scipy.ndimage import minimum_filter1d, uniform_filter1d


NORMALIZATION_MODES = (
    "None",
    "Min-Max [0, 1]",
    "Max Intensity = 100%",
    "Z-Score",
)


def compute_visual_intensity(
    raw_intensity: np.ndarray,
    mode: str = "None",
) -> Tuple[np.ndarray, str]:
    """Return a normalized copy of `raw_intensity` for display only.

    Parameters
    ----------
    raw_intensity : np.ndarray
        Original immutable detector intensity array.
    mode : str
        One of `NORMALIZATION_MODES`.

    Returns
    -------
    Tuple[np.ndarray, str]
        (transformed_intensity_array, y_axis_unit_suffix)
    """
    arr = np.asarray(raw_intensity, dtype=np.float64)
    if arr.size == 0 or mode == "None":
        return (arr.copy(), "")

    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return (arr.copy(), "")

    if mode == "Min-Max [0, 1]":
        ymin = float(np.min(finite))
        ymax = float(np.max(finite))
        span = ymax - ymin
        if span <= 1e-15:
            return (np.zeros_like(arr), "Normalized [0–1]")
        return ((arr - ymin) / span, "Normalized [0–1]")

    if mode == "Max Intensity = 100%":
        ymax = float(np.max(np.abs(finite)))
        if ymax <= 1e-15:
            return (np.zeros_like(arr), "Relative Intensity (%)")
        return ((arr / ymax) * 100.0, "Relative Intensity (%)")

    if mode == "Z-Score":
        mean = float(np.mean(finite))
        std = float(np.std(finite))
        if std <= 1e-15:
            return (np.zeros_like(arr), "Z-Score (σ)")
        return ((arr - mean) / std, "Z-Score (σ)")

    return (arr.copy(), "")


def estimate_visual_baseline(
    retention_time: np.ndarray,
    intensity: np.ndarray,
    window_minutes: float = 1.5,
) -> np.ndarray:
    """Compute a non-destructive rolling minimum-envelope baseline for visual overlay.

    Uses a wide morphological minimum filter followed by a moving-average smoothing
    pass so chromatographic peaks are ignored while slow baseline drift is tracked.
    """
    rt = np.asarray(retention_time, dtype=np.float64)
    y = np.asarray(intensity, dtype=np.float64)
    if y.size < 5:
        return np.zeros_like(y)

    dt_vals = np.diff(rt)
    pos_dt = dt_vals[dt_vals > 0]
    med_dt = float(np.median(pos_dt)) if pos_dt.size > 0 else 0.001667
    win_pts = max(15, int(round(window_minutes / max(med_dt, 1e-6))))
    win_pts = min(win_pts, max(5, y.size // 2))

    eroded = minimum_filter1d(y, size=win_pts, mode="nearest")
    smoothed = uniform_filter1d(eroded, size=win_pts, mode="nearest")
    return smoothed
