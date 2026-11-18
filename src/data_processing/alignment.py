"""Chromatogram retention-time alignment module (Manual and Cross-Correlation).

Strictly preserves original retention-time and intensity arrays; computes and stores
additive retention-time shifts (`t_corrected = t_raw + shift_min`) separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple
import numpy as np
from scipy.signal import correlate, correlation_lags

from src.data_processing.chromatogram import Chromatogram
from src.data_processing.normalization import estimate_visual_baseline


ALIGNMENT_MODES = (
    "No alignment",
    "Manual alignment",
    "Automatic alignment (Cross-Correlation)",
)


@dataclass
class AlignmentEstimate:
    """Proposed or applied alignment parameters for a single chromatogram."""

    trace_id: str
    sample_name: str
    file_name: str
    detector: str
    reference_trace_id: str
    reference_sample_name: str
    estimated_shift_min: float
    correlation_score: float
    region_start_min: float
    region_end_min: float
    max_shift_min: float
    is_reference: bool = False
    warning: str = ""


def estimate_cross_correlation_shift(
    reference: Chromatogram,
    target: Chromatogram,
    region_start_min: Optional[float] = None,
    region_end_min: Optional[float] = None,
    max_shift_min: float = 1.5,
    detrend_baseline: bool = True,
) -> AlignmentEstimate:
    """Estimate the optimal retention-time shift (in minutes) to align `target` to `reference`.

    When the returned `estimated_shift_min` is added to `target.raw_retention_time`,
    features in `target` align with corresponding features in `reference`:
        `t_corrected = t_raw + estimated_shift_min`
    """
    if reference.trace_id == target.trace_id:
        t_start = float(region_start_min) if region_start_min is not None else reference.rt_min
        t_end = float(region_end_min) if region_end_min is not None else reference.rt_max
        return AlignmentEstimate(
            trace_id=target.trace_id,
            sample_name=target.sample_display_name,
            file_name=target.file_name,
            detector=target.assigned_detector,
            reference_trace_id=reference.trace_id,
            reference_sample_name=reference.sample_display_name,
            estimated_shift_min=0.0,
            correlation_score=1.0,
            region_start_min=t_start,
            region_end_min=t_end,
            max_shift_min=max_shift_min,
            is_reference=True,
        )

    # Determine overlapping retention-time window
    overlap_start = max(reference.rt_min, target.rt_min)
    overlap_end = min(reference.rt_max, target.rt_max)

    t_start = max(overlap_start, float(region_start_min)) if region_start_min is not None else overlap_start
    t_end = min(overlap_end, float(region_end_min)) if region_end_min is not None else overlap_end

    warning_msg = ""
    if reference.assigned_detector != target.assigned_detector:
        warning_msg = (
            f"Detector mismatch ({target.assigned_detector} aligned against {reference.assigned_detector})"
        )

    if t_end <= t_start or (t_end - t_start) < 0.05:
        return AlignmentEstimate(
            trace_id=target.trace_id,
            sample_name=target.sample_display_name,
            file_name=target.file_name,
            detector=target.assigned_detector,
            reference_trace_id=reference.trace_id,
            reference_sample_name=reference.sample_display_name,
            estimated_shift_min=0.0,
            correlation_score=0.0,
            region_start_min=t_start,
            region_end_min=t_end,
            max_shift_min=max_shift_min,
            warning="Insufficient overlapping retention-time region for alignment.",
        )

    # Choose uniform time grid step from the finer of the two chromatograms
    dt_ref = reference.sampling_interval_min or 0.001667
    dt_tgt = target.sampling_interval_min or 0.001667
    dt = min(dt_ref, dt_tgt)
    dt = max(dt, 1e-5)

    n_grid = int(round((t_end - t_start) / dt)) + 1
    if n_grid < 10:
        n_grid = 200
        dt = (t_end - t_start) / (n_grid - 1)
    grid_t = np.linspace(t_start, t_end, n_grid)

    # Interpolate both signals onto the uniform region grid (raw arrays are untouched)
    y_ref = np.interp(grid_t, reference.raw_retention_time, reference.raw_intensity)
    y_tgt = np.interp(grid_t, target.raw_retention_time, target.raw_intensity)

    # Internal baseline detrending and zero-meaning for cross-correlation only
    if detrend_baseline and n_grid >= 30:
        y_ref_proc = y_ref - estimate_visual_baseline(grid_t, y_ref, window_minutes=1.0)
        y_tgt_proc = y_tgt - estimate_visual_baseline(grid_t, y_tgt, window_minutes=1.0)
    else:
        y_ref_proc = y_ref.copy()
        y_tgt_proc = y_tgt.copy()

    y_ref_proc = y_ref_proc - np.mean(y_ref_proc)
    y_tgt_proc = y_tgt_proc - np.mean(y_tgt_proc)

    norm_ref = float(np.linalg.norm(y_ref_proc))
    norm_tgt = float(np.linalg.norm(y_tgt_proc))
    if norm_ref <= 1e-15 or norm_tgt <= 1e-15:
        return AlignmentEstimate(
            trace_id=target.trace_id,
            sample_name=target.sample_display_name,
            file_name=target.file_name,
            detector=target.assigned_detector,
            reference_trace_id=reference.trace_id,
            reference_sample_name=reference.sample_display_name,
            estimated_shift_min=0.0,
            correlation_score=0.0,
            region_start_min=t_start,
            region_end_min=t_end,
            max_shift_min=max_shift_min,
            warning="Flat or zero-variance signal in selected alignment region.",
        )

    corr_full = correlate(y_ref_proc, y_tgt_proc, mode="full", method="auto")
    corr_norm = corr_full / (norm_ref * norm_tgt)
    lags = correlation_lags(y_ref_proc.size, y_tgt_proc.size, mode="full")
    shifts_min = lags * dt

    # Restrict search to [-max_shift_min, +max_shift_min]
    max_shift = max(float(max_shift_min), dt)
    valid_mask = np.abs(shifts_min) <= max_shift
    if not np.any(valid_mask):
        valid_mask = np.ones_like(shifts_min, dtype=bool)

    valid_indices = np.flatnonzero(valid_mask)
    sub_corr = corr_norm[valid_indices]
    best_local_idx = int(np.argmax(sub_corr))
    best_idx = int(valid_indices[best_local_idx])
    best_score = float(corr_norm[best_idx])

    # Parabolic 3-point sub-sample interpolation when interior to valid window
    fractional_lag = float(lags[best_idx])
    if 0 < best_idx < (corr_norm.size - 1) and (0 < best_local_idx < valid_indices.size - 1):
        c_left = float(corr_norm[best_idx - 1])
        c_mid = float(corr_norm[best_idx])
        c_right = float(corr_norm[best_idx + 1])
        denom = c_left - 2.0 * c_mid + c_right
        if abs(denom) > 1e-12 and c_mid >= c_left and c_mid >= c_right:
            delta = 0.5 * (c_left - c_right) / denom
            if abs(delta) <= 0.5:
                fractional_lag += delta

    estimated_shift = float(np.clip(fractional_lag * dt, -max_shift, max_shift))

    return AlignmentEstimate(
        trace_id=target.trace_id,
        sample_name=target.sample_display_name,
        file_name=target.file_name,
        detector=target.assigned_detector,
        reference_trace_id=reference.trace_id,
        reference_sample_name=reference.sample_display_name,
        estimated_shift_min=estimated_shift,
        correlation_score=best_score,
        region_start_min=t_start,
        region_end_min=t_end,
        max_shift_min=max_shift,
        is_reference=False,
        warning=warning_msg,
    )


def estimate_batch_alignment(
    reference: Chromatogram,
    targets: Sequence[Chromatogram],
    region_start_min: Optional[float] = None,
    region_end_min: Optional[float] = None,
    max_shift_min: float = 1.5,
) -> List[AlignmentEstimate]:
    """Compute proposed cross-correlation shifts for a list of chromatograms."""
    estimates: List[AlignmentEstimate] = []
    for chrom in targets:
        est = estimate_cross_correlation_shift(
            reference=reference,
            target=chrom,
            region_start_min=region_start_min,
            region_end_min=region_end_min,
            max_shift_min=max_shift_min,
        )
        estimates.append(est)
    return estimates


def apply_alignment_estimates(
    chromatograms: Sequence[Chromatogram],
    estimates: Sequence[AlignmentEstimate],
    method_label: str = "Automatic (Cross-Correlation)",
) -> None:
    """Apply a list of approved AlignmentEstimates to their corresponding Chromatograms."""
    est_by_id: Dict[str, AlignmentEstimate] = {e.trace_id: e for e in estimates}
    for chrom in chromatograms:
        if chrom.trace_id in est_by_id:
            est = est_by_id[chrom.trace_id]
            chrom.apply_shift(
                shift_min=est.estimated_shift_min,
                method=method_label,
                reference_id=est.reference_trace_id,
            )


def apply_manual_shifts(
    chromatograms: Sequence[Chromatogram],
    shifts_by_trace_id: Dict[str, float],
    reference_trace_id: Optional[str] = None,
) -> None:
    """Apply user-specified manual retention-time shifts (in minutes)."""
    for chrom in chromatograms:
        if chrom.trace_id == reference_trace_id:
            chrom.apply_shift(0.0, method="Manual (Reference)", reference_id=reference_trace_id)
        elif chrom.trace_id in shifts_by_trace_id:
            chrom.apply_shift(
                shift_min=float(shifts_by_trace_id[chrom.trace_id]),
                method="Manual",
                reference_id=reference_trace_id,
            )


def reset_all_alignments(chromatograms: Sequence[Chromatogram]) -> None:
    """Reset all chromatograms to their original unaligned retention times."""
    for chrom in chromatograms:
        chrom.reset_alignment()
