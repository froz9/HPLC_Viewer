"""Validation routines for imported chromatographic data."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence
import numpy as np


VALID_DETECTORS = ("UV", "PDA", "ELSD", "Other / Unknown")


@dataclass
class ValidationReport:
    """Comprehensive validation report for a chromatographic channel."""

    file_name: str
    channel: str = ""
    is_valid: bool = True
    status: str = "OK"  # "OK", "Warning", "Error"
    is_empty: bool = False
    is_corrupted: bool = False
    missing_rt_count: int = 0
    missing_intensity_count: int = 0
    non_numeric_count: int = 0
    duplicate_time_count: int = 0
    unsorted_time_count: int = 0
    negative_intensity_count: int = 0
    unequal_sampling_interval: bool = False
    ascii_rounding_jitter: bool = False
    detector_ambiguous: bool = False
    sampling_interval_median_min: Optional[float] = None
    sampling_interval_min_min: Optional[float] = None
    sampling_interval_max_min: Optional[float] = None
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    info_notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert validation report to a dictionary suitable for tables or JSON."""
        return {
            "file_name": self.file_name,
            "channel": self.channel,
            "is_valid": self.is_valid,
            "status": self.status,
            "is_empty": self.is_empty,
            "is_corrupted": self.is_corrupted,
            "missing_rt_count": self.missing_rt_count,
            "missing_intensity_count": self.missing_intensity_count,
            "non_numeric_count": self.non_numeric_count,
            "duplicate_time_count": self.duplicate_time_count,
            "unsorted_time_count": self.unsorted_time_count,
            "negative_intensity_count": self.negative_intensity_count,
            "unequal_sampling_interval": self.unequal_sampling_interval,
            "ascii_rounding_jitter": self.ascii_rounding_jitter,
            "detector_ambiguous": self.detector_ambiguous,
            "sampling_interval_median_min": self.sampling_interval_median_min,
            "warnings": list(self.warnings),
            "errors": list(self.errors),
            "info_notes": list(self.info_notes),
        }

    @property
    def summary_text(self) -> str:
        """Concise human-readable status message for the import table."""
        if self.errors:
            return f"Error: {'; '.join(self.errors)}"
        if self.warnings:
            return f"Warning: {'; '.join(self.warnings)}"
        if self.info_notes:
            return f"OK ({'; '.join(self.info_notes)})"
        return "OK"


def validate_chromatogram_arrays(
    file_name: str,
    retention_time: Sequence[Any],
    intensity: Sequence[Any],
    channel: str = "",
    missing_rt_count: int = 0,
    missing_intensity_count: int = 0,
    non_numeric_count: int = 0,
    detector_ambiguous: bool = False,
    relative_dt_tolerance: float = 0.015,
) -> ValidationReport:
    """Validate retention-time and intensity arrays without modifying raw values.

    Parameters
    ----------
    file_name : str
        Source file name.
    retention_time : Sequence[float]
        Extracted numeric retention time observations.
    intensity : Sequence[float]
        Extracted numeric intensity observations.
    channel : str
        Channel identifier.
    missing_rt_count : int
        Number of rows with missing retention-time values encountered during parsing.
    missing_intensity_count : int
        Number of rows with missing intensity values encountered during parsing.
    non_numeric_count : int
        Number of rows with non-numeric tokens encountered during parsing.
    detector_ambiguous : bool
        True if the detector type could not be unambiguously identified.
    relative_dt_tolerance : float
        Relative deviation threshold above which sampling intervals are flagged as unequal
        (allows ~0.5% Empower 7-significant-digit ASCII rounding jitter by default).
    """
    report = ValidationReport(
        file_name=file_name,
        channel=channel,
        missing_rt_count=missing_rt_count,
        missing_intensity_count=missing_intensity_count,
        non_numeric_count=non_numeric_count,
        detector_ambiguous=detector_ambiguous,
    )

    rt_arr = np.asarray(retention_time, dtype=np.float64)
    int_arr = np.asarray(intensity, dtype=np.float64)

    # Check for NaN / Inf inside numeric arrays
    nan_rt = int(np.sum(~np.isfinite(rt_arr)))
    nan_int = int(np.sum(~np.isfinite(int_arr)))
    if nan_rt > 0:
        report.missing_rt_count += nan_rt
    if nan_int > 0:
        report.missing_intensity_count += nan_int

    if report.missing_rt_count > 0:
        report.warnings.append(f"{report.missing_rt_count} missing retention-time value(s)")
    if report.missing_intensity_count > 0:
        report.warnings.append(f"{report.missing_intensity_count} missing intensity value(s)")
    if report.non_numeric_count > 0:
        report.warnings.append(f"{report.non_numeric_count} non-numeric observation row(s) skipped")
    if report.detector_ambiguous:
        report.warnings.append("Detector type ambiguous (assigned 'Other / Unknown')")

    if rt_arr.size == 0 or int_arr.size == 0:
        report.is_empty = True
        report.is_valid = False
        report.status = "Error"
        report.errors.append("No valid chromatographic data points found")
        return report

    if rt_arr.size != int_arr.size:
        report.is_corrupted = True
        report.is_valid = False
        report.status = "Error"
        report.errors.append(
            f"Mismatched array lengths: {rt_arr.size} time points vs {int_arr.size} intensity points"
        )
        return report

    # Filter finite pairs only for interval/sign statistics (raw arrays remain untouched)
    finite_mask = np.isfinite(rt_arr) & np.isfinite(int_arr)
    rt_fin = rt_arr[finite_mask]
    int_fin = int_arr[finite_mask]

    if rt_fin.size == 0:
        report.is_empty = True
        report.is_valid = False
        report.status = "Error"
        report.errors.append("All observations are NaN or infinite")
        return report

    # Negative intensity check (reported as informational/warning, never deleted)
    neg_count = int(np.sum(int_fin < 0.0))
    report.negative_intensity_count = neg_count
    if neg_count > 0:
        min_val = float(np.min(int_fin))
        report.info_notes.append(
            f"{neg_count} negative baseline point(s) preserved (min={min_val:.4g})"
        )

    # Time ordering, duplicates, and sampling interval checks
    if rt_fin.size >= 2:
        dt = np.diff(rt_fin)
        dup_count = int(np.sum(dt == 0.0))
        unsorted_count = int(np.sum(dt < 0.0))
        report.duplicate_time_count = dup_count
        report.unsorted_time_count = unsorted_count

        if dup_count > 0:
            report.warnings.append(
                f"{dup_count} duplicate retention-time point(s) detected (preserved)"
            )
        if unsorted_count > 0:
            report.warnings.append(
                f"Unsorted retention times detected ({unsorted_count} decreasing step(s))"
            )

        pos_dt = dt[dt > 0.0]
        if pos_dt.size > 0:
            med_dt = float(np.median(pos_dt))
            min_dt = float(np.min(pos_dt))
            max_dt = float(np.max(pos_dt))
            report.sampling_interval_median_min = med_dt
            report.sampling_interval_min_min = min_dt
            report.sampling_interval_max_min = max_dt

            if med_dt > 0.0:
                max_abs_dev = float(np.max(np.abs(pos_dt - med_dt)))
                rel_dev = max_abs_dev / med_dt
                if rel_dev > relative_dt_tolerance or unsorted_count > 0 or dup_count > 0:
                    report.unequal_sampling_interval = True
                    report.warnings.append(
                        f"Unequal sampling intervals (min={min_dt:.6g}, median={med_dt:.6g}, max={max_dt:.6g} min)"
                    )
                elif max_abs_dev > 1e-9:
                    report.ascii_rounding_jitter = True

    if report.errors:
        report.status = "Error"
        report.is_valid = False
    elif report.warnings:
        report.status = "Warning"
        report.is_valid = True
    else:
        report.status = "OK"
        report.is_valid = True

    return report
