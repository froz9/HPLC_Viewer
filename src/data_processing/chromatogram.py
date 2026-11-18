"""Core data structures for representing chromatograms and 3D PDA spectra.

Enforces strict separation between immutable raw data, processed retention-time
shifts, and user-adjustable visualization settings.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple
import numpy as np
import pandas as pd

from src.data_import.file_validation import (
    VALID_DETECTORS,
    ValidationReport,
    validate_chromatogram_arrays,
)


_WAVELENGTH_PATTERN = re.compile(
    r"(?P<wl>\d{3}(?:\.\d+)?)\s*nm(?:\s*@\s*(?P<bw>\d+(?:\.\d+)?)\s*nm)?",
    re.IGNORECASE,
)


def identify_detector_and_wavelength(
    channel: str,
    units: str = "",
    is_3d_pda: bool = False,
) -> Tuple[str, Optional[float], Optional[float], bool]:
    """Identify detector category and wavelength from channel and unit metadata.

    Parameters
    ----------
    channel : str
        Channel string from file metadata (e.g. '2998 Ch1 254nm@1.2nm', 'ELSD Signal').
    units : str
        Detector response units (e.g. 'AU', 'LSU', 'mAU').
    is_3d_pda : bool
        Whether the file is a 3D multi-wavelength PDA matrix.

    Returns
    -------
    Tuple[str, Optional[float], Optional[float], bool]
        (detector_type, wavelength_nm, bandwidth_nm, is_ambiguous)
        where detector_type is one of ('UV', 'PDA', 'ELSD', 'Other / Unknown').
    """
    ch_clean = (channel or "").strip()
    un_clean = (units or "").strip()
    ch_upper = ch_clean.upper()
    un_upper = un_clean.upper()

    wavelength_nm: Optional[float] = None
    bandwidth_nm: Optional[float] = None
    match = _WAVELENGTH_PATTERN.search(ch_clean)
    if match:
        wavelength_nm = float(match.group("wl"))
        if match.group("bw"):
            bandwidth_nm = float(match.group("bw"))

    # 1. ELSD detection
    if "ELSD" in ch_upper or un_upper == "LSU" or "LIGHT SCATTER" in ch_upper:
        return ("ELSD", None, None, False)

    # 2. 3D PDA spectrum detection
    if is_3d_pda:
        return ("PDA", wavelength_nm, bandwidth_nm, False)

    # 3. Explicit PDA channel mention
    if "PDA" in ch_upper or "DAD" in ch_upper or ch_clean == "2998":
        return ("PDA", wavelength_nm, bandwidth_nm, False)

    # 4. Waters 2998 / 996 extracted single-wavelength channel or explicit UV/TUV
    if (
        "2998" in ch_upper
        or "2487" in ch_upper
        or "2489" in ch_upper
        or "TUV" in ch_upper
        or "UV" in ch_upper
        or un_upper in ("AU", "MAU")
        or wavelength_nm is not None
    ):
        return ("UV", wavelength_nm, bandwidth_nm, False)

    # 5. Ambiguous detector: do not silently infer
    return ("Other / Unknown", wavelength_nm, bandwidth_nm, True)


def build_trace_id(file_name: str, channel: str, wavelength_nm: Optional[float] = None) -> str:
    """Create a deterministic unique identifier for an imported chromatogram trace."""
    wl_token = f"{wavelength_nm:.2f}" if wavelength_nm is not None else "none"
    raw_key = f"{file_name}::{channel}::{wl_token}"
    digest = hashlib.sha1(raw_key.encode("utf-8")).hexdigest()[:10]
    safe_file = re.sub(r"[^A-Za-z0-9_.-]", "_", file_name)
    return f"{safe_file}_{digest}"


class Chromatogram:
    """Represents a single chromatographic channel with immutable raw observations."""

    def __init__(
        self,
        file_name: str,
        retention_time: Sequence[float],
        intensity: Sequence[float],
        sample_name: str = "",
        injection_id: str = "1",
        detector: str = "Other / Unknown",
        channel: str = "",
        units: str = "",
        wavelength_nm: Optional[float] = None,
        bandwidth_nm: Optional[float] = None,
        raw_metadata: Optional[Dict[str, str]] = None,
        validation_report: Optional[ValidationReport] = None,
        color: str = "#0072B2",
        trace_id: Optional[str] = None,
        is_from_3d_pda: bool = False,
    ) -> None:
        self.file_name: str = file_name
        self.original_sample_name: str = sample_name if sample_name else file_name
        self.sample_display_name: str = self.original_sample_name
        self.injection_id: str = str(injection_id) if injection_id else "1"

        if detector not in VALID_DETECTORS:
            detector = "Other / Unknown"
        self.original_detector: str = detector
        self.assigned_detector: str = detector

        self.channel: str = channel
        self.units: str = units
        self.wavelength_nm: Optional[float] = wavelength_nm
        self.bandwidth_nm: Optional[float] = bandwidth_nm
        self.raw_metadata: Dict[str, str] = dict(raw_metadata) if raw_metadata else {}
        self.is_from_3d_pda: bool = is_from_3d_pda

        self.trace_id: str = trace_id or build_trace_id(file_name, channel, wavelength_nm)

        # Immutable raw arrays
        rt_arr = np.array(retention_time, dtype=np.float64, copy=True)
        int_arr = np.array(intensity, dtype=np.float64, copy=True)
        rt_arr.flags.writeable = False
        int_arr.flags.writeable = False
        self._raw_retention_time: np.ndarray = rt_arr
        self._raw_intensity: np.ndarray = int_arr

        # Validation report
        if validation_report is None:
            validation_report = validate_chromatogram_arrays(
                file_name=file_name,
                retention_time=self._raw_retention_time,
                intensity=self._raw_intensity,
                channel=channel,
                detector_ambiguous=(detector == "Other / Unknown"),
            )
        self.validation_report: ValidationReport = validation_report

        # Processed alignment state (separated from raw data)
        self.time_shift_min: float = 0.0
        self.is_aligned: bool = False
        self.alignment_method: str = "None"
        self.alignment_reference_id: Optional[str] = None

        # Visualization settings
        self.color: str = color
        self.visible: bool = True

    @property
    def raw_retention_time(self) -> np.ndarray:
        """Read-only view of the original retention times (in minutes)."""
        return self._raw_retention_time

    @property
    def raw_intensity(self) -> np.ndarray:
        """Read-only view of the original detector intensities."""
        return self._raw_intensity

    @property
    def corrected_retention_time(self) -> np.ndarray:
        """Retention-time array with alignment shift applied (`t_raw + time_shift_min`)."""
        if not self.is_aligned or self.time_shift_min == 0.0:
            return self._raw_retention_time
        shifted = self._raw_retention_time + float(self.time_shift_min)
        shifted.flags.writeable = False
        return shifted

    @property
    def active_retention_time(self) -> np.ndarray:
        """Return corrected retention time when aligned, otherwise raw retention time."""
        return self.corrected_retention_time if self.is_aligned else self._raw_retention_time

    def apply_shift(
        self,
        shift_min: float,
        method: str = "Manual",
        reference_id: Optional[str] = None,
    ) -> None:
        """Apply a retention-time shift (in minutes) without modifying raw data."""
        self.time_shift_min = float(shift_min)
        self.is_aligned = True
        self.alignment_method = method
        self.alignment_reference_id = reference_id

    def reset_alignment(self) -> None:
        """Revert chromatogram to its original unaligned state."""
        self.time_shift_min = 0.0
        self.is_aligned = False
        self.alignment_method = "None"
        self.alignment_reference_id = None

    @property
    def n_points(self) -> int:
        """Number of data points in the chromatogram."""
        return int(self._raw_retention_time.size)

    @property
    def rt_min(self) -> float:
        """Minimum original retention time (minutes)."""
        if self._raw_retention_time.size == 0:
            return 0.0
        return float(np.nanmin(self._raw_retention_time))

    @property
    def rt_max(self) -> float:
        """Maximum original retention time (minutes)."""
        if self._raw_retention_time.size == 0:
            return 0.0
        return float(np.nanmax(self._raw_retention_time))

    @property
    def rt_range_str(self) -> str:
        """Formatted retention-time range string."""
        if self._raw_retention_time.size == 0:
            return "N/A"
        return f"{self.rt_min:.2f} – {self.rt_max:.2f} min"

    @property
    def sampling_interval_min(self) -> Optional[float]:
        """Median sampling interval in minutes."""
        return self.validation_report.sampling_interval_median_min

    @property
    def sampling_interval_str(self) -> str:
        """Formatted sampling interval (in min and Hz)."""
        dt = self.sampling_interval_min
        if dt is None or dt <= 0:
            return "N/A"
        hz = 1.0 / (dt * 60.0)
        return f"{dt:.6f} min ({hz:.1f} Hz)"

    @property
    def wavelength_str(self) -> str:
        """Formatted wavelength string when available."""
        if self.wavelength_nm is None:
            return "—"
        if self.bandwidth_nm is not None:
            return f"{self.wavelength_nm:g} nm (@{self.bandwidth_nm:g} nm)"
        return f"{self.wavelength_nm:g} nm"

    @property
    def import_status(self) -> str:
        """Status string from the validation report."""
        return self.validation_report.summary_text

    @property
    def display_label(self) -> str:
        """Human-readable label for legends and selectors."""
        det_info = self.assigned_detector
        if self.wavelength_nm is not None:
            det_info = f"{self.assigned_detector} {self.wavelength_nm:g} nm"
        return f"{self.sample_display_name} [{det_info}]"

    def to_summary_dict(self) -> Dict[str, Any]:
        """Return a summary dictionary matching Section 4.3 import table columns."""
        return {
            "trace_id": self.trace_id,
            "File Name": self.file_name,
            "Sample Name": self.sample_display_name,
            "Original Sample Name": self.original_sample_name,
            "Injection": self.injection_id,
            "Detector": self.assigned_detector,
            "Channel": self.channel or "—",
            "Wavelength (nm)": round(self.wavelength_nm, 2) if self.wavelength_nm is not None else None,
            "Units": self.units or "—",
            "Data Points": self.n_points,
            "RT Range": self.rt_range_str,
            "Sampling Interval": self.sampling_interval_str,
            "Status": self.import_status,
            "Visible": self.visible,
            "Color": self.color,
            "Time Shift (min)": round(self.time_shift_min, 5),
        }

    def to_export_dataframe(self) -> pd.DataFrame:
        """Build a tabular DataFrame for CSV export preserving original and aligned data."""
        data: Dict[str, Any] = {
            "sample_identifier": self.sample_display_name,
            "original_sample_name": self.original_sample_name,
            "file_name": self.file_name,
            "injection_id": self.injection_id,
            "detector": self.assigned_detector,
            "channel": self.channel,
            "wavelength_nm": self.wavelength_nm if self.wavelength_nm is not None else np.nan,
            "units": self.units,
            "original_retention_time_min": self._raw_retention_time,
            "original_intensity": self._raw_intensity,
        }
        if self.is_aligned:
            data["corrected_retention_time_min"] = self.corrected_retention_time
            data["applied_shift_min"] = self.time_shift_min
            data["alignment_method"] = self.alignment_method
        else:
            data["corrected_retention_time_min"] = self._raw_retention_time
            data["applied_shift_min"] = 0.0
            data["alignment_method"] = "None"
        return pd.DataFrame(data)


@dataclass
class PDA3DSpectrum:
    """Stores a 3D Photodiode Array (time x wavelength) dataset from Empower .arw."""

    file_name: str
    sample_name: str
    injection_id: str
    channel: str
    units: str
    raw_metadata: Dict[str, str]
    retention_time: np.ndarray
    wavelengths_nm: np.ndarray
    absorbance_matrix: np.ndarray  # shape: (n_times, n_wavelengths)

    def __post_init__(self) -> None:
        self.retention_time.flags.writeable = False
        self.wavelengths_nm.flags.writeable = False
        self.absorbance_matrix.flags.writeable = False

    @property
    def n_times(self) -> int:
        return int(self.retention_time.size)

    @property
    def n_wavelengths(self) -> int:
        return int(self.wavelengths_nm.size)

    @property
    def wl_min(self) -> float:
        return float(np.min(self.wavelengths_nm)) if self.wavelengths_nm.size else 0.0

    @property
    def wl_max(self) -> float:
        return float(np.max(self.wavelengths_nm)) if self.wavelengths_nm.size else 0.0

    def extract_wavelength(self, target_nm: float, color: str = "#0072B2") -> Chromatogram:
        """Extract the diode channel closest to `target_nm` as a 2D Chromatogram."""
        if self.wavelengths_nm.size == 0:
            raise ValueError(f"3D PDA spectrum '{self.file_name}' contains no wavelengths.")
        idx = int(np.argmin(np.abs(self.wavelengths_nm - float(target_nm))))
        actual_wl = round(float(self.wavelengths_nm[idx]), 2)
        intensity_slice = self.absorbance_matrix[:, idx].copy()
        channel_label = f"{self.channel} PDA {actual_wl:g}nm (target {target_nm:g}nm)"
        report = validate_chromatogram_arrays(
            file_name=self.file_name,
            retention_time=self.retention_time,
            intensity=intensity_slice,
            channel=channel_label,
            detector_ambiguous=False,
        )
        return Chromatogram(
            file_name=self.file_name,
            retention_time=self.retention_time,
            intensity=intensity_slice,
            sample_name=self.sample_name,
            injection_id=self.injection_id,
            detector="PDA",
            channel=channel_label,
            units=self.units or "AU",
            wavelength_nm=actual_wl,
            bandwidth_nm=None,
            raw_metadata=self.raw_metadata,
            validation_report=report,
            color=color,
            is_from_3d_pda=True,
        )
