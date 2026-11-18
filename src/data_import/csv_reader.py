"""Parser for CSV and general delimited chromatographic files (.csv)."""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd

from src.data_import.arw_reader import (
    _extract_metadata_from_header_lines,
    _is_float_token,
    _normalize_metadata_fields,
    _strip_quotes,
    detect_and_decode_bytes,
)
from src.data_import.file_validation import validate_chromatogram_arrays
from src.data_processing.chromatogram import (
    Chromatogram,
    identify_detector_and_wavelength,
)
from src.utils.logging_config import get_logger

logger = get_logger(__name__)


class DelimitedImportError(ValueError):
    """Raised when a delimited text/CSV file cannot be parsed unambiguously."""


@dataclass
class DelimitedParseResult:
    """Result returned from parsing a CSV or TXT chromatographic file."""

    file_name: str
    chromatograms: List[Chromatogram] = field(default_factory=list)
    delimiter: str = ","
    encoding: str = "utf-8"
    raw_metadata: Dict[str, str] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return len(self.errors) == 0 and len(self.chromatograms) > 0


def _detect_delimiter(lines: List[str], default_delim: str = ",") -> str:
    """Detect the consistent delimiter across candidate data lines."""
    candidates = ["\t", ",", ";"]
    # Prefer the last 10 non-empty lines (which are in the numeric data section)
    sample_lines = [ln for ln in lines[-20:] if ln.strip() and not ln.strip().startswith("#")]
    if not sample_lines:
        sample_lines = [ln for ln in lines[:20] if ln.strip()]

    best_delim = default_delim
    best_score = -1
    for delim in candidates:
        counts = [ln.count(delim) for ln in sample_lines]
        if counts and min(counts) >= 1 and len(set(counts)) == 1:
            score = counts[0] * 10
            if score > best_score:
                best_score = score
                best_delim = delim
        elif counts and sum(c >= 1 for c in counts) >= max(1, int(0.8 * len(counts))):
            med_c = int(np.median(counts))
            if med_c > best_score:
                best_score = med_c
                best_delim = delim
    return best_delim


def parse_delimited_chromatogram(
    source: Union[str, Path, bytes, io.BytesIO],
    file_name: Optional[str] = None,
    preferred_delimiter: str = ",",
    raise_on_error: bool = False,
) -> DelimitedParseResult:
    """Parse a delimited chromatographic file (.csv or .txt)."""
    if isinstance(source, (str, Path)):
        path = Path(source)
        resolved_name = file_name or path.name
        try:
            raw_bytes = path.read_bytes()
        except OSError as exc:
            msg = f"Unable to read file '{resolved_name}': {exc}"
            if raise_on_error:
                raise DelimitedImportError(msg) from exc
            return DelimitedParseResult(file_name=resolved_name, errors=[msg])
    elif isinstance(source, bytes):
        resolved_name = file_name or "uploaded.csv"
        raw_bytes = source
    else:
        resolved_name = file_name or getattr(source, "name", "uploaded.csv")
        if hasattr(source, "getvalue"):
            raw_bytes = source.getvalue()
        else:
            source.seek(0)
            raw_bytes = source.read()

    if not raw_bytes or not raw_bytes.strip():
        msg = f"File '{resolved_name}' is empty."
        if raise_on_error:
            raise DelimitedImportError(msg)
        return DelimitedParseResult(file_name=resolved_name, errors=[msg])

    if b"\x00" in raw_bytes[:4096] and not (
        raw_bytes.startswith(b"\xff\xfe") or raw_bytes.startswith(b"\xfe\xff")
    ):
        msg = f"File '{resolved_name}' appears to be binary or corrupted (contains NUL bytes)."
        if raise_on_error:
            raise DelimitedImportError(msg)
        return DelimitedParseResult(file_name=resolved_name, errors=[msg])

    normalized_bytes = raw_bytes.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    text, encoding = detect_and_decode_bytes(normalized_bytes)
    all_lines = [ln for ln in text.split("\n") if ln.strip() != ""]

    if not all_lines:
        msg = f"File '{resolved_name}' contains no non-empty lines."
        if raise_on_error:
            raise DelimitedImportError(msg)
        return DelimitedParseResult(file_name=resolved_name, encoding=encoding, errors=[msg])

    delim = _detect_delimiter(all_lines, default_delim=preferred_delimiter)

    # Check if this is an HPLC Viewer exported CSV (with `original_retention_time_min` & `original_intensity`)
    first_NonComment = next((ln for ln in all_lines if not ln.strip().startswith("#")), "")
    first_cols = [_strip_quotes(c).lower() for c in first_NonComment.split(delim)]
    if "original_retention_time_min" in first_cols and "original_intensity" in first_cols:
        return _parse_exported_hplc_csv(
            text=text,
            resolved_name=resolved_name,
            delim=delim,
            encoding=encoding,
            raise_on_error=raise_on_error,
        )

    # Scan for the first numeric data row
    data_start_idx: Optional[int] = None
    for idx, line in enumerate(all_lines):
        if line.strip().startswith("#"):
            continue
        tokens = [_strip_quotes(t) for t in line.split(delim)]
        non_empty = [t for t in tokens if t != ""]
        if len(non_empty) >= 2 and _is_float_token(non_empty[0]) and _is_float_token(non_empty[1]):
            data_start_idx = idx
            break

    if data_start_idx is None:
        msg = (
            f"Format error in '{resolved_name}': Could not locate numeric retention-time "
            "and intensity columns."
        )
        if raise_on_error:
            raise DelimitedImportError(msg)
        return DelimitedParseResult(
            file_name=resolved_name, delimiter=delim, encoding=encoding, errors=[msg]
        )

    # Parse metadata and potential column header row above `data_start_idx`
    header_lines = all_lines[:data_start_idx]
    metadata: Dict[str, str] = {}
    column_headers: Optional[List[str]] = None

    for idx, hline in enumerate(header_lines):
        s = hline.strip()
        if s.startswith("#"):
            s = s.lstrip("#").strip()
        if ":" in s and delim not in s:
            k, _, v = s.partition(":")
            metadata[_strip_quotes(k)] = _strip_quotes(v)
            continue
        toks = [_strip_quotes(t) for t in s.split(delim)]
        # Check if this is the column header row immediately preceding data_start_idx
        if idx == len(header_lines) - 1 and len(toks) >= 2:
            low0 = toks[0].lower()
            if any(w in low0 for w in ("time", "rt", "retention", "min")):
                column_headers = toks
                continue
        if len(toks) == 2 and toks[0] and not _is_float_token(toks[0]):
            metadata[toks[0]] = toks[1]

    # Also support Empower-style 2-row header if delimiter is tab or comma
    if not metadata and len(header_lines) >= 2:
        tab_Normalized = [delim.join(_strip_quotes(c) for c in ln.split(delim)).replace(delim, "\t") for ln in header_lines]
        metadata.update(_extract_metadata_from_header_lines(tab_Normalized))

    sample_name, injection_id, channel, units = _normalize_metadata_fields(
        metadata, resolved_name
    )

    # Parse data lines
    data_lines = all_lines[data_start_idx:]
    times: List[float] = []
    intensities: List[float] = []
    missing_rt = 0
    missing_int = 0
    non_numeric = 0

    for line in data_lines:
        if line.strip().startswith("#"):
            continue
        parts = [_strip_quotes(p) for p in line.split(delim)]
        if len(parts) < 2:
            non_numeric += 1
            continue
        t_str = parts[0].strip()
        i_str = parts[1].strip()
        if t_str == "" and i_str == "":
            continue
        if t_str == "":
            missing_rt += 1
            continue
        if i_str == "":
            missing_int += 1
            continue
        try:
            t_val = float(t_str)
            i_val = float(i_str)
            if not np.isfinite(t_val):
                missing_rt += 1
                continue
            if not np.isfinite(i_val):
                missing_int += 1
                continue
            times.append(t_val)
            intensities.append(i_val)
        except ValueError:
            non_numeric += 1

    if not times:
        msg = f"No valid numeric observations found in '{resolved_name}'."
        if raise_on_error:
            raise DelimitedImportError(msg)
        return DelimitedParseResult(
            file_name=resolved_name,
            delimiter=delim,
            encoding=encoding,
            raw_metadata=metadata,
            errors=[msg],
        )

    if not channel and column_headers and len(column_headers) >= 2:
        col1_hdr = column_headers[1]
        if col1_hdr.lower() not in ("intensity", "signal", "response", "value", "y"):
            channel = col1_hdr

    detector_meta = metadata.get("Detector") or metadata.get("detector") or ""
    if detector_meta.upper() in ("UV", "PDA", "ELSD"):
        _, wavelength_nm, bandwidth_nm, _ = identify_detector_and_wavelength(
            channel=channel or detector_meta, units=units
        )
        detector = detector_meta.upper()
        is_ambiguous = False
    else:
        detector, wavelength_nm, bandwidth_nm, is_ambiguous = identify_detector_and_wavelength(
            channel=channel, units=units
        )

    # Check explicit Wavelength metadata field if present
    wl_meta = metadata.get("Wavelength") or metadata.get("wavelength") or metadata.get("wavelength_nm")
    if wl_meta and wavelength_nm is None:
        try:
            wavelength_nm = float(str(wl_meta).lower().replace("nm", "").strip())
            if detector == "Other / Unknown":
                detector = "UV"
                is_ambiguous = False
        except ValueError:
            pass

    report = validate_chromatogram_arrays(
        file_name=resolved_name,
        retention_time=times,
        intensity=intensities,
        channel=channel,
        missing_rt_count=missing_rt,
        missing_intensity_count=missing_int,
        non_numeric_count=non_numeric,
        detector_ambiguous=is_ambiguous,
    )

    chrom = Chromatogram(
        file_name=resolved_name,
        retention_time=times,
        intensity=intensities,
        sample_name=sample_name,
        injection_id=injection_id,
        detector=detector,
        channel=channel,
        units=units,
        wavelength_nm=wavelength_nm,
        bandwidth_nm=bandwidth_nm,
        raw_metadata=metadata,
        validation_report=report,
    )

    return DelimitedParseResult(
        file_name=resolved_name,
        chromatograms=[chrom],
        delimiter=delim,
        encoding=encoding,
        raw_metadata=metadata,
    )


def _parse_exported_hplc_csv(
    text: str,
    resolved_name: str,
    delim: str,
    encoding: str,
    raise_on_error: bool,
) -> DelimitedParseResult:
    """Parse a CSV file previously exported by HPLC Viewer."""
    try:
        df = pd.read_csv(io.StringIO(text), sep=delim, comment="#")
    except Exception as exc:
        msg = f"Failed to parse exported CSV '{resolved_name}': {exc}"
        if raise_on_error:
            raise DelimitedImportError(msg) from exc
        return DelimitedParseResult(file_name=resolved_name, errors=[msg])

    group_cols = [c for c in ("sample_identifier", "channel", "detector") if c in df.columns]
    if not group_cols:
        group_iter = [(("Sample",), df)]
    else:
        group_iter = list(df.groupby(group_cols, dropna=False, sort=False))

    chroms: List[Chromatogram] = []
    for _, sub_df in group_iter:
        sample_name = (
            str(sub_df["sample_identifier"].iloc[0])
            if "sample_identifier" in sub_df.columns
            else Path(resolved_name).stem
        )
        orig_sample = (
            str(sub_df["original_sample_name"].iloc[0])
            if "original_sample_name" in sub_df.columns
            else sample_name
        )
        inj = (
            str(sub_df["injection_id"].iloc[0])
            if "injection_id" in sub_df.columns
            else "1"
        )
        det = (
            str(sub_df["detector"].iloc[0])
            if "detector" in sub_df.columns
            else "Other / Unknown"
        )
        ch = (
            str(sub_df["channel"].iloc[0])
            if "channel" in sub_df.columns and pd.notna(sub_df["channel"].iloc[0])
            else ""
        )
        un = (
            str(sub_df["units"].iloc[0])
            if "units" in sub_df.columns and pd.notna(sub_df["units"].iloc[0])
            else ""
        )
        wl: Optional[float] = None
        if "wavelength_nm" in sub_df.columns and pd.notna(sub_df["wavelength_nm"].iloc[0]):
            try:
                wl = float(sub_df["wavelength_nm"].iloc[0])
            except ValueError:
                wl = None

        rt = sub_df["original_retention_time_min"].to_numpy(dtype=np.float64)
        inten = sub_df["original_intensity"].to_numpy(dtype=np.float64)

        chrom = Chromatogram(
            file_name=resolved_name,
            retention_time=rt,
            intensity=inten,
            sample_name=orig_sample,
            injection_id=inj,
            detector=det,
            channel=ch,
            units=un,
            wavelength_nm=wl,
        )
        chrom.sample_display_name = sample_name
        if "applied_shift_min" in sub_df.columns:
            shift_val = float(sub_df["applied_shift_min"].iloc[0])
            method_val = (
                str(sub_df["alignment_method"].iloc[0])
                if "alignment_method" in sub_df.columns
                else "Imported"
            )
            if abs(shift_val) > 0.0 and method_val != "None":
                chrom.apply_shift(shift_val, method=method_val)
        chroms.append(chrom)

    return DelimitedParseResult(
        file_name=resolved_name,
        chromatograms=chroms,
        delimiter=delim,
        encoding=encoding,
    )


def read_csv_chromatogram(
    source: Union[str, Path, bytes, io.BytesIO],
    file_name: Optional[str] = None,
    raise_on_error: bool = False,
) -> DelimitedParseResult:
    """Read a CSV chromatographic file (.csv)."""
    return parse_delimited_chromatogram(
        source=source,
        file_name=file_name,
        preferred_delimiter=",",
        raise_on_error=raise_on_error,
    )
