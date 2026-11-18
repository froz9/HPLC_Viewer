"""Parser for Waters Empower 3 ASCII Raw Data (.arw) exports.

Supports both 2D single-channel exports (ELSD, extracted UV/PDA wavelengths)
and 3D Photodiode Array (PDA) multi-wavelength spectral matrix exports,
with variable metadata header lengths and CR/CRLF/LF line endings.
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union
import numpy as np
import pandas as pd

from src.data_import.file_validation import (
    ValidationReport,
    validate_chromatogram_arrays,
)
from src.data_processing.chromatogram import (
    Chromatogram,
    PDA3DSpectrum,
    identify_detector_and_wavelength,
)
from src.utils.logging_config import get_logger

logger = get_logger(__name__)


class ARWImportError(ValueError):
    """Raised when an .arw file has an unrecognizable or corrupted format."""


@dataclass
class ARWParseResult:
    """Result returned from parsing an Empower .arw file."""

    file_name: str
    chromatograms: List[Chromatogram] = field(default_factory=list)
    pda_spectrum: Optional[PDA3DSpectrum] = None
    encoding: str = "ascii"
    is_3d_pda: bool = False
    raw_metadata: Dict[str, str] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return len(self.errors) == 0 and len(self.chromatograms) > 0


def detect_and_decode_bytes(raw_bytes: bytes) -> Tuple[str, str]:
    """Decode raw bytes safely, returning (decoded_text, encoding_name)."""
    if not raw_bytes:
        return ("", "ascii")
    if raw_bytes.startswith(b"\xef\xbb\xbf"):
        return (raw_bytes.decode("utf-8-sig"), "utf-8-sig")
    if raw_bytes.startswith(b"\xff\xfe") or raw_bytes.startswith(b"\xfe\xff"):
        return (raw_bytes.decode("utf-16"), "utf-16")
    try:
        return (raw_bytes.decode("ascii"), "ascii")
    except UnicodeDecodeError:
        pass
    try:
        return (raw_bytes.decode("utf-8"), "utf-8")
    except UnicodeDecodeError:
        return (raw_bytes.decode("cp1252", errors="replace"), "cp1252")


def _strip_quotes(token: str) -> str:
    """Strip surrounding whitespace and double quotes from a metadata token."""
    t = token.strip()
    if len(t) >= 2 and t[0] == '"' and t[-1] == '"':
        return t[1:-1].strip()
    return t


def _is_float_token(token: str) -> bool:
    """Return True if `token` can be parsed as a finite float."""
    s = token.strip()
    if not s:
        return False
    try:
        val = float(s)
        return bool(np.isfinite(val))
    except ValueError:
        return False


def _is_2d_numeric_row(tokens: List[str]) -> bool:
    """Return True if tokens represent a 2-column numeric (time, intensity) row."""
    non_empty = [t.strip() for t in tokens if t.strip() != ""]
    if len(non_empty) != 2:
        return False
    return _is_float_token(non_empty[0]) and _is_float_token(non_empty[1])


def _extract_metadata_from_header_lines(header_lines: List[str]) -> Dict[str, str]:
    """Extract metadata key-value pairs from variable-length Empower header rows."""
    metadata: Dict[str, str] = {}
    if not header_lines:
        return metadata

    tokenized = [[_strip_quotes(col) for col in line.split("\t")] for line in header_lines]

    header_keywords = {
        "samplename",
        "sample name",
        "injection",
        "channel",
        "det. units",
        "detector units",
        "vial",
        "date acquired",
    }
    # Case 1: Paired tabular metadata rows (e.g., Row i = headers, Row i+1 = values)
    i = 0
    while i < len(tokenized):
        row = tokenized[i]
        row_lower = [c.lower() for c in row]
        keyword_hits = sum(1 for c in row_lower if c in header_keywords)
        if keyword_hits >= 2 and i + 1 < len(tokenized) and len(tokenized[i + 1]) == len(row):
            vals = tokenized[i + 1]
            for k, v in zip(row, vals):
                if k:
                    metadata[k] = v
            i += 2
            continue
        # Case 2: 2-column Key-Value metadata line (e.g. "SampleName"\t"Diox.MeOH")
        if len(row) == 2 and row[0] and not _is_float_token(row[0]):
            metadata[row[0]] = row[1]
        elif len(row) == 1 and ":" in row[0]:
            k, _, v = row[0].partition(":")
            if k.strip() and not _is_float_token(k.strip()):
                metadata[_strip_quotes(k)] = _strip_quotes(v)
        i += 1

    return metadata


def _normalize_metadata_fields(metadata: Dict[str, str], fallback_name: str) -> Tuple[str, str, str, str]:
    """Return (sample_name, injection_id, channel, units) from parsed metadata dictionary."""
    lower_map = {k.lower().strip(): v for k, v in metadata.items()}

    sample_name = (
        lower_map.get("samplename")
        or lower_map.get("sample name")
        or lower_map.get("sample")
        or Path(fallback_name).stem
    )
    injection_id = (
        lower_map.get("injection")
        or lower_map.get("inj")
        or lower_map.get("injection id")
        or "1"
    )
    channel = (
        lower_map.get("channel")
        or lower_map.get("detector")
        or lower_map.get("signal")
        or ""
    )
    units = (
        lower_map.get("det. units")
        or lower_map.get("detector units")
        or lower_map.get("units")
        or ""
    )
    return (sample_name.strip(), injection_id.strip(), channel.strip(), units.strip())


def read_arw(
    source: Union[str, Path, bytes, io.BytesIO],
    file_name: Optional[str] = None,
    target_pda_wavelengths_nm: Optional[Sequence[float]] = (254.0, 365.0),
    raise_on_error: bool = False,
) -> ARWParseResult:
    """Parse a Waters Empower 3 .arw file (2D single-channel or 3D PDA export).

    Parameters
    ----------
    source : str, Path, bytes, or BytesIO-like object
        File path or in-memory uploaded content.
    file_name : str, optional
        Display file name when `source` is bytes or a stream.
    target_pda_wavelengths_nm : Sequence[float], optional
        Default wavelengths (in nm) to extract as 2D Chromatogram traces when a
        3D PDA matrix file is encountered. Defaults to (254.0, 365.0).
    raise_on_error : bool
        If True, raise `ARWImportError` on invalid/corrupted files instead of
        returning an `ARWParseResult` with `.errors`.
    """
    if isinstance(source, (str, Path)):
        path = Path(source)
        resolved_name = file_name or path.name
        try:
            raw_bytes = path.read_bytes()
        except OSError as exc:
            msg = f"Unable to read file '{resolved_name}': {exc}"
            if raise_on_error:
                raise ARWImportError(msg) from exc
            return ARWParseResult(file_name=resolved_name, errors=[msg])
    elif isinstance(source, bytes):
        resolved_name = file_name or "uploaded.arw"
        raw_bytes = source
    else:
        resolved_name = file_name or getattr(source, "name", "uploaded.arw")
        if hasattr(source, "getvalue"):
            raw_bytes = source.getvalue()
        else:
            source.seek(0)
            raw_bytes = source.read()

    if not raw_bytes or not raw_bytes.strip():
        msg = f"File '{resolved_name}' is empty."
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(file_name=resolved_name, errors=[msg])

    # Check for binary / non-text corruption (null bytes)
    if b"\x00" in raw_bytes[:4096] and not (
        raw_bytes.startswith(b"\xff\xfe") or raw_bytes.startswith(b"\xfe\xff")
    ):
        msg = f"File '{resolved_name}' appears to be a binary or corrupted file (contains NUL bytes)."
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(file_name=resolved_name, errors=[msg])

    # Normalize CR/CRLF line endings at the byte level for fast scanning
    normalized_bytes = raw_bytes.replace(b"\r\n", b"\n").replace(b"\r", b"\n")

    # Inspect the first ~64 KB to decode header lines and determine file type
    head_chunk = normalized_bytes[:65536]
    head_text, encoding = detect_and_decode_bytes(head_chunk)
    head_lines = [line for line in head_text.split("\n") if line.strip() != ""]

    if not head_lines:
        msg = f"File '{resolved_name}' contains no readable lines."
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(file_name=resolved_name, errors=[msg])

    # Scan header lines to detect either a 3D PDA Wavelength header or 2D numeric start
    is_3d_pda = False
    wavelength_line_idx: Optional[int] = None
    data_start_line_idx: Optional[int] = None

    for idx, line in enumerate(head_lines[:50]):
        tokens = line.split("\t")
        first_tok = _strip_quotes(tokens[0]).lower()
        if first_tok == "wavelength" and len(tokens) >= 3:
            # Verify subsequent tokens on this line are numeric wavelengths
            wl_candidates = [t.strip() for t in tokens[1:] if t.strip() != ""]
            if wl_candidates and all(_is_float_token(w) for w in wl_candidates[:5]):
                is_3d_pda = True
                wavelength_line_idx = idx
                # Check if the next non-empty line is "Time"
                if idx + 1 < len(head_lines):
                    next_toks = [t.strip() for t in head_lines[idx + 1].split("\t") if t.strip()]
                    if len(next_toks) == 1 and _strip_quotes(next_toks[0]).lower() == "time":
                        data_start_line_idx = idx + 2
                    else:
                        data_start_line_idx = idx + 1
                break

        if _is_2d_numeric_row(tokens):
            data_start_line_idx = idx
            break

    if data_start_line_idx is None:
        msg = (
            f"Format error in '{resolved_name}': Could not identify a valid chromatographic "
            "data section (expected 2-column retention-time/intensity rows or 3D PDA Wavelength matrix)."
        )
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(file_name=resolved_name, encoding=encoding, errors=[msg])

    # Extract metadata from the header lines preceding the data section (or wavelength line)
    meta_end_idx = wavelength_line_idx if is_3d_pda and wavelength_line_idx is not None else data_start_line_idx
    header_lines = head_lines[:meta_end_idx]
    raw_metadata = _extract_metadata_from_header_lines(header_lines)
    sample_name, injection_id, channel, units = _normalize_metadata_fields(
        raw_metadata, resolved_name
    )

    if is_3d_pda and wavelength_line_idx is not None:
        return _parse_3d_pda_arw(
            normalized_bytes=normalized_bytes,
            resolved_name=resolved_name,
            encoding=encoding,
            raw_metadata=raw_metadata,
            sample_name=sample_name,
            injection_id=injection_id,
            channel=channel,
            units=units,
            wavelength_line_str=head_lines[wavelength_line_idx],
            data_start_line_idx=data_start_line_idx,
            target_wavelengths_nm=target_pda_wavelengths_nm,
            raise_on_error=raise_on_error,
        )

    return _parse_2d_arw(
        normalized_bytes=normalized_bytes,
        resolved_name=resolved_name,
        encoding=encoding,
        raw_metadata=raw_metadata,
        sample_name=sample_name,
        injection_id=injection_id,
        channel=channel,
        units=units,
        data_start_line_idx=data_start_line_idx,
        raise_on_error=raise_on_error,
    )


def _parse_2d_arw(
    normalized_bytes: bytes,
    resolved_name: str,
    encoding: str,
    raw_metadata: Dict[str, str],
    sample_name: str,
    injection_id: str,
    channel: str,
    units: str,
    data_start_line_idx: int,
    raise_on_error: bool,
) -> ARWParseResult:
    """Parse a 2D single-channel .arw file into a Chromatogram."""
    full_text, _ = detect_and_decode_bytes(normalized_bytes)
    non_empty_lines = [line for line in full_text.split("\n") if line.strip() != ""]
    data_lines = non_empty_lines[data_start_line_idx:]

    if not data_lines:
        msg = f"File '{resolved_name}' contains metadata headers but no chromatographic data rows."
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(
            file_name=resolved_name,
            encoding=encoding,
            raw_metadata=raw_metadata,
            errors=[msg],
        )

    times: List[float] = []
    intensities: List[float] = []
    missing_rt = 0
    missing_int = 0
    non_numeric = 0
    wrong_col_count = 0

    for line in data_lines:
        parts = line.split("\t")
        if len(parts) < 2:
            wrong_col_count += 1
            non_numeric += 1
            continue
        # Check if there are extra non-empty columns (ambiguous multi-column file)
        extra = [p.strip() for p in parts[2:] if p.strip() != ""]
        if extra:
            wrong_col_count += 1
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

    if wrong_col_count > 0 and len(times) == 0:
        msg = (
            f"Ambiguous or corrupted column structure in '{resolved_name}': "
            f"{wrong_col_count} row(s) did not contain 2 tab-separated columns."
        )
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(
            file_name=resolved_name,
            encoding=encoding,
            raw_metadata=raw_metadata,
            errors=[msg],
        )

    if len(times) == 0:
        msg = f"No valid numeric retention-time/intensity pairs found in '{resolved_name}'."
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(
            file_name=resolved_name,
            encoding=encoding,
            raw_metadata=raw_metadata,
            errors=[msg],
        )

    detector, wavelength_nm, bandwidth_nm, is_ambiguous = identify_detector_and_wavelength(
        channel=channel, units=units, is_3d_pda=False
    )

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
        raw_metadata=raw_metadata,
        validation_report=report,
        is_from_3d_pda=False,
    )

    logger.info(
        "Imported 2D .arw '%s' (Sample='%s', Detector='%s', Channel='%s', Points=%d)",
        resolved_name,
        sample_name,
        detector,
        channel,
        chrom.n_points,
    )

    return ARWParseResult(
        file_name=resolved_name,
        chromatograms=[chrom],
        pda_spectrum=None,
        encoding=encoding,
        is_3d_pda=False,
        raw_metadata=raw_metadata,
    )


def _parse_3d_pda_arw(
    normalized_bytes: bytes,
    resolved_name: str,
    encoding: str,
    raw_metadata: Dict[str, str],
    sample_name: str,
    injection_id: str,
    channel: str,
    units: str,
    wavelength_line_str: str,
    data_start_line_idx: int,
    target_wavelengths_nm: Optional[Sequence[float]],
    raise_on_error: bool,
) -> ARWParseResult:
    """Parse a 3D PDA multi-wavelength .arw file using fast vectorized C parsing."""
    wl_tokens = [t.strip() for t in wavelength_line_str.split("\t")[1:] if t.strip() != ""]
    try:
        wavelengths_nm = np.array([float(w) for w in wl_tokens], dtype=np.float64)
    except ValueError as exc:
        msg = f"Invalid wavelength header in 3D PDA file '{resolved_name}': {exc}"
        if raise_on_error:
            raise ARWImportError(msg) from exc
        return ARWParseResult(file_name=resolved_name, errors=[msg])

    # Advance byte offset past the `data_start_line_idx` non-empty header lines
    pos = 0
    non_empty_seen = 0
    n_bytes = len(normalized_bytes)
    while pos < n_bytes and non_empty_seen < data_start_line_idx:
        nxt = normalized_bytes.find(b"\n", pos)
        if nxt == -1:
            pos = n_bytes
            break
        line_slice = normalized_bytes[pos:nxt].strip()
        if line_slice:
            non_empty_seen += 1
        pos = nxt + 1

    body_bytes = normalized_bytes[pos:].strip()
    if not body_bytes:
        msg = f"3D PDA file '{resolved_name}' contains headers but no time-series matrix."
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(file_name=resolved_name, errors=[msg])

    try:
        df = pd.read_csv(
            io.BytesIO(body_bytes),
            sep="\t",
            header=None,
            engine="c",
            dtype=np.float64,
        )
    except Exception as exc:
        msg = f"Failed to parse 3D PDA numeric matrix in '{resolved_name}': {exc}"
        if raise_on_error:
            raise ARWImportError(msg) from exc
        return ARWParseResult(file_name=resolved_name, errors=[msg])

    expected_cols = 1 + wavelengths_nm.size
    if df.shape[1] != expected_cols:
        msg = (
            f"Column count mismatch in 3D PDA file '{resolved_name}': "
            f"expected {expected_cols} columns (1 Time + {wavelengths_nm.size} Wavelengths), "
            f"found {df.shape[1]}."
        )
        if raise_on_error:
            raise ARWImportError(msg)
        return ARWParseResult(file_name=resolved_name, errors=[msg])

    mat = df.to_numpy(dtype=np.float64, copy=False)
    rt_arr = mat[:, 0].copy()
    abs_mat = mat[:, 1:].copy()

    pda_spectrum = PDA3DSpectrum(
        file_name=resolved_name,
        sample_name=sample_name,
        injection_id=injection_id,
        channel=channel or "2998",
        units=units or "AU",
        raw_metadata=raw_metadata,
        retention_time=rt_arr,
        wavelengths_nm=wavelengths_nm,
        absorbance_matrix=abs_mat,
    )

    extracted_chroms: List[Chromatogram] = []
    targets = list(target_wavelengths_nm) if target_wavelengths_nm else [254.0]
    seen_indices = set()
    for target_wl in targets:
        idx = int(np.argmin(np.abs(wavelengths_nm - float(target_wl))))
        if idx not in seen_indices:
            seen_indices.add(idx)
            extracted_chroms.append(pda_spectrum.extract_wavelength(float(target_wl)))

    logger.info(
        "Imported 3D PDA .arw '%s' (Sample='%s', %d time points x %d wavelengths [%.1f-%.1f nm])",
        resolved_name,
        sample_name,
        pda_spectrum.n_times,
        pda_spectrum.n_wavelengths,
        pda_spectrum.wl_min,
        pda_spectrum.wl_max,
    )

    return ARWParseResult(
        file_name=resolved_name,
        chromatograms=extracted_chroms,
        pda_spectrum=pda_spectrum,
        encoding=encoding,
        is_3d_pda=True,
        raw_metadata=raw_metadata,
    )
