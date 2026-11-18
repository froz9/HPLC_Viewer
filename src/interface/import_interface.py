"""File import caching, batch import progress, and interactive import summary table."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import pandas as pd
import streamlit as st

from src.data_import.arw_reader import ARWParseResult, read_arw
from src.data_import.csv_reader import DelimitedParseResult, read_csv_chromatogram
from src.data_import.file_validation import VALID_DETECTORS
from src.data_import.txt_reader import read_txt_chromatogram
from src.data_processing.chromatogram import Chromatogram, PDA3DSpectrum
from src.visualization.color_manager import assign_colors_to_chromatograms


@st.cache_data(show_spinner=False, max_entries=128)
def cached_parse_file_bytes(
    file_name: str,
    file_hash: str,
    raw_bytes: bytes,
    target_pda_wls: Tuple[float, ...] = (254.0, 365.0),
) -> UnionResult:
    """Parse a chromatographic file once and cache the result by content digest."""
    _ = file_hash  # Used as cache key discriminator
    ext = Path(file_name).suffix.lower()
    if ext == ".arw":
        return read_arw(
            raw_bytes,
            file_name=file_name,
            target_pda_wavelengths_nm=target_pda_wls,
            raise_on_error=False,
        )
    if ext == ".csv":
        return read_csv_chromatogram(raw_bytes, file_name=file_name, raise_on_error=False)
    if ext == ".txt":
        return read_txt_chromatogram(raw_bytes, file_name=file_name, raise_on_error=False)

    # Fallback: sniff content for .arw vs delimited
    if b"SampleName" in raw_bytes[:512] and b"\t" in raw_bytes[:512]:
        return read_arw(
            raw_bytes,
            file_name=file_name,
            target_pda_wavelengths_nm=target_pda_wls,
            raise_on_error=False,
        )
    return read_csv_chromatogram(raw_bytes, file_name=file_name, raise_on_error=False)


UnionResult = ARWParseResult | DelimitedParseResult


def _disambiguate_sample_display_name(
    new_chrom: Chromatogram,
    existing_chroms: Dict[str, Chromatogram],
) -> None:
    """Append file number suffix if a different file injection shares the same SampleName + Channel."""
    for existing in existing_chroms.values():
        if (
            existing.file_name != new_chrom.file_name
            and existing.original_sample_name == new_chrom.original_sample_name
            and existing.channel == new_chrom.channel
        ):
            m = re.search(r"(\d+)(?=\.[^.]+$)", new_chrom.file_name)
            suffix = m.group(1) if m else Path(new_chrom.file_name).stem
            new_chrom.sample_display_name = f"{new_chrom.original_sample_name} (#{suffix})"
            return


def ingest_files_into_session(
    items: Sequence[Tuple[str, bytes]],
    target_pda_wls: Tuple[float, ...] = (254.0, 365.0),
) -> Tuple[int, int]:
    """Parse a batch of `(file_name, raw_bytes)` items and add them to `st.session_state`."""
    if not items:
        return (0, 0)

    chrom_map: Dict[str, Chromatogram] = st.session_state.setdefault("chromatograms", {})
    pda_map: Dict[str, PDA3DSpectrum] = st.session_state.setdefault("pda_spectra", {})
    import_errors: Dict[str, List[str]] = st.session_state.setdefault("import_errors", {})
    processed_keys: set = st.session_state.setdefault("processed_file_keys", set())

    added_count = 0
    error_count = 0
    total = len(items)
    progress = st.progress(0.0, text="Importing chromatographic files...") if total > 1 else None

    for idx, (fname, raw_bytes) in enumerate(items, start=1):
        digest = hashlib.sha1(raw_bytes[:65536] + str(len(raw_bytes)).encode()).hexdigest()[:16]
        session_key = f"{fname}::{digest}"
        if session_key in processed_keys:
            if progress is not None:
                progress.progress(idx / total, text=f"Skipped already imported {fname}")
            continue

        res = cached_parse_file_bytes(
            file_name=fname,
            file_hash=digest,
            raw_bytes=raw_bytes,
            target_pda_wls=target_pda_wls,
        )
        processed_keys.add(session_key)

        if not res.is_valid:
            import_errors[fname] = list(res.errors)
            error_count += 1
        else:
            import_errors.pop(fname, None)
            if isinstance(res, ARWParseResult) and res.pda_spectrum is not None:
                pda_map[fname] = res.pda_spectrum
            for chrom in res.chromatograms:
                # Clone so session edits do not mutate cached object
                fresh_chrom = Chromatogram(
                    file_name=chrom.file_name,
                    retention_time=chrom.raw_retention_time,
                    intensity=chrom.raw_intensity,
                    sample_name=chrom.original_sample_name,
                    injection_id=chrom.injection_id,
                    detector=chrom.original_detector,
                    channel=chrom.channel,
                    units=chrom.units,
                    wavelength_nm=chrom.wavelength_nm,
                    bandwidth_nm=chrom.bandwidth_nm,
                    raw_metadata=chrom.raw_metadata,
                    validation_report=chrom.validation_report,
                    trace_id=chrom.trace_id,
                    is_from_3d_pda=chrom.is_from_3d_pda,
                )
                if fresh_chrom.trace_id not in chrom_map:
                    _disambiguate_sample_display_name(fresh_chrom, chrom_map)
                    chrom_map[fresh_chrom.trace_id] = fresh_chrom
                    added_count += 1

        if progress is not None:
            progress.progress(idx / total, text=f"Processed {fname} ({idx}/{total})")

    if progress is not None:
        progress.empty()

    if added_count > 0:
        palette = st.session_state.get("palette_name", "Okabe-Ito (Colorblind Safe)")
        group_by_sample = st.session_state.get("group_colors_by_sample", True)
        assign_colors_to_chromatograms(
            list(chrom_map.values()),
            palette_name=palette,
            group_by_sample_name=group_by_sample,
        )

    return (added_count, error_count)


def render_import_summary_and_table(chromatograms: List[Chromatogram]) -> None:
    """Render the KPI summary bar, editable sample/channel table, and validation report."""
    import_errors: Dict[str, List[str]] = st.session_state.get("import_errors", {})
    if import_errors:
        for fname, errs in import_errors.items():
            st.error(f"**Import Error — `{fname}`**: {'; '.join(errs)}")

    if not chromatograms:
        return

    unique_samples = {c.sample_display_name for c in chromatograms}
    detectors_present = sorted({c.assigned_detector for c in chromatograms})
    aligned_count = sum(1 for c in chromatograms if c.is_aligned)
    warning_count = sum(1 for c in chromatograms if c.validation_report.status != "OK")

    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("Chromatogram Channels", len(chromatograms))
    col2.metric("Unique Samples", len(unique_samples))
    col3.metric("Detectors", ", ".join(detectors_present))
    col4.metric("Aligned Traces", f"{aligned_count} / {len(chromatograms)}")
    col5.metric(
        "Validation Status",
        "All OK" if warning_count == 0 else f"{warning_count} Warning(s)",
    )

    with st.expander("📋 Imported Chromatograms & Sample Assignment Table (Editable)", expanded=True):
        st.caption(
            "Edit **Visible**, **Sample Name**, **Detector**, or **Color** directly in the table below. "
            "Original experimental files and raw metadata are never modified."
        )
        rows = [c.to_summary_dict() for c in chromatograms]
        df = pd.DataFrame(rows)

        display_cols = [
            "Visible",
            "Color",
            "File Name",
            "Sample Name",
            "Injection",
            "Detector",
            "Channel",
            "Wavelength (nm)",
            "Units",
            "Data Points",
            "RT Range",
            "Sampling Interval",
            "Time Shift (min)",
            "Status",
        ]

        edited_df = st.data_editor(
            df[display_cols],
            key="main_import_table_editor",
            width="stretch",
            hide_index=True,
            column_config={
                "Visible": st.column_config.CheckboxColumn(
                    "Plot",
                    help="Show or hide this chromatogram channel in the visualization panel.",
                    default=True,
                ),
                "Color": st.column_config.TextColumn(
                    "Color (Hex)",
                    help="Hex color code (e.g. #0072B2) for this chromatogram.",
                ),
                "File Name": st.column_config.TextColumn("File Name", disabled=True),
                "Sample Name": st.column_config.TextColumn(
                    "Sample Name",
                    help="Editable display name for legends and reports (does not alter source file).",
                ),
                "Injection": st.column_config.TextColumn("Inj.", disabled=True, width="small"),
                "Detector": st.column_config.SelectboxColumn(
                    "Detector",
                    options=list(VALID_DETECTORS),
                    required=True,
                    help="Assigned detector type (UV, PDA, ELSD, or Other / Unknown).",
                ),
                "Channel": st.column_config.TextColumn("Channel", disabled=True),
                "Wavelength (nm)": st.column_config.NumberColumn(
                    "Wavelength (nm)", format="%.1f", disabled=True
                ),
                "Units": st.column_config.TextColumn("Units", disabled=True, width="small"),
                "Data Points": st.column_config.NumberColumn("Points", disabled=True),
                "RT Range": st.column_config.TextColumn("RT Range", disabled=True),
                "Sampling Interval": st.column_config.TextColumn("Sampling Interval", disabled=True),
                "Time Shift (min)": st.column_config.NumberColumn(
                    "Shift (min)", format="%+.4f", disabled=True
                ),
                "Status": st.column_config.TextColumn("Import Status", disabled=True),
            },
        )

        # Apply user edits back to Chromatogram instances
        for idx, chrom in enumerate(chromatograms):
            row = edited_df.iloc[idx]
            chrom.visible = bool(row["Visible"])
            new_name = str(row["Sample Name"]).strip()
            if new_name:
                chrom.sample_display_name = new_name
            new_det = str(row["Detector"])
            if new_det in VALID_DETECTORS:
                chrom.assigned_detector = new_det
            new_color = str(row["Color"]).strip()
            if new_color.startswith("#") and len(new_color) in (4, 7):
                chrom.color = new_color

    with st.expander("🔍 Data Validation & Quality Report", expanded=(warning_count > 0)):
        val_rows = []
        for c in chromatograms:
            r = c.validation_report
            val_rows.append(
                {
                    "File Name": c.file_name,
                    "Sample": c.sample_display_name,
                    "Channel": c.channel,
                    "Status": r.status,
                    "Missing RT": r.missing_rt_count,
                    "Missing Int.": r.missing_intensity_count,
                    "Non-Numeric": r.non_numeric_count,
                    "Duplicate RT": r.duplicate_time_count,
                    "Unsorted RT": r.unsorted_time_count,
                    "Negative Baseline Pts": r.negative_intensity_count,
                    "Unequal Δt": "Yes" if r.unequal_sampling_interval else "No",
                    "Notes": r.summary_text,
                }
            )
        st.dataframe(pd.DataFrame(val_rows), width="stretch", hide_index=True)
