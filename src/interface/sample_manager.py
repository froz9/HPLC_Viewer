"""Sidebar Sample Manager for selection, renaming, detector assignment, color control, and 3D PDA extraction."""

from __future__ import annotations

from typing import Dict, List
import streamlit as st

from src.data_import.file_validation import VALID_DETECTORS
from src.data_processing.chromatogram import Chromatogram, PDA3DSpectrum
from src.visualization.color_manager import (
    COLOR_PALETTES,
    DEFAULT_PALETTE_NAME,
    reset_colors_to_palette,
)


def render_sidebar_sample_manager(
    chrom_map: Dict[str, Chromatogram],
    pda_map: Dict[str, PDA3DSpectrum],
) -> None:
    """Render Section B (Samples) inside the Streamlit sidebar."""
    chromatograms = list(chrom_map.values())
    if not chromatograms:
        st.caption("Import chromatographic files above to manage samples.")
        return

    # 1. Quick bulk selection buttons
    col_all, col_none = st.columns(2)
    if col_all.button("Select All", width="stretch", help="Display all imported chromatograms."):
        for c in chromatograms:
            c.visible = True
        st.rerun()
    if col_none.button("Deselect All", width="stretch", help="Hide all imported chromatograms."):
        for c in chromatograms:
            c.visible = False
        st.rerun()

    # 2. Color Palette & Reset
    palette_names = list(COLOR_PALETTES.keys())
    current_palette = st.session_state.get("palette_name", DEFAULT_PALETTE_NAME)
    selected_palette = st.selectbox(
        "Color Palette",
        options=palette_names,
        index=palette_names.index(current_palette) if current_palette in palette_names else 0,
        help="Select a scientific color palette for automatic trace color assignment.",
    )
    st.session_state["palette_name"] = selected_palette

    group_by_sample = st.checkbox(
        "Match colors by sample name",
        value=st.session_state.get("group_colors_by_sample", True),
        help="When enabled, UV and ELSD channels from the same sample share the same color.",
    )
    st.session_state["group_colors_by_sample"] = group_by_sample

    if st.button("Reset Colors to Palette", width="stretch"):
        reset_colors_to_palette(
            chromatograms,
            palette_name=selected_palette,
            group_by_sample_name=group_by_sample,
        )
        st.rerun()

    # 3. Individual Sample / Trace Editor & Removal
    with st.expander("✏️ Edit or Remove Individual Trace", expanded=False):
        trace_options = {c.trace_id: f"{c.display_label} ({c.file_name})" for c in chromatograms}
        selected_tid = st.selectbox(
            "Select Chromatogram",
            options=list(trace_options.keys()),
            format_func=lambda tid: trace_options[tid],
            key="sidebar_edit_trace_selector",
        )
        if selected_tid and selected_tid in chrom_map:
            target = chrom_map[selected_tid]
            new_vis = st.checkbox(
                "Visible in Plot",
                value=target.visible,
                key=f"sb_vis_{selected_tid}",
            )
            target.visible = new_vis

            new_name = st.text_input(
                "Sample Display Name",
                value=target.sample_display_name,
                key=f"sb_name_{selected_tid}",
                help="Rename sample for plots and exports without modifying original data.",
            )
            if new_name.strip():
                target.sample_display_name = new_name.strip()

            det_idx = (
                list(VALID_DETECTORS).index(target.assigned_detector)
                if target.assigned_detector in VALID_DETECTORS
                else 0
            )
            new_det = st.selectbox(
                "Assigned Detector",
                options=list(VALID_DETECTORS),
                index=det_idx,
                key=f"sb_det_{selected_tid}",
            )
            target.assigned_detector = new_det

            new_color = st.color_picker(
                "Trace Color",
                value=target.color if target.color.startswith("#") and len(target.color) == 7 else "#0072B2",
                key=f"sb_col_{selected_tid}",
            )
            target.color = new_color

            if st.button(
                "🗑️ Remove This Chromatogram",
                key=f"sb_rm_{selected_tid}",
                width="stretch",
            ):
                chrom_map.pop(selected_tid, None)
                st.rerun()

    # 4. 3D PDA Wavelength Extractor (when 3D PDA files are loaded)
    if pda_map:
        with st.expander("🌈 Extract 3D PDA Wavelength", expanded=False):
            pda_file = st.selectbox(
                "3D PDA File",
                options=list(pda_map.keys()),
                key="sb_pda_file_select",
            )
            spec = pda_map[pda_file]
            target_wl = st.number_input(
                f"Target Wavelength ({spec.wl_min:.1f}–{spec.wl_max:.1f} nm)",
                min_value=float(spec.wl_min),
                max_value=float(spec.wl_max),
                value=min(max(280.0, float(spec.wl_min)), float(spec.wl_max)),
                step=1.0,
                key="sb_pda_target_wl",
            )
            if st.button("Extract Wavelength Channel", width="stretch"):
                extracted = spec.extract_wavelength(float(target_wl))
                chrom_map[extracted.trace_id] = extracted
                reset_colors_to_palette(
                    list(chrom_map.values()),
                    palette_name=selected_palette,
                    group_by_sample_name=group_by_sample,
                )
                st.success(f"Extracted {extracted.wavelength_nm:g} nm from {pda_file}")
                st.rerun()
