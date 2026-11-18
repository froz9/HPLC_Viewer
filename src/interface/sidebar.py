"""Sidebar configuration panel for HPLC Viewer (Sections A, B, C, and D)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import streamlit as st

from src.data_import.file_validation import VALID_DETECTORS
from src.data_processing.alignment import (
    ALIGNMENT_MODES,
    apply_alignment_estimates,
    apply_manual_shifts,
    estimate_batch_alignment,
    reset_all_alignments,
)
from src.data_processing.chromatogram import Chromatogram, PDA3DSpectrum
from src.data_processing.normalization import NORMALIZATION_MODES
from src.interface.import_interface import ingest_files_into_session
from src.interface.sample_manager import render_sidebar_sample_manager
from src.utils.session_manager import restore_session_config
from src.visualization.chromatogram_plot import VISUALIZATION_MODES, PlotConfig


RAW_DATA_DIR = Path(__file__).resolve().parents[2] / "raw_data"


def render_sidebar() -> Tuple[PlotConfig, Dict[str, Any]]:
    """Render the complete scientific configuration sidebar and return active PlotConfig & alignment state."""
    chrom_map: Dict[str, Chromatogram] = st.session_state.setdefault("chromatograms", {})
    pda_map: Dict[str, PDA3DSpectrum] = st.session_state.setdefault("pda_spectra", {})

    with st.sidebar:
        st.title("🔬 HPLC Viewer")
        st.caption("Analytical Chromatogram Visualization & Alignment")

        # =====================================================================
        # SECTION A: IMPORT
        # =====================================================================
        st.header("A. Import")
        uploaded_files = st.file_uploader(
            "Drag and drop chromatographic files",
            type=["arw", "csv", "txt"],
            accept_multiple_files=True,
            help="Supported formats: Waters Empower 3 ASCII (.arw), CSV (.csv), and Tab-delimited (.txt).",
            key=f"uploader_{st.session_state.get('uploader_epoch', 0)}",
        )

        if uploaded_files:
            batch = [(uf.name, uf.getvalue()) for uf in uploaded_files]
            ingest_files_into_session(batch)

        # Convenient local workspace loader for `raw_data/`
        if RAW_DATA_DIR.exists():
            local_arws = sorted(RAW_DATA_DIR.glob("*.arw"))
            if local_arws:
                with st.expander("📂 Load from Local `raw_data/` Folder", expanded=(len(chrom_map) == 0)):
                    non_3d_files = [f for f in local_arws if f.stat().st_size < 5_000_000]
                    if st.button(
                        "Load Example Set (3 Samples: UV + ELSD)",
                        width="stretch",
                        help="Load 2D UV (254/365 nm) and ELSD files for 01006-18-1, 01006-18-2, and 01006-18-3.",
                    ):
                        sample_names = [
                            "HPLC_RAW_EXPORT1276.arw",
                            "HPLC_RAW_EXPORT1278.arw",
                            "HPLC_RAW_EXPORT1279.arw",
                            "HPLC_RAW_EXPORT1282.arw",
                            "HPLC_RAW_EXPORT1284.arw",
                            "HPLC_RAW_EXPORT1285.arw",
                            "HPLC_RAW_EXPORT1288.arw",
                            "HPLC_RAW_EXPORT1290.arw",
                            "HPLC_RAW_EXPORT1291.arw",
                        ]
                        items = [
                            (fname, (RAW_DATA_DIR / fname).read_bytes())
                            for fname in sample_names
                            if (RAW_DATA_DIR / fname).exists()
                        ]
                        ingest_files_into_session(items)
                        st.rerun()

                    selected_local = st.multiselect(
                        "Select specific `raw_data/` files",
                        options=[f.name for f in local_arws],
                        default=[],
                        help="Choose any 2D or 3D (.arw) files from raw_data/ to import into the session.",
                    )
                    if selected_local and st.button("Import Selected `raw_data/` Files", width="stretch"):
                        items = [
                            (fname, (RAW_DATA_DIR / fname).read_bytes())
                            for fname in selected_local
                        ]
                        ingest_files_into_session(items)
                        st.rerun()

        # Reload saved session configuration (.json)
        with st.expander("⚙️ Reload Session Config (.json)", expanded=False):
            cfg_file = st.file_uploader(
                "Upload saved session JSON",
                type=["json"],
                key="session_json_uploader",
                help="Restore visualization, color, detector, and alignment settings.",
            )
            if cfg_file is not None and st.button("Apply Session Configuration", width="stretch"):
                try:
                    content = cfg_file.getvalue().decode("utf-8")
                    plot_cfg, align_cfg, exp_cfg, pal_name, matched = restore_session_config(
                        content, list(chrom_map.values())
                    )
                    st.session_state["restored_plot_config"] = plot_cfg
                    st.session_state["alignment_settings"] = align_cfg
                    st.session_state["export_settings"] = exp_cfg
                    st.session_state["palette_name"] = pal_name
                    st.success(f"Restored configuration ({matched} sample trace(s) matched).")
                    st.rerun()
                except Exception as exc:
                    st.error(f"Could not restore session JSON: {exc}")

        if st.button(
            "🧹 Clear Current Session",
            width="stretch",
            help="Remove all imported chromatograms, alignment previews, and session data.",
        ):
            st.session_state["chromatograms"] = {}
            st.session_state["pda_spectra"] = {}
            st.session_state["import_errors"] = {}
            st.session_state["processed_file_keys"] = set()
            st.session_state["proposed_alignment"] = None
            st.session_state["uploader_epoch"] = st.session_state.get("uploader_epoch", 0) + 1
            st.rerun()

        st.divider()

        # =====================================================================
        # SECTION B: SAMPLES
        # =====================================================================
        st.header("B. Samples")
        render_sidebar_sample_manager(chrom_map, pda_map)

        st.divider()

        # =====================================================================
        # SECTION C: VISUALIZATION
        # =====================================================================
        st.header("C. Visualization")
        chromatograms = list(chrom_map.values())
        restored_cfg: Optional[PlotConfig] = st.session_state.get("restored_plot_config")

        vis_mode = st.radio(
            "Plot Arrangement Mode",
            options=list(VISUALIZATION_MODES),
            index=(
                list(VISUALIZATION_MODES).index(restored_cfg.mode)
                if restored_cfg and restored_cfg.mode in VISUALIZATION_MODES
                else 0
            ),
            help="Overlay: shared axes per detector. Stacked: vertically separated panels. Individual: single trace.",
        )

        individual_trace_id: Optional[str] = None
        if vis_mode == "Individual" and chromatograms:
            trace_opts = {c.trace_id: f"{c.display_label} ({c.file_name})" for c in chromatograms}
            individual_trace_id = st.selectbox(
                "Select Individual Chromatogram",
                options=list(trace_opts.keys()),
                format_func=lambda tid: trace_opts[tid],
            )

        present_detectors = sorted({c.assigned_detector for c in chromatograms}) or list(VALID_DETECTORS)
        active_detectors = st.multiselect(
            "Active Detectors",
            options=list(VALID_DETECTORS),
            default=present_detectors,
            help="UV/PDA and ELSD signals are automatically plotted in separate synchronized panels.",
        )

        # Wavelength filter when UV/PDA wavelengths exist
        available_wls = sorted(
            {round(c.wavelength_nm, 1) for c in chromatograms if c.wavelength_nm is not None}
        )
        selected_wavelengths: Optional[List[float]] = None
        if available_wls:
            selected_wavelengths = st.multiselect(
                "UV / PDA Wavelengths (nm)",
                options=available_wls,
                default=available_wls,
                help="Filter displayed UV/PDA channels by wavelength (nm).",
            )

        # Axis range controls
        custom_x = st.checkbox(
            "Custom Retention-Time (X) Range",
            value=False,
            help="Specify explicit start and end retention time (minutes).",
        )
        x_range: Optional[Tuple[float, float]] = None
        if custom_x:
            global_min_t = min((c.rt_min for c in chromatograms), default=0.0)
            global_max_t = max((c.rt_max for c in chromatograms), default=32.0)
            xc1, xc2 = st.columns(2)
            xmin = xc1.number_input("X Min (min)", value=float(global_min_t), step=0.5)
            xmax = xc2.number_input("X Max (min)", value=float(global_max_t), step=0.5)
            if xmax > xmin:
                x_range = (float(xmin), float(xmax))

        custom_y = st.checkbox(
            "Custom Intensity (Y) Range",
            value=False,
            help="Specify explicit Y-axis ranges independently for UV/PDA and ELSD panels.",
        )
        y_range_uv: Optional[Tuple[float, float]] = None
        y_range_elsd: Optional[Tuple[float, float]] = None
        if custom_y:
            st.caption("UV / PDA Y-Range (AU)")
            uc1, uc2 = st.columns(2)
            u_min = uc1.number_input("UV Y Min", value=-0.05, step=0.05, format="%.4f")
            u_max = uc2.number_input("UV Y Max", value=1.0, step=0.1, format="%.4f")
            if u_max > u_min:
                y_range_uv = (float(u_min), float(u_max))

            st.caption("ELSD Y-Range (LSU)")
            ec1, ec2 = st.columns(2)
            e_min = ec1.number_input("ELSD Y Min", value=-5.0, step=1.0, format="%.2f")
            e_max = ec2.number_input("ELSD Y Max", value=100.0, step=5.0, format="%.2f")
            if e_max > e_min:
                y_range_elsd = (float(e_min), float(e_max))

        line_width = st.slider(
            "Line Thickness (px)",
            min_value=0.5,
            max_value=5.0,
            value=float(restored_cfg.line_width) if restored_cfg else 1.8,
            step=0.1,
        )
        opacity = st.slider(
            "Line Opacity",
            min_value=0.15,
            max_value=1.0,
            value=float(restored_cfg.opacity) if restored_cfg else 0.90,
            step=0.05,
        )

        stacked_spacing = 0.04
        if vis_mode == "Stacked":
            stacked_spacing = st.slider(
                "Stacked Panel Vertical Spacing",
                min_value=0.01,
                max_value=0.15,
                value=0.04,
                step=0.01,
            )

        col_g1, col_g2 = st.columns(2)
        show_grid = col_g1.checkbox("Show Grid", value=True)
        show_legend = col_g2.checkbox("Show Legend", value=True)
        sync_x_axes = st.checkbox(
            "Synchronize X-Axis Zoom Across Panels",
            value=True,
            help="Keep retention-time axes synchronized between UV/PDA and ELSD subplots.",
        )

        # Data Processing Options (Section 8)
        st.subheader("Data Processing (Visual Only)")
        show_baseline = st.checkbox(
            "Display Estimated Baseline",
            value=False,
            help="Overlay a non-destructive rolling minimum baseline envelope.",
        )
        norm_mode = st.selectbox(
            "Visual Intensity Normalization",
            options=list(NORMALIZATION_MODES),
            index=0,
            help="Optional display scaling. Raw intensity values are never modified.",
        )

        st.divider()

        # =====================================================================
        # SECTION D: ALIGNMENT
        # =====================================================================
        st.header("D. Alignment")
        align_mode = st.selectbox(
            "Alignment Mode",
            options=list(ALIGNMENT_MODES),
            index=0,
            help="Correct retention-time drift without overwriting original retention times.",
        )

        alignment_state: Dict[str, Any] = {
            "mode": align_mode,
            "reference_id": None,
            "region_start_min": None,
            "region_end_min": None,
            "max_shift_min": 1.5,
        }

        if not chromatograms:
            st.caption("Import chromatograms to configure alignment.")
        elif align_mode == "No alignment":
            if any(c.is_aligned for c in chromatograms):
                if st.button("Revert All to Original Retention Times", width="stretch"):
                    reset_all_alignments(chromatograms)
                    st.session_state["proposed_alignment"] = None
                    st.rerun()
        else:
            trace_opts = {c.trace_id: f"{c.display_label} ({c.file_name})" for c in chromatograms}
            ref_id = st.selectbox(
                "Reference Chromatogram",
                options=list(trace_opts.keys()),
                format_func=lambda tid: trace_opts[tid],
                help="Chromatogram against which target retention times are aligned.",
            )
            alignment_state["reference_id"] = ref_id
            ref_chrom = chrom_map[ref_id]

            if align_mode == "Automatic alignment (Cross-Correlation)":
                st.caption(
                    f"Aligning traces with **{ref_chrom.assigned_detector}** detector "
                    + (f"({ref_chrom.wavelength_nm:g} nm)" if ref_chrom.wavelength_nm else "")
                )
                ac1, ac2 = st.columns(2)
                reg_start = ac1.number_input(
                    "Region Start (min)",
                    min_value=0.0,
                    max_value=float(ref_chrom.rt_max),
                    value=float(ref_chrom.rt_min),
                    step=0.5,
                    help="Start of retention-time window used for cross-correlation.",
                )
                reg_end = ac2.number_input(
                    "Region End (min)",
                    min_value=0.0,
                    max_value=float(ref_chrom.rt_max),
                    value=float(ref_chrom.rt_max),
                    step=0.5,
                    help="End of retention-time window used for cross-correlation.",
                )
                max_shift = st.number_input(
                    "Max Permitted Shift (± min)",
                    min_value=0.01,
                    max_value=10.0,
                    value=1.5,
                    step=0.1,
                    help="Maximum retention-time shift allowed during cross-correlation search.",
                )
                same_det_only = st.checkbox(
                    "Align matching detector/wavelength channels only",
                    value=True,
                    help="Recommended (Section 7.3): Align chromatograms sharing the reference detector type.",
                )

                alignment_state["region_start_min"] = float(reg_start)
                alignment_state["region_end_min"] = float(reg_end)
                alignment_state["max_shift_min"] = float(max_shift)

                if same_det_only:
                    targets = [
                        c
                        for c in chromatograms
                        if c.assigned_detector == ref_chrom.assigned_detector
                        and (
                            ref_chrom.wavelength_nm is None
                            or c.wavelength_nm is None
                            or abs(c.wavelength_nm - ref_chrom.wavelength_nm) <= 0.5
                        )
                    ]
                else:
                    targets = list(chromatograms)

                if st.button("🔍 Preview Proposed Alignment", width="stretch"):
                    proposals = estimate_batch_alignment(
                        reference=ref_chrom,
                        targets=targets,
                        region_start_min=float(reg_start),
                        region_end_min=float(reg_end),
                        max_shift_min=float(max_shift),
                    )
                    st.session_state["proposed_alignment"] = proposals
                    st.rerun()

                proposals = st.session_state.get("proposed_alignment")
                if proposals:
                    b_col1, b_col2 = st.columns(2)
                    if b_col1.button("✅ Accept & Apply", width="stretch"):
                        apply_alignment_estimates(chromatograms, proposals)
                        st.session_state["proposed_alignment"] = None
                        st.success("Applied automatic alignment shifts.")
                        st.rerun()
                    if b_col2.button("❌ Reject Preview", width="stretch"):
                        st.session_state["proposed_alignment"] = None
                        st.rerun()

            elif align_mode == "Manual alignment":
                st.caption("Specify additive retention-time shifts Δt (minutes) for each sample:")
                manual_shifts: Dict[str, float] = {}
                for c in chromatograms:
                    if c.trace_id == ref_id:
                        st.text_input(
                            f"{c.display_label} (Reference)",
                            value="0.0000 min",
                            disabled=True,
                            key=f"man_ref_{c.trace_id}",
                        )
                        manual_shifts[c.trace_id] = 0.0
                    else:
                        val = st.number_input(
                            f"Δt (min): {c.display_label}",
                            min_value=-15.0,
                            max_value=15.0,
                            value=float(c.time_shift_min),
                            step=0.01,
                            format="%.4f",
                            key=f"man_shift_{c.trace_id}",
                        )
                        manual_shifts[c.trace_id] = float(val)

                if st.button("✅ Apply Manual Shifts", width="stretch"):
                    apply_manual_shifts(chromatograms, manual_shifts, reference_trace_id=ref_id)
                    st.session_state["proposed_alignment"] = None
                    st.success("Applied manual retention-time shifts.")
                    st.rerun()

            if st.button("↩️ Reset All Alignments", width="stretch"):
                reset_all_alignments(chromatograms)
                st.session_state["proposed_alignment"] = None
                st.rerun()

    exp_settings = st.session_state.get("export_settings", {})
    plot_config = PlotConfig(
        mode=vis_mode,
        active_detectors=active_detectors,
        selected_wavelengths=selected_wavelengths,
        individual_trace_id=individual_trace_id,
        x_range=x_range,
        y_range_uv=y_range_uv,
        y_range_elsd=y_range_elsd,
        line_width=line_width,
        opacity=opacity,
        show_grid=show_grid,
        show_legend=show_legend,
        sync_x_axes=sync_x_axes,
        stacked_spacing=stacked_spacing,
        show_baseline=show_baseline,
        normalization_mode=norm_mode,
        use_aligned_time=True,
        title=exp_settings.get("title", ""),
        x_axis_label=exp_settings.get("x_axis_label", "Retention Time (min)"),
        y_axis_label_uv=exp_settings.get("y_axis_label_uv", "Absorbance (AU)"),
        y_axis_label_elsd=exp_settings.get("y_axis_label_elsd", "ELSD Signal (LSU)"),
        font_size=int(exp_settings.get("font_size", 13)),
        legend_position=exp_settings.get("legend_position", "right"),
    )
    return (plot_config, alignment_state)
