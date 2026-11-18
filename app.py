"""Main Streamlit entry point for the HPLC Chromatogram Viewer application."""

from __future__ import annotations

from typing import Any, Dict, List
import pandas as pd
import streamlit as st

from src.data_processing.alignment import AlignmentEstimate
from src.data_processing.chromatogram import Chromatogram
from src.interface.import_interface import render_import_summary_and_table
from src.interface.sidebar import render_sidebar
from src.utils.logging_config import get_logger
from src.utils.session_manager import serialize_session_config
from src.visualization.chromatogram_plot import (
    PlotConfig,
    create_alignment_comparison_figure,
    create_chromatogram_figure,
    filter_chromatograms_for_plot,
)
from src.visualization.figure_export import (
    SUPPORTED_FIGURE_FORMATS,
    export_chromatograms_to_csv_bytes,
    export_figure_bytes,
)

logger = get_logger("hplc_viewer.app")


def _render_alignment_panel(
    chromatograms: List[Chromatogram],
    alignment_state: Dict[str, Any],
    line_width: float,
) -> None:
    """Render Section 3 of the main panel: Alignment Status and Before/After preview."""
    proposals: List[AlignmentEstimate] | None = st.session_state.get("proposed_alignment")
    aligned_traces = [c for c in chromatograms if c.is_aligned]

    if not proposals and not aligned_traces and alignment_state.get("mode") == "No alignment":
        return

    st.subheader("3. Alignment Status & Before / After Comparison")

    if proposals:
        st.info(
            "**Alignment Preview Active** — Review the proposed cross-correlation shifts and "
            "Before/After chromatograms below, then click **Accept & Apply** or **Reject Preview** in the sidebar."
        )
        prop_rows = [
            {
                "Sample Name": p.sample_name,
                "File Name": p.file_name,
                "Detector": p.detector,
                "Role": "Reference" if p.is_reference else "Target",
                "Proposed Shift Δt (min)": round(p.estimated_shift_min, 5),
                "Proposed Shift (sec)": round(p.estimated_shift_min * 60.0, 2),
                "Correlation Score": round(p.correlation_score, 4),
                "Alignment Window (min)": f"{p.region_start_min:.2f} – {p.region_end_min:.2f}",
                "Notes": p.warning or "OK",
            }
            for p in proposals
        ]
        st.dataframe(pd.DataFrame(prop_rows), width="stretch", hide_index=True)

        target_ids = {p.trace_id for p in proposals}
        preview_traces = [c for c in chromatograms if c.trace_id in target_ids]
        comp_fig = create_alignment_comparison_figure(
            chromatograms=preview_traces,
            proposed_estimates=proposals,
            region_start_min=alignment_state.get("region_start_min"),
            region_end_min=alignment_state.get("region_end_min"),
            line_width=line_width,
        )
        st.plotly_chart(
            comp_fig,
            width="stretch",
            config={"scrollZoom": True, "displaylogo": False},
        )
    elif aligned_traces:
        status_rows = [
            {
                "Sample Name": c.sample_display_name,
                "File Name": c.file_name,
                "Detector": c.assigned_detector,
                "Channel": c.channel,
                "Applied Shift Δt (min)": round(c.time_shift_min, 5),
                "Applied Shift (sec)": round(c.time_shift_min * 60.0, 2),
                "Method": c.alignment_method,
            }
            for c in chromatograms
            if c.is_aligned
        ]
        st.dataframe(pd.DataFrame(status_rows), width="stretch", hide_index=True)

        with st.expander("📈 View Before vs. After Alignment Chromatograms", expanded=True):
            comp_fig = create_alignment_comparison_figure(
                chromatograms=aligned_traces,
                proposed_estimates=None,
                region_start_min=alignment_state.get("region_start_min"),
                region_end_min=alignment_state.get("region_end_min"),
                line_width=line_width,
            )
            st.plotly_chart(
                comp_fig,
                width="stretch",
                config={"scrollZoom": True, "displaylogo": False},
            )


def _render_export_panel(
    chromatograms: List[Chromatogram],
    plot_config: PlotConfig,
    alignment_state: Dict[str, Any],
) -> None:
    """Render Section 4 of the main panel: Publication Figure, CSV Data, and Session Config Exports."""
    st.subheader("4. Export Controls (Figures, Data & Session Configuration)")
    exp_settings: Dict[str, Any] = st.session_state.setdefault(
        "export_settings",
        {
            "figure_format": "png",
            "width_px": 1200,
            "height_px": 750,
            "dpi_scale": 2,
            "font_size": 13,
            "legend_position": "right",
            "title": "",
            "x_axis_label": "Retention Time (min)",
            "y_axis_label_uv": "Absorbance (AU)",
            "y_axis_label_elsd": "ELSD Signal (LSU)",
        },
    )

    tab_fig, tab_csv, tab_cfg = st.tabs(
        ["🖼️ Publication Figure Export", "📄 Chromatographic Data (CSV)", "⚙️ Session Configuration (.json)"]
    )

    with tab_fig:
        c1, c2, c3, c4 = st.columns(4)
        fig_title = c1.text_input(
            "Figure Title",
            value=str(exp_settings.get("title", "")),
            placeholder="Optional figure title",
        )
        x_label = c2.text_input(
            "X-Axis Label",
            value=str(exp_settings.get("x_axis_label", "Retention Time (min)")),
        )
        y_label_uv = c3.text_input(
            "UV / PDA Y-Axis Label",
            value=str(exp_settings.get("y_axis_label_uv", "Absorbance (AU)")),
        )
        y_label_elsd = c4.text_input(
            "ELSD Y-Axis Label",
            value=str(exp_settings.get("y_axis_label_elsd", "ELSD Signal (LSU)")),
        )

        d1, d2, d3, d4, d5, d6 = st.columns(6)
        fmt = d1.selectbox(
            "Format",
            options=list(SUPPORTED_FIGURE_FORMATS),
            index=list(SUPPORTED_FIGURE_FORMATS).index(exp_settings.get("figure_format", "png")),
            format_func=lambda s: s.upper(),
        )
        width_px = int(
            d2.number_input("Width (px)", min_value=400, max_value=4000, value=int(exp_settings.get("width_px", 1200)), step=50)
        )
        height_px = int(
            d3.number_input("Height (px)", min_value=300, max_value=3000, value=int(exp_settings.get("height_px", 750)), step=50)
        )
        dpi_scale = int(
            d4.selectbox(
                "Resolution Scale",
                options=[1, 2, 3, 4],
                index=1,
                format_func=lambda s: f"{s}x (~{s*100} DPI)",
            )
        )
        font_size = int(
            d5.number_input("Font Size (pt)", min_value=8, max_value=24, value=int(exp_settings.get("font_size", 13)), step=1)
        )
        legend_pos = d6.selectbox(
            "Legend Position",
            options=["right", "top", "bottom", "inside"],
            index=["right", "top", "bottom", "inside"].index(exp_settings.get("legend_position", "right")),
        )

        exp_settings.update(
            {
                "figure_format": fmt,
                "width_px": width_px,
                "height_px": height_px,
                "dpi_scale": dpi_scale,
                "font_size": font_size,
                "legend_position": legend_pos,
                "title": fig_title,
                "x_axis_label": x_label,
                "y_axis_label_uv": y_label_uv,
                "y_axis_label_elsd": y_label_elsd,
            }
        )

        export_plot_cfg = PlotConfig(
            mode=plot_config.mode,
            active_detectors=plot_config.active_detectors,
            selected_wavelengths=plot_config.selected_wavelengths,
            individual_trace_id=plot_config.individual_trace_id,
            x_range=plot_config.x_range,
            y_range_uv=plot_config.y_range_uv,
            y_range_elsd=plot_config.y_range_elsd,
            line_width=plot_config.line_width,
            opacity=plot_config.opacity,
            show_grid=plot_config.show_grid,
            show_legend=plot_config.show_legend,
            sync_x_axes=plot_config.sync_x_axes,
            stacked_spacing=plot_config.stacked_spacing,
            show_baseline=plot_config.show_baseline,
            normalization_mode=plot_config.normalization_mode,
            use_aligned_time=plot_config.use_aligned_time,
            title=fig_title,
            x_axis_label=x_label,
            y_axis_label_uv=y_label_uv,
            y_axis_label_elsd=y_label_elsd,
            font_size=font_size,
            legend_position=legend_pos,
            height_px=height_px,
            width_px=width_px,
        )

        mime_map = {
            "png": "image/png",
            "svg": "image/svg+xml",
            "pdf": "application/pdf",
            "html": "text/html",
        }
        export_fig = create_chromatogram_figure(chromatograms, export_plot_cfg)
        fig_bytes = export_figure_bytes(
            export_fig,
            fmt=fmt,
            width_px=width_px,
            height_px=height_px,
            scale=dpi_scale,
        )
        st.download_button(
            label=f"⬇️ Download Publication Figure ({fmt.upper()})",
            data=fig_bytes,
            file_name=f"hplc_chromatogram.{fmt}",
            mime=mime_map[fmt],
        )

    with tab_csv:
        st.caption(
            "Exports original retention time (`original_retention_time_min`), original intensity (`original_intensity`), "
            "and corrected retention time (`corrected_retention_time_min`) as separate columns. "
            "Original experimental observations are never overwritten."
        )
        scope = st.radio(
            "Export Scope",
            options=["Currently Plotted (Visible) Chromatograms", "All Imported Chromatograms"],
            horizontal=True,
        )
        target_traces = (
            filter_chromatograms_for_plot(chromatograms, plot_config)
            if scope.startswith("Currently")
            else chromatograms
        )
        csv_bytes = export_chromatograms_to_csv_bytes(target_traces)
        st.download_button(
            label=f"⬇️ Download Chromatographic Data CSV ({len(target_traces)} channel(s))",
            data=csv_bytes,
            file_name="hplc_chromatograms_export.csv",
            mime="text/csv",
            disabled=(len(target_traces) == 0),
        )

    with tab_cfg:
        st.caption(
            "Save sample display names, colors, detector assignments, visibility, plot settings, "
            "and alignment parameters to a lightweight `.json` file (no experimental data arrays embedded)."
        )
        palette_name = st.session_state.get("palette_name", "Okabe-Ito (Colorblind Safe)")
        session_json = serialize_session_config(
            chromatograms=chromatograms,
            plot_config=plot_config,
            alignment_settings=alignment_state,
            export_settings=exp_settings,
            palette_name=palette_name,
        )
        st.download_button(
            label="⬇️ Download Session Configuration (.json)",
            data=session_json.encode("utf-8"),
            file_name="hplc_viewer_session.json",
            mime="application/json",
        )


def main() -> None:
    """Run the HPLC Chromatogram Viewer Streamlit application."""
    st.set_page_config(
        page_title="HPLC Viewer — Chromatogram Analysis & Alignment",
        page_icon="🔬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    plot_config, alignment_state = render_sidebar()
    chrom_map: Dict[str, Chromatogram] = st.session_state.get("chromatograms", {})
    chromatograms = list(chrom_map.values())

    st.markdown(
        "### HPLC Chromatogram Viewer\n"
        "Local scientific environment for importing, inspecting, comparing, and aligning "
        "Waters Empower 3 (`.arw`), `.csv`, and `.txt` chromatographic data."
    )

    if not chromatograms and not st.session_state.get("import_errors"):
        st.info(
            "👈 **Get Started:** Drag and drop Waters Empower `.arw`, `.csv`, or `.txt` files into "
            "**Sidebar → A. Import**, or click **Load Example Set (3 Samples: UV + ELSD)** under "
            "**Load from Local `raw_data/` Folder**."
        )
        return

    # 1. Summary of imported files & editable table
    st.subheader("1. Summary of Imported Files")
    render_import_summary_and_table(chromatograms)

    if not chromatograms:
        return

    # 2. Chromatogram Visualization
    st.subheader("2. Chromatogram Visualization")
    fig = create_chromatogram_figure(chromatograms, plot_config)
    st.plotly_chart(
        fig,
        width="stretch",
        config={
            "scrollZoom": True,
            "displaylogo": False,
            "modeBarButtonsToAdd": ["v1hovermode", "toggleSpikelines"],
        },
    )

    # 3. Alignment Status & Before/After Preview
    _render_alignment_panel(chromatograms, alignment_state, line_width=plot_config.line_width)

    st.divider()

    # 4. Export Controls
    _render_export_panel(chromatograms, plot_config, alignment_state)


if __name__ == "__main__":
    main()
