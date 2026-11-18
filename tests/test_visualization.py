"""Tests for chromatogram visualization, color management, figure/CSV export, and session save/restore."""

import json
from pathlib import Path
import numpy as np
import pandas as pd

from src.data_import.arw_reader import read_arw
from src.data_import.csv_reader import read_csv_chromatogram
from src.utils.session_manager import restore_session_config, serialize_session_config
from src.visualization.chromatogram_plot import (
    PlotConfig,
    create_alignment_comparison_figure,
    create_chromatogram_figure,
)
from src.visualization.color_manager import (
    assign_colors_to_chromatograms,
    reset_colors_to_palette,
)
from src.visualization.figure_export import (
    export_chromatograms_to_csv_bytes,
    export_figure_bytes,
)


RAW_DATA_DIR = Path(__file__).resolve().parents[1] / "raw_data"


def _load_sample_triplet():
    """Load real ELSD, UV 254nm, and UV 365nm chromatograms for sample 01006-18-1."""
    c_elsd = read_arw(RAW_DATA_DIR / "HPLC_RAW_EXPORT1276.arw").chromatograms[0]
    c_254 = read_arw(RAW_DATA_DIR / "HPLC_RAW_EXPORT1278.arw").chromatograms[0]
    c_365 = read_arw(RAW_DATA_DIR / "HPLC_RAW_EXPORT1279.arw").chromatograms[0]
    return [c_elsd, c_254, c_365]


def test_detector_separation_and_wavelength_filtering_in_overlay() -> None:
    """Verify UV and ELSD are plotted in separate synchronized panels and wavelengths filter cleanly."""
    traces = _load_sample_triplet()
    assign_colors_to_chromatograms(traces, group_by_sample_name=False)

    # Both UV and ELSD active -> 2 subplot rows (y, y2), 3 traces total
    cfg_both = PlotConfig(mode="Overlay", active_detectors=["UV", "PDA", "ELSD"])
    fig_both = create_chromatogram_figure(traces, cfg_both)
    assert len(fig_both.data) == 3
    y_axes_used = {trace.yaxis for trace in fig_both.data}
    assert y_axes_used == {"y", "y2"}

    # Filter to only 254 nm -> 1 UV trace + 1 ELSD trace = 2 traces
    cfg_254_only = PlotConfig(
        mode="Overlay",
        active_detectors=["UV", "PDA", "ELSD"],
        selected_wavelengths=[254.0],
    )
    fig_254 = create_chromatogram_figure(traces, cfg_254_only)
    assert len(fig_254.data) == 2

    # Hide ELSD -> 1 subplot row (y only)
    cfg_uv_only = PlotConfig(mode="Overlay", active_detectors=["UV", "PDA"])
    fig_uv = create_chromatogram_figure(traces, cfg_uv_only)
    assert len(fig_uv.data) == 2
    assert {trace.yaxis for trace in fig_uv.data} == {"y"}


def test_stacked_and_individual_modes_and_color_consistency() -> None:
    """Verify Stacked and Individual plotting modes preserve consistent colors and retention times."""
    traces = _load_sample_triplet()
    assign_colors_to_chromatograms(traces, palette_name="Nature / NPG", group_by_sample_name=False)
    initial_colors = [c.color for c in traces]

    # Stacked mode creates 3 vertically separated subplot rows
    cfg_stacked = PlotConfig(mode="Stacked", show_baseline=True)
    fig_stacked = create_chromatogram_figure(traces, cfg_stacked)
    # 3 main traces + 3 baseline traces = 6 traces across y, y2, y3
    assert len(fig_stacked.data) == 6
    assert {t.yaxis for t in fig_stacked.data} == {"y", "y2", "y3"}
    main_colors = [t.line.color for t in fig_stacked.data if t.showlegend is not False]
    assert main_colors == initial_colors

    # Individual mode displays only the requested trace
    cfg_indiv = PlotConfig(mode="Individual", individual_trace_id=traces[1].trace_id)
    fig_indiv = create_chromatogram_figure(traces, cfg_indiv)
    assert len(fig_indiv.data) == 1
    assert fig_indiv.data[0].line.color == initial_colors[1]
    assert np.array_equal(fig_indiv.data[0].x, traces[1].raw_retention_time)

    # Reset colors to default Okabe-Ito palette
    reset_colors_to_palette(traces, group_by_sample_name=True)
    assert len({c.color for c in traces}) == 1  # Same sample name -> consistent color


def test_aligned_retention_time_axis_and_comparison_plot() -> None:
    """Verify aligned chromatograms plot corrected RT while preserving raw RT."""
    traces = _load_sample_triplet()
    c_254 = traces[1]
    c_254.apply_shift(-0.15, method="Manual")

    fig_aligned = create_chromatogram_figure(
        [c_254], PlotConfig(mode="Overlay", use_aligned_time=True)
    )
    assert np.allclose(fig_aligned.data[0].x, c_254.raw_retention_time - 0.15)

    fig_raw = create_chromatogram_figure(
        [c_254], PlotConfig(mode="Overlay", use_aligned_time=False)
    )
    assert np.allclose(fig_raw.data[0].x, c_254.raw_retention_time)

    comp_fig = create_alignment_comparison_figure([c_254], region_start_min=2.0, region_end_min=15.0)
    assert len(comp_fig.data) == 2


def test_csv_figure_and_session_exports() -> None:
    """Verify CSV data export, static/HTML figure export, and session JSON save/restore."""
    traces = _load_sample_triplet()
    traces[1].sample_display_name = "Renamed_Sample_1"
    traces[1].apply_shift(0.08, method="Manual")

    # 1. CSV Export preserves raw RT, raw intensity, and separate corrected RT
    csv_bytes = export_chromatograms_to_csv_bytes([traces[1]])
    assert len(csv_bytes) > 1000
    df = pd.read_csv(pd.io.common.BytesIO(csv_bytes))
    assert "original_retention_time_min" in df.columns
    assert "original_intensity" in df.columns
    assert "corrected_retention_time_min" in df.columns
    assert np.allclose(df["original_retention_time_min"], traces[1].raw_retention_time)
    assert np.allclose(df["original_intensity"], traces[1].raw_intensity)
    assert np.allclose(
        df["corrected_retention_time_min"],
        traces[1].raw_retention_time + 0.08,
    )

    # Round-trip re-import of exported CSV
    reimported = read_csv_chromatogram(csv_bytes, file_name="exported.csv")
    assert reimported.is_valid
    assert reimported.chromatograms[0].sample_display_name == "Renamed_Sample_1"
    assert np.isclose(reimported.chromatograms[0].time_shift_min, 0.08)

    # 2. Figure Export (PNG, SVG, PDF, HTML)
    fig = create_chromatogram_figure([traces[1]], PlotConfig(title="Export Test"))
    for fmt in ("png", "svg", "pdf", "html"):
        out = export_figure_bytes(fig, fmt=fmt, width_px=800, height_px=500, scale=1)
        assert len(out) > 500

    # 3. Session configuration save & restore (must NOT embed raw numerical arrays)
    plot_cfg = PlotConfig(mode="Stacked", line_width=2.4, title="Saved Session")
    session_json = serialize_session_config(
        chromatograms=traces,
        plot_config=plot_cfg,
        alignment_settings={"mode": "Manual alignment", "max_shift_min": 2.0},
        export_settings={"figure_format": "svg", "width_px": 1400, "height_px": 800, "dpi_scale": 2},
    )
    parsed_json = json.loads(session_json)
    assert "raw_retention_time" not in session_json
    assert "raw_intensity" not in session_json
    assert len(parsed_json["samples"]) == 3

    # Modify traces and restore from session_json
    traces[1].sample_display_name = "TempName"
    traces[1].reset_alignment()
    restored_cfg, align_cfg, exp_cfg, pal, matched = restore_session_config(session_json, traces)
    assert matched == 3
    assert traces[1].sample_display_name == "Renamed_Sample_1"
    assert traces[1].is_aligned
    assert np.isclose(traces[1].time_shift_min, 0.08)
    assert restored_cfg.mode == "Stacked"
    assert np.isclose(restored_cfg.line_width, 2.4)
    assert exp_cfg["figure_format"] == "svg"


def test_streamlit_app_end_to_end() -> None:
    """Run headless Streamlit AppTest on app.py to verify full UI lifecycle."""
    from streamlit.testing.v1 import AppTest

    app_path = str(Path(__file__).resolve().parents[1] / "app.py")
    at = AppTest.from_file(app_path, default_timeout=20)
    at.run()
    assert not at.exception

    # Click 'Load Example Set (3 Samples: UV + ELSD)'
    load_btn = next(b for b in at.button if "Load Example Set" in b.label)
    load_btn.click().run(timeout=20)
    assert not at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics.get("Chromatogram Channels") == "9"
    assert metrics.get("Unique Samples") == "3"

    # Click 'Clear Current Session'
    clear_btn = next(b for b in at.button if "Clear Current Session" in b.label)
    clear_btn.click().run(timeout=20)
    assert not at.exception
    assert len(at.metric) == 0

