"""Session configuration serialization and state management.

Saves and restores sample display names, colors, detector assignments, selection
visibility, visualization settings, alignment parameters, and export settings
WITHOUT embedding raw experimental data arrays.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from src.data_import.file_validation import VALID_DETECTORS
from src.data_processing.chromatogram import Chromatogram
from src.visualization.chromatogram_plot import PlotConfig


DEFAULT_SETTINGS_PATH = Path(__file__).resolve().parents[2] / "config" / "default_settings.json"


def load_default_settings(path: Optional[Path] = None) -> Dict[str, Any]:
    """Load default settings from `config/default_settings.json`."""
    target = path or DEFAULT_SETTINGS_PATH
    if target.exists():
        try:
            return json.loads(target.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {
        "application": {"name": "HPLC Viewer", "version": "1.0.0"},
        "visualization": {},
        "alignment": {},
        "export": {},
    }


def serialize_session_config(
    chromatograms: Sequence[Chromatogram],
    plot_config: PlotConfig,
    alignment_settings: Dict[str, Any],
    export_settings: Dict[str, Any],
    palette_name: str = "Okabe-Ito (Colorblind Safe)",
) -> str:
    """Serialize session settings and sample metadata to JSON (excluding raw data arrays)."""
    samples_payload: List[Dict[str, Any]] = []
    for chrom in chromatograms:
        samples_payload.append(
            {
                "trace_id": chrom.trace_id,
                "file_name": chrom.file_name,
                "original_sample_name": chrom.original_sample_name,
                "sample_display_name": chrom.sample_display_name,
                "injection_id": chrom.injection_id,
                "channel": chrom.channel,
                "wavelength_nm": chrom.wavelength_nm,
                "assigned_detector": chrom.assigned_detector,
                "color": chrom.color,
                "visible": chrom.visible,
                "is_aligned": chrom.is_aligned,
                "time_shift_min": chrom.time_shift_min,
                "alignment_method": chrom.alignment_method,
                "alignment_reference_id": chrom.alignment_reference_id,
            }
        )

    payload: Dict[str, Any] = {
        "application": {
            "name": "HPLC Viewer",
            "version": "1.0.0",
            "schema": "hplc_viewer_session_v1",
        },
        "samples": samples_payload,
        "visualization": {
            "mode": plot_config.mode,
            "active_detectors": list(plot_config.active_detectors),
            "selected_wavelengths": (
                list(plot_config.selected_wavelengths)
                if plot_config.selected_wavelengths is not None
                else None
            ),
            "individual_trace_id": plot_config.individual_trace_id,
            "x_range": list(plot_config.x_range) if plot_config.x_range else None,
            "y_range_uv": list(plot_config.y_range_uv) if plot_config.y_range_uv else None,
            "y_range_elsd": list(plot_config.y_range_elsd) if plot_config.y_range_elsd else None,
            "line_width": plot_config.line_width,
            "opacity": plot_config.opacity,
            "show_grid": plot_config.show_grid,
            "show_legend": plot_config.show_legend,
            "sync_x_axes": plot_config.sync_x_axes,
            "stacked_spacing": plot_config.stacked_spacing,
            "show_baseline": plot_config.show_baseline,
            "normalization_mode": plot_config.normalization_mode,
            "use_aligned_time": plot_config.use_aligned_time,
            "palette": palette_name,
        },
        "alignment": dict(alignment_settings),
        "export": {
            "figure_format": export_settings.get("figure_format", "png"),
            "width_px": export_settings.get("width_px", 1200),
            "height_px": export_settings.get("height_px", 750),
            "dpi_scale": export_settings.get("dpi_scale", 2),
            "font_size": plot_config.font_size,
            "legend_position": plot_config.legend_position,
            "title": plot_config.title,
            "x_axis_label": plot_config.x_axis_label,
            "y_axis_label_uv": plot_config.y_axis_label_uv,
            "y_axis_label_elsd": plot_config.y_axis_label_elsd,
        },
    }
    return json.dumps(payload, indent=2)


def restore_session_config(
    json_content: str,
    chromatograms: Sequence[Chromatogram],
) -> Tuple[PlotConfig, Dict[str, Any], Dict[str, Any], str, int]:
    """Apply saved session JSON to currently loaded chromatograms and return settings.

    Returns
    -------
    Tuple[PlotConfig, Dict[str, Any], Dict[str, Any], str, int]
        (restored_plot_config, alignment_settings, export_settings, palette_name, matched_sample_count)
    """
    data = json.loads(json_content)
    if not isinstance(data, dict):
        raise ValueError("Invalid session configuration JSON: root must be a JSON object.")

    samples_list = data.get("samples", [])
    by_trace_id = {s.get("trace_id"): s for s in samples_list if isinstance(s, dict) and s.get("trace_id")}
    by_file_ch = {
        (s.get("file_name"), s.get("channel")): s
        for s in samples_list
        if isinstance(s, dict) and s.get("file_name")
    }

    matched_count = 0
    for chrom in chromatograms:
        saved = by_trace_id.get(chrom.trace_id) or by_file_ch.get((chrom.file_name, chrom.channel))
        if saved is None:
            continue
        matched_count += 1
        if saved.get("sample_display_name"):
            chrom.sample_display_name = str(saved["sample_display_name"])
        det = saved.get("assigned_detector")
        if det in VALID_DETECTORS:
            chrom.assigned_detector = str(det)
        if saved.get("color"):
            chrom.color = str(saved["color"])
        if "visible" in saved:
            chrom.visible = bool(saved["visible"])
        if saved.get("is_aligned"):
            chrom.apply_shift(
                shift_min=float(saved.get("time_shift_min", 0.0)),
                method=str(saved.get("alignment_method", "Restored")),
                reference_id=saved.get("alignment_reference_id"),
            )
        else:
            chrom.reset_alignment()

    vis = data.get("visualization", {})
    exp = data.get("export", {})
    align = data.get("alignment", {})

    x_rng = tuple(vis["x_range"]) if isinstance(vis.get("x_range"), list) and len(vis["x_range"]) == 2 else None
    y_uv = tuple(vis["y_range_uv"]) if isinstance(vis.get("y_range_uv"), list) and len(vis["y_range_uv"]) == 2 else None
    y_el = tuple(vis["y_range_elsd"]) if isinstance(vis.get("y_range_elsd"), list) and len(vis["y_range_elsd"]) == 2 else None

    plot_cfg = PlotConfig(
        mode=str(vis.get("mode", "Overlay")),
        active_detectors=list(vis.get("active_detectors", ["UV", "PDA", "ELSD", "Other / Unknown"])),
        selected_wavelengths=vis.get("selected_wavelengths"),
        individual_trace_id=vis.get("individual_trace_id"),
        x_range=x_rng,  # type: ignore[arg-type]
        y_range_uv=y_uv,  # type: ignore[arg-type]
        y_range_elsd=y_el,  # type: ignore[arg-type]
        line_width=float(vis.get("line_width", 1.8)),
        opacity=float(vis.get("opacity", 0.90)),
        show_grid=bool(vis.get("show_grid", True)),
        show_legend=bool(vis.get("show_legend", True)),
        sync_x_axes=bool(vis.get("sync_x_axes", True)),
        stacked_spacing=float(vis.get("stacked_spacing", 0.04)),
        show_baseline=bool(vis.get("show_baseline", False)),
        normalization_mode=str(vis.get("normalization_mode", "None")),
        use_aligned_time=bool(vis.get("use_aligned_time", True)),
        title=str(exp.get("title", "")),
        x_axis_label=str(exp.get("x_axis_label", "Retention Time (min)")),
        y_axis_label_uv=str(exp.get("y_axis_label_uv", "Absorbance (AU)")),
        y_axis_label_elsd=str(exp.get("y_axis_label_elsd", "ELSD Signal (LSU)")),
        font_size=int(exp.get("font_size", 13)),
        legend_position=str(exp.get("legend_position", "right")),
    )
    palette_name = str(vis.get("palette", "Okabe-Ito (Colorblind Safe)"))
    return (plot_cfg, align, exp, palette_name, matched_count)
