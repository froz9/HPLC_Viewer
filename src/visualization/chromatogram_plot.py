"""Interactive Plotly chromatogram plotting engine.

Implements Overlay, Stacked, and Individual visualization modes, synchronized
separate UV/PDA and ELSD detector panels, optional baseline overlay, optional
display normalization, and Before/After alignment comparison plots.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.data_processing.alignment import AlignmentEstimate
from src.data_processing.chromatogram import Chromatogram
from src.data_processing.normalization import (
    compute_visual_intensity,
    estimate_visual_baseline,
)


VISUALIZATION_MODES = ("Overlay", "Stacked", "Individual")


@dataclass
class PlotConfig:
    """Configuration parameters for chromatogram figure generation."""

    mode: str = "Overlay"  # "Overlay", "Stacked", "Individual"
    active_detectors: List[str] = field(
        default_factory=lambda: ["UV", "PDA", "ELSD", "Other / Unknown"]
    )
    selected_wavelengths: Optional[List[float]] = None
    individual_trace_id: Optional[str] = None
    x_range: Optional[Tuple[float, float]] = None
    y_range_uv: Optional[Tuple[float, float]] = None
    y_range_elsd: Optional[Tuple[float, float]] = None
    line_width: float = 1.8
    opacity: float = 0.90
    show_grid: bool = True
    show_legend: bool = True
    sync_x_axes: bool = True
    stacked_spacing: float = 0.04
    show_baseline: bool = False
    normalization_mode: str = "None"
    use_aligned_time: bool = True
    title: str = ""
    x_axis_label: str = "Retention Time (min)"
    y_axis_label_uv: str = "Absorbance (AU)"
    y_axis_label_elsd: str = "ELSD Signal (LSU)"
    font_size: int = 13
    legend_position: str = "right"  # "right", "top", "bottom", "inside"
    height_px: Optional[int] = None
    width_px: Optional[int] = None


def filter_chromatograms_for_plot(
    chromatograms: Sequence[Chromatogram],
    config: PlotConfig,
) -> List[Chromatogram]:
    """Filter chromatograms by visibility, active detectors, wavelength, and mode."""
    visible_traces = [c for c in chromatograms if c.visible]

    if config.mode == "Individual" and config.individual_trace_id:
        indiv = [c for c in chromatograms if c.trace_id == config.individual_trace_id]
        if indiv:
            return indiv

    filtered: List[Chromatogram] = []
    for chrom in visible_traces:
        if chrom.assigned_detector not in config.active_detectors:
            continue
        if chrom.assigned_detector in ("UV", "PDA") and config.selected_wavelengths is not None:
            if chrom.wavelength_nm is not None:
                # Match within 0.5 nm tolerance
                if not any(abs(chrom.wavelength_nm - wl) <= 0.5 for wl in config.selected_wavelengths):
                    continue
        filtered.append(chrom)

    if config.mode == "Individual" and filtered:
        return [filtered[0]]

    return filtered


def _line_dash_for_wavelength(chrom: Chromatogram, has_multiple_uv_wls: bool) -> str:
    """Assign distinct line dash style when multiple UV wavelengths of a sample are overlaid."""
    if not has_multiple_uv_wls or chrom.wavelength_nm is None:
        return "solid"
    wl = chrom.wavelength_nm
    if abs(wl - 254.0) < 2.0:
        return "solid"
    if abs(wl - 365.0) < 2.0:
        return "dash"
    if abs(wl - 280.0) < 2.0:
        return "dot"
    if abs(wl - 210.0) < 2.0:
        return "dashdot"
    return "solid"


def _apply_legend_layout(fig: go.Figure, position: str, show_legend: bool, font_size: int) -> None:
    """Configure Plotly legend placement."""
    if not show_legend:
        fig.update_layout(showlegend=False)
        return

    legend_cfg = dict(font=dict(size=max(10, font_size - 1)))
    if position == "top":
        legend_cfg.update(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.0)
    elif position == "bottom":
        legend_cfg.update(orientation="h", yanchor="top", y=-0.18, xanchor="left", x=0.0)
    elif position == "inside":
        legend_cfg.update(
            orientation="v",
            yanchor="top",
            y=0.98,
            xanchor="right",
            x=0.98,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="rgba(180,180,180,0.6)",
            borderwidth=1,
        )
    else:  # "right"
        legend_cfg.update(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02)

    fig.update_layout(showlegend=True, legend=legend_cfg)


def create_chromatogram_figure(
    chromatograms: Sequence[Chromatogram],
    config: Optional[PlotConfig] = None,
) -> go.Figure:
    """Build an interactive Plotly figure for the selected chromatograms and visualization mode."""
    if config is None:
        config = PlotConfig()

    traces = filter_chromatograms_for_plot(chromatograms, config)
    if not traces:
        fig = go.Figure()
        fig.update_layout(
            template="plotly_white",
            height=420,
            annotations=[
                dict(
                    text="No chromatograms match the current sample/detector selection.",
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    font=dict(size=15, color="#555555"),
                )
            ],
        )
        return fig

    if config.mode == "Stacked":
        return _create_stacked_figure(traces, config)

    # Mode 1 (Overlay) and Mode 3 (Individual): separate panels by detector family
    uv_traces = [c for c in traces if c.assigned_detector in ("UV", "PDA")]
    elsd_traces = [c for c in traces if c.assigned_detector == "ELSD"]
    other_traces = [c for c in traces if c.assigned_detector == "Other / Unknown"]

    panels: List[Tuple[str, str, List[Chromatogram], Optional[Tuple[float, float]]]] = []
    if uv_traces:
        panels.append(("UV / PDA Detector", config.y_axis_label_uv, uv_traces, config.y_range_uv))
    if elsd_traces:
        panels.append(("ELSD Detector", config.y_axis_label_elsd, elsd_traces, config.y_range_elsd))
    if other_traces:
        panels.append(("Other / Unknown Detector", "Signal Intensity", other_traces, None))

    n_rows = len(panels)
    default_height = 500 if n_rows == 1 else 380 * n_rows
    fig_height = config.height_px or default_height

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        shared_xaxes=config.sync_x_axes if n_rows > 1 else False,
        vertical_spacing=0.10 if n_rows > 1 else 0.0,
        subplot_titles=[p[0] for p in panels] if n_rows > 1 else None,
    )

    uv_wls = {c.wavelength_nm for c in uv_traces if c.wavelength_nm is not None}
    has_multiple_uv_wls = len(uv_wls) > 1

    for row_idx, (panel_title, base_ylabel, panel_traces, y_range) in enumerate(panels, start=1):
        norm_suffix = ""
        for chrom in panel_traces:
            rt = chrom.active_retention_time if config.use_aligned_time else chrom.raw_retention_time
            y_vis, norm_suffix = compute_visual_intensity(
                chrom.raw_intensity, mode=config.normalization_mode
            )
            dash_style = _line_dash_for_wavelength(chrom, has_multiple_uv_wls)
            shift_note = (
                f"<br>RT Shift: {chrom.time_shift_min:+.4f} min ({chrom.alignment_method})"
                if (config.use_aligned_time and chrom.is_aligned)
                else ""
            )
            wl_note = f" ({chrom.wavelength_nm:g} nm)" if chrom.wavelength_nm is not None else ""
            unit_str = norm_suffix or chrom.units or ""

            fig.add_trace(
                go.Scatter(
                    x=rt,
                    y=y_vis,
                    mode="lines",
                    name=chrom.display_label,
                    legendgroup=chrom.trace_id,
                    line=dict(
                        color=chrom.color,
                        width=config.line_width,
                        dash=dash_style,
                    ),
                    opacity=config.opacity,
                    hovertemplate=(
                        f"<b>{chrom.sample_display_name}</b><br>"
                        f"Detector: {chrom.assigned_detector}{wl_note}<br>"
                        "Retention Time: %{x:.4f} min<br>"
                        f"Intensity: %{{y:.5g}} {unit_str}"
                        f"{shift_note}<extra></extra>"
                    ),
                ),
                row=row_idx,
                col=1,
            )

            if config.show_baseline:
                base_y = estimate_visual_baseline(rt, y_vis)
                fig.add_trace(
                    go.Scatter(
                        x=rt,
                        y=base_y,
                        mode="lines",
                        name=f"{chrom.sample_display_name} (Baseline)",
                        legendgroup=chrom.trace_id,
                        showlegend=False,
                        line=dict(color=chrom.color, width=1.0, dash="dot"),
                        opacity=max(0.35, config.opacity * 0.55),
                        hoverinfo="skip",
                    ),
                    row=row_idx,
                    col=1,
                )

        ylabel = f"{base_ylabel} [{norm_suffix}]" if norm_suffix else base_ylabel
        fig.update_yaxes(
            title_text=ylabel,
            showgrid=config.show_grid,
            gridcolor="rgba(200,200,200,0.45)",
            zeroline=True,
            zerolinecolor="rgba(160,160,160,0.5)",
            range=list(y_range) if y_range is not None else None,
            row=row_idx,
            col=1,
        )
        fig.update_xaxes(
            showgrid=config.show_grid,
            gridcolor="rgba(200,200,200,0.45)",
            range=list(config.x_range) if config.x_range is not None else None,
            row=row_idx,
            col=1,
        )

    # Label bottom X-axis (and all X-axes if not shared)
    for r in range(1, n_rows + 1):
        if r == n_rows or not config.sync_x_axes:
            fig.update_xaxes(title_text=config.x_axis_label, row=r, col=1)

    fig.update_layout(
        template="plotly_white",
        title=dict(text=config.title, x=0.01, xanchor="left") if config.title else None,
        height=fig_height,
        width=config.width_px,
        font=dict(family="Arial, sans-serif", size=config.font_size, color="#1F2937"),
        hovermode="closest",
        margin=dict(l=70, r=40, t=55 if config.title or n_rows > 1 else 30, b=55),
    )
    _apply_legend_layout(fig, config.legend_position, config.show_legend, config.font_size)
    return fig


def _create_stacked_figure(
    traces: Sequence[Chromatogram],
    config: PlotConfig,
) -> go.Figure:
    """Create a vertically stacked multi-panel figure with one chromatogram per subplot row."""
    n_rows = len(traces)
    spacing = float(np.clip(config.stacked_spacing, 0.005, max(0.01, 0.85 / max(1, n_rows - 1))))
    default_height = max(380, min(2400, 210 * n_rows))
    fig_height = config.height_px or default_height

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        shared_xaxes=config.sync_x_axes,
        vertical_spacing=spacing if n_rows > 1 else 0.0,
        subplot_titles=[
            f"{c.sample_display_name} — {c.assigned_detector}"
            + (f" ({c.wavelength_nm:g} nm)" if c.wavelength_nm is not None else "")
            for c in traces
        ],
    )

    for idx, chrom in enumerate(traces, start=1):
        rt = chrom.active_retention_time if config.use_aligned_time else chrom.raw_retention_time
        y_vis, norm_suffix = compute_visual_intensity(
            chrom.raw_intensity, mode=config.normalization_mode
        )
        unit_str = norm_suffix or chrom.units or ("AU" if chrom.assigned_detector in ("UV", "PDA") else "LSU")

        fig.add_trace(
            go.Scatter(
                x=rt,
                y=y_vis,
                mode="lines",
                name=chrom.display_label,
                legendgroup=chrom.trace_id,
                line=dict(color=chrom.color, width=config.line_width),
                opacity=config.opacity,
                hovertemplate=(
                    f"<b>{chrom.sample_display_name}</b><br>"
                    "Retention Time: %{x:.4f} min<br>"
                    f"Intensity: %{{y:.5g}} {unit_str}<extra></extra>"
                ),
            ),
            row=idx,
            col=1,
        )

        if config.show_baseline:
            base_y = estimate_visual_baseline(rt, y_vis)
            fig.add_trace(
                go.Scatter(
                    x=rt,
                    y=base_y,
                    mode="lines",
                    name=f"{chrom.sample_display_name} (Baseline)",
                    showlegend=False,
                    line=dict(color=chrom.color, width=1.0, dash="dot"),
                    opacity=0.5,
                    hoverinfo="skip",
                ),
                row=idx,
                col=1,
            )

        y_range = config.y_range_uv if chrom.assigned_detector in ("UV", "PDA") else config.y_range_elsd
        fig.update_yaxes(
            title_text=unit_str,
            showgrid=config.show_grid,
            gridcolor="rgba(200,200,200,0.45)",
            range=list(y_range) if y_range is not None else None,
            row=idx,
            col=1,
        )
        fig.update_xaxes(
            showgrid=config.show_grid,
            gridcolor="rgba(200,200,200,0.45)",
            range=list(config.x_range) if config.x_range is not None else None,
            row=idx,
            col=1,
        )

    fig.update_xaxes(title_text=config.x_axis_label, row=n_rows, col=1)
    fig.update_layout(
        template="plotly_white",
        title=dict(text=config.title, x=0.01, xanchor="left") if config.title else None,
        height=fig_height,
        width=config.width_px,
        font=dict(family="Arial, sans-serif", size=config.font_size, color="#1F2937"),
        hovermode="closest",
        margin=dict(l=70, r=40, t=55, b=55),
    )
    _apply_legend_layout(fig, config.legend_position, config.show_legend, config.font_size)
    return fig


def create_alignment_comparison_figure(
    chromatograms: Sequence[Chromatogram],
    proposed_estimates: Optional[Sequence[AlignmentEstimate]] = None,
    region_start_min: Optional[float] = None,
    region_end_min: Optional[float] = None,
    line_width: float = 1.6,
) -> go.Figure:
    """Create a 2-panel Before vs. After alignment comparison plot."""
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.10,
        subplot_titles=(
            "Before Alignment (Original Retention Times)",
            "After Alignment (Corrected Retention Times)",
        ),
    )

    shift_lookup = {}
    if proposed_estimates is not None:
        shift_lookup = {e.trace_id: e.estimated_shift_min for e in proposed_estimates}

    for chrom in chromatograms:
        # Row 1: Original retention time
        fig.add_trace(
            go.Scatter(
                x=chrom.raw_retention_time,
                y=chrom.raw_intensity,
                mode="lines",
                name=chrom.display_label,
                legendgroup=chrom.trace_id,
                showlegend=True,
                line=dict(color=chrom.color, width=line_width),
                hovertemplate=(
                    f"<b>{chrom.sample_display_name} (Original)</b><br>"
                    "RT: %{x:.4f} min<br>Intensity: %{y:.5g}<extra></extra>"
                ),
            ),
            row=1,
            col=1,
        )

        # Row 2: Proposed shift (if previewing) or currently applied shift
        shift_val = shift_lookup.get(chrom.trace_id, chrom.time_shift_min if chrom.is_aligned else 0.0)
        corrected_rt = chrom.raw_retention_time + float(shift_val)
        fig.add_trace(
            go.Scatter(
                x=corrected_rt,
                y=chrom.raw_intensity,
                mode="lines",
                name=f"{chrom.display_label} ({shift_val:+.3f} min)",
                legendgroup=chrom.trace_id,
                showlegend=False,
                line=dict(color=chrom.color, width=line_width),
                hovertemplate=(
                    f"<b>{chrom.sample_display_name} (Aligned {shift_val:+.4f} min)</b><br>"
                    "Corrected RT: %{x:.4f} min<br>Intensity: %{y:.5g}<extra></extra>"
                ),
            ),
            row=2,
            col=1,
        )

    if region_start_min is not None and region_end_min is not None and region_end_min > region_start_min:
        for r in (1, 2):
            fig.add_vrect(
                x0=region_start_min,
                x1=region_end_min,
                fillcolor="rgba(0, 114, 178, 0.08)",
                line_width=1,
                line_dash="dash",
                line_color="rgba(0, 114, 178, 0.5)",
                row=r,
                col=1,
            )

    fig.update_yaxes(title_text="Intensity", row=1, col=1)
    fig.update_yaxes(title_text="Intensity", row=2, col=1)
    fig.update_xaxes(title_text="Retention Time (min)", row=2, col=1)
    fig.update_layout(
        template="plotly_white",
        height=620,
        font=dict(family="Arial, sans-serif", size=12),
        margin=dict(l=70, r=40, t=50, b=50),
    )
    return fig
