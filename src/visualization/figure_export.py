"""Publication-quality figure and tabular data export utilities."""

from __future__ import annotations

import io
from typing import Sequence
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go

from src.data_processing.chromatogram import Chromatogram
from src.utils.logging_config import get_logger

matplotlib.use("Agg")
logger = get_logger(__name__)

SUPPORTED_FIGURE_FORMATS = ("png", "svg", "pdf", "html")


def export_figure_bytes(
    fig: go.Figure,
    fmt: str = "png",
    width_px: int = 1200,
    height_px: int = 750,
    scale: int = 2,
) -> bytes:
    """Export a Plotly chromatogram figure to PNG, SVG, PDF, or interactive HTML bytes.

    Uses Kaleido for static formats with a Matplotlib fallback if Kaleido/browser
    rendering is unavailable on the host system.
    """
    fmt_clean = fmt.lower().strip().lstrip(".")
    if fmt_clean not in SUPPORTED_FIGURE_FORMATS:
        raise ValueError(f"Unsupported export format '{fmt}'. Choose from {SUPPORTED_FIGURE_FORMATS}.")

    if fmt_clean == "html":
        html_str = fig.to_html(
            include_plotlyjs=True,
            full_html=True,
            config={"scrollZoom": True, "displaylogo": False},
        )
        return html_str.encode("utf-8")

    try:
        return fig.to_image(
            format=fmt_clean,
            width=int(width_px),
            height=int(height_px),
            scale=int(scale),
        )
    except Exception as exc:
        logger.warning(
            "Kaleido static export failed (%s); falling back to Matplotlib renderer.", exc
        )
        return _export_via_matplotlib_fallback(
            fig=fig,
            fmt=fmt_clean,
            width_px=width_px,
            height_px=height_px,
            scale=scale,
        )


def _export_via_matplotlib_fallback(
    fig: go.Figure,
    fmt: str,
    width_px: int,
    height_px: int,
    scale: int,
) -> bytes:
    """Render Plotly Scatter traces with Matplotlib if Kaleido is unavailable."""
    dpi = 100 * max(1, int(scale))
    w_in = max(4.0, width_px / 100.0)
    h_in = max(3.0, height_px / 100.0)

    mpl_fig, ax = plt.subplots(figsize=(w_in, h_in), dpi=dpi)
    for trace in fig.data:
        if getattr(trace, "visible", True) is False:
            continue
        x = getattr(trace, "x", None)
        y = getattr(trace, "y", None)
        if x is None or y is None:
            continue
        line_dict = getattr(trace, "line", {}) or {}
        color = getattr(line_dict, "color", None) or "#0072B2"
        width = getattr(line_dict, "width", 1.5) or 1.5
        label = getattr(trace, "name", "")
        show_leg = getattr(trace, "showlegend", True)
        ax.plot(
            x,
            y,
            color=color,
            linewidth=width,
            label=label if show_leg else "_nolegend_",
            alpha=getattr(trace, "opacity", 0.9) or 0.9,
        )

    title_obj = getattr(fig.layout, "title", None)
    title_text = getattr(title_obj, "text", "") if title_obj else ""
    if title_text:
        ax.set_title(title_text)
    xaxis_title = getattr(getattr(fig.layout, "xaxis", None), "title", None)
    ax.set_xlabel(getattr(xaxis_title, "text", None) or "Retention Time (min)")
    yaxis_title = getattr(getattr(fig.layout, "yaxis", None), "title", None)
    ax.set_ylabel(getattr(yaxis_title, "text", None) or "Signal Intensity")
    ax.grid(True, alpha=0.35)
    if getattr(fig.layout, "showlegend", True) and any(
        getattr(t, "showlegend", True) for t in fig.data
    ):
        ax.legend(loc="best", frameon=True)

    mpl_fig.tight_layout()
    buf = io.BytesIO()
    mpl_fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight")
    plt.close(mpl_fig)
    return buf.getvalue()


def export_chromatograms_to_csv_bytes(chromatograms: Sequence[Chromatogram]) -> bytes:
    """Export selected chromatograms to CSV bytes preserving raw and corrected retention times."""
    if not chromatograms:
        return b""
    frames = [c.to_export_dataframe() for c in chromatograms]
    combined = pd.concat(frames, ignore_index=True)
    buf = io.StringIO()
    combined.to_csv(buf, index=False, float_format="%.7g")
    return buf.getvalue().encode("utf-8")
