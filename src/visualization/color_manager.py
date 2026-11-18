"""Color palette management for scientific chromatogram visualization."""

from __future__ import annotations

from typing import Dict, List, Sequence
from src.data_processing.chromatogram import Chromatogram


COLOR_PALETTES: Dict[str, List[str]] = {
    "Okabe-Ito (Colorblind Safe)": [
        "#0072B2",  # Blue
        "#D55E00",  # Vermilion
        "#009E73",  # Bluish green
        "#CC79A7",  # Reddish purple
        "#E69F00",  # Orange
        "#56B4E9",  # Sky blue
        "#332288",  # Indigo
        "#F0E442",  # Yellow
        "#882255",  # Wine
        "#44AA99",  # Teal
    ],
    "Nature / NPG": [
        "#E64B35",
        "#4DBBD5",
        "#00A087",
        "#3C5488",
        "#F39B7F",
        "#8491B4",
        "#91D1C2",
        "#DC0000",
        "#7E6148",
        "#B09C85",
    ],
    "Tableau 10": [
        "#1F77B4",
        "#FF7F0E",
        "#2CA02C",
        "#D62728",
        "#9467BD",
        "#8C564B",
        "#E377C2",
        "#7F7F7F",
        "#BCBD22",
        "#17BECF",
    ],
    "Viridis Discrete": [
        "#440154",
        "#482878",
        "#3E4989",
        "#31688E",
        "#26828E",
        "#1F9E89",
        "#35B779",
        "#6DCD59",
        "#B4DE2C",
        "#FDE725",
    ],
    "Classic Analytical": [
        "#0B3C5D",
        "#B82601",
        "#1C6B2A",
        "#6B2D5C",
        "#D9822B",
        "#2A7B88",
        "#4A4A4A",
        "#8C1D40",
        "#3D5A80",
        "#556B2F",
    ],
}

DEFAULT_PALETTE_NAME = "Okabe-Ito (Colorblind Safe)"


def get_palette_colors(palette_name: str = DEFAULT_PALETTE_NAME) -> List[str]:
    """Return the hex color list for a named palette."""
    return list(COLOR_PALETTES.get(palette_name, COLOR_PALETTES[DEFAULT_PALETTE_NAME]))


def assign_colors_to_chromatograms(
    chromatograms: Sequence[Chromatogram],
    palette_name: str = DEFAULT_PALETTE_NAME,
    group_by_sample_name: bool = True,
) -> None:
    """Assign palette colors consistently across chromatograms.

    When `group_by_sample_name` is True, channels belonging to the same sample
    receive the same color so UV and ELSD panels share consistent sample colors.
    """
    colors = get_palette_colors(palette_name)
    if not colors:
        return

    if group_by_sample_name:
        sample_color_map: Dict[str, str] = {}
        for chrom in chromatograms:
            key = chrom.sample_display_name
            if key not in sample_color_map:
                sample_color_map[key] = colors[len(sample_color_map) % len(colors)]
            chrom.color = sample_color_map[key]
    else:
        for idx, chrom in enumerate(chromatograms):
            chrom.color = colors[idx % len(colors)]


def reset_colors_to_palette(
    chromatograms: Sequence[Chromatogram],
    palette_name: str = DEFAULT_PALETTE_NAME,
    group_by_sample_name: bool = True,
) -> None:
    """Reset all chromatogram colors to the specified palette."""
    assign_colors_to_chromatograms(
        chromatograms=chromatograms,
        palette_name=palette_name,
        group_by_sample_name=group_by_sample_name,
    )
