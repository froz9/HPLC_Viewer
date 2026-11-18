"""Tests for manual and automatic cross-correlation chromatogram alignment."""

from pathlib import Path
import numpy as np
import pytest

from src.data_import.arw_reader import read_arw
from src.data_processing.alignment import (
    apply_alignment_estimates,
    apply_manual_shifts,
    estimate_batch_alignment,
    estimate_cross_correlation_shift,
    reset_all_alignments,
)
from src.data_processing.chromatogram import Chromatogram


RAW_DATA_DIR = Path(__file__).resolve().parents[1] / "raw_data"


def _make_synthetic_chromatogram(
    name: str,
    shift_min: float = 0.0,
    noise_std: float = 0.0,
    baseline_offset: float = 0.0,
    asymmetric: bool = False,
    seed: int = 42,
) -> Chromatogram:
    """Generate a synthetic multi-peak chromatogram where peaks are shifted by `shift_min`."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 20.0, 4001)  # dt = 0.005 min
    # Evaluate underlying profile at (t - shift_min) so features appear at t_peak + shift_min
    tau = t - shift_min
    y = (
        1.2 * np.exp(-0.5 * ((tau - 5.0) / 0.12) ** 2)
        + 0.7 * np.exp(-0.5 * ((tau - 9.5) / 0.18) ** 2)
        + 1.8 * np.exp(-0.5 * ((tau - 14.2) / 0.15) ** 2)
    )
    if asymmetric:
        # Add tailing component
        y += 0.4 * np.where(tau > 5.0, np.exp(-(tau - 5.0) / 0.35), 0.0)

    y = y + baseline_offset
    if noise_std > 0.0:
        y = y + rng.normal(0.0, noise_std, size=t.size)

    return Chromatogram(
        file_name=f"{name}.csv",
        retention_time=t,
        intensity=y,
        sample_name=name,
        detector="UV",
        channel="UV 254nm",
        units="AU",
        wavelength_nm=254.0,
    )


@pytest.mark.parametrize("true_delay_min", [0.15, -0.22, 0.037, -0.45, 0.0])
def test_cross_correlation_recovers_known_shifts_clean_and_noisy(true_delay_min: float) -> None:
    """Verify automatic alignment estimates the corrective shift (-true_delay_min) within tolerance."""
    ref = _make_synthetic_chromatogram("Ref", shift_min=0.0, noise_std=0.005, seed=10)
    # Target peaks elute at t_ref + true_delay_min, so corrective shift is -true_delay_min
    target = _make_synthetic_chromatogram(
        "Target",
        shift_min=true_delay_min,
        noise_std=0.02,
        baseline_offset=-0.35,
        asymmetric=True,
        seed=20,
    )

    est = estimate_cross_correlation_shift(
        reference=ref,
        target=target,
        region_start_min=2.0,
        region_end_min=18.0,
        max_shift_min=1.0,
    )
    expected_correction = -true_delay_min
    assert abs(est.estimated_shift_min - expected_correction) < 0.01, (
        f"Expected {expected_correction:.4f}, got {est.estimated_shift_min:.4f}"
    )
    assert est.correlation_score > 0.85


def test_alignment_region_selection_ignores_out_of_window_artifacts() -> None:
    """Verify that specifying an alignment window isolates the target peak from solvent front artifacts."""
    t = np.linspace(0.0, 15.0, 3001)
    # Reference has analytical peak at 8.0 min and a solvent front spike at 1.0 min
    y_ref = np.exp(-0.5 * ((t - 8.0) / 0.12) ** 2) + 5.0 * np.exp(-0.5 * ((t - 1.0) / 0.08) ** 2)
    # Target has analytical peak shifted to 8.25 min (+0.25 min delay), but solvent front still at 1.0 min
    y_tgt = np.exp(-0.5 * ((t - 8.25) / 0.12) ** 2) + 5.0 * np.exp(-0.5 * ((t - 1.0) / 0.08) ** 2)

    ref = Chromatogram("ref.csv", t, y_ref, sample_name="Ref", detector="UV")
    tgt = Chromatogram("tgt.csv", t, y_tgt, sample_name="Tgt", detector="UV")

    # Aligning specifically in the [6.0, 10.0] min region must recover -0.25 min
    est = estimate_cross_correlation_shift(
        reference=ref,
        target=tgt,
        region_start_min=6.0,
        region_end_min=10.0,
        max_shift_min=1.0,
    )
    assert abs(est.estimated_shift_min - (-0.25)) < 0.005


def test_real_empower_arw_alignment_and_immutability() -> None:
    """Verify shift estimation and raw data immutability using real Empower .arw files."""
    ref_res = read_arw(RAW_DATA_DIR / "HPLC_RAW_EXPORT1278.arw")  # 01006-18-1 UV 254nm
    assert ref_res.is_valid
    ref = ref_res.chromatograms[0]

    # Create a shifted version of the real chromatogram (+0.12 min delay)
    known_delay = 0.12
    shifted_real = Chromatogram(
        file_name="shifted_1278.arw",
        retention_time=ref.raw_retention_time,
        intensity=np.interp(
            ref.raw_retention_time - known_delay,
            ref.raw_retention_time,
            ref.raw_intensity,
        ),
        sample_name="01006-18-1_shifted",
        detector="UV",
        channel=ref.channel,
        units=ref.units,
        wavelength_nm=ref.wavelength_nm,
    )

    orig_rt_copy = shifted_real.raw_retention_time.copy()
    orig_int_copy = shifted_real.raw_intensity.copy()

    estimates = estimate_batch_alignment(
        reference=ref,
        targets=[ref, shifted_real],
        region_start_min=2.0,
        region_end_min=28.0,
        max_shift_min=1.0,
    )
    assert len(estimates) == 2
    assert estimates[0].is_reference
    assert estimates[0].estimated_shift_min == 0.0
    assert abs(estimates[1].estimated_shift_min - (-known_delay)) < 0.005

    # Apply alignment and check that raw data is unchanged while corrected RT is shifted
    apply_alignment_estimates([ref, shifted_real], estimates)
    assert shifted_real.is_aligned
    assert np.array_equal(shifted_real.raw_retention_time, orig_rt_copy)
    assert np.array_equal(shifted_real.raw_intensity, orig_int_copy)
    assert np.allclose(
        shifted_real.corrected_retention_time,
        orig_rt_copy + estimates[1].estimated_shift_min,
    )

    # Manual alignment override and reset
    apply_manual_shifts([ref, shifted_real], {shifted_real.trace_id: -0.25}, reference_trace_id=ref.trace_id)
    assert np.isclose(shifted_real.time_shift_min, -0.25)
    assert shifted_real.alignment_method == "Manual"

    reset_all_alignments([ref, shifted_real])
    assert not shifted_real.is_aligned
    assert shifted_real.time_shift_min == 0.0
    assert np.array_equal(shifted_real.active_retention_time, orig_rt_copy)
