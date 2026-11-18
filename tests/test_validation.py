"""Tests for chromatographic data validation (Section 4.5)."""

import numpy as np

from src.data_import.arw_reader import read_arw
from src.data_import.file_validation import validate_chromatogram_arrays


def test_validation_preserves_negative_and_duplicate_points() -> None:
    """Verify negative intensities and duplicate timestamps are flagged but never deleted."""
    rt = [0.0, 0.1, 0.1, 0.2, 0.3]
    intensity = [-0.05, 0.10, 0.12, -0.02, 0.50]

    report = validate_chromatogram_arrays(
        file_name="dup_neg.arw",
        retention_time=rt,
        intensity=intensity,
        channel="UV 254nm",
    )
    assert report.is_valid
    assert report.status == "Warning"
    assert report.negative_intensity_count == 2
    assert report.duplicate_time_count == 1
    assert report.unequal_sampling_interval


def test_validation_detects_unsorted_retention_times() -> None:
    """Verify unsorted retention times are detected and reported."""
    rt = [0.0, 0.2, 0.15, 0.3, 0.4]
    intensity = [0.1, 0.2, 0.3, 0.2, 0.1]

    report = validate_chromatogram_arrays(
        file_name="unsorted.arw",
        retention_time=rt,
        intensity=intensity,
    )
    assert report.unsorted_time_count == 1
    assert report.unequal_sampling_interval
    assert report.status == "Warning"
    assert any("Unsorted" in w for w in report.warnings)


def test_validation_detects_unequal_sampling_intervals_vs_empower_rounding() -> None:
    """Distinguish Empower 7-sig-fig ASCII rounding from genuinely unequal intervals."""
    # Empower 7-significant-digit timestamps around 31.99 min (step 1/600 min)
    rt_empower = [31.99333, 31.995, 31.99667, 31.99833, 32.0]
    inten = [0.01, 0.01, 0.01, 0.01, 0.01]
    rep_emp = validate_chromatogram_arrays("emp.arw", rt_empower, inten)
    assert rep_emp.status == "OK"
    assert not rep_emp.unequal_sampling_interval
    assert rep_emp.ascii_rounding_jitter

    # Genuinely unequal sampling intervals
    rt_unequal = [0.0, 0.1, 0.2, 0.5, 0.6]
    rep_uneq = validate_chromatogram_arrays("uneq.arw", rt_unequal, inten)
    assert rep_uneq.unequal_sampling_interval
    assert rep_uneq.status == "Warning"


def test_validation_missing_and_non_numeric_in_arw() -> None:
    """Verify missing RT, missing intensity, and non-numeric lines are counted and reported."""
    raw = (
        b'"SampleName"\t"Injection"\t"Channel"\t"Det. Units"\r'
        b'"TestSample"\t1\t"ELSD Signal"\t"LSU"\r'
        b"0.0\t1.0\r"
        b"\t2.0\r"          # missing RT
        b"0.2\t\r"          # missing intensity
        b"0.3\tbad_val\r"   # non-numeric
        b"0.4\t3.5\r"
    )
    res = read_arw(raw, file_name="partial.arw")
    assert res.is_valid
    chrom = res.chromatograms[0]
    assert chrom.n_points == 2
    rep = chrom.validation_report
    assert rep.missing_rt_count == 1
    assert rep.missing_intensity_count == 1
    assert rep.non_numeric_count == 1
    assert rep.status == "Warning"


def test_validation_empty_and_mismatched_arrays() -> None:
    """Verify empty and mismatched arrays return status='Error'."""
    rep_empty = validate_chromatogram_arrays("empty.csv", [], [])
    assert not rep_empty.is_valid
    assert rep_empty.is_empty
    assert rep_empty.status == "Error"

    rep_mismatch = validate_chromatogram_arrays("mismatch.csv", [0.0, 1.0], [1.0])
    assert not rep_mismatch.is_valid
    assert rep_mismatch.is_corrupted
    assert rep_mismatch.status == "Error"
