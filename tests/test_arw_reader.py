"""Tests for the Waters Empower 3 .arw parser on real and synthetic files."""

from pathlib import Path
import numpy as np
import pytest

from src.data_import.arw_reader import ARWImportError, read_arw


RAW_DATA_DIR = Path(__file__).resolve().parents[1] / "raw_data"


def test_read_real_2d_elsd_arw() -> None:
    """Verify exact retention-time, intensity, and metadata extraction on a real ELSD .arw file."""
    arw_path = RAW_DATA_DIR / "HPLC_RAW_EXPORT1264.arw"
    result = read_arw(arw_path)
    assert result.is_valid
    assert not result.is_3d_pda
    assert len(result.chromatograms) == 1

    chrom = result.chromatograms[0]
    assert chrom.file_name == "HPLC_RAW_EXPORT1264.arw"
    assert chrom.original_sample_name == "Diox.MeOH"
    assert chrom.injection_id == "1"
    assert chrom.assigned_detector == "ELSD"
    assert chrom.channel == "ELSD Signal"
    assert chrom.units == "LSU"
    assert chrom.wavelength_nm is None
    assert chrom.n_points == 6001
    assert np.isclose(chrom.rt_min, 0.0)
    assert np.isclose(chrom.rt_max, 10.0)
    assert np.isclose(chrom.raw_intensity[0], 0.5277824)
    assert np.isclose(chrom.raw_intensity[-1], -0.009641647)

    # Ensure raw arrays are immutable
    with pytest.raises(ValueError):
        chrom.raw_retention_time[0] = 999.0
    with pytest.raises(ValueError):
        chrom.raw_intensity[0] = 999.0


def test_read_real_2d_uv_arw_channels() -> None:
    """Verify UV detector and wavelength extraction across 254nm, 365nm, and 280nm real exports."""
    res_254 = read_arw(RAW_DATA_DIR / "HPLC_RAW_EXPORT1266.arw")
    assert res_254.is_valid
    c_254 = res_254.chromatograms[0]
    assert c_254.assigned_detector == "UV"
    assert c_254.wavelength_nm == 254.0
    assert c_254.bandwidth_nm == 1.2
    assert c_254.units == "AU"
    assert c_254.n_points == 6001
    assert np.isclose(c_254.raw_intensity[0], -8.24e-5)

    res_365 = read_arw(RAW_DATA_DIR / "HPLC_RAW_EXPORT1267.arw")
    c_365 = res_365.chromatograms[0]
    assert c_365.assigned_detector == "UV"
    assert c_365.wavelength_nm == 365.0
    assert c_365.bandwidth_nm == 1.2

    res_280 = read_arw(RAW_DATA_DIR / "HPLC_RAW_EXPORT1351.arw")
    c_280 = res_280.chromatograms[0]
    assert c_280.original_sample_name == "Shutdown"
    assert c_280.assigned_detector == "UV"
    assert c_280.wavelength_nm == 280.0
    assert c_280.n_points == 12061


def test_read_real_3d_pda_arw() -> None:
    """Verify 3D PDA spectral matrix parsing and wavelength extraction on real 16MB file."""
    arw_path = RAW_DATA_DIR / "HPLC_RAW_EXPORT1265.arw"
    result = read_arw(arw_path, target_pda_wavelengths_nm=(254.0, 365.0))
    assert result.is_valid
    assert result.is_3d_pda
    assert result.pda_spectrum is not None
    assert result.pda_spectrum.n_times == 6001
    assert result.pda_spectrum.n_wavelengths == 261
    assert np.isclose(result.pda_spectrum.wl_min, 190.1094)
    assert np.isclose(result.pda_spectrum.wl_max, 500.2398)
    assert len(result.chromatograms) == 2
    assert result.chromatograms[0].assigned_detector == "PDA"
    assert result.chromatograms[0].is_from_3d_pda

    # Extract an arbitrary custom wavelength from the 3D PDA spectrum
    custom_trace = result.pda_spectrum.extract_wavelength(210.0)
    assert custom_trace.assigned_detector == "PDA"
    assert abs(custom_trace.wavelength_nm - 210.1) < 1.0
    assert custom_trace.n_points == 6001


def test_all_48_real_2d_arw_files_parse_cleanly() -> None:
    """Validate that all 48 2D .arw files in raw_data/ parse with zero errors."""
    all_files = sorted(RAW_DATA_DIR.glob("*.arw"))
    assert len(all_files) == 64

    count_2d = 0
    for fpath in all_files:
        # Skip the 16 large 3D files in the bulk loop to keep test suite fast (<2s)
        if fpath.stat().st_size > 5_000_000:
            continue
        res = read_arw(fpath)
        assert res.is_valid, f"Failed on {fpath.name}: {res.errors}"
        chrom = res.chromatograms[0]
        assert chrom.validation_report.status == "OK"
        assert not chrom.validation_report.unequal_sampling_interval
        count_2d += 1

    assert count_2d == 48


def test_variable_header_lengths_and_malformed_arw() -> None:
    """Verify handling of different metadata header lengths, empty files, and corrupted inputs."""
    # 4-line key-value metadata header before numeric section
    custom_arw = (
        b'"SampleName"\t"CustomSample"\r'
        b'"Injection"\t"3"\r'
        b'"Channel"\t"UV 214nm"\r'
        b'"Det. Units"\t"AU"\r'
        b"0.0\t0.001\r0.5\t0.025\r1.0\t0.010\r"
    )
    res = read_arw(custom_arw, file_name="custom.arw")
    assert res.is_valid
    chrom = res.chromatograms[0]
    assert chrom.original_sample_name == "CustomSample"
    assert chrom.injection_id == "3"
    assert chrom.assigned_detector == "UV"
    assert chrom.wavelength_nm == 214.0
    assert chrom.n_points == 3

    # Empty file
    empty_res = read_arw(b"", file_name="empty.arw")
    assert not empty_res.is_valid
    assert "empty" in empty_res.errors[0].lower()

    # Binary corrupted file
    with pytest.raises(ARWImportError):
        read_arw(b"\x00\x01\x02\x03\xff", file_name="corrupt.arw", raise_on_error=True)

    # Ambiguous detector channel is marked 'Other / Unknown' and triggers warning
    ambig_arw = (
        b'"SampleName"\t"Injection"\t"Channel"\t"Det. Units"\r'
        b'"UnkSample"\t1\t"Aux_Channel_A"\t"mV"\r'
        b"0.0\t1.0\r0.1\t1.5\r0.2\t1.2\r"
    )
    ambig_res = read_arw(ambig_arw, file_name="ambig.arw")
    assert ambig_res.is_valid
    ambig_chrom = ambig_res.chromatograms[0]
    assert ambig_chrom.assigned_detector == "Other / Unknown"
    assert ambig_chrom.validation_report.detector_ambiguous
