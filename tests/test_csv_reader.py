"""Tests for CSV and TXT chromatographic readers across delimiters and formats."""

import numpy as np
import pytest

from src.data_import.csv_reader import DelimitedImportError, read_csv_chromatogram
from src.data_import.txt_reader import read_txt_chromatogram


def test_read_comma_delimited_csv_with_metadata() -> None:
    """Test CSV import with metadata header lines and column labels."""
    csv_bytes = (
        b"# SampleName: Standard_Mix_A\n"
        b"# Injection: 2\n"
        b"# Detector: UV\n"
        b"# Wavelength: 254 nm\n"
        b"# Units: AU\n"
        b"Retention Time (min),2998 Ch1 254nm\n"
        b"0.0,0.0012\n"
        b"0.1,0.0154\n"
        b"0.2,0.0450\n"
        b"0.3,0.0120\n"
    )
    res = read_csv_chromatogram(csv_bytes, file_name="std_mix.csv")
    assert res.is_valid
    assert res.delimiter == ","
    chrom = res.chromatograms[0]
    assert chrom.original_sample_name == "Standard_Mix_A"
    assert chrom.injection_id == "2"
    assert chrom.assigned_detector == "UV"
    assert chrom.wavelength_nm == 254.0
    assert chrom.units == "AU"
    assert chrom.n_points == 4
    assert np.allclose(chrom.raw_retention_time, [0.0, 0.1, 0.2, 0.3])
    assert np.allclose(chrom.raw_intensity, [0.0012, 0.0154, 0.0450, 0.0120])


def test_read_semicolon_delimited_csv() -> None:
    """Verify automatic detection of semicolon delimiters."""
    csv_bytes = (
        b"SampleName;Caffeic_Acid\n"
        b"Channel;ELSD Signal\n"
        b"Units;LSU\n"
        b"Time;Signal\n"
        b"0.0;-0.05\n"
        b"0.5;12.4\n"
        b"1.0;45.8\n"
    )
    res = read_csv_chromatogram(csv_bytes, file_name="semi.csv")
    assert res.is_valid
    assert res.delimiter == ";"
    chrom = res.chromatograms[0]
    assert chrom.original_sample_name == "Caffeic_Acid"
    assert chrom.assigned_detector == "ELSD"
    assert chrom.units == "LSU"
    assert chrom.n_points == 3


def test_read_tab_delimited_txt() -> None:
    """Verify TXT reader on tab-delimited files with CR/LF line endings."""
    txt_bytes = (
        b"SampleName\tExtract_01\r\n"
        b"Channel\tPDA 365nm\r\n"
        b"Units\tAU\r\n"
        b"0.0\t0.000\r\n"
        b"0.2\t0.125\r\n"
        b"0.4\t0.040\r\n"
    )
    res = read_txt_chromatogram(txt_bytes, file_name="extract.txt")
    assert res.is_valid
    assert res.delimiter == "\t"
    chrom = res.chromatograms[0]
    assert chrom.original_sample_name == "Extract_01"
    assert chrom.assigned_detector == "PDA"
    assert chrom.wavelength_nm == 365.0
    assert chrom.n_points == 3


def test_malformed_csv_raises_or_reports_error() -> None:
    """Verify clear error reporting on non-numeric or empty CSV files."""
    bad_csv = b"HeaderA,HeaderB\nfoo,bar\nbaz,qux\n"
    res = read_csv_chromatogram(bad_csv, file_name="bad.csv")
    assert not res.is_valid
    assert len(res.errors) == 1

    with pytest.raises(DelimitedImportError):
        read_csv_chromatogram(bad_csv, file_name="bad.csv", raise_on_error=True)
