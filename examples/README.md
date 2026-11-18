# Example Chromatographic Files

This directory contains example `.csv` and `.txt` chromatographic files demonstrating the supported delimited text formats in addition to the Waters Empower 3 `.arw` files located in `../raw_data/`.

## Files Included

1. **`example_chromatogram.csv`**:
   - Comma-delimited UV chromatogram (`254 nm`, `AU`) with optional `# Key: Value` header metadata lines.
2. **`example_chromatogram.txt`**:
   - Tab-delimited ELSD chromatogram (`LSU`) with key-value metadata rows.
3. **`../raw_data/*.arw`**:
   - 64 real Waters Empower 3 ASCII exports across 16 injections (ELSD 2D, UV 254 nm 2D, UV 365/280 nm 2D, and 3D PDA spectral matrices).
   - See `../docs/arw_format_report.md` for the complete technical format analysis.
