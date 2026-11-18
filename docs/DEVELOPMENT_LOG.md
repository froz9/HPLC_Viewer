# HPLC Viewer — Development Log & Architectural Decisions

## Phase 1: Empirical Data Format Investigation (`raw_data/*.arw`)

- **Action**: Inspected all 64 Waters Empower 3 `.arw` files in `raw_data/` (`HPLC_RAW_EXPORT1264.arw` to `1359.arw`).
- **Key Findings**:
  1. **Encoding**: Pure 7-bit ASCII across all 64 files.
  2. **Line Endings**: Strictly Carriage Return (`\r`, `0x0D`) with zero Line Feed (`\n`, `0x0A`) characters. Standard `BytesIO.readline()` or `io.StringIO` without newline normalization fails on these files by reading the entire file as a single line.
  3. **Two `.arw` Export Layouts**:
     - **2D Single-Channel Exports** (`ELSD Signal`, `2998 Ch1 254nm@1.2nm`, `2998 Ch2 365nm@1.2nm`, `2998 Ch2 280nm@1.2nm`): 2 metadata rows (`"SampleName"\t"Injection"\t"Channel"\t"Det. Units"` + values) followed immediately by 2-column tab-delimited numeric rows (`time_min \t intensity`).
     - **3D PDA Spectral Matrix Exports** (`Channel = "2998"`, ~16–53.7 MB): 2 metadata rows, followed by a `Wavelength\t 190.1094\t...` row (170 to 261 wavelengths) and a `Time` row, followed by `(1 + N_wavelengths)` tab-separated columns.
  4. **7-Significant-Figure Rounding in Empower**: At retention times $t \ge 10\text{ min}$, Empower rounds exported timestamps to 5 decimal places (`31.99667`, `31.99833`, `32`), causing adjacent $\Delta t$ steps to alternate between `0.00166` and `0.00167` min ($\sim 0.42\%$ relative jitter around the true $1/600\text{ min} = 10\text{ Hz}$ interval).
  5. **Duplicate `SampleName` Across Injections**: `Blank_Medium` (Injection 1) appears in two separate file groups (`1270–1273` and `1356–1359`).
- **Deliverable**: Documented in [`docs/arw_format_report.md`](arw_format_report.md).

---

## Phase 2: Data Import Engine & Validation (`src/data_import/`)

- **Architectural Decisions**:
  1. **Byte-Level Line Normalization**: All file readers normalize `\r\n` and `\r` to `\n` in memory before header scanning or C-engine tabular parsing.
  2. **Variable Header Scanning**: Rather than hardcoding a 2-line offset, `arw_reader.py` dynamically scans header rows until it detects either the `Wavelength` 3D header or the start of numeric 2-column observations, supporting both horizontal tabular headers and vertical key-value headers.
  3. **Fast Vectorized 3D PDA Parsing**: For 50+ MB 3D PDA files, `arw_reader.py` locates the byte offset of the numeric matrix and passes the normalized byte buffer directly to `pandas.read_csv(..., engine="c", dtype=np.float64)`, parsing 5 million floats in `< 0.3 s`. It stores an immutable `PDA3DSpectrum` and extracts default 2D channels (`254 nm`, `365 nm`) while allowing users to extract any custom wavelength on demand in the UI.
  4. **Immutable Raw Arrays**: `Chromatogram` sets `flags.writeable = False` on `_raw_retention_time` and `_raw_intensity` NumPy arrays so raw experimental observations can never be accidentally mutated in-place.
  5. **Scientifically Calibrated Validation**: `file_validation.py` checks all 9 required integrity conditions (missing RT, missing intensity, non-numeric tokens, duplicate timestamps, unsorted retention times, negative intensities, unequal sampling intervals, empty files, corrupted files) while distinguishing Empower's 7-significant-figure rounding jitter from true sampling interval anomalies.

---

## Phase 3: Interactive Visualization & Color Management (`src/visualization/`)

- **Architectural Decisions**:
  1. **Strict Separation of UV/PDA and ELSD Axes**: In Overlay mode, `create_chromatogram_figure()` automatically separates absorbance traces (`UV` / `PDA` in `AU`) and light-scattering traces (`ELSD` in `LSU`) into two vertically stacked, X-axis-synchronized subplots so incompatible physical units are never combined on a shared Y-axis.
  2. **Three Visualization Modes**: Implemented `Overlay` (shared axes per detector family), `Stacked` (vertically separated subplot per chromatogram with adjustable spacing and independent Y-scales), and `Individual` (single trace view).
  3. **Consistent Color Management**: `color_manager.py` provides 5 scientific palettes (defaulting to colorblind-safe Okabe-Ito) and groups colors by sample name by default so a sample's UV and ELSD traces share the same color across panels.

---

## Phase 4: Chromatogram Alignment (`src/data_processing/alignment.py`)

- **Architectural Decisions**:
  1. **Non-Destructive Shift Representation**: Alignment stores only `time_shift_min`, `alignment_method`, and `alignment_reference_id` on `Chromatogram`. `corrected_retention_time` is computed on the fly as `raw_retention_time + time_shift_min`.
  2. **Cross-Correlation with Baseline Detrending & Sub-Sample Refinement**: `estimate_cross_correlation_shift()` interpolates the reference and target signals onto a common uniform grid within the user-selected `[region_start_min, region_end_min]` window, detrends slow baseline drift on an internal working copy (never modifying raw intensity), computes normalized cross-correlation bounded by `[-max_shift_min, +max_shift_min]`, and applies 3-point parabolic interpolation around the correlation peak for sub-sample precision.
  3. **Preview Before Commit**: Users can preview proposed shifts in a dedicated 2-panel Before/After comparison plot before accepting or rejecting the alignment.

---

## Phase 5 & 6: Export, Session Management, Testing & Packaging

- **Architectural Decisions**:
  1. **Dual Static Figure Export Engine**: `figure_export.py` uses Plotly + Kaleido for high-DPI PNG, SVG, and PDF generation, with an automatic Matplotlib fallback if a target machine lacks a browser/Kaleido runtime.
  2. **Lightweight Session JSON**: `session_manager.py` serializes sample display names, colors, detector assignments, visibility, plot configuration, and alignment shifts without embedding raw numerical arrays.
  3. **Automated Test Suite**: Created 25+ automated tests across `tests/test_arw_reader.py`, `tests/test_csv_reader.py`, `tests/test_validation.py`, `tests/test_alignment.py`, and `tests/test_visualization.py`.
