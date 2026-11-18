# HPLC Viewer — Local Analytical Chromatogram Visualization & Alignment

**HPLC Viewer** is a local scientific desktop application designed for analytical chemistry and natural-products research laboratories to import, inspect, compare, align, and export high-performance liquid chromatography (HPLC) chromatograms without requiring a Waters Empower installation or an internet connection.

Built specifically around **Waters Empower 3** exports with **UV / Photodiode Array (PDA, Waters 2998)** and **Evaporative Light Scattering (ELSD)** detectors, the application prioritizes data integrity, exact retention-time representation, and reproducible transformations.

---

## 1. Core Scientific & Engineering Principles

1. **Zero Modification of Source Files**: Original experimental `.arw`, `.csv`, and `.txt` files are read strictly in-memory and are never altered on disk.
2. **Immutable Raw Observations**: Underlying raw retention-time and intensity arrays are locked as read-only (`writeable = False`) NumPy arrays upon import.
3. **Strict Separation of Detector Units**: UV/PDA absorbance signals (`AU`) and ELSD light-scattering signals (`LSU`) are never combined on a shared Y-axis. When both detectors are active, they are rendered in vertically synchronized panels sharing the retention-time axis.
4. **Non-Destructive Retention-Time Alignment**: Alignment stores an additive shift $\Delta t$ (`corrected_retention_time = raw_retention_time + time_shift_min`) separately from `raw_retention_time`, allowing instant before/after comparison and full reversal.
5. **100% Local Operation**: All processing runs locally on `localhost` with telemetry disabled (`gatherUsageStats = false`).

---

## 2. Project Architecture

```text
HPLC_Viewer/
├── app.py                          # Main Streamlit application entry point
├── launch_hplc_viewer.bat          # One-click Windows launch script
├── requirements.txt                # Explicit Python package dependencies
├── pytest.ini                      # Pytest configuration
├── README.md                       # Documentation, usage, and packaging guide
├── .streamlit/
│   └── config.toml                 # Local-only server & scientific theme config
├── config/
│   ├── default_settings.json       # Default visualization, alignment & export settings
│   └── example_session.json        # Example saved session configuration
├── docs/
│   ├── arw_format_report.md        # Technical report on Waters Empower 3 .arw format
│   └── DEVELOPMENT_LOG.md          # Development log & architectural decisions
├── examples/
│   ├── README.md                   # Guide to example files
│   ├── example_chromatogram.csv    # Example comma-delimited UV chromatogram
│   └── example_chromatogram.txt    # Example tab-delimited ELSD chromatogram
├── raw_data/                       # 64 Waters Empower 3 .arw files (16 injections)
├── logs/                           # Local application logs
├── src/
│   ├── data_import/
│   │   ├── arw_reader.py           # Empower 3 2D & 3D PDA .arw parser
│   │   ├── csv_reader.py           # Comma/semicolon CSV parser & round-trip loader
│   │   ├── txt_reader.py           # Tab-delimited TXT parser
│   │   └── file_validation.py      # 9-point chromatographic data validator
│   ├── data_processing/
│   │   ├── chromatogram.py         # Immutable Chromatogram & PDA3DSpectrum models
│   │   ├── alignment.py            # Manual & cross-correlation RT alignment engine
│   │   └── normalization.py        # Non-destructive visual scaling & baseline envelope
│   ├── visualization/
│   │   ├── chromatogram_plot.py    # Interactive Plotly Overlay, Stacked & Individual plots
│   │   ├── color_manager.py        # Scientific colorblind-safe palettes
│   │   └── figure_export.py        # PNG, SVG, PDF (Kaleido + Matplotlib fallback), HTML & CSV export
│   ├── interface/
│   │   ├── sidebar.py              # Sidebar Sections A (Import), B (Samples), C (Vis), D (Alignment)
│   │   ├── sample_manager.py       # Sample selection, renaming, colors & 3D PDA wavelength extractor
│   │   └── import_interface.py     # Cached file ingestion & interactive editable sample table
│   └── utils/
│       ├── logging_config.py       # Structured file & console logging
│       └── session_manager.py      # JSON session state serialization & restoration
└── tests/
    ├── test_arw_reader.py          # Real & synthetic .arw import tests
    ├── test_csv_reader.py          # CSV & TXT reader tests
    ├── test_validation.py          # Data validation & integrity tests
    ├── test_alignment.py           # Cross-correlation & manual alignment tests
    └── test_visualization.py       # Plotting, color management & export tests
```

---

## 3. Installation Instructions

### Prerequisites
- **Python 3.11+** (tested on Python 3.11–3.13 on Windows).

### Setup
Open PowerShell or Command Prompt in the `HPLC_Viewer` folder:

```powershell
# Optional: create and activate a virtual environment
python -m venv .venv
.venv\Scripts\activate

# Install required dependencies
pip install -r requirements.txt
```

---

## 4. Launching the Application

### Option A: Windows Double-Click Launcher
Double-click **`launch_hplc_viewer.bat`** in Windows Explorer, or run:
```powershell
.\launch_hplc_viewer.bat
```

### Option B: Command Line
```powershell
python -m streamlit run app.py --server.address=localhost --browser.gatherUsageStats=false
```
The application will open automatically in your local browser at `http://localhost:8501`.

---

## 5. Usage Guide

### A. Importing Chromatograms (`Sidebar → A. Import`)
- **Drag and Drop**: Drag one or more `.arw`, `.csv`, or `.txt` files into the uploader.
- **Local `raw_data/` Workspace Loader**: Expand **Load from Local `raw_data/` Folder** to load the pre-configured 3-sample UV+ELSD example set with one click, or pick specific `.arw` files from `raw_data/`.
- **Interactive Sample Table**: In **Main Panel → 1. Summary of Imported Files**, inspect file metadata, sampling interval, retention-time range, and validation status. Double-click any cell in **Plot**, **Color**, **Sample Name**, or **Detector** to customize display settings without modifying the source file.

### B. Managing Samples & 3D PDA Wavelengths (`Sidebar → B. Samples`)
- Use **Select All** / **Deselect All** to toggle visibility in bulk.
- Choose from 5 scientific color palettes (**Okabe-Ito Colorblind Safe**, **Nature / NPG**, **Tableau 10**, **Viridis Discrete**, **Classic Analytical**) and toggle **Match colors by sample name** so UV and ELSD traces from the same sample share identical colors.
- When a **3D PDA `.arw` file** (`Channel = "2998"`) is imported, expand **🌈 Extract 3D PDA Wavelength** in the sidebar to extract any diode wavelength (`190.1–500.2 nm`) as an independent chromatographic trace.

### C. Visualization Modes (`Sidebar → C. Visualization`)
- **Overlay**: Plots multiple chromatograms on shared axes per detector family (Upper panel: UV/PDA in `AU`; Lower panel: ELSD in `LSU`).
- **Stacked**: Displays each chromatogram in its own vertically separated subplot panel with independent Y-axis scaling and adjustable vertical spacing.
- **Individual**: Inspects a single selected chromatogram channel in isolation.
- Filter active detectors (`UV`, `PDA`, `ELSD`) and individual UV/PDA wavelengths (e.g., `254 nm`, `365 nm`, `280 nm`), customize X/Y axis bounds, line thickness, opacity, and optional visual baseline or intensity normalization.

### D. Retention-Time Alignment (`Sidebar → D. Alignment`)
- **Automatic Alignment (Cross-Correlation)**:
  1. Select a **Reference Chromatogram**.
  2. Specify the **Region Start (min)** and **Region End (min)** window containing the diagnostic peak(s) and set **Max Permitted Shift (± min)**.
  3. Click **🔍 Preview Proposed Alignment** to inspect the estimated shifts, correlation scores, and the **Before vs. After Alignment** comparison plot in the main panel.
  4. Click **✅ Accept & Apply** to apply the shifts or **❌ Reject Preview** to discard them.
- **Manual Alignment**: Enter exact additive shifts $\Delta t$ (in minutes) per sample and click **✅ Apply Manual Shifts**.
- Click **↩️ Reset All Alignments** at any time to revert to original retention times.

### E. Exporting Figures, Data & Sessions (`Main Panel → 4. Export Controls`)
- **Publication Figures**: Customize title, axis labels, width/height (px), DPI scale (`1x`–`4x`), font size, and legend placement, and download as **PNG**, **SVG**, **PDF**, or **Interactive HTML**.
- **Aligned Data (CSV)**: Export tabular CSV files containing both `original_retention_time_min`, `original_intensity`, and `corrected_retention_time_min` (plus `applied_shift_min` and sample/detector metadata).
- **Session Configuration (`.json`)**: Save and reload all sample names, colors, detector assignments, alignment shifts, and plot settings without embedding raw experimental arrays.

---

## 6. Supported Formats & Empower `.arw` Technical Summary

Full technical details are documented in [`docs/arw_format_report.md`](docs/arw_format_report.md).

- **Waters Empower 3 ASCII Raw Data (`.arw`)**:
  - Supports **CR-only (`\r`)**, **CRLF (`\r\n`)**, and **LF (`\n`)** line endings.
  - Supports both **2D single-channel exports** (`ELSD Signal`, `2998 Ch1 254nm@1.2nm`, `2998 Ch2 365nm@1.2nm`, `2998 Ch2 280nm@1.2nm`) and **3D PDA spectral matrix exports** (`Channel = "2998"` with `Wavelength` and `Time` header rows).
  - Handles variable metadata header lengths (horizontal tabular headers or key-value rows).
- **Delimited Text Files (`.csv`, `.txt`)**:
  - Automatic delimiter detection (`,`, `;`, `\t`), optional metadata headers, and round-trip re-import of CSV files exported by HPLC Viewer.

---

## 7. Supported vs. Unsupported Features (Version 1.0)

| Category | Supported in Version 1.0 | Unsupported in Version 1.0 (Future Modules) |
| :--- | :--- | :--- |
| **File Formats** | Waters Empower 3 `.arw` (2D & 3D PDA), `.csv`, `.txt` | Proprietary binary Empower database files (`.raw` folders) |
| **Detectors** | UV, PDA (multi-wavelength & extracted), ELSD, Other/Unknown | MS (`m/z` mass spectra fragmentation viewers) |
| **Visualization** | Overlay, Stacked, Individual; synchronized UV & ELSD panels | 3D surface/contour heatmaps of PDA matrices |
| **Alignment** | Linear retention-time shifting (Manual & Cross-Correlation) | Nonlinear retention-time warping (COW / DTW) |
| **Signal Processing** | Optional visual baseline envelope & display normalization | Automated peak integration, peak picking & compound ID |
| **Export** | PNG, SVG, PDF, Interactive HTML, CSV (raw + aligned), JSON session | Direct write-back to Waters Empower database |

---

## 8. Running Automated Tests

Run the complete test suite with `pytest`:

```powershell
pytest -v
```

---

## 9. Troubleshooting

- **Large 3D PDA `.arw` files (>50 MB)**: The default upload limit is configured to `400 MB` in `.streamlit/config.toml`. Importing a 53 MB 3D PDA `.arw` file takes ~0.3 seconds and extracts `254 nm` and `365 nm` by default; use **Sidebar → B. Samples → Extract 3D PDA Wavelength** to extract additional wavelengths.
- **Static PNG/PDF Export**: Uses `kaleido` (with local Chrome/Edge) by default and automatically falls back to `matplotlib` if a browser runtime is unavailable.
- **Ambiguous Detector Warning**: If an imported `.csv` or `.txt` file does not specify `UV`, `PDA`, `ELSD`, `AU`, `LSU`, or a wavelength in `nm`, it is safely assigned `Other / Unknown` rather than guessing. You can assign its detector in the editable table or sidebar.

---

## 10. Packaging as a Standalone Windows Desktop Executable

To distribute **HPLC Viewer** to laboratories without Python installed, you can bundle it using **PyInstaller** (already available in the environment):

1. Create a wrapper entry script `run_desktop.py`:
   ```python
   import os
   import sys
   from streamlit.web import cli as stcli

   if __name__ == "__main__":
       base_dir = os.path.dirname(os.path.abspath(__file__))
       app_path = os.path.join(base_dir, "app.py")
       sys.argv = [
           "streamlit",
           "run",
           app_path,
           "--global.developmentMode=false",
           "--server.address=localhost",
           "--browser.gatherUsageStats=false",
       ]
       sys.exit(stcli.main())
   ```
2. Build the one-folder executable distribution:
   ```powershell
   pyinstaller --noconfirm --onedir --windowed ^
       --name "HPLC_Viewer" ^
       --add-data "app.py;." ^
       --add-data "src;src" ^
       --add-data "config;config" ^
       --add-data ".streamlit;.streamlit" ^
       --copy-metadata streamlit ^
       run_desktop.py
   ```
3. Distribute the resulting `dist/HPLC_Viewer/` folder to target Windows workstations.
