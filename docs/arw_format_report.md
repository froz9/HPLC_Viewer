# Technical Report: Waters Empower 3 ASCII Raw Data (`.arw`) Format Investigation

## 1. Executive Summary

An empirical investigation was conducted on all **64 Waters Empower 3 `.arw` files** located in the `raw_data/` directory (`HPLC_RAW_EXPORT1264.arw` through `HPLC_RAW_EXPORT1359.arw`), representing **16 chromatographic injections** across 4 exported channels per injection.

This report documents the exact byte-level encoding, line termination, delimiter structure, metadata headers, data payload layout, time and signal representations, differences between 2D single-channel exports (ELSD and extracted UV wavelengths) and 3D Photodiode Array (PDA) spectral matrix exports, and potential parsing complications.

---

## 2. Inventory of Inspected `.arw` Files

The 64 files correspond to 16 injections (4 files per injection: 1 ELSD 2D channel, 1 3D PDA spectral matrix, and 2 extracted 2D UV/PDA wavelength channels):

| File Range | Sample Name | Injection | Run Time (min) | Data Points | Exported Channels |
| :--- | :--- | :---: | :---: | :---: | :--- |
| `1264–1267` | `Diox.MeOH` | 1 | `0.0 – 10.0` | 6,001 | `ELSD Signal`, `2998` (3D), `2998 Ch1 254nm@1.2nm`, `2998 Ch2 365nm@1.2nm` |
| `1270–1273` | `Blank_Medium` | 1 | `0.0 – 32.0` | 19,201 | `ELSD Signal`, `2998` (3D), `2998 Ch1 254nm@1.2nm`, `2998 Ch2 365nm@1.2nm` |
| `1276–1345` | `01006-18-1` to `01006-18-12` (12 samples) | 1 | `0.0 – 32.0` | 19,201 | `ELSD Signal`, `2998` (3D), `2998 Ch1 254nm@1.2nm`, `2998 Ch2 365nm@1.2nm` |
| `1348–1351` | `Shutdown` | 1 | `0.0 – 20.1` | 12,061 | `ELSD Signal`, `2998` (3D, 170 WLs), `2998 Ch1 254nm@1.2nm`, `2998 Ch2 280nm@1.2nm` |
| `1356–1359` | `Blank_Medium` | 1 | `0.0 – 32.0` | 19,201 | `ELSD Signal`, `2998` (3D), `2998 Ch1 254nm@1.2nm`, `2998 Ch2 365nm@1.2nm` |

---

## 3. File Encoding and Line Termination

1. **Character Encoding**:
   - All 64 files contain strictly 7-bit **ASCII** bytes (`0x00–0x7F`; zero bytes `> 127` across all 64 files).
   - To accommodate sample names from other laboratories that may include accented characters or symbols (e.g., `µ`, `°`, Latin characters), the parser should attempt `utf-8` (with `utf-8-sig` BOM handling) and fall back cleanly to `cp1252` / `latin-1`.

2. **Line Termination (Critical Finding)**:
   - Every `.arw` file uses **Carriage Return only (`\r`, ASCII `0x0D`)** as its line terminator.
   - There are **zero Line Feed (`\n`, `0x0A`) bytes** in any of the 64 files.
   - Every file ends with a trailing `\r`.
   - **Implication for Python / Streamlit**: Standard binary `BytesIO.readline()` or `io.StringIO(text)` (which defaults to `newline='\n'`) will **fail to split lines** on `\r`-delimited files and read the entire file (up to 53.7 MB) as a single line. The parser must normalize `\r\n` and `\r` to `\n` (or split explicitly on `\r` / universal newlines) before parsing headers and numeric tables.

---

## 4. Delimiters and Quoting

- **Column Delimiter**: Tab character (`\t`, ASCII `0x09`).
- **String Quoting**: Metadata header names and string metadata values are enclosed in double quotes (`"`), while numeric metadata values (such as `Injection` = `1`) and all numeric data values are unquoted.
- **Leading Whitespace in 3D PDA Files**: In 3D PDA exports, the first column after `Wavelength\t` (and the first intensity column after `<time>\t` on each data row) is padded with two leading spaces (e.g., `Wavelength\t  190.1094\t191.2837...` and `0.00166665\t  -0.0350025\t...`). Whitespace must be stripped when parsing token values.

---

## 5. Metadata Structure

All inspected `.arw` files begin with a 2-line metadata table:

- **Row 0 (Metadata Column Headers)**:
  ```text
  "SampleName"	"Injection"	"Channel"	"Det. Units"
  ```
- **Row 1 (Metadata Values)**:
  - ELSD example (`HPLC_RAW_EXPORT1264.arw`):
    ```text
    "Diox.MeOH"	1	"ELSD Signal"	"LSU"
    ```
  - 2D Extracted UV/PDA Channel example (`HPLC_RAW_EXPORT1266.arw`):
    ```text
    "Diox.MeOH"	1	"2998 Ch1 254nm@1.2nm"	"AU"
    ```
  - 3D PDA Matrix example (`HPLC_RAW_EXPORT1265.arw`):
    ```text
    "Diox.MeOH"	1	"2998"	"AU"
    ```

> **Note on Variable Header Lengths**: While these 64 exports use a 2-row metadata block (or 4 rows including `Wavelength` and `Time` headers in 3D PDA files), Empower export methods can be configured with different sets of metadata columns or transposed key-value header rows. Therefore, the parser dynamically scans header rows until it locates the data section (either a `Wavelength\t...` + `Time` 3D header or consecutive numeric tab-separated rows).

---

## 6. Data Payload Structure

Two distinct `.arw` data structures exist in the dataset:

### 6.1 Type A: 2D Single-Channel Chromatogram (`.arw`, ~107–363 KB)

Immediately following the metadata lines (starting at Row 2), each row contains exactly **2 tab-separated numeric columns** with no column header line (`Time` / `Intensity` headers are omitted):

```text
0	-8.24e-005
0.001666667	-0.0001348
0.003333333	-0.0002279
0.005	-0.000324
```

- **Column 0**: Retention time in **minutes** (`float64`).
- **Column 1**: Detector response intensity (`float64`, formatted in standard decimal or 3-digit exponent scientific notation such as `-8.24e-005`).

### 6.2 Type B: 3D PDA Multi-Wavelength Matrix (`.arw`, ~16.0–53.8 MB)

When the entireWaters 2998 PDA 3D channel (`Channel = "2998"`) is exported:

- **Row 2**: Tab-separated wavelength header row starting with literal `Wavelength`, followed by $N_{\lambda}$ floating-point wavelengths in **nm**:
  - Standard methods (`1265–1343`, `1357`): **261 wavelengths** from `190.1094 nm` to `500.2398 nm` ($\Delta\lambda \approx 1.188\text{ nm}$).
  - `Shutdown` method (`1349`): **170 wavelengths** from `199.5079 nm` to `399.495 nm`.
- **Row 3**: Single token `Time`.
- **Row 4 to End**: Tab-separated numeric rows of length $1 + N_{\lambda}$, where Column 0 is retention time (minutes) and Columns $1 \dots N_{\lambda}$ are absorbance values (`AU`) at each corresponding diode wavelength:
  ```text
  Wavelength	  190.1094	191.2837	...	500.2398
  Time
  0	  0	0	...	0
  0.00166665	  -0.0350025	0.0064186	...	-4.85e-005
  ```

---

## 7. Time Units, Sampling Rate, and Precision Artifacts

1. **Time Units**: Retention time is expressed in **minutes** (`0.0` to `10.0`, `20.1`, or `32.0` min).
2. **Sampling Rate**: All 64 files are sampled at **10 Hz** ($\Delta t = \frac{1}{600}\text{ min} \approx 0.001666667\text{ min} = 0.1\text{ s}$).
3. **Significant-Digit Rounding Artifact**:
   - Empower formats exported floating-point numbers to **7 significant digits**.
   - Consequently, for retention times $t \ge 10.0\text{ min}$, timestamps are rounded to 5 decimal places (e.g., `31.995`, `31.99667`, `31.99833`, `32`).
   - First differences $\Delta t_i = t_{i+1} - t_i$ alternate between `0.00166` and `0.00167` ($\pm 7 \times 10^{-6}\text{ min}$, or $\sim 0.42\%$ relative rounding jitter).
   - Furthermore, in 3D PDA files (`1271`, `1277`, etc.), Empower steps by `0.00166691` min, reaching `32.00469` min at sample 19,201 instead of `32.00000` min.
   - **Validation Rule**: Sampling interval uniformity validation must account for 7-significant-figure ASCII rounding jitter ($\le 1\%$ relative deviation from median $\Delta t$) so that normal Empower exports are not falsely flagged as having corrupted/unequal sampling intervals, while still reporting true step changes or missing time points.

---

## 8. Differences Between UV/PDA and ELSD Files

| Property | UV / PDA Extracted Channels | 3D PDA Matrix | ELSD Channels |
| :--- | :--- | :--- | :--- |
| **Channel String** | `"2998 Ch1 254nm@1.2nm"`, `"2998 Ch2 365nm@1.2nm"`, `"2998 Ch2 280nm@1.2nm"` | `"2998"` | `"ELSD Signal"` |
| **Detector Units** | `"AU"` (Absorbance Units) | `"AU"` (Absorbance Units) | `"LSU"` (Light Scattering Units) |
| **Wavelength Info** | Embedded in `Channel` string (e.g., `254 nm`, bandwidth `1.2 nm`) | Explicit `Wavelength` row (`190.1–500.2 nm`) | None (evaporative light scattering is wavelength-independent) |
| **Signal Range in Dataset** | `-0.018` to `+3.90 AU` | `-0.4` to `+3.9 AU` | `-8.84` to `+230.22 LSU` |
| **Baseline / Negative Values** | Small negative baseline drift (`-10^{-5}` to `-10^{-2} AU`) is normal | Normal baseline noise | Negative baseline offset (`-0.1` to `-8.8 LSU`) is common and physically normal after autozero |

---

## 9. Potential Complications and Parser Design Requirements

1. **CR-Only (`\r`) Line Endings**: Must be handled at the byte/stream level before line iteration or CSV parsing.
2. **Duplicate Sample Names Across Injections**: `Blank_Medium` (Injection 1) appears twice in the sequence (`1270–1273` and `1356–1359`). The application must preserve the exact `SampleName` metadata while assigning each imported file/channel a unique internal identifier and disambiguated label when needed.
3. **Multiple Channels per Sample Injection**: Because a single sample (`01006-18-1`) has separate files for `ELSD`, `254 nm`, `365 nm`, and `3D PDA`, each file/wavelength represents an independent chromatographic trace.
4. **Large 3D PDA Files (~53 MB each)**:
   - Parsing must use vectorized C-engine reading (`pandas.read_csv` after byte-level `\r` $\to$ `\n` normalization, which parses 53 MB in `< 0.3 s`).
   - Rather than flooding the sample list with 261 separate traces per 3D file, importing a 3D PDA `.arw` file should extract user-requested target wavelength(s) (defaulting to `254 nm` and `365 nm`, or user-specified wavelengths in the UI) while keeping the available wavelength list accessible so users can extract any additional wavelength channel on demand.
5. **Negative Intensities**: Every single file in `raw_data/` contains negative baseline values (ranging from 46 to 16,033 points per chromatogram). As required by the specification, negative intensities must **never** be deleted or clipped; their count and minimum value are reported transparently in the validation summary.
