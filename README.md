# ERG Wavelet Analysis Toolkit

This repository contains MATLAB code for analyzing and visualizing Electroretinography (ERG) signals using advanced wavelet transforms. The code implements the methods described in "A Practical Introduction to Wavelet Analysis in Electroretinography" by Shwetar et al., providing researchers and clinicians with tools to extract time-frequency insights from ERG recordings.

---

## Features
- **File Support**: Handles multiple file formats (.xlsx, .mat, .csv, .asc)
- **Signal Preprocessing**: Automatic detrending to remove DC offset
- **Amplitude Spectrum**: Reveals frequency components in ERG signals
- **Continuous Wavelet Transform (CWT)**: Produces detailed scalograms visualizing energy across time and frequency
- **Discrete Wavelet Transform (DWT)**: Displays decomposition energy heatmaps with discrete frequency bands
- **Comparison Mode**: Analyze two signals simultaneously (e.g., normal vs. pathological ERGs)
- **Comprehensive Visualization**: Publication-quality 4-panel plots for each signal

---

## Getting Started

### Prerequisites

- **MATLAB**: Ensure MATLAB is installed along with the following toolboxes:
  - Signal Processing Toolbox
  - Wavelet Toolbox

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yshwetar/ERG-Wavelet-Analysis.git
   ```
2. Navigate to the project directory in MATLAB.

---

## Usage

1. **Run the Main Script**:
   ```matlab
   run('ERG_WT_Analysis.m')
   ```

2. **Input Selection**:
   - The script will prompt you to select an ERG data file
   - You'll have the option to load a second file for comparison
   - If no file is selected, the script will run with demo data

3. **Data Formats**:
   - **Time Series**: The script automatically detects time columns when available
   - **Sampling Rate**: Either detected from time column or user-prompted
   - **Support**: Works with .xlsx, .mat, .csv, and .asc formats

4. **Handling Different Signal Lengths**:
   - When comparing two signals of different lengths, you can choose to:
     - Truncate both signals to the length of the shorter one
     - Pad the shorter signal with zeros to match the longer one

---

## Analysis Methods

The toolkit implements the wavelet analysis methods described in the manuscript:

### 1. Amplitude Spectrum
### 2. Continuous Wavelet Transform (CWT)
### 3. Discrete Wavelet Transform (DWT)

---

## Output Visualization

The script generates a comprehensive figure with 4 panels for each signal:

1. **Time-Domain Plot**: Shows the original ERG waveform
2. **Amplitude Spectrum**: Displays frequency components
3. **CWT Scalogram**: Shows energy across continuous time-frequency space
4. **DWT Heatmap**: Shows energy across discrete frequency bands

Additionally, the console output provides:
- Signal length and duration
- Amplitude ranges
- Dominant frequency
- Sampling rate

---

## Case Study Application

As demonstrated in the manuscript, this toolkit can be used to analyze:

- **Standard ISCEV ERG Recordings**:
  - LA Flicker response
  - LA 3.0 cd·s/m²
  - DA 3.0 cd·s/m²
  - DA 0.01 cd·s/m²

- **Normal vs. Pathological Comparison**: 
  - The toolkit is particularly useful for comparing normal ERGs with those from patients with conditions like congenital stationary night blindness (CSNB)
  - Wavelet methods can reveal subtle differences in time-frequency characteristics that may not be apparent in standard time-domain analysis

---

## Citation

If you use this toolkit in your research, please cite:
- **Paper**: Shwetar Y, Lalush D, McAnany J, Jeffrey B, Haendel M. "A Practical Introduction to Wavelet Analysis in Electroretinography" (Add journal and DOI once published)
- **Code**: GitHub repository: https://github.com/yshwetar/ERG-Wavelet-Analysis

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## Contact

For questions or suggestions, contact **Yousif Shwetar** at [yousif_shwetar@med.unc.edu].
