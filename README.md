
# ERG Signal Analysis and Visualization

This repository contains MATLAB code for analyzing and visualizing Electroretinography (ERG) signals. The code applies advanced signal processing techniques to extract meaningful insights and generates publication-quality plots. Designed for research and education, it is adaptable and user-friendly.

---

## Features
- **Signal Preprocessing**: Includes detrending and low-pass filtering.
- **Fourier Transform (FT)**: Power Spectrum computation and visualization.
- **Short-Time Fourier Transform (STFT)**: Generates spectrograms for time-frequency analysis.
- **Continuous Wavelet Transform (CWT)**: Produces scalograms for detailed time-frequency insights.
- **Discrete Wavelet Transform (DWT)**: Displays energy heatmaps across decomposition levels.
- **Customizable Outputs**: Flexible layout and parameter adjustment for tailored results.

---

## Getting Started

### Prerequisites

- **MATLAB**: Ensure MATLAB is installed along with the following toolboxes:
  - Signal Processing Toolbox
  - Wavelet Toolbox

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/yshwetar/ERG-Signal-Analysis-For-Users.git
   ```
2. Navigate to the project directory in MATLAB.

---

## File Structure

### Input Data

The input data must be in `.mat` or `.csv` format:
- **Time Column**: Represents the time vector (not used directly in analysis).
- **Signal Column(s)**: Contains the ERG signals for processing.
- **Example**:
  ```csv
  time,signal
  0.000,0.123
  0.001,0.234
  ...
  ```

### Test File

The file `0001.mat` is included as a test file. It is part of the **PERG-IOBA Dataset** from the original publication available at:
- PhysioNet: [https://physionet.org/content/perg-ioba-dataset/1.0.0/](https://physionet.org/content/perg-ioba-dataset/1.0.0/)
- Nature Scientific Data: [https://www.nature.com/articles/s41597-024-03857-1](https://www.nature.com/articles/s41597-024-03857-1)

For more information, please refer to the publication.

---

## Scripts and Associated Figures

The provided MATLAB scripts correspond to the figures in the publication:
1. **`PS_STFT.m`**: Power Spectrum and Short-Time Fourier Transform.
   - *Figure*: Manuscript Figure 1.
2. **`CWT.m`**: Continuous Wavelet Transform.
   - *Figure*: Manuscript Figure 2.
3. **`DWT.m`**: Discrete Wavelet Transform.
   - *Figure*: Manuscript Figure 3.
4. **`FT.m`**: Fourier Transform analysis.
   - *Figure*: Supplementary Figure 2.
5. **`Windowing.m`**: Exploring windowing techniques for STFT.
   - *Figure*: Supplementary Figure 3.
6. **`Technical_Example.m`**: Demonstrates all methods in a single analysis.
   - *Figure*: Supplementary Figure 1.
7. **`Publication_Code.m`**: Integrates all methods for reproducibility.
8. **`User_Code.m`**: A template for users to analyze custom ERG data.

---

## Usage

1. **Prepare the Input File**:
   - Update `base_path` in the script with the path to your `.mat` or `.csv` file.

2. **Customize Parameters**:
   - Modify parameters such as `sampling_rate`, `filter_order`, and `cutoff_frequency`.

3. **Run the Script**:
   - Execute the desired script in MATLAB:
     ```matlab
     load('0001.mat'); % Load sample ERG data
     run('CWT.m');     % Run Continuous Wavelet Transform
     ```

4. **View Outputs**:
   - The scripts generate plots such as power spectra, spectrograms, scalograms, and decomposition heatmaps.

---

## Outputs

The scripts generate visualizations including:
1. **Filtered ERG Signal**: Time-domain plot of preprocessed signals.
2. **Power Spectrum**: Frequency-domain visualization.
3. **STFT Spectrogram**: Temporal frequency analysis.
4. **CWT Scalogram**: Rich time-frequency representations.
5. **DWT Heatmap**: Hierarchical decomposition energy visualization.

---

## Citation

If you use this repository, please cite:
- **Paper**: *Signal Processing in Electroretinography* (Add full citation once published).
- **Authors**: Yousif Shwetar, David Lalush, Brett Jeffrey, Melissa Haendel.
- **DOI**: (Add once available).

For the test file `0001.mat` and '0019.mat', please also see and cite the PERG-IOBA Dataset if necessary:
- **Dataset**: [https://physionet.org/content/perg-ioba-dataset/1.0.0/](https://physionet.org/content/perg-ioba-dataset/1.0.0/)
- **Original Paper**: [https://www.nature.com/articles/s41597-024-03857-1](https://www.nature.com/articles/s41597-024-03857-1)

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## Contact

For questions or suggestions, contact **Yousif Shwetar** at [yousif_shwetar@med.unc.edu].
