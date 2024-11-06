
# Photobleach Analysis Tool

This repository contains a MATLAB-based tool for analyzing photobleaching events in fluorescence microscopy data. The tool provides functionalities to load, analyze, and visualize photobleaching dynamics, with various customizable parameters and tools for identifying tethered molecules and integrating intensity of individual molecules in Cy3 and Cy5 channels.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [1. Generate Reference Set](#1-generate-reference-set)
  - [2. Generate Input Parameter File](#2-generate-input-parameter-file)
  - [3. Identify Tethered Molecules and Drift Correction](#3-identify-tethered-molecules-and-drift-correction)
  - [4. Identify Photobleaching Events](#4-identify-photobleaching-events)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites

1. **MATLAB Toolboxes**:
   - Image Processing Toolbox
   - Statistics and Machine Learning Toolbox

2. **External Packages**:
   - [Fast PSF Fitting](https://github.com/scstein/FastPsfFitting)
   - [Bio-formats plugin (bfmatlab)](https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/index.html)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/morenoat/Photobleach_analysis.git
   ```
2. Add the `photobleach_scripts` folder and experimental data folder to the MATLAB path:
   ```matlab
   addpath(genpath(pwd));
   ```

## Usage

###   Open MATLAB live script:   Analyze\_photobleaching 

-This script contains all the information and excutable sections to analyze datasets.

### 1. Generate Reference Set

This script uses reference images to identify pairs of molecules in Cy3 and Cy5 channels for calibration between the two. The output is a MATLAB structure containing reference points.

```matlab
Generate_reference_set(parms_output_folder_name, ID_molecules, plot_each_image, last_round_plot, path_to_data);
```

- **Parameters**:
  - `path_to_data`: Path to folder with images for calibration.
  - `parms_output_folder_name`: Output folder name for saving the reference data.
  - `ID_molecules`: Set `true` to identify molecules in images.
  - `plot_each_image`: Set `true` to visualize each image.

### 2. Generate Input Parameter File

This script configures parameters related to the camera detector, excitation frames, and image acquisition timing. It generates a `parms_in.txt` file for use in subsequent steps.

```matlab
generate_parms(number_cy3_excitation_frames, number_cy5_excitation_frames, exposuret, delay, starting_frame, parms_output_folder_name);
```

- **Parameters**:
  - `number_cy3_excitation_frames` and `number_cy5_excitation_frames`: Frames per cycle for Cy3 and Cy5.
  - `exposuret`: Laser exposure time in seconds.
  - `delay`: Delay time between exposures.
  - `starting_frame`: First frame for analysis.
  - `parms_output_folder_name`: Output folder for parameter file.

### 3. Identify Tethered Molecules and Drift Correction

This step uses the generated reference set to integrate molecule intensity over time and identify regions of interest (ROI) for Cy3 and Cy5 channels.

```matlab
define_movie_segments(parms_output_folder_name, output_name_data_integration, input_name_data_integration, name_reference_file);
```

- **Parameters**:
  - `input_name_data_integration`: Name of the movie file.
  - `output_name_data_integration`: Output folder for integrated data.
  - `name_reference_file`: Reference file created in Step 1.
  - `parms_output_folder_name`: Folder containing parameter file.

### 4. Identify Photobleaching Events

To analyze photobleaching events, launch the GUI using the following command:

```matlab
app = ScanForBleach();
```

This GUI provides an interactive interface for analyzing and identifying photobleaching events based on user-defined parameters.

## Contributing

Contributions are welcome! Please fork this repository and submit a pull request with your improvements.

## License

This project is licensed under the MIT License.
