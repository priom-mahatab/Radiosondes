# Radiosonde Data Processing and Skew-T Diagram Generation

This project processes radiosonde data files, extracts relevant atmospheric data, and generates Skew-T diagrams for meteorological analysis. The script selects specific files based on their timestamps, removes duplicate pressure values, and computes thermodynamic parameters using the MetPy library.

There are three functions, one for files at 11 UTC, another for 17 UTC, and the last one for overlapping the two SkewTs from the same date.

## Functions

- `skewT_generator_11(input_folder, output_folder)`: Generates SkewT, Hodograph, and a Stats section for files with sonde data from 11 UTC.
- `skewT_generator_17(input_folder, output_folder)`: Generates SkewT, Hodograph, and a Stats section for files with sonde data from 17 UTC.
- `overlapping_skewTs(input_folder, output_folder)`: Generates SkewT by overlapping the data from the same days.

## Features

- **Automated File Selection**: Processes files matching specific timestamps (11XX or 17XX hours).
- **Data Cleaning**: Removes duplicate pressure values to ensure accurate calculations.
- **Skew-T Diagram Generation**: Plots temperature, dew point, and other fiducial lines.
- **CAPE/CIN Calculations**: Computes Convective Available Potential Energy (CAPE) and Convective Inhibition (CIN) using MetPy.
- **Legend Positioning**: Ensures the legend is always positioned within the top-right corner of the plot.

## Prerequisites

- Python 3.x
- Required Libraries:
  - numpy
  - pandas
  - metpy
  - matplotlib
  - jupyter notebook

## Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/priom-mahatab/Physics-Research.git
   ```
2. Install dependencies:
   ```sh
   pip install numpy pandas metpy matplotlib jupyter
   ```

## Usage

Run the script to process radiosonde data and generate Skew-T diagrams:

```sh
jupyter notebook "TRACER M1.ipynb"
```

### Handling Errors

- **Duplicate Pressure Warnings**: The script automatically removes duplicate pressure values.
- **MetPy Calculation Errors**: If a calculation fails, the script gracefully skips the affected computation.

## Output

- Processed Skew-T diagrams are saved in the `output/` directory with unique filenames based on the input data.

## Notes

- Ensure the input CSV files follow the naming convention: `housondewnpnM1.b1.YYYYMMDD.HHMMSS.custom.csv`.
- Adjust the legend positioning in the plot by modifying:
  ```python
  skew.ax.legend(loc='upper right', bbox_to_anchor=(0.85, 1), fontsize=10, frameon=True)
  ```

## License

This project is licensed under the MIT License.

