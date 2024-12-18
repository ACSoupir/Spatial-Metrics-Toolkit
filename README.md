# Spatial Metrics Toolkit

## Table of Contents
1. [Overview](#overview)
2. [Input Data Requirements](#input-data-requirements)
3. [Preparing Your Data](#preparing-your-data)
4. [Usage Instructions](#usage-instructions)
5. [Output](#output)
6. [Installation](#installation)
7. [License](#license)

---

## Overview

Analyzing large spatial datasets in bioinformatics can quickly become computationally demanding, especially when traditional parallelization approaches like R's parallel `package` (the `mclapply` function) are used. These methods often duplicate the entire dataset across multiple cores, leading to excessive memory usage and potential RAM limitations on standard systems. The **Spatial Metrics Toolkit** is designed to overcome these challenges by adopting an efficient on-disk data processing strategy. Instead of passing large datasets between cores, the toolkit keeps data on disk and reads it in as needed, ensuring minimal memory overhead even for massive datasets.

This approach not only makes the toolkit more memory-efficient but also allows it to scale seamlessly across a range of hardware setups, from personal laptops to high-performance computing (HPC) clusters. When deployed on an HPC environment, the toolkit leverages the power of distributed computing by integrating with SLURM job scheduling. This enables researchers to distribute tasks across multiple nodes, maximizing computational resources and reducing analysis time significantly.

Whether you're working with single-node setups or large multi-node clusters, the **Spatial Metrics Toolkit** provides a flexible and scalable solution for extracting meaningful insights from spatial proteomics and transcriptomics data.

### Key Features
- **Automated Spatial Metric Calculation**: Effortlessly compute a wide range of spatial summary metrics for large and complex tissue microenvironment datasets.
- **Scalable Computing**: Supports SLURM for distributed computing on HPC clusters, enabling efficient scaling from single-node to multi-node environments.
- **Interactive Visualizations**: Generate detailed plots to visualize spatial distributions, clustering behavior, and more, enhancing data interpretation.
- **Flexible Data Input**: Works seamlessly with standard CSV data formats, making it easy to integrate with existing workflows.
- **Customizable Metrics**: Add new spatial summary metrics directly via the config.yml file and specify file paths, empowering users to tailor the toolkit to their unique analysis needs.

---

## Input Data Requirements

The input data must be a CSV file with the following required columns:

-   Sample Identifier (like `sample_id` or `id`)
-   **`x`**: The x-coordinate of the cell or point in the tissue.
-   **`y`**: The y-coordinate of the cell or point in the tissue.
-   **Marker Columns**: Binary columns (1/0) indicating the presence (1) or absence (0) of specific markers for each cell. Each marker should have its own column with a descriptive name (e.g., `CD3+`, `CD20+`).
-   **`compartment`** (optional): A column indicating the tissue compartment or region for each cell (e.g., `tumor`, `stroma`, `lymph`).

### Example Input Data

| x      | y      | CD3 | CD20 | CD68 | compartment |
|--------|--------|-----|------|------|-------------|
| 100.23 | 200.45 | 1   | 0    | 0    | tumor       |
| 150.67 | 250.89 | 0   | 1    | 0    | stroma      |
| 300.11 | 400.56 | 0   | 0    | 1    | lymph       |

### Notes:

-   The `x` and `y` columns must be in the same units and coordinate system for all points in the dataset.
-   Marker columns should use consistent naming conventions and binary values only.
-   The `compartment` column is optional but recommended for stratified analyses.

---

## Preparing Your Data

1.  **Check your data for missing values:** Ensure there are no missing values in the `x`, `y`, or marker columns. Missing values can cause errors during processing.
2.  **Convert marker positivity to binary values:** If your data contains continuous values for marker intensity, binarize them into 1 (positive) and 0 (negative) based on your chosen threshold.
3.  **Format tissue compartments (optional):** If applicable, include a column indicating tissue compartments. Use consistent and meaningful labels.
4.  **Save as CSV:** Ensure your file is saved in CSV format with a `.csv` extension.

**Tip**: Use tools like R, Python, or Excel for preprocessing. Scripts for common tasks are included in the `rpgm/` directory.
**Tip**: Use tools like R, Python, or Excel for preprocessing. Scripts for common tasks are included in the `rpgm/` directory.

### Setting Up the Configuration File

Before running the toolkit, ensure that the `config.yml` file is properly configured to match your dataset and analysis requirements. The configuration file is used to specify critical parameters for data processing and computation:

- **Column Names**: Set the correct column names for:
  - **`x_value`**: The column name for the x-coordinate of the cells (e.g., `x`).
  - **`y_value`**: The column name for the y-coordinate of the cells (e.g., `y`).
  - **`compartment`** (optional): The column name for tissue compartments, such as `tumor`, `stroma`, or `lymph`. If not applicable, this can be left undefined.

- **Metrics**: Specify the spatial metrics to calculate, such as Ripley’s K, G-function, or DBSCAN clustering, under the `metrics` parameter. For example:
  ```yaml
  metrics:
    - kest
    - gest
    - dbscan
  ```
- **Markers**:Define the marker columns used to identify cell types under the markers section. Example:
  ```yaml
  markers:
    - CD3
    - CD20
    - CD68

  ```
- **SLURM**: Enable or disable SLURM-based parallelization with the `slurm` parameter:
  - `TRUE` to use SLURM on an HPC.
  - `FALSE` to fall back to single-node parallel processing using the `parallel` package.

---

## Usage Instructions

### Running the Toolkit
1. Place your formatted CSV file in the `data/per-cell/` directory.
2. Configure the YAML file (`config.yml`) to specify analysis parameters.

#### Example Command
```bash
Rscript main.R --yaml config.yml --cores 4
```
- `--yaml`: Path to the YAML configuration file.
- `--cores`: Number of CPU cores to use for parallel processing.

---

## Output

The toolkit generates:
1. **Spatial Summary Metrics**: CSV files with computed metrics.
2. **Visualizations**: Plots showing spatial distribution and clustering.

### Example Folder Structure
```
output/
├── metrics/
|   ├──kest/
|   |   ├──TMA1_[3,B].tif.csv.gz
├── figures/
|   ├──barplot/
│   │   ├──marker_distribution.pdf
|   ├──metrics/
│   │   ├──kest/
│   │   │   ├──TMA1_[3,B].tif.pdf
|   ├──point/
│   │   ├──TMA1_[3,B].tif.pdf
```

---

## Installation

### Prerequisites
- **R (≥4.0)** with required packages installed (e.g., `parallel`, `yaml`). 
- SLURM installed for HPC usage (optional).

### Installation Steps
1. Clone the repository:
   ```bash
   git clone <repository_url>
   ```
2. Install required R packages:

    Packages should automatically be installed when first running the toolkit. It does help to manually install packages if there are errors. Have noted that system libraries are not always available on HPCs and requires communication with administrators for installing. If manually installing, can be done with:
   
   ```R
   install.packages(c("yaml", "parallel"))
   ```

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

