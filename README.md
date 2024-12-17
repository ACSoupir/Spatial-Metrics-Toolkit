# Spatial Metrics Toolkit

## Overview

The Spatial Metrics Toolkit is designed to automate the computation of spatial summary metrics for tissue microenvironment data. This tool supports scalable and flexible analysis, enabling researchers to extract meaningful insights from spatial proteomics and transcriptomics data.

Within the config file, the specification of `slurm` as true or false will tell the Toolkit what kind of backend to use. If set to `TRUE`, the Toolkit will attempt to find the location of SLURM on the system and if present, will spawn jobs using it. If not able to locate it, then the Toolkit will default back to using the `parallel` package. The implementation of SLURM will allow for jobs to be distributed across nodes on a multi-node HPC while `parallel` has the limitation of working on a single node.

## Input Data Format

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

## Preparing Your Data

1.  **Check your data for missing values:** Ensure there are no missing values in the `x`, `y`, or marker columns. Missing values can cause errors during processing.
2.  **Convert marker positivity to binary values:** If your data contains continuous values for marker intensity, binarize them into 1 (positive) and 0 (negative) based on your chosen threshold.
3.  **Format tissue compartments (optional):** If applicable, include a column indicating tissue compartments. Use consistent and meaningful labels.
4.  **Save as CSV:** Ensure your file is saved in CSV format with a `.csv` extension.

## Using the Toolkit

Once your data is formatted correctly: 1. Place your CSV file in the `data/per-cell/` directory specified by the toolkit. 2. Run the tool using the appropriate command or script (refer to the [Usage Instructions](#usage-instructions) below). 3. Outputs will be saved in a structured folder with summary metrics and visualizations.

## Usage Instructions

1.  **Other Parameters in YAML**

    -   `metrics`: Choose spatial summary metrics [Ripley's K, G-function, DBSCAN] for analyzing point pattern distributions and clustering behavior. 🔍
        -   `- dbscan`: DBSCAN: Density-based spatial clustering algorithm that groups points with many neighbors within a radius ε while detecting noise, enabling discovery of arbitrary-shaped clusters in spatial data. 📊
        -   `- kest`: Estimates Ripley's K-function to analyze spatial point patterns by measuring inter-point distances and clustering intensity at multiple scales within any window shape. 📊
        -   `- gest`: Estimates nearest neighbor distance distribution G(r) to analyze spatial point patterns by measuring point-to-point proximity relationships in any window shape. 📊
    -   `variables`: Use this to specify variables for analyses.
        -   `markers`: columns that are used to identify the cell types
        -   `x_value`/`y_value`: column specifying the location of centroid of cells
    -   `slurm`: `TRUE`/`FALSE` for whether or not to use SLURM
        -   If set to `TRUE`, will look to ensure SLURM is available on the system. If not available, will fall back to single-node parallel processing

2.  **Command-line Execution**\
    Use the following command to run the toolkit:

    ``` bash
    Rscript main.R --yaml config.yml --cores 4
    ```

    `--yaml`: path to the yaml configuration file. 🛠️

    `--cores`: number of CPU cores to use for parallel processing of data.

## Output

The toolkit generates:

-   **Spatial Summary Metrics:** CSV files containing computed spatial metrics.

-   **Visualizations:** Plots showing spatial distribution and summary analyses.
