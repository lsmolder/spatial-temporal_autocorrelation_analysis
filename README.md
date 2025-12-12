# Autocorrelation Analysis Script

This script calculates Temporal Autocorrelation (TA) and Spatial Autocorrelation (SA) for brain imaging data, following the methods described in Shinn et al. (2023).

## Requirements

- Python 3
- See `requirements.txt` for dependencies.

## Usage

```bash
python calculate_autocorrelation.py --fc_matrix <path_to_fc_matrix.csv> --timeseries <path_to_timeseries.csv> --atlas <path_to_atlas.nii.gz> --output <output_directory>
```

### Arguments

- `--fc_matrix`: Path to the Functional Connectivity matrix CSV.
    - Format: First row and first column should be region labels.
- `--timeseries`: Path to the fMRI time series CSV.
    - Format: First row should be region labels. Data starts from the second row (no time column).
- `--atlas`: Path to the NIfTI atlas file (`.nii` or `.nii.gz`) defining the brain regions.
- `--output`: Directory where results will be saved.

## Outputs

The script generates the following files in the output directory:

1.  `results_summary.csv`: Contains the Global TA, SA-lambda (length scale), and SA-infinity (asymptote).
2.  `results_by_region.csv`: Contains the TA value for each region, along with the extracted (x, y, z) coordinates.
3.  `sa_curve_fit_data.csv`: Contains the data used for the SA exponential decay fit (distance bins, mean correlation, and fitted model), useful for plotting.

## Methodology

### Temporal Autocorrelation (TA)
Calculated as the lag-1 autocorrelation of the time series for each region. Global TA is the mean of all regional TAs.

### Spatial Autocorrelation (SA)
Calculated by fitting an exponential decay model to the relationship between Functional Connectivity (correlation) and Euclidean distance between regions.

Model: $y = SA_{inf} + (1 - SA_{inf}) \times e^{-x / SA_{\lambda}}$

- $SA_{\lambda}$: Length scale (rate of decay).
- $SA_{inf}$: Asymptote (baseline correlation at long distances).
