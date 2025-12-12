import argparse
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.optimize import curve_fit
from scipy.spatial.distance import pdist, squareform
import os
import sys

def load_and_validate_data(fc_path, ts_path, atlas_path):
    """
    Loads the FC Matrix, Timeseries, and Atlas.
    Validates that they have matching dimensions and region labels.
    """
    print(f"Loading FC Matrix from: {fc_path}")
    # Load FC Matrix (assumes labels in first row and first column)
    fc_df = pd.read_csv(fc_path, header=0, index_col=0)
    
    # Handle NaNs in FC matrix: replace with 0
    # Common issue: diagonals or missing connections might be NaN
    if fc_df.isnull().values.any():
        print(f"  Warning: FC Matrix contains {fc_df.isnull().sum().sum()} NaNs. Replacing with 0.")
        fc_df = fc_df.fillna(0)
        
    fc_matrix = fc_df.values
    fc_labels = fc_df.index.astype(float).astype(int).values
    
    print(f"  Shape: {fc_matrix.shape}")
    print(f"  Found {len(fc_labels)} regions in FC Matrix.")

    print(f"Loading Timeseries from: {ts_path}")
    # Load Timeseries (assumes labels in first row, data starts from row 1)
    # Note: User said "first row is labels", "first column is actual datapoints"
    # This implies header=None, and we manually split
    ts_raw = pd.read_csv(ts_path, header=None)
    
    # Extract labels from the first row
    ts_labels_raw = ts_raw.iloc[0, :].values
    ts_labels = ts_labels_raw.astype(float).astype(int)
    
    # Extract data (all rows after the first)
    ts_data = ts_raw.iloc[1:, :].values
    
    print(f"  Shape: {ts_data.shape} (Timepoints x Regions)")
    print(f"  Found {len(ts_labels)} regions in Timeseries.")

    # Validate Consistency
    if len(fc_labels) != len(ts_labels):
        print("\nERROR: Mismatch in number of regions!")
        print(f"  FC Matrix: {len(fc_labels)}")
        print(f"  Timeseries: {len(ts_labels)}")
        sys.exit(1)
        
    # Check if labels match (using set comparison for robustness against ordering)
    if set(fc_labels) != set(ts_labels):
        print("\nERROR: Mismatch in region labels!")
        print(f"  FC Labels (first 5): {fc_labels[:5]}")
        print(f"  TS Labels (first 5): {ts_labels[:5]}")
        sys.exit(1)
        
    print("Loading Atlas...")
    atlas_img = nib.load(atlas_path)
    atlas_data = atlas_img.get_fdata()
    atlas_affine = atlas_img.affine
    
    unique_atlas_labels = np.unique(atlas_data)
    unique_atlas_labels = unique_atlas_labels[unique_atlas_labels > 0].astype(int) # Exclude 0 (background)
    
    print(f"  Found {len(unique_atlas_labels)} unique regions in Atlas.")
    
    # Validate Atlas against CSVs
    # We allow the atlas to have MORE regions (e.g., whole brain atlas vs masked data)
    # But we must ensure all regions in our CSVs exist in the Atlas
    missing_in_atlas = set(fc_labels) - set(unique_atlas_labels)
    if missing_in_atlas:
        print("\nERROR: The following regions found in CSVs are missing from the Atlas:")
        print(list(missing_in_atlas)[:10]) # Print first 10
        sys.exit(1)
        
    print("\n✅ Data validation passed! All inputs are consistent.")
    
    return fc_matrix, ts_data, fc_labels, atlas_data, atlas_affine

def extract_centroids(region_labels, atlas_data, atlas_affine):
    """
    Extracts the (x, y, z) Center of Mass (CoM) for each region label.
    Converts from Voxel coordinates to World coordinates (mm) using the affine.
    """
    print("\nExtracting region centroids...")
    centroids = []
    
    from scipy.ndimage import center_of_mass
    
    for label in region_labels:
        # Create a boolean mask for the current region
        mask = (atlas_data == label)
        
        if not np.any(mask):
            print(f"  Warning: Region {label} not found in atlas volume (should have been caught by validation).")
            centroids.append([np.nan, np.nan, np.nan])
            continue
            
        # Calculate center of mass in voxel coordinates
        # center_of_mass returns (z, y, x) or (x, y, z) depending on data layout? 
        # It returns indices in the order of dimensions.
        com_voxel = center_of_mass(mask)
        
        # Convert to world coordinates using the affine matrix
        # nibabel's apply_affine works with (x, y, z)
        com_world = nib.affines.apply_affine(atlas_affine, com_voxel)
        
        centroids.append(com_world)
        
    return np.array(centroids)

def calculate_ta(timeseries_data, region_labels):
    """
    Calculates Temporal Autocorrelation (TA) for each region.
    
    METHODOLOGY:
    -------------
    TA is defined as the lag-1 autocorrelation of the fMRI time series.
    
    Formula:
      TA = Correlation( X(t), X(t+1) )
      
    Where X(t) is the signal at time t, and X(t+1) is the signal shifted by 1 TR.
    
    Steps:
    1. For each region (column in timeseries_data):
    2. Create two vectors: 
       - v1: data from index 0 to N-2
       - v2: data from index 1 to N-1
    3. Calculate the Pearson correlation coefficient between v1 and v2.
    4. Global TA is the mean of all regional TA values.
    
    Ref: Shinn et al. 2023
    """
    print("\nCalculating Temporal Autocorrelation (TA)...")
    
    ta_values = []
    
    # Iterate over each region (column)
    for i in range(timeseries_data.shape[1]):
        signal = timeseries_data[:, i]
        
        # Calculate lag-1 autocorrelation
        # np.corrcoef returns a 2x2 matrix, we want the off-diagonal element [0, 1]
        ta = np.corrcoef(signal[:-1], signal[1:])[0, 1]
        ta_values.append(ta)
        
    ta_values = np.array(ta_values)
    global_ta = np.nanmean(ta_values)
    
    print(f"  Global TA (Mean Lag-1 Autocorr): {global_ta:.4f}")
    
    # Create a DataFrame for per-region results
    ta_df = pd.DataFrame({
        'Region_Label': region_labels,
        'TA_Lag1': ta_values
    })
    
    return global_ta, ta_df

def calculate_sa(fc_matrix, centroids):
    """
    Calculates Spatial Autocorrelation (SA) parameters.
    
    METHODOLOGY:
    -------------
    SA describes how functional connectivity (FC) decays as the physical distance 
    between brain regions increases.
    
    Model: Exponential Decay
      y = SA_inf + (1 - SA_inf) * exp( -x / SA_lambda )
      
      Where:
      - y:  Functional Connectivity (Correlation)
      - x:  Euclidean Distance (mm)
      - SA_lambda: Length scale (rate of decay)
      - SA_inf:    Asymptote (baseline correlation at long distances)
    
    Steps:
    1. Calculate pairwise Euclidean distances between all region centroids.
    2. Extract the upper triangle of the FC matrix (unique pairs).
    3. Extract the corresponding upper triangle of the Distance matrix.
    4. Bin the data by distance (1mm bins) to reduce noise and weight equally.
       - Calculate the mean FC for all pairs falling within each 1mm bin.
    5. Fit the exponential decay curve to these binned values using `scipy.optimize.curve_fit`.
    
    Ref: Shinn et al. 2023
    """
    print("\nCalculating Spatial Autocorrelation (SA)...")
    
    # 1. Calculate Euclidean Distance Matrix
    print("  Computing pairwise distances...")
    dist_matrix = squareform(pdist(centroids))
    
    # 2. Extract Upper Triangles (excluding diagonal)
    n_regions = fc_matrix.shape[0]
    triu_indices = np.triu_indices(n_regions, k=1)
    
    flat_fc = fc_matrix[triu_indices]
    flat_dist = dist_matrix[triu_indices]
    
    # 3. Binning (1mm bins)
    print("  Binning data by distance (1mm bins)...")
    max_dist = int(np.ceil(flat_dist.max()))
    bins = np.arange(0, max_dist + 1, 1)
    
    bin_centers = []
    bin_means = []
    
    for i in range(len(bins) - 1):
        # Find pairs with distance in current bin range
        mask = (flat_dist >= bins[i]) & (flat_dist < bins[i+1])
        
        if np.any(mask):
            # Calculate mean correlation, filtering out NaNs in FC data if any remain
            valid_corrs = flat_fc[mask]
            valid_corrs = valid_corrs[~np.isnan(valid_corrs)]
            
            if len(valid_corrs) > 0:
                mean_corr = np.mean(valid_corrs)
                bin_centers.append((bins[i] + bins[i+1]) / 2)
                bin_means.append(mean_corr)
            
    bin_centers = np.array(bin_centers)
    bin_means = np.array(bin_means)
    
    # Filter out any bins that might still be NaN (just in case)
    valid_bins = ~np.isnan(bin_means)
    bin_centers = bin_centers[valid_bins]
    bin_means = bin_means[valid_bins]
    
    print(f"  Data points for fitting: {len(bin_centers)} bins")
    
    # Check if we have enough points to fit
    if len(bin_centers) < 3:
        print("  ERROR: Not enough data points (bins) to fit curve. Need at least 3.")
        return np.nan, np.nan, bin_centers, bin_means, None
    
    # 4. Define Exponential Decay Function
    def exponential_decay(x, sa_lambda, sa_inf):
        # Formula: A + (1-A) * exp(-x/L)
        # where A = sa_inf, L = sa_lambda
        return sa_inf + (1 - sa_inf) * np.exp(-x / sa_lambda)
    
    # 5. Fit Curve
    print("  Fitting exponential decay model...")
    try:
        # Initial guess (p0): lambda=10mm, inf=0.1
        # Bounds: lambda (0, 200), inf (-1, 1)
        popt, pcov = curve_fit(
            exponential_decay, 
            bin_centers, 
            bin_means,
            p0=[10.0, 0.1],
            bounds=([0.001, -1.0], [200.0, 1.0]),
            maxfev=5000
        )
        
        sa_lambda, sa_inf = popt
        print(f"  SA-lambda (Length Scale): {sa_lambda:.4f} mm")
        print(f"  SA-inf (Asymptote): {sa_inf:.4f}")
        
        return sa_lambda, sa_inf, bin_centers, bin_means, exponential_decay(bin_centers, *popt)
        
    except Exception as e:
        print(f"  ERROR: Curve fitting failed: {e}")
        return np.nan, np.nan, bin_centers, bin_means, None

def main():
    parser = argparse.ArgumentParser(description="Calculate Temporal and Spatial Autocorrelation (Shinn et al. 2023)")
    parser.add_argument("--fc_matrix", required=True, help="Path to FC Matrix CSV")
    parser.add_argument("--timeseries", required=True, help="Path to Timeseries CSV")
    parser.add_argument("--atlas", required=True, help="Path to Atlas NIfTI file")
    parser.add_argument("--output", required=True, help="Output directory")
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        
    # 1. Load and Validate
    fc_matrix, ts_data, region_labels, atlas_data, atlas_affine = load_and_validate_data(
        args.fc_matrix, args.timeseries, args.atlas
    )
    
    # 2. Extract Centroids
    centroids = extract_centroids(region_labels, atlas_data, atlas_affine)
    
    # 3. Calculate TA
    global_ta, ta_df = calculate_ta(ts_data, region_labels)
    
    # Add coordinates to the region dataframe
    ta_df['x'] = centroids[:, 0]
    ta_df['y'] = centroids[:, 1]
    ta_df['z'] = centroids[:, 2]
    
    # 4. Calculate SA
    sa_lambda, sa_inf, bin_x, bin_y, fit_y = calculate_sa(fc_matrix, centroids)
    
    # 5. Save Results
    print("\nSaving results...")
    
    # Summary CSV
    summary_df = pd.DataFrame({
        'Metric': ['Global_TA', 'SA_Lambda', 'SA_Infinity'],
        'Value': [global_ta, sa_lambda, sa_inf]
    })
    summary_path = os.path.join(args.output, 'results_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"  Summary saved to: {summary_path}")
    
    # By-Region CSV
    region_path = os.path.join(args.output, 'results_by_region.csv')
    ta_df.to_csv(region_path, index=False)
    print(f"  Region stats saved to: {region_path}")
    
    # Optional: Save curve fit data for plotting later
    if bin_x is not None and len(bin_x) > 0:
        curve_df = pd.DataFrame({
            'Distance_mm': bin_x,
            'Mean_Correlation': bin_y,
            'Fitted_Model': fit_y if fit_y is not None else [np.nan]*len(bin_x)
        })
        curve_path = os.path.join(args.output, 'sa_curve_fit_data.csv')
        curve_df.to_csv(curve_path, index=False)
        print(f"  Curve fit data saved to: {curve_path}")
    
    print("\nAnalysis Complete! 🚀")

if __name__ == "__main__":
    main()