"""
Chromatin Scaling Analysis Pipeline
Author: Amit Kumar Mondal
Description: Parses Hi-C files, normalizes contact matrices, and calculates 
             multi-scale polymer scaling exponents to analyze chromatin topology.
"""

import cooler
import matplotlib.pyplot as plt
import numpy as np


def analyze_polymer_scaling(cool_file_path, region="chr1:0-10000000"):
  """Loads balanced matrix and calculates multi-scale polymer scaling exponents."""
  c = cooler.Cooler(cool_file_path)
  resolution = c.binsize

  # Fetch Knight-Ruiz normalized matrix
  print(f"[*] Fetching balanced matrix for region: {region}...")
  mat_balanced = c.matrix(balance="KR").fetch(region)

  # Extract diagonals to calculate contact probability vs. genomic distance
  n = mat_balanced.shape[0]
  distances = []
  average_contacts = []

  for d in range(1, n):
    diagonal_values = np.diag(mat_balanced, k=d)
    # Safely filter out infinite values and NaNs
    clean_diagonal = diagonal_values[np.isfinite(diagonal_values)]

    if len(clean_diagonal) > 0:
      mean_val = np.mean(clean_diagonal)
      if mean_val > 0:
        distances.append(d * resolution)
        average_contacts.append(mean_val)

  # Convert list variables to NumPy arrays for mathematical slicing
  distances_arr = np.array(distances)
  log_s = np.log10(distances_arr)
  log_P = np.log10(average_contacts)

  # Global linear fit
  global_slope, global_intercept = np.polyfit(log_s, log_P, 1)

  # --- Multi-scale segmentation analysis ---
  # Short Range (Fractal Window: <= 2.5 Mb)
  short_mask = distances_arr <= 2500000
  log_s_short, log_P_short = log_s[short_mask], log_P[short_mask]
  slope_short = (
      np.polyfit(log_s_short, log_P_short, 1)[0]
      if len(log_s_short) > 1
      else None
  )

  # Long Range (Polymer thermodynamic limit: >= 3.0 Mb)
  long_mask = distances_arr >= 3000000
  log_s_long, log_P_long = log_s[long_mask], log_P[long_mask]
  slope_long = (
      np.polyfit(log_s_long, log_P_long, 1)[0] if len(log_s_long) > 1 else None
  )

  # Pack results
  results = {
      "matrix": mat_balanced,
      "distances": distances_arr,
      "log_s": log_s,
      "log_P": log_P,
      "global_slope": global_slope,
      "global_intercept": global_intercept,
      "short_slope": slope_short,
      "long_slope": slope_long,
  }

  return results


def plot_results(results, region_label="chr1:0-10Mb"):
  """Plots both the spatial matrix heatmap and the log-log scaling curves."""
  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

  # Plot 1: Balanced Hi-C Heatmap
  im = ax1.imshow(np.log1p(results["matrix"]), cmap="YlOrRd")
  fig.colorbar(im, ax=ax1, label="log(Contacts + 1)")
  ax1.set_title(f"Normalized Hi-C Contact Map ({region_label})")

  # Plot 2: Log-Log Scaling Decay
  ax2.scatter(
      results["log_s"],
      results["log_P"],
      color="darkgreen",
      alpha=0.7,
      label="KR-Balanced data",
  )
  fit_line = (
      results["global_slope"] * results["log_s"] + results["global_intercept"]
  )
  ax2.plot(
      results["log_s"],
      fit_line,
      color="black",
      linestyle="--",
      label=f"Balanced Global Fit (Alpha = {results['global_slope']:.2f})",
  )

  ax2.set_xlabel("log10(Genomic Distance in bp)")
  ax2.set_ylabel("log10(Contact Probability)")
  ax2.set_title("Biophysical Scaling Law Decay")
  ax2.legend()
  ax2.grid(True, linestyle=":", alpha=0.6)

  plt.tight_layout()
  plt.show()


def run_pipeline(hic_path, chrom="1", resolution=500000, region_end=10000000):
  """Orchestrates the entire chromatin analysis pipeline."""
  import os

  cool_output = "temp_data.cool"

  # Convert formats if not already done
  if not os.path.exists(cool_output):
    print(f"[*] Starting conversion of {hic_path} at resolution {resolution}...")
    cmd = f"hic2cool convert {hic_path} {cool_output} -r {resolution}"
    os.system(cmd)

  c = cooler.Cooler(cool_output)
  first_chrom = c.chromnames[0]
  target_region = f"{first_chrom}:0-{region_end}"

  results = analyze_polymer_scaling(cool_output, target_region)

  print("\n=== BIOPHYSICAL ANALYSIS COMPLETE ===")
  print(f"Global Scaling Exponent (Alpha): {results['global_slope']:.2f}")
  print(f"Short-Range Slope (500kb-2.5Mb): {results['short_slope']:.2f}")
  print(f"Long-Range Slope (3.0Mb-10Mb):   {results['long_slope']:.2f}")

  plot_results(results, target_region)
