#!/usr/bin/env python3
"""
Plot distribution samples from the simulated annealing distribution tests.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Configure matplotlib to use serif fonts
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman']
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

output_ext = '.pdf'  # for output plots

# Read the CSV file
csv_file = 'distribution_samples.csv'
if not Path(csv_file).exists():
    print(f"Error: {csv_file} not found!")
    print("Please run 'fpm test distribution_tests' first to generate the data.")
    exit(1)

df = pd.read_csv(csv_file)

# Get distribution names (all columns except 'sample')
dist_names = [col for col in df.columns if col != 'sample']
n_dists = len(dist_names)

# Convert all distribution columns to numeric, handling any formatting issues
for col in dist_names:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Drop any rows with NaN values that might have resulted from conversion
df = df.dropna()

print(f"Loaded {len(df)} valid samples for {n_dists} distributions")
print()

# Create a figure with subplots for each distribution
fig, axes = plt.subplots(1, 5, figsize=(20, 8))
fig.suptitle('Distribution Samples from Simulated Annealing Perturbations',
             fontsize=16, fontweight='bold')

axes = axes.flatten()

# Plot histogram for each distribution
for i, dist_name in enumerate(dist_names):
    ax = axes[i]
    data = df[dist_name].to_numpy()

    # Create histogram
    n, bins, patches = ax.hist(data, bins=50, density=True, alpha=0.7,
                                color=f'C{i}', edgecolor='black', linewidth=0.5)

    # Add statistics text
    mean_val = np.mean(data)
    std_val = np.std(data)
    min_val = np.min(data)
    max_val = np.max(data)

    stats_text = f'Mean: {mean_val:.2f}\nStd: {std_val:.2f}\nMin: {min_val:.2f}\nMax: {max_val:.2f}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_title(dist_name.replace('_', ' ').title(), fontweight='bold')
    ax.set_xlabel('Value')
    ax.set_ylabel('Probability Density')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'distribution_histograms{output_ext}', dpi=300, bbox_inches='tight')
print(f"Saved histogram plot to: distribution_histograms{output_ext}")

# Create a second figure showing time series (samples over time)
fig2, axes2 = plt.subplots(1, 5, figsize=(20, 8))
fig2.suptitle('Distribution Samples Time History',
              fontsize=16, fontweight='bold')

axes2 = axes2.flatten()

for i, dist_name in enumerate(dist_names):
    ax = axes2[i]
    data = df[dist_name].to_numpy()
    sample_nums = df['sample'].to_numpy()

    # Plot time series
    ax.plot(sample_nums, data, alpha=0.6, linewidth=0.5, color=f'C{i}')
    ax.axhline(y=np.mean(data), color='red', linestyle='--',
               linewidth=1.5, label=f'Mean: {np.mean(data):.2f}')

    ax.set_title(dist_name.replace('_', ' ').title(), fontweight='bold')
    ax.set_xlabel('Sample Number')
    ax.set_ylabel('Value')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([-5.5, 5.5])  # Consistent y-axis for comparison

plt.tight_layout()
plt.savefig(f'distribution_timeseries{output_ext}', dpi=300, bbox_inches='tight')
print(f"Saved time series plot to: distribution_timeseries{output_ext}")

# Create a comparison plot - all distributions on one plot
fig3, ax3 = plt.subplots(figsize=(14, 8))

for i, dist_name in enumerate(dist_names):
    data = df[dist_name].to_numpy()
    ax3.hist(data, bins=50, alpha=0.5, label=dist_name.replace('_', ' '),
             density=True, histtype='stepfilled')

ax3.set_title('All Distributions Comparison', fontsize=16, fontweight='bold')
ax3.set_xlabel('Value', fontsize=12)
ax3.set_ylabel('Probability Density', fontsize=12)
ax3.legend(loc='upper right', fontsize=10)
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'distribution_comparison{output_ext}', dpi=300, bbox_inches='tight')
print(f"Saved comparison plot to: distribution_comparison{output_ext}")

# Create box plot comparison
fig4, ax4 = plt.subplots(figsize=(14, 6))

data_for_boxplot = [df[dist_name].to_numpy() for dist_name in dist_names]
labels_for_boxplot = [name.replace('_', '\n') for name in dist_names]

bp = ax4.boxplot(data_for_boxplot, labels=labels_for_boxplot,
                  patch_artist=True, showmeans=True)

# Color the boxes
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(f'C{i}')
    patch.set_alpha(0.6)

ax4.set_title('Distribution Comparison - Box Plots', fontsize=16, fontweight='bold')
ax4.set_ylabel('Value', fontsize=12)
ax4.grid(True, alpha=0.3, axis='y')
ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

plt.tight_layout()
plt.savefig(f'distribution_boxplots{output_ext}', dpi=300, bbox_inches='tight')
print(f"Saved box plot to: distribution_boxplots{output_ext}")

# Print summary statistics
print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)
for dist_name in dist_names:
    data = df[dist_name].to_numpy()
    print(f"\n{dist_name.upper()}:")
    print(f"  Mean:     {np.mean(data):8.4f}")
    print(f"  Std Dev:  {np.std(data):8.4f}")
    print(f"  Min:      {np.min(data):8.4f}")
    print(f"  Max:      {np.max(data):8.4f}")
    print(f"  Median:   {np.median(data):8.4f}")
    print(f"  Q1:       {np.percentile(data, 25):8.4f}")
    print(f"  Q3:       {np.percentile(data, 75):8.4f}")

print("\n" + "="*60)
print("Plots saved successfully!")
print("="*60)
