#!/usr/bin/env python3
"""
Nucleosome Periodicity Analysis
=================================
Reproduces the fragment length periodicity analysis from Snyder et al. (2016) Cell.
Shows the ~167bp peak and ~10.4bp periodicity reflecting nucleosome positioning.
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.signal import savgol_filter
from scipy.fft import fft, fftfreq


# ─────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────

def load_fragment_lengths(tsv_dir, condition=None):
    """Load fragment lengths from TSV files, optionally filtering by condition."""
    files = glob.glob(os.path.join(tsv_dir, "*.tsv"))
    if not files:
        raise FileNotFoundError(f"No TSV files found in {tsv_dir}")

    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        if condition:
            df = df[df['condition'] == condition]
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(combined):,} fragments ({condition or 'all conditions'})")
    return combined


# ─────────────────────────────────────────────
# 2. COMPUTE FRAGMENT LENGTH DISTRIBUTION
# ─────────────────────────────────────────────

def compute_distribution(lengths, min_len=100, max_len=220):
    """Compute normalized fragment length distribution at 1bp resolution."""
    bins = np.arange(min_len, max_len + 1)
    counts, _ = np.histogram(lengths, bins=np.arange(min_len, max_len + 2))
    # Normalize to density
    density = counts / counts.sum()
    return bins, density


def smooth_distribution(density, window=11, polyorder=3):
    """Apply Savitzky-Golay smoothing to reveal periodicity."""
    return savgol_filter(density, window_length=window, polyorder=polyorder)


# ─────────────────────────────────────────────
# 3. FFT PERIODICITY ANALYSIS
# ─────────────────────────────────────────────

def compute_periodicity(density, min_len=100):
    """Use FFT to identify dominant periodicities in fragment length distribution."""
    # Subtract mean to remove DC component
    signal = density - density.mean()
    
    # Compute FFT
    n = len(signal)
    yf = np.abs(fft(signal))[:n//2]
    xf = fftfreq(n, d=1)[:n//2]  # frequency in cycles per bp
    
    # Convert to period (bp)
    with np.errstate(divide='ignore'):
        periods = 1.0 / xf[1:]  # skip DC component
    power = yf[1:]
    
    return periods, power


# ─────────────────────────────────────────────
# 4. PLOTS
# ─────────────────────────────────────────────

def plot_periodicity(tsv_dir, outdir):
    """
    Reproduce Figure 1B from Snyder et al. (2016).
    Shows fragment length distribution with nucleosome periodicity.
    """
    # Load data by condition
    healthy_df = load_fragment_lengths(tsv_dir, condition='healthy')
    cancer_df  = load_fragment_lengths(tsv_dir, condition='cancer')

    healthy_lengths = healthy_df['fragment_length'].values
    cancer_lengths  = cancer_df['fragment_length'].values

    # Compute distributions at 1bp resolution
    bins_h, density_h = compute_distribution(healthy_lengths)
    bins_c, density_c = compute_distribution(cancer_lengths)

    # Smooth
    smooth_h = smooth_distribution(density_h)
    smooth_c = smooth_distribution(density_c)

    # ── Figure 1: Fragment length distributions with periodicity ──
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    colors = {'healthy': '#2196F3', 'cancer': '#F44336'}

    for ax, bins, density, smooth, label, color in [
        (axes[0], bins_h, density_h, smooth_h, 'Healthy', colors['healthy']),
        (axes[1], bins_c, density_c, smooth_c, 'Cancer',  colors['cancer'])
    ]:
        # Raw distribution
        ax.bar(bins, density, width=1, alpha=0.3, color=color, label='Raw')
        # Smoothed overlay
        ax.plot(bins, smooth, color=color, lw=2, label='Smoothed')

        # Mark 167bp peak
        ax.axvline(x=167, color='gray', linestyle='--', alpha=0.7, label='167bp (chromatosome)')

        # Mark nucleosome periodicity peaks (~10bp spacing from 167bp)
        for peak in [157, 147, 137, 127, 117, 107]:
            ax.axvline(x=peak, color='green', linestyle=':', alpha=0.4, linewidth=0.8)
        for peak in [177, 187, 197, 207, 217]:
            ax.axvline(x=peak, color='green', linestyle=':', alpha=0.4, linewidth=0.8)

        ax.set_xlabel('Fragment Length (bp)', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title(f'cfDNA Fragment Length Distribution — {label} (n={len(healthy_lengths if label=="Healthy" else cancer_lengths):,} fragments)', fontsize=13)
        ax.legend(fontsize=10)
        ax.set_xlim(100, 220)

    plt.suptitle('Nucleosome Periodicity in cfDNA Fragment Lengths\nReproducing Snyder et al. (2016) Figure 1B', fontsize=14, y=1.01)
    plt.tight_layout()
    path = os.path.join(outdir, 'nucleosome_periodicity.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"Periodicity plot saved to {path}")
    plt.close()


def plot_fft(tsv_dir, outdir):
    """Plot FFT power spectrum to quantify periodicity."""
    healthy_df = load_fragment_lengths(tsv_dir, condition='healthy')
    cancer_df  = load_fragment_lengths(tsv_dir, condition='cancer')

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, df, label, color in [
        (axes[0], healthy_df, 'Healthy', '#2196F3'),
        (axes[1], cancer_df,  'Cancer',  '#F44336')
    ]:
        bins, density = compute_distribution(df['fragment_length'].values)
        periods, power = compute_periodicity(density)

        # Focus on biologically relevant periods (5-30bp)
        mask = (periods >= 5) & (periods <= 30)
        ax.plot(periods[mask], power[mask], color=color, lw=2)
        ax.axvline(x=10.4, color='gray', linestyle='--', alpha=0.7, label='10.4bp (helical pitch)')
        ax.set_xlabel('Period (bp)', fontsize=12)
        ax.set_ylabel('FFT Power', fontsize=12)
        ax.set_title(f'FFT Power Spectrum — {label}', fontsize=13)
        ax.legend(fontsize=10)

    plt.suptitle('Nucleosome Helical Pitch (~10.4bp) in cfDNA', fontsize=14)
    plt.tight_layout()
    path = os.path.join(outdir, 'fft_periodicity.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"FFT plot saved to {path}")
    plt.close()


def plot_comparison(tsv_dir, outdir):
    """Overlay healthy vs cancer at 1bp resolution."""
    healthy_df = load_fragment_lengths(tsv_dir, condition='healthy')
    cancer_df  = load_fragment_lengths(tsv_dir, condition='cancer')

    bins_h, density_h = compute_distribution(healthy_df['fragment_length'].values)
    bins_c, density_c = compute_distribution(cancer_df['fragment_length'].values)

    smooth_h = smooth_distribution(density_h)
    smooth_c = smooth_distribution(density_c)

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(bins_h, smooth_h, color='#2196F3', lw=2.5, label='Healthy')
    ax.plot(bins_c, smooth_c, color='#F44336', lw=2.5, label='Cancer')
    ax.fill_between(bins_h, smooth_h, smooth_c,
                    where=(smooth_h > smooth_c),
                    alpha=0.15, color='#2196F3', label='Healthy enriched')
    ax.fill_between(bins_h, smooth_h, smooth_c,
                    where=(smooth_c > smooth_h),
                    alpha=0.15, color='#F44336', label='Cancer enriched')

    ax.axvline(x=167, color='gray', linestyle='--', alpha=0.7, label='167bp peak')
    ax.set_xlabel('Fragment Length (bp)', fontsize=13)
    ax.set_ylabel('Density', fontsize=13)
    ax.set_title('cfDNA Fragment Length: Healthy vs Cancer\n1bp Resolution with Smoothing', fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(100, 220)

    plt.tight_layout()
    path = os.path.join(outdir, 'periodicity_comparison.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"Comparison plot saved to {path}")
    plt.close()


# ─────────────────────────────────────────────
# 5. MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Nucleosome periodicity analysis of cfDNA fragment lengths'
    )
    parser.add_argument('--tsv_dir', required=True,
                        help='Directory containing fragment length TSV files')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for plots')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Generating nucleosome periodicity plots...")
    plot_periodicity(args.tsv_dir, args.outdir)
    plot_fft(args.tsv_dir, args.outdir)
    plot_comparison(args.tsv_dir, args.outdir)
    print("\nDone!")


if __name__ == '__main__':
    main()