#!/usr/bin/env python3
"""
cfDNA Fragment Length ML Classifier
=====================================
Trains a classifier to distinguish cancer from healthy samples
based on fragment length distributions from shallow WGS cfDNA data.

Based on: Snyder et al. (2016) Cell. GSE71378.
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import roc_curve, auc, classification_report
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# ─────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────

def load_fragment_data(tsv_dir):
    """Load all fragment length TSV files into a single DataFrame."""
    files = glob.glob(os.path.join(tsv_dir, "*.tsv"))
    if not files:
        raise FileNotFoundError(f"No TSV files found in {tsv_dir}")
    
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        dfs.append(df)
    
    return pd.concat(dfs, ignore_index=True)


# ─────────────────────────────────────────────
# 2. FEATURE ENGINEERING
# ─────────────────────────────────────────────

def engineer_features(df):
    """
    Convert raw fragment lengths into per-sample features.
    
    Features:
    - short_ratio: fraction of fragments in 100-150bp (enriched in ctDNA)
    - long_ratio: fraction of fragments in 151-220bp (healthy cfDNA)
    - short_long_ratio: ratio of short to long fragments (key cancer biomarker)
    - mean_length: mean fragment length
    - std_length: standard deviation of fragment length
    - median_length: median fragment length
    """
    records = []
    
    for sample, group in df.groupby('sample'):
        lengths = group['fragment_length'].values
        condition = group['condition'].iloc[0]
        
        short = np.sum((lengths >= 100) & (lengths <= 150)) / len(lengths)
        long  = np.sum((lengths >= 151) & (lengths <= 220)) / len(lengths)
        
        records.append({
            'sample':           sample,
            'condition':        condition,
            'short_ratio':      short,
            'long_ratio':       long,
            'short_long_ratio': short / long if long > 0 else 0,
            'mean_length':      np.mean(lengths),
            'std_length':       np.std(lengths),
            'median_length':    np.median(lengths)
        })
    
    return pd.DataFrame(records)


# ─────────────────────────────────────────────
# 3. TRAIN AND EVALUATE
# ─────────────────────────────────────────────

def train_and_evaluate(features_df, outdir):
    """Train classifier and evaluate with leave-one-out cross validation."""

    # Filter to only healthy and cancer for binary classification
    features_df = features_df[features_df['condition'].isin(['healthy', 'cancer'])].copy()
    print(f"\nFiltered to {len(features_df)} samples (healthy + cancer only)")

    # Deduplicate — IH02/SRR2130051 and IH03/SRR2130052 are the same samples
    features_df = features_df[~features_df['sample'].isin(['IC17_liver', 'IC26_prostate', 'IH02_healthy', 'IH03_healthy'])].copy()

    feature_cols = [
        'short_ratio', 'long_ratio', 'short_long_ratio',
        'mean_length', 'std_length', 'median_length'
    ]

    X = features_df[feature_cols].values
    y = (features_df['condition'] == 'cancer').astype(int).values
    samples = features_df['sample'].values

    # Scale features and use pipeline for clean cross-validation
    clf = Pipeline([
        ('scaler', StandardScaler()),
        ('clf', LogisticRegression(random_state=42, max_iter=1000))
    ])

    # Leave-one-out CV (appropriate for small sample sizes)
    loo = LeaveOneOut()
    y_prob = cross_val_predict(clf, X, y, cv=loo, method='predict_proba')[:, 1]
    y_pred = (y_prob >= 0.5).astype(int)

    # Print classification report
    print("\n=== Classification Report ===")
    print(classification_report(y, y_pred, target_names=['healthy', 'cancer']))

    # Save predictions
    results_df = pd.DataFrame({
        'sample':     samples,
        'condition':  features_df['condition'].values,
        'true_label': y,
        'pred_prob':  y_prob,
        'pred_label': y_pred
    })
    results_df.to_csv(os.path.join(outdir, 'predictions.csv'), index=False)
    print(f"\nPredictions saved to {outdir}/predictions.csv")

    return y, y_prob, features_df


# ─────────────────────────────────────────────
# 4. PLOTS
# ─────────────────────────────────────────────

def plot_fragment_distributions(df, outdir):
    """Plot fragment length distributions for cancer vs healthy."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = {'healthy': '#2196F3', 'cancer': '#F44336', 'other': '#4CAF50'}
    
    for condition, group in df.groupby('condition'):
        ax.hist(
            group['fragment_length'],
            bins=60, alpha=0.6, density=True,
            color=colors[condition],
            label=condition.capitalize()
        )
    
    ax.axvspan(100, 150, alpha=0.1, color='red', label='Short fragment zone (100-150bp)')
    ax.set_xlabel('Fragment Length (bp)', fontsize=13)
    ax.set_ylabel('Density', fontsize=13)
    ax.set_title('cfDNA Fragment Length Distributions: Cancer vs Healthy', fontsize=14)
    ax.legend(fontsize=11)
    
    plt.tight_layout()
    path = os.path.join(outdir, 'fragment_length_distributions.png')
    plt.savefig(path, dpi=150)
    print(f"Fragment distribution plot saved to {path}")
    plt.close()


def plot_roc_curve(y_true, y_prob, outdir):
    """Plot ROC curve."""
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc_auc = auc(fpr, tpr)
    
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(fpr, tpr, color='#E91E63', lw=2,
            label=f'ROC curve (AUC = {roc_auc:.2f})')
    ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
    ax.set_xlabel('False Positive Rate', fontsize=13)
    ax.set_ylabel('True Positive Rate', fontsize=13)
    ax.set_title('ROC Curve: Cancer vs Healthy Classification', fontsize=14)
    ax.legend(fontsize=11)
    
    plt.tight_layout()
    path = os.path.join(outdir, 'roc_curve.png')
    plt.savefig(path, dpi=150)
    print(f"ROC curve saved to {path}")
    plt.close()


def plot_feature_distributions(features_df, outdir):
    """Plot key features by condition."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    colors = {'healthy': '#2196F3', 'cancer': '#F44336', 'other': '#4CAF50'}
    features = ['short_long_ratio', 'mean_length', 'short_ratio']
    titles = ['Short/Long Ratio', 'Mean Fragment Length', 'Short Fragment Ratio']
    
    for ax, feat, title in zip(axes, features, titles):
        for condition, group in features_df.groupby('condition'):
            ax.bar(condition, group[feat].mean(),
                   color=colors[condition], alpha=0.8)
            ax.set_title(title, fontsize=12)
            ax.set_ylabel('Value', fontsize=11)
    
    plt.suptitle('Key Features by Condition', fontsize=14)
    plt.tight_layout()
    path = os.path.join(outdir, 'feature_distributions.png')
    plt.savefig(path, dpi=150)
    print(f"Feature distributions plot saved to {path}")
    plt.close()


# ─────────────────────────────────────────────
# 5. MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='cfDNA fragment length ML classifier'
    )
    parser.add_argument('--tsv_dir',  required=True,
                        help='Directory containing fragment length TSV files')
    parser.add_argument('--outdir',   required=True,
                        help='Output directory for plots and results')
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print("Loading fragment length data...")
    df = load_fragment_data(args.tsv_dir)
    print(f"Loaded {len(df):,} fragment length measurements across {df['sample'].nunique()} samples")
    
    print("\nEngineering features...")
    features_df = engineer_features(df)
    print(features_df[['sample', 'condition', 'short_long_ratio', 'mean_length']].to_string())
    
    print("\nTraining classifier...")
    y_true, y_prob, features_df = train_and_evaluate(features_df, args.outdir)
    
    print("\nGenerating plots...")
    plot_fragment_distributions(df, args.outdir)
    plot_roc_curve(y_true, y_prob, args.outdir)
    plot_feature_distributions(features_df, args.outdir)
    
    print("\nDone!")


if __name__ == '__main__':
    main()