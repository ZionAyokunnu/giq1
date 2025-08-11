#!/usr/bin/env python3
"""
Linearised BUSCO dotplotter:
- Loads two BUSCO TSVs (busco_id, sequence, gene_start, gene_end, strand, ...)
- Infers or ingests chromosome sizes
- Linearises per-chromosome coordinates across the whole genome
- Draws a dotplot with chromosome boundary lines and labels
- Robust when there are 0 common genes (saves an informative plot instead of crashing)

Usage:
python3 dotplotter.py \
  agora_ancestral_genome.tsv \
  giq_ancestral_genome.tsv \
  agora_vs_giq_comparison.png \
  --name1 "AGORA Ancestral" \
  --name2 "GIQ Ancestral" \
  --linearize \
  --sizes1 optional_sizes_genome1.csv \
  --sizes2 optional_sizes_genome2.csv

CSV for sizes (optional) should have columns:
  chromosome, chromsome_size_b  (typo accepted), or chromosome_size_b
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Loading + utilities
# -----------------------------
def load_busco_tsv(file_path: str, name: str) -> pd.DataFrame:
    print(f"Loading {name} from: {file_path}")
    df = pd.read_csv(file_path, sep='\t')
    print(f"  Loaded {len(df)} genes")

    required_cols = ['busco_id', 'sequence', 'gene_start', 'gene_end']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name}: missing required columns: {missing}")

    # Sort within each chromosome by start
    df = df.sort_values(['sequence', 'gene_start']).reset_index(drop=True)

    # Add simple per-chromosome rank position (backup plotting index)
    df['rank_in_chr'] = df.groupby('sequence').cumcount()

    # Report per-chromosome counts
    for chrom, n in df['sequence'].value_counts().sort_index().items():
        print(f"  {chrom}: {n} genes")

    return df


def _pick_size_column(df: pd.DataFrame) -> str:
    for col in ['chromsome_size_b', 'chromosome_size_b', 'size', 'size_b']:
        if col in df.columns:
            return col
    raise ValueError("Chromosome sizes file must contain one of: "
                     "chromsome_size_b, chromosome_size_b, size, size_b")


def _natural_sort(seq_list):
    # Basic natural sort for names like CM085462.1, chr1, chr10, etc.
    import re
    def key(s):
        return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]
    return sorted(seq_list, key=key)


def infer_chromosome_sizes_from_busco(df: pd.DataFrame) -> pd.DataFrame:
    """
    Infer chromosome sizes from BUSCO gene_end maxima per 'sequence'.
    This needs no external sheet.
    """
    g = df.groupby('sequence')['gene_end'].max().reset_index()
    g.columns = ['chromosome', 'chromsome_size_b']  # tolerate the common typo
    return g


def load_chromosome_sizes_csv(path: str) -> pd.DataFrame:
    """
    Load a CSV/TSV with columns:
      chromosome, chromsome_size_b  (or chromosome_size_b / size_b / size)
    Any extra columns are ignored.
    """
    sep = '\t' if str(path).lower().endswith(('.tsv', '.tab')) else ','
    df = pd.read_csv(path, sep=sep)
    if 'chromosome' not in df.columns:
        raise ValueError("Sizes CSV must contain a 'chromosome' column.")
    size_col = _pick_size_column(df)
    out = df[['chromosome', size_col]].copy()
    out.columns = ['chromosome', 'chromsome_size_b']
    return out


def build_offsets(size_df: pd.DataFrame) -> dict:
    """
    Compute cumulative offsets per chromosome in natural order.
    Returns: dict {chromosome -> cumulative_offset_in_bases}
    """
    size_df = size_df.copy()
    size_df['chromosome'] = size_df['chromosome'].astype(str)
    size_df = size_df.set_index('chromosome').loc[_natural_sort(size_df['chromosome'].unique())].reset_index()

    offsets = {}
    cum = 0
    for _, row in size_df.iterrows():
        offsets[row['chromosome']] = cum
        cum += int(row['chromsome_size_b'])
    return offsets


# -----------------------------
# Linearisation helpers (your requested simple style)
# -----------------------------
def linearise_coordinates(df: pd.DataFrame, size_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert per-chromosome coordinates to genome-linearised coordinates.
    Assumes:
      df columns: sequence, gene_start, gene_end
      size_df columns: chromosome, chromsome_size_b
    """
    offsets = build_offsets(size_df)
    out = df.copy()
    out['linear_start'] = out.apply(lambda r: offsets.get(str(r['sequence']), 0) + int(r['gene_start']), axis=1)
    out['linear_end']   = out.apply(lambda r: offsets.get(str(r['sequence']), 0) + int(r['gene_end']),   axis=1)
    return out, size_df


def add_chromosome_boundaries(ax, size_df: pd.DataFrame, axis='x'):
    """
    Add boundary lines at cumulative chromosome ends for the given axis.
    """
    size_df = size_df.copy()
    size_df['chromosome'] = size_df['chromosome'].astype(str)
    size_df = size_df.set_index('chromosome').loc[_natural_sort(size_df['chromosome'].unique())].reset_index()

    cum = 0
    for i, (_, row) in enumerate(size_df.iterrows()):
        cum += int(row['chromsome_size_b'])
        if i < len(size_df) - 1:
            if axis == 'x':
                ax.axvline(x=cum, linewidth=0.8, alpha=0.7, color='grey')
            else:
                ax.axhline(y=cum, linewidth=0.8, alpha=0.7, color='grey')


def add_chromosome_labels(ax, size_df_x: pd.DataFrame, size_df_y: pd.DataFrame):
    """
    Add chromosome labels as secondary axes (top for X, right for Y).
    """
    # X labels
    s1 = size_df_x.copy()
    s1['chromosome'] = s1['chromosome'].astype(str)
    s1 = s1.set_index('chromosome').loc[_natural_sort(s1['chromosome'].unique())].reset_index()

    x_ticks, x_labels, cum = [], [], 0
    for _, r in s1.iterrows():
        mid = cum + int(r['chromsome_size_b']) / 2
        x_ticks.append(mid)
        x_labels.append(r['chromosome'])
        cum += int(r['chromsome_size_b'])

    # Y labels
    s2 = size_df_y.copy()
    s2['chromosome'] = s2['chromosome'].astype(str)
    s2 = s2.set_index('chromosome').loc[_natural_sort(s2['chromosome'].unique())].reset_index()

    y_ticks, y_labels, cum = [], [], 0
    for _, r in s2.iterrows():
        mid = cum + int(r['chromsome_size_b']) / 2
        y_ticks.append(mid)
        y_labels.append(r['chromosome'])
        cum += int(r['chromsome_size_b'])

    ax_top = ax.secondary_xaxis('top')
    ax_top.set_xticks(x_ticks)
    ax_top.set_xticklabels(x_labels, rotation=45, ha='left', fontsize=10)
    ax_top.tick_params(length=0)

    ax_right = ax.secondary_yaxis('right')
    ax_right.set_yticks(y_ticks)
    ax_right.set_yticklabels(y_labels, fontsize=10)
    ax_right.tick_params(length=0)


# -----------------------------
# Plotting
# -----------------------------
def create_linear_dotplot(genome1_df: pd.DataFrame,
                          genome2_df: pd.DataFrame,
                          genome1_name: str,
                          genome2_name: str,
                          output_path: str,
                          linearize: bool = False,
                          sizes1: pd.DataFrame | None = None,
                          sizes2: pd.DataFrame | None = None):
    """
    Build comparison plot. If linearize=True, uses linearised coordinates and draws
    chromosome boundaries/labels.
    """

    # Common BUSCO ids
    common = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"Common genes found: {len(common)}")

    # Early path for 0 commons (save info plot, avoid crash)
    if len(common) == 0:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.axis('off')
        ax.text(0.03, 0.95,
                f'NO COMMON GENES FOUND\n\n'
                f'{genome1_name}: {len(genome1_df)} genes, {genome1_df["sequence"].nunique()} chromosomes\n'
                f'{genome2_name}: {len(genome2_df)} genes, {genome2_df["sequence"].nunique()} chromosomes\n\n'
                f'Check BUSCO ID formats and extraction.',
                transform=ax.transAxes, va='top', fontsize=12, color='red')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved informative plot to: {output_path}")
        return {
            'total_genes_genome1': len(genome1_df),
            'total_genes_genome2': len(genome2_df),
            'common_genes': 0,
            'correlation': 0.0,
            'chromosomes_genome1': genome1_df['sequence'].nunique(),
            'chromosomes_genome2': genome2_df['sequence'].nunique()
        }

    # Filter to common genes in order of appearance
    g1 = genome1_df[genome1_df['busco_id'].isin(common)].copy()
    g2 = genome2_df[genome2_df['busco_id'].isin(common)].copy()

    # Derive plotting coordinates
    if linearize:
        # Get or infer chromosome sizes
        s1 = sizes1 if sizes1 is not None else infer_chromosome_sizes_from_busco(genome1_df)
        s2 = sizes2 if sizes2 is not None else infer_chromosome_sizes_from_busco(genome2_df)

        g1_lin, s1_used = linearise_coordinates(g1, s1)
        g2_lin, s2_used = linearise_coordinates(g2, s2)

        # Map by BUSCO id → linear_start
        x_pos = dict(zip(g1_lin['busco_id'], g1_lin['linear_start']))
        y_pos = dict(zip(g2_lin['busco_id'], g2_lin['linear_start']))
        coords = [(x_pos[b], y_pos[b]) for b in common if b in x_pos and b in y_pos]

        xs = np.array([c[0] for c in coords], dtype=float)
        ys = np.array([c[1] for c in coords], dtype=float)

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        ax = axes[0, 0]
        ax.scatter(xs, ys, s=10, alpha=0.7, edgecolors='none')

        # Diagonal + labels
        m = max(xs.max(), ys.max())
        ax.plot([0, m], [0, m], 'r--', alpha=0.4, linewidth=1.5, label='y=x')
        ax.set_xlabel(f'{genome1_name} (linearised, Mb)')
        ax.set_ylabel(f'{genome2_name} (linearised, Mb)')
        ax.set_title(f'BUSCO Linearised Gene Order: {genome1_name} vs {genome2_name}')

        # Format axes as Mb
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f}'))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y/1e6:.0f}'))

        # Boundaries + labels
        add_chromosome_boundaries(ax, s1_used, axis='x')
        add_chromosome_boundaries(ax, s2_used, axis='y')
        add_chromosome_labels(ax, s1_used, s2_used)

        # Side panels: counts per chromosome
        ax_chr1 = axes[0, 1]
        g1_counts = g1.groupby('sequence').size().reindex(_natural_sort(g1['sequence'].unique()), fill_value=0)
        ax_chr1.bar(range(len(g1_counts)), g1_counts.values, alpha=0.7)
        ax_chr1.set_xticks(range(len(g1_counts)))
        ax_chr1.set_xticklabels(g1_counts.index, rotation=45, ha='right')
        ax_chr1.set_title(f'{genome1_name} BUSCOs / chr')
        ax_chr1.set_ylabel('Count')

        ax_chr2 = axes[1, 0]
        g2_counts = g2.groupby('sequence').size().reindex(_natural_sort(g2['sequence'].unique()), fill_value=0)
        ax_chr2.bar(range(len(g2_counts)), g2_counts.values, alpha=0.7)
        ax_chr2.set_xticks(range(len(g2_counts)))
        ax_chr2.set_xticklabels(g2_counts.index, rotation=45, ha='right')
        ax_chr2.set_title(f'{genome2_name} BUSCOs / chr')
        ax_chr2.set_ylabel('Count')

        # Stats panel
        ax_stats = axes[1, 1]
        ax_stats.axis('off')
        corr = np.corrcoef(xs, ys)[0, 1] if len(xs) > 1 else 0.0
        txt = (
            f"COMPARISON SUMMARY\n\n"
            f"Total common BUSCOs: {len(common)}\n"
            f"{genome1_name}: {len(genome1_df)} genes, {genome1_df['sequence'].nunique()} chromosomes\n"
            f"{genome2_name}: {len(genome2_df)} genes, {genome2_df['sequence'].nunique()} chromosomes\n\n"
            f"Linear correlation: {corr:.3f}\n"
        )
        ax_stats.text(0.02, 0.98, txt, va='top', family='monospace')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Dotplot saved to: {output_path}")

        return {
            'total_genes_genome1': len(genome1_df),
            'total_genes_genome2': len(genome2_df),
            'common_genes': len(common),
            'correlation': float(corr),
            'chromosomes_genome1': genome1_df['sequence'].nunique(),
            'chromosomes_genome2': genome2_df['sequence'].nunique()
        }

    # Non-linearised fallback: rank-based positions (original approach)
    g1_positions = dict(zip(g1['busco_id'], g1['rank_in_chr']))
    g2_positions = dict(zip(g2['busco_id'], g2['rank_in_chr']))

    xs = np.array([g1_positions[b] for b in common if b in g1_positions], dtype=float)
    ys = np.array([g2_positions[b] for b in common if b in g2_positions], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax = axes[0, 0]
    ax.scatter(xs, ys, s=10, alpha=0.7, edgecolors='none')

    m = max(xs.max(), ys.max())
    ax.plot([0, m], [0, m], 'r--', alpha=0.4, linewidth=1.5, label='y=x')
    ax.set_xlabel(f'{genome1_name} rank')
    ax.set_ylabel(f'{genome2_name} rank')
    ax.set_title(f'BUSCO Gene Order (rank): {genome1_name} vs {genome2_name}')

    # Side panels: counts per chromosome
    ax_chr1 = axes[0, 1]
    g1_counts = g1['sequence'].value_counts().sort_index()
    ax_chr1.bar(range(len(g1_counts)), g1_counts.values, alpha=0.7)
    ax_chr1.set_xticks(range(len(g1_counts)))
    ax_chr1.set_xticklabels(g1_counts.index, rotation=45, ha='right')
    ax_chr1.set_title(f'{genome1_name} BUSCOs / chr')
    ax_chr1.set_ylabel('Count')

    ax_chr2 = axes[1, 0]
    g2_counts = g2['sequence'].value_counts().sort_index()
    ax_chr2.bar(range(len(g2_counts)), g2_counts.values, alpha=0.7)
    ax_chr2.set_xticks(range(len(g2_counts)))
    ax_chr2.set_xticklabels(g2_counts.index, rotation=45, ha='right')
    ax_chr2.set_title(f'{genome2_name} BUSCOs / chr')
    ax_chr2.set_ylabel('Count')

    # Stats panel
    ax_stats = axes[1, 1]
    ax_stats.axis('off')
    corr = np.corrcoef(xs, ys)[0, 1] if len(xs) > 1 else 0.0
    txt = (
        f"COMPARISON SUMMARY\n\n"
        f"Total common BUSCOs: {len(common)}\n"
        f"{genome1_name}: {len(genome1_df)} genes, {genome1_df['sequence'].nunique()} chromosomes\n"
        f"{genome2_name}: {len(genome2_df)} genes, {genome2_df['sequence'].nunique()} chromosomes\n\n"
        f"Rank correlation: {corr:.3f}\n"
    )
    ax_stats.text(0.02, 0.98, txt, va='top', family='monospace')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Dotplot saved to: {output_path}")

    return {
        'total_genes_genome1': len(genome1_df),
        'total_genes_genome2': len(genome2_df),
        'common_genes': len(common),
        'correlation': float(corr),
        'chromosomes_genome1': genome1_df['sequence'].nunique(),
        'chromosomes_genome2': genome2_df['sequence'].nunique()
    }


# -----------------------------
# CLI
# -----------------------------
def main():
    p = argparse.ArgumentParser(description='Linearised BUSCO dotplot comparison')
    p.add_argument('genome1_tsv')
    p.add_argument('genome2_tsv')
    p.add_argument('output_plot')
    p.add_argument('--name1', default='Genome1')
    p.add_argument('--name2', default='Genome2')
    p.add_argument('--linearize', action='store_true', help='Use genome-linear coordinates and draw chr boundaries')
    p.add_argument('--sizes1', help='Optional CSV/TSV for genome1 sizes (chromosome, chromsome_size_b)')
    p.add_argument('--sizes2', help='Optional CSV/TSV for genome2 sizes (chromosome, chromsome_size_b)')
    args = p.parse_args()

    # Load TSVs
    g1 = load_busco_tsv(args.genome1_tsv, args.name1)
    g2 = load_busco_tsv(args.genome2_tsv, args.name2)

    # Optional size tables
    s1 = load_chromosome_sizes_csv(args.sizes1) if args.sizes1 else None
    s2 = load_chromosome_sizes_csv(args.sizes2) if args.sizes2 else None

    stats = create_linear_dotplot(
        g1, g2, args.name1, args.name2, args.output_plot,
        linearize=args.linearize, sizes1=s1, sizes2=s2
    )

    print("\n" + "="*50)
    print("DOTPLOT COMPARISON COMPLETE")
    print("="*50)
    print(f"Common genes: {stats['common_genes']}")
    print(f"Correlation: {stats['correlation']:.3f}")
    print(f"Plot saved to: {args.output_plot}")


if __name__ == "__main__":
    main()

