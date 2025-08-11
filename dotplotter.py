#!/usr/bin/env python3
"""
Independent linear BUSCO dotplotter for comparing gene orders between two BUSCO TSV files.
"""

"""
script:

python3 busco_dotplotter.py \
  agora_ancestral_genome.tsv \
  giq_ancestral_genome.tsv \
  agora_vs_giq_comparison.png \
  --name1 "AGORA Ancestral" \
  --name2 "GIQ Ancestral"
  --linear-output agora_vs_giq_linearised.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import seaborn as sns

def load_busco_tsv(file_path: str, name: str):
    """Load and process BUSCO TSV file"""
    print(f"Loading {name} from: {file_path}")
    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"  Loaded {len(df)} genes")
        required_cols = ['busco_id', 'sequence', 'gene_start', 'gene_end']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        df = df.sort_values(['sequence', 'gene_start'])
        df['linear_position'] = 0
        cumulative_pos = 0
        for chrom in sorted(df['sequence'].unique()):
            chrom_mask = df['sequence'] == chrom
            chrom_genes = df[chrom_mask]
            positions = np.arange(cumulative_pos, cumulative_pos + len(chrom_genes))
            df.loc[chrom_mask, 'linear_position'] = positions
            cumulative_pos += len(chrom_genes)
            print(f"  {chrom}: {len(chrom_genes)} genes")
        return df
    except Exception as e:
        print(f"Error loading {name}: {e}")
        raise

def create_linear_dotplot(genome1_df, genome2_df, genome1_name, genome2_name, output_path):
    """Create linear dotplot comparing two genomes - Fixed to handle empty common genes"""
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"Common genes found: {len(common_genes)}")
    if len(common_genes) == 0:
        print("Warning: No common genes found between the two genomes!")
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        ax_main = axes[0, 0]
        ax_main.text(0.5, 0.5, 'NO COMMON GENES FOUND\nCheck BUSCO ID formats',
                    ha='center', va='center', transform=ax_main.transAxes,
                    fontsize=16, color='red', weight='bold')
        ax_main.set_xlabel(f'{genome1_name} Linear Position')
        ax_main.set_ylabel(f'{genome2_name} Linear Position')
        ax_main.set_title(f'Linear Gene Order Comparison\n{genome1_name} vs {genome2_name}')
        ax_chr1 = axes[0, 1]; ax_chr1.axis('off')
        sample1 = list(genome1_df['busco_id'].head(10))
        ax_chr1.text(0.05, 0.95, f'{genome1_name} Sample IDs:\n' + '\n'.join(sample1),
                    transform=ax_chr1.transAxes, fontsize=10, verticalalignment='top')
        ax_chr2 = axes[1, 0]; ax_chr2.axis('off')
        sample2 = list(genome2_df['busco_id'].head(10))
        ax_chr2.text(0.05, 0.95, f'{genome2_name} Sample IDs:\n' + '\n'.join(sample2),
                    transform=ax_chr2.transAxes, fontsize=10, verticalalignment='top')
        ax_stats = axes[1, 1]; ax_stats.axis('off')
        stats_text = f"""
ERROR: NO COMMON GENES FOUND

{genome1_name}: {len(genome1_df)} genes
{genome2_name}: {len(genome2_df)} genes

Possible causes:
• Different BUSCO ID formats
• Different gene sets
• Parsing errors

Check that both files contain 
valid BUSCO IDs like:
100019at7147, 77167at7147, etc.
"""
        ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes, fontsize=12, verticalalignment='top', color='red')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Error plot saved to: {output_path}")
        return {
            'total_genes_genome1': len(genome1_df),
            'total_genes_genome2': len(genome2_df),
            'common_genes': 0,
            'correlation': 0.0,
            'chromosomes_genome1': len(genome1_df['sequence'].unique()),
            'chromosomes_genome2': len(genome2_df['sequence'].unique())
        }
    g1_common = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2_common = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()
    g1_positions = dict(zip(g1_common['busco_id'], g1_common['linear_position']))
    g2_positions = dict(zip(g2_common['busco_id'], g2_common['linear_position']))
    plot_data = []
    for gene in common_genes:
        if gene in g1_positions and gene in g2_positions:
            plot_data.append({'gene': gene,'genome1_pos': g1_positions[gene],'genome2_pos': g2_positions[gene]})
    plot_df = pd.DataFrame(plot_data)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    ax_main = axes[0, 0]
    ax_main.scatter(plot_df['genome1_pos'], plot_df['genome2_pos'],
                    alpha=0.7, s=60, c='blue', edgecolors='darkblue', linewidth=0.5)
    ax_main.set_xlabel(f'{genome1_name} Linear Position')
    ax_main.set_ylabel(f'{genome2_name} Linear Position')
    ax_main.set_title(f'Linear Gene Order Comparison\n{genome1_name} vs {genome2_name}')
    ax_main.grid(True, alpha=0.3)
    max_pos = max(plot_df['genome1_pos'].max(), plot_df['genome2_pos'].max())
    ax_main.plot([0, max_pos], [0, max_pos], 'r--', alpha=0.5, linewidth=2, label='Perfect correlation')
    ax_main.legend()
    correlation = plot_df['genome1_pos'].corr(plot_df['genome2_pos'])
    ax_main.text(0.05, 0.95, f'Correlation: {correlation:.3f}',
                transform=ax_main.transAxes, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    ax_chr1 = axes[0, 1]
    chr_counts1 = g1_common['sequence'].value_counts()
    ax_chr1.bar(range(len(chr_counts1)), chr_counts1.values, alpha=0.7, color='lightblue')
    ax_chr1.set_xticks(range(len(chr_counts1)))
    ax_chr1.set_xticklabels(chr_counts1.index, rotation=45, ha='right')
    ax_chr1.set_ylabel('Gene Count')
    ax_chr1.set_title(f'{genome1_name} Chromosome Distribution')
    ax_chr2 = axes[1, 0]
    chr_counts2 = g2_common['sequence'].value_counts()
    ax_chr2.bar(range(len(chr_counts2)), chr_counts2.values, alpha=0.7, color='lightgreen')
    ax_chr2.set_xticks(range(len(chr_counts2)))
    ax_chr2.set_xticklabels(chr_counts2.index, rotation=45, ha='right')
    ax_chr2.set_ylabel('Gene Count')
    ax_chr2.set_title(f'{genome2_name} Chromosome Distribution')
    ax_stats = axes[1, 1]; ax_stats.axis('off')
    stats_text = f"""
COMPARISON SUMMARY

Total genes compared: {len(plot_data)}
{genome1_name}: {len(genome1_df)} genes, {len(genome1_df['sequence'].unique())} chromosomes
{genome2_name}: {len(genome2_df)} genes, {len(genome2_df['sequence'].unique())} chromosomes

Linear correlation: {correlation:.3f}

INTERPRETATION:
• Correlation > 0.7: High synteny conservation
• Correlation 0.3-0.7: Moderate synteny
• Correlation < 0.3: Low synteny / high rearrangement

COMMON GENES BY CHROMOSOME:
"""
    for chr1 in sorted(g1_common['sequence'].unique()):
        chr1_genes = set(g1_common[g1_common['sequence'] == chr1]['busco_id'])
        for chr2 in sorted(g2_common['sequence'].unique()):
            chr2_genes = set(g2_common[g2_common['sequence'] == chr2]['busco_id'])
            overlap = len(chr1_genes & chr2_genes)
            if overlap > 0:
                stats_text += f"\n{chr1} ↔ {chr2}: {overlap} genes"
    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes,
                 fontsize=9, verticalalignment='top', fontfamily='monospace')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Dotplot saved to: {output_path}")
    return {
        'total_genes_genome1': len(genome1_df),
        'total_genes_genome2': len(genome2_df),
        'common_genes': len(common_genes),
        'correlation': correlation,
        'chromosomes_genome1': len(genome1_df['sequence'].unique()),
        'chromosomes_genome2': len(genome2_df['sequence'].unique())
    }

# -------------------------------
# NEW FUNCTION (single addition)
# -------------------------------
def create_linearised_dotplot_single(genome1_df, genome2_df, genome1_name, genome2_name, linear_output_path):
    """
    Create a SEPARATE PNG: linearised genome comparison using Google Sheet chromosome sizes.
    Does NOT alter the original plot or CLI behavior.
    """
    # 1) Intersect BUSCOs
    common_genes = set(genome1_df['busco_id']) & set(genome2_df['busco_id'])
    print(f"[Linearised] Common genes found: {len(common_genes)}")

    # 2) Filter to common
    g1 = genome1_df[genome1_df['busco_id'].isin(common_genes)].copy()
    g2 = genome2_df[genome2_df['busco_id'].isin(common_genes)].copy()

    # 3) Handle no-overlap case: still produce a PNG
    if len(g1) == 0 or len(g2) == 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, 'NO COMMON GENES FOUND\n(Linearised view)',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=16, color='red', weight='bold')
        ax.set_axis_off()
        fig.savefig(linear_output_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"[Linearised] Error plot saved to: {linear_output_path}")
        return

    # 4) Fetch chromosome sizes from Google Sheet
    sheet_url = "https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825"
    csv_url = sheet_url.replace('/edit?gid=', '/export?format=csv&gid=')
    try:
        sheet_df = pd.read_csv(csv_url)
        # Expect columns: species, chromosome, chromsome_size_b
        if not {'chromosome', 'chromsome_size_b'}.issubset(set(sheet_df.columns)):
            print("[Linearised] Sheet columns missing; falling back to infer sizes.")
            sheet_df = None
    except Exception as e:
        print(f"[Linearised] Failed to load sheet: {e}; falling back to infer sizes.")
        sheet_df = None

    # 5) Build size dicts per genome
    def build_sizes(df, sheet):
        sizes = {}
        if sheet is not None:
            # Map by chromosome name
            size_map = dict(zip(sheet['chromosome'].astype(str), sheet['chromsome_size_b'].astype(int)))
            for chrom in df['sequence'].unique():
                if str(chrom) in size_map:
                    sizes[chrom] = int(size_map[str(chrom)])
        # Fallback for any missing chromosomes: infer from gene_end
        missing = [c for c in df['sequence'].unique() if c not in sizes]
        if missing:
            inferred = df.groupby('sequence')['gene_end'].max().astype(int).to_dict()
            for c in missing:
                sizes[c] = inferred.get(c, 0)
        # Order by chromosome label (stable, minimal assumption)
        sizes = dict(sorted(sizes.items(), key=lambda x: str(x[0])))
        return sizes

    sizes1 = build_sizes(g1, sheet_df)
    sizes2 = build_sizes(g2, sheet_df)

    # 6) Offsets
    def offsets_from_sizes(sizes):
        offs, cum = {}, 0
        for chrom, sz in sizes.items():
            offs[chrom] = cum
            cum += int(sz)
        return offs

    off1 = offsets_from_sizes(sizes1)
    off2 = offsets_from_sizes(sizes2)

    # 7) Linear positions
    g1['start'] = g1['gene_start'].astype(int)
    g2['start'] = g2['gene_start'].astype(int)
    g1['linear_bp'] = g1.apply(lambda r: off1.get(r['sequence'], 0) + r['start'], axis=1)
    g2['linear_bp'] = g2.apply(lambda r: off2.get(r['sequence'], 0) + r['start'], axis=1)

    # 8) Build plotting table
    p1 = dict(zip(g1['busco_id'], g1['linear_bp']))
    p2 = dict(zip(g2['busco_id'], g2['linear_bp']))
    pts = [{'x': p1[b], 'y': p2[b]} for b in common_genes if b in p1 and b in p2]
    plot_df = pd.DataFrame(pts)

    # 9) Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    if len(plot_df) > 0:
        ax.scatter(plot_df['x'], plot_df['y'], alpha=0.6, s=30, c='red', edgecolors='darkred', linewidth=0.3)
        # Boundaries
        for chrom, off in list(off1.items())[1:]:
            ax.axvline(off, color='grey', linewidth=0.8, alpha=0.5, linestyle='--')
        for chrom, off in list(off2.items())[1:]:
            ax.axhline(off, color='grey', linewidth=0.8, alpha=0.5, linestyle='--')
        # Mb formatting
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f}'))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, p: f'{y/1e6:.0f}'))
        ax.set_xlabel(f'{genome1_name} Linearised Position (Mb)')
        ax.set_ylabel(f'{genome2_name} Linearised Position (Mb)')
        ax.set_title('Linearised Genome Comparison (separate PNG)')
        ax.grid(True, alpha=0.2)
        # Correlation
        if len(plot_df) > 1:
            R = plot_df['x'].corr(plot_df['y'])
            ax.text(0.05, 0.95, f'Linear R = {R:.3f}', transform=ax.transAxes,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.8))
    else:
        ax.text(0.5, 0.5, 'NO PLOT POINTS AFTER LINEARISATION',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=14, color='red', weight='bold')
        ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(linear_output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"[Linearised] Plot saved to: {linear_output_path}")

def main():
    parser = argparse.ArgumentParser(description='Create linear BUSCO dotplot comparison')
    parser.add_argument('genome1_tsv', help='First genome BUSCO TSV file')
    parser.add_argument('genome2_tsv', help='Second genome BUSCO TSV file')
    parser.add_argument('output_plot', help='Output plot file (PNG/PDF)')
    parser.add_argument('--name1', default='Genome1', help='Name for first genome')
    parser.add_argument('--name2', default='Genome2', help='Name for second genome')
    # NEW: optional separate PNG for linearised view
    parser.add_argument('--linear-output', dest='linear_output', default=None,
                        help='Optional: path to save a SEPARATE PNG with linearised genome comparison')

    args = parser.parse_args()
    try:
        genome1_df = load_busco_tsv(args.genome1_tsv, args.name1)
        genome2_df = load_busco_tsv(args.genome2_tsv, args.name2)

        # Original plot (unchanged behavior)
        stats = create_linear_dotplot(genome1_df, genome2_df,
                                      args.name1, args.name2,
                                      args.output_plot)

        print("\n" + "="*50)
        print("DOTPLOT COMPARISON COMPLETE")
        print("="*50)
        print(f"Common genes: {stats['common_genes']}")
        print(f"Linear correlation: {stats['correlation']:.3f}")
        print(f"Plot saved to: {args.output_plot}")

        # NEW: optional second PNG for linearised genome comparison
        if args.linear_output:
            create_linearised_dotplot_single(genome1_df, genome2_df,
                                             args.name1, args.name2,
                                             args.linear_output)
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
