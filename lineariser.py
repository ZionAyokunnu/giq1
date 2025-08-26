#!/usr/bin/env python3
"""
Independent script to linearize multi-chromosome BUSCO genome to single artificial chromosome.
Creates TSV output compatible with ancestral genome formats.

Usage:
python3 linearize_genome.py input_busco.tsv output_linearized.tsv [--chr-name chr1] [--bin-size 20000]

python3 lineariser.py \
    /Users/zionayokunnu/Documents/Giq/compare/root_giq_ancestral_ordinal.tsv \
    compare/root_giq_anc_ordinal_linearized.tsv \
    --chr-name chr1 \
    --bin-size 288000
"""

import pandas as pd
import argparse
import sys
from pathlib import Path


def load_busco_file(file_path):
    """Load BUSCO file handling both formats"""
    print(f"Loading BUSCO file: {file_path}")
    
    # Check if file has comment headers
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
    
    if first_line.startswith('#'):
        # BUSCO format with comments
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None, 
                         names=['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length'])
        
        # Convert positions to numeric
        df['gene_start'] = pd.to_numeric(df['gene_start'], errors='coerce')
        df['gene_end'] = pd.to_numeric(df['gene_end'], errors='coerce')
        
        # Filter to complete genes only
        df = df[df['status'] == 'Complete'].copy()
        print(f"  Filtered to {len(df)} complete genes")
    else:
        # Clean TSV format
        df = pd.read_csv(file_path, sep='\t')
        df['gene_start'] = pd.to_numeric(df['gene_start'], errors='coerce')
        df['gene_end'] = pd.to_numeric(df['gene_end'], errors='coerce')
    
    print(f"  Loaded {len(df)} genes from {len(df['sequence'].unique())} chromosomes")
    
    # Show chromosome distribution
    chr_counts = df['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        print(f"    {chr_name}: {count} genes")
    
    return df


def create_linearized_single_chromosome_genome(multi_chr_df, artificial_chr_name="chr1", bin_size=20000):
    """
    Convert multi-chromosome genome to single artificial chromosome
    while preserving gene order within each chromosome.
    
    Args:
        multi_chr_df: DataFrame with multiple chromosomes
        artificial_chr_name: Name for the artificial chromosome
        bin_size: Size of each gene bin in bp
        
    Returns:
        linearized_df: DataFrame with single artificial chromosome
    """
    
    print(f"\nLinearizing to single chromosome: {artificial_chr_name}")
    print(f"Using bin size: {bin_size:,} bp")
    
    linearized_df = multi_chr_df.copy()
    
    # Sort by chromosome name and position to preserve natural order
    linearized_df = linearized_df.sort_values(['sequence', 'gene_start'])
    print(f"  Sorted {len(linearized_df)} genes by chromosome and position")
    
    # Create sequential artificial positions
    artificial_positions = []
    artificial_ends = []
    cumulative_pos = 0
    
    # Track original chromosome for reference
    original_chromosomes = []
    
    current_chr = None
    chr_gene_count = 0
    
    for idx, row in linearized_df.iterrows():
        # Track chromosome transitions
        if row['sequence'] != current_chr:
            if current_chr is not None:
                print(f"    {current_chr}: {chr_gene_count} genes → positions {start_pos:,} to {cumulative_pos-bin_size:,}")
            current_chr = row['sequence']
            chr_gene_count = 0
            start_pos = cumulative_pos
        
        artificial_positions.append(cumulative_pos)
        artificial_ends.append(cumulative_pos + bin_size)
        original_chromosomes.append(row['sequence'])
        
        cumulative_pos += bin_size
        chr_gene_count += 1
    
    # Print final chromosome
    if current_chr is not None:
        print(f"    {current_chr}: {chr_gene_count} genes → positions {start_pos:,} to {cumulative_pos-bin_size:,}")
    
    # Update to single artificial chromosome
    linearized_df['original_sequence'] = linearized_df['sequence']  # Keep original for reference
    linearized_df['sequence'] = artificial_chr_name
    linearized_df['gene_start'] = artificial_positions
    linearized_df['gene_end'] = artificial_ends
    
    # Add linear position for compatibility
    linearized_df['linear_position'] = range(len(linearized_df))
    
    print(f"  Created linearized genome: {len(linearized_df)} genes")
    print(f"  Total span: 0 to {cumulative_pos-bin_size:,} bp ({(cumulative_pos-bin_size)/1e6:.1f} Mb)")
    
    return linearized_df


def save_linearized_tsv(linearized_df, output_path):
    """Save linearized genome in TSV format compatible with ancestral genomes"""
    
    print(f"\nSaving linearized genome to: {output_path}")
    
    # Create output in standard BUSCO format
    output_df = linearized_df[['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length']].copy()
    
    # Ensure status is Complete
    output_df['status'] = 'Complete'
    
    # Save with tab separation
    output_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"  Saved {len(output_df)} genes")
    print(f"  Format: Standard BUSCO TSV")
    print(f"  Chromosome: {output_df['sequence'].iloc[0]}")
    print(f"  Position range: {output_df['gene_start'].min():,} to {output_df['gene_end'].max():,}")
    
    return output_path


def create_comparison_summary(original_df, linearized_df, output_dir):
    """Create summary comparison between original and linearized genomes"""
    
    summary_path = Path(output_dir) / "linearization_summary.txt"
    
    summary_text = f"""
GENOME LINEARIZATION SUMMARY
{"="*50}

ORIGINAL GENOME:
Total genes: {len(original_df):,}
Chromosomes: {len(original_df['sequence'].unique())}
Chromosome distribution:
"""
    
    chr_counts = original_df['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        summary_text += f"  {chr_name}: {count} genes\n"
    
    summary_text += f"""

LINEARIZED GENOME:
Total genes: {len(linearized_df):,}
Chromosomes: {len(linearized_df['sequence'].unique())} ({linearized_df['sequence'].iloc[0]})
Position range: {linearized_df['gene_start'].min():,} to {linearized_df['gene_end'].max():,} bp
Total span: {(linearized_df['gene_end'].max()/1e6):.1f} Mb

VALIDATION:
Gene count preserved: {len(original_df) == len(linearized_df)}
All genes present: {set(original_df['busco_id']) == set(linearized_df['busco_id'])}

NEXT STEPS:
1. Test linearization: Plot original vs linearized (expect perfect diagonal)
2. Test segmentation: Use linearized as "AGORA-type" genome
3. Compare results: Should match direct genome-to-genome comparison
"""
    
    with open(summary_path, 'w') as f:
        f.write(summary_text)
    
    print(f"Summary saved to: {summary_path}")
    return summary_path


def main():
    parser = argparse.ArgumentParser(description='Linearize multi-chromosome BUSCO genome to single artificial chromosome')
    parser.add_argument('input_busco', help='Input BUSCO TSV file')
    parser.add_argument('output_linearized', help='Output linearized TSV file')
    parser.add_argument('--chr-name', default='chr1', help='Name for artificial chromosome (default: chr1)')
    parser.add_argument('--bin-size', type=int, default=20000, help='Gene bin size in bp (default: 20000)')
    
    args = parser.parse_args()
    
    try:
        print("GENOME LINEARIZATION TOOL")
        print("="*50)
        
        # Load original genome
        original_df = load_busco_file(args.input_busco)
        
        # Check if already single chromosome
        if len(original_df['sequence'].unique()) == 1:
            print(f"\nWarning: Input genome already has single chromosome ({original_df['sequence'].iloc[0]})")
            print("Linearization may not be necessary, but proceeding...")
        
        # Create linearized version
        linearized_df = create_linearized_single_chromosome_genome(
            original_df, 
            args.chr_name, 
            args.bin_size
        )
        
        # Save linearized genome
        output_path = save_linearized_tsv(linearized_df, args.output_linearized)
        
        # Create summary
        output_dir = Path(args.output_linearized).parent
        summary_path = create_comparison_summary(original_df, linearized_df, output_dir)
        
        print("\n" + "="*50)
        print("LINEARIZATION COMPLETED SUCCESSFULLY")
        print("="*50)
        print(f"Original: {len(original_df)} genes, {len(original_df['sequence'].unique())} chromosomes")
        print(f"Linearized: {len(linearized_df)} genes, 1 chromosome ({args.chr_name})")
        print(f"Output: {output_path}")
        print(f"Summary: {summary_path}")
        
        print(f"\nTest commands:")
        print(f"# Test 1: Plot original vs linearized (expect diagonal)")
        print(f"python3 busco_dotplotter.py {args.input_busco} {args.output_linearized} test1_original_vs_linearized")
        print(f"")
        print(f"# Test 2: Use linearized as AGORA-type genome")
        print(f"python3 busco_dotplotter.py {args.output_linearized} [other_genome.tsv] test2_linearized_vs_other")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
    