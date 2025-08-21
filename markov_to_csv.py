#!/usr/bin/env python3
"""
Extract GIQ Markov profiles to separate BUSCO TSV files for positional and ordinal conservation.
Creates two independent ancestral genomes for different types of analysis.

Usage:
python3 markov_to_csv.py \
    comparison_output/stages/6a_positional_profile.csv \
    comparison_output/stages/6b_ordinal_profile.csv \
    compare/root_giq_ancestral_positional.tsv \
    compare/root_giq_ancestral_ordinal.tsv \
    --min-genes 400 \
    --freq-threshold 0.5 \
    --bin-size 20 \
    --rank-spacing 20
"""

import pandas as pd
import argparse
from pathlib import Path

def extract_positional_profile(positional_csv, output_tsv, min_genes=500, freq_threshold=0.8, bin_size_kb=40):
    """
    Extract positional profile to BUSCO TSV format.
    
    Args:
        positional_csv: Path to positional profile CSV
        output_tsv: Path to output BUSCO TSV
        min_genes: Minimum genes per chromosome
        freq_threshold: Minimum frequency threshold (0.8 = 80% of genomes)
        bin_size_kb: Bin size in kb for position calculation
        
    Returns:
        busco_df: BUSCO-formatted DataFrame
    """
    
    print("=" * 60)
    print("EXTRACTING POSITIONAL PROFILE")
    print("=" * 60)
    
    # Read positional profile
    print(f"Reading: {positional_csv}")
    df = pd.read_csv(positional_csv)
    print(f"Loaded {len(df)} positional entries")
    
    # Check what columns we have
    print(f"Available columns: {list(df.columns)}")
    
    # Extract chromosome and bin information (handle different column names)
    if 'bin_or_rank_id' in df.columns:
        bin_col = 'bin_or_rank_id'
    elif 'bin_id' in df.columns:
        bin_col = 'bin_id'
    else:
        print("Error: Cannot find bin ID column")
        print(f"Available columns: {list(df.columns)}")
        return None
    
    df['chromosome'] = df[bin_col].str.extract(r'(.+?)_bin_')[0]
    df['bin_number'] = df[bin_col].str.extract(r'_bin_(\d+)')[0].astype(int)
    
    # Calculate frequency ratio if possible
    if 'genome_count' in df.columns and 'total_genomes' in df.columns:
        df['frequency_ratio'] = df['genome_count'] / df['total_genomes']
        print(f"Frequency threshold: {freq_threshold} ({freq_threshold*100}% of genomes)")
        high_freq = df[df['frequency_ratio'] >= freq_threshold].copy()
        print(f"High-frequency entries: {len(high_freq)}/{len(df)} ({len(high_freq)/len(df)*100:.1f}%)")
    else:
        print("No frequency data found - using all entries")
        high_freq = df.copy()
        high_freq['frequency_ratio'] = 1.0
    
    # Filter chromosomes by gene count
    genes_per_chr = high_freq.groupby('chromosome')['busco_id'].nunique().sort_values(ascending=False)
    major_chromosomes = genes_per_chr[genes_per_chr >= min_genes].index.tolist()
    
    print(f"Chromosome filtering (min {min_genes} genes):")
    print(f"  Total chromosomes: {len(genes_per_chr)}")
    print(f"  Major chromosomes: {len(major_chromosomes)}")
    
    for chr_name in major_chromosomes:
        print(f"    {chr_name}: {genes_per_chr[chr_name]} genes")
    
    # Filter to major chromosomes
    major_data = high_freq[high_freq['chromosome'].isin(major_chromosomes)].copy()
    print(f"Filtered to {len(major_data)} entries on major chromosomes")
    
    # For each gene on each chromosome, find the best bin (highest overlap)
    print("Finding best positional assignments...")
    best_assignments = []
    
    for chromosome in sorted(major_chromosomes):
        chr_data = major_data[major_data['chromosome'] == chromosome]
        
        for busco_id in chr_data['busco_id'].unique():
            gene_entries = chr_data[chr_data['busco_id'] == busco_id]
            
            # Find entry with highest overlap percentage
            if 'overlap_or_frequency' in gene_entries.columns:
                best_entry = gene_entries.loc[gene_entries['overlap_or_frequency'].idxmax()]
                overlap_score = best_entry['overlap_or_frequency']
            else:
                best_entry = gene_entries.iloc[0]
                overlap_score = 100.0
            
            best_assignments.append({
                'busco_id': busco_id,
                'chromosome': chromosome,
                'bin_number': best_entry['bin_number'],
                'gene_position': best_entry.get('gene_position', best_entry['bin_number'] * bin_size_kb * 1000),
                'overlap_score': overlap_score,
                'frequency_ratio': best_entry['frequency_ratio'],
                'genome_count': best_entry.get('genome_count', 1),
                'total_genomes': best_entry.get('total_genomes', 1)
            })
    
    best_df = pd.DataFrame(best_assignments)
    print(f"Best assignments: {len(best_df)} genes")
    
    # Convert to BUSCO format
    print("Converting to BUSCO format...")
    busco_entries = []
    
    for chromosome in sorted(best_df['chromosome'].unique()):
        chr_genes = best_df[best_df['chromosome'] == chromosome].copy()
        chr_genes = chr_genes.sort_values('bin_number')  # Sort by position
        
        print(f"  {chromosome}: {len(chr_genes)} genes")
        
        for _, gene in chr_genes.iterrows():
            start_pos = int(gene['bin_number'] * bin_size_kb * 1000)
            end_pos = start_pos + (bin_size_kb * 1000)
            
            busco_entries.append({
                'busco_id': gene['busco_id'],
                'status': 'Complete',
                'sequence': chromosome,
                'gene_start': start_pos,
                'gene_end': end_pos,
                'strand': '+',
                'score': gene['overlap_score'],
                'length': bin_size_kb * 1000,
                'frequency': gene['frequency_ratio'],
                'genome_count': f"{gene['genome_count']}/{gene['total_genomes']}"
            })
    
    busco_df = pd.DataFrame(busco_entries)
    busco_df = busco_df.sort_values(['sequence', 'gene_start'])
    
    # Save to file
    print(f"Saving to: {output_tsv}")
    busco_df.to_csv(output_tsv, sep='\t', index=False)
    
    # Summary
    print(f"\nPositional Profile Summary:")
    print(f"  Total genes: {len(busco_df)}")
    print(f"  Chromosomes: {len(busco_df['sequence'].unique())}")
    print(f"  Average frequency: {busco_df['frequency'].mean():.2f}")
    
    chr_counts = busco_df['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        avg_freq = busco_df[busco_df['sequence'] == chr_name]['frequency'].mean()
        print(f"    {chr_name}: {count} genes (avg freq: {avg_freq:.2f})")
    
    return busco_df

def extract_ordinal_profile(ordinal_csv, output_tsv, min_genes=500, freq_threshold=0.8, rank_spacing_kb=40):
    """
    Extract ordinal profile to BUSCO TSV format.
    
    Args:
        ordinal_csv: Path to ordinal profile CSV
        output_tsv: Path to output BUSCO TSV
        min_genes: Minimum genes per chromosome
        freq_threshold: Minimum frequency threshold
        rank_spacing_kb: Spacing between ranks in kb
        
    Returns:
        busco_df: BUSCO-formatted DataFrame
    """
    
    print("=" * 60)
    print("EXTRACTING ORDINAL PROFILE")
    print("=" * 60)
    
    # Read ordinal profile
    print(f"Reading: {ordinal_csv}")
    df = pd.read_csv(ordinal_csv)
    print(f"Loaded {len(df)} ordinal entries")
    
    # Check what columns we have
    print(f"Available columns: {list(df.columns)}")
    
    # Extract chromosome and rank information (handle different column names)
    if 'bin_or_rank_id' in df.columns:
        rank_col = 'bin_or_rank_id'
    elif 'rank_id' in df.columns:
        rank_col = 'rank_id'
    else:
        print("Error: Cannot find rank ID column")
        print(f"Available columns: {list(df.columns)}")
        return None
    
    df['chromosome'] = df[rank_col].str.extract(r'(.+?)_rank_')[0]
    df['rank_number'] = df[rank_col].str.extract(r'_rank_(\d+)')[0].astype(int)
    
    # Calculate frequency ratio if possible
    if 'genome_count' in df.columns and 'total_genomes' in df.columns:
        df['frequency_ratio'] = df['genome_count'] / df['total_genomes']
        print(f"Frequency threshold: {freq_threshold} ({freq_threshold*100}% of genomes)")
        high_freq = df[df['frequency_ratio'] >= freq_threshold].copy()
        print(f"High-frequency entries: {len(high_freq)}/{len(df)} ({len(high_freq)/len(df)*100:.1f}%)")
    else:
        print("No frequency data found - using all entries")
        high_freq = df.copy()
        high_freq['frequency_ratio'] = 1.0
    
    # Filter chromosomes by gene count
    genes_per_chr = high_freq.groupby('chromosome')['busco_id'].nunique().sort_values(ascending=False)
    major_chromosomes = genes_per_chr[genes_per_chr >= min_genes].index.tolist()
    
    print(f"Chromosome filtering (min {min_genes} genes):")
    print(f"  Total chromosomes: {len(genes_per_chr)}")
    print(f"  Major chromosomes: {len(major_chromosomes)}")
    
    for chr_name in major_chromosomes:
        print(f"    {chr_name}: {genes_per_chr[chr_name]} genes")
    
    # Filter to major chromosomes
    major_data = high_freq[high_freq['chromosome'].isin(major_chromosomes)].copy()
    print(f"Filtered to {len(major_data)} entries on major chromosomes")
    
    # For ordinal, typically one rank per gene per chromosome
    print("Processing ordinal assignments...")
    ordinal_assignments = []
    
    for chromosome in sorted(major_chromosomes):
        chr_data = major_data[major_data['chromosome'] == chromosome]
        
        for busco_id in chr_data['busco_id'].unique():
            gene_entries = chr_data[chr_data['busco_id'] == busco_id]
            
            # For ordinal, take the highest frequency entry if multiple
            if len(gene_entries) > 1:
                best_entry = gene_entries.loc[gene_entries['frequency_ratio'].idxmax()]
            else:
                best_entry = gene_entries.iloc[0]
            
            ordinal_assignments.append({
                'busco_id': busco_id,
                'chromosome': chromosome,
                'rank_number': best_entry['rank_number'],
                'gene_position': best_entry.get('gene_position', best_entry['rank_number'] * rank_spacing_kb * 1000),
                'frequency_ratio': best_entry['frequency_ratio'],
                'genome_count': best_entry.get('genome_count', 1),
                'total_genomes': best_entry.get('total_genomes', 1)
            })
    
    ordinal_df = pd.DataFrame(ordinal_assignments)
    print(f"Ordinal assignments: {len(ordinal_df)} genes")
    
    # Convert to BUSCO format
    print("Converting to BUSCO format...")
    busco_entries = []
    
    for chromosome in sorted(ordinal_df['chromosome'].unique()):
        chr_genes = ordinal_df[ordinal_df['chromosome'] == chromosome].copy()
        chr_genes = chr_genes.sort_values('rank_number')  # Sort by rank order
        
        print(f"  {chromosome}: {len(chr_genes)} genes")
        
        for _, gene in chr_genes.iterrows():
            start_pos = int(gene['rank_number'] * rank_spacing_kb * 1000)
            end_pos = start_pos + (rank_spacing_kb * 1000)
            
            busco_entries.append({
                'busco_id': gene['busco_id'],
                'status': 'Complete',
                'sequence': chromosome,
                'gene_start': start_pos,
                'gene_end': end_pos,
                'strand': '+',
                'score': gene['frequency_ratio'] * 100,  # Convert to percentage
                'length': rank_spacing_kb * 1000,
                'frequency': gene['frequency_ratio'],
                'genome_count': f"{gene['genome_count']}/{gene['total_genomes']}",
                'rank': gene['rank_number']
            })
    
    busco_df = pd.DataFrame(busco_entries)
    busco_df = busco_df.sort_values(['sequence', 'rank'])
    
    # Save to file
    print(f"Saving to: {output_tsv}")
    busco_df.to_csv(output_tsv, sep='\t', index=False)
    
    # Summary
    print(f"\nOrdinal Profile Summary:")
    print(f"  Total genes: {len(busco_df)}")
    print(f"  Chromosomes: {len(busco_df['sequence'].unique())}")
    print(f"  Average frequency: {busco_df['frequency'].mean():.2f}")
    
    chr_counts = busco_df['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        avg_freq = busco_df[busco_df['sequence'] == chr_name]['frequency'].mean()
        print(f"    {chr_name}: {count} genes (avg freq: {avg_freq:.2f})")
    
    return busco_df

def main():
    parser = argparse.ArgumentParser(description='Extract GIQ profiles to separate positional and ordinal TSV files')
    parser.add_argument('positional_csv', help='Positional profile CSV (6a_positional_profile.csv)')
    parser.add_argument('ordinal_csv', help='Ordinal profile CSV (6b_ordinal_profile.csv)')
    parser.add_argument('positional_tsv', help='Output positional BUSCO TSV file')
    parser.add_argument('ordinal_tsv', help='Output ordinal BUSCO TSV file')
    parser.add_argument('--min-genes', type=int, default=500, help='Minimum genes per chromosome (default: 500)')
    parser.add_argument('--freq-threshold', type=float, default=0.8, help='Minimum frequency threshold (default: 0.8)')
    parser.add_argument('--bin-size', type=int, default=40, help='Bin size in kb for positional (default: 40)')
    parser.add_argument('--rank-spacing', type=int, default=40, help='Rank spacing in kb for ordinal (default: 40)')
    
    args = parser.parse_args()
    
    try:
        print("EXTRACTING GIQ PROFILES TO SEPARATE TSV FILES")
        print("=" * 80)
        
        # Extract positional profile
        print("\n1. EXTRACTING POSITIONAL PROFILE")
        positional_df = extract_positional_profile(
            args.positional_csv,
            args.positional_tsv,
            args.min_genes,
            args.freq_threshold,
            args.bin_size
        )
        
        print("\n" + "="*80)
        
        # Extract ordinal profile  
        print("\n2. EXTRACTING ORDINAL PROFILE")
        ordinal_df = extract_ordinal_profile(
            args.ordinal_csv,
            args.ordinal_tsv,
            args.min_genes,
            args.freq_threshold,
            args.rank_spacing
        )
        
        # Final summary
        print("\n" + "="*80)
        print("EXTRACTION COMPLETE")
        print("="*80)
        print(f"Positional ancestral genome: {args.positional_tsv}")
        print(f"  {len(positional_df)} genes across {len(positional_df['sequence'].unique())} chromosomes")
        print(f"Ordinal ancestral genome: {args.ordinal_tsv}")
        print(f"  {len(ordinal_df)} genes across {len(ordinal_df['sequence'].unique())} chromosomes")
        
        # Compare overlap
        pos_genes = set(positional_df['busco_id'])
        ord_genes = set(ordinal_df['busco_id'])
        overlap = pos_genes & ord_genes
        
        print(f"\nGene set comparison:")
        print(f"  Positional only: {len(pos_genes - ord_genes)}")
        print(f"  Ordinal only: {len(ord_genes - pos_genes)}")
        print(f"  Both profiles: {len(overlap)}")
        print(f"  Total unique: {len(pos_genes | ord_genes)}")
        
        print("✅ Dual profile extraction completed successfully!")
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())