#!/usr/bin/env python3
"""
Extract GIQ Markov profile highest probability positions into standard BUSCO TSV format.
Now includes chromosome consolidation by name across species.
"""


"""
script:  

python3 markov_to_csv.py \
    comparison_output/stages/5_markov_profile.csv \
    compare/root_giq_ancestral_genome2.tsv \
    --consolidation-tsv giq_consolidation_details.tsv

"""

import pandas as pd
import json
import argparse
from pathlib import Path


def consolidate_giq_by_chromosome_name(profile_df):
    """
    Consolidate chromosomes by their names across species.
    For each gene, select the position with highest probability within each chromosome.
    
    Args:
        profile_df: DataFrame with columns [bin_id, busco_id, average_percentage, ...]
        
    Returns:
        consolidated_df: DataFrame with best positions per gene per chromosome
    """
    
    print("Consolidating GIQ chromosomes by name...")
    
    # Extract chromosome names from bin_ids
    profile_df['chromosome'] = profile_df['bin_id'].str.extract(r'(.+?)_bin_')[0]
    profile_df['bin_number'] = profile_df['bin_id'].str.extract(r'_bin_(\d+)')[0].astype(int)
    
    print(f"Found chromosomes: {sorted(profile_df['chromosome'].unique())}")
    
    # Group by chromosome and gene, then find best position
    consolidated_genes = []
    
    for chromosome in sorted(profile_df['chromosome'].unique()):
        chr_data = profile_df[profile_df['chromosome'] == chromosome]
        print(f"  Processing {chromosome}: {len(chr_data)} bin entries")
        
        # For each gene in this chromosome, find the bin with highest probability
        for busco_id in chr_data['busco_id'].unique():
            gene_bins = chr_data[chr_data['busco_id'] == busco_id]
            
            # Find bin with highest average_percentage
            best_bin = gene_bins.loc[gene_bins['average_percentage'].idxmax()]
            
            consolidated_genes.append({
                'busco_id': busco_id,
                'chromosome': chromosome,
                'bin_number': best_bin['bin_number'],
                'average_percentage': best_bin['average_percentage'],
                'genome_frequency': best_bin['genome_frequency'],
                'genome_count': best_bin['genome_count'],
                'total_genomes': best_bin['total_genomes']
            })
    
    consolidated_df = pd.DataFrame(consolidated_genes)
    
    print(f"Consolidated to {len(consolidated_df)} gene positions across {len(consolidated_df['chromosome'].unique())} chromosomes")
    
    # Print chromosome summary
    chr_counts = consolidated_df['chromosome'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        print(f"  {chr_name}: {count} genes")
    
    return consolidated_df


def extract_giq_profile_to_busco(profile_csv_path: str, output_tsv_path: str):
    """
    Convert GIQ Markov profile to BUSCO TSV format using consolidated chromosomes.
    
    Args:
        profile_csv_path: Path to GIQ markov_profile.csv (from stages/5_markov_profile.csv)
        output_tsv_path: Path to output BUSCO TSV file
    """
    
    # Read the GIQ profile
    print(f"Reading GIQ profile from: {profile_csv_path}")
    profile_df = pd.read_csv(profile_csv_path)
    
    print(f"Loaded {len(profile_df)} gene-bin entries")
    
    # Consolidate chromosomes by name
    consolidated_df = consolidate_giq_by_chromosome_name(profile_df)
    
    # Convert to BUSCO TSV format
    busco_entries = []
    
    # Process each chromosome separately to maintain gene order
    for chromosome in sorted(consolidated_df['chromosome'].unique()):
        chr_genes = consolidated_df[consolidated_df['chromosome'] == chromosome].copy()
        
        # Sort genes by bin number (position within chromosome)
        chr_genes = chr_genes.sort_values('bin_number')
        
        print(f"Converting {chromosome}: {len(chr_genes)} genes")
        
        # Convert bin positions to genomic coordinates
        for _, gene in chr_genes.iterrows():
            # Assuming 100kb bins as default
            bin_size_kb = 100
            start_pos = gene['bin_number'] * bin_size_kb * 1000
            end_pos = start_pos + bin_size_kb * 1000
            
            busco_entries.append({
                'busco_id': gene['busco_id'],
                'status': 'Complete',
                'sequence': chromosome,  # Use consolidated chromosome name
                'gene_start': start_pos,
                'gene_end': end_pos,
                'strand': '+',
                'score': gene['average_percentage'],
                'length': bin_size_kb * 1000
            })
    
    # Create DataFrame and save
    busco_df = pd.DataFrame(busco_entries)
    busco_df = busco_df.sort_values(['sequence', 'gene_start'])
    
    # Save to TSV (standard BUSCO format)
    busco_df.to_csv(output_tsv_path, sep='\t', index=False, header=True)
    
    print(f"Extracted {len(busco_df)} genes to: {output_tsv_path}")
    
    # Print final summary
    print("\nGIQ Consolidated Profile Summary:")
    print(f"  Total genes: {len(busco_df)}")
    print(f"  Chromosomes: {sorted(busco_df['sequence'].unique())}")
    print(f"  Position range: {busco_df['gene_start'].min():,} - {busco_df['gene_end'].max():,}")
    
    # Show chromosome distribution
    chr_counts = busco_df['sequence'].value_counts().sort_index()
    print("  Gene distribution:")
    for chr_name, count in chr_counts.items():
        print(f"    {chr_name}: {count} genes")
    
    return busco_df


def extract_giq_alternative_format(profile_json_path: str, output_tsv_path: str):
    """
    Alternative: Extract from the main markov_profile.json file with consolidation.
    
    Args:
        profile_json_path: Path to markov_profile.json
        output_tsv_path: Path to output BUSCO TSV file
    """
    
    print(f"Reading GIQ profile from JSON: {profile_json_path}")
    
    with open(profile_json_path, 'r') as f:
        profile_data = json.load(f)
    
    markov_profile = profile_data['markov_profile']
    
    # Convert JSON to DataFrame format for consolidation
    profile_entries = []
    
    for bin_id, genes_data in markov_profile.items():
        for busco_id, gene_info in genes_data.items():
            profile_entries.append({
                'bin_id': bin_id,
                'busco_id': busco_id,
                'average_percentage': gene_info['average_percentage'],
                'genome_frequency': gene_info['genome_frequency'],
                'genome_count': gene_info.get('genome_count', 1),
                'total_genomes': gene_info.get('total_genomes', 1)
            })
    
    profile_df = pd.DataFrame(profile_entries)
    
    # Use the same consolidation logic
    consolidated_df = consolidate_giq_by_chromosome_name(profile_df)
    
    # Convert to BUSCO format
    busco_entries = []
    
    for chromosome in sorted(consolidated_df['chromosome'].unique()):
        chr_genes = consolidated_df[consolidated_df['chromosome'] == chromosome].copy()
        chr_genes = chr_genes.sort_values('bin_number')
        
        print(f"Converting {chromosome}: {len(chr_genes)} genes")
        
        for _, gene in chr_genes.iterrows():
            # Get bin size from config if available
            bin_size_kb = profile_data.get('config', {}).get('position_bin_size_kb', 100)
            start_pos = gene['bin_number'] * bin_size_kb * 1000
            end_pos = start_pos + bin_size_kb * 1000
            
            busco_entries.append({
                'busco_id': gene['busco_id'],
                'status': 'Complete',
                'sequence': chromosome,
                'gene_start': start_pos,
                'gene_end': end_pos,
                'strand': '+',
                'score': gene['average_percentage'],
                'length': bin_size_kb * 1000
            })
    
    # Create DataFrame and save
    busco_df = pd.DataFrame(busco_entries)
    busco_df = busco_df.sort_values(['sequence', 'gene_start'])
    
    busco_df.to_csv(output_tsv_path, sep='\t', index=False, header=True)
    
    print(f"Extracted {len(busco_df)} genes to: {output_tsv_path}")
    
    # Print summary
    print("\nGIQ Consolidated Profile Summary:")
    print(f"  Total genes: {len(busco_df)}")
    print(f"  Chromosomes: {sorted(busco_df['sequence'].unique())}")
    
    chr_counts = busco_df['sequence'].value_counts().sort_index()
    print("  Gene distribution:")
    for chr_name, count in chr_counts.items():
        print(f"    {chr_name}: {count} genes")
    
    return busco_df


def main():
    parser = argparse.ArgumentParser(description='Extract GIQ Markov profile to BUSCO TSV with chromosome consolidation')
    parser.add_argument('input_file', help='Input file (CSV or JSON)')
    parser.add_argument('output_tsv', help='Output BUSCO TSV file')
    parser.add_argument('--format', choices=['csv', 'json'], default='csv',
                       help='Input format: csv (stages/5_markov_profile.csv) or json (markov_profile.json)')
    parser.add_argument('--consolidation-tsv', help='Optional: Output consolidation details to TSV for debugging')
    
    args = parser.parse_args()
    
    try:
        if args.format == 'csv':
            busco_df = extract_giq_profile_to_busco(args.input_file, args.output_tsv)
        else:
            busco_df = extract_giq_alternative_format(args.input_file, args.output_tsv)
        
        # Save consolidation details if requested
        if args.consolidation_tsv and 'consolidated_df' in locals():
            print(f"Saving consolidation details to: {args.consolidation_tsv}")
            consolidated_df.to_csv(args.consolidation_tsv, sep='\t', index=False)
            print(f"  Consolidation details: {len(consolidated_df)} gene-chromosome pairs")
        
        print("✅ GIQ profile extraction with consolidation completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())