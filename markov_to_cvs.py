#!/usr/bin/env python3
"""
Extract GIQ Markov profile highest probability positions into standard BUSCO TSV format.
"""


"""
script:
python3 markov_to_cvs.py \
  giq_profile_results/stages/5_markov_profile.csv \
  giq_ancestral_genome.tsv \
  --format csv

"""

import pandas as pd
import json
import argparse
from pathlib import Path


def extract_giq_profile_to_busco(profile_csv_path: str, output_tsv_path: str):
    """
    Convert GIQ Markov profile to BUSCO TSV format using highest probability positions.
    
    Args:
        profile_csv_path: Path to GIQ markov_profile.csv (from stages/5_markov_profile.csv)
        output_tsv_path: Path to output BUSCO TSV file
    """
    
    # Read the GIQ profile
    print(f"Reading GIQ profile from: {profile_csv_path}")
    profile_df = pd.read_csv(profile_csv_path)
    
    print(f"Loaded {len(profile_df)} gene-bin entries")
    
    # Find highest probability position for each gene
    highest_prob_positions = []
    
    for busco_id in profile_df['busco_id'].unique():
        gene_data = profile_df[profile_df['busco_id'] == busco_id]
        
        # Find the bin with highest average percentage for this gene
        best_entry = gene_data.loc[gene_data['average_percentage'].idxmax()]
        
        # Extract bin information
        bin_id = best_entry['bin_id']
        chromosome = bin_id.split('_bin_')[0] if '_bin_' in bin_id else 'chr1'
        bin_number = int(bin_id.split('_bin_')[1]) if '_bin_' in bin_id else 0
        
        # Convert bin to genomic coordinates (approximate)
        # Assuming 100kb bins as default
        bin_size_kb = 100
        start_pos = bin_number * bin_size_kb * 1000
        end_pos = start_pos + bin_size_kb * 1000
        
        highest_prob_positions.append({
            'busco_id': busco_id,
            'status': 'Complete',
            'sequence': chromosome,
            'gene_start': start_pos,
            'gene_end': end_pos,
            'strand': '+',
            'score': best_entry['average_percentage'],
            'length': bin_size_kb * 1000,
            'bin_number': bin_number,
            'genome_frequency': best_entry['genome_frequency'],
            'probability': best_entry['average_percentage'] / 100.0
        })
    
    # Convert to DataFrame and sort by position
    busco_df = pd.DataFrame(highest_prob_positions)
    busco_df = busco_df.sort_values(['sequence', 'gene_start'])
    
    # Create standard BUSCO TSV format
    busco_tsv = busco_df[['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length']]
    
    # Save to TSV
    busco_tsv.to_csv(output_tsv_path, sep='\t', index=False, header=True)
    
    print(f"Extracted {len(busco_tsv)} genes to: {output_tsv_path}")
    
    # Print summary
    print("\nGIQ Profile Summary:")
    print(f"  Total genes: {len(busco_tsv)}")
    print(f"  Chromosomes: {sorted(busco_tsv['sequence'].unique())}")
    print(f"  Position range: {busco_tsv['gene_start'].min():,} - {busco_tsv['gene_end'].max():,}")
    print(f"  Average probability: {busco_df['probability'].mean():.3f}")
    
    return busco_tsv




def extract_giq_alternative_format(profile_json_path: str, output_tsv_path: str):
    """
    Alternative: Extract from the main markov_profile.json file.
    
    Args:
        profile_json_path: Path to markov_profile.json
        output_tsv_path: Path to output BUSCO TSV file
    """
    
    print(f"Reading GIQ profile from JSON: {profile_json_path}")
    
    with open(profile_json_path, 'r') as f:
        profile_data = json.load(f)
    
    markov_profile = profile_data['markov_profile']
    
    # Extract highest probability positions
    gene_best_positions = {}
    
    for bin_id, genes_data in markov_profile.items():
        chromosome = bin_id.split('_bin_')[0] if '_bin_' in bin_id else 'chr1'
        bin_number = int(bin_id.split('_bin_')[1]) if '_bin_' in bin_id else 0
        
        for busco_id, gene_info in genes_data.items():
            avg_percentage = gene_info['average_percentage']
            
            if busco_id not in gene_best_positions or avg_percentage > gene_best_positions[busco_id]['probability']:
                
                # Convert bin to genomic coordinates
                bin_size_kb = profile_data.get('config', {}).get('position_bin_size_kb', 100)
                start_pos = bin_number * bin_size_kb * 1000
                end_pos = start_pos + bin_size_kb * 1000
                
                gene_best_positions[busco_id] = {
                    'busco_id': busco_id,
                    'status': 'Complete',
                    'sequence': chromosome,
                    'gene_start': start_pos,
                    'gene_end': end_pos,
                    'strand': '+',
                    'score': avg_percentage,
                    'length': bin_size_kb * 1000,
                    'bin_number': bin_number,
                    'genome_frequency': gene_info['genome_frequency'],
                    'probability': avg_percentage / 100.0
                }
    
    # Convert to DataFrame and sort
    busco_df = pd.DataFrame(list(gene_best_positions.values()))
    busco_df = busco_df.sort_values(['sequence', 'gene_start'])
    
    # Create standard BUSCO TSV format
    busco_tsv = busco_df[['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length']]
    
    # Save to TSV
    busco_tsv.to_csv(output_tsv_path, sep='\t', index=False, header=True)
    
    print(f"Extracted {len(busco_tsv)} genes to: {output_tsv_path}")
    
    # Print summary
    print("\nGIQ Profile Summary:")
    print(f"  Total genes: {len(busco_tsv)}")
    print(f"  Chromosomes: {sorted(busco_tsv['sequence'].unique())}")
    print(f"  Position range: {busco_tsv['gene_start'].min():,} - {busco_tsv['gene_end'].max():,}")
    print(f"  Average probability: {busco_df['probability'].mean():.3f}")
    
    return busco_tsv




def main():
    parser = argparse.ArgumentParser(description='Extract GIQ Markov profile to BUSCO TSV')
    parser.add_argument('input_file', help='Input file (CSV or JSON)')
    parser.add_argument('output_tsv', help='Output BUSCO TSV file')
    parser.add_argument('--format', choices=['csv', 'json'], default='csv',
                       help='Input format: csv (stages/5_markov_profile.csv) or json (markov_profile.json)')
    
    args = parser.parse_args()
    
    try:
        if args.format == 'csv':
            extract_giq_profile_to_busco(args.input_file, args.output_tsv)
        else:
            extract_giq_alternative_format(args.input_file, args.output_tsv)
        
        print("✅ GIQ profile extraction completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0




if __name__ == "__main__":
    import sys
    sys.exit(main())



