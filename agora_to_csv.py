#!/usr/bin/env python3
"""
Extract AGORA ancestral genome into standard BUSCO TSV format.
"""
"""
Script:
python3 agora_to_csv.py \
  /Users/za7/Documents/giq/agora_results/ancGenome.553100.00.list.bz2 \
  agora_ancestral_genome2.tsv

"""

import pandas as pd
import argparse
import bz2
from pathlib import Path


def extract_agora_ancestral_genome_to_busco(agora_file: str, output_tsv_path: str):
    """
    Convert AGORA ancestral genome to BUSCO TSV format.
    Fixed to properly parse AGORA format.
    """
    
    print(f"Reading AGORA ancestral genome from: {agora_file}")
    
    # Read the bzipped file
    if agora_file.endswith('.bz2'):
        with bz2.open(agora_file, 'rt') as f:
            lines = f.readlines()
    else:
        with open(agora_file, 'r') as f:
            lines = f.readlines()
    
    print(f"Read {len(lines)} lines from AGORA file")
    
    ancestral_genes = []
    current_chromosome = "chr1"  
    gene_length = 1000  
    
    for line_num, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        parts = line.split('\t')
        
        # Handle both AGORA formats
        if len(parts) >= 6:
            # First format: 0  0  1  1  552100.00.1  Species|BUSCO_ID  Species|BUSCO_ID  ...
            ancestral_id = parts[4]  # 552100.00.1
            gene_entries = parts[5:]  # All the Species|BUSCO_ID entries
        elif len(parts) >= 2:
            # Second format: 552100.00.1  Species|BUSCO_ID  Species|BUSCO_ID  ...
            ancestral_id = parts[0]  # 552100.00.1
            gene_entries = parts[1:]  # All the Species|BUSCO_ID entries
        else:
            print(f"Warning: Skipping malformed line {line_num}: {line}")
            continue
        
        # Extract BUSCO IDs from the gene entries
        busco_ids = []
        for entry in gene_entries:
            if '|' in entry:
                # Extract BUSCO ID from "Species|BUSCO_ID" format
                busco_id = entry.split('|')[-1]  # Take the part after the last |
                if busco_id and 'at7147' in busco_id:  # Validate it's a real BUSCO ID
                    busco_ids.append(busco_id)
        
        # Use the first valid BUSCO ID we find (they should all be the same)
        if busco_ids:
            busco_id = busco_ids[0]  # Take first one
            
            start_pos = line_num * gene_length
            end_pos = start_pos + gene_length
            
            ancestral_genes.append({
                'busco_id': busco_id,
                'status': 'Complete',
                'sequence': current_chromosome,
                'gene_start': start_pos,
                'gene_end': end_pos,
                'strand': '+',
                'score': 100.0,
                'length': gene_length,
                'ancestral_position': line_num,
                'ancestral_id': ancestral_id,
                'original_entries': gene_entries
            })
        else:
            print(f"Warning: No valid BUSCO ID found in line {line_num}: {line[:100]}...")
    
    if not ancestral_genes:
        raise ValueError("No valid BUSCO genes could be parsed from the AGORA file")
    
    # Convert to DataFrame
    busco_df = pd.DataFrame(ancestral_genes)
    
    # Create standard BUSCO TSV format
    busco_tsv = busco_df[['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length']]
    
    # Save to TSV
    busco_tsv.to_csv(output_tsv_path, sep='\t', index=False, header=True)
    
    print(f"Extracted {len(busco_tsv)} ancestral genes to: {output_tsv_path}")
    
    # Print summary with actual BUSCO IDs
    print("\nAGORA Ancestral Genome Summary:")
    print(f"  Total genes: {len(busco_tsv)}")
    print(f"  Chromosomes: {sorted(busco_tsv['sequence'].unique())}")
    print(f"  Position range: {busco_tsv['gene_start'].min():,} - {busco_tsv['gene_end'].max():,}")
    
    # Show first few genes with actual BUSCO IDs
    print(f"\nFirst 10 BUSCO IDs found:")
    for i, busco_id in enumerate(busco_tsv['busco_id'].head(10)):
        print(f"  {i}: {busco_id}")
    
    return busco_tsv




def main():
    parser = argparse.ArgumentParser(description='Extract AGORA ancestral genome to BUSCO TSV')
    parser.add_argument('agora_file', help='AGORA ancestral genome file (.bz2 or plain text)')
    parser.add_argument('output_tsv', help='Output BUSCO TSV file')
    
    args = parser.parse_args()
    
    try:
        extract_agora_ancestral_genome_to_busco(args.agora_file, args.output_tsv)
        print("✅ AGORA ancestral genome extraction completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0




if __name__ == "__main__":
    import sys
    sys.exit(main())



