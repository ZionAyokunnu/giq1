#!/usr/bin/env python3
"""
Extract AGORA ancestral genome into standard BUSCO TSV format.
"""
"""
Script:
python3 agora_to_csv.py \
  /Users/za7/Documents/giq/agora_results/ancGenome.552100.00.list.bz2 \
  agora_ancestral_genome.tsv

"""

import pandas as pd
import argparse
import bz2
from pathlib import Path


def extract_agora_ancestral_genome_to_busco(agora_file: str, output_tsv_path: str):
    """
    Convert AGORA ancestral genome to BUSCO TSV format.
    
    Args:
        agora_file: Path to AGORA ancGenome.*.list.bz2 file
        output_tsv_path: Path to output BUSCO TSV file
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
    
    # Parse AGORA ancestral genome format
    ancestral_genes = []
    current_chromosome = "chr1"  # Default chromosome for ancestral genome
    position = 0
    gene_length = 1000  # Default gene length for ancestral genes
    
    for line_num, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        
        # AGORA ancestral genome format varies, let's handle common formats
        parts = line.split('\t')
        
        if len(parts) >= 5:
            # Format: ancestor_id, block_size, gene_ids, orientations, weights
            ancestor_id = parts[0]
            block_size = int(parts[1])
            gene_ids = parts[2].split()
            orientations = parts[3].split() if parts[3] else ['+'] * len(gene_ids)
            
            for i, gene_id in enumerate(gene_ids):
                strand = '+' if orientations[i] == '1' else '-'
                
                # Extract BUSCO ID from gene name
                if '|' in gene_id:
                    busco_id = gene_id.split('|')[-1]
                else:
                    busco_id = gene_id
                
                start_pos = position * gene_length
                end_pos = start_pos + gene_length
                
                ancestral_genes.append({
                    'busco_id': busco_id,
                    'status': 'Complete',
                    'sequence': current_chromosome,
                    'gene_start': start_pos,
                    'gene_end': end_pos,
                    'strand': strand,
                    'score': 100.0,  # Perfect score for ancestral reconstruction
                    'length': gene_length,
                    'ancestral_position': position,
                    'block_size': block_size,
                    'original_gene_id': gene_id
                })
                
                position += 1
        
        elif len(parts) == 1:
            # Simple format: just gene IDs, one per line
            gene_id = parts[0]
            
            # Extract BUSCO ID
            if '|' in gene_id:
                busco_id = gene_id.split('|')[-1]
            else:
                busco_id = gene_id
            
            start_pos = position * gene_length
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
                'ancestral_position': position,
                'block_size': 1,
                'original_gene_id': gene_id
            })
            
            position += 1
    
    if not ancestral_genes:
        print("Warning: No genes found. Trying alternative parsing...")
        
        # Alternative: try parsing as simple gene list
        for line_num, line in enumerate(lines):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Treat each line as a gene ID
            gene_id = line
            busco_id = gene_id.split('|')[-1] if '|' in gene_id else gene_id
            
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
                'block_size': 1,
                'original_gene_id': gene_id
            })
    
    # Convert to DataFrame
    busco_df = pd.DataFrame(ancestral_genes)
    
    if len(busco_df) == 0:
        raise ValueError("No genes could be parsed from the AGORA file")
    
    # Create standard BUSCO TSV format
    busco_tsv = busco_df[['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length']]
    
    # Save to TSV
    busco_tsv.to_csv(output_tsv_path, sep='\t', index=False, header=True)
    
    print(f"Extracted {len(busco_tsv)} ancestral genes to: {output_tsv_path}")
    
    # Print summary
    print("\nAGORA Ancestral Genome Summary:")
    print(f"  Total genes: {len(busco_tsv)}")
    print(f"  Chromosomes: {sorted(busco_tsv['sequence'].unique())}")
    print(f"  Position range: {busco_tsv['gene_start'].min():,} - {busco_tsv['gene_end'].max():,}")
    
    # Show first few genes
    print(f"\nFirst 5 genes:")
    for i, row in busco_tsv.head().iterrows():
        print(f"  {row['busco_id']}: {row['gene_start']:,} - {row['gene_end']:,}")
    
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



