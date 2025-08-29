#!/usr/bin/env python3
"""
Extract AGORA ancestral genome into standard BUSCO TSV format.
Handles both standard and merged chromosome-aware formats.



python3 /Users/zionayokunnu/Documents/Giq/agora_to_csv.py \
    /Users/zionayokunnu/Documents/Giq/homology_results/merged_ancGenome.504100.00.list.bz2 \
    /Users/zionayokunnu/Documents/Giq/compare/root_chr-aware_agora_ancestral_genome.tsv \
    --merged \
    --gene-length 288000


python3 /Users/zionayokunnu/Documents/Giq/agora_to_csv.py \
    /Users/zionayokunnu/Documents/Giq/agora_results/ancGenome.504100.00.list.bz2 \
    /Users/zionayokunnu/Documents/Giq/compare/root_agora_ancestral_genome.tsv \
    --gene-length 288000
    
    
    
    
python3 /Users/zionayokunnu/Documents/Giq/agora_to_csv.py \
    /Users/zionayokunnu/Documents/Giq/testingAgora/merged_ancGenome.504100.00.list.bz2 \
    /Users/zionayokunnu/Documents/Giq/testingAgora/root_chr-aware_agora_ancestral_genome.tsv \
    --merged \
    --gene-length 288000
    
python3 /Users/zionayokunnu/Documents/Giq/agora_to_csv.py \
    /Users/zionayokunnu/Documents/Giq/testingAgora/ancGenome.504100.00.list.bz2 \
    /Users/zionayokunnu/Documents/Giq/testingAgora/root_agora_ancestral_genome.tsv \
    --gene-length 288000
    
"""

import pandas as pd
import argparse
import bz2
from pathlib import Path


def extract_agora_ancestral_genome_to_busco(agora_file: str, output_tsv_path: str, is_merged: bool = False, gene_length: int = 288000):
    """
    Convert AGORA ancestral genome to BUSCO TSV format.
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
    chromosome_boundaries = []
    current_chromosome_id = 1
    current_chromosome = f"chr{current_chromosome_id}"
    genes_in_current_chr = 0
    
    previous_chromosome_group = None
    running_position_offset = 0
    previous_agora_chr_id = None
    max_position_in_current_agora_chr = 0
    
    processed_lines = 0
    for line_num, line in enumerate(lines, 1):
        
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        parts = line.split('\t')
        # parts = line.split() 
        
        if is_merged:
            # Merged format: block start end strand ancestral_id gene_entries... chromosome_group
            if len(parts) >= 6:
                agora_chr_id = int(parts[0])
                agora_start = int(parts[1])
                agora_end = int(parts[2]) 
                chromosome_group = parts[5]  # Last field
                field5_parts = parts[4].split()  # Split the gene field by spaces
                ancestral_id = field5_parts[0]  # First part is ancestral ID
                gene_entries = field5_parts[1:]  # Rest are gene entries Exclude last column
                        
                # if line_num <= 10 or chromosome_group != previous_chromosome_group:
                #     print(f"DEBUG Line {line_num}: agora_chr={agora_chr_id}, chr_group={chromosome_group}, agora_pos={agora_start}-{agora_end}")
                
            
            else:
                continue
            
            
            max_position_in_current_agora_chr = max(max_position_in_current_agora_chr, agora_end)
            
                # Check for biological chromosome change
            if chromosome_group != previous_chromosome_group:
                # New biological chromosome - reset everything
                current_chromosome = chromosome_group
                running_position_offset = 0
                previous_agora_chr_id = agora_chr_id
                max_position_in_current_agora_chr = agora_end
                previous_chromosome_group = chromosome_group
            
            # Check for AGORA chromosome change within same biological chromosome  
            elif agora_chr_id != previous_agora_chr_id:
                # Same biological chromosome, new AGORA chromosome
                # Add offset to continue coordinates from where previous AGORA chr ended
                running_position_offset += max_position_in_current_agora_chr
                previous_agora_chr_id = agora_chr_id
                max_position_in_current_agora_chr = agora_end
            
            # # Calculate continuous coordinates preserving AGORA order
            # start_pos = (running_position_offset + agora_start) * gene_length
            # end_pos = (running_position_offset + agora_end) * gene_length
            
            
            
            try:
                # Calculate continuous coordinates preserving AGORA order
                start_pos = (running_position_offset + agora_start) * gene_length
                end_pos = (running_position_offset + agora_end) * gene_length
            except Exception as e:
                print(f"ERROR at line {line_num}: {e}")
                print(f"  agora_start={agora_start}, agora_end={agora_end}")
                print(f"  running_offset={running_position_offset}")
                continue


        else:
            
            parts = line.split()
            # Standard format: block start end strand ancestral_id gene_entries...
            if len(parts) >= 6:
                agora_chromosome_id = int(parts[0])  # AGORA's chromosome assignment
                ancestral_id = parts[4]  # 504100.00.1
                gene_entries = parts[5:]  # All the Species|BUSCO_ID entries
            else:
                continue
                
            # Check if we've moved to a new chromosome
            if agora_chromosome_id != current_chromosome_id:
                if genes_in_current_chr > 0:
                    chromosome_boundaries.append({
                        'chromosome': current_chromosome,
                        'gene_count': genes_in_current_chr,
                        'end_line': line_num - 1
                    })
                
                current_chromosome_id = agora_chromosome_id
                current_chromosome = f"chr{current_chromosome_id}"
                genes_in_current_chr = 0
                
         
            # Original artificial spacing logic
            position_in_chr = genes_in_current_chr
            start_pos = position_in_chr * gene_length
            end_pos = start_pos + gene_length

            
        processed_lines += 1
        
        # Extract BUSCO IDs from the gene entries
        busco_ids = []
        for entry in gene_entries:
            if '|' in entry:
                # Extract BUSCO ID from "Species|BUSCO_ID" format
                busco_id = entry.split('|')[-1]  # Take the part after the last |
                # if busco_id and 'at' in busco_id.lower():  # Validate it's a real BUSCO ID
                #     busco_ids.append(busco_id)
                if busco_id:  # Remove the 'at' validation
                    busco_ids.append(busco_id.strip())  # Add strip() for whitespace

            # Debug: Print what we're processing
            # print(f"DEBUG: Processing {len(gene_entries)} gene entries")
            for i, entry in enumerate(gene_entries[:3]):
                # print(f"DEBUG: Entry {i}: '{entry}'")
                if '|' in entry:
                    busco_id = entry.split('|')[-1]
                    # print(f"DEBUG: Extracted BUSCO ID: '{busco_id}', has 'at': {'at' in busco_id.lower()}")
            
            # This debug code runs for EVERY entry, not just |entries
            for i, entry in enumerate(gene_entries[:3]):
                if '|' in entry:
                    busco_id = entry.split('|')[-1]
        
        if not busco_ids and ancestral_id:
            # Fallback: create synthetic BUSCO ID when no gene mappings exist
            synthetic_busco_id = f"agora_{line_num:06d}_{ancestral_id.split('.')[-1]}"
            busco_ids = [synthetic_busco_id]
            
        # Use the first valid BUSCO ID we find (they should all be the same)
        if busco_ids:
            busco_id = busco_ids[0]  # Take first one

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
                'original_entries': '|'.join(gene_entries)
            })
            
            if not is_merged:
                genes_in_current_chr += 1

    # Record the final chromosome for standard format
    if not is_merged and genes_in_current_chr > 0:
        chromosome_boundaries.append({
            'chromosome': current_chromosome,
            'gene_count': genes_in_current_chr,
            'end_line': len(lines)
        })
        
    print(f"Processed {processed_lines} lines total")
    
    if not ancestral_genes:
        raise ValueError("No valid BUSCO genes could be parsed from the AGORA file")

    # Convert to DataFrame
    busco_df = pd.DataFrame(ancestral_genes)

    # Create standard BUSCO TSV format
    busco_tsv = busco_df[['busco_id', 'status', 'sequence', 'gene_start', 'gene_end', 'strand', 'score', 'length']]

    # Save to TSV
    busco_tsv.to_csv(output_tsv_path, sep='\t', index=False, header=True)

    print(f"Extracted {len(busco_tsv)} ancestral genes to: {output_tsv_path}")

    # Print summary
    print(f"\nAGORA Ancestral Genome Summary:")
    print(f"  Total genes: {len(busco_tsv)}")
    # print(f"  Chromosomes: {sorted(busco_tsv['sequence'].unique())}")
    print(f"  Position range: {busco_tsv['gene_start'].min():,} - {busco_tsv['gene_end'].max():,}")

    # Show chromosome distribution
    # print(f"\nChromosome distribution:")
    # chr_counts = busco_tsv['sequence'].value_counts().sort_index()
    # for chr_name, count in chr_counts.items():
    #     print(f"  {chr_name}: {count} genes")
        
    return busco_tsv


def main():
    parser = argparse.ArgumentParser(description='Extract AGORA ancestral genome to BUSCO TSV')
    parser.add_argument('agora_file', help='AGORA ancestral genome file (.bz2 or plain text)')
    parser.add_argument('output_tsv', help='Output BUSCO TSV file')
    parser.add_argument('--merged', action='store_true', help='Input is merged chromosome-aware format')
    parser.add_argument('--gene-length', type=int, default=288000, help='Gene length for spacing (default: 288000)')

    args = parser.parse_args()

    try:
        extract_agora_ancestral_genome_to_busco(args.agora_file, args.output_tsv, args.merged, args.gene_length)
        print("AGORA ancestral genome extraction completed successfully!")

    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())