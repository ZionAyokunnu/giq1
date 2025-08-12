#!/usr/bin/env python3
"""
Extract AGORA ancestral genome into standard BUSCO TSV format.
Maps ancestral genes back to original chromosome names using input BUSCO files.
"""
"""
Script:
python3 agora_to_csv.py \
  /Users/za7/Documents/giq/agora_results/ancGenome.553100.00.list.bz2 \
  /Users/za7/Documents/Bibionidae/busco-tables/bibio_marci.tsv \
  /Users/za7/Documents/Bibionidae/busco-tables/plecia_longiforceps.tsv \
  /Users/za7/Documents/Bibionidae/busco-tables/dilophus_febrilis.tsv \
  agora_ancestral_genome.tsv

"""

import pandas as pd
import argparse
import bz2
from pathlib import Path
from collections import Counter


def load_busco_file(busco_file_path: str, species_name: str):
    """Load a BUSCO TSV file and return a dictionary mapping BUSCO_ID -> chromosome"""
    print(f"Loading {species_name} BUSCO file: {busco_file_path}")
    
    try:
        # Read the file to find the header line
        with open(busco_file_path, 'r') as f:
            lines = f.readlines()
        
        # Find the header line (starts with # and contains column names)
        header_line = None
        data_start_idx = 0
        
        for i, line in enumerate(lines):
            if line.startswith('#') and ('Busco id' in line or 'busco_id' in line):
                header_line = line.strip('# \n').split('\t')
                data_start_idx = i + 1
                break
        
        if header_line is None:
            # Try to read normally if no commented header found
            df = pd.read_csv(busco_file_path, sep='\t', comment='#')
        else:
            # Use the found header and read from data start
            data_lines = lines[data_start_idx:]
            # Create a temporary file-like object
            from io import StringIO
            data_content = ''.join(data_lines)
            df = pd.read_csv(StringIO(data_content), sep='\t', names=header_line)
        
        print(f"  Loaded {len(df)} genes from {species_name}")
        
        # Clean up column names (remove extra spaces)
        df.columns = df.columns.str.strip()
        
        # Handle different possible column names
        busco_id_col = None
        status_col = None
        sequence_col = None
        
        for col in df.columns:
            col_lower = col.lower()
            if 'busco' in col_lower and 'id' in col_lower:
                busco_id_col = col
            elif col_lower in ['status', 'state']:
                status_col = col
            elif col_lower in ['sequence', 'scaffold', 'chromosome', 'contig']:
                sequence_col = col
        
        if not all([busco_id_col, status_col, sequence_col]):
            print(f"  Available columns: {list(df.columns)}")
            raise ValueError(f"Could not find required columns. Found: busco_id={busco_id_col}, status={status_col}, sequence={sequence_col}")
        
        print(f"  Using columns: {busco_id_col} | {status_col} | {sequence_col}")
        
        # Create mapping: busco_id -> chromosome
        busco_to_chr = {}
        complete_genes = df[df[status_col] == 'Complete']
        
        for _, row in complete_genes.iterrows():
            busco_id = row[busco_id_col]
            chromosome = row[sequence_col]  # chromosome name
            busco_to_chr[busco_id] = chromosome
        
        print(f"  Found {len(busco_to_chr)} complete BUSCO genes in {species_name}")
        return busco_to_chr
        
    except Exception as e:
        print(f"Error loading {species_name} BUSCO file: {e}")
        raise


def resolve_chromosome_assignment(busco_id: str, species_mappings: dict):
    """
    Resolve which chromosome a BUSCO gene should be assigned to based on input species.
    
    Args:
        busco_id: The BUSCO gene ID
        species_mappings: Dict of {species_name: {busco_id: chromosome}}
    
    Returns:
        tuple: (assigned_chromosome, confidence_metrics)
    """
    
    # Collect chromosomes this BUSCO appears on across species
    chromosomes = []
    species_with_gene = []
    
    for species_name, busco_mapping in species_mappings.items():
        if busco_id in busco_mapping:
            chromosome = busco_mapping[busco_id]
            chromosomes.append(chromosome)
            species_with_gene.append(species_name)
    
    if not chromosomes:
        return None, {'status': 'not_found', 'species_count': 0, 'chromosomes': []}
    
    # Count chromosome occurrences
    chr_counts = Counter(chromosomes)
    most_common_chr, max_count = chr_counts.most_common(1)[0]
    
    # Determine confidence/conflict status
    total_species = len(chromosomes)
    unique_chromosomes = len(chr_counts)
    
    if unique_chromosomes == 1:
        # All species agree
        status = 'unanimous'
    elif max_count > 1:
        # Majority consensus (2+ species on same chromosome)
        status = 'majority_consensus'
    else:
        # All different chromosomes, arbitrary choice
        status = 'conflict_arbitrary'
    
    metrics = {
        'status': status,
        'species_count': total_species,
        'chromosome_votes': dict(chr_counts),
        'chromosomes': chromosomes,
        'species_with_gene': species_with_gene,
        'confidence': max_count / total_species
    }
    
    return most_common_chr, metrics


def extract_agora_ancestral_genome_to_busco(agora_file: str, busco_files: list, output_tsv_path: str):
    """
    Convert AGORA ancestral genome to BUSCO TSV format.
    Maps genes back to original chromosome names using input BUSCO files.
    """
    
    print(f"Reading AGORA ancestral genome from: {agora_file}")
    
    # Load all input BUSCO files
    species_mappings = {}
    for busco_file in busco_files:
        # Extract species name from file path
        species_name = Path(busco_file).stem.replace('_busco', '').replace('.tsv', '')
        species_mappings[species_name] = load_busco_file(busco_file, species_name)
    
    # Read the AGORA bzipped file
    if agora_file.endswith('.bz2'):
        with bz2.open(agora_file, 'rt') as f:
            lines = f.readlines()
    else:
        with open(agora_file, 'r') as f:
            lines = f.readlines()
    
    print(f"Read {len(lines)} lines from AGORA file")
    
    ancestral_genes = []
    gene_length = 1000000  # 1Mb spacing between genes
    
    # Track chromosome positions for realistic coordinates
    chromosome_positions = {}
    
    # Track metrics for chromosome assignment
    assignment_metrics = {
        'unanimous': 0,
        'majority_consensus': 0,
        'conflict_arbitrary': 0,
        'not_found': 0
    }
    conflict_details = []
    
    for line_num, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        # Split on tabs first (should give us exactly 5 fields)
        parts = line.split('\t')
        
        if len(parts) != 5:
            print(f"Warning: Expected 5 tab-separated fields, got {len(parts)} in line {line_num}: {line[:100]}...")
            continue
        
        # Parse the first 4 tab-separated fields
        try:
            agora_chromosome_id = int(parts[0])  # AGORA's reconstructed chromosome ID
            position_in_chr = int(parts[1])  # Position within AGORA chromosome
            gene_order = int(parts[2])  # Gene order
            orientation = int(parts[3])  # Orientation (1 = +, -1 = -)
        except ValueError as e:
            print(f"Warning: Could not parse numeric fields in line {line_num}: {e}")
            continue
        
        # The 5th field contains space-separated entries: ancestral_id + gene entries
        fifth_field_parts = parts[4].split()
        
        if len(fifth_field_parts) < 1:
            print(f"Warning: No entries in 5th field of line {line_num}")
            continue
        
        ancestral_id = fifth_field_parts[0]  # First space-separated entry
        gene_entries = fifth_field_parts[1:]  # Remaining space-separated entries
        
        # Determine strand from orientation
        strand = '+' if orientation == 1 else '-'
        
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
            
            # Resolve chromosome assignment using input BUSCO files
            assigned_chromosome, metrics = resolve_chromosome_assignment(busco_id, species_mappings)
            
            # Track assignment metrics
            assignment_metrics[metrics['status']] += 1
            
            # Record conflicts for reporting
            if metrics['status'] in ['conflict_arbitrary', 'majority_consensus']:
                conflict_details.append({
                    'busco_id': busco_id,
                    'status': metrics['status'],
                    'chromosomes': metrics['chromosomes'],
                    'species': metrics['species_with_gene'],
                    'assigned_to': assigned_chromosome,
                    'confidence': metrics['confidence']
                })
            
            if assigned_chromosome is None:
                print(f"Warning: Could not resolve chromosome for {busco_id}, using AGORA assignment chr{agora_chromosome_id}")
                assigned_chromosome = f"chr{agora_chromosome_id}"
            
            # Track position within each chromosome
            if assigned_chromosome not in chromosome_positions:
                chromosome_positions[assigned_chromosome] = 0
            
            # Calculate realistic genomic coordinates
            start_pos = chromosome_positions[assigned_chromosome]
            end_pos = start_pos + gene_length
            chromosome_positions[assigned_chromosome] += gene_length
            
            ancestral_genes.append({
                'busco_id': busco_id,
                'status': 'Complete',
                'sequence': assigned_chromosome,  # Now uses resolved chromosome name
                'gene_start': start_pos,
                'gene_end': end_pos,
                'strand': strand,
                'score': 100.0,
                'length': gene_length,
                'ancestral_position': line_num,
                'ancestral_id': ancestral_id,
                'agora_chromosome_id': agora_chromosome_id,
                'position_in_chr': position_in_chr,
                'gene_order': gene_order,
                'assignment_status': metrics['status'],
                'assignment_confidence': metrics.get('confidence', 0.0),
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
    
    # Print detailed summary with chromosome assignment metrics
    print("\n" + "="*70)
    print("AGORA ANCESTRAL GENOME SUMMARY")
    print("="*70)
    print(f"Total genes: {len(busco_tsv)}")
    print(f"Chromosomes: {sorted(busco_tsv['sequence'].unique())}")
    print(f"Position range: {busco_tsv['gene_start'].min():,} - {busco_tsv['gene_end'].max():,}")
    
    # Chromosome assignment metrics
    print("\nCHROMOSOME ASSIGNMENT METRICS:")
    print("-" * 40)
    total_genes = len(ancestral_genes)
    for status, count in assignment_metrics.items():
        percentage = (count / total_genes) * 100
        print(f"  {status.replace('_', ' ').title()}: {count} ({percentage:.1f}%)")
    
    # Show chromosome distribution
    print(f"\nCHROMOSOME DISTRIBUTION:")
    print("-" * 30)
    chr_counts = busco_tsv['sequence'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        print(f"  {chr_name}: {count} genes")
    
    # Show conflicts if any
    if conflict_details:
        print(f"\nCONFLICT DETAILS (first 10):")
        print("-" * 40)
        for i, conflict in enumerate(conflict_details[:10]):
            print(f"  {conflict['busco_id']}: {conflict['chromosomes']} → {conflict['assigned_to']} "
                  f"({conflict['status']}, confidence: {conflict['confidence']:.2f})")
        
        if len(conflict_details) > 10:
            print(f"  ... and {len(conflict_details) - 10} more conflicts")
    
    # Show first few genes with chromosome assignments
    print(f"\nFIRST 10 GENES:")
    print("-" * 30)
    for i, (_, row) in enumerate(busco_tsv.head(10).iterrows()):
        assignment_info = busco_df.iloc[i]
        print(f"  {row['sequence']} | {row['busco_id']} | {row['strand']} | "
              f"{row['gene_start']:,}-{row['gene_end']:,} | {assignment_info['assignment_status']}")
    
    return busco_tsv


def main():
    parser = argparse.ArgumentParser(description='Extract AGORA ancestral genome to BUSCO TSV with chromosome mapping')
    parser.add_argument('agora_file', help='AGORA ancestral genome file (.bz2 or plain text)')
    parser.add_argument('busco_files', nargs='+', help='Original BUSCO TSV files used as AGORA input')
    parser.add_argument('output_tsv', help='Output BUSCO TSV file')
    
    args = parser.parse_args()
    
    try:
        extract_agora_ancestral_genome_to_busco(args.agora_file, args.busco_files, args.output_tsv)
        print("\n✅ AGORA ancestral genome extraction completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())