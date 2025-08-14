#!/usr/bin/env python3
"""
Chromosome Alignment Parser for GIQ Pipeline
Dynamically parses RagTag AGP files to create chromosome mappings
"""


import os
import sys
import json
import subprocess
from collections import defaultdict
from pathlib import Path




def parse_fasta_headers(fasta_path):
   """
   Extract chromosome names from FASTA file headers
   Returns: set of chromosome names
   """
   chromosomes = set()
  
   if not os.path.exists(fasta_path):
       raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
  
   with open(fasta_path, 'r') as f:
       for line in f:
           if line.startswith('>'):
               # Extract chromosome name (first part of header)
               chr_name = line.strip()[1:].split()[0]
               chromosomes.add(chr_name)
  
   return chromosomes




def parse_agp_alignments(agp_file_path):
   """
   Parse AGP file to extract chromosome alignments
   Returns: dict of {query_chr: [(ref_chr, alignment_length), ...]}
   """
   alignments = defaultdict(list)
  
   if not os.path.exists(agp_file_path):
       raise FileNotFoundError(f"AGP file not found: {agp_file_path}")
  
   with open(agp_file_path, 'r') as f:
       for line in f:
           # Skip comments and headers
           if line.startswith('#') or line.startswith('##'):
               continue
          
           parts = line.strip().split('\t')
           if len(parts) < 9:
               continue
          
           # Only process actual sequence alignments (W = sequence, U = gap)
           if parts[4] == 'W':
               ref_chr = parts[0].replace('_RagTag', '')  # Remove RagTag suffix
               query_chr = parts[5]
               start_pos = int(parts[1])
               end_pos = int(parts[2])
               alignment_length = end_pos - start_pos + 1
              
               alignments[query_chr].append((ref_chr, alignment_length))
  
   return alignments




def create_chromosome_mappings(agp_file_path, reference_fasta_path, query_fasta_path):
   """
   Create chromosome mappings between query and reference genomes
  
   Args:
       agp_file_path: Path to RagTag AGP file
       reference_fasta_path: Path to reference genome FASTA
       query_fasta_path: Path to query genome FASTA
  
   Returns:
       dict: {
           'query_to_standard': {query_chr: standard_name, ...},
           'reference_to_standard': {ref_chr: standard_name, ...}
       }
   """
  
   # Parse FASTA files to get chromosome lists
   ref_chromosomes = parse_fasta_headers(reference_fasta_path)
   query_chromosomes = parse_fasta_headers(query_fasta_path)
  
   print(f"Reference chromosomes: {len(ref_chromosomes)}")
   print(f"Query chromosomes: {len(query_chromosomes)}")
  
   # Parse AGP alignments
   alignments = parse_agp_alignments(agp_file_path)
  
   # Create mapping from query chromosomes to reference chromosomes
   # Use the largest alignment for primary mapping
   query_to_ref = {}
   for query_chr, alignment_list in alignments.items():
       if alignment_list:
           # Sort by alignment length (descending) and take the largest
           largest_alignment = max(alignment_list, key=lambda x: x[1])
           query_to_ref[query_chr] = largest_alignment[0]
  
   # Group reference chromosomes and assign standard names
   ref_to_standard = {}
   query_to_standard = {}
  
   # Get unique reference chromosomes and sort them
   unique_ref_chrs = sorted(set(query_to_ref.values()))
  
   # Assign standard chromosome names (chr1, chr2, etc.)
   for i, ref_chr in enumerate(unique_ref_chrs, 1):
       standard_name = f"chr{i}"
       ref_to_standard[ref_chr] = standard_name
  
   # Map query chromosomes to standard names
   ref_chr_usage = defaultdict(list)  # Track which query chrs map to each ref chr
  
   for query_chr, ref_chr in query_to_ref.items():
       ref_chr_usage[ref_chr].append(query_chr)
  
   # Assign standard names, handling multiple query chrs per reference chr
   for ref_chr, query_chr_list in ref_chr_usage.items():
       standard_base = ref_to_standard[ref_chr]
      
       if len(query_chr_list) == 1:
           # Single mapping - use base name
           query_to_standard[query_chr_list[0]] = standard_base
       else:
           # Multiple mappings - sort by alignment length and assign suffixes
           query_alignments = []
           for query_chr in query_chr_list:
               # Get total alignment length for this query chromosome
               total_length = sum(length for ref, length in alignments[query_chr] if ref == ref_chr)
               query_alignments.append((query_chr, total_length))
          
           # Sort by alignment length (descending)
           query_alignments.sort(key=lambda x: x[1], reverse=True)
          
           # Assign names: largest gets base name, others get suffixes
           for i, (query_chr, _) in enumerate(query_alignments):
               if i == 0:
                   query_to_standard[query_chr] = standard_base
               else:
                   suffix = chr(ord('b') + i - 1)  # b, c, d, ...
                   query_to_standard[query_chr] = f"{standard_base}{suffix}"
  
   # Add reference chromosomes that don't appear in AGP but exist in FASTA
   for ref_chr in ref_chromosomes:
       if ref_chr not in ref_to_standard:
           # Assign next available chromosome number
           next_num = len(unique_ref_chrs) + 1
           ref_to_standard[ref_chr] = f"chr{next_num}"
           unique_ref_chrs.append(ref_chr)
  
   print(f"Created mappings for {len(query_to_standard)} query chromosomes")
   print(f"Created mappings for {len(ref_to_standard)} reference chromosomes")
  
   return {
       'query_to_standard': query_to_standard,
       'reference_to_standard': ref_to_standard,
       'alignment_stats': {
           'total_alignments': sum(len(aligns) for aligns in alignments.values()),
           'mapped_chromosomes': len(query_to_standard)
       }
   }




def standardize_chromosome_name(mappings, genome_type, chr_name):
   """
   Convert chromosome name to standardized name using mappings
  
   Args:
       mappings: Output from create_chromosome_mappings()
       genome_type: 'reference' or 'query'
       chr_name: Original chromosome name
  
   Returns:
       str: Standardized chromosome name (chr1, chr2, etc.)
   """
   if genome_type == 'reference':
       return mappings['reference_to_standard'].get(chr_name, chr_name)
   elif genome_type == 'query':
       return mappings['query_to_standard'].get(chr_name, chr_name)
   else:
       raise ValueError("genome_type must be 'reference' or 'query'")




def print_mapping_summary(mappings):
   """
   Print a summary of chromosome mappings
   """
   print("\n=== Chromosome Mapping Summary ===")
  
   print("\nReference to Standard Mappings:")
   for ref_chr, std_name in sorted(mappings['reference_to_standard'].items()):
       print(f"  {ref_chr} → {std_name}")
  
   print("\nQuery to Standard Mappings:")
   for query_chr, std_name in sorted(mappings['query_to_standard'].items()):
       print(f"  {query_chr} → {std_name}")
  
   print(f"\nAlignment Statistics:")
   stats = mappings['alignment_stats']
   print(f"  Total alignments: {stats['total_alignments']}")
   print(f"  Mapped chromosomes: {stats['mapped_chromosomes']}")




def run_ragtag_scaffold(reference_fasta, query_fasta, output_dir):
   """
   Run RagTag scaffold command and return AGP file path
  
   Args:
       reference_fasta: Path to reference genome FASTA
       query_fasta: Path to query genome FASTA 
       output_dir: Directory for RagTag output
  
   Returns:
       str: Path to generated AGP file
   """
   # Create output directory if it doesn't exist
   os.makedirs(output_dir, exist_ok=True)
  
   # Run RagTag scaffold command
   cmd = [
       'ragtag.py', 'scaffold',
       reference_fasta, query_fasta,
       '-o', output_dir
   ]
  
   print(f"Running RagTag: {' '.join(cmd)}")
  
   try:
       result = subprocess.run(cmd, capture_output=True, text=True, check=True)
       print(f"RagTag completed successfully for {os.path.basename(query_fasta)}")
   except subprocess.CalledProcessError as e:
       raise RuntimeError(f"RagTag failed for {query_fasta}: {e.stderr}")
  
   # Return path to AGP file
   agp_file = os.path.join(output_dir, 'ragtag.scaffold.agp')
   if not os.path.exists(agp_file):
       raise FileNotFoundError(f"Expected AGP file not found: {agp_file}")
  
   return agp_file




def create_mappings_from_fastas(reference_fasta, query_fasta_list, output_base_dir):
   """
   Complete pipeline: Run RagTag for multiple genomes and create chromosome mappings
  
   Args:
       reference_fasta: Path to reference genome FASTA (e.g., Bibio_marci.fna)
       query_fasta_list: List of paths to query genome FASTAs
       output_base_dir: Base directory for all RagTag outputs
  
   Returns:
       dict: {genome_name: mappings_dict, ...} for all genomes
   """
   all_mappings = {}
  
   print(f"=== Running RagTag Pipeline ===")
   print(f"Reference: {reference_fasta}")
   print(f"Query genomes: {len(query_fasta_list)}")
   print(f"Output directory: {output_base_dir}")
  
   # Create base output directory
   os.makedirs(output_base_dir, exist_ok=True)
  
   for query_fasta in query_fasta_list:
       # Extract genome name from filename
       genome_name = os.path.basename(query_fasta).replace('.fna', '').replace('.fasta', '')
       output_dir = os.path.join(output_base_dir, f"{genome_name}_vs_reference")
      
       print(f"\n--- Processing {genome_name} ---")
      
       try:
           # Run RagTag scaffold
           agp_file = run_ragtag_scaffold(reference_fasta, query_fasta, output_dir)
          
           # Create chromosome mappings
           mappings = create_chromosome_mappings(agp_file, reference_fasta, query_fasta)
           all_mappings[genome_name] = mappings
          
           print(f"✓ Successfully created mappings for {genome_name}")
          
       except Exception as e:
           print(f"✗ Failed to process {genome_name}: {e}")
           continue
  
   print(f"\n=== Pipeline Complete ===")
   print(f"Successfully processed {len(all_mappings)} genomes")
  
   return all_mappings




def get_unified_chromosome_mappings(all_mappings):
   """
   Create unified mappings across all genomes using reference genome as standard
  
   Args:
       all_mappings: Output from create_mappings_from_fastas()
  
   Returns:
       dict: Unified mappings for use in GIQ pipeline
   """
   unified_mappings = {
       'genome_mappings': {},  # {genome_name: {chr_name: standard_name}}
       'standard_chromosomes': set()  # All standard chromosome names
   }
  
   for genome_name, mappings in all_mappings.items():
       # Store mappings for this genome
       unified_mappings['genome_mappings'][genome_name] = mappings['query_to_standard'].copy()
      
       # Add to set of all standard chromosomes
       unified_mappings['standard_chromosomes'].update(mappings['query_to_standard'].values())
  
   # Also include reference mappings (if we have them)
   if all_mappings:
       # Get reference mappings from first genome
       first_mappings = next(iter(all_mappings.values()))
       reference_mappings = first_mappings['reference_to_standard']
      
       # Add reference genome mappings (extract reference name from function usage)
       unified_mappings['reference_mappings'] = reference_mappings
       unified_mappings['standard_chromosomes'].update(reference_mappings.values())
  
   return unified_mappings




def save_unified_mappings(unified_mappings, output_dir):
   """
   Save unified mappings to JSON file for use in build-profile
  
   Args:
       unified_mappings: Output from get_unified_chromosome_mappings()
       output_dir: Directory to save the mappings file
  
   Returns:
       str: Path to saved mappings file
   """
   mappings_file = os.path.join(output_dir, 'chromosome_mappings.json')
  
   # Convert sets to lists for JSON serialization
   serializable_mappings = {
       'genome_mappings': unified_mappings['genome_mappings'],
       'standard_chromosomes': list(unified_mappings['standard_chromosomes'])
   }
  
   if 'reference_mappings' in unified_mappings:
       serializable_mappings['reference_mappings'] = unified_mappings['reference_mappings']
  
   with open(mappings_file, 'w') as f:
       json.dump(serializable_mappings, f, indent=2)
  
   print(f"Chromosome mappings saved to: {mappings_file}")
   return mappings_file




def load_chromosome_mappings(mappings_file):
   """
   Load chromosome mappings from JSON file
  
   Args:
       mappings_file: Path to chromosome mappings JSON file
  
   Returns:
       dict: Unified mappings dictionary
   """
   with open(mappings_file, 'r') as f:
       mappings = json.load(f)
  
   # Convert standard_chromosomes back to set
   mappings['standard_chromosomes'] = set(mappings['standard_chromosomes'])
  
   return mappings




def standardize_chromosome_name_unified(unified_mappings, genome_name, chr_name):
   """
   Standardize chromosome name using unified mappings from pipeline
  
   Args:
       unified_mappings: Output from get_unified_chromosome_mappings()
       genome_name: Name of genome (e.g., 'Dilophus_febrilis')
       chr_name: Original chromosome name
  
   Returns:
       str: Standard chromosome name
   """
   genome_mappings = unified_mappings['genome_mappings'].get(genome_name, {})
   return genome_mappings.get(chr_name, chr_name)




def align_chromosomes_command(reference_fasta, query_fastas, output_dir):
   """
   Main command function for chromosome alignment
  
   Args:
       reference_fasta: Path to reference genome FASTA
       query_fastas: List of paths to query genome FASTAs
       output_dir: Output directory for alignment results
  
   Returns:
       str: Path to chromosome mappings file
   """
   print("=" * 60)
   print("CHROMOSOME ALIGNMENT PIPELINE")
   print("=" * 60)
  
   # Create output directory structure
   os.makedirs(output_dir, exist_ok=True)
   ragtag_dir = os.path.join(output_dir, 'ragtag_outputs')
   stages_dir = os.path.join(output_dir, 'stages')
   os.makedirs(stages_dir, exist_ok=True)
  
   # Run RagTag pipeline
   all_mappings = create_mappings_from_fastas(reference_fasta, query_fastas, ragtag_dir)
  
   # Create unified mappings
   unified_mappings = get_unified_chromosome_mappings(all_mappings)
  
   # Save mappings file for build-profile
   mappings_file = save_unified_mappings(unified_mappings, output_dir)
  
   # Save stage data for clarity
   import json
  
   # Stage 1: Save individual genome mappings
   for genome_name, mappings in all_mappings.items():
       stage_file = os.path.join(stages_dir, f'1_mappings_{genome_name}.json')
       with open(stage_file, 'w') as f:
           json.dump({
               'query_to_standard': mappings['query_to_standard'],
               'reference_to_standard': mappings['reference_to_standard'],
               'alignment_stats': mappings['alignment_stats']
           }, f, indent=2)
       print(f"Saved individual mappings: {stage_file}")
  
   # Stage 2: Save unified mappings
   unified_stage_file = os.path.join(stages_dir, '2_unified_mappings.json')
   with open(unified_stage_file, 'w') as f:
       json.dump({
           'genome_mappings': unified_mappings['genome_mappings'],
           'standard_chromosomes': list(unified_mappings['standard_chromosomes']),
           'reference_mappings': unified_mappings.get('reference_mappings', {})
       }, f, indent=2)
   print(f"Saved unified mappings: {unified_stage_file}")
  
   # Stage 3: Save alignment summary
   summary = {
       'reference_genome': os.path.basename(reference_fasta),
       'query_genomes': [os.path.basename(f) for f in query_fastas],
       'total_genomes': len(query_fastas) + 1,
       'standard_chromosomes': list(unified_mappings['standard_chromosomes']),
       'mappings_file': mappings_file,
       'ragtag_outputs': ragtag_dir
   }
  
   summary_file = os.path.join(output_dir, 'alignment_summary.json')
   with open(summary_file, 'w') as f:
       json.dump(summary, f, indent=2)
  
   print("\n" + "=" * 60)
   print("ALIGNMENT PIPELINE COMPLETE")
   print("=" * 60)
   print(f"Reference genome: {os.path.basename(reference_fasta)}")
   print(f"Query genomes: {len(query_fastas)}")
   print(f"Standard chromosomes: {len(unified_mappings['standard_chromosomes'])}")
   print(f"Mappings file: {mappings_file}")
   print(f"Summary file: {summary_file}")
   print("=" * 60)
  
   return mappings_file


