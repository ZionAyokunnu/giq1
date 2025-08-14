#!/usr/bin/env python3
"""
Multi-Directional Chromosome Alignment System for GIQ Pipeline
- Runs RagTag between ALL pairs of genomes (no single reference)
- Finds homologous chromosome groups across all genomes
- Assigns hierarchical standard names: chr1, chr2... then chr11a, chr12a... then chr21b, chr22b...
"""

import os
import sys
import json
import subprocess
from collections import defaultdict, deque
from pathlib import Path
import itertools


def normalize_genome_name(genome_name):
    """Ensure consistent lowercase naming"""
    return genome_name.lower()


def get_genome_name_from_fasta(fasta_path):
    """Extract genome name from FASTA file path"""
    return normalize_genome_name(
        os.path.basename(fasta_path).replace('.fna', '').replace('.fasta', '')
    )


def parse_fasta_headers(fasta_path):
    """Extract chromosome names from FASTA file headers"""
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
    """Parse AGP file to extract chromosome alignments"""
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
    """Create chromosome mappings between query and reference genomes"""
    
    # Parse FASTA files to get chromosome lists
    ref_chromosomes = parse_fasta_headers(reference_fasta_path)
    query_chromosomes = parse_fasta_headers(query_fasta_path)
    
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
    
    return {
        'query_to_reference': query_to_ref,
        'alignment_details': alignments,
        'stats': {
            'total_alignments': sum(len(aligns) for aligns in alignments.values()),
            'mapped_chromosomes': len(query_to_ref),
            'reference_chromosomes': len(ref_chromosomes),
            'query_chromosomes': len(query_chromosomes)
        }
    }


def run_ragtag_scaffold(reference_fasta, query_fasta, output_dir):
    """Run RagTag scaffold command and return AGP file path"""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Run RagTag scaffold command
    cmd = [
        'ragtag.py', 'scaffold',
        reference_fasta, query_fasta,
        '-o', output_dir
    ]
    
    print(f"  Running: ragtag.py scaffold {os.path.basename(reference_fasta)} {os.path.basename(query_fasta)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"RagTag failed: {e.stderr}")
    
    # Return path to AGP file
    agp_file = os.path.join(output_dir, 'ragtag.scaffold.agp')
    if not os.path.exists(agp_file):
        raise FileNotFoundError(f"Expected AGP file not found: {agp_file}")
    
    return agp_file


def run_all_pairwise_alignments(genome_fastas, output_base_dir):
    """Run RagTag between ALL pairs of genomes"""
    
    print("=" * 60)
    print("RUNNING PAIRWISE ALIGNMENTS")
    print("=" * 60)
    
    alignments = {}
    total_pairs = len(genome_fastas) * (len(genome_fastas) - 1)
    current_pair = 0
    
    for ref_fasta in genome_fastas:
        ref_name = get_genome_name_from_fasta(ref_fasta)
        
        for query_fasta in genome_fastas:
            query_name = get_genome_name_from_fasta(query_fasta)
            
            # Skip self-alignment
            if ref_name == query_name:
                continue
            
            current_pair += 1
            print(f"\n[{current_pair}/{total_pairs}] {query_name} vs {ref_name}")
            
            # Create output directory for this alignment
            output_dir = os.path.join(output_base_dir, 'pairwise_alignments', f"{query_name}_vs_{ref_name}")
            
            try:
                # Run RagTag scaffold
                agp_file = run_ragtag_scaffold(ref_fasta, query_fasta, output_dir)
                
                # Create chromosome mappings
                mappings = create_chromosome_mappings(agp_file, ref_fasta, query_fasta)
                
                alignments[f"{query_name}_vs_{ref_name}"] = {
                    'agp_file': agp_file,
                    'mappings': mappings,
                    'query_genome': query_name,
                    'reference_genome': ref_name
                }
                
                print(f"    ✓ Mapped {mappings['stats']['mapped_chromosomes']} chromosomes")
                
            except Exception as e:
                print(f"    ✗ Failed: {e}")
                continue
    
    print(f"\n✓ Completed {len(alignments)} pairwise alignments")

##debug

    print(f"\nDEBUG: All alignments summary:")
    for key, data in alignments.items():
        print(f"  {key}: {len(data['mappings']['query_to_reference'])} mappings")
        print(f"    Sample: {list(data['mappings']['query_to_reference'].items())[:3]}")
        
        
    return alignments
    
def find_connected_components(graph):
    """Find connected components in homology graph using BFS with debugging"""
    visited = set()
    components = []
    
    print(f"DEBUG: Graph has {len(graph)} nodes")
    print(f"DEBUG: Sample nodes: {list(graph.keys())[:5]}")
    
    for node in graph:
        if node not in visited:
            # Start new component
            component = set()
            queue = deque([node])
            
            while queue:
                current = queue.popleft()
                if current not in visited:
                    visited.add(current)
                    component.add(current)
                    
                    # DEBUG: Print current node and its neighbors
                    neighbors = graph[current]
                    if len(component) <= 3:  # Only print for small components
                        print(f"DEBUG: Node {current} has {len(neighbors)} neighbors")
                    
                    # Add all neighbors to queue
                    for neighbor in neighbors:
                        if neighbor not in visited:
                            queue.append(neighbor)
            
            if component:
                components.append(component)
                # DEBUG: Print component info
                genomes_in_component = set(genome for genome, chr in component)
                print(f"DEBUG: Found component with {len(component)} chromosomes from {len(genomes_in_component)} genomes")
                if len(component) <= 10:  # Show details for small components
                    print(f"  Component: {sorted(component)}")
    
    print(f"DEBUG: Total components found: {len(components)}")
    return components


def find_all_homologous_groups(all_pairwise_alignments):
    """Find homologous chromosome groups across all genomes"""
    
    print("\n" + "=" * 60)
    print("DETECTING HOMOLOGOUS CHROMOSOME GROUPS")
    print("=" * 60)
    
    # Build homology graph
    homology_graph = defaultdict(set)
    
    # For each pairwise alignment, add homologous pairs
    alignment_count = 0
    for alignment_key, alignment_data in all_pairwise_alignments.items():
        query_name = alignment_data['query_genome']
        ref_name = alignment_data['reference_genome']
        mappings = alignment_data['mappings']['query_to_reference']
        
        for query_chr, ref_chr in mappings.items():
            # Add bidirectional homology
            homology_graph[(query_name, query_chr)].add((ref_name, ref_chr))
            homology_graph[(ref_name, ref_chr)].add((query_name, query_chr))
            alignment_count += 1
    
    print(f"Built homology graph with {alignment_count} chromosome pairs")
    
    # Find connected components (homologous groups)
    homologous_groups = find_connected_components(homology_graph)
    
    # Sort groups by size (number of genomes) and then by total chromosomes
    homologous_groups.sort(key=lambda g: (-len(set(genome for genome, chr in g)), -len(g)))
    
    print(f"Found {len(homologous_groups)} homologous chromosome groups:")
    for i, group in enumerate(homologous_groups[:10]):  # Show first 10
        genomes_in_group = set(genome for genome, chr in group)
        print(f"  Group {i+1}: {len(genomes_in_group)} genomes, {len(group)} chromosomes")
        if len(group) <= 6:  # Show details for small groups
            for genome, chr in sorted(group):
                print(f"    {genome}:{chr}")
    
    if len(homologous_groups) > 10:
        print(f"  ... and {len(homologous_groups) - 10} more groups")
    
    
    
    ##debug
    
    
    # Add to find_all_homologous_groups():
    print(f"DEBUG: Sample homologous groups:")
    for i, group in enumerate(homologous_groups[:5]):
        print(f"  Group {i}: {group}")
        
        
        return homologous_groups


def get_total_gene_content(group, genome_gene_counts=None):
    """Estimate total gene content for a homologous group (placeholder)"""
    # This is a placeholder - in real implementation, you'd count BUSCO genes
    # For now, just use group size as a proxy
    return len(group)


def assign_hierarchical_standard_names(homologous_groups, total_genomes):
    """Assign hierarchical standard names: chr1, chr2... then chr11a, chr12a... then chr21b, chr22b..."""
    
    print("\n" + "=" * 60)
    print("ASSIGNING HIERARCHICAL STANDARD NAMES")
    print("=" * 60)
    
    standard_names = {}
    chr_counter = 1
    
    # Group by number of genomes represented
    groups_by_genome_count = defaultdict(list)
    for group in homologous_groups:
        genome_count = len(set(genome for genome, chr in group))
        groups_by_genome_count[genome_count].append(group)
    
    # Sort within each genome count by total size
    for genome_count in groups_by_genome_count:
        groups_by_genome_count[genome_count].sort(
            key=lambda g: (-get_total_gene_content(g), -len(g))
        )
    
    # Level 0: Complete groups (all genomes present)
    if total_genomes in groups_by_genome_count:
        complete_groups = groups_by_genome_count[total_genomes]
        print(f"Level 0 (Complete): {len(complete_groups)} groups with all {total_genomes} genomes")
        
        ##Debug
        
        for i, group in enumerate(complete_groups[:3]):  # Show first 3 groups
            print(f"  Group {i}: {group}")
            print(f"    Genomes in group: {set(genome for genome, chr in group)}")
            print(f"    Size: {len(group)} chromosomes")
            
        
        
        for group in complete_groups:
            standard_name = f"chr{chr_counter}"
            for genome_name, chr_name in group:
                standard_names[(genome_name, chr_name)] = standard_name
            print(f"  {standard_name}: {len(group)} chromosomes")
            chr_counter += 1
    
    # Level 1: (n-1) genome groups
    if (total_genomes - 1) in groups_by_genome_count:
        n_minus_1_groups = groups_by_genome_count[total_genomes - 1]
        print(f"Level 1 (n-1): {len(n_minus_1_groups)} groups with {total_genomes-1} genomes")
        
        for group in n_minus_1_groups:
            standard_name = f"chr{chr_counter}a"
            for genome_name, chr_name in group:
                standard_names[(genome_name, chr_name)] = standard_name
            print(f"  {standard_name}: {len(group)} chromosomes")
            chr_counter += 1
    
    # Level 2: (n-2) genome groups
    if (total_genomes - 2) in groups_by_genome_count and total_genomes >= 3:
        n_minus_2_groups = groups_by_genome_count[total_genomes - 2]
        print(f"Level 2 (n-2): {len(n_minus_2_groups)} groups with {total_genomes-2} genomes")
        
        for group in n_minus_2_groups:
            standard_name = f"chr{chr_counter}b"
            for genome_name, chr_name in group:
                standard_names[(genome_name, chr_name)] = standard_name
            print(f"  {standard_name}: {len(group)} chromosomes")
            chr_counter += 1
    
    # Continue pattern for smaller groups
    suffixes = ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
    suffix_index = 0
    
    for genome_count in range(total_genomes - 3, 0, -1):
        if genome_count in groups_by_genome_count and suffix_index < len(suffixes):
            suffix = suffixes[suffix_index]
            groups = groups_by_genome_count[genome_count]
            print(f"Level {suffix}: {len(groups)} groups with {genome_count} genomes")
            
            for group in groups:
                standard_name = f"chr{chr_counter}{suffix}"
                for genome_name, chr_name in group:
                    standard_names[(genome_name, chr_name)] = standard_name
                print(f"  {standard_name}: {len(group)} chromosomes")
                chr_counter += 1
            
            suffix_index += 1
    
    # Remaining singletons (genome-specific chromosomes)
    singleton_counter = chr_counter
    for genome_count in [1]:
        if genome_count in groups_by_genome_count:
            singletons = groups_by_genome_count[genome_count]
            print(f"Singletons: {len(singletons)} genome-specific chromosomes")
            
            for group in singletons:
                # Group should have exactly one chromosome
                if len(group) == 1:
                    genome_name, chr_name = list(group)[0]
                    standard_name = f"chr{singleton_counter}_singleton_{genome_name}"
                    standard_names[(genome_name, chr_name)] = standard_name
                    singleton_counter += 1
    
    print(f"\nAssigned standard names to {len(standard_names)} chromosomes")
    return standard_names


def create_unified_mappings_multi_reference(all_pairwise_alignments, genome_fastas):
    """Create unified mappings from multi-directional alignments"""
    
    total_genomes = len(genome_fastas)
    
    # Find homologous groups
    homologous_groups = find_all_homologous_groups(all_pairwise_alignments)
    
    # Assign hierarchical standard names
    standard_names = assign_hierarchical_standard_names(homologous_groups, total_genomes)
    
    # Build unified mappings
    unified_mappings = {
        'genome_mappings': defaultdict(dict),
        'standard_chromosomes': set(),
        'homology_levels': {
            'complete': [],        # chr1, chr2, chr3...
            'n_minus_1': [],       # chr11a, chr12a, chr13a...
            'n_minus_2': [],       # chr21b, chr22b, chr23b...
            'other_levels': [],    # chr31c, chr41d, etc.
            'singletons': []       # genome-specific chromosomes
        },
        'statistics': {
            'total_genomes': total_genomes,
            'total_groups': len(homologous_groups),
            'total_mapped_chromosomes': len(standard_names)
        }
    }
    
    # Populate mappings and categorize by level
    for (genome_name, chr_name), standard_name in standard_names.items():
        
        ##Debug
        
        if standard_name in ['chr1', 'chr2']:  # Debug first few
          print(f"DEBUG: Adding {standard_name} for {genome_name}:{chr_name}")
        
        ##
        
        unified_mappings['genome_mappings'][genome_name][chr_name] = standard_name
        unified_mappings['standard_chromosomes'].add(standard_name)
        
        # Categorize by level
        if 'singleton' in standard_name:
            unified_mappings['homology_levels']['singletons'].append(standard_name)
        elif standard_name.endswith('a'):
            unified_mappings['homology_levels']['n_minus_1'].append(standard_name)
        elif standard_name.endswith('b'):
            unified_mappings['homology_levels']['n_minus_2'].append(standard_name)
        elif any(c.isalpha() for c in standard_name[-1:]):
            # Other letter suffixes
            unified_mappings['homology_levels']['other_levels'].append(standard_name)
        else:
            # No suffix - complete groups
            unified_mappings['homology_levels']['complete'].append(standard_name)
    
    # Add any unmapped chromosomes from original FASTA files
    for fasta_path in genome_fastas:
        genome_name = get_genome_name_from_fasta(fasta_path)
        all_chromosomes = parse_fasta_headers(fasta_path)
        mapped_chromosomes = set(unified_mappings['genome_mappings'][genome_name].keys())
        
        unmapped = all_chromosomes - mapped_chromosomes
        if unmapped:
            print(f"Adding {len(unmapped)} unmapped chromosomes for {genome_name}")
            for chr_name in sorted(unmapped):
                singleton_name = f"chr{len(standard_names)+1}_singleton_{genome_name}"
                unified_mappings['genome_mappings'][genome_name][chr_name] = singleton_name
                unified_mappings['standard_chromosomes'].add(singleton_name)
                unified_mappings['homology_levels']['singletons'].append(singleton_name)
    
    return unified_mappings


def save_unified_mappings(unified_mappings, output_dir):
    """Save unified mappings to JSON file for use in build-profile"""
    mappings_file = os.path.join(output_dir, 'chromosome_mappings.json')
    
    # Convert sets to lists for JSON serialization
    serializable_mappings = {
        'genome_mappings': dict(unified_mappings['genome_mappings']),
        'standard_chromosomes': sorted(list(unified_mappings['standard_chromosomes'])),
        'homology_levels': {
            level: sorted(chrs) for level, chrs in unified_mappings['homology_levels'].items()
        },
        'statistics': unified_mappings['statistics']
    }
    
    with open(mappings_file, 'w') as f:
        json.dump(serializable_mappings, f, indent=2)
    
    print(f"\nChromosome mappings saved to: {mappings_file}")
    return mappings_file


def load_chromosome_mappings(mappings_file):
    """Load chromosome mappings from JSON file"""
    with open(mappings_file, 'r') as f:
        mappings = json.load(f)
    
    # Convert standard_chromosomes back to set for compatibility
    if 'standard_chromosomes' in mappings and isinstance(mappings['standard_chromosomes'], list):
        mappings['standard_chromosomes'] = set(mappings['standard_chromosomes'])
    
    return mappings


def standardize_chromosome_name_unified(unified_mappings, genome_name, chr_name):
    """
    Standardize chromosome name using unified mappings from multi-directional pipeline
    
    Args:
        unified_mappings: Mappings dictionary from load_chromosome_mappings()
        genome_name: Name of the genome 
        chr_name: Original chromosome name
        
    Returns:
        str: Standardized chromosome name
    """
    # Normalize input genome name to lowercase
    normalized_genome_name = normalize_genome_name(genome_name)
    
    genome_mappings_dict = unified_mappings['genome_mappings']
    
    # Try normalized name first
    if normalized_genome_name in genome_mappings_dict:
        return genome_mappings_dict[normalized_genome_name].get(chr_name, chr_name)
    
    # Fallback: try exact match (for backward compatibility)  
    if genome_name in genome_mappings_dict:
        return genome_mappings_dict[genome_name].get(chr_name, chr_name)
    
    # If no mapping found, return original name
    print(f"WARNING: No mapping found for genome '{genome_name}' (normalized: '{normalized_genome_name}')")
    print(f"Available genomes: {list(genome_mappings_dict.keys())}")
    return chr_name


def multi_reference_alignment_pipeline(genome_fastas, output_dir):
    """Complete multi-directional alignment pipeline"""
    
    print("=" * 80)
    print("MULTI-DIRECTIONAL CHROMOSOME ALIGNMENT PIPELINE")
    print("=" * 80)
    print(f"Input genomes: {len(genome_fastas)}")
    for fasta in genome_fastas:
        print(f"  - {get_genome_name_from_fasta(fasta)}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory structure
    os.makedirs(output_dir, exist_ok=True)
    stages_dir = os.path.join(output_dir, 'stages')
    os.makedirs(stages_dir, exist_ok=True)
    
    # Step 1: Run all pairwise alignments
    print("\nSTEP 1: Running pairwise RagTag alignments...")
    all_alignments = run_all_pairwise_alignments(genome_fastas, output_dir)
    
    # Save stage 1 data
    stage1_file = os.path.join(stages_dir, '1_pairwise_alignments.json')
    with open(stage1_file, 'w') as f:
        # Serialize alignment data (excluding file paths for cleaner JSON)
        serializable_alignments = {}
        for key, data in all_alignments.items():
            serializable_alignments[key] = {
                'query_genome': data['query_genome'],
                'reference_genome': data['reference_genome'],
                'mappings': data['mappings'],
                'agp_file': data['agp_file']
            }
        json.dump(serializable_alignments, f, indent=2)
    print(f"Saved pairwise alignments: {stage1_file}")
    
    # Step 2: Create unified mappings
    print("\nSTEP 2: Creating unified multi-reference mappings...")
    unified_mappings = create_unified_mappings_multi_reference(all_alignments, genome_fastas)
    
    # Step 3: Save results
    print("\nSTEP 3: Saving mappings...")
    mappings_file = save_unified_mappings(unified_mappings, output_dir)
    
    # Save stage 2 data  
    stage2_file = os.path.join(stages_dir, '2_unified_mappings.json')
    with open(stage2_file, 'w') as f:
        json.dump({
            'genome_mappings': dict(unified_mappings['genome_mappings']),
            'standard_chromosomes': sorted(list(unified_mappings['standard_chromosomes'])),
            'homology_levels': unified_mappings['homology_levels'],
            'statistics': unified_mappings['statistics']
        }, f, indent=2)
    print(f"Saved unified mappings: {stage2_file}")
    
    # Save summary
    summary = {
        'pipeline_type': 'multi_directional',
        'input_genomes': [get_genome_name_from_fasta(f) for f in genome_fastas],
        'total_genomes': len(genome_fastas),
        'pairwise_alignments': len(all_alignments),
        'standard_chromosomes': len(unified_mappings['standard_chromosomes']),
        'homology_levels': {
            level: len(chrs) for level, chrs in unified_mappings['homology_levels'].items()
        },
        'mappings_file': mappings_file,
        'output_directory': output_dir
    }
    
    summary_file = os.path.join(output_dir, 'alignment_summary.json')
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n" + "=" * 80)
    print("MULTI-DIRECTIONAL ALIGNMENT PIPELINE COMPLETE")
    print("=" * 80)
    print(f"Total genomes: {len(genome_fastas)}")
    print(f"Pairwise alignments: {len(all_alignments)}")
    print(f"Standard chromosomes: {len(unified_mappings['standard_chromosomes'])}")
    print("\nHierarchy breakdown:")
    for level, chrs in unified_mappings['homology_levels'].items():
        print(f"  {level}: {len(chrs)} chromosomes")
    print(f"\nMappings file: {mappings_file}")
    print(f"Summary file: {summary_file}")
    print("=" * 80)
    
    return mappings_file