#!/usr/bin/env python3

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


def parse_agp_alignments(agp_file_path, min_alignment_length=5000000):  # 5MB threshold
    """Parse AGP file to extract chromosome alignments with contig filtering"""
    alignments = defaultdict(list)
    
    if not os.path.exists(agp_file_path):
        raise FileNotFoundError(f"AGP file not found: {agp_file_path}")
    
    filtered_count = 0
    total_count = 0
    
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
                
                total_count += 1
                
                # Filter out small contigs
                if alignment_length >= min_alignment_length:
                    alignments[query_chr].append((ref_chr, alignment_length))
                else:
                    filtered_count += 1
    
    print(f"  Filtered out {filtered_count}/{total_count} small alignments (< {min_alignment_length/1000000:.1f}MB)")
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
        'ragtagg.py', 'scaffold',
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


def find_all_homologous_groups(all_pairwise_alignments=None, genome_fastas=None, busco_file_mapping=None):
    """Find homologous chromosome groups - with RagTag or gene content fallback"""
    
    print("\n" + "=" * 60)
    print("DETECTING HOMOLOGOUS CHROMOSOME GROUPS")
    print("=" * 60)
    
    # Initialize empty homology graph
    homology_graph = defaultdict(set)
    alignment_count = 0
    
    # Method 1: RagTag-based (if available and successful)
    if all_pairwise_alignments:
        print("Using RagTag alignments...")
        processed_pairs = set()
        
        for alignment_key, alignment_data in all_pairwise_alignments.items():
            query_name = alignment_data['query_genome']
            ref_name = alignment_data['reference_genome']
            mappings = alignment_data['mappings']['query_to_reference']
            
            for query_chr, ref_chr in mappings.items():
                pair = tuple(sorted([(query_name, query_chr), (ref_name, ref_chr)]))
                
                if pair not in processed_pairs:
                    homology_graph[(query_name, query_chr)].add((ref_name, ref_chr))
                    homology_graph[(ref_name, ref_chr)].add((query_name, query_chr))
                    processed_pairs.add(pair)
                    alignment_count += 1
        
        print(f"Built homology graph with {alignment_count} RagTag chromosome pairs")
    
    # Method 2: Gene content fallback (if RagTag failed or as enhancement)
    if genome_fastas and busco_file_mapping:
        if alignment_count == 0:
            print("No RagTag data - using gene content as PRIMARY method")
            homology_graph = build_gene_content_homology_primary(homology_graph, genome_fastas, busco_file_mapping)
        else:
            print("Enhancing RagTag with gene content...")
            homology_graph = enhance_homology_with_gene_content(homology_graph, genome_fastas, busco_file_mapping)
    
    # Filter by gene content
    if busco_file_mapping:
        homology_graph = filter_chromosomes_by_gene_count(homology_graph, busco_file_mapping, min_genes_threshold=300)
    
    if len(homology_graph) == 0:
        raise ValueError("No homology data available - need either RagTag alignments or BUSCO files")
    
    # Find connected components (rest unchanged)
    homologous_groups = find_connected_components(homology_graph)
    if busco_file_mapping:
        homology_graph = filter_chromosomes_by_gene_count(homology_graph, busco_file_mapping, min_genes_threshold=300)
        
    # NEW: Enhance with gene content
    if genome_fastas and busco_file_mapping:
        homology_graph = enhance_homology_with_gene_content(homology_graph, genome_fastas, busco_file_mapping)
    
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
    
    if len(homologous_groups) > 400:
        print(f"  ... and {len(homologous_groups) - 10} more groups")
    
    return homologous_groups

def build_gene_content_homology_primary(homology_graph, genome_fastas, busco_file_mapping):
    """Build homology graph based purely on gene content (fallback method)"""
    print("Building homology graph from gene content...")
    
    # Load BUSCO data
    busco_data = {}
    for fasta_path in genome_fastas:
        genome_name = get_genome_name_from_fasta(fasta_path)
        if genome_name in busco_file_mapping:
            try:
                from core import parse_busco_table, filter_busco_genes
                from config import CONFIG
                
                busco_df = parse_busco_table(busco_file_mapping[genome_name], CONFIG)
                filtered_df = filter_busco_genes(busco_df, CONFIG)
                busco_data[genome_name] = filtered_df
                print(f"  Loaded {len(filtered_df)} genes for {genome_name}")
            except Exception as e:
                print(f"  Warning: Could not load BUSCO data for {genome_name}: {e}")
    
    # Build homology from gene content
    connection_count = 0
    gene_content_threshold = 50  # Lower threshold for primary method
    
    genome_names = list(busco_data.keys())
    for i in range(len(genome_names)):
        for j in range(i + 1, len(genome_names)):
            genome1, genome2 = genome_names[i], genome_names[j]
            df1, df2 = busco_data[genome1], busco_data[genome2]
            
            chrs1 = df1['sequence'].unique()
            chrs2 = df2['sequence'].unique()
            
            # Score all chromosome pairs
            chromosome_scores = []
            for chr1 in chrs1:
                for chr2 in chrs2:
                    genes1 = set(df1[df1['sequence'] == chr1]['busco_id'])
                    genes2 = set(df2[df2['sequence'] == chr2]['busco_id'])
                    overlap = len(genes1 & genes2)
                    
                    if overlap >= gene_content_threshold:
                        chromosome_scores.append((chr1, chr2, overlap))
            
            # Add best matches to homology graph
            chromosome_scores.sort(key=lambda x: x[2], reverse=True)
            used_chr1, used_chr2 = set(), set()
            
            for chr1, chr2, overlap in chromosome_scores:
                if chr1 not in used_chr1 and chr2 not in used_chr2:
                    node1 = (genome1, chr1)
                    node2 = (genome2, chr2)
                    
                    homology_graph[node1].add(node2)
                    homology_graph[node2].add(node1)
                    connection_count += 1
                    print(f"  Primary: {genome1}:{chr1} ↔ {genome2}:{chr2} ({overlap} shared genes)")
                    
                    used_chr1.add(chr1)
                    used_chr2.add(chr2)
    
    print(f"Built primary gene content homology with {connection_count} connections")
    return homology_graph

def enhance_homology_with_gene_content(homology_graph, genome_fastas, busco_file_mapping=None):
    """Enhance RagTag homology with gene content validation"""
    if not busco_file_mapping:
        print("No BUSCO data provided - skipping gene content enhancement")
        return homology_graph
    
    print("Enhancing homology with gene content analysis...")
    
    # Load BUSCO data for each genome
    busco_data = {}
    for fasta_path in genome_fastas:
        genome_name = get_genome_name_from_fasta(fasta_path)
        if genome_name in busco_file_mapping:
            busco_file = busco_file_mapping[genome_name]
            try:
                # Parse BUSCO file (you'd need to import this from core)
                from core import parse_busco_table, filter_busco_genes
                from config import CONFIG
                
                busco_df = parse_busco_table(busco_file, CONFIG)
                filtered_df = filter_busco_genes(busco_df, CONFIG)
                busco_data[genome_name] = filtered_df
                print(f"  Loaded {len(filtered_df)} genes for {genome_name}")
            except Exception as e:
                print(f"  Warning: Could not load BUSCO data for {genome_name}: {e}")
    
    if len(busco_data) < 2:
        print("Insufficient BUSCO data - skipping enhancement")
        return homology_graph
    
    # Find high-confidence gene-content-based connections
    enhanced_connections = 0
    gene_content_threshold = 400  # Minimum shared genes
    
    genome_names = list(busco_data.keys())
    for i in range(len(genome_names)):
        for j in range(i + 1, len(genome_names)):
            genome1, genome2 = genome_names[i], genome_names[j]
            df1, df2 = busco_data[genome1], busco_data[genome2]
            
            # Get all chromosomes for each genome
            chrs1 = df1['sequence'].unique()
            chrs2 = df2['sequence'].unique()
            
            # Score all chromosome pairs by gene overlap
            chromosome_scores = []
            for chr1 in chrs1:
                for chr2 in chrs2:
                    genes1 = set(df1[df1['sequence'] == chr1]['busco_id'])
                    genes2 = set(df2[df2['sequence'] == chr2]['busco_id'])
                    overlap = len(genes1 & genes2)
                    
                    if overlap >= gene_content_threshold:
                        chromosome_scores.append((chr1, chr2, overlap))
            
            # Sort by overlap and take best matches
            chromosome_scores.sort(key=lambda x: x[2], reverse=True)
            
            # Add high-confidence connections to homology graph
            used_chr1, used_chr2 = set(), set()
            for chr1, chr2, overlap in chromosome_scores:
                if chr1 not in used_chr1 and chr2 not in used_chr2:
                    # Check if this connection already exists from RagTag
                    node1 = (genome1, chr1)
                    node2 = (genome2, chr2)
                    
                    if node2 not in homology_graph[node1]:
                        # New connection found by gene content
                        homology_graph[node1].add(node2)
                        homology_graph[node2].add(node1)
                        enhanced_connections += 1
                        print(f"  Enhanced: {genome1}:{chr1} ↔ {genome2}:{chr2} ({overlap} shared genes)")
                    
                    used_chr1.add(chr1)
                    used_chr2.add(chr2)
    
    print(f"Enhanced homology graph with {enhanced_connections} additional connections")
    return homology_graph

def filter_chromosomes_by_gene_count(homology_graph, busco_file_mapping, min_genes_threshold=400):
    """Simple filtering: remove chromosomes with < min_genes_threshold BUSCO genes"""
    if not busco_file_mapping:
        return homology_graph
    
    print(f"Filtering chromosomes with < {min_genes_threshold} BUSCO genes...")
    
    # Load gene counts per chromosome
    chr_gene_counts = {}
    for genome_name, busco_file in busco_file_mapping.items():
        try:
            from core import parse_busco_table, filter_busco_genes
            from config import CONFIG
            
            busco_df = parse_busco_table(busco_file, CONFIG)
            filtered_df = filter_busco_genes(busco_df, CONFIG)
            
            # Count genes per chromosome
            gene_counts = filtered_df['sequence'].value_counts()
            for chr_name, count in gene_counts.items():
                chr_gene_counts[(genome_name, chr_name)] = count
        except Exception as e:
            print(f"Warning: Could not process {genome_name}: {e}")
    
    # Filter homology graph
    filtered_graph = defaultdict(set)
    filtered_count = 0
    
    for node, neighbors in homology_graph.items():
        # Check if this chromosome has enough genes
        if chr_gene_counts.get(node, 0) >= min_genes_threshold:
            # Keep this chromosome and its valid connections
            for neighbor in neighbors:
                if chr_gene_counts.get(neighbor, 0) >= min_genes_threshold:
                    filtered_graph[node].add(neighbor)
        else:
            filtered_count += 1
    
    print(f"Filtered out {filtered_count} chromosomes with insufficient genes")
    return filtered_graph

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


def create_unified_mappings_multi_reference(all_pairwise_alignments, genome_fastas, busco_file_mapping=None):
    """Create unified mappings from multi-directional alignments"""
    
    total_genomes = len(genome_fastas)
    
    # Find homologous groups
    homologous_groups = find_all_homologous_groups(all_pairwise_alignments, genome_fastas, busco_file_mapping)
    
    # Assign hierarchical standard names
    standard_names = assign_hierarchical_standard_names(homologous_groups, total_genomes)
    
    # Build unified mappings
    unified_mappings = {
        'genome_mappings': defaultdict(dict),
        'standard_chromosomes': set(),
        'homology_levels': {
            'complete': set(),        # ← Change to sets to avoid duplicates
            'n_minus_1': set(),       
            'n_minus_2': set(),       
            'other_levels': set(),    
            'singletons': set()       
        },
        'statistics': {
            'total_genomes': total_genomes,
            'total_groups': len(homologous_groups),
            'total_mapped_chromosomes': len(standard_names)
        }
    }
    
    # Populate mappings and categorize by level
    for (genome_name, chr_name), standard_name in standard_names.items():
        
        unified_mappings['genome_mappings'][genome_name][chr_name] = standard_name
        unified_mappings['standard_chromosomes'].add(standard_name)
        
        # Categorize by level (add to sets, not lists)
        if 'singleton' in standard_name:
            unified_mappings['homology_levels']['singletons'].add(standard_name)
        elif standard_name.endswith('a'):
            unified_mappings['homology_levels']['n_minus_1'].add(standard_name)
        elif standard_name.endswith('b'):
            unified_mappings['homology_levels']['n_minus_2'].add(standard_name)
        elif any(c.isalpha() for c in standard_name[-1:]):
            unified_mappings['homology_levels']['other_levels'].add(standard_name)
        else:
            unified_mappings['homology_levels']['complete'].add(standard_name)
    
    # Convert sets back to sorted lists for JSON compatibility
    for level in unified_mappings['homology_levels']:
        unified_mappings['homology_levels'][level] = sorted(list(unified_mappings['homology_levels'][level]))
    
    # Add any unmapped chromosomes from original FASTA files
    for fasta_path in genome_fastas:
        genome_name = get_genome_name_from_fasta(fasta_path)
        all_chromosomes = parse_fasta_headers(fasta_path)
        mapped_chromosomes = set(unified_mappings['genome_mappings'][genome_name].keys())
        
        unmapped = all_chromosomes - mapped_chromosomes
        if unmapped:
            print(f"Processing {len(unmapped)} unmapped chromosomes for {genome_name}")
            
            # Filter unmapped chromosomes by gene content
            if busco_file_mapping and genome_name in busco_file_mapping:
                unmapped = filter_unmapped_by_gene_count(unmapped, genome_name, busco_file_mapping[genome_name])
            
            print(f"Adding {len(unmapped)} unmapped chromosomes after filtering")
            for chr_name in sorted(unmapped):
                singleton_name = f"chr{len(standard_names)+1}_singleton_{genome_name}"
                unified_mappings['genome_mappings'][genome_name][chr_name] = singleton_name
                unified_mappings['standard_chromosomes'].add(singleton_name)
                unified_mappings['homology_levels']['singletons'].append(singleton_name)
    
    return unified_mappings


def filter_unmapped_by_gene_count(unmapped_chromosomes, genome_name, busco_file, min_genes=400):
    """Filter unmapped chromosomes by gene content"""
    try:
        from core import parse_busco_table, filter_busco_genes
        from config import CONFIG
        
        busco_df = parse_busco_table(busco_file, CONFIG)
        filtered_df = filter_busco_genes(busco_df, CONFIG)
        
        # Count genes per chromosome
        gene_counts = filtered_df['sequence'].value_counts()
        
        # Filter chromosomes
        filtered_chromosomes = set()
        for chr_name in unmapped_chromosomes:
            gene_count = gene_counts.get(chr_name, 0)
            if gene_count >= min_genes:
                filtered_chromosomes.add(chr_name)
            else:
                print(f"  Filtered out {chr_name}: {gene_count} genes < {min_genes}")
        
        print(f"  Kept {len(filtered_chromosomes)}/{len(unmapped_chromosomes)} unmapped chromosomes")
        return filtered_chromosomes
        
    except Exception as e:
        print(f"  Warning: Could not filter unmapped chromosomes for {genome_name}: {e}")
        return unmapped_chromosomes


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

