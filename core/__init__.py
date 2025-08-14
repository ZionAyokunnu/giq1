"""Core analysis modules for genome inversion analyser"""

from .quality_assessment import (
    assess_assembly_quality,
)

from .busco_processor import (
    parse_busco_table,
    filter_busco_genes,
    detect_flips,
    correct_strand_orientation,
    process_multiple_genomes,
    
)

from .profile import (
    process_genomes_binning,
    group_genomes_by_chromosome
    
)

from .query_reference import (
    extract_gene_distribution,
    calculate_distribution_stats,
    calculate_bit_scores,
    compare_query_genome_to_profile,
)

from .query_movement import (
    extract_gene_ranges,
    extract_current_ranges,
    filter_same_chromosome,
    calculate_gene_movement,
    analyse_query_movements,
    extract_gene_distribution,
    get_movement_summary,
)

from .reverse import (
    extract_movement_sequence,
    detect_adjacency_inversions,
    detect_extended,
    detect_flip_in_pattern,
    apply_adjacency_inversion,
    apply_flip_inversion,
    iterative_detection,
    get_permutable_positions,
    calculate_position_probability,
    evaluate_inversion_step_probability,
    generate_gene_specific_pathways,
    generate_inversion_steps_for_gene,
    iterative_detection_gene_specific,
    evaluate_pathway_steps,
    get_gene_final_position,
    track_all_gene_final_positions,
    probability_weighted_inversion_analysis,
    check_events_iteration
)

from .chr_aligner import (
    normalize_genome_name,
    parse_fasta_headers,
    parse_agp_alignments,
    create_chromosome_mappings,
    run_ragtag_scaffold,
    run_all_pairwise_alignments,
    save_unified_mappings,
    find_connected_components,
    find_all_homologous_groups,
    get_total_gene_content,
    assign_hierarchical_standard_names,
    create_unified_mappings_multi_reference,
    save_unified_mappings,
    multi_reference_alignment_pipeline,
    get_genome_name_from_fasta,
    load_chromosome_mappings,
    standardize_chromosome_name_unified,
    
    
)

__all__ = [
    'assess_assembly_quality',
    'parse_busco_table',
    'filter_busco_genes',
    'detect_flips',
    'correct_strand_orientation',
    'process_multiple_genomes',
    'process_genomes_binning',
    'extract_gene_distribution',
    'calculate_distribution_stats',
    'calculate_bit_scores',
    'compare_query_genome_to_profile',
    'group_genomes_by_chromosome'
    

]