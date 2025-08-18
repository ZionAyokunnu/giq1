"""Core analysis modules for genome inversion analyser"""

from .quality_assessment import (
    assess_assembly_quality,
)

from .busco_processor import (
    parse_busco_table,
    filter_busco_genes,
    detect_inversions
)

from .chr_aligner import (
    parse_fasta_headers,
    parse_agp_alignments,
    create_chromosome_mappings,
    run_ragtag_scaffold,
    run_pairwise_ragtag,
    prepare_diagonal_synteny_data,
    select_optimal_chromosome_pairs
)
__all__ = [
    'assess_assembly_quality',
    'parse_busco_table',
    'filter_busco_genes',
    'detect_inversions'
]