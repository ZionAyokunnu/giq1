"""
Quality assessment module for the Genome Inversion analyser1

"""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def assess_assembly_quality(fasta_path, busco_df, config):
    """Comprehensive assembly quality assessment"""
    
    quality_metrics = {}
    # BUSCO completeness
    if len(busco_df) > 0:
        complete_buscos = len(busco_df[busco_df['status'] == 'Complete'])
        fragmented_buscos = len(busco_df[busco_df['status'] == 'Fragmented'])
        missing_buscos = len(busco_df[busco_df['status'] == 'Missing'])
        duplicated_buscos = len(busco_df[busco_df['status'] == 'Duplicated'])
        total_buscos = len(busco_df)
        
        completeness = complete_buscos / total_buscos if total_buscos > 0 else 0
        fragmentation = fragmented_buscos / total_buscos if total_buscos > 0 else 0
        duplication = duplicated_buscos / total_buscos if total_buscos > 0 else 0
        
        quality_metrics = {
            'busco_completeness': completeness,
            'busco_fragmentation': fragmentation,
            'busco_duplication': duplication,
            'busco_missing': missing_buscos / total_buscos if total_buscos > 0 else 0
        }

    
    
    logger.info(f"  Assembly quality: {quality_metrics}")
    
    return {
        'metrics': quality_metrics,
    }