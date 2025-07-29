"""
File utilities for the Genome Inversion analyser1
"""

import pandas as pd
from pathlib import Path


def create_output_directory(config):
    """Create output directory structure"""
    base_dir = Path(config.get('base_output_dir', 'results'))
    base_dir.mkdir(exist_ok=True)
    
    subdirs = ['plots', 'data', 'reports', 'debug', 'cache']
    for subdir in subdirs:
        (base_dir / subdir).mkdir(exist_ok=True)
    
    return base_dir


def standardize_sequence_id(id_str):
    """Standardize sequence IDs to ensure consistency across tools"""
    if pd.isna(id_str) or not id_str:
        return None
    
    clean_id = str(id_str).strip()
    clean_id = clean_id.replace('>', '').split()[0]
    clean_id = clean_id.replace('|', '_').replace(':', '_')
    
    return clean_id