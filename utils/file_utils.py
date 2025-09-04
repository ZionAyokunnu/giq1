

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
