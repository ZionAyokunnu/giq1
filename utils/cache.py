"""
Caching utilities for the Genome Inversion analyser
Handles caching of alignment results to avoid recomputation and save the life of my cute little Mac
"""

import hashlib
import pickle
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def generate_cache_key(first_busco_df, second_busco_df, config):
    """Generate a cache key based on input data and configuration"""
    key_config = {
        #Empty till further computation
    }

    return


    return None