"""
Genome Inversion Quantifier - giq2
Modular toolkit for detecting and analysing chromosomal inversions
"""

__version__ = "1.0.0"
__author__ = "Zion Ayokunnu"
__supervisors__ = "Kamil Jaron, Sasha, Arif"
__description__ = "A comprehensive toolkit for genome synteny and inversion analysis"

from . import config
from . import utils
from . import core
from . import visualisation

try:
    from . import registry
except ImportError:
    registry = None

    detect_inversion = None

__all__ = [
    'config',
    'utils', 
    'core',
    'visualisation'
]

if registry is not None:
    __all__.append('registry')

if detect_inversion is not None:
    __all__.append('detect_inversion')