# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .stokespy import StokesCube, MagVectorCube

try:
    from .version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ['StokesCube', 'MagVectorCube']
