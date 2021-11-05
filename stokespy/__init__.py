# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .stokespy import StokesCube, StokesParamCube, StokesParamMap, StokesProfile, MagVectorCube 

from .instload import get_HMI_data, get_SP_data

try:
    from .version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ['StokesCube', 'MagVectorCube', 'get_HMI_data', 'get_SP_data']
