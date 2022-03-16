# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .stokespy import StokesCube, StokesParamCube, StokesParamMap, StokesProfile, MagVectorCube 

from .instload import load_HinodeSP_stokes, load_HinodeSP_magvec, load_HMI_stokes, load_HMI_magvec

try:
    from .version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ['StokesCube', 'MagVectorCube', 'load_HinodeSP_stokes', 'load_HinodeSP_magvec', 'load_HMI_stokes', 'load_HMI_magvec']
