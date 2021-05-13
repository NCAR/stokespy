from .stokespy import StokesCube

try:
    from .version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ['StokesCube']
