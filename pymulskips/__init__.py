"""
pymulskips
"""

__author__ = "Gaetano Calogero, Giuseppe Fisicaro, Ioannis Deretzis, Antonino La Magna"
# __license__ = "MPL-2.0"

import warnings
from . import setuprun
from . import analyze
from . import sockets

# import "process" submodule only if dolfin is available
import warnings
try:
    import cantera
    from . import process
except ImportError:
    # cantera is not installed
    warnings.warn("The 'process' submodule could not be imported because the 'cantera' package is not installed. The others were imported successfully.")
    pass

# import "io" submodule only if dolfin is available
try:
    import dolfin
    from . import io
except ImportError:
    # dolfin is not installed
    warnings.warn("The 'io' submodule could not be imported because the 'dolfin' package is not installed or there is some other problem. The others were imported successfully.")
    pass
