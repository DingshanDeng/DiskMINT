"""
 DiskMINT model
 Based RADMC3D, and LIME
 Dingshan Deng
 dingshandeng@arizona.edu
"""

from importlib.metadata import version, PackageNotFoundError
from . import constants
# from . import wrapper_model_parameters
# from . import wrapper_disk_density

# -----------------------------------------------------------
# Version Handling: Reads from installed package metadata
# -----------------------------------------------------------
try:
    __version__ = version("diskmint")
except PackageNotFoundError:
    # This happens if you run the code from source without installing it
    # (e.g. just running python inside the folder)
    __version__ = "unknown (source)"

__author__ = "Dingshan Deng"
__copyright__ = "Copyright (C) 2023-2025 Dingshan Deng"
__all__ = ["constants", "model", "disk_density", "execute", "dustopac", "modelgrid"]