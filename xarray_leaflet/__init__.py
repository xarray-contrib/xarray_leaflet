"""Top-level package for xarray-leaflet."""

__author__ = """David Brochart"""
__email__ = 'david.brochart@gmail.com'
__version__ = '0.1.10'

from .xarray_leaflet import LeafletMap
from .server_extension import _jupyter_server_extension_paths, _load_jupyter_server_extension, _jupyter_nbextension_paths


load_jupyter_server_extension = _load_jupyter_server_extension
