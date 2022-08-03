__version__ = "0.1.16"

from .xarray_leaflet import LeafletMap  # noqa
from .server_extension import _jupyter_server_extension_paths  # noqa
from .server_extension import _load_jupyter_server_extension
from .server_extension import _jupyter_nbextension_paths  # noqa


load_jupyter_server_extension = _load_jupyter_server_extension
