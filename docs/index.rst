Welcome to xarray-leaflet's documentation!
==========================================

**xarray-leaflet** is an **xarray** extension for tiled map plotting. In the
world of web mapping, tiles are chunks of the image you see on your screen,
usually of 256x256 pixels. They are requested by e.g. **Leaflet** when they are
visible in your map view. This allows to handle only parts of the map, which
could be otherwise huge in its totality. It also allows to e.g. reveal more
details as you zoom in.

Chunks also sound familiar to xarray users, because xarray can make use of Dask
and Zarr under the hood to handle parts of data that wouldn't fit in memory
otherwise. The idea with xarray-leaflet is to generate tiles using xarray as
they are requested by Leaflet when the user interacts with the map, either by
dragging it or by zooming. xarray can process only the chunks that are needed to
render each tile, do some aggregation operation on them (e.g. coarsen), or
whatever transformation you want to extract data that is relevant to the current
map view.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   getting_started
   usage
   modules
   contributing
   history

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
