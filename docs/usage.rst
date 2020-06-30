=====
Usage
=====

xarray-leaflet tries hard to make it possible to visualize big data. It will
only use the necessary data to create every tile of the map, and if your data
is stored in a chunked format (e.g. with `Zarr
<https://zarr.readthedocs.io>`_), it will never hold the whole DataArray in
memory, but only parts of it as it is needed by each tile. It will also take
advantage of `Dask <https://dask.org>`_'s lazy evaluation system to generate
each tile in parallel.

Projections
===========

xarray-leaflet currently accepts DataArrays that are in the ``EPSG:4326``
projection (aka ``WGS84``).  This means that the ``x`` and ``y`` coordinates
are expressed in degrees (longitudes and latitudes).  Furthermore, they must be
on a regular grid, e.g. there is always 0.1 between two coordinates (although
it could be different for ``x`` and ``y``).  xarray-leaflet's primary Leaflet
projection is Web Mercator, or ``EPSG:3857``. This is Leaflet's default
projection, and most base maps that you can find on the Internet support this
projection.  This means that xarray-leaflet will be much more efficient when
using this projection, basically because it can slice easily.  However, it can
also use a custom Leaflet projection, e.g. ``EPSG:3413`` which is centered on
the north pole. But to do so it will have to work on the whole DataArray, which
might not fit in memory if it is too big.

Static or dynamic mode
======================

xarray-leaflet has two modes, one is for static maps and the other for dynamic
maps.  By dynamic we mean that you won't see the same thing depending on the
current map view.  For instance, when you zoom in you might want to show
details that would otherwise be hidden.  So the map really adapts to where you
look on the Earth The downside of dynamic maps is that you will see more
flickering as you interact with the map, because tiles have to be refreshed as
soon as you drag or zoom.

On the contrary, static maps won't be as reactive with respect to the location.
Of course they will refine as you zoom in, but they won't change dramatically.
On the other hand, they are rendered really smoothly.

Data pipeline
=============

An automatic data pipeline is provided, which uses default transformations to
show your data.  However, you will probably want to adapt these transformations
to fit your needs.  The data goes through a number of functions that are
chained, one function producing data for the next one.  These functions are:

- ``transform0``: it operates on the input DataArray as a whole, even on the
  parts that are not visible on the map.  This is often where you want to
  normalize your data to keep the values between 0 and 1.  The default
  transformation does exactly that.

- ``transform1``: it operates on the part of the data that is visible on the
  map. Often, this is where you want to apply dynamic transformations.

- ``transform2``: it applies to the data contained in each Leaflet tile before
  reprojection. Reprojection needs your data to fit into memory, so you may
  want to downsample your data and keep approximately the same number of points
  as there are in a tile (256 x 256). The default transformation does just
  that.

- ``transform3``: it applies to the data in each tile after reprojection. This
  is the last opportunity to transform your data before it is saved to a PNG
  file and sent to the browser. That's often when you want to apply styling.

Appart from this data flow, you can also pass a `matplotlib
<https://matplotlib.org/tutorials/colors/colormaps.html>`_ colormap function,
as well as a `rasterio
<https://rasterio.readthedocs.io/en/latest/api/rasterio.warp.html>`_ resampling
method.
