[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/davidbrochart/xarray_leaflet/master?filepath=examples%2Fintroduction.ipynb)

# xarray-leaflet: an xarray extension for tiled map plotting

[xarray](http://xarray.pydata.org) and [Leaflet](https://leafletjs.com) share this ability to work with fragments of data, xarray through Dask's chunks, and Leaflet through map tiles. In the end this is really the same concept, so it was a natural thing to make them work together.

Fortunately xarray is written in Python, and we happen to have a great Python binding for Leaflet, [ipyleaflet](https://ipyleaflet.readthedocs.io).

xarray-leaflet uses ipyleaflet as a plotting backend for data arrays. It generates map tiles on the fly, possibly using Dask's lazy evaluation system and Zarr's chunked data storage, and serves them through the Jupyter server, allowing for big data visualization.

See the example notebooks:
- [examples/introduction.ipynb](https://github.com/davidbrochart/xarray_leaflet/blob/master/examples/introduction.ipynb) to get started.
- [examples/dynamic.ipynb](https://github.com/davidbrochart/xarray_leaflet/blob/master/examples/dynamic.ipynb) for more advances visualizations, and in particular dynamic maps.

## How does it compare to other visualization libraries?

xarray-leaflet doesn't try to reinvent the wheel. It stands on the shoulders of giants: xarray, Dask, Jupyter widgets, Leaflet. By combining this software stack, it opens up new possibilities while being a relatively small library.

## Installation

Using conda:

```bash
conda install -c conda-forge xarray_leaflet
```

Using pip:

```bash
pip install xarray_leaflet
```

By default xarray-leaflet will generate tiles in temporary directories. For dynamic maps it will create a new directory each time you interact with the map, either by dragging or zooming. This is because there is a direct mapping between the tile directory and the URL where tiles are served. Since for dynamic maps, tiles should not be cached by the browser, the URL has to constantly change. These temporary directories are not automatically cleaned up at the moment, so you might want to do it yourself. In Unix-like systems they are under `/tmp/xarray_leaflet_*`.
