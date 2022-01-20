[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/davidbrochart/xarray_leaflet/main?urlpath=lab%2Ftree%2Fexamples%2Fintroduction.ipynb)
[![Build Status](https://github.com/davidbrochart/xarray_leaflet/workflows/UI%20Tests/badge.svg)](https://github.com/davidbrochart/xarray_leaflet/actions)

# xarray-leaflet: an xarray extension for tiled map plotting

[xarray](http://xarray.pydata.org) and [Leaflet](https://leafletjs.com) share this ability to work with fragments of data, xarray through Dask's chunks, and Leaflet through map tiles. In the end this is really the same concept, so it was a natural thing to make them work together.

Fortunately xarray is written in Python, and we happen to have a great Python binding for Leaflet, [ipyleaflet](https://ipyleaflet.readthedocs.io).

xarray-leaflet uses ipyleaflet as a plotting backend for data arrays. It generates map tiles on the fly, possibly using Dask's lazy evaluation system and Zarr's chunked data storage, and serves them through the Jupyter server, allowing for big data visualization.

See the example notebooks:
- [examples/introduction.ipynb](https://github.com/davidbrochart/xarray_leaflet/blob/master/examples/introduction.ipynb) to get started.
- [examples/dynamic.ipynb](https://github.com/davidbrochart/xarray_leaflet/blob/master/examples/dynamic.ipynb) for more advanced visualizations, and in particular dynamic maps.
- [examples/custom_projection.ipynb](https://github.com/davidbrochart/xarray_leaflet/blob/master/examples/custom_projection.ipynb) for an example of non-Mercator projection.

## How does it compare to other visualization libraries?

xarray-leaflet doesn't try to reinvent the wheel. It stands on the shoulders of giants: xarray, Jupyter widgets, Leaflet. By combining this software stack, it opens up new possibilities while being a relatively small library.

## Installation

Using conda:

```bash
conda install -c conda-forge xarray_leaflet
```

Using pip:

```bash
pip install xarray_leaflet
```

## Using xarray-leaflet with Voila

To work with xarray-leaflet, Voila has to be launched with the following command:

```bash
jupyter server --ServerApp.open_browser=True --ServerApp.default_url="voila/render/path_to_notebook.ipynb"
```
