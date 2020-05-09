# xarray-leaflet: an xarray extension for tiled map plotting

[xarray](http://xarray.pydata.org) and [Leaflet](https://leafletjs.com) share this ability to work with fragments of data, xarray through Dask's chunks, and Leaflet through map tiles. In the end this is really the same concept, so it was a natural thing to make them work together.

Fortunately xarray is written in Python, and we happen to have a great Python binding for Leaflet, [ipyleaflet](https://ipyleaflet.readthedocs.io).

`xarray-leaflet` uses `ipyleaflet` as a plotting backend for data arrays. It generates map tiles on the fly, possibly using Dask's lazy evaluation system and Zarr's chunked data storage, and serves them through the Jupyter server, allowing for big data visualization.

See [examples/introduction.ipynb](https://github.com/davidbrochart/xarray_leaflet/blob/master/examples/introduction.ipynb) for an example notebook.

## Installation

```bash
pip install xarray_leaflet
```

`xarray-leaflet` currently only works in the classic Jupyter Notebook (not JupyterLab).
