"""Main module."""
import os
from tempfile import mkdtemp

import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import mercantile
from ipyleaflet import LocalTileLayer
from traitlets import observe

from .transform import passthrough, normalize, coarsen
from .utils import get_webmercator, write_image, reindex, get_bbox_tiles, get_transform


@xr.register_dataarray_accessor('leaflet')
class LeafletMap:
    """A xarray.DataArray extension for tiled map plotting, based on (ipy)leaflet.

    """
    def __init__(self, da):
        self._da = da

    def plot(self, m, y_dim='latitude', x_dim='longitude',
             transform0=normalize,
             transform1=passthrough,
             transform2=coarsen(),
             transform3=passthrough,
             colormap=plt.cm.inferno,
             persist=True,
             dynamic=False,
             tile_dir=None):
        """Display an array as an interactive map.

        Assumes that the pixels are given on a regular grid
        (fixed spacing in x and y).

        Parameters
        ----------
        m : ipyleaflet.Map
            The map on while to show the layer
        y_dim : str, optional
            Name of the y dimension/coordinate
            (default: 'latitude').
        x_dim : str, optional
            Name of the x dimension/coordinate
            (default: 'longitude').
        transform0 : function, optional
            Transformation over the visible DataArray
        transform1 : function, optional
            Transformation over the tiles before reprojection
        transform2 : function, optional
            Transformation over the tiles before saving to PNG
        persist : bool, optional
            Whether to keep the tile files (True) or not (False)
        dynamic : bool, optional
            Whether the map is dynamic (True) or not (False). If True then the
            tiles will refreshed each time the map is dragged or zoomed.
        tile_dir : str, optional
            The path to the tile directory (must be absolute)

        Returns
        -------
        l : ipyleaflet.LocalTileLayer
            A handler to the layer that is added to the map
        """

        da = self._da
        var_dims = da.dims
        if set(var_dims) != set([y_dim, x_dim]):
            raise ValueError(
                "Invalid dimensions in DataArray: "
                "should include only {}, found {}."
                .format((y_dim, x_dim), var_dims)
            )

        da = da.rename({y_dim: 'y', x_dim: 'x'})

        # ensure latitudes are descending
        if np.any(np.diff(da.y.values) >= 0):
            da = da.sel(y=slice(None, None, -1))

        da, transform0_args = get_transform(transform0(da))

        # infer grid specifications (assume a rectangular grid)
        y = da.y.values
        x = da.x.values

        x_left = float(x.min())
        x_right = float(x.max())
        y_lower = float(y.min())
        y_upper = float(y.max())

        dx = float((x_right - x_left) / (x.size - 1))
        dy = float((y_upper - y_lower) / (y.size - 1))

        map_started = False
        l = LocalTileLayer()

        def main(change):
            nonlocal da, map_started, tile_dir, persist
            if not map_started and len(m.bounds) > 0:
                map_started = True

                url = m.window_url
                if url.endswith('/lab'):
                    # we are in JupyterLab
                    base_url = url[:-4]
                else:
                    if '/notebooks/' in url:
                        # we are in classical Notebook
                        i = url.rfind('/notebooks/')
                        base_url = url[:i]
                    else:
                        # we are in Voila
                        # TODO: make it work when Voila uses Jupyter server's ExtensionApp
                        base_url = url.rstrip('/')

                if dynamic:
                    persist = True
                    tile_dir = None

                tile_path = tile_dir or mkdtemp(prefix='xarray_leaflet_')
                url = base_url + '/xarray_leaflet' + tile_path + '/{z}/{x}/{y}.png'
                l.path = url

                def get_tiles(change):
                    nonlocal url, tile_path
                    if dynamic:
                        new_tile_path = mkdtemp(prefix='xarray_leaflet_')
                        new_url = base_url + '/xarray_leaflet' + new_tile_path + '/{z}/{x}/{y}.png'
                        if l in m.layers:
                            m.remove_layer(l)
                    ((south, west), (north, east)) = change['new']
                    tiles = list(mercantile.tiles(west, south, east, north, m.zoom))
                    if dynamic:
                        da_visible = da.sel(y=slice(north, south), x=slice(west, east))
                    else:
                        bbox = get_bbox_tiles(tiles)
                        da_visible = da.sel(y=slice(bbox.north, bbox.south), x=slice(bbox.west, bbox.east))
                    # check if we have visible data
                    if 0 not in da_visible.shape:
                        da_visible, transform1_args = get_transform(transform1(da_visible, *transform0_args))
                    if dynamic:
                        tile_path = new_tile_path
                        url = new_url
                    for tile in tiles:
                        path = f'{tile_path}/{int(tile.z)}/{tile.x}/{tile.y}.png'
                        # if static map, check if we already have the tile
                        # if dynamic map, new tiles are always created
                        if dynamic or not os.path.exists(path):
                            bbox = mercantile.bounds(tile)
                            da_tile = da_visible.sel(y=slice(bbox.north, bbox.south), x=slice(bbox.west, bbox.east))
                            # check if we have data for this tile
                            if 0 not in da_tile.shape:
                                da_tile = reindex(da_tile, bbox, dx, dy)
                                da_tile, transform2_args = get_transform(transform2(da_tile, *transform1_args))
                                np_tile = get_webmercator(da_tile.values, bbox.west, bbox.north, dx, dy)
                                np_tile, transform3_args = get_transform(transform3(np_tile, *transform2_args))
                                np_tile = colormap(np_tile)
                                write_image(np_tile, path, persist)
                    if dynamic:
                        l.path = url
                        m.add_layer(l)
                        l.redraw()

                get_tiles({'new': m.bounds})
                m.observe(get_tiles, names='bounds')
                if not dynamic:
                    m.add_layer(l)

        main(None)
        m.observe(main, names='bounds')
        return l
