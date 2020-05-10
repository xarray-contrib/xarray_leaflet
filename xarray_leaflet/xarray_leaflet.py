"""Main module."""
import os
from tempfile import mkdtemp

import xarray as xr
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from PIL import Image
import numpy as np
from affine import Affine
import mercantile
from ipyleaflet import LocalTileLayer
from traitlets import observe


def no_transform(array, *args):
    return array


@xr.register_dataarray_accessor('leaflet')
class LeafletMap:
    """A xarray.DataArray extension for tiled map plotting, based on (ipy)leaflet.

    """
    def __init__(self, da):
        self._da = da

    def plot(self, m, lat_dim='latitude', lon_dim='longitude',
             transform0=no_transform,
             transform1=no_transform,
             transform2=no_transform,
             persist=True,
             dynamic=False,
             tile_dir=None):
        """Display a lat/lon array as an interactive map.

        Assumes that the pixels are given on a regular grid
        (fixed spacing in latitude and longitude).

        Parameters
        ----------
        m : ipyleaflet.Map
            The map on while to show the layer
        lat_dim : str, optional
            Name of the latitude dimension/coordinate
            (default: 'latitude').
        lon_dim : str, optional
            Name of the longitude dimension/coordinate
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
        if dynamic:
            persist = True
            tile_dir = None

        tile_root_dir = './xarray_leaflet_tiles'
        os.makedirs(tile_root_dir, exist_ok=True)
        tile_path = tile_dir or mkdtemp(dir=tile_root_dir)
        url = tile_path + '/{z}/{x}/{y}.png'
        l = LocalTileLayer(path=url)

        da = self._da.copy()
        var_dims = da.dims

        if set(var_dims) != set([lat_dim, lon_dim]):
            raise ValueError(
                "Invalid dimensions in DataArray: "
                "should include only {}, found {}."
                .format((lat_dim, lon_dim), var_dims)
            )

        # ensure latitudes are descending
        if np.any(np.diff(da[lat_dim].values) >= 0):
            da = da.sel(**{lat_dim: slice(None, None, -1)})

        # infer grid specifications (assume a rectangular grid)
        lat = da[lat_dim].values
        lon = da[lon_dim].values

        lon_left = float(lon.min())
        lon_right = float(lon.max())
        lat_lower = float(lat.min())
        lat_upper = float(lat.max())

        dx = float((lon_right - lon_left) / (lon.size - 1))
        dy = float((lat_upper - lat_lower) / (lat.size - 1))

        def get_tiles(change):
            nonlocal url, tile_path
            if dynamic:
                new_tile_path = mkdtemp(dir=tile_root_dir)
                new_url = new_tile_path + '/{z}/{x}/{y}.png'
                if l in m.layers:
                    m.remove_layer(l)
            ((south, west), (north, east)) = change['new']
            tiles = list(mercantile.tiles(west, south, east, north, m.zoom))
            if dynamic:
                da_visible = da.sel(**{lat_dim: slice(north, south), lon_dim: slice(west, east)})
            else:
                bbox = get_bbox_tiles(tiles)
                da_visible = da.sel(**{lat_dim: slice(bbox.north, bbox.south), lon_dim: slice(bbox.west, bbox.east)})
            # check if we have visible data
            if 0 not in da_visible.shape:
                da_visible, transform1_args = get_transform(transform0(da_visible))
            if dynamic:
                for tile in tiles:
                    path = f'{tile_path}/{int(tile.z)}/{tile.x}/{tile.y}.png'
                    empty_image(path)
                    with open(path[:-4] + '.done', 'wt') as f:
                        f.write('delete')
                tile_path = new_tile_path
                url = new_url
            for tile in tiles:
                path = f'{tile_path}/{int(tile.z)}/{tile.x}/{tile.y}.png'
                # check if we already have the tile
                if dynamic or not os.path.exists(path):
                    bbox = mercantile.bounds(tile)
                    da_tile = da_visible.sel(**{lat_dim: slice(bbox.north, bbox.south), lon_dim: slice(bbox.west, bbox.east)})
                    # check if we have data for this tile
                    if 0 in da_tile.shape:
                        empty_image(path)
                    else:
                        da_tile = reindex(da_tile, lon_dim, lat_dim, bbox, dx, dy)
                        da_tile, transform2_args = get_transform(transform1(da_tile, *transform1_args))
                        np_tile = get_webmercator(da_tile.values, bbox.west, bbox.north, dx, dy)
                        np_tile, transform3_args = get_transform(transform2(np_tile, *transform2_args))
                        write_image(np_tile, path)
                    with open(path[:-4] + '.done', 'wt') as f:
                        if persist:
                            f.write('keep')
                        else:
                            f.write('delete')
            if dynamic:
                l.path = url
                m.add_layer(l)
                l.redraw()

        get_tiles({'new': m.bounds})
        m.observe(get_tiles, names='bounds')
        if not dynamic:
            m.add_layer(l)
        return l


def get_webmercator(source, west, north, dx, dy):
    with rasterio.Env():
        rows, cols = source.shape
        east = west + dx * cols
        south = north - dy * rows
        src_transform = Affine(dx, 0, west, 0, -dy, north)
        src_crs = {'init': 'EPSG:4326'}

        dst_crs = {'init': 'EPSG:3857'}
        width = height = 256
        dst_transform, new_width, new_height = calculate_default_transform(src_crs, dst_crs, cols, rows, bottom=south, left=east, top=north, right=west, dst_width=width, dst_height=height)
        assert width == new_width
        assert height == new_height

        destination = np.zeros((height, width))

        reproject(
            source,
            destination,
            src_transform=src_transform,
            src_crs=src_crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest)

    return destination


def write_image(data_array, path):
    im = Image.fromarray(np.uint8(data_array*255))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    im.save(path)


def empty_image(path):
    im = Image.new('RGB', (256, 256))
    im.putalpha(256)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    im.save(path)


def reindex(da_tile, lon_dim, lat_dim, bbox, dx, dy):
    do_reindex = False
    lons = da_tile[lon_dim].values
    # check if we have data on the left
    if lons[0] - dx < bbox.west:
        # check if we have data on the right
        if not (lons[-1] + dx > bbox.east):
            # we need to pad on the right
            do_reindex = True
            lons = np.arange(lons[0], bbox.east, dx)
    else:
        # we need to pad on the left
        do_reindex = True
        lons = np.arange(lons[-1], bbox.west, -dx)[::-1]
    lats = da_tile[lat_dim].values
    # check if we have data at the top
    if lats[0] + dy > bbox.north:
        # check if we have data at the bottom
        if not (lats[-1] - dy < bbox.south):
            # we need to pad at the bottom
            do_reindex = True
            lats = np.arange(lats[0], bbox.south, -dy)
    else:
        # we need to pad at the top
        do_reindex = True
        lats = np.arange(lats[-1], bbox.north, dy)[::-1]
    if do_reindex:
        da_tile = da_tile.reindex(**{lat_dim: lats, lon_dim: lons}, method='nearest', tolerance=dx/4)
    return da_tile


def get_bbox_tiles(tiles):
    north = east = -float('inf')
    south = west = float('inf')
    for tile in tiles:
        bbox = mercantile.bounds(tile)
        north = max(north, bbox.north)
        south = min(south, bbox.south)
        west = min(west, bbox.west)
        east = max(east, bbox.east)
    bbox_tiles = mercantile.LngLatBbox(west, south, east, north)
    return bbox_tiles


def get_transform(result):
    if type(result) == tuple:
        array, *args = result
    else:
        array = result
        args = []
    return array, args
