import os
import asyncio
from tempfile import TemporaryDirectory

import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import mercantile
from ipyleaflet import LocalTileLayer, WidgetControl, DrawControl
from ipyspin import Spinner
from traitlets import HasTraits, Bool, observe
from rasterio.warp import Resampling

from .transform import passthrough, normalize, coarsen
from .utils import reproject_custom, reproject_not_custom, write_image, get_bbox_tiles, get_transform, wait_for_change


@xr.register_dataarray_accessor('leaflet')
class LeafletMap(HasTraits):
    """A xarray.DataArray extension for tiled map plotting, based on (ipy)leaflet.
    """

    map_ready = Bool(False)

    @observe('map_ready')
    def _map_ready_changed(self, change):
        self._start()


    def __init__(self, da):
        self._da = da
        self._da_selected = None

    def plot(self,
             m,
             x_dim='x',
             y_dim='y',
             fit_bounds=True,
             rgb_dim=None,
             transform0=None,
             transform1=passthrough,
             transform2=coarsen(),
             transform3=passthrough,
             colormap=None,
             persist=True,
             dynamic=False,
             tile_dir=None,
             tile_height=256,
             tile_width=256,
             resampling=Resampling.nearest,
             debug_output=None):
        """Display an array as an interactive map.

        Assumes that the pixels are given on a regular grid
        (fixed spacing in x and y).

        Parameters
        ----------
        m : ipyleaflet.Map
            The map on while to show the layer.
        y_dim : str, optional
            Name of the y dimension/coordinate
            (default: 'y').
        x_dim : str, optional
            Name of the x dimension/coordinate
            (default: 'x').
        rgb_dim : str, optional
            Name of the RGB dimension/coordinate
            (default: None).
        transform0 : function, optional
            Transformation over the whole DataArray.
        transform1 : function, optional
            Transformation over the visible DataArray.
        transform2 : function, optional
            Transformation over a tile before reprojection.
        transform3 : function, optional
            Transformation over a tile before saving to PNG.
        colormap : function, optional
            The colormap function to use for the tile PNG
            (default: matplotlib.pyplot.cm.inferno).
        persist : bool, optional
            Whether to keep the tile files (True) or not (False).
        dynamic : bool, optional
            Whether the map is dynamic (True) or not (False). If True then the
            tiles will refreshed each time the map is dragged or zoomed.
        tile_dir : str, optional
            The path to the tile directory (must be absolute).
        tile_height : int, optional
            The heiht of a tile in pixels (default: 256).
        tile_width : int, optional
            The width of a tile in pixels (default: 256).
        resampling : int, optional
            The resampling method to use, see rasterio.warp.reproject
            (default: Resampling.nearest).

        Returns
        -------
        l : ipyleaflet.LocalTileLayer
            A handler to the layer that is added to the map.
        """

        self.debug_output = debug_output

        if 'proj4def' in m.crs:
            # it's a custom projection
            if dynamic:
                raise RuntimeError('Dynamic maps are only supported for Web Mercator (EPSG:3857), not {}'.format(m.crs))
            self.dst_crs = m.crs['proj4def']
            self.web_mercator = False
            self.custom_proj = True
        elif m.crs['name'].startswith('EPSG'):
            epsg = m.crs['name'][4:]
            if dynamic and epsg != '3857':
                raise RuntimeError('Dynamic maps are only supported for Web Mercator (EPSG:3857), not {}'.format(m.crs))
            self.dst_crs = 'EPSG:' + epsg
            self.web_mercator = epsg == '3857'
            self.custom_proj = False
        else:
            raise RuntimeError('Unsupported map projection: {}'.format(m.crs))

        self.nodata = self._da.rio.nodata
        var_dims = self._da.dims
        expected_dims = [y_dim, x_dim]
        if rgb_dim is not None:
            expected_dims.append(rgb_dim)
        if set(var_dims) != set(expected_dims):
            raise ValueError(
                "Invalid dimensions in DataArray: "
                "should include only {}, found {}."
                .format(tuple(expected_dims), var_dims)
            )

        if rgb_dim is not None and colormap is not None:
            raise ValueError(
                "Cannot have a RGB dimension and a "
                "colormap at the same time."
            )
        elif rgb_dim is None:
            if colormap is None:
                colormap = plt.cm.inferno
            if transform0 is None:
                transform0 = normalize
        else:
            # there is a RGB dimension
            if transform0 is None:
                transform0 = passthrough

        self.resampling = resampling
        self.tile_dir = tile_dir
        self.persist = persist
        self.attrs = self._da.attrs
        self.m = m
        self.dynamic = dynamic
        self.tile_width = tile_width
        self.tile_height = tile_height
        self.transform0 = transform0
        self.transform1 = transform1
        self.transform2 = transform2
        self.transform3 = transform3
        self.colormap = colormap
        if self.dynamic:
            self.persist = False
            self.tile_dir = None

        self._da = self._da.rename({y_dim: 'y', x_dim: 'x'})
        if rgb_dim is None:
            self.is_rgb = False
        else:
            self.is_rgb = True
            self._da = self._da.rename({rgb_dim: 'rgb'})

        # ensure latitudes are descending
        if np.any(np.diff(self._da.y.values) >= 0):
            self._da = self._da.sel(y=slice(None, None, -1))

        # infer grid specifications (assume a rectangular grid)
        y = self._da.y.values
        x = self._da.x.values

        self.x_left = float(x.min())
        self.x_right = float(x.max())
        self.y_lower = float(y.min())
        self.y_upper = float(y.max())

        self.dx = float((self.x_right - self.x_left) / (x.size - 1))
        self.dy = float((self.y_upper - self.y_lower) / (y.size - 1))

        if fit_bounds:
            asyncio.ensure_future(self.async_fit_bounds())
        else:
            asyncio.ensure_future(self.async_wait_for_bounds())

        self.l = LocalTileLayer()
        if self._da.name is not None:
            self.l.name = self._da.name

        self._da_notransform = self._da

        self.spinner = Spinner()
        self.spinner.radius = 5
        self.spinner.length = 3
        self.spinner.width = 5
        self.spinner.lines = 8
        self.spinner.color = '#000000'
        self.spinner.layout.height = '30px'
        self.spinner.layout.width = '30px'
        self.spinner_control = WidgetControl(widget=self.spinner, position='bottomright')

        return self.l


    def select(self, draw_control=None):
        if draw_control is None:
            self._draw_control = DrawControl()
            self._draw_control.polygon = {}
            self._draw_control.polyline = {}
            self._draw_control.circlemarker = {}
            self._draw_control.rectangle = {
                'shapeOptions': {
                    'fillOpacity': 0.5
                }
            }
        else:
            self._draw_control = draw_control
        self._draw_control.on_draw(self._get_selection)
        self.m.add_control(self._draw_control)


    def unselect(self):
        self.m.remove_control(self._draw_control)


    def get_selection(self):
        return self._da_selected


    def _get_selection(self, *args, **kwargs):
        if self._draw_control.last_draw['geometry'] is not None:
            lonlat = self._draw_control.last_draw['geometry']['coordinates'][0]
            lats = [ll[1] for ll in lonlat]
            lons = [ll[0] for ll in lonlat]
            lt0, lt1 = min(lats), max(lats)
            ln0, ln1 = min(lons), max(lons)
            self._da_selected = self._da_notransform.sel(y=slice(lt1, lt0), x=slice(ln0, ln1))


    def _start(self):
        self.m.add_control(self.spinner_control)
        self._da, self.transform0_args = get_transform(self.transform0(self._da, debug_output=self.debug_output))

        self.url = self.m.window_url
        if self.url.endswith('/lab'):
            # we are in JupyterLab
            self.base_url = self.url[:-4]
        else:
            if '/notebooks/' in self.url:
                # we are in classical Notebook
                i = self.url.rfind('/notebooks/')
                self.base_url = self.url[:i]
            elif '/voila/' in self.url:
                # we are in Voila
                i = self.url.rfind('/voila/')
                self.base_url = self.url[:i]

        if self.tile_dir is None:
            self.tile_temp_dir = TemporaryDirectory(prefix='xarray_leaflet_')
            self.tile_path = self.tile_temp_dir.name
        else:
            self.tile_path = self.tile_dir
        self.url = self.base_url + '/xarray_leaflet' + self.tile_path + '/{z}/{x}/{y}.png'
        self.l.path = self.url

        self.m.remove_control(self.spinner_control)
        self._get_tiles()
        self.m.observe(self._get_tiles, names='pixel_bounds')
        if not self.dynamic:
            self.m.add_layer(self.l)


    def _get_tiles(self, change=None):
        with self.debug_output:
            self.m.add_control(self.spinner_control)
            if self.dynamic:
                self.tile_temp_dir.cleanup()
                self.tile_temp_dir = TemporaryDirectory(prefix='xarray_leaflet_')
                new_tile_path = self.tile_temp_dir.name
                new_url = self.base_url + '/xarray_leaflet' + new_tile_path + '/{z}/{x}/{y}.png'
                if self.l in self.m.layers:
                    self.m.remove_layer(self.l)

            (left, top), (right, bottom) = self.m.pixel_bounds
            (south, west), (north, east) = self.m.bounds
            z = int(self.m.zoom)  # TODO: support non-integer zoom levels?
            if self.custom_proj:
                resolution = self.m.crs['resolutions'][z]

            if self.web_mercator:
                tiles = list(mercantile.tiles(west, south, east, north, z))
            else:
                x0, x1 = int(left) // self.tile_width, int(right) // self.tile_width + 1
                y0, y1 = int(top) // self.tile_height, int(bottom) // self.tile_height + 1
                tiles = [mercantile.Tile(x, y, z) for x in range(x0, x1) for y in range(y0, y1)]

            if self.dynamic:
                # dynamic maps are redrawn at each interaction with the map
                # so we can take exactly the corresponding slice in the original data
                da_visible = self._da.sel(y=slice(north, south), x=slice(west, east))
            elif self.web_mercator:
                # for static web mercator maps we can't redraw a tile once it has been (partly) displayed,
                # so we must slice the original data on tile boundaries
                bbox = get_bbox_tiles(tiles)
                # take one more source data point to avoid glitches
                da_visible = self._da.sel(y=slice(bbox.north + self.dy, bbox.south - self.dy), x=slice(bbox.west - self.dx, bbox.east + self.dx))
            else:
                # it's a custom projection or not web mercator, the visible tiles don't translate easily
                # to a slice of the original data, so we keep everything
                # TODO: slice the data for EPSG3395, EPSG4326, Earth, Base and Simple
                da_visible = self._da

            # check if we have some data to show
            if 0 not in da_visible.shape:
                da_visible, transform1_args = get_transform(self.transform1(da_visible, *self.transform0_args, debug_output=self.debug_output))

            if self.dynamic:
                self.tile_path = new_tile_path
                self.url = new_url

            for tile in tiles:
                x, y, z = tile
                path = f'{self.tile_path}/{z}/{x}/{y}.png'
                # if static map, check if we already have the tile
                # if dynamic map, new tiles are always created
                if self.dynamic or not os.path.exists(path):
                    if self.web_mercator:
                        bbox = mercantile.bounds(tile)
                        xy_bbox = mercantile.xy_bounds(tile)
                        x_pix = (xy_bbox.right - xy_bbox.left) / self.tile_width
                        y_pix = (xy_bbox.top - xy_bbox.bottom) / self.tile_height
                        # take one more source data point to avoid glitches
                        da_tile = da_visible.sel(y=slice(bbox.north + self.dy, bbox.south - self.dy), x=slice(bbox.west - self.dx, bbox.east + self.dx))
                    else:
                        da_tile = da_visible
                    # check if we have data for this tile
                    if 0 in da_tile.shape:
                        write_image(path, None, self.persist)
                    else:
                        da_tile.attrs = self.attrs
                        da_tile, transform2_args = get_transform(self.transform2(da_tile, tile_width=self.tile_width, tile_height=self.tile_height, debug_output=self.debug_output), *transform1_args)
                        # reproject each RGB component if needed
                        # TODO: must be doable with xarray.apply_ufunc
                        if self.is_rgb:
                            das = [da_tile.isel(rgb=i) for i in range(3)]
                        else:
                            das = [da_tile]
                        for i in range(len(das)):
                            das[i] = das[i].rio.write_nodata(self.nodata)
                            if self.custom_proj:
                                das[i] = reproject_custom(das[i], self.dst_crs, x, y, z, resolution, resolution, self.tile_width, self.tile_height, self.resampling)
                            else:
                                das[i] = reproject_not_custom(das[i], self.dst_crs, xy_bbox.left, xy_bbox.top, x_pix, y_pix, self.tile_width, self.tile_height, self.resampling)
                            das[i], transform3_args = get_transform(self.transform3(das[i], *transform2_args, debug_output=self.debug_output))
                        if self.is_rgb:
                            alpha = np.where(das[0]==self._da.rio.nodata, 0, 255)
                            das.append(alpha)
                            da_tile = np.stack(das, axis=2)
                            write_image(path, da_tile, self.persist)
                        else:
                            da_tile = self.colormap(das[0])
                            write_image(path, da_tile*255, self.persist)

            if self.dynamic:
                self.l.path = self.url
                self.m.add_layer(self.l)
                self.l.redraw()

            self.m.remove_control(self.spinner_control)


    async def async_wait_for_bounds(self):
        if len(self.m.bounds) == 0:
            await wait_for_change(self.m, 'bounds')
        self.map_ready = True


    async def async_fit_bounds(self):
        center = self.y_lower + (self.y_upper - self.y_lower) / 2, self.x_left + (self.x_right - self.x_left) / 2
        if center != self.m.center:
            self.m.center = center
            await wait_for_change(self.m, 'bounds')
        zoomed_out = False
        # zoom out
        while True:
            if self.m.zoom <= 1:
                break
            if self.debug_output is not None:
                with self.debug_output:
                    print('Zooming out')
            (south, west), (north, east) = self.m.bounds
            if south > self.y_lower or north < self.y_upper or west > self.x_left or east < self.x_right:
                self.m.zoom = self.m.zoom - 1
                await wait_for_change(self.m, 'bounds')
                zoomed_out = True
            else:
                break
        if not zoomed_out:
            # zoom in
            while True:
                if self.debug_output is not None:
                    with self.debug_output:
                        print('Zooming in')
                (south, west), (north, east) = self.m.bounds
                if south < self.y_lower and north > self.y_upper and west < self.x_left and east > self.x_right:
                    self.m.zoom = self.m.zoom + 1
                    await wait_for_change(self.m, 'bounds')
                else:
                    self.m.zoom = self.m.zoom - 1
                    await wait_for_change(self.m, 'bounds')
                    break
        self.map_ready = True
