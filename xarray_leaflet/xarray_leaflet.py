import asyncio
import os
import tempfile
from typing import Callable, Optional

import geopandas as gpd
import matplotlib as mpl
import mercantile
import numpy as np
import pandas as pd
import xarray as xr
from ipyleaflet import DrawControl, LocalTileLayer, WidgetControl
from ipyspin import Spinner
from IPython.display import Image, display
from ipyurl import Url
from ipywidgets import Output
from matplotlib import pyplot as plt
from rasterio.warp import Resampling
from traitlets import Bool, HasTraits, observe

from .transform import coarsen, normalize, passthrough
from .utils import debug  # noqa
from .utils import (
    get_bbox_tiles,
    get_transform,
    reproject_custom,
    reproject_not_custom,
    wait_for_change,
    write_image,
)
from .vector import Zvect


class Leaflet(HasTraits):

    is_vector: bool

    map_ready = Bool(False)

    @observe("map_ready")
    def _map_ready_changed(self, change):
        self._start()

    def plot(
        self,
        m,
        *,
        # raster or vector options:
        get_base_url: Optional[Callable] = None,
        dynamic: Optional[bool] = None,
        persist: bool = True,
        tile_dir=None,
        tile_height: int = 256,
        tile_width: int = 256,
        layer: Callable = LocalTileLayer,
        transform3=passthrough,
        colormap=None,
        # raster-only options:
        x_dim="x",
        y_dim="y",
        fit_bounds=True,
        rgb_dim=None,
        transform0=None,
        transform1=passthrough,
        transform2=coarsen(),
        colorbar_position="topright",
        resampling=Resampling.nearest,
        # vector-only options:
        measurement: Optional[str] = None,
        visible_callback: Optional[Callable] = None,
        rasterize_function: Optional[Callable] = None,
    ):
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
        fit_bounds: bool, optional
            Set the map to fit the bounds of the array (default True).
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
            (default: matplotlib.pyplot.cm.viridis).
        colorbar_position : str, optional
            Where to show the colorbar (default: "topright").
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
        get_base_url: callable, optional
            A function taking the window URL and returning the base URL to use.
        measurement: str, optional
            The geocube measurement.
        visible_callback: callable, optional
            A callable taking the following arguments:
            - the ipyleaflet.Map
            - the xarray.DataArray of the visible region
            - the mercantile.LngLatBbox of the visible region

            and returning True if the layer should be shown, False otherwise.
        rasterize_function: callable, optional
            A callable passed to make_geocube. Defaults to:
            partial(rasterize_image, merge_alg=MergeAlg.add, all_touched=True)

        Returns
        -------
        layer : ipyleaflet.LocalTileLayer
            A handler to the layer that is added to the map.
        """

        self.layer = layer()

        self.transform0 = transform0
        self.transform1 = transform1
        self.transform2 = transform2
        self.transform3 = transform3

        if self.is_vector:
            # source is a GeoDataFrame (vector)
            self.visible_callback = visible_callback
            if measurement is None:
                raise RuntimeError("You must provide a 'measurement'.")
            if dynamic is None:
                dynamic = True
            if not dynamic:
                self.vmin = self._df[measurement].min()
                self.vmax = self._df[measurement].max()
            self.measurement = measurement
            zarr_temp_dir = tempfile.TemporaryDirectory(prefix="xarray_leaflet_zarr_")
            self.zvect = Zvect(
                self._df,
                measurement,
                rasterize_function,
                tile_width,
                tile_height,
                zarr_temp_dir.name,
            )
            if colormap is None:
                colormap = plt.cm.viridis
        else:
            # source is a DataArray (raster)
            if dynamic is None:
                dynamic = False
            if "proj4def" in m.crs:
                # it's a custom projection
                if dynamic:
                    raise RuntimeError(
                        f"Dynamic maps are only supported for Web Mercator (EPSG:3857), not {m.crs}"
                    )
                self.dst_crs = m.crs["proj4def"]
                self.web_mercator = False
                self.custom_proj = True
            elif m.crs["name"].startswith("EPSG"):
                epsg = m.crs["name"][4:]
                if dynamic and epsg != "3857":
                    raise RuntimeError(
                        f"Dynamic maps are only supported for Web Mercator (EPSG:3857), not {m.crs}"
                    )
                self.dst_crs = "EPSG:" + epsg
                self.web_mercator = epsg == "3857"
                self.custom_proj = False
            else:
                raise RuntimeError(f"Unsupported map projection: {m.crs}")

            self.nodata = self._da.rio.nodata
            var_dims = self._da.dims
            expected_dims = [y_dim, x_dim]
            if rgb_dim is not None:
                expected_dims.append(rgb_dim)
            if set(var_dims) != set(expected_dims):
                raise ValueError(
                    "Invalid dimensions in DataArray: "
                    "should include only {}, found {}.".format(
                        tuple(expected_dims), var_dims
                    )
                )

            if rgb_dim is not None and colormap is not None:
                raise ValueError(
                    "Cannot have a RGB dimension and a " "colormap at the same time."
                )
            elif rgb_dim is None:
                if colormap is None:
                    colormap = plt.cm.viridis
                if transform0 is None:
                    transform0 = normalize
            else:
                # there is a RGB dimension
                if transform0 is None:
                    transform0 = passthrough

            self.attrs = self._da.attrs
            self._da = self._da.rename({y_dim: "y", x_dim: "x"})
            if rgb_dim is None:
                self.is_rgb = False
            else:
                self.is_rgb = True
                self._da = self._da.rename({rgb_dim: "rgb"})

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

            if self._da.name is not None:
                self.layer.name = self._da.name

            self._da_notransform = self._da

        self.resampling = resampling
        self.dynamic = dynamic
        self.tile_dir = tile_dir
        self.persist = persist
        self.m = m
        self.tile_width = tile_width
        self.tile_height = tile_height
        self.colormap = colormap
        self.colorbar = None
        self.colorbar_position = colorbar_position

        if dynamic:
            self.persist = False
            self.tile_dir = None

        if get_base_url is None:
            self.base_url = None
            self.url_widget = Url()
            display(self.url_widget)
        else:
            self.base_url = get_base_url(self.m.window_url)

        if fit_bounds and not self.is_vector:
            asyncio.ensure_future(self.async_fit_bounds())
        else:
            asyncio.ensure_future(self.async_wait_for_bounds())

        self.spinner = Spinner()
        self.spinner.radius = 5
        self.spinner.length = 3
        self.spinner.width = 5
        self.spinner.lines = 8
        self.spinner.color = "#000000"
        self.spinner.layout.height = "30px"
        self.spinner.layout.width = "30px"
        self.spinner_control = WidgetControl(
            widget=self.spinner, position="bottomright"
        )

        return self.layer

    def select(self, draw_control=None):
        if draw_control is None:
            self._draw_control = DrawControl()
            self._draw_control.polygon = {}
            self._draw_control.polyline = {}
            self._draw_control.circlemarker = {}
            self._draw_control.rectangle = {"shapeOptions": {"fillOpacity": 0.5}}
        else:
            self._draw_control = draw_control
        self._draw_control.on_draw(self._get_selection)
        self.m.add_control(self._draw_control)

    def unselect(self):
        self.m.remove_control(self._draw_control)

    def get_selection(self):
        return self._da_selected

    def _get_selection(self, *args, **kwargs):
        if self._draw_control.last_draw["geometry"] is not None:
            lonlat = self._draw_control.last_draw["geometry"]["coordinates"][0]
            lats = [ll[1] for ll in lonlat]
            lons = [ll[0] for ll in lonlat]
            lt0, lt1 = min(lats), max(lats)
            ln0, ln1 = min(lons), max(lons)
            self._da_selected = self._da_notransform.sel(
                y=slice(lt1, lt0), x=slice(ln0, ln1)
            )

    def _start(self):
        self.m.add_control(self.spinner_control)
        if not self.is_vector:
            self._da, self.transform0_args = get_transform(self.transform0(self._da))
        else:
            self.layer.name = self.measurement

        if self.tile_dir is None:
            self.tile_temp_dir = tempfile.TemporaryDirectory(prefix="xarray_leaflet_")
            self.tile_path = self.tile_temp_dir.name
        else:
            self.tile_path = self.tile_dir
        self.url = (
            self.base_url + "/xarray_leaflet/" + self.tile_path + "/{z}/{x}/{y}.png"
        )
        self.layer.path = self.url

        self.m.remove_control(self.spinner_control)
        if not self.is_vector:
            get_tiles = self._get_raster_tiles
        elif self._df is not None:
            get_tiles = self._get_vector_tiles
        else:
            raise RuntimeError("No DataArray or GeoDataFrame provided.")
        get_tiles()
        self.m.observe(get_tiles, names="pixel_bounds")
        if not self.dynamic:
            if not self.is_vector:
                self._show_colorbar(self._da_notransform)
        self.m.add_layer(self.layer)

    def _show_colorbar(self, da):
        if self.colorbar_position and self.colormap is not None:
            vmin = da.min().values
            vmax = da.max().values
            fig = plt.figure(figsize=(8, 3))
            ax = fig.add_axes([0.05, 0.8, 0.5, 0.07])
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            cbar = mpl.colorbar.ColorbarBase(  # noqa
                ax, cmap=self.colormap, norm=norm, orientation="horizontal"
            )
            f = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
            output = Output()
            try:
                plt.savefig(f.name, bbox_inches="tight")
                with output:
                    display(Image(filename=f.name))
            finally:
                os.unlink(f.name)
                f.close()
            self.colorbar = WidgetControl(
                widget=output, position=self.colorbar_position, transparent_bg=True
            )
            self.m.add_control(self.colorbar)
            plt.close()

    def _get_vector_tiles(self, change=None):
        self.m.add_control(self.spinner_control)

        # visible bounds
        (south, west), (north, east) = self.m.bounds
        z = int(self.m.zoom)
        tiles = mercantile.tiles(west, south, east, north, z)

        if self.dynamic:
            # get DataArray for the visible map
            llbbox = mercantile.LngLatBbox(west, south, east, north)
            da_visible = self.zvect.get_da_llbbox(llbbox, z)
            # check if we must show the layer
            if self.visible_callback and not self.visible_callback(
                self.m, da_visible, llbbox
            ):
                self.m.remove_control(self.spinner_control)
                return
            if da_visible is None:
                vmin = vmax = 0
            else:
                vmin = da_visible.min()
                vmax = da_visible.max()
        else:
            vmin = self.vmin
            vmax = self.vmax
            da_visible_computed = False

        for tile in tiles:
            x, y, z = tile
            path = f"{self.tile_path}/{z}/{x}/{y}.png"
            if self.dynamic or not os.path.exists(path):
                if not self.dynamic and not da_visible_computed:
                    # get DataArray for the visible map
                    llbbox = mercantile.LngLatBbox(west, south, east, north)
                    da_visible = self.zvect.get_da_llbbox(llbbox, z)
                    da_visible_computed = True
                if da_visible is None:
                    da_tile = None
                else:
                    da_tile = self.zvect.get_da(z)
                    if da_tile is not None:
                        da_tile = da_tile[
                            y * self.tile_height : (y + 1) * self.tile_height,  # noqa
                            x * self.tile_width : (x + 1) * self.tile_width,  # noqa
                        ]
                        if np.isnan(da_tile).all():
                            da_tile = None
                if da_tile is None:
                    write_image(path, None)
                else:
                    da_tile = self.transform3(da_tile)
                    # normalize
                    da_tile = (da_tile - vmin) / (vmax - vmin)
                    da_tile = self.colormap(da_tile)
                    write_image(path, da_tile * 255)

        if self.dynamic:
            self.layer.redraw()

        self.m.remove_control(self.spinner_control)

    def _get_raster_tiles(self, change=None):
        self.m.add_control(self.spinner_control)

        (left, top), (right, bottom) = self.m.pixel_bounds
        (south, west), (north, east) = self.m.bounds
        z = int(self.m.zoom)  # TODO: support non-integer zoom levels?
        if self.custom_proj:
            resolution = self.m.crs["resolutions"][z]

        if self.web_mercator:
            tiles = list(mercantile.tiles(west, south, east, north, z))
        else:
            x0, x1 = int(left) // self.tile_width, int(right) // self.tile_width + 1
            y0, y1 = int(top) // self.tile_height, int(bottom) // self.tile_height + 1
            tiles = [
                mercantile.Tile(x, y, z) for x in range(x0, x1) for y in range(y0, y1)
            ]

        if self.dynamic:
            # dynamic maps are redrawn at each interaction with the map
            # so we can take exactly the corresponding slice in the original data
            da_visible = self._da.sel(y=slice(north, south), x=slice(west, east))
        elif self.web_mercator:
            # for static web mercator maps we can't redraw a tile once it has been displayed
            # (even partially) so we must slice the original data on tile boundaries
            bbox = get_bbox_tiles(tiles)
            # take one more source data point to avoid glitches
            da_visible = self._da.sel(
                y=slice(bbox.north + self.dy, bbox.south - self.dy),
                x=slice(bbox.west - self.dx, bbox.east + self.dx),
            )
        else:
            # it's a custom projection or not web mercator, the visible tiles don't translate easily
            # to a slice of the original data, so we keep everything
            # TODO: slice the data for EPSG3395, EPSG4326, Earth, Base and Simple
            da_visible = self._da

        # check if we have some data to show
        if 0 not in da_visible.shape:
            da_visible, transform1_args = get_transform(
                self.transform1(da_visible, *self.transform0_args)
            )

        for tile in tiles:
            x, y, z = tile
            path = f"{self.tile_path}/{z}/{x}/{y}.png"
            # if static map, check if we already have the tile
            # if dynamic map, new tiles are always created
            if self.dynamic or not os.path.exists(path):
                if self.web_mercator:
                    bbox = mercantile.bounds(tile)
                    xy_bbox = mercantile.xy_bounds(tile)
                    x_pix = (xy_bbox.right - xy_bbox.left) / self.tile_width
                    y_pix = (xy_bbox.top - xy_bbox.bottom) / self.tile_height
                    # take one more source data point to avoid glitches
                    da_tile = da_visible.sel(
                        y=slice(bbox.north + self.dy, bbox.south - self.dy),
                        x=slice(bbox.west - self.dx, bbox.east + self.dx),
                    )
                else:
                    da_tile = da_visible
                # check if we have data for this tile
                if 0 in da_tile.shape:
                    write_image(path, None)
                else:
                    da_tile.attrs = self.attrs
                    da_tile, transform2_args = get_transform(
                        self.transform2(
                            da_tile,
                            tile_width=self.tile_width,
                            tile_height=self.tile_height,
                        ),
                        *transform1_args,
                    )
                    # reproject each RGB component if needed
                    # TODO: must be doable with xarray.apply_ufunc
                    if self.is_rgb:
                        das = [da_tile.isel(rgb=i) for i in range(3)]
                    else:
                        das = [da_tile]
                    for i in range(len(das)):
                        das[i] = das[i].rio.write_nodata(self.nodata)
                        if self.custom_proj:
                            das[i] = reproject_custom(
                                das[i],
                                self.dst_crs,
                                x,
                                y,
                                z,
                                resolution,
                                resolution,
                                self.tile_width,
                                self.tile_height,
                                self.resampling,
                            )
                        else:
                            das[i] = reproject_not_custom(
                                das[i],
                                self.dst_crs,
                                xy_bbox.left,
                                xy_bbox.top,
                                x_pix,
                                y_pix,
                                self.tile_width,
                                self.tile_height,
                                self.resampling,
                            )
                        das[i], transform3_args = get_transform(
                            self.transform3(das[i], *transform2_args)
                        )
                    if self.is_rgb:
                        alpha = np.where(das[0] == self._da.rio.nodata, 0, 255)
                        das.append(alpha)
                        da_tile = np.stack(das, axis=2)
                        write_image(path, da_tile)
                    else:
                        da_tile = self.colormap(das[0])
                        write_image(path, da_tile * 255)

        if self.dynamic:
            if self.colorbar in self.m.controls:
                self.m.remove_control(self.colorbar)
            self._show_colorbar(
                self._da_notransform.sel(y=slice(north, south), x=slice(west, east))
            )
            self.layer.redraw()

        self.m.remove_control(self.spinner_control)

    async def async_wait_for_bounds(self):
        if len(self.m.bounds) == 0:
            await wait_for_change(self.m, "bounds")
        if self.base_url is None:
            self.base_url = (await self.url_widget.get_url()).rstrip("/")
        self.map_ready = True

    async def async_fit_bounds(self):
        center = (
            self.y_lower + (self.y_upper - self.y_lower) / 2,
            self.x_left + (self.x_right - self.x_left) / 2,
        )
        if center != self.m.center:
            self.m.center = center
            await wait_for_change(self.m, "bounds")
        zoomed_out = False
        # zoom out
        while True:
            if self.m.zoom <= 1:
                break
            (south, west), (north, east) = self.m.bounds
            if (
                south > self.y_lower
                or north < self.y_upper
                or west > self.x_left
                or east < self.x_right
            ):
                self.m.zoom = self.m.zoom - 1
                await wait_for_change(self.m, "bounds")
                zoomed_out = True
            else:
                break
        if not zoomed_out:
            # zoom in
            while True:
                (south, west), (north, east) = self.m.bounds
                if (
                    south < self.y_lower
                    and north > self.y_upper
                    and west < self.x_left
                    and east > self.x_right
                ):
                    self.m.zoom = self.m.zoom + 1
                    await wait_for_change(self.m, "bounds")
                else:
                    self.m.zoom = self.m.zoom - 1
                    await wait_for_change(self.m, "bounds")
                    break
        if self.base_url is None:
            self.base_url = (await self.url_widget.get_url()).rstrip("/")
        self.map_ready = True


@xr.register_dataarray_accessor("leaflet")
class DataArrayLeaflet(Leaflet):
    """A DataArraye extension for tiled map plotting, based on (ipy)leaflet."""

    def __init__(self, da: xr.DataArray):
        self._da = da
        self._da_selected = None
        self.is_vector = False


@pd.api.extensions.register_dataframe_accessor("leaflet")
class GeoDataFrameLeaflet(Leaflet):
    """A GeoDataFrame extension for tiled map plotting, based on (ipy)leaflet."""

    def __init__(self, df: gpd.GeoDataFrame):
        self._df = df
        self.is_vector = True
