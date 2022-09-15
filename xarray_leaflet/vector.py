import json
import math
from functools import partial
from pathlib import Path
from typing import Callable, Optional

import mercantile
import numpy as np
import xarray as xr
import zarr
from geocube.api.core import make_geocube
from geocube.rasterize import rasterize_image
from geopandas import GeoDataFrame
from rasterio.enums import MergeAlg
from shapely.geometry import box, mapping

from .utils import debug  # noqa


class Zvect:
    def __init__(
        self,
        df: GeoDataFrame,
        measurement: str,
        rasterize_function: Optional[Callable],
        width: int,
        height: int,
        root_path: str = "",
    ):
        # reproject to Web Mercator
        self.df = df.to_crs(epsg=3857)
        self.measurement = measurement
        self.rasterize_function = rasterize_function or partial(
            rasterize_image, merge_alg=MergeAlg.add, all_touched=True
        )
        self.width = width
        self.height = height
        self.zzarr = Zzarr(root_path, width, height)
        self.tiles = []

    def get_da_tile(self, tile: mercantile.Tile) -> Optional[xr.DataArray]:
        xy_bounds = mercantile.xy_bounds(tile)
        dx = (xy_bounds.right - xy_bounds.left) / self.width
        dy = (xy_bounds.top - xy_bounds.bottom) / self.height
        # take one more pixel to avoid glitches
        bbox = box(
            xy_bounds.left - dx,
            xy_bounds.bottom - dy,
            xy_bounds.right + dx,
            xy_bounds.top + dy,
        )
        geom = json.dumps(
            {**mapping(bbox), "crs": {"properties": {"name": "EPSG:3857"}}}
        )
        df_tile = self.df.clip(bbox.bounds)
        if df_tile.empty:
            return None
        ds_tile = make_geocube(
            vector_data=df_tile,
            resolution=(-dy, dx),
            measurements=[self.measurement],
            rasterize_function=self.rasterize_function,
            fill=0,
            geom=geom,
        )
        # remove added pixels
        da_tile = ds_tile[self.measurement][1:-1, 1:-1]
        return da_tile

    def get_da_llbbox(
        self, bbox: mercantile.LngLatBbox, z: int
    ) -> Optional[xr.DataArray]:
        tiles = mercantile.tiles(*bbox, z)
        all_none = True
        for tile in tiles:
            if tile in self.tiles:
                all_none = False
            else:
                da_tile = self.get_da_tile(tile)
                if da_tile is not None:
                    all_none = False
                    self.zzarr.write_to_zarr(tile, da_tile.values)
                self.tiles.append(tile)
        if all_none:
            return None
        da = self.get_da(z)
        y0, x0 = deg2idx(bbox.north, bbox.west, z, self.height, self.width, math.floor)
        y1, x1 = deg2idx(bbox.south, bbox.east, z, self.height, self.width, math.ceil)
        return da[y0:y1, x0:x1]

    def get_da(self, z: int) -> xr.DataArray:
        return self.zzarr.get_ds(z)["da"]


class Zzarr:
    def __init__(self, root_path: str, width: int, height: int):
        self.root_path = Path(root_path)
        self.width = width
        self.height = height
        self.z = None

    def open_zarr(self, mode: str, z: int) -> zarr.Array:
        path = self.root_path / str(z)
        array = zarr.open(
            path / "da",
            mode=mode,
            shape=(2**z * self.height, 2**z * self.width),
            chunks=(self.height, self.width),
            dtype="<f8",
        )
        if mode == "w":
            # write Dataset to zarr
            (path / ".zattrs").write_text(json.dumps(dict()))
            zarray = json.loads((path / "da" / ".zarray").read_text())
            zattrs = dict(_ARRAY_DIMENSIONS=["y", "x"])
            zgroup = dict(zarr_format=zarray["zarr_format"])
            (path / ".zgroup").write_text(json.dumps(zgroup))
            (path / ".zmetadata").write_text(
                json.dumps(
                    dict(
                        metadata={
                            ".zattrs": {},
                            ".zgroup": zgroup,
                            "da/.zarray": zarray,
                            "da/.zattrs": zattrs,
                        },
                        zarr_consolidated_format=1,
                    )
                )
            )
            (path / "da" / ".zattrs").write_text(json.dumps(zattrs))
        return array

    def write_to_zarr(self, tile: mercantile.Tile, data: np.ndarray):
        x, y, z = tile
        path = self.root_path / str(z)
        if path.exists():
            mode = "a"
        else:
            mode = "w"
        array = self.open_zarr(mode, z)
        array[
            y * self.height : (y + 1) * self.height,  # noqa
            x * self.width : (x + 1) * self.width,  # noqa
        ] = data

    def get_ds(self, z: int) -> xr.Dataset:
        path = self.root_path / str(z)
        if z != self.z:
            self.ds_z = xr.open_zarr(path)
            self.z = z
        return self.ds_z


def deg2idx(lat_deg, lon_deg, zoom, height, width, round_fun):
    lat_rad = math.radians(lat_deg)
    n = 2**zoom
    xtile = round_fun(((lon_deg + 180) % 360) / 360 * n * width)
    ytile = round_fun((1 - math.asinh(math.tan(lat_rad)) / math.pi) / 2 * n * height)
    return ytile, xtile
