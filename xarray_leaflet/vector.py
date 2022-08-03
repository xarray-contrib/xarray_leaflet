import json
from functools import partial
from pathlib import Path
from typing import Optional

import mercantile
import numpy as np
import pyproj
import xarray as xr
import zarr
from geocube.api.core import make_geocube
from geocube.rasterize import rasterize_image
from geopandas import GeoDataFrame
from rasterio.enums import MergeAlg
from shapely.geometry import box, mapping
from shapely.ops import transform

from .utils import debug  # noqa


class Zvect:
    def __init__(
        self,
        df: GeoDataFrame,
        measurement: str,
        width: int,
        height: int,
        root_path: str = "",
    ):
        # reproject to Web Mercator
        self.df = df.to_crs(epsg=3857)
        self.measurement = measurement
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
            rasterize_function=partial(
                rasterize_image, merge_alg=MergeAlg.add, all_touched=True
            ),
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
        project = pyproj.Transformer.from_crs(
            pyproj.CRS("EPSG:4326"), pyproj.CRS("EPSG:3857"), always_xy=True
        ).transform
        b = box(*bbox)
        polygon = transform(project, b)
        left, bottom, right, top = polygon.bounds
        return self.zzarr.get_ds(z)["da"].sel(
            x=slice(left, right), y=slice(top, bottom)
        )

    def get_da(self, z: int) -> xr.DataArray:
        return self.zzarr.get_ds(z)["da"]


class Zzarr:
    def __init__(self, root_path: str, width: int, height: int):
        self.root_path = Path(root_path)
        self.width = width
        self.height = height
        self.ds = {}

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
            mi, ma = mercantile.minmax(z)
            ul = mercantile.xy_bounds(mi, mi, z)
            lr = mercantile.xy_bounds(ma, ma, z)
            bbox = mercantile.Bbox(ul.left, lr.bottom, lr.right, ul.top)
            x = zarr.open(
                path / "x",
                mode="w",
                shape=(2**z * self.width,),
                chunks=(2**z * self.width,),
                dtype="<f8",
            )
            x[:] = np.linspace(bbox.left, bbox.right, 2**z * self.width)
            x_zattrs = dict(_ARRAY_DIMENSIONS=["x"])
            (path / "x" / ".zattrs").write_text(json.dumps(x_zattrs))
            y = zarr.open(
                path / "y",
                mode="w",
                shape=(2**z * self.height,),
                chunks=(2**z * self.height,),
                dtype="<f8",
            )
            y[:] = np.linspace(bbox.top, bbox.bottom, 2**z * self.height)
            x_zarray = json.loads((path / "x" / ".zarray").read_text())
            y_zarray = json.loads((path / "y" / ".zarray").read_text())
            y_zattrs = dict(_ARRAY_DIMENSIONS=["y"])
            (path / "y" / ".zattrs").write_text(json.dumps(y_zattrs))
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
                            "x/.zarray": x_zarray,
                            "x/.zattrs": x_zattrs,
                            "y/.zarray": y_zarray,
                            "y/.zattrs": y_zattrs,
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
        self.array = self.open_zarr(mode, z)
        self.array[
            y * self.height : (y + 1) * self.height,  # noqa
            x * self.width : (x + 1) * self.width,  # noqa
        ] = data

    def get_ds(self, z: int) -> xr.Dataset:
        path = self.root_path / str(z)
        if z not in self.ds:
            self.ds[z] = xr.open_zarr(path)
        return self.ds[z]
