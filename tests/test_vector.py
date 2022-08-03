import filecmp
from tempfile import TemporaryDirectory

import geopandas as gpd
import mercantile
from shapely.geometry import box
from xarray_leaflet.vector import Zvect

from .utils import save_fig


HEIGHT = WIDTH = 256


def test_da_tile_full():
    x, y, z = 3, 4, 5
    bounds = mercantile.bounds(x, y, z)
    b = box(bounds.west, bounds.south, bounds.east, bounds.north)
    df = gpd.GeoDataFrame(geometry=gpd.GeoSeries([b]))
    df.set_crs(epsg=4326, inplace=True)
    measurement = "mask"
    df[measurement] = 1
    tile = mercantile.bounding_tile(*df.to_crs(epsg=4326).geometry.total_bounds)
    assert tile.x == x
    assert tile.y == y
    assert tile.z == z
    zvect = Zvect(df, measurement, WIDTH, HEIGHT)
    da_tile = zvect.get_da_tile(tile)
    assert len(da_tile.x) == WIDTH
    assert len(da_tile.y) == HEIGHT
    expected, result = save_fig(da_tile, "test_da_tile_full.png")
    assert filecmp.cmp(expected, result)


def test_da_tile_nybb():
    path_to_data = gpd.datasets.get_path("nybb")
    df = gpd.read_file(path_to_data)
    measurement = "mask"
    df[measurement] = 1
    tile = mercantile.bounding_tile(*df.to_crs(epsg=4326).geometry.total_bounds)
    zvect = Zvect(df, measurement, WIDTH, HEIGHT)
    da_tile = zvect.get_da_tile(tile)
    assert len(da_tile.x) == WIDTH
    assert len(da_tile.y) == HEIGHT
    expected, result = save_fig(da_tile, "test_da_tile_nybb.png")
    assert filecmp.cmp(expected, result)


def test_get_da_llbbox():
    path_to_data = gpd.datasets.get_path("nybb")
    df = gpd.read_file(path_to_data)
    measurement = "mask"
    df[measurement] = 1
    bounds = df.to_crs(epsg=4326).geometry.total_bounds
    tile = mercantile.bounding_tile(*bounds)
    llbbox = mercantile.LngLatBbox(*bounds)
    with TemporaryDirectory() as tmpdirname:
        zvect = Zvect(df, measurement, WIDTH, HEIGHT, tmpdirname)
        da = zvect.get_da_llbbox(llbbox, tile.z + 1)
        expected, result = save_fig(da, "test_get_da_llbbox.png")
    assert filecmp.cmp(expected, result)
