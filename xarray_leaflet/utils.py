import os
import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from PIL import Image
import mercantile
from affine import Affine


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


def write_image(data_array, path, persist):
    im = Image.fromarray(np.uint8(data_array*255))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    im.save(path)
    write_done_file(path, persist)


def write_done_file(png_path, persist):
    with open(png_path[:-4] + '.done', 'wt') as f:
        if persist:
            f.write('keep')
        else:
            f.write('delete')


def reindex(da_tile, bbox, dx, dy):
    do_reindex = False
    lons = da_tile.x.values
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
    lats = da_tile.y.values
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
        da_tile = da_tile.reindex(y=lats, x=lons, method='nearest', tolerance=dx/4)
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
