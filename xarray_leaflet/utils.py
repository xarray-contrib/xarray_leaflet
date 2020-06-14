import os
import numpy as np
from PIL import Image
import mercantile
from affine import Affine


def reproject_custom(source, dst_crs, x0, y0, z, resolution, width, height):
    a = resolution
    b = 0
    c = (x0 - 2 ** z) * width * resolution
    d = 0
    e = -resolution
    f = (2 ** z - y0) * height * resolution
    dst_affine = Affine(a, b, c, d, e, f)
    destination = source.rio.reproject(dst_crs, dst_affine_width_height=(dst_affine, width, height))
    return destination


def reproject_not_custom(source, dst_crs, width, height):
    destination = source.rio.reproject(dst_crs, shape=(height, width))
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
