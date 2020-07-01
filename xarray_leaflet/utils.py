import os
import numpy as np
from PIL import Image
import mercantile
from affine import Affine


def reproject_custom(source, dst_crs, x0, y0, z, x_res, y_res, width, height, resampling):
    a = x_res
    b = 0
    c = (x0 - 2 ** z) * width * x_res
    d = 0
    e = -y_res
    f = (2 ** z - y0) * height * y_res
    dst_affine = Affine(a, b, c, d, e, f)
    destination = source.rio.reproject(dst_crs, transform=dst_affine, shape=(height, width), resampling=resampling)
    return destination


def reproject_not_custom(source, dst_crs, x0, y0, x_res, y_res, width, height, resampling):
    a = x_res
    b = 0
    c = x0
    d = 0
    e = -y_res
    f = y0
    dst_affine = Affine(a, b, c, d, e, f)
    destination = source.rio.reproject(dst_crs, transform=dst_affine, shape=(height, width), resampling=resampling)
    return destination


def write_image(path, data, persist):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if data is None:
        open(path, 'wb').close()
    else:
        im = Image.fromarray(np.uint8(data))
        im.save(path)
    write_done_file(path, persist)


def write_done_file(png_path, persist):
    with open(png_path[:-4] + '.done', 'wt') as f:
        if persist:
            f.write('keep')
        else:
            f.write('delete')


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
