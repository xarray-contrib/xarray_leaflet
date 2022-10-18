import asyncio
import os

import mercantile
import numpy as np
from affine import Affine
from PIL import Image


def reproject_custom(
    source, dst_crs, x0, y0, z, x_res, y_res, width, height, resampling
):
    a = x_res
    b = 0
    c = (x0 - 2**z) * width * x_res
    d = 0
    e = -y_res
    f = (2**z - y0) * height * y_res
    dst_affine = Affine(a, b, c, d, e, f)
    destination = source.rio.reproject(
        dst_crs, transform=dst_affine, shape=(height, width), resampling=resampling
    )
    return destination


def reproject_not_custom(
    source, dst_crs, x0, y0, x_res, y_res, width, height, resampling
):
    a = x_res
    b = 0
    c = x0
    d = 0
    e = -y_res
    f = y0
    dst_affine = Affine(a, b, c, d, e, f)
    destination = source.rio.reproject(
        dst_crs, transform=dst_affine, shape=(height, width), resampling=resampling
    )
    return destination


def write_image(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    status_path = path[:-4] + ".status"
    with open(status_path, "wt") as f:
        f.write("computing")
    if data is None:
        open(path, "wb").close()
    else:
        im = Image.fromarray(np.uint8(data))
        im.save(path)
    with open(status_path, "wt") as f:
        f.write("done")


def get_bbox_tiles(tiles):
    north = east = -float("inf")
    south = west = float("inf")
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


def wait_for_change(widget, value):
    future = asyncio.Future()

    def get_value(change):
        future.set_result(change.new)
        widget.unobserve(get_value, value)

    widget.observe(get_value, value)
    return future


def lng_to_180(lng: float) -> float:
    lng1 = lng + 180
    lng2 = lng1 % 360
    if lng1 < 0:
        lng2 = 360 - lng2
    lng2 -= 180
    return lng2


def debug(message: str):
    with open("debug.txt", "a") as f:
        f.write(f"{message}\n")
