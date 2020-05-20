import numpy as np
import xarray as xr


def passthrough(array, *args):
    return array


def normalize(array, *args):
    vmin = np.min(array).values
    vmax = np.max(array).values
    array = (array - vmin) / (vmax - vmin)
    return array


def coarsen(agg_func=xr.core.rolling.DataArrayCoarsen.mean):
    def _(array):
        ny, nx = array.shape
        wx = nx // (256 // 2)
        wy = ny // (256 // 2)
        dim = {}
        if wx > 1:
            dim['x'] = wx
        if wy > 1:
            dim['y'] = wy
        array = array.coarsen(**dim, boundary='pad')
        array = agg_func(array)
        return array
    return _
