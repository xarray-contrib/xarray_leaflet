import warnings
import numpy as np
import xarray as xr


def passthrough(array, *args, **kwargs):
    return array


def normalize(array, *args, **kwargs):
    vmin = np.min(array).values
    vmax = np.max(array).values
    array = (array - vmin) / (vmax - vmin)
    return array


def coarsen(agg_func=xr.core.rolling.DataArrayCoarsen.mean):
    def _(array, *args, **kwargs):
        tile_width = kwargs['tile_width']
        tile_height = kwargs['tile_height']
        ny, nx = array.shape
        wx = nx // (tile_width * 2)
        wy = ny // (tile_height * 2)
        dim = {}
        if wx > 1:
            dim['x'] = wx
        if wy > 1:
            dim['y'] = wy
        array = array.coarsen(**dim, boundary='pad')
        # ignore "mean of empty slice" warning in np.nanmean
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            array = agg_func(array)
        return array
    return _
