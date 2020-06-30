===============
Getting started
===============

In order to show a DataArray (say, ``da``) on a Leaflet map, all you need to do is create a map:

.. code::

    import xarray_leaflet
    from ipyleaflet import Map

    m = Map()
    m

And pass it to the plotting function, which returns a layer:

.. code::

    l = da.leaflet.plot(m)

You can then interact with the layer and e.g. control its opacity:

.. code::

    l.interact(opacity=(0., 1.)

For more advanced examples, please have a look at the `Notebooks <https://github.com/davidbrochart/xarray_leaflet/tree/master/examples>`_.
