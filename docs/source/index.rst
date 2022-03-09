Welcome to MagCoordPy's documentation!
======================================

To install:

.. code-block:: console

   pip install magcoordpy


.. toctree::
   :maxdepth: 2
   :caption: Contents:


coord_transforms.geodetic2cd
============================

Example usage
-------------

.. code-block:: python

   from magcoordpy import coord_transforms
   long_geo = np.arange(-180, 190, 10)
   lati_geo = np.zeros(len(long_array))
   alti_geo = np.zeros(len(long_array))
   lat_cd, lon_cd, r_cd = coord_transforms.geodetic2cd(lati_geo, long_geo, alti_geo, year=2021.0)


.. autofunction:: coord_transforms.geodetic2cd
   :members:


References
==========

Laundal, K.M., Richmond, A.D. Magnetic Coordinate Systems. Space Sci Rev 206, 27–59 (2017). <https://doi.org/10.1007/s11214-016-0275-y>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
