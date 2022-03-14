Orthogonal Magnetic Coordinate Systems
======================================

.. automodule:: coord_transforms
    :members:
    :show-inheritance:


Example usage
~~~~~~~~~~~~~

.. code-block:: python

   from magcoordpy import coord_transforms
   long_geo = np.arange(-180, 190, 10)
   lati_geo = np.zeros(len(long_array))
   alti_geo = np.zeros(len(long_array))
   lat_cd, lon_cd, r_cd = coord_transforms.geodetic2cd(lati_geo, long_geo, alti_geo, year=2021.0)
