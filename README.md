# MagCoordPy

![test-main](https://github.com/giorgiosavastano/magcoordpy/actions/workflows/python-test-main.yml/badge.svg)
![coverage-main](https://img.shields.io/codecov/c/github/giorgiosavastano/magcoordpy)
![license](https://img.shields.io/github/license/giorgiosavastano/magcoordpy)

A python package for working with magnetic coordinate transformations.
The documentation is available at https://magcoordpy.readthedocs.io/en/latest/.

Installation
------------

    pip install magcoordpy

Example usage
-------------

    from magcoordpy import coord_transforms
    long_geo = np.arange(-180, 190, 10)
    lati_geo = np.zeros(len(long_array))
    alti_geo = np.zeros(len(long_array))
    lat_cd, lon_cd, r_cd = coord_transforms.geodetic2cd(lati_geo, long_geo, alti_geo,
    													year=2021.0)

It includes the following functions (not exhaustive list):

* geodetic2cd
* cd2geodetic


Authors:

- Giorgio Savastano (<giorgiosavastano@gmail.com>)

Please use github issues to make bug reports and request new functionality. Contributions are always welcome.
