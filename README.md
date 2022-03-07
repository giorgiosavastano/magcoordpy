# MagCoordPy

![test-main](https://github.com/giorgiosavastano/magcoordpy/actions/workflows/python-test-main.yml/badge.svg)
![coverage-main](https://img.shields.io/codecov/c/github/giorgiosavastano/magcoordpy)
![license](https://img.shields.io/github/license/giorgiosavastano/magcoordpy)

A python package for working with geomagnetic coordinates.
The documentation is available at https://magcoordpy.readthedocs.io/en/latest/.

Installation
------------

    pip install magcoordpy

Example usage
-------------

    import magcoordpy as mcp
    long_geo = np.arange(-180, 190, 10)
    lati_geo = np.zeros(len(long_array))
    alti_geo = np.zeros(len(long_array))
    lat_cd, long_cd, r_cd = mcp.geodetic2cd(
            lati_array, long_array, alti_array, year=2021.0
        )


Authors:

- Giorgio Savastano (<giorgiosavastano@gmail.com>)

Please use github issues to make bug reports and request new functionality. Contributions are always welcome.
