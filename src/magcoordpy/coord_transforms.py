import logging
from math import pi

import numpy as np
import pymap3d

import magcoordpy.constants as CONSTS
import magcoordpy.helpers as helpers
import magcoordpy.igrf_coeffs as igrf_coeffs

logger = logging.getLogger("spwx-coord")


def geodetic2ecef(geod_lati_deg, geod_long_deg, geod_alti_m):
    """
    Conversion from Geodetic (lat, lon, alt) to geocentric Cartesian (x, y, z) coordinates.
    Author: Giorgio Savastano (giorgio.savastano@spire.com)

    Parameters
    ----------
    geod_long_deg : np.ndarray
        array containing geodetic longitudes in degrees
    geod_lati_deg : np.ndarray
        array containing geodetic longitudes in degrees
    geod_alti_m : np.ndarray
        array containing geodetic altitude in meters

    Returns
    -------
    tuple : tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    geod_long = np.deg2rad(geod_long_deg)
    geod_lati = np.deg2rad(geod_lati_deg)

    rho = CONSTS.RE_M * (1 - CONSTS.E_SQ * (np.sin(geod_lati) ** 2)) ** (-0.5)

    x = (rho + geod_alti_m) * np.cos(geod_lati) * np.cos(geod_long)
    y = (rho + geod_alti_m) * np.cos(geod_lati) * np.sin(geod_long)
    z = (rho + geod_alti_m - CONSTS.E_SQ * rho) * np.sin(geod_lati)

    return x, y, z


def compute_pos_centered_dipole_north_pole(year=2021.0):
    """
    Function that returns the position of the centered dipole (CD) northern pole using first three Gauss coefficients.
    Coordinates of the CD northern pole are returned in Geocentric Spherical (colat, long) coordinates [rad].
    Author: Giorgio Savastano (giorgio.savastano@spire.com)


    Parameters
    ----------
        year : float
            year for the computation

    Returns
    -------
        tuple : tuple[float, float]
    """
    g, h = igrf_coeffs.get_coeffs(CONSTS.GH, year)

    if len(g) == 0:
        raise ValueError(" IGRF coefficients not correctly extracted")

    # first three Gauss coefficients (n = 1)
    g_n1_m0 = g[1][0]
    g_n1_m1 = g[1][1]
    h_n1_m1 = h[1][1]

    b0 = (g_n1_m0 ** 2 + g_n1_m1 ** 2 + h_n1_m1 ** 2) ** 0.5

    logger.debug(f" B0 = {b0}")

    # geocentric colatitude of the CD pole coordinates in NH
    geocentric_colat_n = np.arccos(-g_n1_m0 / b0)

    # geocentric longitude of the CD pole coordinates in NH
    geocentric_long_n = np.arctan2(h_n1_m1, g_n1_m1) - pi

    logger.info(
        f" Geocentric coordinates of the CD pole coordinates in {year}:"
        f" Colatitude: {np.round(np.rad2deg(geocentric_colat_n), 2)}"
        f" Longitude: {np.round(np.rad2deg(geocentric_long_n), 2)}"
    )

    return geocentric_colat_n, geocentric_long_n


def mlon2mlt(mlon_cd, datetime, ssheight=50 * 6371):
    """
    Computes the magnetic local time at the specified magnetic longitude and UT.
    Author: Giorgio Savastano (giorgio.savastano@spire.com) adapted from apexpy.

    Parameters
    ----------
    mlon_cd : array_like
        Centered Dipole (CD) magnetic longitude
    datetime : :class:`datetime.datetime`
        Date and time
    ssheight : float, optional
        Altitude in km to use for converting the subsolar point from
        geographic to magnetic coordinates. A high altitude is used
        to ensure the subsolar point is mapped to high latitudes, which
        prevents the South-Atlantic Anomaly (SAA) from influencing the MLT.
    Returns
    -------
    mlt : ndarray or float
        Magnetic local time [0, 24)
    Notes
    -----
    To compute the MLT, we find the apex longitude of the subsolar point at
    the given time. Then the MLT of the given point will be computed from
    the separation in magnetic longitude from this point (1 hour = 15
    degrees).
    """
    sslat_geo, sslon_geo = helpers.subsol(datetime)
    sslat_cd, sslon_cd, ssalt_cd = geodetic2cd(
        float(sslat_geo), float(sslon_geo), float(ssheight)
    )

    # np.float64 will ensure lists are converted to arrays
    return (180 + np.float64(mlon_cd) - sslon_cd) / 15 % 24


def ecef2eccdf(x_geoc, y_geoc, z_geoc, year=2021.0):
    """
    Transformation from Geocentric Cartesian (x, y, z) to Centered Dipole (CD) Cartesian (x, y, z) coordinates.
    Author: Giorgio Savastano (giorgio.savastano@spire.com)

    Parameters
    ----------
    x_geoc : np.ndarray
        array containing x coordinates
    y_geoc : np.ndarray
        array containing y coordinates
    z_geoc : np.ndarray
        array containing z coordinates
    year : float
        year for computing the IGRF Gaus coefficients

    Returns
    -------
    tuple : tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    if type(x_geoc) == list:
        logger.info(" Converting list to np.ndarrays.")
        x_geoc = np.asarray(x_geoc)
        y_geoc = np.asarray(y_geoc)
        z_geoc = np.asarray(z_geoc)
    elif type(x_geoc) != np.ndarray:
        logger.info(f" Converting {type(x_geoc)} to np.ndarrays.")
        x_geoc = np.asarray([x_geoc])
        y_geoc = np.asarray([y_geoc])
        z_geoc = np.asarray([z_geoc])

    # position of the centered dipole (CD) northern pole
    geoc_colat_n, geoc_long_n = compute_pos_centered_dipole_north_pole(year=year)

    # Compute the base vectors of the CD coordinate system
    # the transpose CD base vectors form the rotation matrix

    # m = unit vector anti-parallel to the dipole axis
    # this correspond to the z axis of the CD coordinate system
    # (Cartesian z axis aligns with the dipole axis)
    m = np.asarray(
        [
            np.sin(geoc_colat_n) * np.cos(geoc_long_n),
            np.sin(geoc_colat_n) * np.sin(geoc_long_n),
            np.cos(geoc_colat_n),
        ]
    )

    # x axis of the CD coordinate system
    xx = np.asarray(
        [
            np.cos(geoc_colat_n) * np.cos(geoc_long_n),
            np.cos(geoc_colat_n) * np.sin(geoc_long_n),
            -np.sin(geoc_colat_n),
        ]
    )

    # y axis of the CD coordinate system
    yy = np.asarray([-np.sin(geoc_long_n), np.cos(geoc_long_n), 0])
    # Rotation matrix
    R_geoc_cd = np.asarray([xx, yy, m])

    matrix_arrays = np.asarray([x_geoc, y_geoc, z_geoc]).T

    # multiplication between the rotation matrix and the array of vectors
    # I had to invert the order of the multiplication to take advantage of numpy vectorization
    x_cd_arr, y_cd_arr, z_cd_arr = (matrix_arrays @ R_geoc_cd.T).T

    return x_cd_arr, y_cd_arr, z_cd_arr


def ecef2spherical(x, y, z, decimals=2):
    """
    Transformation from geocentric Cartesian (x, y, z) to Spherical (colat, long, r) coordinates.
    Author: Giorgio Savastano (giorgio.savastano@spire.com)

    Parameters
    ----------
    x : np.ndarray
        array containing x coordinates
    y : np.ndarray
        array containing y coordinates
    z : np.ndarray
        array containing z coordinates
    decimals : int, optional
            Number of decimal places to round to (default: 2).  If
            decimals is negative, it specifies the number of positions to
            the left of the decimal point

    Returns
    -------
    tuple : tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    # geocentric radial distance 𝑟
    r = np.sqrt(x * x + y * y + z * z)
    # the colatitude 𝜃
    colat = np.round(np.rad2deg(np.arccos(z / r)), decimals)
    # the east longitude 𝜙
    long = np.round(np.rad2deg(np.arctan2(y, x)), decimals)

    return colat, long, np.round(r, decimals)


def spherical2ecef(colat_geoc_arr, long_geoc_arr, radial_dist_geoc_arr):
    """
    Transformation from Geocentric (GEO) Spherical (colat, lon, r) to Cartesian (x,y,z).
    Author: Giorgio Savastano (giorgio.savastano@spire.com)

    Args:
        colat_geoc_arr: np.ndarray
        array containing colatitude component of Geocentric coordinates in rad
        long_geoc_arr:
        array containing longitude component of Geocentric coordinates in rad
        radial_dist_geoc_arr:
        array containing radial distance component of Geocentric coordinates in m

    Returns:

    """
    x_geoc_arr = radial_dist_geoc_arr * np.sin(colat_geoc_arr) * np.cos(long_geoc_arr)
    y_geoc_arr = radial_dist_geoc_arr * np.sin(colat_geoc_arr) * np.sin(long_geoc_arr)
    z_geoc_arr = radial_dist_geoc_arr * np.cos(colat_geoc_arr)

    return x_geoc_arr, y_geoc_arr, z_geoc_arr


def geodetic2cd(
    gglat_deg_array, gglon_deg_array, ggalt_km_array, decimals=2, year=2021.0
):
    """Transformation from Geodetic to Centered Dipole (CD).

    Author: Giorgio Savastano (giorgio.savastano@spire.com)

    Parameters
    ----------
    gglon_deg_array : np.ndarray
        array containing geodetic longitude values in degrees
    gglat_deg_array : np.ndarray
        array containing geodetic latitude values in degrees
    ggalt_km_array : np.ndarray
        array containing geodetic altitude values in km
    decimals : int, default=2
        Number of decimal places to round to. If
        decimals is negative, it specifies the number of positions to
        the left of the decimal point.
    year : float, default=2021.0
        year for computing the IGRF Gauss coefficients

    Returns
    -------
    tuple : tuple[np.ndarray, np.ndarray, np.ndarray]
        CD lat, lon, alt arrays
    """
    if type(gglon_deg_array) == list:
        logger.info(" Converting list to np.ndarrays.")
        gglon_deg_array = np.asarray(gglon_deg_array)
        gglat_deg_array = np.asarray(gglat_deg_array)
        ggalt_km_array = np.asarray(ggalt_km_array)
    elif type(gglon_deg_array) != np.ndarray:
        logger.info(f" Converting {type(gglon_deg_array)} to np.ndarrays.")
        gglon_deg_array = np.asarray([gglon_deg_array])
        gglat_deg_array = np.asarray([gglat_deg_array])
        ggalt_km_array = np.asarray([ggalt_km_array])

    x_geoc, y_geoc, z_geoc = pymap3d.geodetic2ecef(
        gglat_deg_array, gglon_deg_array, ggalt_km_array * 1000.0
    )

    x_cd, y_cd, z_cd = ecef2eccdf(x_geoc, y_geoc, z_geoc, year=year)

    colat_cd, long_cd, r_cd = ecef2spherical(x_cd, y_cd, z_cd)
    lat_cd = np.round(90 - colat_cd, decimals)
    alt_cd = np.round(r_cd - CONSTS.RE_M, decimals)

    return lat_cd, long_cd, alt_cd


def eccdf2ecef(x_cd, y_cd, z_cd, year=2021.0):
    """
    Transformation from Centered Dipole (CD) Cartesian (x, y, z) to Geocentric Cartesian (x, y, z) coordinates.
    Author: Giorgio Savastano (giorgio.savastano@spire.com)

    Parameters
    ----------
    x_cd : np.ndarray
        array containing x component of CD coordinates
    y_cd : np.ndarray
        array containing y component of CD coordinates
    z_cd : np.ndarray
        array containing z component of CD coordinates
    year : float, default=2021.0
        year for computing the IGRF Gaus coefficients

    Returns
    -------
    tuple : tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    if type(x_cd) == list:
        logger.info(" Converting list to np.ndarrays.")
        x_cd = np.asarray(x_cd)
        y_cd = np.asarray(y_cd)
        z_cd = np.asarray(z_cd)
    elif type(x_cd) != np.ndarray:
        logger.info(f" Converting {type(x_cd)} to np.ndarrays.")
        x_cd = np.asarray([x_cd])
        y_cd = np.asarray([y_cd])
        z_cd = np.asarray([z_cd])

    # position of the centered dipole (CD) northern pole
    geoc_colat_n, geoc_long_n = compute_pos_centered_dipole_north_pole(year=year)

    # Compute the base vectors of the CD coordinate system
    # the transpose CD base vectors form the rotation matrix

    # m = unit vector anti-parallel to the dipole axis
    # this correspond to the z axis of the CD coordinate system
    # (Cartesian z axis aligns with the dipole axis)
    m = np.asarray(
        [
            np.sin(geoc_colat_n) * np.cos(geoc_long_n),
            np.sin(geoc_colat_n) * np.sin(geoc_long_n),
            np.cos(geoc_colat_n),
        ]
    )

    # x axis of the CD coordinate system
    xx = np.asarray(
        [
            np.cos(geoc_colat_n) * np.cos(geoc_long_n),
            np.cos(geoc_colat_n) * np.sin(geoc_long_n),
            -np.sin(geoc_colat_n),
        ]
    )

    # y axis of the CD coordinate system
    yy = np.asarray([-np.sin(geoc_long_n), np.cos(geoc_long_n), 0])
    # Rotation matrix is the transpose since the vectors are orthonormal
    R_cd_geoc = np.asarray([xx, yy, m]).T

    matrix_arrays = np.asarray([x_cd, y_cd, z_cd]).T

    # multiplication between the rotation matrix and the array of vectors
    # I had to invert the order of the multiplication to take advantage of numpy vectorization
    x_geoc_arr, y_geoc_arr, z_geoc_arr = (matrix_arrays @ R_cd_geoc.T).T

    return x_geoc_arr, y_geoc_arr, z_geoc_arr


def cd2geodetic(lat_cd_arr, lon_cd_arr, alt_cd_arr, decimals=3, year=2021.0):
    """Transformation from Centered Dipole (CD) to Geodetic.

    Author: Giorgio Savastano (giorgio.savastano@spire.com)

    Parameters
    ----------
    lat_cd_arr : np.ndarray
        array containing latitude component of CD coordinates in degrees
    lon_cd_arr : np.ndarray
        array containing longitude component of CD coordinates in degrees
    alt_cd_arr : np.ndarray
        array containing altitude (N.B. not the radial distance)
        component of CD coordinates in meters
    decimals : int, , default=3
        Number of decimal places to round to. If
        decimals is negative, it specifies the number of positions to
        the left of the decimal point
    year : float, default=2021.0
        year for computing the IGRF Gauss coefficients

    Returns
    -------
    tuple : tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    # convert lat to colat
    colat_cd_arr = 90 - lat_cd_arr
    rad_distance = alt_cd_arr + CONSTS.RE_M

    x_cd_arr, y_cd_arr, z_cd_arr = spherical2ecef(
        np.deg2rad(colat_cd_arr), np.deg2rad(lon_cd_arr), rad_distance
    )

    x_geoc_arr, y_geoc_arr, z_geoc_arr = eccdf2ecef(
        x_cd_arr, y_cd_arr, z_cd_arr, year=year
    )

    lat_geoc, long_geoc, alt_geoc = pymap3d.ecef2geodetic(
        x_geoc_arr, y_geoc_arr, z_geoc_arr
    )

    return (
        np.round(lat_geoc, decimals),
        np.round(long_geoc, decimals),
        np.round(alt_geoc, decimals),
    )
