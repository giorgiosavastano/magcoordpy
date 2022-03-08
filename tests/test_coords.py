import datetime as dt
from pathlib import Path
from unittest import TestCase

import numpy as np
import pytest

from magcoordpy.coord_transforms import (
    cd2geodetic,
    compute_pos_centered_dipole_north_pole,
    ecef2eccdf,
    geodetic2cd,
    geodetic2ecef,
    mlon2mlt,
)

path = Path(__file__)


class TestCoords(TestCase):
    def setUp(self):
        pass

    @pytest.mark.unit
    def test_geodetic2cd(self):
        """
        Using World Data Center for Geomagnetism, Kyoto (http://wdc.kugi.kyoto-u.ac.jp/igrf/gggm/) for validation.
        Author: Giorgio Savastano (giorgio.savastano@spire.com)
        """
        long_array = np.arange(-180, 190, 10)
        lati_array = np.zeros(len(long_array))
        alti_array = np.zeros(len(long_array))
        lat_arr, long_arr, r_arr = geodetic2cd(
            lati_array, long_array, alti_array, year=2021.0
        )
        self.assertEqual(lat_arr[0], -2.78)
        self.assertEqual(long_arr[0], -107.11)
        self.assertEqual(lat_arr[1], -1.19)
        self.assertEqual(long_arr[1], -97.23)
        self.assertEqual(lat_arr[2], 0.43)
        self.assertEqual(long_arr[2], -87.36)
        self.assertEqual(lat_arr[3], 2.04)
        self.assertEqual(long_arr[3], -77.49)
        self.assertEqual(lat_arr[4], 3.59)
        self.assertEqual(long_arr[4], -67.60)
        self.assertEqual(lat_arr[5], 5.04)
        self.assertEqual(long_arr[5], -57.68)
        self.assertEqual(lat_arr[18], 2.78)
        self.assertEqual(long_arr[18], 72.89)
        self.assertEqual(lat_arr[26], -8.31)
        self.assertEqual(long_arr[26], 152.36)
        self.assertEqual(lat_arr[34], -5.66)
        self.assertEqual(long_arr[34], -126.96)
        self.assertEqual(lat_arr[36], -2.78)
        self.assertEqual(long_arr[36], -107.11)

        lati_array = np.arange(-90, 100, 10)
        long_array = -np.ones(len(lati_array)) * 70
        alti_array = np.zeros(len(long_array))
        lat_arr, long_arr, r_arr = geodetic2cd(
            lati_array, long_array, alti_array, year=2021.0
        )
        self.assertEqual(lat_arr[0], -80.64)
        self.assertEqual(long_arr[0], 0)
        self.assertEqual(lat_arr[1], -70.58)
        self.assertEqual(long_arr[1], 1.40)
        self.assertEqual(lat_arr[2], -60.52)
        self.assertEqual(long_arr[2], 1.87)
        self.assertEqual(lat_arr[3], -50.48)
        self.assertEqual(long_arr[3], 2.11)
        self.assertEqual(lat_arr[4], -40.46)
        self.assertEqual(long_arr[4], 2.27)
        self.assertEqual(lat_arr[9], 9.35)
        self.assertEqual(long_arr[9], 2.71)
        self.assertEqual(lat_arr[10], 19.29)
        self.assertEqual(long_arr[10], 2.79)
        self.assertEqual(lat_arr[11], 29.23)
        self.assertEqual(long_arr[11], 2.88)
        self.assertEqual(lat_arr[17], 89.16)
        self.assertEqual(long_arr[17], 33.93)
        self.assertEqual(lat_arr[18], 80.64)
        self.assertEqual(long_arr[18], 180.00)

        # test list as input
        lat_arr, long_arr, r_arr = geodetic2cd([80.0], [110.0], [10.123], year=2021.0)
        self.assertEqual(lat_arr[0], 70.58)
        self.assertEqual(long_arr[0], -178.60)

        # test float as input
        lat_arr, long_arr, r_arr = geodetic2cd(80.0, 110.0, 10.123, year=2021.0)
        self.assertEqual(lat_arr[0], 70.58)
        self.assertEqual(long_arr[0], -178.60)

        # test int as input
        lat_arr, long_arr, r_arr = geodetic2cd(80, 110, 10, year=2021.0)
        self.assertEqual(lat_arr[0], 70.58)
        self.assertEqual(long_arr[0], -178.60)

    @pytest.mark.unit
    def test_ecef2eccdf(self):
        x, y, z = geodetic2ecef(80, 10, 20)

        x_cd, y_cd, z_cd = ecef2eccdf([x], [y], [z])
        self.assertEqual(type(x_cd), np.ndarray)
        self.assertEqual(type(y_cd), np.ndarray)
        self.assertEqual(type(z_cd), np.ndarray)

        x_cd, y_cd, z_cd = ecef2eccdf(int(x), int(y), int(z))
        self.assertEqual(type(x_cd), np.ndarray)
        self.assertEqual(type(y_cd), np.ndarray)
        self.assertEqual(type(z_cd), np.ndarray)

        x_cd, y_cd, z_cd = ecef2eccdf(float(x), float(y), float(z))
        self.assertEqual(type(x_cd), np.ndarray)
        self.assertEqual(type(y_cd), np.ndarray)
        self.assertEqual(type(z_cd), np.ndarray)

    @pytest.mark.unit
    def test_position_centered_dipole_northern_pole(self):
        """
        Test that the code return the right position of the centered dipole northern pole.
        Validation with Table I (up to 2000) in https://www.spenvis.oma.be/help/background/magfield/cd.html.
        Author: Giorgio Savastano (giorgio.savastano@spire.com)
        """
        colat_rad, lon_rad = compute_pos_centered_dipole_north_pole(2015)
        colat_deg = round(np.rad2deg(colat_rad), 2)
        long_deg = round(np.rad2deg(lon_rad), 2)
        self.assertEqual(colat_deg, 9.69)
        self.assertEqual(long_deg, -72.61)

        colat_rad, lon_rad = compute_pos_centered_dipole_north_pole(2000)
        colat_deg = np.rad2deg(colat_rad)
        long_deg = np.rad2deg(lon_rad)
        lat_deg_comp = round(90 - colat_deg, 2)
        lon_deg_comp = round(360 + long_deg, 2)
        self.assertEqual(lat_deg_comp, 79.54)
        self.assertEqual(lon_deg_comp, 288.43)

    @pytest.mark.unit
    def test_cd2geodetic(self):

        # check value of the magnetic equator
        long_array = np.arange(-180, 190, 10)
        lati_array = np.zeros(len(long_array))
        alti_array = np.zeros(len(long_array))

        lat_geo_arr, lon_geo_arr, r_arr = cd2geodetic(
            lati_array, long_array, alti_array, decimals=3, year=2021.0
        )
        self.assertEqual(lat_geo_arr[0], 9.423)
        self.assertEqual(lon_geo_arr[0], 107.328)
        self.assertEqual(lat_geo_arr[1], 9.279)
        self.assertEqual(lon_geo_arr[1], 117.460)
        self.assertEqual(lat_geo_arr[8], 1.629)
        self.assertEqual(lon_geo_arr[8], -172.542)
        self.assertEqual(lat_geo_arr[13], -6.042)
        self.assertEqual(lon_geo_arr[13], -123.050)
        self.assertEqual(lat_geo_arr[18], -9.423)
        self.assertEqual(lon_geo_arr[18], -72.672)
        self.assertEqual(lat_geo_arr[19], -9.279)
        self.assertEqual(lon_geo_arr[19], -62.540)
        self.assertEqual(lat_geo_arr[24], -4.696)
        self.assertEqual(lon_geo_arr[24], -12.341)
        self.assertEqual(lat_geo_arr[30], 4.696)
        self.assertEqual(lon_geo_arr[30], 46.996)
        self.assertEqual(lat_geo_arr[36], 9.423)
        self.assertEqual(lon_geo_arr[36], 107.328)

    @pytest.mark.unit
    def test_mlon2mlt(self):
        """
        Comparison with apexpy results.
        Returns
        -------

        """
        dat = dt.datetime.fromisoformat("2021-12-09 12:00:00")
        mlt_cd = np.round(mlon2mlt(10.12, dat), 2)
        self.assertEqual(mlt_cd, 8.18)

        mlt_cd = np.round(mlon2mlt(-10.12, dat), 2)
        self.assertEqual(mlt_cd, 6.83)

        long_array = np.arange(-180, 190, 10)
        mlt_cd = np.round(mlon2mlt(long_array, dat), 2)

        apexpy_results = np.array(
            [
                19.51,
                20.17,
                20.84,
                21.51,
                22.17,
                22.84,
                23.51,
                0.17,
                0.84,
                1.51,
                2.17,
                2.84,
                3.51,
                4.17,
                4.84,
                5.51,
                6.17,
                6.84,
                7.51,
                8.17,
                8.84,
                9.51,
                10.17,
                10.84,
                11.51,
                12.17,
                12.84,
                13.51,
                14.17,
                14.84,
                15.51,
                16.17,
                16.84,
                17.51,
                18.17,
                18.84,
                19.51,
            ]
        )

        self.assertTrue((mlt_cd == apexpy_results).all())

    def tearDown(self):
        pass
