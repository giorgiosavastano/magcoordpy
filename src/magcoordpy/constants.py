import os

from pymap3d.ellipsoid import Ellipsoid

from magcoordpy.igrf_coeffs import load_coeffs

igrf_file = os.path.dirname(os.path.abspath(__file__)) + "/data/igrf13coeffs.txt"
GH = load_coeffs(igrf_file)
# The World Geodetic System 1984 (WGS84) datum
earth_ellipsoid = Ellipsoid("wgs84")
# equatorial radius of ellipsoid
RE_M = earth_ellipsoid.semimajor_axis
# eccentricity squared of the ellipsoid
E_SQ = earth_ellipsoid.eccentricity ** 2
