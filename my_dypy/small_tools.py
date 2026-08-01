# -*- coding:utf-8 -*-
import math
from datetime import timedelta

import cartopy
import copy
import fiona
import numpy as np
import numpy.ma as ma
import scipy as sp
from PIL import Image, ImageDraw
from fiona import crs
from pyproj import Geod

from dypy.constants import rdv, L, R, cp, t_zero, g
from dypy.intergrid import Intergrid

__all__ = [
    "CrossSection",
    "interpolate",
    "rotate_points",
    "interval",
    "dewpoint",
    "esat",
    "moist_lapse",
    "equivalent_pot_temp",
    "mask_polygon",
    "potential_temperature",
    "great_circle_distance",
    "geodetic_distance",
    "read_shapefile",
]


class NotAllRequiredVarsException(Exception):
    pass


class CrossSection(object):
    """Create a cross-section

    return the distance, the pressure levels,
    and the variables.

    Parameters
    ----------

    variables : dictionary
        A dictionary of the variables with the following structure:
        {'name': np.array}.
        If version is rotated: need to contain as least rlon, rlat
        if version is regular: need to contain at least lon, lat
    coo : list or numpy.ndarray
        If `coo` is a list of coordinates, it is used of the starting and
        end points: [(startlon, startlat), (endlon, endlat)]
        If `coo` is a numpy.ndarray it is used as the cross-section points.
        coo need then to be similar to :
        np.array([[10, 45], [11, 46], [12, 47]])
    pressure : np.array
        pressure coordinate as 1D array
    int2p : bool, optional (default: False)
        True if variables need to be interpolated on pressure levels,
        requires p in the variables, the levels are given by `pressure`
    int2z : bool, optional (default: False)
        True if variables need to be interpolated on heights levels,
        requires z in the variables, the levels are given by `pressure`
        may require `flip`=True
    flip : bool, optional (default: False)
        True if variables need to be flip vertically,
        the first dimension of 3D array will be flipped
    pollon : float, optional (default: -170)
        pole_longitude of the rotated grid
    pollat : float, optional (default: 43)
        pole_latitude of the rotated grid
    nbre : int, optional (default: 10000)
        nbre of points along the cross-section
    order: {1, 2, 3, 4}, optional (default: 1)
        order of interpolation see for details
    version: string (defautl rotated)
        type of grid used as input (rotated or regular)
    modification by LJA: add possibility to return grid
        """

    def __init__(
        self,
        variables,
        coo,
        pressure,
        int2p=False,
        int2z=False,
        flip=False,
        pollon=-170,
        pollat=43,
        nbre=1000,
        order=1,
        return_grid=False,
        version="rotated",
    ):

        variables = copy.deepcopy(variables)
        # test the input
        if version == "rotated":
            required = {"rlon", "rlat"}
        elif version == "regular":
            required = {"lon", "lat"}
        else:
            msg = "only rotated or regular are supported"
            raise NotImplementedError(msg)

        if int2p:
            required.add("p")
        if int2z:
            required.add("z")

        diff = required.difference(set(variables.keys()))

        if diff:
            msg = "{} not in variables dictionary".format(diff)
            raise NotAllRequiredVarsException(msg)

        if version == "rotated":
            if return_grid:
                rlon = variables['rlon']
                rlat = variables['rlat']
            else:
                rlon = variables.pop('rlon')
                rlat = variables.pop('rlat')
        elif version == "regular":
            if return_grid:
                lon = variables['lon']
                lat = variables['lat']
            else:
                lon = variables.pop('lon')
                lat = variables.pop('lat')

        self.pressure = pressure
        self.nz = pressure.size
        self.order = order

        if type(coo) in [list, tuple]:
            self.coo = coo
            self.nbre = nbre
            # find the coordinate of the cross-section
            if version == "rotated":
                (self.lon, self.lat, crlon, crlat, self.distance) = find_cross_points(
                    coo, nbre, pollon, pollat
                )
            elif version == "regular":
                self.lon, self.lat, self.distance = find_cross_points(coo, nbre)
                crlon, crlat = self.lon, self.lat
        elif type(coo) == np.ndarray:
            self.lon = coo[:, 0]
            self.lat = coo[:, 1]
            self.nbre = coo.shape[0]
            if version == "rotated":
                crlon, crlat = rotate_points(pollon, pollat, self.lon, self.lat)
            elif version == "regular":
                crlon, crlat = self.lon, self.lat
            self.distance = geodetic_distance(
                coo[0, 0], coo[0, 1], coo[-1, 0], coo[-1, 1]
            )
        else:
            msg = "<coo> should be of type list, tuple or nump.ndarray"
            raise TypeError(msg)

        self.distances = np.linspace(0, self.distance, self.nbre)

        # get the information for intergrid
        if version == "rotated":
            self._lo = np.array([rlat.min(), rlon.min()])
            self._hi = np.array([rlat.max(), rlon.max()])
        elif version == "regular":
            self._lo = np.array([lat.min(), lon.min()])
            self._hi = np.array([lat.max(), lon.max()])
        self._query_points = [[lat, lon] for lat, lon in zip(crlat, crlon)]

        if int2p or int2z:
            if int2p:
                grid = variables["p"]
            if int2z:
                grid = variables["z"]
            if flip and grid.ndim == 3:
                grid = grid[::-1, ...]
            for var, data in variables.items():
                if data.ndim > 2:
                    if flip:
                        data = data[::-1, ...]
                    dataint = interpolate(data, grid, pressure)
                    variables[var] = dataint

        variables = self._int2cross(variables)
        self.__dict__.update(variables)

    def _int2cross(self, variables):
        for var, data in variables.items():
            if data.squeeze().ndim > 2:
                datacross = np.zeros((self.nz, self.nbre))
                for i, datah in enumerate(data.squeeze()):
                    intfunc = Intergrid(
                        datah, lo=self._lo, hi=self._hi, verbose=False, order=self.order
                    )
                    datacross[i, :] = intfunc.at(self._query_points)
            else:
                datacross = np.zeros((self.nbre,))
                intfunc = Intergrid(
                    data.squeeze(),
                    lo=self._lo,
                    hi=self._hi,
                    verbose=False,
                    order=self.order,
                )
                datacross[:] = intfunc.at(self._query_points)
            variables[var] = datacross
        return variables


def find_cross_points(coo, nbre, pole_longitude=None, pole_latitude=None):
    """ give nbre points along a great circle between coo[0] and coo[1]
        iand rotated the points

        modification by LJA: add the starting point and end point to the cross points
    """

    g = Geod(ellps="WGS84")
    cross_points = g.npts(coo[0][0], coo[0][1], coo[1][0], coo[1][1], nbre - 2)
    cross_points.insert(0,tuple([coo[0][0], coo[0][1]]))
    cross_points.append(tuple([coo[1][0], coo[1][1]]))

    lat = np.array([point[1] for point in cross_points])
    lon = np.array([point[0] for point in cross_points])
    #distance = great_circle_distance(coo[0][0], coo[0][1], coo[1][0], coo[1][1])
    distance = geodetic_distance(cross_points[0][0], cross_points[0][1], cross_points[-1][0], cross_points[-1][1])
    if (pole_latitude is None) and (pole_latitude is None):
        return lon, lat, distance
    else:
        rlon, rlat = rotate_points(pole_longitude, pole_latitude, lon, lat)
        return lon, lat, rlon, rlat, distance


def great_circle_distance(lon1, lat1, lon2, lat2):
    """Calculate the great circle distance between two points

    based on : https://gist.github.com/gabesmed/1826175


    Parameters
    ----------

    lon1: float
        longitude of the starting point
    lat1: float
        latitude of the starting point
    lon2: float
        longitude of the ending point
    lat2: float
        latitude of the ending point

    Returns
    -------

    distance (km): float

    Examples
    --------

    >>> great_circle_distance(0, 55, 8, 45.5)
    1199.3240879770135
    """

    EARTH_CIRCUMFERENCE = 6378137  # earth circumference in meters

    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.cos(
        math.radians(lat1)
    ) * math.cos(math.radians(lat2)) * math.sin(dLon / 2) * math.sin(dLon / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = EARTH_CIRCUMFERENCE * c

    return distance / 1000

def geodetic_distance(lon1, lat1, lon2, lat2):
    """

    added by LJA

    Calculate the WGS84 distance between two points

    based on : https://pyproj4.github.io/pyproj/stable/api/geod.html


    Parameters
    ----------

    lon1: float
        longitude of the starting point
    lat1: float
        latitude of the starting point
    lon2: float
        longitude of the ending point
    lat2: float
        latitude of the ending point

    input can be number or numpy array

    Returns
    -------

    distance (km): np.float

    Examples
    --------

    >>> geodetic_distance(0, 55, 8, 45.5)
    1199.1306595506408
    """

    if type(lon1) != np.ndarray:
        lon1 = np.array(lon1)
        lon2 = np.array(lon2)
        lat1 = np.array(lat1)
        lat2 = np.array(lat2)

    g = Geod(ellps="WGS84")
    distance = g.inv(lon1,lat1,lon2,lat2)[-1]

    return distance / 1000


def interpolate(data, grid, interplevels):
    """interpolate `data` on `grid` for given `interplevels`

    Interpolate the `data` array at every given levels (`interplevels`) using
    the `grid` array as reference.

    The `grid` array need to be increasing along the first axis.
    Therefore ERA-Interim pressure must be flip: p[::-1, ...], but not
    COSMO pressure since its level 0 is located at the top of the atmosphere.

    Parameters
    ----------
    data : array (nz, nlat, nlon)
        data to interpolate
    grid : array (nz, nlat, nlon)
        grid use to perform the interpolation
    interplevels: list, array
        list of the new vertical levels, in the same unit as the grid

    Returns
    -------
    interpolated array: array (len(interplevels), nlat, nlon)

    Examples
    --------

    >>> print(qv.shape)
    (60, 181, 361)
    >>> print(p.shape)
    (60, 181, 361)
    >>> levels = np.arange(200, 1050, 50)
    (17,)
    >>> qv_int = interpolate(qv, p, levels)
    >>> print(qv_int.shape)
    (17, 181, 361)

    """
    data = data.squeeze()
    grid = grid.squeeze()
    shape = list(data.shape)
    if (data.ndim > 3) | (grid.ndim > 3):
        message = "data and grid need to be 3d array"
        raise IndexError(message)

    try:
        nintlev = len(interplevels)
    except Exception:
        interplevels = [interplevels]
        nintlev = len(interplevels)
    shape[-3] = nintlev

    outdata = np.ones(shape) * np.nan
    if nintlev > 20:
        for idx, _ in np.ndenumerate(data[0]):
            column = grid[:, idx[0], idx[1]]
            column_GRID = data[:, idx[0], idx[1]]

            value = np.interp(
                interplevels, column, column_GRID, left=np.nan, right=np.nan
            )
            outdata[:, idx[0], idx[1]] = value[:]
    else:
        for j, intlevel in enumerate(interplevels):
            for lev in range(grid.shape[0]):
                cond1 = grid[lev, :, :] > intlevel
                cond2 = grid[lev - 1, :, :] < intlevel
                right = np.where(cond1 & cond2)
                if right[0].size > 0:
                    sabove = grid[lev, right[0], right[1]]
                    sbelow = grid[lev - 1, right[0], right[1]]
                    dabove = data[lev, right[0], right[1]]
                    dbelow = data[lev - 1, right[0], right[1]]
                    result = (intlevel - sbelow) / (sabove - sbelow) * (
                        dabove - dbelow
                    ) + dbelow
                    outdata[j, right[0], right[1]] = result
    return outdata


def rotate_points(pole_longitude, pole_latitude, lon, lat, direction="n2r"):
    """Rotate lon, lat from/to a rotated system

    Parameters
    ----------
    pole_longitude: float
        longitudinal coordinate of the rotated pole
    pole_latitude: float
        latitudinal coordinate of the rotated pole
    lon: array (1d)
        longitudinal coordinates to rotate
    lat: array (1d)
        latitudinal coordinates to rotate
    direction: string, optional
        direction of the rotation;
        n2r: from non-rotated to rotated (default)
        r2n: from rotated to non-rotated

    Returns
    -------
    rlon: array
    rlat: array
    """
    lon = np.array(lon)
    lat = np.array(lat)
    rotatedgrid = cartopy.crs.RotatedPole(
        pole_longitude=pole_longitude, pole_latitude=pole_latitude
    )
    standard_grid = cartopy.crs.Geodetic()

    if direction == "n2r":
        rotated_points = rotatedgrid.transform_points(standard_grid, lon, lat)
    elif direction == "r2n":
        rotated_points = standard_grid.transform_points(rotatedgrid, lon, lat)

    rlon, rlat, _ = rotated_points.T
    return rlon, rlat


def interval(starting_date, ending_date, delta=timedelta(hours=1)):
    curr = starting_date
    while curr < ending_date:
        yield curr
        curr += delta


def dewpoint(p, qv):
    """Calculate the dew point temperature

    following (eq.8):
    Lawrence, M.G., 2005. The Relationship between Relative Humidity
    and the Dewpoint Temperature in Moist Air:
    A Simple Conversion and Applications.
    Bulletin of the American Meteorological Society 86,
    225–233. doi:10.1175/BAMS-86-2-225

    Parameters
    ----------
    p: float, array
        pressure in Pa
    qv: float, array
        specific humidity in kg/kg

    Returns
    -------
    dewpoint (in °C): float, array
    """
    B1 = 243.04  # °C
    C1 = 610.94  # Pa
    A1 = 17.625

    e = p * qv / (rdv + qv)
    td = (B1 * np.log(e / C1)) / (A1 - np.log(e / C1))

    return td


def esat(t):
    """
    Calculate the saturation vapor pressure for t in °C

    Following eq. 6 of
    Lawrence, M.G., 2005. The Relationship between Relative Humidity
    and the Dewpoint Temperature in Moist Air:
    A Simple Conversion and Applications.
    Bulletin of the American Meteorological Society 86,
    225–233. doi:10.1175/BAMS-86-2-225

    Parameters
    ----------
    t: float, numpy.array
        temperature in K

    Examples
    --------

    >>> esat(0)
    610.94000000000005

    """
    C1 = 610.94  # Pa
    A1 = 17.625  #
    B1 = 243.03  # °C

    return C1 * np.exp((A1 * t) / (B1 + t))


def esat_ice(t):
    """calculate the saturation vapor pressure over ice for t in K

    Follows COSMO, line 2468 in src_gscp.f90
    constant values are from data_constants.f90

    Parameters
    ----------
    t: float, numpy.array
        temperature in K

    Returns
    -------
    esat_ice: float, numpy array
        saturation vapor pressure over ice in Pa

    Examples
    --------

    >>> esat_ice(273.15)
    610.27696635637164


    """
    b1 = 610.78
    b2i = 21.8745584
    b3 = 273.16
    b4i = 7.66
    return b1 * np.exp(b2i * (t - b3) / (t - b4i))


def moist_lapse(t, p):
    """ Calculates moist adiabatic lapse rate

    Note: We calculate dT/dp, not dT/dz
    See formula 3.16 in Rogers&Yau for dT/dz, but this must be combined with
    the dry adiabatic lapse rate (gamma = g/cp) and the
    inverse of the hydrostatic equation (dz/dp = -RT/pg)

    Parameters
    ----------
    t : float, array
        temperature in degree Celsius
    p : float, array
        temperature in Pa

    Returns
    -------
    moist adiabatic lapse rate: float, array
    """

    a = 2.0 / 7.0
    b = rdv * L * L / (R * cp)
    c = a * L / R
    e = esat(t)
    wsat = rdv * e / (p - e)  # Rogers&Yau 2.18
    numer = a * (t + t_zero) + c * wsat
    denom = p * (1 + b * wsat / ((t + t_zero) ** 2))

    return numer / denom


def read_shapefile(shapefile):
    """Return the geometry (list) of each shape and the projection

    Parameters
    ----------
    shapefile: string
        path to a shapefile

    Returns
    -------
    geometries : list of array (number of points, 2)
    projection : projection string (PROJ4.strings)
    """
    with fiona.open(shapefile) as shapef:
        geometries = [
            np.array(shape["geometry"]["coordinates"]).squeeze() for shape in shapef
        ]
        projection = crs.to_string(shapef.crs)

    return geometries, projection


def mask_polygon(polygon, array, lon, lat):
    """Mask value outside polygon, work only for regular grid

    Parameters
    ----------
    polygon: array (X x 2)
        coordinates of the polygon
    array: array (nlon x nlat)
        array to mask
    lon: array (nlon x nlat)
        regular coordinates of the array
    lat: array (nlon x nlat)
        regular coordinates of the array

    Returns
    -------
    masked array: MaskedArray (nlon x nlat)


    Examples
    --------

    >>> import numpy as np
    >>> polygon = np.array([[5, 5],
    >>>                     [10, 5],
    >>>                     [10, 10],
    >>>                     [5, 10],
    >>>                     [5, 5]])
    >>> array = np.zeros((20, 20))
    >>> lon = np.arange(0, 20, dtype=float)
    >>> lat = np.arange(0, 20, dtype=float)
    >>> lons, lats = np.meshgrid(lon, lat)
    >>> clipped = mask_polygon(polygon, array, lons, lats)

    """

    xsize = array.shape[-1]
    ysize = array.shape[-2]

    lon0 = lon[0, 0]
    lat0 = lat[0, 0]
    lon1 = lon[-1, -1]
    lat1 = lat[-1, -1]
    nx, ny = lon.T.shape
    dx = (lon1 - lon0) / nx
    dy = (lat1 - lat0) / ny

    def return_indices(x, y):
        i = int((x - lon0) / dx)
        j = int((y - lat0) / dy)
        return i, j

    # transform the polygon coordinates into pixels coordinates
    pixels = [return_indices(x, y) for x, y in polygon]

    # rasterize of the polygon
    raster_poly = Image.new("L", (xsize, ysize), 1)
    rasterize = ImageDraw.Draw(raster_poly)
    rasterize.polygon(pixels, 0)
    mask = sp.frombuffer(raster_poly.tobytes(), "b")
    mask.shape = raster_poly.im.size[1], raster_poly.im.size[0]

    # clip the raster

    if array.ndim >= 3:
        mask.shape = (1,) + mask.shape

    mask = np.broadcast_arrays(mask, array)[0]

    clipped = ma.masked_where(mask != 0, array)

    return clipped


def potential_temperature(t, p, p0=1000):
    """Calculate the potential temperature

    Parameters
    ----------
    t: array
        temperature in K
    p: array
        pressure in hPa
    p0: float, optional
        pressure at sea level in hPa

    Returns
    -------
    potential temperature: array
    """
    return t * (p0 / p) ** (2 / 7.0)


def virtual_temperature(t, qv):
    """
    Based on http://glossary.ametsoc.org/wiki/Virtual_potential_temperature
    """
    return t * (1 + 0.61 * qv)


def brunt_vaisala(theta, z):
    """
    Based on http://glossary.ametsoc.org/wiki/Brunt-väisälä_frequency
    """
    dthdz = np.gradient(theta)[0] / np.gradient(z)[0]
    return np.sqrt((g / theta) * dthdz)


def equivalent_pot_temp(t, p, qv):
    """Return the equivalent potential temperature in K

    computation of equivalent potential temperature
    according to Bolton (1980) except constant 0.1998 (4.805
    according to Bolton).

    Parameters
    ----------
    t : np.array, np.float
        Temperature in K
    p : np.array, np.float
        Pressure in hPa
    qv : np.array, np.float
        Specific humidity in kg/kg

    Returns
    -------
    THE (equivalent potential temperature in K) : (np.array, np.float)

    Examples
    --------

    >>> t = 16 + 273
    >>> qv = 0.01156
    >>> p = 850
    >>> equivalent_pot_temp(t, p, qv)
    337.59612858187029
    """
    r = 1 / ((1 / qv) - 1)
    e = 100 * p * qv / (0.622 + 0.378 * qv)
    tl = 2840 / (3.5 * np.log(t) - np.log(e) - 0.1998) + 55
    the = t * (1000 / p) ** (0.2854 * (1 - 0.28 * r))
    the *= np.exp(((3.376 / tl) - 0.00254) * 1e3 * r * (1 + 0.81 * r))
    return the


def path_along_domain_boundary(lon, lat, nbnd=0, bnd_nan=False, fld=None):
    """Construct a path parallel to the domain boundary.

    Slice the lon and lat arrays. The distance to the boundary is given in
    grid points as <nbnd>.

    Returns
    -------
     x coordinates of path: array
     y coordinates of path: array
    """
    if not bnd_nan:
        nx, ny = lon.shape
        xs, xe, ys, ye = 0, nx, 0, ny
    else:
        if fld is None:
            raise ValueError("must pass fld for bnd_nan=T")
        xs, xe, ys, ye = locate_domain_in_nans(fld)

    xs += nbnd
    xe -= nbnd
    ys += nbnd
    ye -= nbnd

    pxs = (
        lon[xs:xe, ys].tolist()
        + lon[xe - 1, ys:ye].tolist()  # west
        + lon[xe:xs:-1, ye - 1].tolist()  # north
        + lon[xs, ye:ys:-1].tolist()  # east  # south
    )
    pys = (
        lat[xs:xe, ys].tolist()
        + lat[xe - 1, ys:ye].tolist()  # west
        + lat[xe:xs:-1, ye - 1].tolist()  # north
        + lat[xs, ye:ys:-1].tolist()  # east  # south
    )

    if (pxs[-1], pys[-1]) != (pxs[0], pys[0]):
        pxs.append(pxs[0])
        pys.append(pys[0])

    return pxs, pys


def locate_domain_in_nans(fld):
    """Measure the width of a rectangular NaN boundary zone around a field."""

    nx, ny = fld.shape
    iref, jref = int(nx / 2), int(ny / 2)

    # Left
    for i in range(nx):
        if not np.isnan(fld[i, jref]):
            bxs = i
            break
    else:
        raise Exception("all NaN at y={}".format(jref))

    # Right
    for i in range(nx - 1, -1, -1):
        if not np.isnan(fld[i, jref]):
            bxe = i + 1
            break

    # Bottom
    for j in range(ny):
        if not np.isnan(fld[iref, j]):
            bys = j
            break
    else:
        raise Exception("all NaN at x={}".format(iref))

    # Top
    for j in range(ny - 1, -1, -1):
        if not np.isnan(fld[iref, j]):
            bye = j + 1
            break

    return bxs, bxe, bys, bye


def return_erainterim_pressure(ps):
    """Return the pressure in hPa

    Parameters
    ----------
    ps: np.ndarray (2D)
        surface pressure in hPa

    Returns
    -------
    pressure in hPa (60 x ps.shape): np.ndarray
    """
    aklay = np.array(
        [
            [
                [
                    0.00000000e00,
                    3.68387103e-02,
                    3.66284877e-01,
                    1.38141561e00,
                    3.38863707e00,
                    6.61347675e00,
                    1.12063723e01,
                    1.72484627e01,
                    2.47573814e01,
                    3.36930466e01,
                    4.39634552e01,
                    5.54304657e01,
                    6.79156036e01,
                    8.12057953e01,
                    9.50592804e01,
                    1.09211288e02,
                    1.23379883e02,
                    1.37271729e02,
                    1.50587860e02,
                    1.63029587e02,
                    1.74304138e02,
                    1.84130539e02,
                    1.92245392e02,
                    1.98408646e02,
                    2.02409409e02,
                    2.04071716e02,
                    2.03260330e02,
                    1.99886566e02,
                    1.93914017e02,
                    1.85364380e02,
                    1.74323273e02,
                    1.60996384e02,
                    1.45775635e02,
                    1.29263840e02,
                    1.12267853e02,
                    9.57058945e01,
                    8.03584366e01,
                    6.66232605e01,
                    5.46236343e01,
                    4.43349915e01,
                    3.57835655e01,
                    2.88815498e01,
                    2.33108120e01,
                    1.88145714e01,
                    1.51855755e01,
                    1.22565460e01,
                    9.89247513e00,
                    7.98439217e00,
                    6.44434547e00,
                    5.20134640e00,
                    4.19295025e00,
                    3.36233878e00,
                    2.66637444e00,
                    2.07681704e00,
                    1.57533824e00,
                    1.15060127e00,
                    7.96423793e-01,
                    5.10365665e-01,
                    2.92126685e-01,
                    9.99999940e-02,
                ]
            ]
        ],
        dtype=sp.float32,
    ).T
    bklay = np.array(
        [
            [
                [
                    9.98815060e-01,
                    9.95824814e-01,
                    9.91144776e-01,
                    9.83966410e-01,
                    9.73653972e-01,
                    9.59733367e-01,
                    9.41880941e-01,
                    9.19912100e-01,
                    8.93770397e-01,
                    8.63515913e-01,
                    8.29314172e-01,
                    7.91424990e-01,
                    7.50191212e-01,
                    7.06027210e-01,
                    6.59408033e-01,
                    6.10857964e-01,
                    5.60939193e-01,
                    5.10240734e-01,
                    4.59367245e-01,
                    4.08927530e-01,
                    3.59523505e-01,
                    3.11738938e-01,
                    2.66128063e-01,
                    2.23204523e-01,
                    1.83430105e-01,
                    1.47203416e-01,
                    1.14848614e-01,
                    8.66042674e-02,
                    6.26121163e-02,
                    4.29057851e-02,
                    2.73995195e-02,
                    1.59103926e-02,
                    8.11201334e-03,
                    3.44813662e-03,
                    1.13827549e-03,
                    2.68609205e-04,
                    3.79117482e-05,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                    0.00000000e00,
                ]
            ]
        ],
        dtype=sp.float32,
    ).T
    ps = ps.reshape((1,) + ps.shape)
    return aklay + ps * bklay
