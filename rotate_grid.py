# -*- coding: utf-8 -*-

import numpy as np
import cartopy

def rotate_points(pole_longitude,pole_latitude,lon,lat,direction,hemisphere):

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
    hemisphere: string
        NH: Northern hemisphere
        SH: Southern hemisphere

    Returns
    -------
    rlon: array
    rlat: array
    """
    if hemisphere == 'NH':
        lon = np.array(lon)
        lat = np.array(lat)
        rotatedgrid = cartopy.crs.RotatedPole(
            pole_longitude=pole_longitude,
            pole_latitude=pole_latitude
        )
        standard_grid = cartopy.crs.Geodetic()

        if direction == 'n2r':
            rotated_points = rotatedgrid.transform_points(standard_grid, lon, lat)
        elif direction == 'r2n':
            rotated_points = standard_grid.transform_points(rotatedgrid, lon, lat)

        rlon, rlat, _ = rotated_points.T

    if hemisphere == 'SH':
        lon = (lon*pi)/180; # Convert degrees to radians
        lat = (lat*pi)/180;

        SP_lon = pole_longitude - 180; # Calculate SP coordinates out of NP coordinates
        SP_lat = -pole_latitude ;

        theta = 90+SP_lat;  # Rotation around y-axis
        phi = SP_lon;       # Rotation around z-axis

        theta = (theta*pi)/180;
        phi = (phi*pi)/180;    # Convert degrees to radians

        x = cos(lon)*cos(lat); # Convert from spherical to cartesian coordinates
        y = sin(lon)*cos(lat);
        z = sin(lat);

        if direction == 'n2r':    # Nonrotated -> Rotated

            x_new = cos(theta)*cos(phi)*x + cos(theta)*sin(phi)*y + sin(theta)*z;
            y_new = -sin(phi)*x + cos(phi)*y;
            z_new = -sin(theta)*cos(phi)*x - sin(theta)*sin(phi)*y + cos(theta)*z;

        elif direction == 'r2n':  # Rotated -> Nonrotated

            phi = -phi;
            theta = -theta;

            x_new = cos(theta)*cos(phi)*x + sin(phi)*y + sin(theta)*cos(phi)*z;
            y_new = -cos(theta)*sin(phi)*x + cos(phi)*y - sin(theta)*sin(phi)*z;
            z_new = -sin(theta)*x + cos(theta)*z;

        lon_new = atan2(y_new,x_new); # Convert cartesian back to spherical coordinates
        lat_new = asin(z_new);

        rlon = (lon_new*180)/pi;      # Convert radians back to degrees
        rlat = (lat_new*180)/pi;

        if rlon > 0:
            rlon = rlon - 180

        elif rlon < 0:
            rlon = rlon + 180
   
    return rlon, rlat

#def rotate_points(pole_longitude,pole_latitude,lon,lat,direction):
#
#    #constants
#    zrpi18  = 57.2957795
#    zpir18  = 0.0174532925
#
#    #general settings
#    zsinpol = np.sin(zpir18*pole_latitude)
#    zcospol = np.cos(zpir18*pole_latitude)
#    zlampol = zpir18*pole_longitude
#    zphis    = zpir18*lat
#    zlams    = lon
#    zlams[zlams > 180] = zlams[zlams > 180] - 360
#    zlams    = zpir18*zlams
#
#    if direction == 'r2n':
#
#        #transform latitude
#        zarg    = zcospol*np.cos(zphis)*np.cos(zlams) + zsinpol*np.sin(zphis)
#        phstoph = zrpi18*np.arcsin(zarg)
#
#        #transform longitude
#        zarg1   = np.sin(zlampol)*(-zsinpol*np.cos(zlams)*np.cos(zphis)  + 
#                  zcospol*np.sin(zphis)) - np.cos(zlampol)*np.sin(zlams)*np.cos(zphis)
#        zarg2   = np.cos(zlampol)*(-zsinpol*np.cos(zlams)*np.cos(zphis)  + 
#                  zcospol*np.sin(zphis)) + np.sin(zlampol)*np.sin(zlams)*np.cos(zphis)
#        lmstolm = zrpi18*np.arctan2(zarg1, zarg2)
#
#        return lmstolm, phstoph
#
#
#    elif direction == 'n2r':
#
#        #transform latitude
#        zarg    = zcospol*np.cos(lat)*np.cos(lon-zlampol) + zsinpol*np.sin(lat)
#        phtophs = zrpi18*np.arcsin(zarg)
#
#        #transform longitude
#        zarg1   = - np.sin(lon-zlampol)*np.cos(lat)
#        zarg2   = - zsinpol*np.cos(lat)*np.cos(lon-zlampol)+zcospol*np.sin(lat)
#        lmtolms = zrpi18*np.arctan2(zarg1,zarg2)
#
#        return lmtolms, phtophs


#---------------------------------------------------------
# Converstion of CH coords to WGS latlon
# Source:
# https://github.com/ValentinMinder/Swisstopo-WGS84-LV03/blob/master/scripts/py/wgs84_ch1903.py
#---------------------------------------------------------
# Convert CH y/x to WGS lat
def CHtoWGSlat(y, x):
    # Axiliary values (% Bern)
    y_aux = (y - 2600000) / 1000000
    x_aux = (x - 1200000) / 1000000
    lat = (16.9023892 + (3.238272 * x_aux)) + \
            - (0.270978 * pow(y_aux, 2)) + \
            - (0.002528 * pow(x_aux, 2)) + \
            - (0.0447 * pow(y_aux, 2) * x_aux) + \
            - (0.0140 * pow(x_aux, 3))
    # Unit 10000" to 1" and convert seconds to degrees (dec)
    lat = (lat * 100) / 36
    return lat

# Convert CH y/x to WGS long
def CHtoWGSlng(y, x):
    # Axiliary values (% Bern)
    y_aux = (y - 2600000) / 1000000
    x_aux = (x - 1200000) / 1000000
    lng = (2.6779094 + (4.728982 * y_aux) + \
            + (0.791484 * y_aux * x_aux) + \
            + (0.1306 * y_aux * pow(x_aux, 2))) + \
            - (0.0436 * pow(y_aux, 3))
    # Unit 10000" to 1" and convert seconds to degrees (dec)
    lng = (lng * 100) / 36
    return lng

