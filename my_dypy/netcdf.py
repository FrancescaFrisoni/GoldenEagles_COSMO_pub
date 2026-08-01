# coding:utf-8

"""
Module to read and write netcdf file.

Interface to netCDF4

"""

import sys
from collections import OrderedDict

import netCDF4
import xarray as xr
import numpy as np
import os


def _check_var(nvariables, variables):
    """ Check if <variables> are contained in <nvariables>
        return True/False, difference
    """
    nv = set(nvariables)
    v = set(variables)
    ok = set(nv).issuperset(v)
    return ok, v.difference(nv)


def read_var(filename, variables, index=slice(None), **kwargs):
    """ Extract variables from a netCDF

        Raise an IOerror if the file[s] is not found;
        netcdf time are return as datetime array if possible

    Parameters
    ----------
    filename: string
        path to a netCDF file
    variables: list or string
        variables to read as a list of string;
        a single variable can also be given as a string
    index: slice, optional
        A slice object; for example np.s_[0, :, 10]
    kwargs: keyword arguments
        A list of arguments to pass to the MFDataset class.
        For example to aggregate the files on the <time> dimension
        use aggdim='time'

    Returns
    -------
    list of numpy array: list
    """
    filename = [filename] if type(filename) is not list else filename
    if type(variables) is not list:
        variables = [variables]

    # test if filename exists
    for f in filename:
        if "http://" in f:
            continue
        if not os.path.isfile(f):
            raise IOError("{} was not found".format(f))

    try:
        #with xr.open_mfdataset(filename, **kwargs) as ncfile:
        with netCDF4.Dataset(filename, **kwargs) as ncfile:
            vararray = _extract_from_netcdf(ncfile, variables, index)
    except OSError as err:
        if len(filename) == 1:
            with netCDF4.Dataset(filename[0], **kwargs) as ncfile:
                vararray = _extract_from_netcdf(ncfile, variables, index)
        else:
            raise err
    return vararray


def _extract_from_netcdf(ncfile, variables, index):
    vararray = []
    nvariables = list(ncfile.variables.keys())
    ok, diffvar = _check_var(nvariables, variables)
    if not ok:
        err = "{} not found in file\n. Available :{}".format(
            ",".join(diffvar), ",".join(nvariables)
        )
        raise Exception(err)

    for var in variables:
        if var == "time":
            time = ncfile.variables[var]
            try:
                vardata = netCDF4.num2date(time[:], units=time.units)
            except (AttributeError, ValueError):
                vardata = time[:]
        else:
            ndim = ncfile.variables[var].ndim
            if   (ndim == 2) and (index != slice(None)):
                nindex = index[3:]
            elif (ndim == 3) and (index != slice(None)):
                nindex = index[2:]
            elif (ndim == 4) and (index != slice(None)):
                nindex = index[1:]
            else:
                nindex = index
            vardata = ncfile.variables[var][nindex]
        vararray.append(vardata.squeeze())
    return vararray


def read_var_bbox(filename, variables, bbox, lon="lon", lat="lat", return_index=False):
    """Read var only in the bbox

    Similar to read_var but read only the data in the given bounding box.

    Parameters
    ----------
    filename: list or string
        filename(s) where the data is located
    variables: list or string
        variable(s) name of the data to read from file
    bbox: list
        coordinates of the bounding box as:
        minlon, maxlon, minlat, maxlat
    lon: string or np.ndarray, default lon
        name of the 2D longitude array in the data;
        or 2D longitude array
    lat: string or np.ndarray, default lat
        name of the 2D latitude array in the data:
        or 2D latitude array
    return_index: boolean, default False
        If True return the index used to reduce the variable data, with the
         dimension (time, height, lon, lat)

    Returns
    -------
    bbox_lon: numpy array
        lon restricted to the bbox
    bbox_lat: numpy array
        lat restricted to the bbox
    bbox_dta: numpy arra
        data restricted to the bbox

    Adaption of LJA: Enable reading of COSMO-1E files
    Necessary was: Extension of index dimension --> ndim == 4 exception in read_var

    """
    if type(lon) is str:
        lon, lat = read_var(filename, [lon, lat])
    xll, xur, yll, yur = bbox
    xindex, yindex = np.where((lon > xll) & (lon < xur) & (lat > yll) & (lat < yur))
    xmin, xmax = xindex.min(), xindex.max() + 1
    ymin, ymax = yindex.min(), yindex.max() + 1
    index = np.s_[:,:,:, xmin:xmax, ymin:ymax]
    vardata = read_var(filename, variables, index=index)
    returndata = [lon[xmin:xmax, ymin:ymax], lat[xmin:xmax, ymin:ymax]]
    returndata.extend(vardata)
    if return_index:
        returndata.extend([index])
    return returndata

def read_var_rbbox(filename, variables, rbbox, lon="x_1", lat="y_1", return_index=False):
    """Read var only in the bbox

    Similar to read_var_bbox but takes rotated coordinates as input

    Parameters
    ----------
    filename: list or string
        filename(s) where the data is located
    variables: list or string
        variable(s) name of the data to read from file
    rbbox: list
        rotated coordinates of the bounding box as:
        minlon, maxlon, minlat, maxlat
    lon: string or np.ndarray, default lon
        name of the 2D longitude array in the data;
        or 2D longitude array
    lat: string or np.ndarray, default lat
        name of the 2D latitude array in the data:
        or 2D latitude array
    return_index: boolean, default False
        If True return the index used to reduce the variable data, with the
         dimension (time, height, lon, lat)

    Returns
    -------
    bbox_lon: numpy array
        lon restricted to the bbox
    bbox_lat: numpy array
        lat restricted to the bbox
    bbox_dta: numpy arra
        data restricted to the bbox

    """
    if type(lon) is str:
        lon, lat = read_var(filename, [lon, lat])
    xll, xur, yll, yur = rbbox
    xmin = np.where(lon >=  xll)[0].min()
    xmax = np.where(lon <=  xur)[0].max()
    ymin = np.where(lat >=  yll)[0].min()
    ymax = np.where(lat <=  yur)[0].max()
    index = np.s_[:,:,:, ymin:ymax+1, xmin:xmax+1]
    vardata = read_var(filename, variables, index=index)
    returndata = [lon[xmin:xmax+1], lat[ymin:ymax+1]]
    returndata.extend(vardata)
    if return_index:
        returndata.extend([index])
    return returndata

def read_gattributes(filename):
    """ Read global attributes from a netCDF

    Parameters
    ----------
    filename: string
        path to a netcdf file

    Returns
    -------
    global attributes: dictionary
    """
    with netCDF4.Dataset(filename) as ncfile:
        gattributes = ncfile.__dict__
    return gattributes


def read_dimensions(filename):
    """ Read dimensions

    Unlimited dimensions size are return as None

    Parameters
    ----------
    filename: string

    Returns
    -------
    a dictionary with name as key and the dimension as value: dictionary
    """
    dimensions = OrderedDict()
    with netCDF4.Dataset(filename) as ncfile:
        for name, dim in ncfile.dimensions.items():
            dimensions[name] = None if dim.isunlimited() else dim.size
    return dimensions


def read_variables(filename):
    """ Read variables

    Parameters
    ----------
    filename: string

    Returns
    -------
    a list of variables names: list
    """
    with netCDF4.Dataset(filename) as ncfile:
        variables = list(ncfile.variables.keys())
    return variables


def read_var_attributes(filename, var):
    """ Return the attributes of the variables as a dictionary"""
    with netCDF4.Dataset(filename) as ncfile:
        return OrderedDict(
            (n, getattr(ncfile.variables[var], n))
            for n in ncfile.variables[var].ncattrs()
        )


def create(outname, dimensions, gattributes, format="NETCDF3_CLASSIC"):
    """ create a netCDF """
    ncfile = netCDF4.Dataset(outname, "w", format=format)
    for dimname, size in dimensions.items():
        ncfile.createDimension(dimname, size)
    ncfile.setncatts(gattributes)
    return ncfile


def addvar(ncfile, varname, vardata, dimensions, attributes=None):
    """ Add a variable to an opened netcdf file"""
    var = ncfile.createVariable(varname, vardata.dtype, dimensions)
    if type(attributes) in [dict, OrderedDict]:
        var.setncatts(attributes)
    try:
        var[:] = vardata
    except Exception as e:
        sys.stderr.write("{0} in addvar()\n".format(e))
        sys.exit(1)


def addvar_to_file(filename, varname, vardata, dimensions, attributes=None):
    """Add a variable to a netcdf file"""
    with netCDF4.Dataset(filename, "r+") as ncfile:
        addvar(ncfile, varname, vardata, dimensions, attributes)
