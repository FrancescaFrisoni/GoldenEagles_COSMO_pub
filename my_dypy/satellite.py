# encoding: utf8

# Copyright (c) 2012 Martin Raspaud <martin.raspaud@smhi.se>
# Copyright (c) 2015 Nicolas Piaget <nicolas.piaget@env.ethz.ch>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

""" Collections of tools to simplify the use of pytroll
"""
import mpop
from path import Path
from pyproj import Proj
from pyresample import geometry
from jinja2 import Template
from datetime import datetime


def coord2area(proj, min_lat, max_lat, min_lon, max_lon, resolution, fromstring):
    """Return the parameters need for an AreaDefinition (pyresample)
    """
    lat_0 = (min_lat + max_lat) / 2.0
    lon_0 = (min_lon + max_lon) / 2.0

    if fromstring:
        p = Proj(proj)

    else:
        p = Proj(proj=proj, lat_0=lat_0, lon_0=lon_0, ellps="WGS84")

    x_ur, y_ur = p(max_lon, max_lat)
    x_lr, y_lr = p(max_lon, min_lat)
    x_ll, y_ll = p(min_lon, min_lat)
    x_ul, y_ul = p(min_lon, max_lat)

    area_extent = (min(x_ll, x_ul), min(y_lr, y_ll), max(x_ur, x_lr), max(y_ul, y_ur))
    xsize = int((area_extent[2] - area_extent[0]) / resolution)
    ysize = int((area_extent[3] - area_extent[1]) / resolution)

    return lon_0, lat_0, xsize, ysize, area_extent


def coord2areadef(domain, resolution, proj="eqc", name="default", fromstring=False):
    """Create an AreaDefinition as defined by pyresample

        Args:
            domain (list):      [min_lon, max_lon, min_lat, max_lat]
            resolution (int):   resolution in m
            proj (string):      Projection name or epsg:<number> [default:eqc]
            name (string):      Name of the AreaDefinition
            fromstring (bool):  If True, the proj parameters is used as init
                                string for Proj(<proj>) else:
                                Proj(proj=<proj>, lat_0=lat_0, lon_0=0,
                                    ellps='WGS84')
                                where lat_0 and lon_0 are coordinates of
                                the middle point of <domain>

        return a geometry.AreaDefinition object

    """
    if proj == "cyl":
        proj = "eqc"

    min_lon, max_lon, min_lat, max_lat = domain
    lon_0, lat_0, xsize, ysize, area_extent = coord2area(
        proj, min_lat, max_lat, min_lon, max_lon, resolution, fromstring
    )

    proj_id = "{}_{}_{}".format(proj, lon_0, lat_0)
    if fromstring:
        proj_dict = dict([e.lstrip("+").split("=") for e in proj.split(" ") if e != ""])
    else:
        proj_dict = {"proj": proj, "lat_0": lat_0, "lon_0": lon_0, "ellps": "WGS84"}

    return geometry.AreaDefinition(
        name, name, proj_id, proj_dict, xsize, ysize, area_extent
    )


def prepare_config(
    casefolder,
    sats=[{"name": "Meteosat", "number": 9}],
    seviri="seviri",
    safnwc="safnwc",
    overwrite=True,
    templatedir="/net/litho/atmosdyn/dypy/config/",
):
    """
    Write configuration file in:
        <casefolder>/.config/
    and set the configuration path of mpop to this directory

    The data are excpected to be in <seviri> (relative) to the
    <casefolder> for seviri data and in <safnwc> for SAFNWC products.

    If <overwrite> is not set to False, the configuration template will be
    overwritten and the old ones save up to three time
    """
    casefolder = Path(casefolder)
    templatedir = Path(templatedir)
    if not casefolder.isdir():
        raise NotADirectoryError(str(casefolder))

    configdir = casefolder.joinpath(".config/")
    seviri_path = casefolder.joinpath(seviri)
    safnwc_path = casefolder.joinpath(safnwc)

    if overwrite:
        configdir.makedirs_p()

        # backup the last three configurations
        backup = configdir.joinpath("backup", datetime.now().isoformat())
        if len(configdir.files()) > 0:
            backup.makedirs_p()
            backup_dirs = backup.parent.listdir()
            if len(backup_dirs) > 3:
                sorted(backup_dirs)[0].rmtree_p()
            for f in configdir.files():
                f.move(backup)

        # copy a simple mpop file
        mpopfile = templatedir.joinpath("mpop.cfg")
        mpopfile.copyfile(configdir.joinpath("mpop.cfg"))

        # create an empty areas.def file
        areadef = configdir.joinpath("areas.def")
        areadef.touch()

        # prepare the config file for each satellite
        for sat in sats:
            filename = "{name}-{number}.cfg".format(
                name=sat["name"].capitalize(), number=sat["number"]
            )
            templatefile = templatedir.joinpath(filename)
            outfile = configdir.joinpath(filename)
            _config_from_template(
                templatefile, outfile, seviri_path=seviri_path, safnwc_path=safnwc_path
            )

    mpop.CONFIG_PATH = configdir


def _config_from_template(templatefile, outfile, **kwargs):
    """
    Fill a template and write it in the configdir

    """
    with open(templatefile, "r") as ifile:
        template = Template(ifile.read())
    with open(outfile, "w") as f:
        f.write(template.render(**kwargs))
