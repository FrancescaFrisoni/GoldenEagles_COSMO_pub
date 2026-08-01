# -*- coding:utf-8 -*-

import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.pyplot import subplot
from lagranto.plotting import CartoFigure, plot_trajs

__all__ = ["CartoFigure", "Mapfigure", "plot_trajs", "SkewT"]


class Mapfigure:
    """
        Class based on Basemap with additional functionality
        such as plot_trajectories
    """

    def __init__(
        self,
        resolution="i",
        projection="cyl",
        domain=None,
        lon=None,
        lat=None,
        basemap=None,
        **kwargs
    ):

        if basemap is None:
            if (domain is None) & (lon is not None):
                domain = [lon.min(), lon.max(), lat.min(), lat.max()]
            if projection == "eqc":
                projection = "cyl"
            if domain is not None:
                kwargs["llcrnrlon"] = domain[0]
                kwargs["urcrnrlon"] = domain[1]
                kwargs["llcrnrlat"] = domain[2]
                kwargs["urcrnrlat"] = domain[3]
            kwargs["resolution"] = resolution
            kwargs["projection"] = projection
            self.m = Basemap(**kwargs)
        else:
            self.m = basemap
        if lon is not None:
            self.x, self.y = self.m(lon, lat)

    def __getattr__(self, item):
        return getattr(self.m, item)

    def __call__(self, *args, **kwargs):
        return self.m(*args, **kwargs)

    def __dir__(self):
        return self.m.__dir__() + ["drawmap", "plot_traj"]

    def drawmap(
        self,
        continent=False,
        nbrem=5,
        nbrep=5,
        coastargs={},
        countryargs={},
        meridiansargs={},
        parallelsargs={},
    ):
        """
        draw basic features on the map
        nbrem: interval bewteen meridians
        nbrep: interval between parallels
        """
        self.drawcoastlines(**coastargs)
        self.drawcountries(**countryargs)
        merid = np.arange(-180, 180, nbrem)
        parall = np.arange(-90, 90, nbrep)
        self.drawmeridians(merid, labels=[0, 0, 0, 1], **meridiansargs)
        self.drawparallels(parall, labels=[1, 0, 0, 0], **parallelsargs)
        if continent:
            self.fillcontinents(color="lightgrey")

    def plot_traj(self, trajs, variable, cmap="Spectral", levels=None, **kwargs):
        """Plot trajectories on a map

        Usage:
            m = Mapfigure(domain=[5, 15, 40, 50])
            m.drawmap()
            m.plot_traj(trajs, 'QV')
        """
        if self.ax is None:
            self.ax = subplot()
        lon, lat = trajs["lon"].copy(), trajs["lat"].copy()
        lonnanvalues = np.isnan(lon)
        latnanvalues = np.isnan(lat)
        x, y = self.m(trajs["lon"], trajs["lat"])
        # m() does not account for nan values
        x[lonnanvalues] = np.nan
        y[latnanvalues] = np.nan
        trajs["lon"] = x
        trajs["lat"] = y
        lc = plot_trajs(self.ax, trajs, variable, cmap=cmap, levels=levels, **kwargs)
        trajs["lon"] = lon
        trajs["lat"] = lat
        return lc


class SkewT:
    """
    Create a skewT-logP diagramm
    """

    # Private attributes
    SKEWNESS = 37.5
    T_ZERO = 273.15
    # P_bot is used to define the skewness of the plot
    P_BOT = 100000
    L = 2.501e6  # latent heat of vaporization
    R = 287.04  # gas constant air
    RV = 461.5  # gas constant vapor
    EPS = R / RV

    CP = 1005.0
    CV = 718.0
    KAPPA = (CP - CV) / CP
    G = 9.81

    # constants used to calculate moist adiabatic lapse rate
    # See formula 3.16 in Rogers&Yau
    A = 2.0 / 7.0
    B = EPS * L * L / (R * CP)
    C = A * L / R

    def __init__(self, ax, prange={"pbot": 1000.0, "ptop": 100.0, "dp": 1.0}):
        """ initalize a skewT instance """

        self.pbot = prange["pbot"] * 100.0
        self.ptop = prange["ptop"] * 100.0
        self.dp = prange["dp"] * 100.0
        self.ax = ax
        # Defines the ranges of the plot
        self.plevs = np.arange(self.pbot, self.ptop - 1, -self.dp)
        self._isotherms()
        self._isobars()
        self._dry_adiabats()
        self._moist_adiabats()
        # self._mixing_ratio()

    def _skewnessTerm(self, P):
        return self.SKEWNESS * np.log(self.P_BOT / P)

    def _isotherms(self):
        for temp in np.arange(-100, 50, 10):
            self.ax.semilogy(
                temp + self._skewnessTerm(self.plevs),
                self.plevs,
                basey=np.e,
                color="blue",
                linestyle=("solid" if temp == 0 else "dashed"),
                linewidth=0.5,
            )

    def _isobars(self):
        for n in np.arange(self.P_BOT, self.ptop - 1, -10 ** 4):
            self.ax.plot([-40, 50], [n, n], color="black", linewidth=0.5)

    def _mixing_ratio(self):
        rdv = 0.622
        B1 = 243.04  # °C
        C1 = 610.94  # Pa
        A1 = 17.625
        t = np.arange(-30, 50, 10)
        m = np.zeros((self.plevs.size, t.size))
        for i, temp in enumerate(t):
            es = C1 * np.exp(A1 * temp / (B1 + temp))
            m[:, i] = rdv * es / (self.plevs - es)
        t, p = np.meshgrid(t, self.plevs)
        self.ax.contour(t, p, m)

    def _dry_adiabats(self):
        for tk in self.T_ZERO + np.arange(-30, 210, 10):
            dry_adiabat = tk * (self.plevs / self.P_BOT) ** self.KAPPA - (
                self.T_ZERO + self._skewnessTerm(self.plevs)
            )
            self.ax.semilogy(
                dry_adiabat,
                self.plevs,
                basey=np.e,
                color="brown",
                linestyle="dashed",
                linewidth=0.5,
            )

    def _moist_adiabats(self):
        ps = [p for p in self.plevs if p <= self.P_BOT]
        tlevels = np.concatenate(
            (np.arange(-40.0, 10.1, 5.0), np.arange(12.5, 45.1, 2.5))
        )
        for temp in tlevels:
            moist_adiabat = []
            for p in ps:
                temp -= self.dp * self.gamma_s(temp, p)
                moist_adiabat.append(temp + self._skewnessTerm(p))
            self.ax.semilogy(
                moist_adiabat,
                ps,
                basey=np.e,
                color="green",
                linestyle="dashed",
                linewidth=0.5,
            )

    def plot_data(self, p, T, color="black", style="solid"):
        self.ax.semilogy(
            T + self._skewnessTerm(p * 100),
            p * 100,
            basey=np.e,
            color=(color),
            linestyle=(style),
            linewidth=1.5,
        )
        self.ax.axis([-40, 50, self.pbot, self.ptop])
        self.ax.set_xlabel(r"Temperature ($^{\circ}$ C)")
        xticks = np.arange(-40, 51, 5)
        self.ax.set_xticks(
            xticks, ["" if tick % 10 != 0 else str(tick) for tick in xticks]
        )
        self.ax.set_ylabel("Pressure (hPa)")
        yticks = np.arange(self.pbot, self.ptop - 1, -10 ** 4)
        self.ax.set_yticks(yticks)
        self.ax.set_yticklabels(["{:2.0f}".format(label) for label in yticks / 100])

    def plot_windsbarbs(self, p, u, v, offset=40):
        x = p.copy()
        x[:] = offset
        ax2 = self.ax.twinx()
        ax2.barbs(x[::2], p[::2] * 100, u[::2], v[::2])
        ax2.set_yscale("log", basey=np.e)
        yticks = np.arange(self.pbot, self.ptop - 1, -10 ** 4)
        ax2.yaxis.set_ticks(yticks)
        ax2.set_ylim([self.pbot, self.ptop])
        ax2.set_xlim([-40, 50])

    def es(self, T):
        """Returns saturation vapor pressure (Pascal) at temperature T (Celsius)
        Formula 2.17 in Rogers&Yau"""
        return 611.2 * np.exp(17.67 * T / (T + 243.5))

    def gamma_s(self, T, p):
        """Calculates moist adiabatic lapse rate for T (Celsius) and p (Pa)
        Note: We calculate dT/dp, not dT/dz
        See formula 3.16 in Rogers&Yau for dT/dz,
        but this must be combined with
        the dry adiabatic lapse rate (gamma = g/cp) and the
        inverse of the hydrostatic equation (dz/dp = -RT/pg)"""
        esat = self.es(T)
        wsat = self.EPS * esat / (p - esat)  # Rogers&Yau 2.18
        numer = self.A * (T + self.T_ZERO) + self.C * wsat
        denom = p * (1 + self.B * wsat / ((T + self.T_ZERO) ** 2))
        return numer / denom  # Rogers&Yau 3.16
