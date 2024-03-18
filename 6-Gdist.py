import matplotlib.pyplot as plt
import datetime as dt
import sys
sys.path.extend(["py/","pyrt/"])
from fetch2D import *
from fetch import FetchModel
import util
import mplstyle
import copy

from scipy.io import readsav
import cartopy.crs as ccrs
import cartopy
import matplotlib.ticker as mticker
import numpy as np
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.dates as mdates
import os

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import matplotlib as mpl
mpl.rcParams.update({"xtick.labelsize": 5, "ytick.labelsize":5, "font.size":5})

params = [
    {
        "name": "e",
        "unit_multiplier": 1.,
        "unit": "cc",
        "mol_to_cc": "T",
        "interpolate": {"scale": "log", "type": "linear"},
    },
]
date = dt.datetime(2017,8,21,18,30)

dat = readsav("dataset/waccmx_20170821eclipse_dne.sav")


mpl.rcParams.update({"xtick.labelsize":6, "ytick.labelsize":6, "font.size":6})
fig = plt.figure(figsize=(4, 2), dpi=240)
ax = fig.add_subplot(
            111,
            projection=ccrs.PlateCarree()
        )
kwargs = {}
kwargs["edgecolor"] = "k"
kwargs["facecolor"] = "none"
kwargs["linewidth"] = 0.3
feature = cartopy.feature.NaturalEarthFeature(
            "physical", "coastline", "50m", **kwargs
        )
ax.add_feature(feature, **kwargs)
ax.set_extent([-140, -40, 20, 70], crs=ccrs.PlateCarree())
plt_lons = np.arange(-180, 181, 30)
mark_lons = np.arange(-180, 180, 30)
plt_lats = np.arange(20, 90, 15)
gl = ax.gridlines(ccrs.PlateCarree(), linewidth=0.5, draw_labels=True,)
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlocator = mticker.FixedLocator(plt_lons)
gl.ylocator = mticker.FixedLocator(plt_lats)
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.n_steps = 90
gl.ylabel_style = {'size':8, 'color': 'b'}
gl.xlabel_style = {"size":8, 'color': 'b'}

glat, glon = np.meshgrid(dat["glat"], dat["glon"])
dne = np.copy(dat["dne"])
dne[dne<-1e4] = np.nan

import matplotlib.colors as colors
XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), glon, glat)
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
im = ax.contourf(
    X, Y, dne.T, transform=ccrs.PlateCarree(), cmap="jet_r", alpha=0.9,
    norm=colors.Normalize(vmin=-1e4, vmax=4e4), **kwargs
)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
cb = fig.colorbar(im, ax=ax, cax=cax)
cb.set_label(r"$GC^p$, /cc")
#cb.set_lim(-1e4, 4e4)


#cb = fig.colorbar(im, ax=ax, cax=cax)
#cb.set_label(r"Occultation, $\mathcal{O}$")
ax.text(0.05, 1.05, "21 August 2017, 18:30 UT", transform=ax.transAxes, ha="left", va="center")
ax.text(0.95, 1.05, "Coord: Geo", transform=ax.transAxes, ha="right", va="center")


fig.subplots_adjust(wspace=0.5, hspace=0.3)
fig.savefig("dataset/figures/globe_distribution.png", bbox_inches="tight")