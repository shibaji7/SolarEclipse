import matplotlib.pyplot as plt
import datetime as dt

import cartopy.crs as ccrs
import cartopy
import matplotlib.ticker as mticker
import numpy as np
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.dates as mdates

import pysolar
import pytz

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import matplotlib as mpl
mpl.rcParams.update({"xtick.labelsize": 5, "ytick.labelsize":5, "font.size":5})

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

ax.text(0.05, 1.05, "21 August 2017, 18:30 UT", transform=ax.transAxes, ha="left", va="center")
ax.text(0.95, 1.05, "Coord: Geo", transform=ax.transAxes, ha="right", va="center")

date = dt.datetime(2017, 8, 27, 18, 30, 0, tzinfo=pytz.utc)
glats, glons = np.meshgrid(
    np.arange(20, 72), np.arange(-140, -38)
)
sza = np.zeros_like(glats)
for i in range(glats.shape[0]):
    for j in range(glats.shape[1]):
        sza[i,j] = 90 - pysolar.solar.get_altitude(glats[i,j], glons[i,j], date)

XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), glons, glats)
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]

im = ax.contourf(
    X, Y, np.cos(np.deg2rad(sza)), 
    transform=ccrs.PlateCarree(), 
    cmap="gray_r", alpha=0.5,
    levels=[0,.3,.6,.7,.85,1],
    vmin=0, vmax=1, **kwargs
)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
cb = fig.colorbar(im, ax=ax, cax=cax)
cb.set_label(r"$\cos(\chi)$, $\chi$ is SZA in deg ($^\circ$)")
ax.text(0.05, 0.95, "(a)", transform=ax.transAxes, ha="left", va="center")

center, lonx = 41, -99
lat_low, lat_up = center-7, center+7
XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), 
                                          np.array([[lonx, lonx, lonx]]), 
                                          np.array([[center, lat_low, lat_up]]))
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
ax.scatter(X, Y, s=10, color="magenta", marker=".", fc="None", lw=0.3)

fig.subplots_adjust(wspace=0.5, hspace=0.3)
fig.savefig("dataset/figures/zenith_angle.png", bbox_inches="tight")

cos_szas = [
    np.cos(np.deg2rad(90 - pysolar.solar.get_altitude(center, lonx, date))),
    np.cos(np.deg2rad(90 - pysolar.solar.get_altitude(lat_low, lonx, date))),
    np.cos(np.deg2rad(90 - pysolar.solar.get_altitude(lat_up, lonx, date))),
]

import sys
sys.path.extend(["py/","pyrt/"])
from fetch import DiffWACCMX
stn_low, stn_up = "lat_low", "lat_up"
params = [
    {
        "name": "amb_diff",
        "unit_multiplier": 1.0,
        "unit": "cc",
        "interpolate": {"scale": "linear", "type": "linear"},
    },
    {
        "name": "dfield",
        "unit_multiplier": 1.0,
        "unit": "cc",
        "interpolate": {"scale": "linear", "type": "linear"},
    },
]
dw_low,dw_up =  (
    DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=stn_low,
    ),
    DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=stn_up,
    )
)
dct_low, dct_up = (
    dw_low.diffential_difference_data(),
    dw_up.diffential_difference_data()
)
del dw_low, dw_up

print(dct_low.keys())
fig, axes = plt.subplots(dpi=200, figsize=(6, 5), nrows=2, ncols=2, sharex="all", sharey="all")
axes = axes.ravel()

for ax in axes:
    ax = axes[0]
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.set_xlim([
        dt.datetime(2017,8,21,15),
        dt.datetime(2017,8,21,21)
    ])
    ax.set_ylim(-3,3)

ax = axes[0]
ax.set_ylabel(
    r"$\delta(D_{\alpha})$, $cm^{-3}s^{-1}$"
)
ax.plot(
    dct_low["amb_diff_150"].time, 
    dct_low["amb_diff_150"].d_dif,
    color="r", lw=0.5, ls="-", 
    label=r"$\theta-8$"
)
ax.plot(
    dct_up["amb_diff_150"].time, 
    dct_up["amb_diff_150"].d_dif,
    color="k", lw=0.5, ls="-", 
    label=r"$\theta+8$"
)
ax.legend(loc=1)
ax.text(0.05,1.05,"h=150 km", ha="left", va="center", transform=ax.transAxes)
ax = axes[1]
ax.plot(
    dct_low["amb_diff_240"].time, 
    dct_low["amb_diff_240"].d_dif,
    color="r", lw=0.5, ls="-", 
)
ax.plot(
    dct_up["amb_diff_240"].time, 
    dct_up["amb_diff_240"].d_dif,
    color="k", lw=0.5, ls="-",
)
ax.text(0.05,1.05,"h=240 km", ha="left", va="center", transform=ax.transAxes)

ax = axes[2]
ax.set_ylabel(
    r"$\delta(D_{\vec{E}\times\vec{B}})$, $cm^{-3}s^{-1}$"
)
ax.plot(
    dct_low["dfield_150"].time, 
    dct_low["dfield_150"].d_dif,
    color="r", lw=0.5, ls="-", 
)
ax.plot(
    dct_up["dfield_150"].time, 
    dct_up["dfield_150"].d_dif,
    color="k", lw=0.5, ls="-", 
)
ax.set_xlabel("UT")
ax = axes[3]
ax.plot(
    dct_low["dfield_240"].time, 
    dct_low["dfield_240"].d_dif,
    color="r", lw=0.5, ls="-", 
)
ax.plot(
    dct_up["dfield_240"].time, 
    dct_up["dfield_240"].d_dif,
    color="k", lw=0.5, ls="-",
)
ax.set_xlabel("UT")


fig.subplots_adjust(hspace=0.1, wspace=0.1)
fig.savefig("dataset/figures/latitudinal_D_compare.png", bbox_inches="tight")