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
lat_boulder, lon_boulder = 40.0150, -105.2705
lat_lusk, lon_lusk = 40.0150, -105.2705
XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), 
                                          np.array([[lonx, lonx, lonx]]), 
                                          np.array([[center, lat_low, lat_up]]))
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
ax.scatter(X, Y, s=10, color="magenta", marker=".", fc="None", lw=0.3)

XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), 
                                          np.array([[lon_boulder, lon_lusk]]), 
                                          np.array([[lat_boulder, lat_lusk]]))
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
ax.scatter(X, Y, s=1, color="b", marker="D", fc="None", lw=0.3)

fig.subplots_adjust(wspace=0.5, hspace=0.3)
fig.savefig("dataset/figures/zenith_angle.png", bbox_inches="tight")

cos_szas = [
    np.cos(np.deg2rad(90 - pysolar.solar.get_altitude(center, lonx, date))),
    np.cos(np.deg2rad(90 - pysolar.solar.get_altitude(lat_low, lonx, date))),
    np.cos(np.deg2rad(90 - pysolar.solar.get_altitude(lat_up, lonx, date))),
]

import sys
sys.path.extend(["py/","pyrt/"])
from fetch import FetchModel, pconst, DiffWACCMX
import aacgmv2


mpl.rcParams.update({"xtick.labelsize": 5, "ytick.labelsize":5, "font.size":5})

mpl.rcParams.update({"xtick.labelsize":10, "ytick.labelsize":10, "font.size":10})
# fig = plt.figure()

k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
stns = ["lat_low", "center", "lat_up"]
params = [
    {
        "name": "e",
        "unit_multiplier": 1.,
        "unit": "cc",
        "mol_to_cc": "T",
        "interpolate": {"scale": "log", "type": "linear"},
    },
]
fm = FetchModel(
    "dataset/29jan2022/ECL_1100_0100_all.nc",
    params=params,
)
fig, axes = plt.subplots(dpi=240, figsize=(6, 5), nrows=2, ncols=1, sharex="all")
for ax in axes:
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.set_xlim([
        dt.datetime(2017,8,21,15),
        dt.datetime(2017,8,21,21)
    ])

axes[0].set_ylim([2, 8])
axes[0].set_ylabel(
    r"$foF_2$, MHz"
)
axes[1].set_ylim([2, 5])
axes[1].set_ylabel(
    r"$foF_1$, MHz"
)
axes[1].set_xlabel("UT")
values_150, values_240 = [], []
for stn, c in zip(stns, ["darkred", "darkblue", "darkgreen"]):
    fm.run_height_interpolate(stn)
    loc = fm.__fetch_latlon_by_station__(stn)
    olat, olon, _ = aacgmv2.get_aacgm_coord(loc["lat"], loc["lon"], 300, dt.datetime(2017,8,21,18))

    ttime0 = fm.dataset[0]["interpol_value"]["time"]
    ne = fm.dataset[0]["interpol_value"]["value"]
    f_240, f_150 = (
        np.sqrt(ne[:, np.argmin(np.abs(fm.intp_height-240))]*1e6*k)/(2*np.pi),
        np.sqrt(ne[:, np.argmin(np.abs(fm.intp_height-150))]*1e6*k)/(2*np.pi)
    )
    axes[0].plot(ttime0, f_240/1e6, color=c, ls="-", lw=0.8, label=r"$\theta=%.1f$"%olat)
    axes[1].plot(ttime0, f_150/1e6, color=c, ls="-", lw=0.8,)
    axes[0].legend(loc=1)
    values_150.append(f_150[0])
    values_240.append(f_240[0])
print((np.max(values_150)-np.min(values_150))/np.min(values_150))
print((np.max(values_240)-np.min(values_240))/np.min(values_240))
fig.subplots_adjust(hspace=0.1, wspace=0.1)
fig.savefig("dataset/figures/latitudinal_D_compare.png", bbox_inches="tight")
del fm

fig, axes = plt.subplots(dpi=240, figsize=(6, 5), nrows=2, ncols=1, sharex="all")
for ax in axes:
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.set_xlim([
        dt.datetime(2017,8,21,15),
        dt.datetime(2017,8,21,21)
    ])
axes[0].set_ylabel(
    r"$\delta(p-l)$, $cm^{-3} s^{-1}$"
)
# axes[1].set_ylim([2, 5])
axes[1].set_ylabel(
    r"$\sum\delta(D_x)$, $cm^{-3} s^{-1}$"
)
axes[1].set_xlabel("UT")

params = [
    {
        "name": "Op_CHMP",
        "unit_multiplier": 1.0,
        "unit": "cc",
        "interpolate": {"scale": "linear", "type": "linear"},
    },
    {
        "name": "Op_CHML",
        "unit_multiplier": 1.0,
        "unit": "cc",
        "interpolate": {"scale": "linear", "type": "linear"},
    },
    {
        "name": "amb_diff",
        "unit_multiplier": 1.0,
        "unit": "cc",
        "interpolate": {"scale": "linear", "type": "linear"},
    },
    {
        "name": "dwind",
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

ts_pl, ts_dm = [], []

for stn, c in zip(stns, ["darkred", "darkblue", "darkgreen"]):
    dw = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=stn,
    )
    dct = dw.diffential_difference_data(
        kind="TS", hs=[150, 240]
    )
    loc = dw.eclipse.__fetch_latlon_by_station__(stn)
    olat, olon, _ = aacgmv2.get_aacgm_coord(loc["lat"], loc["lon"], 300, dt.datetime(2017,8,21,18))
    h = "240"
    axes[0].plot(
        dct["Op_CHML_"+h].time, 
        dct["Op_CHMP_"+h].d_dif-dct["Op_CHML_"+h].d_dif, 
        color=c, ls="-", lw=0.8, 
        label=r"$\theta=%.1f$"%olat
    )
    ts_pl.append(dct["Op_CHMP_"+h].d_dif-dct["Op_CHML_"+h].d_dif)
    axes[1].plot(
        dct["amb_diff_"+h].time, 
        dct["amb_diff_"+h].d_dif+dct["dwind_"+h].d_dif+dct["dfield_"+h].d_dif, 
        color=c, ls="-", lw=0.8,
    )
    ts_dm.append(dct["amb_diff_"+h].d_dif+dct["dwind_"+h].d_dif+dct["dfield_"+h].d_dif)
    del dw, dct
axes[0].legend(loc=1)

ts_pl, ts_dm = (np.array(ts_pl), np.array(ts_dm))

dts_pl, dts_dm = (
    np.abs(ts_pl[0,:]-ts_pl[0,2]),
    np.abs(ts_dm[0,:]-ts_dm[0,2])
)
print(
    np.sqrt(((ts_pl[0] - ts_pl[2]) ** 2).mean()),
    np.sqrt(((ts_dm[0] - ts_dm[2]) ** 2).mean())
)
axes[0].text(
    0.05, 0.9, "RMSD: %.1f"%np.sqrt(((ts_pl[0] - ts_pl[2]) ** 2).mean()),
    ha="left", va="center",
    transform=axes[0].transAxes
)
axes[1].text(
    0.05, 0.9, "RMSD: %.1f"%np.sqrt(((ts_dm[0] - ts_dm[2]) ** 2).mean()),
    ha="left", va="center",
    transform=axes[1].transAxes
)
fig.subplots_adjust(hspace=0.1, wspace=0.1)
fig.savefig("dataset/figures/latitudinal_Dx_compare.png", bbox_inches="tight")

