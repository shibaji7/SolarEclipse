# G-Condition Latitudeinal Dependence Study
import matplotlib.pyplot as plt
import datetime as dt
import sys
sys.path.append("py/")
import mplstyle
from fetch import *
import utils
import cartopy.crs as ccrs
import cartopy
import matplotlib.ticker as mticker
import numpy as np
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.dates as mdates

import aacgmv2

import netCDF4 as nc
center, lonx = 41, -99
lat_low, lat_up = center-7, center+7

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import matplotlib as mpl
mpl.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize":15, "font.size":15})

def fetch_occl_data(d=dt.datetime(2017,8,21,18)):
    folder = "dataset/Mark/euv/%s_%dkm_171_1.nc"
    file = folder % (d.strftime("%Y%m%d%H%M%S"), 150)
    
    ds = nc.Dataset(file)
    of = ds.variables["of"][:]
    lat, lon = ds.variables["glat"][:], ds.variables["glon"][:]
    glat, glon = np.meshgrid(lat, lon)
    return lat, lon, glat, glon, of

lat, lon, glat, glon, of = fetch_occl_data()

params = [
    {
        "name": "e",
        "unit_multiplier": 1.,
        "unit": "cc",
        "mol_to_cc": "T",
        "interpolate": {"scale": "log", "type": "linear"},
    },
]

dw0 = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=None,
        loc={"lat": center, "lon":lonx},
    )
k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
e_eclipse0 = dw0.eclipse.dataset[0]["interpol_value"]["value"]
ttime0 = dw0.eclipse.dataset[0]["interpol_value"]["time"]
f0_240, f0_150 = (
    np.sqrt(e_eclipse0[:, np.argmin(np.abs(dw0.eclipse.intp_height-240))]*1e6*k)/(2*np.pi), 
    np.sqrt(e_eclipse0[:, np.argmin(np.abs(dw0.eclipse.intp_height-150))]*1e6*k)/(2*np.pi)
)
f_dif0 = (f0_240-f0_150)
f_dif0_max = (abs(np.min(f_dif0))*(2*np.pi))**2/k/1e6
f_dif0_dur = 5*np.sum(f_dif0 <= 0, axis=0)
del dw0

dw1 = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=None,
        loc={"lat": lat_low, "lon":lonx},
    )
k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
e_eclipse1 = dw1.eclipse.dataset[0]["interpol_value"]["value"]
ttime1 = dw1.eclipse.dataset[0]["interpol_value"]["time"]
f1_240, f1_150 = (
    np.sqrt(e_eclipse1[:, np.argmin(np.abs(dw1.eclipse.intp_height-240))]*1e6*k)/(2*np.pi), 
    np.sqrt(e_eclipse1[:, np.argmin(np.abs(dw1.eclipse.intp_height-150))]*1e6*k)/(2*np.pi)
)
f_dif1 = (f1_240-f1_150)
f_dif1_max = (abs(np.min(f_dif1))*(2*np.pi))**2/k/1e6
f_dif1_dur = 5*np.sum(f_dif1 <= 0, axis=0)
del dw1

dw2 = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=None,
        loc={"lat": lat_up, "lon":lonx},
    )
k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
e_eclipse2 = dw2.eclipse.dataset[0]["interpol_value"]["value"]
ttime2 = dw2.eclipse.dataset[0]["interpol_value"]["time"]
f2_240, f2_150 = (
    np.sqrt(e_eclipse2[:, np.argmin(np.abs(dw2.eclipse.intp_height-240))]*1e6*k)/(2*np.pi), 
    np.sqrt(e_eclipse2[:, np.argmin(np.abs(dw2.eclipse.intp_height-150))]*1e6*k)/(2*np.pi)
)
f_dif2 = (f2_240-f2_150)
f_dif2_max = (abs(np.min(f_dif2))*(2*np.pi))**2/k/1e6
f_dif2_dur = 5*np.sum(f_dif2 <= 0, axis=0)
del dw2

mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
#%matplotlib inline
fig = plt.figure(figsize=(10, 5), dpi=300)
ax = fig.add_subplot(
            221,
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
gl.ylabel_style = {'size':10, 'color': 'b'}
gl.xlabel_style = {"size":10, 'color': 'b',}

XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), glon, glat)
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]

im = ax.contourf(
    X, Y, of.T, transform=ccrs.PlateCarree(), cmap="gray", alpha=0.5,
    levels=[0,.2,.4,.6,.8,1.],
    vmin=0, vmax=1, **kwargs
)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
cb = fig.colorbar(im, ax=ax, cax=cax)
cb.set_label(r"Occultation, $\mathcal{O}$")
ax.text(0.05, 0.95, "(a)", transform=ax.transAxes, ha="left", va="center")

XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), 
                                          np.array([[lonx, lonx, lonx]]), 
                                          np.array([[center, lat_low, lat_up]]))
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
ax.scatter(X, Y, s=10, color="magenta", marker=".", fc="None", lw=0.3)

ax = fig.add_subplot(222)
ax.text(0.05, 0.95, "(b)", transform=ax.transAxes, ha="left", va="center")
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
ax.set_xlabel("UT", fontdict={"size":12})
ax.set_ylabel(
    r"$f_0$, MHz", fontdict={"size":12}
)
a = of[np.argmin(np.abs(lat-center)), np.argmin(np.abs(lon-lonx))]
olat, olon, _ = aacgmv2.get_aacgm_coord(center, lonx, 300, dt.datetime(2017,8,21,18))
ax.text(0.95, 0.98, r"$\theta=%d^{\circ},\phi=%d^{\circ}$"%(olat,olon)+\
        "\n"+r"$\mathcal{O}=%.2f,\Delta T_{GC}=%d min, GC^p=%.1f$"%(a,f_dif0_dur,f_dif0_max)+" el/cc", 
        transform=ax.transAxes, ha="right", va="top", fontdict={"size":10})
ax.plot(ttime0, f0_150/1e6, "ro", ms=0.4, ls="None")
ax.plot(ttime0, f0_240/1e6, "bs", ms=0.4, ls="None")
ax.set_ylim(2,7)
ax.set_xlim(dt.datetime(2017,8,21,16),dt.datetime(2017,8,21,22))

ax = fig.add_subplot(223)
fig.subplots_adjust(wspace=0.5)
ax.text(0.05, 0.95, "(c)", transform=ax.transAxes, ha="left", va="center")
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
ax.set_xlabel("UT", fontdict={"size":12})
ax.set_ylabel(
    r"$f_0$, MHz", fontdict={"size":12}
)
ax.plot(ttime1, f1_150/1e6, "ro", ms=0.4, ls="None")
ax.plot(ttime1, f1_240/1e6, "bs", ms=0.4, ls="None")
ax.set_ylim(2,7)
a = of[np.argmin(np.abs(lat-lat_low)), np.argmin(np.abs(lon-lonx))]
ax.set_xlim(dt.datetime(2017,8,21,16),dt.datetime(2017,8,21,22))
olat, olon, _ = aacgmv2.get_aacgm_coord(lat_low, lonx, 300, dt.datetime(2017,8,21,18))
ax.text(0.95, 0.98, 
        r"$\theta=%d^{\circ},\phi=%d^{\circ}$"%(olat,olon)+\
        "\n"+r"$\mathcal{O}=%.2f,\Delta T_{GC}=%d min, GC^p=%.1f$"%(a,f_dif1_dur,f_dif1_max)+" el/cc", 
        transform=ax.transAxes, ha="right", va="top", fontdict={"size":10})

ax = fig.add_subplot(224)
fig.subplots_adjust(wspace=0.5)
ax.text(0.05, 0.95, "(d)", transform=ax.transAxes, ha="left", va="center", fontdict={"size":10})
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
ax.set_xlabel("UT", fontdict={"size":12})
ax.set_ylabel(
    r"$f_0$, MHz", fontdict={"size":12}
)
a = of[np.argmin(np.abs(lat-lat_up)), np.argmin(np.abs(lon-lonx))]
ax.plot(ttime2, f2_150/1e6, "ro", ms=0.4, ls="None")
ax.plot(ttime2, f2_240/1e6, "bs", ms=0.4, ls="None")
ax.set_ylim(2,7)
ax.set_xlim(dt.datetime(2017,8,21,16),dt.datetime(2017,8,21,22))
olat, olon, _ = aacgmv2.get_aacgm_coord(lat_up, lonx, 300, dt.datetime(2017,8,21,18))
ax.text(0.95, 0.98, r"$\theta=%d^{\circ},\phi=%d^{\circ}$"%(olat,olon)+\
        "\n"+r"$\mathcal{O}=%.2f,\Delta T_{GC}=%d min, GC^p=%.1f$"%(a,f_dif2_dur,f_dif2_max)+" el/cc", 
        transform=ax.transAxes, ha="right", va="top", fontdict={"size":10})

fig.subplots_adjust(wspace=0.5, hspace=0.3)
fig.savefig("dataset/figures/latitude_distribution.png", bbox_inches="tight")