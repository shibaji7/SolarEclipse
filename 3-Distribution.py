# G-Condition Study
import matplotlib.pyplot as plt
import datetime as dt
import sys
sys.path.append("py/")
from fetch import *
import utils
import mplstyle

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

def get_gridded_parameters(q, xparam="lat", yparam="lon", zparam="gc_del_t"):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[ [xparam, yparam, zparam] ]
    plotParamDF[xparam] = plotParamDF[xparam].tolist()
    plotParamDF[yparam] = plotParamDF[yparam].tolist()
    plotParamDF = plotParamDF.groupby( [xparam, yparam] ).mean().reset_index()
    plotParamDF = plotParamDF[ [xparam, yparam, zparam] ].pivot( xparam, yparam )
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y  = np.meshgrid( x, y )
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
            np.isnan(plotParamDF[zparam].values),
            plotParamDF[zparam].values)
    return X,Y,Z

def fetch_occl_data(d):
    folder = "dataset/Mark/euv/%s_%dkm_171_1.nc"
    file = folder % (d.strftime("%Y%m%d%H%M%S"), 150)
    
    ds = nc.Dataset(file)
    of = ds.variables["of"][:]
    lat, lon = ds.variables["glat"][:], ds.variables["glon"][:]
    glat, glon = np.meshgrid(lat, lon)
    return lat, lon, glat, glon, of

date = dt.datetime(2017,8,21,18)
lat, lon, glat, glon, of = fetch_occl_data(date)

params = [
    {
        "name": "e",
        "unit_multiplier": 1.,
        "unit": "cc",
        "mol_to_cc": "T",
        "interpolate": {"scale": "log", "type": "linear"},
    },
]

lon_list, lat_list = (
    np.arange(-150, -60, 1), 
    np.arange(20, 80, 1)
)

f_dif_max = np.nan*np.zeros((len(lat_list), len(lon_list)))
f_dif_dur = np.nan*np.zeros((len(lat_list), len(lon_list)))

for i, lat in enumerate(lat_list):
    for j, lon in enumerate(lon_list):
        print(len(lat_list)*len(lon_list),"----", (len(lon_list)*i+j))
        file = f"dataset/ds/{int(lat)}_{int(lon)}.txt"
        if not os.path.exists(file):
            dw0 = DiffWACCMX(
                    "dataset/29jan2022/ECL_1100_0100_all.nc",
                    "dataset/29jan2022/BGC_1100_0100_all.nc",
                    params=params,
                    stn=None,
                    loc={"lat": lat, "lon":lon},
                )
            k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
            e_eclipse0 = dw0.eclipse.dataset[0]["interpol_value"]["value"]
            ttime0 = dw0.eclipse.dataset[0]["interpol_value"]["time"]
            f0_240, f0_150 = (
                np.sqrt(e_eclipse0[:, np.argmin(np.abs(dw0.eclipse.intp_height-240))]*1e6*k)/(2*np.pi), 
                np.sqrt(e_eclipse0[:, np.argmin(np.abs(dw0.eclipse.intp_height-150))]*1e6*k)/(2*np.pi)
            )
            f_dif = (f0_240-f0_150)
            f_dif_max[i,j] = (abs(np.min(f_dif))*(2*np.pi))**2/k/1e6
            f_dif_dur[i,j] = 5*np.sum(f_dif <= 0, axis=0)
            o = np.array([
                (abs(np.min(f_dif))*(2*np.pi))**2/k/1e6, 5*np.sum(f_dif <= 0, axis=0)
            ])
            np.savetxt(file, o, fmt='%.18e')
            del dw0
import glob
import pandas as pd
files = glob.glob("dataset/ds/*.txt")
obj = []
for f in files:
    lat, lon = (
        int(f.split("/")[-1].replace(".txt", "").split("_")[0]),
        int(f.split("/")[-1].replace(".txt", "").split("_")[1])
    )
    dat = np.loadtxt(f)
    obj.append(dict(
        lat=lat, lon=lon, 
        gc_del_t=dat[1] if dat[1] > 0 else np.nan, 
        gc_del_n=dat[0] if dat[1] > 0 else np.nan
    ))
obj = pd.DataFrame.from_records(obj)
print(obj.head())

mpl.rcParams.update({"xtick.labelsize":8, "ytick.labelsize":8, "font.size":8})
#%matplotlib inline
fig = plt.figure(figsize=(3, 4.5), dpi=240)
ax = fig.add_subplot(
            311,
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

ax = fig.add_subplot(
            312,
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
Lat, Lon, Val = get_gridded_parameters(obj)
Val[Val==np.nan] = 0
from scipy.ndimage import gaussian_filter
Val = gaussian_filter(Val, sigma=0.5)
XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), Lon, Lat)
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
im = ax.contourf(
    X, Y, Val.T, transform=ccrs.PlateCarree(), cmap="gray_r", alpha=0.5,
    #levels=[0,.2,.4,.6,.8,1.],
    vmin=0, vmax=150, **kwargs
)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
cb = fig.colorbar(im, ax=ax, cax=cax)
cb.set_label(r"$\Delta T_{GC}$, minutes")
ax.text(0.05, 0.95, "(b)", transform=ax.transAxes, ha="left", va="center")

ax = fig.add_subplot(
            313,
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
Lat, Lon, Val = get_gridded_parameters(obj, zparam="gc_del_n")
Val[Val==np.nan] = 0
from scipy.ndimage import gaussian_filter
Val = gaussian_filter(Val, sigma=0.5)
XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), Lon, Lat)
X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
im = ax.contourf(
    X, Y, Val.T, transform=ccrs.PlateCarree(), cmap="gray_r", alpha=0.5,
    #levels=[0,.2,.4,.6,.8,1.],
    vmin=0, vmax=5000, **kwargs
)
cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
cb = fig.colorbar(im, ax=ax, cax=cax)
cb.set_label(r"$GC^p$, /cc")
ax.text(0.05, 0.95, "(c)", transform=ax.transAxes, ha="left", va="center")

fig.subplots_adjust(wspace=0.5, hspace=0.3)
fig.savefig("dataset/figures/globe_distribution.png", bbox_inches="tight")