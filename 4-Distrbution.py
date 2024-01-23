import matplotlib.pyplot as plt
import datetime as dt
import sys
sys.path.extend(["py/","pyrt/"])
from fetch2D import *
from fetch import FetchModel
import util
import mplstyle
import copy


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
date = dt.datetime(2017,8,21,18)
k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
lats = np.arange(20, 70, 2)
lons = np.arange(-140, -40, 2)
Lats, Lons = np.meshgrid(lats, lons)

# eclipse = Fetch2DModel(
#     "dataset/29jan2022/ECL_1100_0100_all.nc", 
#     copy.copy(params), date, 2, 2, 150
# )
eclipse = FetchModel(
    "dataset/29jan2022/ECL_1100_0100_all.nc",
    copy.copy(params), date
)

for i in range(len(lats)):
    for j in range(len(lons)):
        eclipse.run_height_interpolate(
            loc={"lat": lats[i], "lon":lons[j]},
        )
        print("--")
        o = eclipse.fetch_data()
        dif = (o["e_150"].e - o["e_240"].e).iloc[0]
        if dif>0:
            print(lats[i], lons[j], dif)
        
#         eden_150 = eclipse.run_2D_interpolation(dict(
#             lat=lats[i], lon=lons[j], h=150
#         ))
#         eden_240 = eclipse.run_2D_interpolation(dict(
#             lat=lats[i], lon=lons[j], h=240
#         ))
#         f0_240, f0_150 = (
#                 np.sqrt(eden_240["e"]*1e6*k)/(2*np.pi), 
#                 np.sqrt(eden_150["e"]*1e6*k)/(2*np.pi)
#             )
#         f_dif = (f0_240-f0_150)
#         print(lats[i], lons[j], f_dif)
        #break

# print(eclipse.time[0], eclipse.time[-1])
# timeI = eclipse.time.index(date)
# print(timeI, eclipse.time[timeI])
# eden = eclipse.dataset[0]["value"][timeI, :, :, :]
# latx, lonx = eclipse.latx, eclipse.lonx
# print(np.min(lonx), np.max(lonx))
# latsx = latx[(latx>20) & (latx<70)]
# lonsx = lonx[(lonx>-140) & (lonx<-40)]
# print(np.min(lonsx), np.max(lonsx))
# print(eden.shape, len(latsx), len(lonsx), len(eclipse.latx), len(eclipse.lonx))
# eden = eden[:,(latx>20) & (latx<70),:]
# eden = eden[:,:,(lonx>-140) & (lonx<-40)]
# print(latsx, lonsx, eden.shape)

# f_dif_max = np.nan*np.zeros((len(latsx),len(lonsx)))
# f_dif_dur = np.nan*np.zeros((len(latsx),len(lonsx)))
# n_dif_max = np.nan*np.zeros((len(latsx),len(lonsx)))

# print(eclipse.Z3.shape)
# Z3 = eclipse.Z3[timeI,:,:,:]
# Z3 = Z3[:,(latx>20) & (latx<70),:]
# print(Z3.shape)
# Z3 = Z3[:,:,(lonx>-140) & (lonx<-40)]
# print(len(latsx), len(lonsx), Z3.shape, eden.shape)
# hx = np.arange(100,300)
# for i in range(eden.shape[1]):
#     for j in range(eden.shape[2]):
#         eden_new = util.interpolate_by_altitude(
#             Z3[:,i,j], hx, 
#             eden[:,i,j]
#         )
#         alt_H150_index = hx.tolist().index(150)#np.abs(Z3[:,i,j]-160).argmin()
#         alt_H240_index = hx.tolist().index(240)#np.abs(Z3[:,i,j]-220).argmin()
#         #print(Z3[alt_H150_index,i,j],Z3[alt_H240_index,i,j])
#         f0_240, f0_150 = (
#                 np.sqrt(eden_new[alt_H240_index]*1e6*k)/(2*np.pi), 
#                 np.sqrt(eden_new[alt_H150_index]*1e6*k)/(2*np.pi)
#             )     
#         n0_240, n0_150 = (
#             eden_new[alt_H240_index],
#             eden_new[alt_H150_index]
#         )
#         f_dif = (f0_240-f0_150)
#         n_dif = n0_240-n0_150
#         print(
#                 lonsx[j], latsx[i],
#                 hx[alt_H150_index], hx[alt_H240_index], 
#                 f_dif, n_dif
#             )
#         # if n_dif < 0:
#         #     print(
#         #         lonsx[j], latsx[i],
#         #         Z3[alt_H150_index,i,j], Z3[alt_H240_index,i,j], 
#         #         f_dif, n_dif, eclipse_den[alt_H150_index,i,j], 
#         #         eclipse_den[alt_H240_index,i,j]
#         #     )
#         #     pass
#         f_dif_max[i,j] = (np.min(f_dif))*(2*np.pi)**2/k/1e6
#         f_dif_dur[i,j] = 5*np.sum(f_dif <= 0, axis=0)
#         n_dif_max[i,j] = n_dif
#print(np.nanmax(n_dif_max), np.nanmin(n_dif_max))

# def fetch_occl_data(d):
#     folder = "dataset/Mark/euv/%s_%dkm_171_1.nc"
#     file = folder % (d.strftime("%Y%m%d%H%M%S"), 150)
    
#     ds = nc.Dataset(file)
#     of = ds.variables["of"][:]
#     lat, lon = ds.variables["glat"][:], ds.variables["glon"][:]
#     glat, glon = np.meshgrid(lat, lon)
#     return lat, lon, glat, glon, of
# lat, lon, glat, glon, of = fetch_occl_data(date)

# #print(np.nanmin(of), np.nanmax(of))

# select_latlon = of>0

# mpl.rcParams.update({"xtick.labelsize":8, "ytick.labelsize":8, "font.size":8})
# #%matplotlib inline
# fig = plt.figure(figsize=(3, 4.5), dpi=240)
# ax = fig.add_subplot(
#             311,
#             projection=ccrs.PlateCarree()
#         )
# kwargs = {}
# kwargs["edgecolor"] = "k"
# kwargs["facecolor"] = "none"
# kwargs["linewidth"] = 0.3
# feature = cartopy.feature.NaturalEarthFeature(
#             "physical", "coastline", "50m", **kwargs
#         )
# ax.add_feature(feature, **kwargs)
# ax.set_extent([-140, -40, 20, 70], crs=ccrs.PlateCarree())
# plt_lons = np.arange(-180, 181, 30)
# mark_lons = np.arange(-180, 180, 30)
# plt_lats = np.arange(20, 90, 15)
# gl = ax.gridlines(ccrs.PlateCarree(), linewidth=0.5, draw_labels=True,)
# gl.xlabels_top = False
# gl.ylabels_right = False
# gl.xlocator = mticker.FixedLocator(plt_lons)
# gl.ylocator = mticker.FixedLocator(plt_lats)
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# gl.n_steps = 90
# gl.ylabel_style = {'size':8, 'color': 'b'}
# gl.xlabel_style = {"size":8, 'color': 'b'}

# XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), glon, glat)
# X, Y = XYZ[:, :, 0], XYZ[:, :, 1]

# im = ax.contourf(
#     X, Y, of.T, transform=ccrs.PlateCarree(), cmap="gray", alpha=0.5,
#     levels=[0,.2,.4,.6,.8,1.],
#     vmin=0, vmax=1, **kwargs
# )
# cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
# cb = fig.colorbar(im, ax=ax, cax=cax)
# cb.set_label(r"Occultation, $\mathcal{O}$")
# ax.text(0.05, 0.95, "(a)", transform=ax.transAxes, ha="left", va="center")

# ax = fig.add_subplot(
#             312,
#             projection=ccrs.PlateCarree()
#         )
# kwargs = {}
# kwargs["edgecolor"] = "k"
# kwargs["facecolor"] = "none"
# kwargs["linewidth"] = 0.3
# feature = cartopy.feature.NaturalEarthFeature(
#             "physical", "coastline", "50m", **kwargs
#         )
# ax.add_feature(feature, **kwargs)
# ax.set_extent([-140, -40, 20, 70], crs=ccrs.PlateCarree())
# plt_lons = np.arange(-180, 181, 30)
# mark_lons = np.arange(-180, 180, 30)
# plt_lats = np.arange(20, 90, 15)
# gl = ax.gridlines(ccrs.PlateCarree(), linewidth=0.5, draw_labels=True,)
# gl.xlabels_top = False
# gl.ylabels_right = False
# gl.xlocator = mticker.FixedLocator(plt_lons)
# gl.ylocator = mticker.FixedLocator(plt_lats)
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# gl.n_steps = 90
# gl.ylabel_style = {'size':8, 'color': 'b'}
# gl.xlabel_style = {"size":8, 'color': 'b'}
# Lon, Lat = np.meshgrid(eclipse.lonx, eclipse.latx)
# XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), Lon, Lat)
# f_dif_dur[f_dif_dur<=0] = np.nan
# X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
# im = ax.contourf(
#     X, Y, f_dif_dur, transform=ccrs.PlateCarree(), cmap="gray_r", alpha=0.5,
#     #levels=[0,.2,.4,.6,.8,1.],
#     vmin=0, vmax=150, **kwargs
# )
# cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
# cb = fig.colorbar(im, ax=ax, cax=cax)
# cb.set_label(r"$\Delta T_{GC}$, minutes")
# ax.text(0.05, 0.95, "(b)", transform=ax.transAxes, ha="left", va="center")

# ax = fig.add_subplot(
#             313,
#             projection=ccrs.PlateCarree()
#         )
# kwargs = {}
# kwargs["edgecolor"] = "k"
# kwargs["facecolor"] = "none"
# kwargs["linewidth"] = 0.3
# feature = cartopy.feature.NaturalEarthFeature(
#             "physical", "coastline", "50m", **kwargs
#         )
# ax.add_feature(feature, **kwargs)
# ax.set_extent([-140, -40, 20, 70], crs=ccrs.PlateCarree())
# plt_lons = np.arange(-180, 181, 30)
# mark_lons = np.arange(-180, 180, 30)
# plt_lats = np.arange(20, 90, 15)
# gl = ax.gridlines(ccrs.PlateCarree(), linewidth=0.5, draw_labels=True,)
# gl.xlabels_top = False
# gl.ylabels_right = False
# gl.xlocator = mticker.FixedLocator(plt_lons)
# gl.ylocator = mticker.FixedLocator(plt_lats)
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
# gl.n_steps = 90
# gl.ylabel_style = {'size':8, 'color': 'b'}
# gl.xlabel_style = {"size":8, 'color': 'b'}
# n_dif_max[n_dif_max>=0] = np.nan
# XYZ = ccrs.PlateCarree().transform_points(ccrs.Geodetic(), Lon, Lat)
# X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
# im = ax.contourf(
#     X, Y, n_dif_max, transform=ccrs.PlateCarree(), cmap="gray_r", alpha=0.5,
#     #levels=[0,.2,.4,.6,.8,1.],
#     #vmin=0, vmax=5000, **kwargs
# )
# cax = ax.inset_axes([1.05, 0.1, 0.05, 0.8], transform=ax.transAxes)
# cb = fig.colorbar(im, ax=ax, cax=cax)
# cb.set_label(r"$GC^p$, /cc")
# ax.text(0.05, 0.95, "(c)", transform=ax.transAxes, ha="left", va="center")

# fig.subplots_adjust(wspace=0.5, hspace=0.3)
# fig.savefig("dataset/figures/globe_distribution.png", bbox_inches="tight")