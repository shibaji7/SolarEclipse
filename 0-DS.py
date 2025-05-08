import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


import datetime as dt
import sys
sys.path.append("py/")
from fetch import *
import mplstyle

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import matplotlib as mpl
mpl.rcParams.update({"xtick.labelsize": 15, "ytick.labelsize":15, "font.size":15})
plt.rcParams.update({
    "text.usetex": False,
})

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return

from netCDF4 import Dataset
co_ds_240 = Dataset("database/EOFs/39.992lat_-105.455lon_20170821160000_20170821193000-240.nc")
co_ds_150 = Dataset("database/EOFs/39.992lat_-105.455lon_20170821160000_20170821193000-150.nc")

ls_ds_240 = Dataset("database/EOFs/42.75lat_-104.455lon_20170821160000_20170821193000-240.nc")
ls_ds_150 = Dataset("database/EOFs/42.75lat_-104.455lon_20170821160000_20170821193000-150.nc")

digi = FetchDigisonde("dataset/Mabie/lusk.sav")
of150x = FetchOccult([dt.datetime(2017,8,21,16), 42], 150)
of240x = FetchOccult([dt.datetime(2017,8,21,16), 42], 240)
print(digi.data)

digi_co = FetchDigisonde("dataset/Mabie/boulder.sav")
digi_co.data.head()

headersm = [
    'GDALT', 
    'GDLAT',
    'GLON', 
    'TEC', 
    'NEMAX',
    'HMAX', 
    'POPL', 
    'DPOPL',
    'TI',
    'DTI', 
    'PM',
    'DPM',
    'VI3', 
    'DVI3', 
    'TR',
    'DTR',
    'CCTITR',
    'DCCTITR',
    'NEL',
    'DNEL',
    'NE',
    'DNE'
]
headersj = [
    "NE",
    "TI",
    "TR",
    "VO",
    "PH+",
    "PM",
    "CO",
    "GDALT",
    "RANGE",
    "AZ1",
    "AZ2",
    "EL1",
    "EL2",
    "GDLAT",
    "GLON",
]

of150 = FetchOccult([dt.datetime(2017,8,21,15,15), 72], 150, "mil")
of240 = FetchOccult([dt.datetime(2017,8,21,15,15), 72], 240, "mil")
isr = FetchISR("dataset/ISR/mlh170821m.004.txt", headers=headersm)
ISR = isr.data.copy()
def get_gridded_parameters(q, xparam="DATE", yparam="GDALT", zparam="NE", r=0, rounding=False):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[ [xparam, yparam, zparam] ]
    if rounding:
        plotParamDF.loc[:, xparam] = np.round(plotParamDF[xparam].tolist(), r)
        plotParamDF.loc[:, yparam] = np.round(plotParamDF[yparam].tolist(), r)
    else:
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
X, Y, Z = get_gridded_parameters(ISR)
hmax, Tx = [], []
ISRO = ISR[(ISR.DATE>=dt.datetime(2017,8,21,15)) & (ISR.DATE<=dt.datetime(2017,8,21,21))]
for t in ISRO.DATE.unique():
    Tx.append(t)
    u = ISRO[ISRO.DATE==t]
    hmax.append(u.GDALT.tolist()[np.argmax(u.NE)])

params = [
    {
        "name": "e",
        "unit_multiplier": 1.,
        "unit": "cc",
        "mol_to_cc": "T",
        "interpolate": {"scale": "log", "type": "linear"},
    }
]
dw = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn="lusk",
    )

odf = []
with open("dataset/ISR/MLH_ISR_Ne_gridded_21-Aug-2017.txt", "r") as f:
    lines = f.readlines()
    for ix, line in enumerate(lines[1:]):
        line = list(filter(None, line.replace("\n", "").split(" ")))
        odf.append({"UT": float(line[0]), "alt": float(line[1]), "el": float(line[2])})
odf = pd.DataFrame.from_records(odf)
print(odf.head(1000))

k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
e_eclipse = dw.eclipse.dataset[0]["interpol_value"]["value"]
ttime = dw.eclipse.dataset[0]["interpol_value"]["time"]
f_240, f_150 = (
    np.sqrt(e_eclipse[:, np.argmin(np.abs(dw.eclipse.intp_height-240))]*1e6*k)/(2*np.pi), 
    np.sqrt(e_eclipse[:, np.argmin(np.abs(dw.eclipse.intp_height-150))]*1e6*k)/(2*np.pi)
)
print(e_eclipse.shape)

f_dif = (f_240-f_150)
f_dif_max = (abs(np.min(f_dif))*(2*np.pi))**2/k/1e6
f_dif_dur = 5*np.sum(f_dif <= 0, axis=0)

import copy
fig = plt.figure(dpi=300, figsize=(8,9))
ax = fig.add_subplot(311)
ax.set_ylabel(r"$foF_{1,2}$, MHz")
ax.set_xlabel("UT")
ax.set_ylim(2,6)
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,1)))
ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0,60,5)))
#ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
ax.set_xlim([dt.datetime(2017,8,21,16), dt.datetime(2017,8,21,19, 30)])
ax.plot(digi.data.time, digi.data.fof1, "ro", ms=1.5, ls="None")
ax.plot(digi.data.time, digi.data.fof2, "bo", ms=1.5, ls="None")
ax0 = copy.copy(ax)
ax = ax.twinx()
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,1)))
ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0,60,5)))
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
ax.set_ylim(0, 1.2)
ax.set_ylabel("Obscuration")
# ax.plot(of150x.data.time, of150x.data.of, color="darkred", ls="None", marker="+", ms=4, alpha=0.6)
# ax.plot(of240x.data.time, of240x.data.of, color="darkblue", ls="None", marker="+", ms=4, alpha=0.6)
ax.plot([dt.datetime(2017,8,21,16) + dt.timedelta(minutes=int(i))
         for i in ls_ds_150.variables["time"][:].tolist()[::5]], co_ds_150.variables["171"][:].tolist()[::5], 
        color="darkred", marker="+", ms=4, ls="None", alpha=0.8)
ax.plot([dt.datetime(2017,8,21,16) + dt.timedelta(minutes=int(i))
         for i in ls_ds_240.variables["time"][:].tolist()[::5]], co_ds_240.variables["171"][:].tolist()[::5], 
        color="darkblue", marker="+", ms=4, ls="None", alpha=0.6)
ax.axvline(dt.datetime(2017,8,21,16,15), ls="--", lw=0.5, color="k")
ax.axvline(dt.datetime(2017,8,21,17,45), ls="--", lw=0.5, color="k")
ax.axvline(dt.datetime(2017,8,21,19,15), ls="--", lw=0.5, color="k")
ax.axvspan(dt.datetime(2017,8,21,18,15), dt.datetime(2017,8,21,19,5), alpha=0.1)
ax.text(0.05, 0.9, "(a) Lusk, WY", ha="left", va="center", transform=ax.transAxes)

bounds = list(np.arange(0.5, 3.01, .5))
cmap = plt.cm.jet

from scipy import constants as C
k = C.elementary_charge**2/(C.m_e*C.epsilon_0)

Tx, hmax = [], []
for t in odf.UT.unique():
    Tx.append(t)
    u = odf[odf.UT==t]
    hmax.append(u.alt.tolist()[np.argmax(u.el)])

ax = fig.add_subplot(313)
ax.set_ylabel("Height, km")
ax.set_xlabel("UT")
X, Y, Z = get_gridded_parameters(odf, "UT", "alt", "el")
im = ax.pcolor(X, Y, Z.T/1e11, cmap=cmap, vmax=3, vmin=0)
ax.set_ylim(120, 350)
ax.set_xlim(15, 22)
ax.plot(Tx, hmax, color="k", marker="+", ms=4, markeredgewidth=0.3, ls="None")
ax.text(0.05, 0.9, "(c) MHISR, MA", ha="left", va="center", transform=ax.transAxes, fontdict={"color":"k"})
ax.axvline(17.50, ls="--", lw=0.6, color="k")
ax.axvline(18.75, ls="--", lw=0.6, color="k")
ax.axvline(20.00, ls="--", lw=0.6, color="k")


pos = ax.get_position()
cpos = [pos.x1 + 0.025, pos.y0 + 0.0125,
        0.015, pos.height * 0.7]                # this list defines (left, bottom, width, height
cax = fig.add_axes(cpos)
cb2 = fig.colorbar(im, cax,
                   spacing="uniform",
                   orientation="vertical", 
                   cmap=cmap)
cb2.set_label(r"$N_e$ ($\times 10^{11}$), $m^{-3}$")


ax = fig.add_subplot(312)
ax.set_ylabel(r"$foF_{1,2}$, MHz")
ax.set_xlabel("UT")
ax.set_ylim(2,6)
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,1)))
ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0,60,5)))
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
#ax.plot(ttime, f_150/1e6, "ro", ms=0.4, ls="None")
#ax.plot(ttime, f_240/1e6, "bs", ms=0.4, ls="None")
#co_ds.variables["alt_km"]
ax.plot(digi_co.data.time, digi_co.data.fof1, "ro", ms=1.5, ls="None")
ax.plot(digi_co.data.time, digi_co.data.fof2, "bo", ms=1.5, ls="None")
ax.axvline(dt.datetime(2017,8,21,16,15), ls="--", lw=0.5, color="k")
ax.axvline(dt.datetime(2017,8,21,17,45), ls="--", lw=0.5, color="k")
ax.axvline(dt.datetime(2017,8,21,19,15), ls="--", lw=0.5, color="k")
ax.set_xlim([dt.datetime(2017,8,21,16), dt.datetime(2017,8,21,19, 30)])
ax.text(0.05, 0.9, "(b) Boulder, CO", ha="left", va="center", transform=ax.transAxes)
_ = ax.axvspan(dt.datetime(2017,8,21,18,2), dt.datetime(2017,8,21,18,55), alpha=0.1)
ax = ax.twinx()
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,1)))
ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0,60,5)))
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
ax.set_ylim(0, 1.2)
ax.set_ylabel("Obscuration")
ax.plot([dt.datetime(2017,8,21,16) + dt.timedelta(minutes=int(i))
         for i in co_ds_150.variables["time"][:].tolist()[::5]], co_ds_150.variables["171"][:].tolist()[::5], 
        color="darkred", marker="+", ms=4, ls="None", alpha=0.8)
ax.plot([dt.datetime(2017,8,21,16) + dt.timedelta(minutes=int(i))
         for i in co_ds_240.variables["time"][:].tolist()[::5]], co_ds_240.variables["171"][:].tolist()[::5], 
        color="darkblue", marker="+", ms=4, ls="None", alpha=0.6)

fig.subplots_adjust(wspace=0.3,hspace=0.3)
fig.savefig("dataset/figures/intro_to_gc.png", bbox_inches="tight")
plt.close()


fig = plt.figure(dpi=300, figsize=(8,3))
ax = fig.add_subplot(111)
ax.set_ylabel(r"$f_{150,240}$, MHz")
ax.set_xlabel("UT")
ax.set_ylim(2,6)
ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,1)))
ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0,60,5)))
#ax.xaxis.set_minor_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
ax.set_xlim([dt.datetime(2017,8,21,16), dt.datetime(2017,8,21,19, 30)])
ttime = dw.eclipse.dataset[0]["interpol_value"]["time"]
ax.plot(ttime, f_150*1e-6, "ro", ms=0.8, ls="None")
ax.plot(ttime, f_240*1e-6, "bo", ms=0.8, ls="None")
ax.axvline(dt.datetime(2017,8,21,16,15), ls="--", lw=0.5, color="k")
ax.axvline(dt.datetime(2017,8,21,17,45), ls="--", lw=0.5, color="k")
ax.axvline(dt.datetime(2017,8,21,19,15), ls="--", lw=0.5, color="k")
ax.axvspan(dt.datetime(2017,8,21,18,15), dt.datetime(2017,8,21,19,5), alpha=0.1)
ax.text(0.05, 0.9, "Lusk, WY", ha="left", va="center", transform=ax.transAxes)
ax.text(0.95, 0.9, r"$\Delta T_{GC}$=%d"%f_dif_dur + "\n" + r"$GC^p$=%.1f/cc"%f_dif_max, ha="right", va="center", transform=ax.transAxes)
fig.subplots_adjust(wspace=0.3,hspace=0.3)
fig.savefig("dataset/figures/gc_matrices.png", bbox_inches="tight")
plt.close()

fig = plt.figure(dpi=240, figsize=(5,7))
ax = fig.add_subplot(311)
ax.set_ylabel(r"$foF_{1,2}^{v,d}$, MHz")
# ax.set_xlabel("UT")
ax.set_ylim(2,6)
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
ax.set_xlim([dt.datetime(2017,8,21,16), dt.datetime(2017,8,21,19, 30)])
ax.plot(digi.data.time, digi.data.fof1, "ro", ms=0.4, ls="None")
ax.plot(digi.data.time, digi.data.fof2, "bs", ms=0.4, ls="None")
ax = ax.twinx()
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
ax.set_ylim(0, 1.2)
ax.set_ylabel("Obscuration")
ax.plot(of150x.data.time, of150x.data.of, "ko", ms=0.4, ls="None")
ax.plot(of240x.data.time, of240x.data.of, "ks", ms=0.4, ls="None")
ax.axvline(dt.datetime(2017,8,21,16,15), ls="--", lw=0.2, color="k")
ax.axvline(dt.datetime(2017,8,21,17,50), ls="--", lw=0.2, color="k")
ax.axvline(dt.datetime(2017,8,21,19,15), ls="--", lw=0.2, color="k")
ax.axvspan(dt.datetime(2017,8,21,18,15), dt.datetime(2017,8,21,19,5), alpha=0.1)
ax.text(0.05, 0.9, "(a)", ha="left", va="center", transform=ax.transAxes)

ax = fig.add_subplot(312)
ax.set_ylabel(r"$foF_{1,2}^{v,m}$, MHz")
ax.set_xlabel("UT")
ax.set_ylim(2,6)
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"$%H^{%M}$"))
ax.plot(ttime, f_150/1e6, "ro", ms=0.4, ls="None")
ax.plot(ttime, f_240/1e6, "bs", ms=0.4, ls="None")
ax.axvline(dt.datetime(2017,8,21,16,15), ls="--", lw=0.2, color="k")
ax.axvline(dt.datetime(2017,8,21,17,50), ls="--", lw=0.2, color="k")
ax.axvline(dt.datetime(2017,8,21,19,15), ls="--", lw=0.2, color="k")
ax.set_xlim([dt.datetime(2017,8,21,16), dt.datetime(2017,8,21,19, 30)])
ax.text(0.05, 0.9, "(b)", ha="left", va="center", transform=ax.transAxes)
ax.text(0.95, 0.98, r"$\Delta T_{GC}=%d min, GC^p=%.1f /cc$"%(f_dif_dur,f_dif_max), 
        transform=ax.transAxes, ha="right", va="top", fontdict={"size":8})
_ = ax.axvspan(dt.datetime(2017,8,21,18,15), dt.datetime(2017,8,21,19,5), alpha=0.1)
fig.subplots_adjust(wspace=0.3,hspace=0.3)
fig.savefig("dataset/figures/data_model_comparison.png", bbox_inches="tight")
plt.close()
