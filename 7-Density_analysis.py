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


params = [
    {
        "name": "TElec",
        "unit_multiplier": 1.,
        "unit": "K",
        "interpolate": {"scale": "linear", "type": "linear"},
    },
    {
        "name": "TIon",
        "unit_multiplier": 1.,
        "unit": "K",
        "interpolate": {"scale": "linear", "type": "linear"},
    },
    {
        "name": "UI",
        "unit_multiplier": 1.,
        "unit": "m/s",
        "interpolate": {"scale": "linear", "type": "linear"},
    }
]
dw = DiffWACCMX(
    "dataset/29jan2022/ECL_1100_0100_all.nc",
    "dataset/29jan2022/BGC_1100_0100_all.nc",
    params=params,
    stn="lusk",
)
change_TElec = (
    dw.bgc.dataset[0]["interpol_value"]["value"]-
    dw.eclipse.dataset[0]["interpol_value"]["value"]
)
change_TIon = (
    dw.eclipse.dataset[1]["interpol_value"]["value"]
    -dw.bgc.dataset[1]["interpol_value"]["value"]
)
change_WI = (
    dw.eclipse.dataset[2]["interpol_value"]["value"]-
    dw.bgc.dataset[2]["interpol_value"]["value"]
)
alt = dw.eclipse.intp_height
time = dw.eclipse.dataset[0]["interpol_value"]["time"]


fig, axes = plt.subplots(dpi=200, figsize=(6, 6), nrows=2, ncols=1)
ax = axes[0]
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
ax.set_xlim(dt.datetime(2017,8,21,15), dt.datetime(2017,8,21,22))
ax.set_ylabel("Heights, km")
vlines=[
    dt.datetime(2017, 8, 21, 16, 15),
    dt.datetime(2017, 8, 21, 17, 45),
    dt.datetime(2017, 8, 21, 19, 30),
]
im = ax.pcolor(
    time, alt, change_TElec.T, vmax=100, vmin=-200, cmap="Spectral"
)
cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
cb = fig.colorbar(im, ax=ax, cax=cax)
cb.set_label(r"$T_{Elec}$, K")
for d in vlines:
    ax.axvline(d, ls="--", lw=0.4, color="k")
ax.set_xticks([])
ax.set_ylim(100, 400)

ax = axes[1]
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
ax.set_xlim(dt.datetime(2017,8,21,15), dt.datetime(2017,8,21,22))
ax.set_ylabel("Heights, km")
im = ax.pcolor(
    time, alt, change_TIon.T, vmax=50, vmin=-50, cmap="Spectral"
)
cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
cb = fig.colorbar(im, ax=ax, cax=cax)
cb.set_label(r"$T_{Ion}$, K")
for d in vlines:
    ax.axvline(d, ls="--", lw=0.4, color="k")
ax.set_xlabel("UT")
ax.set_ylim(100, 400)

# ax = axes[2]
# ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
# ax.set_xlim(dt.datetime(2017,8,21,15), dt.datetime(2017,8,21,22))
# ax.set_ylabel("Heights, km")
# im = ax.pcolor(
#     time, alt, change_WI.T, vmax=10, vmin=-10, cmap="Spectral"
# )
# cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
# cb = fig.colorbar(im, ax=ax, cax=cax)
# cb.set_label(r"$W_I$, m/s")
# for d in vlines:
#     ax.axvline(d, ls="--", lw=0.4, color="k")
# ax.set_xlabel("UT")
# ax.set_ylim(100, 400)

fig.subplots_adjust(wspace=0.1, hspace=0.1)
fig.savefig("dataset/figures/temp_vel_dist.png", bbox_inches="tight")