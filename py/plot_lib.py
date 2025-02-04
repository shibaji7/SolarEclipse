#!/usr/bin/env python

"""plot_lib.py: Dedicated for utility functions."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import pickle

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from scipy import array
from scipy.interpolate import interp1d

import mplstyle
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
import matplotlib as mpl

mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize": 12, "font.size": 12})


def smooth(x, window_len=51, window="hanning"):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )
    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    if window == "flat":
        w = numpy.ones(window_len, "d")
    else:
        w = eval("np." + window + "(window_len)")
    y = np.convolve(w / w.sum(), s, mode="valid")
    d = window_len - 1
    y = y[int(d / 2) : -int(d / 2)]
    return y


def extrap1d(x, y, kind="linear"):
    """This method is used to extrapolate 1D paramteres"""
    interpolator = interp1d(x, y, kind=kind)
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        elif x > xs[-1]:
            return ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(list(map(pointwise, array(xs))))

    return ufunclike


def compute_1D_TS(
    axes,
    Hs=[150, 240],
    date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 20)],
    lab=["(a)", "(b)"],
    vlines=[
        dt.datetime(2017, 8, 21, 16, 15),
        dt.datetime(2017, 8, 21, 17, 45),
        dt.datetime(2017, 8, 21, 19, 30),
    ],
    stn="",
    prm="d_dif",
):
    """
    Plot 1D TS
    """
    files = [
        "latest_e.pickle",
        "latest_Op_CHMP.pickle",
        "latest_Op_CHML.pickle",
        "latest_dwind.pickle",
        "latest_dfield.pickle",
        "latest_amb_diff.pickle",
    ]
    params = ["e", "Op_CHMP","Op_CHML", "dwind", "dfield", "amb_diff"]
    idat, jdat = [], []
    for i, f in enumerate(files):
        with open("dataset/" + f, "rb") as handle:
            b = pickle.load(handle)
            time = [dt.datetime.strptime(tx, "%Y-%m-%d %H:%M") for tx in b["time"]]
            hs = b["Hs"]
            o = b[params[i] + ".d_dif"]
            idat.append(o[:, hs.tolist().index(150)])
            jdat.append(o[:, hs.tolist().index(240)])
    dat = np.array([idat, jdat])

    multipliers = [0.1, 1]
    for i, h in enumerate(Hs):
        ax = axes[i]
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
        ax.xaxis.set_major_locator(mdates.HourLocator())
        ax.set_xlabel("UT")
        ax.set_xlim(date_lim)
        ax.set_ylabel(
            r"$\delta [e/\frac{\partial O^+}{\partial t}]$, $cm^{-3}/cm^{-3}s^{-1}$"
        )
        ax.plot(
            time,
            dat[i][0] / 1000,
            "gray",
            lw=1.,
            ls="-",
            label=r"$\delta(e)\times 10^3$",
        )
        m = multipliers[i]
        ax.plot(time, m*(dat[i][1]+dat[i][2]), "r", lw=0.5, ls="-", label=r"$\delta(p-l)$")
        ax.plot(time, dat[i][3], "b", lw=0.5, ls="-", label=r"$\delta(D_{wind})$")
        ax.plot(
            time,
            dat[i][3],
            "k",
            lw=0.5,
            ls="-",
            label=r"$\delta(D_{\vec{E}\times\vec{B}})$",
        )
        ax.plot(
            time, dat[i][5], "darkgreen", lw=0.5, ls="-", label=r"$\delta(D_{\alpha})$"
        )
        ax.plot(
            time, m*(dat[i][1]+dat[i][2]+dat[i][3]+dat[i][4]+dat[i][5]),
            "cyan", lw=0.5, ls="-", label=r"$\sum\delta(\mu)$"
        )
        ax.axhline(0, ls="--", lw=0.3, alpha=0.4, color="k")
        ax.text(
            0.05,
            1.05,
            lab[i] + " h=%d km" % h,
            ha="left",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_ylim(-20, 20)
        if i == 1:
            ax.set_ylabel("")
            ax.set_yticklabels(["", "", "", "", ""])
        if i == 1:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
            ax.text(0.05, 0.9, stn, ha="left", va="center", transform=ax.transAxes)
        for d in vlines:
            ax.axvline(d, ls="--", lw=0.4, color="k")
    return


def compute_2D_TS_F1(
    axes,
    dwx,
    params,
    fig,
    date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 20)],
    lab=["(a)", "(b)"],
    vlines=[
        dt.datetime(2017, 8, 21, 16, 15),
        dt.datetime(2017, 8, 21, 17, 45),
        dt.datetime(2017, 8, 21, 19, 30),
    ],
    stn="",
):
    """
    Plot 2D histograms for the parameters.
    """
    dct = dwx.diffential_difference_2D()
    labels = [r"e^-", r"M^+"]
    time, Hs = dct["time"], dct["Hs"]
    l0, l1 = time.index(vlines[0] + dt.timedelta(minutes=60)), time.index(
        vlines[-1] - dt.timedelta(minutes=60)
    )
    ranges = [8000, 8000]
    dct["Mp.d_dif"] = dct["NOp.d_dif"] + dct["O2p.d_dif"] + dct["O2p.d_dif"]
    for i, p in enumerate(["e", "Mp"]):
        Ts = []
        o = dct[p + ".d_dif"]
        with open("dataset/latest_%s.pickle" % p, "wb") as handle:
            pickle.dump(
                {
                    "time": [tx.strftime("%Y-%m-%d %H:%M") for tx in time[1:]],
                    "Hs": Hs,
                    p + ".d_dif": o,
                },
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        for h in Hs:
            arg = np.argmin(np.abs(o[l0:l1, Hs.tolist().index(h)]))
            Ts.append(time[l0 + arg])
        ax = axes[i]
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
        ax.set_xlabel("UT")
        ax.set_xlim(date_lim)
        ax.set_ylabel("Heights, km")
        im = ax.pcolor(
            time[1:], Hs, o.T, vmax=ranges[i], vmin=-1 * ranges[i], cmap="jet"
        )
        if i == 0:
            Tx = smooth(np.array([(x - Ts[0]).total_seconds() for x in Ts]))
            Tx = [Ts[0] + dt.timedelta(seconds=x) for x in Tx]
            ax.plot(Tx, Hs, "w-", lw=1.5)
            ax.plot(Tx, Hs, "r--", lw=0.8)
        else:
            ax.axvline(vlines[1], color="w", ls="-", lw=1.5)
            ax.axvline(vlines[1], color="m", ls="--", lw=0.8)
        ax.text(
            0.05,
            1.10,
            lab[i] + r" $\delta(%s)$" % labels[i],
            ha="left",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_ylim(100, 300)
        for d in vlines:
            ax.axvline(d, ls="--", lw=0.4, color="k")
        if i == 1:
            cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
            cb = fig.colorbar(im, ax=ax, cax=cax)
            cb.set_label(r"$cm^{-3}$")
            ax.set_ylabel("")
            ax.set_yticklabels(["", "", "", "", ""])
    return


def compute_2D_TS(
    axes,
    dwx,
    params,
    fig,
    date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 20)],
    lab=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"],
    vlines=[
        dt.datetime(2017, 8, 21, 16, 15),
        dt.datetime(2017, 8, 21, 17, 45),
        dt.datetime(2017, 8, 21, 19, 30),
    ],
    stn="",
):
    """
    Plot 2D histograms for the parameters.
    """
    labels = ["p", "l", "p-l", r"D_{wind}", r"D_{\alpha}", r"D_{\vec{E}\times\vec{B}}"]
    dct = dwx.diffential_difference_2D()
    time, Hs = dct["time"], dct["Hs"]
    ranges = [100, 100, 10, 10, 5, 5]
    l0, l1 = time.index(vlines[0] + dt.timedelta(minutes=60)), time.index(
        vlines[-1] - dt.timedelta(minutes=60)
    )
    for i, p in enumerate(
        ["Op_CHMP", "Op_CHML", "Op_CHM", "dwind", "amb_diff", "dfield"]
    ):
        Ts = []
        o = dct[p + ".d_dif"]
        with open("dataset/latest_%s.pickle" % p, "wb") as handle:
            pickle.dump(
                {
                    "time": [tx.strftime("%Y-%m-%d %H:%M") for tx in time[1:]],
                    "Hs": Hs,
                    p + ".d_dif": o,
                },
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
        ax = axes[i]
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
        ax.set_xlabel("UT")
        ax.set_xlim(date_lim)
        ax.set_ylabel("Heights, km")
        im = ax.pcolor(
            time[1:], Hs, o.T, vmax=ranges[i], vmin=-1 * ranges[i], cmap="jet"
        )
        if i < 2:
            Ts = []
            for h in Hs:
                arg = np.argmin(np.abs(o[l0:l1, Hs.tolist().index(h)]))
                Ts.append(time[l0 + arg])
            Tx = smooth(np.array([(x - Ts[0]).total_seconds() for x in Ts]))
            Tx = [Ts[0] + dt.timedelta(seconds=x) for x in Tx]
            ax.plot(Tx, Hs, "w-", lw=1.5)
            ax.plot(Tx, Hs, "m--", lw=0.8)
            # if "Op" in p:
            #    for h in Hs:
            #        arg = np.argmin(np.abs(o[l0:l1, Hs.tolist().index(h)]))
            #        Ts.append(time[l0+arg])
            #    ax.plot(Ts, Hs, "k-", lw=0.8,)
        ax.text(
            0.05,
            1.10,
            lab[i] + r" $\delta(%s)$" % labels[i],
            ha="left",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_ylim(100, 300)
        for d in vlines:
            ax.axvline(d, ls="--", lw=0.4, color="k")
        if np.mod(i, 2) == 1:
            cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
            cb = fig.colorbar(im, ax=ax, cax=cax)
            cb.set_label(r"$cm^{-3}s^{-1}$")
            ax.set_ylabel("")
            ax.set_yticklabels(["", "", "", "", ""])
    return


def create_master_plot(
    stn="lusk",
):
    from fetch import DiffWACCMX

    params = [
        {
            "name": "NOp",
            "unit_multiplier": 1.0,
            "unit": "cc",
            "mol_to_cc": "T",
            "interpolate": {"scale": "log", "type": "linear"},
        },
        {
            "name": "O2p",
            "unit_multiplier": 1.0,
            "unit": "cc",
            "mol_to_cc": "T",
            "interpolate": {"scale": "log", "type": "linear"},
        },
        {
            "name": "N2p",
            "unit_multiplier": 1.0,
            "unit": "cc",
            "mol_to_cc": "T",
            "interpolate": {"scale": "log", "type": "linear"},
        },
        {
            "name": "e",
            "unit_multiplier": 1.0,
            "unit": "cc",
            "mol_to_cc": "T",
            "interpolate": {"scale": "log", "type": "linear"},
        },
    ]
    dw = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=stn,
    )
    fig, axes = plt.subplots(dpi=240, figsize=(6, 2.5), nrows=1, ncols=2)
    axes = axes.ravel()
    compute_2D_TS_F1(axes[:2], dw, params, fig)
    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.savefig("dataset/figures/density_distribution.png", bbox_inches="tight")
    del dw

    fig, axes = plt.subplots(dpi=200, figsize=(6, 8), nrows=3, ncols=2)
    axes = axes.ravel()
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
    dw = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        stn=stn,
    )
    compute_2D_TS(axes, dw, params, fig)
    fig.subplots_adjust(wspace=0.1, hspace=0.4)
    fig.savefig("dataset/figures/difference_distribution.png", bbox_inches="tight")
    del dw

    fig, axes = plt.subplots(dpi=200, figsize=(6, 2.5), nrows=1, ncols=2)
    axes = axes.ravel()
    compute_1D_TS(axes)
    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    fig.savefig("dataset/figures/onedim_distribution.png", bbox_inches="tight")
    return fig
