#!/usr/bin/env python

"""utils.py: Dedicated for utility functions."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from scipy import array
from scipy.interpolate import interp1d
import pickle

plt.style.use(["science", "ieee"])


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


class SummaryPlots(object):
    """
    Create a summary plot.
    """

    def __init__(self, dwx=None, isr=None, digi=None):
        """
        Parameters:
        -----------
        dwx: WACCMX model
        """
        self.dwx = dwx
        self.isr = isr
        self.digi = digi
        return

    def compute_1D_TS(
        self,
        Hs=[150, 240],
        date_lim=[dt.datetime(2017, 8, 21, 14), dt.datetime(2017, 8, 22)],
        lab=["(a)", "(b)"],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
        stn="",
        prm="d_dif",
    ):
        """
        Plot 1D TS
        """
        dct = self.dwx.diffential_difference_data(kind="TS", hs=Hs)
        fig = plt.figure(dpi=300, figsize=(6, 2))
        for i, h in enumerate(Hs):
            oChmP, oChmL = dct["Op_CHMP_" + str(h)], dct["Op_CHML_" + str(h)]
            oWind, oAmb, oField = (
                dct["dwind_" + str(h)],
                dct["amb_diff_" + str(h)],
                dct["dfield_" + str(h)],
            )
            oElec = dct["e_" + str(h)]
            ax = fig.add_subplot(121 + i)
            ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
            ax.set_xlabel("UT")
            ax.set_xlim(date_lim)
            ax.set_ylabel(
                r"$\delta [\frac{\partial O^+}{\partial t}]$, $cm^{-3}s^{-1}$"
            )
            ax.plot(
                oElec.time,
                oElec[prm],
                "k",
                lw=0.5,
                ls="--",
                label=r"$\delta(e)$"
            )
            ax.plot(
                oChmP.time,
                (oChmP[prm] - oChmL[prm]),
                "r",
                lw=0.5,
                ls="-",
                label=r"$\delta(p-l)$",
            )
            ax.plot(
                oWind.time,
                oWind[prm],
                "b",
                lw=0.5,
                ls="-",
                label=r"$\delta(D_{wind})$",
            )
            ax.plot(
                oAmb.time,
                oAmb[prm],
                "k",
                lw=0.5,
                ls="-",
                label=r"$\delta(D_{\alpha})$",
            )
            ax.plot(
                oField.time,
                oField[prm],
                "darkgreen",
                lw=0.5,
                ls="-",
                label=r"$\delta(D_{\vec{E}\times\vec{B}})$",
            )
            #ax.set_ylim(-20, 20)
            ax.axhline(0, ls="--", lw=0.3, alpha=0.4, color="k")
            ax.text(
                0.05,
                1.05,
                lab[i] + " h=%d km" % h,
                ha="left",
                va="center",
                transform=ax.transAxes,
            )
            if i == 0:
                ax.legend(loc=1, fontsize=6)
                ax.text(0.05, 0.9, stn, ha="left", va="center", transform=ax.transAxes)
            for d in vlines:
                ax.axvline(d, ls="--", lw=0.4, color="k")
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        return fig

    def compute_1D_d_TS(
        self,
        Hs=[120, 150, 240, 300],
        date_lim=[dt.datetime(2017, 8, 21, 14), dt.datetime(2017, 8, 22)],
        lab=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
        stn="",
        prm="d_dif",
    ):
        dct = self.dwx.diffential_difference_data(kind="TS", hs=Hs)
        fig = plt.figure(dpi=300, figsize=(6, 4))
        for i, h in enumerate(Hs):
            oOp, oO2p, oNOp, oN2p = (
                dct["Op_" + str(h)],
                dct["O2p_" + str(h)],
                dct["NOp_" + str(h)],
                dct["N2p_" + str(h)],
            )
            ax = fig.add_subplot(221 + i)
            ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
            ax.set_xlabel("UT")
            ax.set_xlim(date_lim)
            ax.set_ylabel(r"$\delta [N^+]$, $cm^{-3}$")
            ax.plot(
                oOp.time,
                oOp[prm],
                "r",
                lw=0.5,
                ls="-",
                label=r"$\delta(O^+)$",
            )
            ax.plot(
                oO2p.time,
                oO2p[prm],
                "b",
                lw=0.5,
                ls="-",
                label=r"$\delta(O_2^+)$",
            )
            ax.plot(
                oN2p.time,
                oN2p[prm],
                "k",
                lw=0.5,
                ls="-",
                label=r"$\delta(N_2^+)$",
            )
            ax.plot(
                oNOp.time,
                oNOp[prm],
                "g",
                lw=0.5,
                ls="-",
                label=r"$\delta(NO^+)$",
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
            #             del_t = (
            #                 oe.time.tolist()[np.argmin(oe[prm])]
            #                 - oOp.time.tolist()[np.argmin(oOp[prm])]
            #             ).total_seconds()
            #             ax.text(
            #                 0.99,
            #                 1.05,
            #                 "Delay: %d s" % (del_t),
            #                 ha="right",
            #                 va="center",
            #                 transform=ax.transAxes,
            #             )
            if i == 0:
                ax.legend(loc=1, fontsize=6)
                ax.text(0.05, 0.9, stn, ha="left", va="center", transform=ax.transAxes)
            for d in vlines:
                ax.axvline(d, ls="--", lw=0.4, color="k")
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        return fig

    def compute_1D_pl_TS(
        self,
        Hs=[120, 150, 240, 300],
        date_lim=[dt.datetime(2017, 8, 21, 14), dt.datetime(2017, 8, 22)],
        lab=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
        stn="",
        prm="d_dif",
    ):
        """
        1D plot of p-l from composition
        """
        dct = self.dwx.diffential_difference_data(kind="TS", hs=Hs)
        fig = plt.figure(dpi=300, figsize=(6, 4))
        for i, h in enumerate(Hs):
            oChmP, oChmL = (dct["Op_CHMP_" + str(h)], dct["Op_CHML_" + str(h)])
            ax = fig.add_subplot(221 + i)
            ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
            ax.set_xlabel("UT")
            ax.set_xlim(date_lim)
            ax.set_ylabel(r"$\delta [p/l]$, $cm^{-3}s^{-1}$")
            ax.plot(
                oChmP.time,
                oChmP[prm],
                "r",
                lw=0.5,
                ls="-",
                label=r"$\delta(p)$",
            )
            ax.plot(
                oChmL.time,
                oChmL[prm],
                "b",
                lw=0.5,
                ls="-",
                label=r"$\delta(l)$",
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
            del_t = (
                oChmP.time.tolist()[np.argmin(oChmP[prm])]
                - oChmL.time.tolist()[np.argmin(oChmL[prm])]
            ).total_seconds()
            ax.text(
                0.99,
                1.05,
                "Delay: %d s" % (del_t),
                ha="right",
                va="center",
                transform=ax.transAxes,
            )
            if i == 0:
                ax.legend(loc=1, fontsize=6)
                ax.text(0.05, 0.9, stn, ha="left", va="center", transform=ax.transAxes)
            for d in vlines:
                ax.axvline(d, ls="--", lw=0.4, color="k")
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        return fig

    def compute_2D_TS(
        self,
        params,
        date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 20)],
        lab=["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
        stn="",
    ):
        """
        Plot 2D histograms for the parameters.
        """
        labels = ["p", "l", "p-l", r"D_{wind}", r"D_{\alpha}", r"D_{\vec{E}\times\vec{B}}"]
        fig = plt.figure(dpi=300, figsize=(6, 6))
        dct = self.dwx.diffential_difference_2D()
        time, Hs = dct["time"], dct["Hs"]
        ranges = [80, 80, 20, 5, 5, 5]
        l0, l1 = time.index(vlines[0] + dt.timedelta(minutes=60)), time.index(
            vlines[-1] - dt.timedelta(minutes=60)
        )
        for i, p in enumerate(["Op_CHMP", "Op_CHML", "Op_CHM", "dwind", "amb_diff", "dfield"]):
            Ts = []
            o = dct[p + ".d_dif"]
            with open("dataset/latest_%s.pickle"%p, "wb") as handle:
                pickle.dump(
                    {
                        "time": [tx.strftime("%Y-%m-%d %H:%M") for tx in time[1:]],
                        "Hs": Hs,
                        p + ".d_dif": o
                    }, 
                    handle, 
                    protocol=pickle.HIGHEST_PROTOCOL
                )
            ax = fig.add_subplot(321 + i)
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
                ax.plot(Tx, Hs, "m-", lw=0.8)
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
            cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
            cb = fig.colorbar(im, ax=ax, cax=cax)
            cb.set_label(r"$cm^{-3}s^{-1}$")
        fig.subplots_adjust(wspace=0.8, hspace=0.8)
        return fig

    def compute_2D_TS_F1(
        self,
        params,
        date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 20)],
        lab=["(a)", "(b)"],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
        stn="",
    ):
        """
        Plot 2D histograms for the parameters.
        """
        labels = [r"e^-", r"M^+"]
        fig = plt.figure(dpi=300, figsize=(6, 2))
        dct = self.dwx.diffential_difference_2D()
        time, Hs = dct["time"], dct["Hs"]
        l0, l1 = time.index(vlines[0] + dt.timedelta(minutes=60)), time.index(
            vlines[-1] - dt.timedelta(minutes=60)
        )
        ranges = [8000, 8000]
        dct["Mp.d_dif"] = dct["NOp.d_dif"] + dct["O2p.d_dif"] + dct["O2p.d_dif"]
        for i, p in enumerate(["e", "Mp"]):
            Ts = []
            o = dct[p + ".d_dif"]
            with open("dataset/latest_%s.pickle"%p, "wb") as handle:
                pickle.dump(
                    {
                        "time": [tx.strftime("%Y-%m-%d %H:%M") for tx in time[1:]],
                        "Hs": Hs,
                        p + ".d_dif": o
                    }, 
                    handle, 
                    protocol=pickle.HIGHEST_PROTOCOL
                )
            for h in Hs:
                arg = np.argmin(np.abs(o[l0:l1, Hs.tolist().index(h)]))
                Ts.append(time[l0 + arg])
            ax = fig.add_subplot(121 + i)
            ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
            ax.set_xlabel("UT")
            ax.set_xlim(date_lim)
            ax.set_ylabel("Heights, km")
            im = ax.pcolor(
                time[1:], Hs, o.T, vmax=ranges[i], vmin=-1 * ranges[i], cmap="jet"
            )
            Tx = smooth(np.array([(x - Ts[0]).total_seconds() for x in Ts]))
            Tx = [Ts[0] + dt.timedelta(seconds=x) for x in Tx]
            ax.plot(Tx, Hs, "m-", lw=0.8)
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
            if i==1:
                cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
                cb = fig.colorbar(im, ax=ax, cax=cax)
                cb.set_label(r"$cm^{-3}$")
                ax.set_ylabel("")
                ax.set_yticklabels(["", "", "", "", ""])
        fig.subplots_adjust(wspace=0.1, hspace=0.6)
        return fig

    def compute_2D_TS_Tne(
        self,
        params,
        date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 20)],
        lab=["(a)", "(b)"],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
        stn="",
    ):
        """
        Plot 2D histograms for the parameters.
        """
        labels = [r"T_n", r"T_e"]
        fig = plt.figure(dpi=300, figsize=(3, 3))
        ax = fig.add_subplot(111)
        dct = self.dwx.diffential_difference_2D()
        time, Hs = dct["time"], dct["Hs"]
        for i, p in enumerate(["T", "TElec"]):
            o = dct[p + ".d_dif"]
            ax = fig.add_subplot(121 + i)
            ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
            ax.set_xlabel("UT")
            ax.set_xlim(date_lim)
            ax.set_ylabel("Heights, km")
            im = ax.pcolor(
                time[1:], Hs, o.T, vmax=ranges[i], vmin=-1 * ranges[i], cmap="jet_r"
            )
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
            cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
            cb = fig.colorbar(im, ax=ax, cax=cax)
            cb.set_label(r"$K$")
        fig.subplots_adjust(wspace=0.8, hspace=0.8)
        return fig

    def compute_1D_del(
        self,
        date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 22)],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
    ):
        labels = [r"e^-", "p", "l"]
        fig = plt.figure(dpi=300, figsize=(6, 3))
        dct = self.dwx.diffential_difference_2D()
        time, Hs = dct["time"], dct["Hs"]
        ax = fig.add_subplot(111)
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
        ax.set_xlabel("UT")
        ax.set_xlim(date_lim)
        ax.set_ylabel("Heights, km")
        l0, l1 = time.index(vlines[0] + dt.timedelta(minutes=60)), time.index(
            vlines[-1] - dt.timedelta(minutes=60)
        )
        l2, l3 = time.index(vlines[1] + dt.timedelta(minutes=60)), time.index(
            vlines[-1] + dt.timedelta(minutes=180)
        )
        colors = ["r", "b", "g", "k"]
        for i, p in enumerate(["e", "Op_CHMP", "Op_CHML"]):
            Ts0, Ts1 = [], []
            o = dct[p + ".d_dif"]
            for h in Hs:
                arg = np.argmin(np.abs(o[l0:l1, Hs.tolist().index(h)]))
                Ts0.append(time[l0 + arg])
                arg = np.argmin(np.abs(o[l2:l3, Hs.tolist().index(h)]))
                Ts1.append(time[l0 + arg])
            Tx = smooth(np.array([(x - Ts0[0]).total_seconds() for x in Ts0]))
            Tx = [Ts0[0] + dt.timedelta(seconds=x) for x in Tx]
            ax.plot(
                Tx,
                Hs,
                color=colors[i],
                ls="-",
                lw=0.6,
                label=r" $\delta(%s)$" % labels[i],
            )
            Tx = smooth(np.array([(x - Ts1[0]).total_seconds() for x in Ts1]))
            Tx = [Ts1[0] + dt.timedelta(seconds=x) for x in Tx]
            ax.plot(Tx, Hs, color=colors[i], ls="--", lw=0.6)
        for d in vlines:
            ax.axvline(d, ls="--", lw=0.4, color="k")
        ax.legend(loc=2)
        ax.set_ylim(100, 300)
        fig.subplots_adjust(wspace=0.5, hspace=0.5)
        return

def compute_1D_TS(
    Hs=[150, 240],
    date_lim=[dt.datetime(2017, 8, 21, 15), dt.datetime(2017, 8, 21, 20)],
    lab=["(i)", "(j)"],
    vlines=[
        dt.datetime(2017, 8, 21, 16),
        dt.datetime(2017, 8, 21, 17, 45),
        dt.datetime(2017, 8, 21, 20),
    ],
    stn="",
    prm="d_dif",
):
    """
    Plot 1D TS
    """
    files = [
        "latest_e.pickle",
        "latest_Op_CHM.pickle",
        "latest_dwind.pickle",
        "latest_dfield.pickle",
        "latest_amb_diff.pickle",
    ]
    params = ["e", "Op_CHM", "dwind", "dfield", "amb_diff"]
    idat, jdat = [], []
    for i, f in enumerate(files):
        with open("dataset/"+f, "rb") as handle:
            b = pickle.load(handle)
            time = [dt.datetime.strptime(tx, "%Y-%m-%d %H:%M") for tx in b["time"]]
            hs = b["Hs"]
            o = b[params[i]+".d_dif"]
            idat.append(o[:, hs.tolist().index(150)])
            jdat.append(o[:, hs.tolist().index(240)])
    dat = np.array([idat, jdat])
    fig = plt.figure(dpi=300, figsize=(6, 2))
    for i, h in enumerate(Hs):
        ax = fig.add_subplot(121 + i)
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
        ax.xaxis.set_major_locator(mdates.HourLocator())
        ax.set_xlabel("UT")
        ax.set_xlim(date_lim)
        ax.set_ylabel(
            r"$\delta [e/\frac{\partial O^+}{\partial t}]$, $cm^{-3}/cm^{-3}s^{-1}$"
        )
        ax.plot(
            time,
            dat[i][0]/1000,
            "gray",
            lw=0.5,
            ls="-",
            label=r"$\delta(e)\times 10^3$"
        )
        ax.plot(
            time,
            dat[i][1],
            "r",
            lw=0.5,
            ls="-",
            label=r"$\delta(p-l)$"
        )
        ax.plot(
            time,
            dat[i][2],
            "b",
            lw=0.5,
            ls="-",
            label=r"$\delta(D_{wind})$"
        )
        ax.plot(
            time,
            dat[i][3],
            "k",
            lw=0.5,
            ls="-",
            label=r"$\delta(D_{\vec{E}\times\vec{B}})$"
        )
        ax.plot(
            time,
            dat[i][4],
            "darkgreen",
            lw=0.5,
            ls="-",
            label=r"$\delta(D_{\alpha})$"
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
        ax.set_ylim(-20,20)
        if i==1: 
            ax.set_ylabel("")
            ax.set_yticklabels(["", "", "", "", ""])
        if i == 1:
            ax.legend(bbox_to_anchor=(1.01, 1.05), fontsize=6)
            ax.text(0.05, 0.9, stn, ha="left", va="center", transform=ax.transAxes)
        for d in vlines:
            ax.axvline(d, ls="--", lw=0.4, color="k")
    fig.subplots_adjust(wspace=0.1, hspace=0.6)
    return fig

    