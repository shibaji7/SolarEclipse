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
        Hs=[120, 150, 240, 300],
        factors=[1, 1, 1, 1],
        date_lim=[dt.datetime(2017, 8, 21, 14), dt.datetime(2017, 8, 21, 22)],
        lab=["(a)", "(b)", "(c)", "(d)"],
        vlines=[
            dt.datetime(2017, 8, 21, 16),
            dt.datetime(2017, 8, 21, 17, 45),
            dt.datetime(2017, 8, 21, 20),
        ],
        stn="",
    ):
        dct = self.dwx.fetch_data(kind="TS", hs=Hs)
        fig = plt.figure(dpi=300, figsize=(6, 4))
        for i, h in enumerate(Hs):
            oChmP, oChmL = dct["Op_CHMP_" + str(h)], dct["Op_CHML_" + str(h)]
            oWind, oAmb, oField = (
                dct["dwind_" + str(h)],
                dct["amb_diff_" + str(h)],
                dct["dfield_" + str(h)],
            )
            ax = fig.add_subplot(221 + i)
            ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
            ax.set_xlabel("UT")
            ax.set_xlim(date_lim)
            ax.set_ylabel(
                r"$\delta [\frac{\partial O^+}{\partial t}]\times %.1f$, $cm^{-3}s^{-1}$"
                % factors[i]
            )
            ax.plot(
                oChmP.time,
                (oChmP.dif - oChmL.dif) * factors[i],
                "r",
                lw=0.5,
                ls="-",
                label=r"$\delta(p-l)$",
            )
            ax.plot(
                oWind.time,
                oWind.dif * factors[i],
                "b",
                lw=0.5,
                ls="-",
                label=r"$\delta(D_{wind})$",
            )
            ax.plot(
                oAmb.time,
                oAmb.dif * factors[i],
                "k",
                lw=0.5,
                ls="-",
                label=r"$\delta(D_{\alpha})$",
            )
            ax.plot(
                oField.time,
                oField.dif * factors[i],
                "darkgreen",
                lw=0.5,
                ls="-",
                label=r"$\delta(D_{\vec{E}\times\vec{B}})$",
            )
            ax.set_ylim(-20, 20)
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

    def compute_2D_TS(self, hs=[120, 150, 240, 300]):
        return
