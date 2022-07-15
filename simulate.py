import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])
import datetime as dt
import sys
sys.path.append("py/")
import os
from fetch2D import *

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
    }
]
times = [dt.datetime(2017,8,21,15)+ dt.timedelta(minutes=i*5) for i in range(7*12)]
h = 200
lat_res = 2
lon_res = 2
p_names = ["Op_CHM.ddif", "amb_diff.ddif", "dwind.ddif", "dfield.ddif"]
colorbar_labels = [r"$\delta(p-l), cm^{-3}s^{-1}$", 
                   r"$\delta(D_{\alpha}), cm^{-3}s^{-1}$",
                   r"$\delta(D_{wind}), cm^{-3}s^{-1}$", 
                   r"$\delta(D_{\vec{E}\times\vec{B}}), cm^{-3}s^{-1}$"]
for time in times:
    fname = "tmp/%02d-%02d.png"%(time.hour, time.minute)
    if not os.path.exists(fname):
        summary_plots_2D(params, time, h, lat_res, lon_res, p_names, 
                         colorbar_labels, fname)