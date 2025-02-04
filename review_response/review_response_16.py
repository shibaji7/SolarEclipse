import sys
sys.path.extend(["py/"])
from fetch import DiffWACCMX

import matplotlib.ticker as mticker
import numpy as np
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.dates as mdates
import datetime as dt

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import matplotlib as mpl
mpl.rcParams.update({"xtick.labelsize": 5, "ytick.labelsize":5, "font.size":5})

mpl.rcParams.update({"xtick.labelsize":6, "ytick.labelsize":6, "font.size":6})


stn = "lusk"
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
dct = dw.diffential_difference_data(
    kind="TS", hs=[150, 240, 300, 350]
)
print(dct.keys())
fig, ax = plt.subplots(dpi=200, figsize=(5, 3), nrows=1, ncols=1)
ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
ax.xaxis.set_major_locator(mdates.HourLocator())
ax.set_xlim([
    dt.datetime(2017,8,21,15),
    dt.datetime(2017,8,21,21)
])
# ax.set_ylim(-3,3)
ax.plot(
    dct["Op_CHML_150"].time, 
    dct["Op_CHMP_150"].d_dif+dct["Op_CHML_150"].d_dif,
    color="r", lw=0.5, ls="-", 
)
# ax.axvline(
#     dct["Op_CHML_150"].time.tolist()[
#         np.argmax(
#             dct["Op_CHML_150"].d_dif
#         )
#     ],
#     ls="--", color="r"
# )
ax.plot(
    dct["Op_CHML_240"].time, 
    dct["Op_CHMP_240"].d_dif+dct["Op_CHML_240"].d_dif,
    color="b", lw=0.5, ls="-", 
)
# ax.axvline(
#     dct["Op_CHML_240"].time.tolist()[
#         np.argmax(
#             dct["Op_CHML_240"].d_dif
#         )
#     ],
#     ls="--", color="b"
# )
# ax.plot(
#     dct["Op_CHML_300"].time, 
#     dct["Op_CHML_300"].d_dif,
#     color="k", lw=0.5, ls="-", 
# )
# ax.axvline(
#     dct["Op_CHML_300"].time.tolist()[
#         np.argmax(
#             dct["Op_CHML_300"].d_dif
#         )
#     ],
#     ls="--", color="k"
)
# ax.plot(
#     dct["Op_CHMP_350"].time, 
#     dct["Op_CHML_350"].d_dif,
#     color="g", lw=0.5, ls="-", 
# )
ax.set_xlabel("UT")
fig.savefig("dataset/figures/height.png", bbox_inches="tight")
