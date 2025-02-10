# from pynasonde.digisonde.dvl import DvlExtractor
from pynasonde.digisonde.digi_plots import SaoSummaryPlots
from pynasonde.digisonde.edp import EdpExtractor
from pynasonde.digisonde.sao import SaoExtractor


import pandas as pd
import numpy as np
import datetime as dt
import sys
sys.path.append("review_response")
from eutils import create_eclipse_path_local, smooth
from rr_helper import (
    calculate_GC_stats_model, estimate_GCparams,
    lay_eclipse_occl
)

import matplotlib.ticker as mticker
import matplotlib.dates as mdates

def download_possible_datasets(stations):
    from pynasonde.webhook import Webhook

    wh = Webhook()
    for stn_code in stations:
        sources = [
            dict(
                uri=f"https://data.ngdc.noaa.gov/instruments/remote-sensing/active/profilers-sounders/ionosonde/mids09/{stn_code}/individual/2017/233/scaled/",
                folder=f"/media/chakras4/Crucial X9/NOAA_Archives/profilers-sounders/ionosonde/mids09/{stn_code}/2017/233/scaled/",
            ),
            dict(
                uri=f"https://data.ngdc.noaa.gov/instruments/remote-sensing/active/profilers-sounders/ionosonde/mids09/{stn_code}/2017-264/image/",
                folder=f"/media/chakras4/Crucial X9/NOAA_Archives/profilers-sounders/ionosonde/mids09/{stn_code}/2017/233/image/",
            ),
        ]
        for source in sources:
            wh.__check_all_sub_folders__(
                source["uri"],
                source["folder"],
                ["SAO", "EDP", "PNG", "MMM", "16C"],
            )
    return


def generate_digisonde_pfh_profiles(
    folders,
    fig_file_name,
    fig_title="",
    draw_local_time=False,
):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    plt.rc('text', usetex=False)
    out = calculate_GC_stats_model(stn_code)
    df = SaoExtractor.load_SAO_files(
        folders=folders, func_name="height_profile",
        n_procs=12,
    )
    df.ed = df.ed / 1e6
    dfsc = SaoExtractor.load_SAO_files(
        folders=folders, func_name="scaled",
        n_procs=12,
    )
    sao_plot = SaoSummaryPlots(
        font_size=12,
        figsize=(6, 4 * 2),
        nrows=2,
        fig_title=fig_title,
        draw_local_time=draw_local_time,
    )
    xlim = [dt.datetime(2017,8,21,16), dt.datetime(2017,8,21,21)]
    ylim = [90, 300]
    gcp = estimate_GCparams(df, ylim, xlim, out["stn"])
    ax, _ = sao_plot.add_TS(
        df,
        zparam="ed",
        prange=[0, 0.3],
        ylim=ylim,
        zparam_lim=10,
        xlim=xlim,
        cbar_label=r"$N_e$,$\times 10^{6}$ /cc",
        plot_type="scatter",
        title=f"Stn Code: {out['stn']}",
        scatter_ms=300,
        xlabel="",  # if i == 2 else "",
        cmap="Blues",
    )
    ax.plot(
        dfsc.datetime, 
        dfsc.hmF1,
        "+", color="darkred", 
        ls="None", ms="5", zorder=4
    )
    lay_eclipse_occl(ax.twinx(), gcp["hmax_t"], gcp["p"])
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax = sao_plot.get_axes(False)
    ax.plot(
        out["time"], out["f1"],
        color="r", ls="-", lw=0.9,
        label="$foF_1$"
    )
    ax.plot(
        out["time"], out["f2"],
        color="b", ls="-", lw=0.9,
        label="$foF_2$"
    )
    lay_eclipse_occl(ax.twinx(), gcp["hmax_t"], gcp["p"])
    ax.legend(loc=4)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.set_xlim(xlim)
    ax.set_ylim(2, 6)
    txt = r"$\theta=%.1f^{\circ}$, $\phi=%.1f^{\circ}$"%(out["olat"], out["olon"])
    txt += "\n"
    txt += r" $\mathcal{O}=%.2f$, $\Delta T_{GC}$=%d min, $GC^p$=%.1f el/cc"%(out["o"], out["f_dif2_dur"], out["n_dif2_max"])
    ax.set_xlabel("Time, UT")
    ax.set_ylabel(r"$foF_{1,2}$, MHz")
    ax.text(
        0.05, 0.05, txt,
        ha="left", va="bottom",
        transform=ax.transAxes
    )
    sao_plot.save(fig_file_name)
    sao_plot.close()
    print(out["local_time"], out["dist"])
    print(gcp["f_dif2_dur"], gcp["f_dif2_max"])
    return


######################################
## Download all dataset 2017
######################################
# download_possible_datasets(["BC840", "AU930", "AL945", "WI937", "EG931"])

## Analyzing the dataset form 2017 Eclipse
stn_code = "EG931"
folders = [
    f"/media/chakras4/Crucial X9/NOAA_Archives/profilers-sounders/ionosonde/mids09/{stn_code}/2017/233/scaled/",
]
generate_digisonde_pfh_profiles(
    folders,
    f"dataset/figures/2017_{stn_code.lower()}_pf.png",
    fig_title="Digisondes / 21 August, 2017",
)