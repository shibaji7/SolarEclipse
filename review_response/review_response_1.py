# from pynasonde.digisonde.dvl import DvlExtractor
from pynasonde.digisonde.digi_plots import SaoSummaryPlots
from pynasonde.digisonde.edp import EdpExtractor
from pynasonde.digisonde.sao import SaoExtractor
from pynasonde.digisonde.digi_utils import get_digisonde_info

import pandas as pd
import numpy as np
import datetime as dt
import sys
sys.path.append("review_response")
from utils import create_eclipse_path_local, smooth


def estimate_GCparams(df, ylim, xlim, stn_code):
    stn_info = get_digisonde_info(stn_code)
    print(stn_info)
    df = df[
        (df.th>=ylim[0]) &
        (df.th<=ylim[1]) &
        (df.datetime>=xlim[0]) &
        (df.datetime<=xlim[1]) 
    ]
    hmax = []
    hmax_t = []
    for x in df.datetime.unique():
        o = df[df.datetime==x]
        if len(o) > 0:
            hmax.append(o.th.tolist()[o.ed.argmax()])
            hmax_t.append(x)
    
    hmax_t = [
        dt.datetime.utcfromtimestamp(
            (x-np.datetime64('1970-01-01T00:00:00'))/np.timedelta64(1, 's')
        ) 
        for x in hmax_t
    ]
    p = create_eclipse_path_local(
        hmax_t, 
        stn_info["LAT"], stn_info["LONG"]
    )
    return dict(hmax=hmax, hmax_t=hmax_t, p=p)

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
    func_name,
    fig_file_name,
    fig_title="",
    draw_local_time=False,
    stns=[],
):
    dfs = [
        SaoExtractor.load_SAO_files(
            folders=[folder], func_name=func_name,
            n_procs=12,
        )
        for folder in folders
    ]
    dfs_scaled = [
        SaoExtractor.load_SAO_files(
            folders=[folder], func_name="scaled",
            n_procs=12,
        )
        for folder in folders
    ]
    N = len(folders)
    sao_plot = SaoSummaryPlots(
        font_size=12,
        figsize=(6, 4 * N),
        nrows=N,
        fig_title=fig_title,
        draw_local_time=draw_local_time,
    )
    xlim = [dt.datetime(2017,8,21,15), dt.datetime(2017,8,21,21)]
    ylim = [90, 300]
    for i in range(N):
        df = dfs[i]
        dfsc = dfs_scaled[i]
        df.ed = df.ed / 1e6
        gcp = estimate_GCparams(df, ylim, xlim, stns[i])
        ax, _ = sao_plot.add_TS(
            df,
            zparam="ed",
            prange=[0, 0.5],
            ylim=ylim,
            zparam_lim=10,
            xlim=xlim,
            cbar_label=r"$N_e$,$\times 10^{6}$ /cc",
            plot_type="scatter",
            title="Stn Code: " + stns[i],
            scatter_ms=300,
            xlabel="Time, UT",  # if i == 2 else "",
            cmap="Blues",
        )
        ax.plot(
            dfsc.datetime, 
            dfsc.hmF1,
            "+", color="darkred", 
            ls="None", ms="5", zorder=4
        )
        axc = ax.twinx()
        axc.plot(
            gcp["hmax_t"],
            1-gcp["p"], ls="-", 
            lw=0.8, color="k"
        )
        axc.set_yticks([])
    sao_plot.save(fig_file_name)
    sao_plot.close()
    return


def generate_digisonde_hNmf_profiles(
    folders,
    func_name,
    fig_file_name,
    fig_title="",
    draw_local_time=False,
    stns=[],
):
    dfs = [
        SaoExtractor.load_SAO_files(folders=[folder], func_name=func_name)
        for folder in folders
    ]
    N = len(folders)
    sao_plot = SaoSummaryPlots(
        font_size=12,
        figsize=(6, 4 * N),
        nrows=N,
        fig_title=fig_title,
        draw_local_time=draw_local_time,
    )

    for i in range(N):
        df = dfs[i]
        ax, _ = sao_plot.plot_TS(
            df,
            left_yparams=["foF1"],
            left_ylim=[1, 15],
            right_ylim=[80, 400],
            right_yparams=["hmF1"],
            title="Stn Code: " + stns[i],
        )

    sao_plot.save(fig_file_name)
    sao_plot.close()
    return


def generate_digisonde_edp_profiles(
    folders,
    func_name,
    fig_file_name,
    fig_title="",
    draw_local_time=False,
    stns=[],
):
    dfs = [
        EdpExtractor.load_EDP_files(folders=[folder], func_name=func_name)
        for folder in folders
    ]
    N = len(folders)
    sao_plot = SaoSummaryPlots(
        font_size=12,
        figsize=(6, 4 * N),
        nrows=N,
        fig_title=fig_title,
        draw_local_time=draw_local_time,
    )

    for i in range(N):
        df = dfs[i]
        print(df.density.min(),df.density.max())
        df.density = df.density/1e12
        ax, _ = sao_plot.add_TS(
            df,
            yparam="height",
            zparam="density",
            prange=[0, 0.5],
            ylim=[90, 600],
            zparam_lim=1e6,
            cbar_label=r"$N_e$,$\times 10^{12}$ /cm",
            plot_type="scatter",
            title="Stn Code: " + stns[i],
            scatter_ms=10,
            xlabel="Time, UT",  # if i == 2 else "",
        )
    sao_plot.save(fig_file_name)
    sao_plot.close()
    return


######################################
## Download all dataset 2017
######################################
# download_possible_datasets(["BC840", "AU930", "AL945", "WI937"])

## Analyzing the dataset form 2017 Eclipse
stn_code = "WI937"
folders = [
    f"/media/chakras4/Crucial X9/NOAA_Archives/profilers-sounders/ionosonde/mids09/{stn_code}/2017/233/scaled/",
]
func_name = "height_profile"
generate_digisonde_pfh_profiles(
    folders,
    func_name,
    f"dataset/figures/2017_{stn_code.lower()}_pf.png",
    fig_title="Digisondes / 21 August, 2017",
    stns=[stn_code],
)

# generate_digisonde_hNmf_profiles(
#     folders,
#     "scaled",
#     f"dataset/figures/2017_{stn_code.lower()}_hNmf.png",
#     fig_title="Digisondes / 21 August, 2017",
#     stns=[stn_code],
# )

# drift_dataset = DvlExtractor.load_DVL_files(folders=[])


# generate_digisonde_edp_profiles(
#     folders,
#     func_name,
#     f"dataset/figures/2017_{stn_code.lower()}_edp.png",
#     fig_title="Digisondes / 21 August, 2017",
#     stns=[stn_code],
# )
