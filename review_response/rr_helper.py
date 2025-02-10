
import sys
sys.path.extend(["py/","pyrt/"])
from fetch import FetchModel, pconst

import aacgmv2
import datetime as dt
import numpy as np

from eutils import create_eclipse_path_local, read_eclispe_path, smooth
from pynasonde.digisonde.digi_utils import get_digisonde_info

def calculate_GC_stats_model(stn, hf1=150, hf2=240):
    k = pconst["q_e"]**2/(pconst["m_e"]*pconst["eps0"])
    params = [
        {
            "name": "e",
            "unit_multiplier": 1.,
            "unit": "cc",
            "mol_to_cc": "T",
            "interpolate": {"scale": "log", "type": "linear"},
        },
    ]
    fm = FetchModel(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        params=params,
    )
    fm.run_height_interpolate(stn)
    loc = fm.__fetch_latlon_by_station__(stn)
    loc["olat"], loc["olon"], _ = aacgmv2.get_aacgm_coord(loc["lat"], loc["lon"], 300, dt.datetime(2017,8,21,18))
    
    efm = fm.dataset[0]["interpol_value"]["value"]
    f2, f1 = (
        np.sqrt(efm[:, np.argmin(np.abs(fm.intp_height-hf2))]*1e6*k)/(2*np.pi), 
        np.sqrt(efm[:, np.argmin(np.abs(fm.intp_height-hf1))]*1e6*k)/(2*np.pi)
    )
    f_dif2 = (f2-f1)
    f_dif2_max = (abs(np.min(f_dif2))*(2*np.pi))**2/(k*1e6)
    f_dif2_dur = 5*np.sum(f_dif2 <= 0, axis=0)
    print(f_dif2_max, f_dif2_dur)

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import matplotlib.dates as mdates

    import matplotlib as mpl
    # mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    plt.rc('text', usetex=False)
    fig = plt.figure(figsize=(10, 5), dpi=300)
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.plot(
        fm.dataset[0]["interpol_value"]["time"],
        f2/1e6, color="b"
    )
    ax.plot(
        fm.dataset[0]["interpol_value"]["time"],
        f1/1e6, color="r"
    )
    ax.axis(
        xmin=dt.datetime(2017,8,21,15), 
        xmax=dt.datetime(2017,8,21,20)
    )
    # fig.savefig(f"review_response/{stn}.png")

    dates = [dt.datetime(2017,8,21,17)+dt.timedelta(minutes=i) for i in range(180)]
    p = create_eclipse_path_local(
        dates,
        loc["lat"], loc["lon"]
    )
    print(">>>>>>>>>>>>>>>>>>>>>",np.nanmin(p), np.nanmax(p))
    t = fm.dataset[0]["interpol_value"]["time"][np.argmax(p)]

    loc["time"] = fm.dataset[0]["interpol_value"]["time"]
    loc["f1"], loc["f2"] = f1/1e6, f2/1e6
    loc["hf1"], loc["hf2"] = hf1, hf2
    loc["n_dif2_max"] = f_dif2_max*3.5
    loc["f_dif2_dur"] = f_dif2_dur
    loc["o"] = 1-np.nanmax(p)
    loc["stn"] = stn
    loc["local_time"] = convert_to_local_time(t, loc["lat"], loc["lon"])
    loc["dist"] = distance_from_center(stn)
    return loc


def estimate_GCparams(df, ylim, xlim, stn_code, hf1=150, hf2=240):
    stn_info = get_digisonde_info(stn_code)
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
    df1, df2 = (
        df[
            (df.th>=hf1-5)
            & (df.th<=hf1+5)
        ],
        df[
            (df.th>=hf2-5)
            & (df.th<=hf2+5)
        ]
    )
    df1 = df1.groupby("datetime").max().reset_index()
    df2 = df2.groupby("datetime").max().reset_index()
    df1.ed, df2.ed = (
        smooth(df1.ed*1e6, 3), smooth(df2.ed*1e6, 3)
    )
    f_dif2 = np.array(df2.ed-df1.ed)
    f_dif2_max = abs(np.nanmin(f_dif2))
    f_dif2_dur = 5*np.sum(f_dif2 <= 0, axis=0)
    print(f_dif2, df1, len(df2))

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import matplotlib.dates as mdates
    import matplotlib as mpl
    # mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
    mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})
    plt.rc('text', usetex=False)
    fig = plt.figure(figsize=(10, 5), dpi=300)
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.plot(
        df2.datetime,
        df2.ed, color="b"
    )
    ax.plot(
        df1.datetime,
        df1.ed, color="r"
    )
    ax.axis(
        xmin=dt.datetime(2017,8,21,15), 
        xmax=dt.datetime(2017,8,21,20)
    )
    fig.savefig(f"review_response/{stn_code}.png")

    return dict(
        hmax=hmax, hmax_t=hmax_t, 
        p=p, f_dif2_max=f_dif2_max,
        f_dif2_dur=f_dif2_dur
    )


def lay_eclipse_occl(ax, time, ocl):
    import matplotlib.ticker as mticker
    import matplotlib.dates as mdates

    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    ax.plot(
        time, 1-ocl, color="k", ls="--", lw=1.2
    )
    ax.set_ylim(0,1)
    ax.set_yticks([])

    return ax


def convert_to_local_time(timestamp, latitude, longitude):
    from timezonefinder import TimezoneFinder
    from datetime import datetime
    import pytz
    tf = TimezoneFinder()
    timezone_str = tf.timezone_at(lng=longitude, lat=latitude)

    if timezone_str:
        timezone = pytz.timezone(timezone_str)
        utc_time = timestamp.replace(tzinfo=pytz.utc)
        local_time = utc_time.replace(tzinfo=pytz.utc).astimezone(timezone)
        return local_time
    else:
        return None  # Timezone not found for given coordinates

def distance_from_center(stn):
    stn_info = get_digisonde_info(stn)
    o = read_eclispe_path(2017)
    o = o[["LatC", "LonC"]]
    o = o[
        (o.LonC>=stn_info["LONG"]-2)
        & (o.LonC<stn_info["LONG"]+2)
    ]
    return o.LatC.mean()-stn_info["LAT"]
