import os
import ephem
import numpy as n
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                   "Lucida Grande", "Verdana"]
import matplotlib as mpl
mpl.rcParams.update({"xtick.labelsize": 12, "ytick.labelsize":12, "font.size":12})

import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
import datetime as dt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import shapefile
#from read_rad import Radar
import pydarn
import numpy as np

figsize = (6,6)
scalesize = 0.4
def magnetic_inclination_angle(lats, lons, date):
    import pyIGRF
    FACT = 180.0 / np.pi
    r = 6371.0
    I = np.zeros_like(lats)*np.nan
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            _, ix, _, _, _, _, _ = pyIGRF.igrf_value(
                lats[i, j], lons[i, j], r, date.year
            )
            I[i,j] = ix
    return I

def overlay_instrument(ax, inst, to, tx=ccrs.PlateCarree(), marker="D", zorder=20, markerColor="b", markerSize=12, 
    fontSize=11, font_color="darkgreen", xOffset=-2, yOffset=-2, annotate=True, mult=1):
    lat, lon = inst["lat"], inst["lon"]
    #lat *= mult
    ax.scatter([lon], [lat], s=markerSize, marker=marker,
        color=markerColor, zorder=zorder, transform=tx, lw=0.8, alpha=0.4)
    if annotate:
        x, y = to.transform_point(lon+xOffset, lat+yOffset*mult, src_crs=tx)
        ax.text(x, y, inst["name"], ha="center", va="center", transform=to,
                     fontdict={"color":font_color, "size":fontSize}, alpha=0.8)
    return

def overlay_radar(ax, rad, to, tx=ccrs.PlateCarree(), marker="D", zorder=2, markerColor="k", markerSize=2, 
    fontSize="small", font_color="darkgreen", xOffset=-5, yOffset=-1.5, annotate=True, mult=1):
    """ Adding the radar location """
    hdw = pydarn.read_hdw_file(rad)
    lat, lon = hdw.geographic.lat, hdw.geographic.lon
    lat *= mult
    ax.scatter([lon], [lat], s=markerSize, marker=marker,
        color=markerColor, zorder=zorder, transform=tx, lw=0.8, alpha=0.4)
    nearby_rad = [["adw", "kod", "cve", "fhe", "wal", "gbr", "pyk", "aze", "sys"],
                ["ade", "ksr", "cvw", "fhw", "bks", "sch", "sto", "azw", "sye"]]
    if annotate:
        if rad in nearby_rad[0]: xOff, ha = -5 if not xOffset else -xOffset, -2
        elif rad in nearby_rad[1]: xOff, ha = 5 if not xOffset else xOffset, -2
        else: xOff, ha = xOffset, -1
        x, y = to.transform_point(lon+xOff, lat+ha, src_crs=tx)
        ax.text(x, y, rad.upper(), ha="center", va="center", transform=to,
                     fontdict={"color":font_color, "size":fontSize}, alpha=0.8)
    return

def overlay_fov(ax, rad, to, tx=ccrs.PlateCarree(), maxGate=60, rangeLimits=None, beamLimits=None,
                    model="IS", fov_dir="front", fovColor=None, fovAlpha=0.2,
                    fovObj=None, zorder=1, lineColor="k", lineWidth=0.6, ls="-"):
    """ Overlay radar FoV """
    from numpy import transpose, ones, concatenate, vstack, shape
    hdw = pydarn.read_hdw_file(rad)
    sgate = 0
    egate = hdw.gates if not maxGate else maxGate
    ebeam = hdw.beams
    if beamLimits is not None: sbeam, ebeam = beamLimits[0], beamLimits[1]
    else: sbeam = 0
    fov = pydarn.Coords.GEOGRAPHIC(hdw.stid)
    latFull, lonFull = fov[0].T, fov[1].T
    xyz = to.transform_points(tx, lonFull, latFull)
    x, y = xyz[:, :, 0], xyz[:, :, 1]
    contour_x = concatenate((x[sbeam, sgate:egate], x[sbeam:ebeam, egate],
                x[ebeam, egate:sgate:-1],
                x[ebeam:sbeam:-1, sgate]))
    contour_y = concatenate((y[sbeam, sgate:egate], y[sbeam:ebeam, egate],
            y[ebeam, egate:sgate:-1],
            y[ebeam:sbeam:-1, sgate]))
    ax.plot(contour_x, contour_y, color=lineColor, zorder=zorder, linewidth=lineWidth, ls=ls, alpha=0.6)
    return

def draw_images(lats, lons, date, p, I, i=0, to=ccrs.Orthographic(-90, -90), cb=True, 
    rads=[], beam_sound={}):
    fig = plt.figure(dpi=180, figsize=figsize)
    ax = plt.axes(projection=to)
    ax.set_global()
    #ax.add_feature(Nightshade(date, alpha=0.3))
    ax.coastlines()
    p = np.ma.masked_invalid(p)
    im = ax.contourf(
        lons,
        lats,
        p.T,
        transform=ccrs.PlateCarree(),
        cmap="gray_r", alpha=0.6,
        levels=[0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0]
    )
    if cb: _add_colorbar(fig, ax, im)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.3, 
        color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.xlocator = mticker.FixedLocator(n.arange(-180,180,60))
    gl.ylocator = mticker.FixedLocator(n.arange(-90,90,30))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    for rad in rads:
        overlay_radar(ax, rad, to)
        overlay_fov(ax, rad, to, beamLimits=[0,16])
        overlay_fov(ax, rad, to, beamLimits=[beam_sound["beam"], beam_sound["beam"]+1], 
            lineColor=beam_sound["color"])
    return ax, fig

def _add_colorbar(fig, ax, im, colormap="gray_r", label=r"Occultation, ($\alpha$)"):
        """Add a colorbar to the right of an axis."""
        pos = ax.get_position()
        cpos = [
            pos.x1 + 0.1,
            pos.y0 + 0.2 * pos.height,
            0.02,
            pos.height * scalesize,
        ]  # this list defines (left, bottom, width, height)
        cax = fig.add_axes(cpos)
        cb2 = fig.colorbar(
            im,
            cax=cax,
            cmap=colormap,
            spacing="uniform",
            orientation="vertical",
        )
        cb2.set_label(label)
        return

def intersection(r0,r1,d,n_s=100):
    A1=n.zeros([n_s,n_s])
    A2=n.zeros([n_s,n_s])
    I=n.zeros([n_s,n_s])
    x=n.linspace(-2.0*r0,2.0*r0,num=n_s)
    y=n.linspace(-2.0*r0,2.0*r0,num=n_s)
    xx,yy=n.meshgrid(x,y)
    A1[n.sqrt((xx+d)**2.0+yy**2.0) < r0]=1.0
    n_sun=n.sum(A1)
    A2[n.sqrt(xx**2.0+yy**2.0) < r1]=1.0
    S=A1+A2
    I[S>1]=1.0
    eclipse=n.sum(I)/n_sun
    return(eclipse)

def get_eclipse(t0,n_t,dt=60.0,alts=n.linspace(0,600e3,num=600),lats=[38.0],lons=[-120]):
    # Location
    obs = ephem.Observer()
    n_alts=len(alts)
    n_lats=len(lats)
    n_lons=len(lons)
    
    p=n.zeros([n_t,n_alts,n_lats,n_lons])
    times=n.arange(n_t)*dt
    dts=[]
    for ti,t in enumerate(times):
        print("Time %1.2f (s)"%(t), t0)
        for ai,alt in enumerate(alts):
            for lai,lat in enumerate(lats):
                for loi,lon in enumerate(lons):
                    #obs.lon, obs.lat = '-1.268304', '51.753101'#'16.02', '78.15' # ESR
                    obs.lon, obs.lat = '%1.2f'%(lon), '%1.2f'%(lat) # ESR
                    obs.elevation=alt
                    obs.date= t0#(ephem.date(ephem.date(t0)+t*ephem.second))
                    sun, moon = ephem.Sun(), ephem.Moon()
                    
                    # Output list
                    results=[]
                    seps=[]
                    sun.compute(obs)
                    moon.compute(obs)
                    r_sun=(sun.size/2.0)/3600.0
                    r_moon=(moon.size/2.0)/3600.0
                    s=n.degrees(ephem.separation((sun.az, sun.alt), (moon.az, moon.alt)))
                    percent_eclipse=0.0
                            
                    if s < (r_moon+r_sun):
#                        print("eclipsed")
                        if s < 1e-3:
                            percent_eclipse=1.0
                        else:
                            percent_eclipse=intersection(r_moon,r_sun,s,n_s=100)

                    if n.degrees(sun.alt) <= r_sun:
                        if n.degrees(sun.alt) <= -r_sun:
                            percent_eclipse=n.nan
                        else:
                            percent_eclipse=1.0-((n.degrees(sun.alt)+r_sun)/(2.0*r_sun))*(1.0-percent_eclipse)
            
                    p[ti,ai,lai,loi]=percent_eclipse
        dts.append(obs.date)
    return(p,times,dts)


def smooth(x,window_len=51,window="hanning"):
    if x.ndim != 1: raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len: raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3: return x
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]: raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == "flat": w = np.ones(window_len,"d")
    else: w = eval("np."+window+"(window_len)")
    y = np.convolve(w/w.sum(),s,mode="valid")
    d = window_len - 1
    y = y[int(d/2):-int(d/2)]
    return y

def load_path(ax, fname = "database/2024-path/center.shp", color="k"):
    shape = shapefile.Reader(fname)
    if len(shape.shapeRecords()) == 1:
        feature = shape.shapeRecords()[0]
        first = feature.shape.__geo_interface__    
        coords = first["coordinates"]
        lats, lons = (
            [d[1] for d in list(coords)],
            [d[0] for d in list(coords)]
        )
        xyz = to.transform_points(ccrs.PlateCarree(), n.array(lons), n.array(lats))
        x, y = xyz[:, 0], xyz[:, 1]
        ax.plot(x, y, ls=":", color=color, lw=0.7, alpha=1.)
    else:
        latsU, lonsU = ([], [])
        latsL, lonsL = ([], [])
        for feature in shape.shapeRecords():
            first = feature.shape.__geo_interface__
            coords = first["coordinates"]
    return

def fetch_event(date, lats, lons, i, to, rads, fname=None, cb=True, beam_sound={}):
    t0 = ephem.date(
        (
            date.year,
            date.month,
            date.day,
            date.hour,
            date.minute,
            date.second,
        )
    )
    p,_,_ = get_eclipse(t0, 
        n_t = 1,
        alts = n.array([100]), 
        lats = lats, 
        lons = lons,
    )
    p[p<.15] = 0
    Lats, Lons = n.meshgrid(lats, lons)
    # Overlay Inc
    #I = magnetic_inclination_angle(Lats, Lons, date)
    ax, fig = draw_images(Lats, Lons, date, p[0,0,:,:], None, i, to, cb=cb, rads=rads, beam_sound=beam_sound)
    overlay_instrument(ax, dict(lat=39.992, lon=-105.269, name="Boulder"), to=to, mult=1)
    overlay_instrument(ax, dict(lat=42.75, lon=-104.455, name="Lusk"), to=to, mult=-1)
    overlay_instrument(ax, dict(lat=42.6, lon=-71.5, name="MHISR"), to=to, mult=-1, markerColor="r")

    # Totality path
    load_path(ax, "database/2024-path/center.shp", "k")
    #load_path(ax, "database/2024-path/ppath01.shp", "r")
    
    ax.set_extent([-120, -60, 10, 80], crs = ccrs.PlateCarree())
    ax.text(
        0.99,
        1.05,
        date.strftime("%H:%M:%S UT") +\
            f" / {beam_sound['freq']} / {beam_sound['beam']}" ,
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontdict={"color":beam_sound["color"]}
    )
    fig.savefig("tmp/%04d.png"%i, bbox_inches="tight")
    return ax, fig


if __name__ == "__main__":
    mode = "mode1"
    beams=[0, 3, 7, 11, 15]
    freqs=["f0","f1","f2","f3"]
    colors=["darkred","darkgreen","darkblue","m"]
    beam_sound = 3
    runtime = int(120*60/beam_sound)
    dates = [
        dt.datetime(2024,4,8,18) + dt.timedelta(seconds=beam_sound*i)
        for i in range(runtime)
    ]
    lats = n.linspace(0,90,num=90*2)
    lons = n.linspace(-180,180,num=91*2)
    to = ccrs.Orthographic(-90, 50)
    for i, date in enumerate(dates):
        # Mode 1 sounding 
        beam_sound = dict(
            beam = beams[np.mod(i,len(beams))],
            freq = freqs[np.mod(int(i/len(beams)),len(freqs))],
            color = colors[np.mod(int(i/len(beams)),len(colors))],
            date=date
        )
        if not os.path.exists("tmp/%04d.png"%i):
            fetch_event(date, lats, lons, i, to, ["fhe"], "database/2024.csv", beam_sound=beam_sound)
            plt.close()
        #break
    cmd = f"ffmpeg -framerate 2 -pattern_type glob -i 'tmp/*.png' -c:v libx264 tmp/{mode}.mp4"
    os.system(cmd)
    cmd = f"ffmpeg -framerate 10 -pattern_type glob -i 'tmp/*.png' tmp/{mode}.gif"
    os.system(cmd)