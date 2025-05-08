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
#from read_rad import Radar
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

def draw_images(lats, lons, date, p, I, i=0, to=ccrs.Orthographic(-90, -90), cb=True):
    fig = plt.figure(dpi=180, figsize=figsize)
    ax = plt.axes(projection=to)
    ax.set_global()
    ax.add_feature(Nightshade(date, alpha=0.3))
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
    return ax, fig

def _add_colorbar(fig, ax, im, colormap="gray_r", label=r"Obscuration, ($\mathcal{O}$)"):
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
        print("Time %1.2f (s)"%(t))
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

def fetch_event(date, lats, lons, i, to, rads, fname=None, cb=True):
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

    ax, fig = draw_images(Lats, Lons, date, p[0,0,:,:], None, i, to, cb=cb)
    overlay_instrument(ax, dict(lat=39.992, lon=-105.269, name="Boulder"), to=to, mult=1)
    overlay_instrument(ax, dict(lat=42.75, lon=-104.455, name="Lusk"), to=to, mult=-1)
    overlay_instrument(ax, dict(lat=42.6, lon=-71.5, name="MHISR"), to=to, mult=-1, markerColor="r")
    overlay_instrument(ax, dict(lat=37.94, lon=284.42 - 360, name="WI937"), to=to, mult=1, markerColor="m")
    overlay_instrument(ax, dict(lat=40.00, lon=254.70 - 360, name="BC840"), to=to, mult=-.2, markerColor="m", markerSize=3)
    overlay_instrument(ax, dict(lat=30.40, lon=262.30 - 360, name="AU930"), to=to, mult=-1, markerColor="m")

    # Totality path
    if fname:
        LatC, LonC = [], []
        LatN, LonN = [], []
        LatS, LonS = [], []
        with open(fname, "r") as f: lines = f.readlines()
        for line in lines:
            line = line.split("  ")
            locN, loc, locS = line[1], line[3], line[2]

            latcomp = -1 if "S" in loc else 1
            loc = loc.split(" ")
            LatC.append(
                latcomp*float(loc[0])+
                float(loc[1].\
                    replace(".","").\
                    replace("N","").\
                    replace("S",""))/1e3
            )
            LonC.append(
                -1*float(loc[2])+
                float(loc[3].\
                    replace(".","").\
                    replace("W",""))/1e3
            )

            locS = locS.split(" ")
            LatS.append(
                latcomp*float(locS[0])+
                float(locS[1].\
                    replace(".","").\
                    replace("N","").\
                    replace("S",""))/1e3
            )
            LonS.append(
                -1*float(locS[2])+
                float(locS[3].\
                    replace(".","").\
                    replace("W",""))/1e3
            )

            locN = locN.split(" ")
            LatN.append(
                latcomp*float(locN[0])+
                float(locN[1].\
                    replace(".","").\
                    replace("N","").\
                    replace("S",""))/1e3
            )
            LonN.append(
                -1*float(locN[2])+
                float(locN[3].\
                    replace(".","").\
                    replace("W",""))/1e3
            )
        LatC, LonC = smooth(np.array(LatC))+0.4, smooth(np.array(LonC))
        LatS, LonS = smooth(np.array(LatS))+0.4, smooth(np.array(LonS))
        LatN, LonN = smooth(np.array(LatN))+0.4, smooth(np.array(LonN))
        xyz = to.transform_points(ccrs.PlateCarree(), n.array(LonC), n.array(LatC))
        x, y = xyz[:, 0], xyz[:, 1]
        ax.plot(x, y, ls=":", color="k", lw=0.7, alpha=1.)
        xyz = to.transform_points(ccrs.PlateCarree(), n.array(LonS), n.array(LatS))
        x, y = xyz[:, 0], xyz[:, 1]
        ax.plot(x, y, ls=":", color="r", lw=0.7, alpha=1.)
        xyz = to.transform_points(ccrs.PlateCarree(), n.array(LonN), n.array(LatN))
        x, y = xyz[:, 0], xyz[:, 1]
        ax.plot(x, y, ls=":", color="r", lw=0.7, alpha=1.)
    ax.set_extent([-120, -60, 25, 55], crs = ccrs.PlateCarree())
    fig.savefig("dataset/figures/%03d.png"%i, bbox_inches="tight")
    return ax, fig


if __name__ == "__main__":
    date = dt.datetime(2017,8,21,18,26,40)
    lats = n.linspace(0,90,num=90*2)
    lons = n.linspace(-180,180,num=91*2)
    i = 10
    to = ccrs.Orthographic(-90, 50)
    fetch_event(date, lats, lons, i, to, [], "database/2017.csv")