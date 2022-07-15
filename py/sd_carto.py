#!/usr/bin/env python

"""sd_carto.py: utility module for Costom Carto py geoaxes to plot data on aacgmv2 coordinates."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys

import cartopy
import matplotlib
import numpy as np
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.projections import register_projection
from shapely.geometry import LineString, MultiLineString, mapping

sys.path.extend(["py/"])


def get_gridded_parameters(
    q, xparam="lon", yparam="lat", zparam="", r=0, rounding=True
):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[[xparam, yparam, zparam]]
    if rounding:
        plotParamDF.loc[:, xparam] = np.round(plotParamDF[xparam].tolist(), r)
        plotParamDF.loc[:, yparam] = np.round(plotParamDF[yparam].tolist(), r)
    else:
        plotParamDF[xparam] = plotParamDF[xparam].tolist()
        plotParamDF[yparam] = plotParamDF[yparam].tolist()
    plotParamDF = plotParamDF.groupby([xparam, yparam]).mean().reset_index()
    plotParamDF = plotParamDF[[xparam, yparam, zparam]].pivot(xparam, yparam)
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y = np.meshgrid(x, y)
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
        np.isnan(plotParamDF[zparam].values), plotParamDF[zparam].values
    )
    return X, Y, Z


class SDCarto(GeoAxes):
    name = "sdcarto"

    def __init__(self, *args, **kwargs):
        self.supported_coords = ["geo", "aacgmv2", "aacgmv2_mlt"]
        if "coords" in kwargs and kwargs["coords"] is not None:
            self.coords = kwargs.pop("coords")
            if self.coords not in self.supported_coords:
                err_str = "Coordinates not supported, choose from : "
                for _n, _sc in enumerate(self.supported_coords):
                    if _n + 1 != len(self.supported_coords):
                        err_str += _sc + ", "
                    else:
                        err_str += _sc
                raise TypeError(err_str)
        else:
            self.coords = "geo"
        if "map_projection" in kwargs and kwargs["map_projection"] is not None:
            self.map_projection = kwargs.pop("map_projection")
        else:
            self.map_projection = cartopy.crs.PlateCarree()
        if "rad" in kwargs and kwargs["rad"] is not None:
            self.rad = kwargs.pop("rad")
        if "plot_date" in kwargs and kwargs["plot_date"] is not None:
            self.plot_date = kwargs.pop("plot_date")
        else:
            if self.coords == "aacgmv2" or self.coords == "aacgmv2_mlt":
                raise TypeError(
                    "Need to provide a date using 'plot_date' keyword for aacgmv2 plotting"
                )
        super().__init__(map_projection=self.map_projection, *args, **kwargs)
        return

    def overaly_coast_lakes(self, resolution="50m", color="black", **kwargs):
        """Overlay AACGM coastlines and lakes"""
        kwargs["edgecolor"] = color
        kwargs["facecolor"] = "none"
        kwargs["linewidth"] = 0.3
        # overaly coastlines
        feature = cartopy.feature.NaturalEarthFeature(
            "physical", "coastline", resolution, **kwargs
        )
        self.add_feature(cartopy.feature.COASTLINE, **kwargs)
        self.add_feature(cartopy.feature.LAKES, **kwargs)
        return

    def coastlines(self, resolution="50m", color="black", **kwargs):
        # details!
        kwargs["edgecolor"] = color
        kwargs["facecolor"] = "none"
        kwargs["linewidth"] = 0.3
        feature = cartopy.feature.NaturalEarthFeature(
            "physical", "coastline", resolution, **kwargs
        )
        return self.add_feature(feature, **kwargs)

    def add_feature(self, feature, **kwargs):
        if "edgecolor" not in kwargs:
            kwargs["edgecolor"] = "black"
        if self.coords == "geo":
            super().add_feature(feature, **kwargs)
        else:
            aacgm_geom = self.get_aacgm_geom(feature)
            aacgm_feature = cartopy.feature.ShapelyFeature(
                aacgm_geom, cartopy.crs.Geodetic(), **kwargs
            )
            kwargs["facecolor"] = "none"
            super().add_feature(aacgm_feature, **kwargs)
        return

    def add_dn_terminator(self, **kwargs):
        """Adding day night terminator"""
        from cartopy.feature.nightshade import Nightshade

        if self.plot_date:
            ns_feature = Nightshade(self.plot_date, alpha=0.2)
            super().add_feature(ns_feature, **kwargs)
        return

    def data_lay(self, df, rf, p_name, to, fm, idx=None, gflg_mask=None):
        """
        Data scatter-plots
        """
        df = df[df.slist <= self.maxGate]
        o = (
            df[df[gflg_mask] == idx]
            if (idx is not None) and (gflg_mask is not None)
            else df.copy()
        )
        if len(o) > 0:
            lons, lats = rf.lonFull, rf.latFull
            Xb, Yg, Px = get_gridded_parameters(
                o, xparam="bmnum", yparam="slist", zparam=p_name, rounding=True
            )
            Xb, Yg = Xb.astype(int), Yg.astype(int)
            lons, lats = lons[Xb.ravel(), Yg.ravel()].reshape(Xb.shape), lats[
                Xb.ravel(), Yg.ravel()
            ].reshape(Xb.shape)
            XYZ = to.transform_points(fm, lons, lats)
            Px = np.ma.masked_invalid(Px)
            return XYZ[:, :, 0], XYZ[:, :, 1], Px
        else:
            return [], [], []

    def overlay_data(
        self,
        df,
        to,
        fm,
        fig,
        p_name="Op_CHM_150.ddif",
        colorbar_label="",
        cmap="RdBu",
        vlim=[-20, 20],
        **kwargs
    ):
        """
        Adding radar data
        df: dataframe object
        """
        cmap = (
            matplotlib.pyplot.get_cmap("jet_r")
            if cmap is None
            else matplotlib.pyplot.get_cmap(cmap)
        )
        Xb, Yg, Px = get_gridded_parameters(
            df, xparam="lon", yparam="lat", zparam=p_name, rounding=True
        )
        XYZ = to.transform_points(fm, Xb, Yg)
        X, Y = XYZ[:, :, 0], XYZ[:, :, 1]
        im = self.pcolor(
            X, Y, Px.T, transform=to, cmap=cmap, vmin=vlim[0], vmax=vlim[1], **kwargs
        )
        cax = self.inset_axes([1.1, 0.1, 0.05, 0.5], transform=self.transAxes)
        cb = fig.colorbar(im, ax=self, cax=cax)
        cb.set_label(colorbar_label)
        return

    #######################
    # Additional
    #######################

    def add_feature(self, feature, **kwargs):
        # Now we"ll set facecolor as None because aacgm doesn"t close
        # continents near equator and it turns into a problem
        if "edgecolor" not in kwargs:
            kwargs["edgecolor"] = "black"
        if "facecolor" in kwargs:
            print(
                "manually setting facecolor keyword to none as aacgm fails for fill! want to know why?? think about equator!"
            )
        kwargs["facecolor"] = "none"
        if self.coords == "geo":
            super().add_feature(feature, **kwargs)
        else:
            aacgm_geom = self.get_aacgm_geom(feature)
            aacgm_feature = cartopy.feature.ShapelyFeature(
                aacgm_geom, cartopy.crs.Geodetic(), **kwargs
            )
            super().add_feature(aacgm_feature, **kwargs)

    def get_aacgm_geom(self, feature, out_height=300.0):
        new_i = []
        # cartopy.feature.COASTLINE
        for _n, i in enumerate(feature.geometries()):
            aa = mapping(i)
            mag_list = []
            geo_coords = aa["coordinates"]
            for _ni, _list in enumerate(geo_coords):
                mlon_check_jump_list = []
                split_mag_list = None
                if len(_list) == 1:
                    _loop_list = _list[0]
                else:
                    _loop_list = _list
                for _ngc, _gc in enumerate(_loop_list):
                    _mc = aacgmv2.get_aacgm_coord(
                        _gc[1], _gc[0], out_height, self.plot_date
                    )
                    if np.isnan(_mc[0]):
                        continue
                    mlon_check_jump_list.append(_mc[1])
                    if self.coords == "aacgmv2":
                        mag_list.append((_mc[1], _mc[0]))
                    else:
                        if _mc[2] * 15.0 > 180.0:
                            mag_list.append((_mc[2] * 15.0 - 360.0, _mc[0]))
                        else:
                            mag_list.append((_mc[2] * 15.0, _mc[0]))
                # check for unwanted jumps
                mlon_check_jump_list = np.array(mlon_check_jump_list)

                jump_arr = np.diff(mlon_check_jump_list)
                bad_inds = np.where(np.abs(jump_arr) > 10.0)[0]
                # delete the range of bad values
                # This is further complicated because
                # in some locations mlon jumps from -177 to +178
                # and this causes jumps in the maps! To deal with
                # this we"ll split arrays of such jumps
                # (these jumps typically have just one bad ind )
                # and make them into two seperate entities (LineStrings)
                # so that shapely will treat them as two seperate boundaries!
                if len(bad_inds) > 0:
                    if len(bad_inds) > 1:
                        mag_list = [
                            i
                            for j, i in enumerate(mag_list)
                            if j - 1 not in np.arange(bad_inds[0], bad_inds[1])
                        ]
                    else:
                        split_mag_list = mag_list[bad_inds[0] + 1 :]
                        mag_list = mag_list[: bad_inds[0] + 1]
                mag_coords = tuple(mag_list)
                if len(mag_list) > 1:
                    new_i.append(mag_coords)
                if split_mag_list is not None:
                    #             print(split_mag_list)
                    if len(split_mag_list) > 1:
                        new_i.append(tuple(split_mag_list))

        aacgm_coast = MultiLineString(new_i)
        return aacgm_coast

    def mark_latitudes(self, lat_arr, lon_location=45, **kwargs):
        """
        mark the latitudes
        Write down the latitudes on the map for labeling!
        we are using this because cartopy doesn"t have a
        label by default for non-rectangular projections!
        """
        if isinstance(lat_arr, list):
            lat_arr = np.array(lat_arr)
        else:
            if not isinstance(lat_arr, np.ndarray):
                raise TypeError("lat_arr must either be a list or numpy array")
        # make an array of lon_location
        lon_location_arr = np.full(lat_arr.shape, lon_location)
        proj_xyz = self.projection.transform_points(
            cartopy.crs.Geodetic(), lon_location_arr, lat_arr
        )
        # plot the lats now!
        out_extent_lats = False
        for _np, _pro in enumerate(proj_xyz[..., :2].tolist()):
            # check if lats are out of extent! if so ignore them
            lat_lim = self.get_extent(crs=cartopy.crs.Geodetic())[2::]
            if (lat_arr[_np] >= min(lat_lim)) and (lat_arr[_np] <= max(lat_lim)):
                self.text(
                    _pro[0], _pro[1], r"$%s^{\circ}$" % str(lat_arr[_np]), **kwargs
                )
            else:
                out_extent_lats = True
        if out_extent_lats:
            print("some lats were out of extent ignored them")

    def mark_longitudes(self, lon_arr=np.arange(-180, 180, 60), **kwargs):
        """
        mark the longitudes
        Write down the longitudes on the map for labeling!
        we are using this because cartopy doesn"t have a
        label by default for non-rectangular projections!
        This is also trickier compared to latitudes!
        """
        if isinstance(lon_arr, list):
            lon_arr = np.array(lon_arr)
        else:
            if not isinstance(lon_arr, np.ndarray):
                raise TypeError("lat_arr must either be a list or numpy array")
        # get the boundaries
        [x1, y1], [x2, y2] = self.viewLim.get_points()
        bound_lim_arr = []
        right_bound = LineString(([-x1, y1], [x2, y2]))
        top_bound = LineString(([x1, -y1], [x2, y2]))
        bottom_bound = LineString(([x1, y1], [x2, -y2]))
        left_bound = LineString(([x1, y1], [-x2, y2]))
        plot_outline = MultiLineString(
            [right_bound, top_bound, bottom_bound, left_bound]
        )
        # get the plot extent, we"ll get an intersection
        # to locate the ticks!
        plot_extent = self.get_extent(cartopy.crs.Geodetic())
        line_constructor = lambda t, n, b: np.vstack(
            (np.zeros(n) + t, np.linspace(b[2], b[3], n))
        ).T
        for t in lon_arr[:-1]:
            xy = line_constructor(t, 30, plot_extent)
            # print(xy)
            proj_xyz = self.projection.transform_points(
                cartopy.crs.Geodetic(), xy[:, 0], xy[:, 1]
            )
            xyt = proj_xyz[..., :2]
            ls = LineString(xyt.tolist())
            locs = plot_outline.intersection(ls)
            if not locs:
                continue
            # we need to get the alignment right
            # so get the boundary closest to the label
            # and plot it!
            closest_bound = min(
                [
                    right_bound.distance(locs),
                    top_bound.distance(locs),
                    bottom_bound.distance(locs),
                    left_bound.distance(locs),
                ]
            )
            if closest_bound == right_bound.distance(locs):
                ha = "left"
                va = "top"
            elif closest_bound == top_bound.distance(locs):
                ha = "left"
                va = "bottom"
            elif closest_bound == bottom_bound.distance(locs):
                ha = "left"
                va = "top"
            else:
                ha = "right"
                va = "top"
            if self.coords == "aacgmv2_mlt":
                marker_text = str(int(t / 15.0))
            else:
                marker_text = r"$%s^{\circ}$" % str(t)
            self.text(
                locs.bounds[0] + 0.05 * locs.bounds[0],
                locs.bounds[1] + 0.05 * locs.bounds[1],
                marker_text,
                ha=ha,
                va=va,
                **kwargs
            )
        return

    def to_aagcm(self, lat, lon):
        if "aacgmv2" in self.coords:
            lat, lon, mlt = aacgmv2.get_aacgm_coord(lat, lon, 300, self.plot_date)
            if self.coords == "aacgmv2_mlt":
                lon = mlt * 15
        return lat, lon

    def to_aagcms(self, lats, lons):
        mlats, mlons = np.zeros_like(lats), np.zeros_like(lats)
        if "aacgmv2" in self.coords:
            for i in range(lats.shape[0]):
                mlats[i, :], mlons[i, :], mlt = aacgmv2.get_aacgm_coord_arr(
                    lats[i, :], lons[i, :], 300, self.plot_date
                )
                if self.coords == "aacgmv2_mlt":
                    mlons[i, :] = mlt * 15
        else:
            mlats, mlons = lats, lons
        return mlats, mlons

    def date_string(self, h, label_style="web"):
        # Set the date and time formats
        dfmt = "%d/%b/%Y" if label_style == "web" else "%d %b %Y,"
        tfmt = "%H:%M"
        date_str = "{:{dd} {tt}} UT, ".format(self.plot_date, dd=dfmt, tt=tfmt)
        date_str += " h=%d km" % h
        self.text(
            0.5,
            0.95,
            date_str,
            ha="center",
            va="center",
            transform=self.transAxes,
            fontsize="medium",
        )
        return


register_projection(SDCarto)
