#!/usr/bin/env python

"""plot_lib.py: module is dedicated to plot and create the movies."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from loguru import logger

sys.path.extend(["py/"])


plt.style.use(["science", "ieee"])


class OverlayDataset(object):
    """
    Ovealay SD and SM data
    """

    def __init__(self, figsize=(8, 8), nrows=2, ncols=2, dpi=150, coords="geo"):
        logger.info("Fan plot initialization")
        self.fig = plt.figure(figsize=figsize, dpi=dpi)
        self.nrows = nrows
        self.ncols = ncols
        self._num_subplots_created = 0
        self.coords = coords
        return

    def _add_axis(self, date, h):
        # from sd_carto import SDCarto
        from sd_carto import SDCarto
        self.frm = ccrs.Geodetic()
        self.to = ccrs.NorthPolarStereo(-90, 90)
        self._num_subplots_created += 1
        logger.info(f"Adding axes to fan plot : {self._num_subplots_created}")
        pos = int(str(self.nrows) + str(self.ncols) + str(self._num_subplots_created))
        ax = self.fig.add_subplot(
            pos,
            projection="sdcarto",
            map_projection=self.to,
            plot_date=date,
            coords=self.coords,
        )
        ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        ax.set_extent([-180, 180, 20, 90], crs=ccrs.PlateCarree())
        plt_lons = np.arange(-180, 181, 30)
        mark_lons = np.arange(-180, 180, 30)
        plt_lats = np.arange(20, 90, 15)
        gl = ax.gridlines(ccrs.PlateCarree(), linewidth=0.5)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        ax.mark_latitudes(plt_lats, fontsize="small", color="darkblue")
        ax.mark_longitudes(plt_lons, fontsize="small", color="darkblue")
        #ax.date_string(h)
        return ax

    def plot_data(self, df, date, h, p_name, colorbar_label):
        """Plot data - overlay on top of map"""
        ax = self._add_axis(date, h)
        ax.overlay_data(
            df,
            self.to,
            self.frm,
            self.fig,
            p_name,
            colorbar_label,
        )
        ax.add_dn_terminator()
        return ax

    def save(self, filepath):
        logger.info(f"Save FoV figure to : {filepath}")
        self.fig.subplots_adjust(wspace=0.7, hspace=0.1)
        self.fig.savefig(filepath, bbox_inches="tight")
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return
