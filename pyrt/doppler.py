#!/usr/bin/env python3

"""doppler.py: simulate python program for Doppler computation along the path of ray"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import numpy as np
import pandas as pd
from scipy.interpolate import interp2d

class Doppler(object):

    def __init__(
        self,
        grange,
        height, 
        ne_eclipse,
        ne_no_eclipse,
    ):
        self.kconst, self.cconst = 80.6, 3e8
        self.grange, self.height, self.ne_eclipse, self.ne_no_eclipse = (
            grange, height, ne_eclipse, ne_no_eclipse
        )
        self.func_eclipse = interp2d(
            self.grange, self.height,
            np.log10(self.ne_eclipse)
        )
        self.func_no_eclipse = interp2d(
            self.grange,
            self.height,
            np.log10(self.ne_no_eclipse)
        )
        return
        
    def compute_doppler_f_from_ray(
        self,
        ray,
        elv,
        frequency
    ):
        k_const, c_const, delt = 80.6, 3e8, 60
        dop, dth, sth = [], [], []
        gi, hi = np.array(ray.grange), np.array(ray.height)
        dth.append(0.)
        sth.append(0.)
        for k in range(len(gi[:-1])):
            dth.append(np.abs(hi[k]-hi[k+1]))
            sth.append(np.sqrt((hi[k]-hi[k+1])**2+(gi[k]-gi[k+1])**2))
        dth, sth = np.array(dth)*1000., np.array(sth)*1000
        for k, g, h in zip(range(len(gi)), gi, hi):
            if h > 50.:
                dne = (10**self.func_eclipse(g,h) - 10**self.func_no_eclipse(g, h))[0]
                #df = (kconst / (cconst * self.frequency*1e6)) * (dne / delt) * (dth[k] / np.cos(np.deg2rad(90.-elvi)))
                df = (k_const / (c_const * frequency*1e6)) * (dne / delt) *1e6 * (dth[k])
                if np.isnan(df):  df = 0.
                dop.append(df)
            else: dop.append(0.)
        ray["dop"] = 1e3*np.array(dop) / np.cos(np.deg2rad(90.-elv))
        ray["sth"] = sth
        ray["dth"] = dth
        return ray

    def compute_doppler_f_from_file(
        self,
        date,
        rad,
        beam,
        elv,
        fold_base = "simulation_results/{date}/{rad}/",
        fname_base = "{ut}.{beam}.elv({e}).csv"
    ):
        folder = fold_base.format(
            date=date.strftime("%Y.%m.%d"),
            rad=rad
        )
        fname = folder + fname_base.format(
            ut=date.strftime("%H%M"),
            beam=beam,
            e=elv if "." in str(elv) else int(elv)
        )
        return