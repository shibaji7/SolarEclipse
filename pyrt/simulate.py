#!/usr/bin/env python3

"""simulate.py: simulate python program for RT"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import argparse
import copy
import datetime as dt
import glob
import os

import numpy as np
import pandas as pd
import pydarn
from dateutil import parser as dparser
from geopy.distance import great_circle as GC
from gitm import GITM
from doppler import Doppler
from loguru import logger
from scipy.io import savemat


def read_params(fname="rt.json"):
    import json
    from types import SimpleNamespace

    with open(fname, "r") as f:
        param = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
    return param


class RayTrace(object):
    """
    Ray trace class to trace all the points
    """

    def __init__(
        self,
        event,
        rad,
        beam,
        cfg,
    ):
        self.beam = beam
        self.event = event
        self.rad = rad
        self.folder = "simulation_results/{dn}/{rad}/".format(
            dn=self.event.strftime("%Y.%m.%d"), rad=self.rad
        )
        os.makedirs(self.folder, exist_ok=True)
        self.cfg = cfg
        self.hdw = pydarn.read_hdw_file(self.rad)
        self._estimate_bearing_()
        return

    def _estimate_bearing_(self):
        """Estimate laitude and logitude bearings"""
        fname = self.folder + f"bearing_{'%02d'%self.beam}.mat"
        bearing = (self.beam - self.hdw.beams) * self.hdw.beam_separation
        lat, lon = (self.hdw.geographic.lat, self.hdw.geographic.lon)
        p = (lat, lon)
        gc = GC(p, p)
        dist = np.linspace(
            0, self.cfg.max_ground_range_km, self.cfg.number_of_ground_step_km
        )

        m = {}
        lats, lons = [], []
        for d in dist:
            x = gc.destination(p, bearing, distance=d)
            lats.append(x[0])
            lons.append(x[1])
        m["dist"], m["lat"], m["lon"] = dist, np.array(lats), np.array(lons)
        (
            m["olat"],
            m["olon"],
            m["rb"],
            m["num_range"],
            m["max_range"],
            m["range_inc"],
        ) = (
            lat,
            lon,
            bearing,
            float(len(dist)),
            float(self.cfg.max_ground_range_km),
            float(dist[1] - dist[0]),
        )
        m["ht"] = np.arange(
            self.cfg.start_height_km,
            self.cfg.end_height_km,
            self.cfg.height_incriment_km,
        ).astype(float)
        m["start_height"], m["height_inc"], m["num_heights"], m["heights"] = (
            float(self.cfg.start_height_km),
            float(self.cfg.height_incriment_km),
            float(
                len(
                    np.arange(
                        self.cfg.start_height_km,
                        self.cfg.end_height_km,
                        self.cfg.height_incriment_km,
                    )
                )
            ),
            np.arange(
                self.cfg.start_height_km,
                self.cfg.end_height_km,
                self.cfg.height_incriment_km,
            ),
        )

        m["freq"], m["tol"], m["nhops"] = (
            float(self.cfg.frequency),
            float(1e-7),
            float(self.cfg.nhops),
        )
        m["elev_s"], m["elev_i"], m["elev_e"] = (
            float(self.cfg.start_elevation),
            float(self.cfg.elevation_inctiment),
            float(self.cfg.end_elevation),
        )
        m["elvs"] = np.arange(
            float(self.cfg.start_elevation),
            float(self.cfg.end_elevation),
            float(self.cfg.elevation_inctiment),
        )
        m["radius_earth"] = self.cfg.radius_earth
        savemat(fname, m)
        self.bearing_object = copy.copy(m)
        return

    def compile(self, density):
        """Compute RT using Pharlap"""
        fname = self.folder + "{date}.{bm}.elv(<elv>).csv".format(
            bm="%02d" % self.beam, date=self.event.strftime("%H%M")
        )
        pwd = os.getcwd() + "/pharlap/pharlap_4.1.3/dat"
        cmd = "export DIR_MODELS_REF_DAT={pwd};\
                cd pharlap;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fname='{fname}';bm={bm};\
                rt_1D;exit;\"".format(
            pwd=pwd,
            ut=self.event.strftime("%Y %m %d %H %M"),
            rad=self.rad,
            dic=self.folder,
            bm=self.beam,
            fname=fname,
        )
        logger.info(f"Running command: {cmd}")
        os.system(cmd)
        self.load_simulated_results(density)
        self.run_doppler_computation()
        return

    def run_doppler_computation(self):
        """
        Run/Compute simulated rays Doppler
        """
        self.dop = Doppler(
            self.bearing_object["dist"],
            self.bearing_object["ht"],
            self.pharlap_simulation["density"],
            self.pharlap_simulation["density"]*1.1,
        )
        self.pharlap_simulation["ray_dop_sim"] = []
        for indx, elv in enumerate(self.bearing_object["elvs"]):
            ray = self.dop.compute_doppler_f_from_ray(
                self.pharlap_simulation["rays"][indx], elv,
                self.bearing_object["freq"]
            )
            self.pharlap_simulation["ray_dop_sim"].append(ray)
        return

    def load_simulated_results(self, density):
        """
        Load simulated rays and e-density
        """
        fname_prev = self.folder + "{dn}_{bm}.mat".format(
                dn=(args.event - dt.timedelta(minutes=1)).strftime("%H.%M"), 
                bm="%02d" % self.beam
            )
        files = glob.glob(
            self.folder
            + "{time}.{beam}.elv(*).csv".format(
                time=self.event.strftime("%H%M"), beam="%02d" % self.beam
            )
        )
        files.sort()
        self.pharlap_simulation = {
            "time": self.event,
            "beam": self.beam,
            "rad": self.rad,
            "rays": [],
            "bearing_object": self.bearing_object,
            "density": density,
        }
        for f in files:
            self.pharlap_simulation["rays"].append(pd.read_csv(f))
        return


def execute_gitm_simulations(
    rtobj,
    cfg,
    args,
):
    """
    Execute GITM simulation for ray tracing
    """
    fname_now = rtobj.folder + "{dn}_{bm}.mat".format(
        dn=args.event.strftime("%H.%M"), bm="%02d" % args.beam
    )
    fname_prev = rtobj.folder + "{dn}_{bm}.mat".format(
        dn=(args.event - dt.timedelta(minutes=30)).strftime("%H.%M"), bm="%02d" % args.beam
    )
    gitm = GITM.create_density_files(
        cfg,
        args.event.year,
        f"dataset/GITM/{args.event.strftime('%Y%m%d')}/",
        "eden",
        args.event,
        rtobj.bearing_object["lat"],
        rtobj.bearing_object["lon"],
        rtobj.bearing_object["ht"],
        fname_now,
        fname_prev
    )
    rtobj.compile(gitm.param)
    import plots
    plots.plot_rays(
        rtobj.folder,
        rtobj.pharlap_simulation,
        "GITM + Eclipse Mask (Mrak et al., 2022)",
        overlay_density=True,
        showrefract=False,
    )
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--model", default="gitm", help="Model name [gitm]"
    )
    parser.add_argument("-r", "--rad", default="cvw", help="Radar code (default cvw)")
    parser.add_argument(
        "-bm", "--beam", default=11, type=int, help="Radar beam (default 11)"
    )
    parser.add_argument(
        "-ev",
        "--event",
        default=dt.datetime(2017, 8, 21, 17, 30),
        help="Event date for simulation [YYYY-mm-ddTHH:MM]",
        type=dparser.isoparse,
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))
    if args.model == "gitm":
        cfg = read_params()
        rtobj = RayTrace(args.event, args.rad, args.beam, cfg)
        execute_gitm_simulations(rtobj, cfg, args)
