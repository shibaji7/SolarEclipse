#!/usr/bin/env python

"""fetch.py: fetch simulated (WACCM-X) data from netcdf files."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys

sys.path.extend(["py/"])

import copy
import datetime as dt

import netCDF4 as nc
import numpy as np
import pandas as pd
import scipy.io as io
import utils
from loguru import logger

pconst = {
    "boltz": 1.38066e-23,  # Boltzmann constant  in Jule K^-1
    "h": 6.626e-34,  # Planks constant  in ergs s
    "c": 2.9979e08,  # in m s^-1
    "avo": 6.023e23,  # avogadro's number
    "Re": 6371.0e3,
    "amu": 1.6605e-27,
    "q_e": 1.602e-19,  # Electron charge in C
    "m_e": 9.109e-31,  # Electron mass in kg
    "g": 9.81,  # Gravitational acceleration on the surface of the Earth
    "eps0": 1e-9 / (36 * np.pi),
    "R": 8.31,  # J mol^-1 K^-1
}


class FetchDigisonde(object):
    """
    Fetch Digisonde data from file.
    """

    def __init__(self, fname, event_start=dt.datetime(2017, 8, 21, 16)):
        """
        Parameters:
        -----------
        fname: File name of the file
        event_start: Start time
        """
        self.fname = fname
        self.event_start = event_start
        self.data = self.__read_dataset__()
        return

    def __read_dataset__(self):
        o = io.readsav(self.fname)
        o["time"] = [
            self.event_start + dt.timedelta(minutes=int(x)) for x in o["timevals"]
        ]
        ox = pd.DataFrame.from_dict(o)
        keys = [
            "foe",
            "fof1",
            "fxf1",
            "fof2",
            "fxf2",
            "ofof1",
            "ofxf1",
            "ofof2",
            "ofxf2",
        ]
        for k in keys:
            ox[k] = ox[k].apply(lambda x: np.nan if x == 0.0 else x)
        return ox


class FetchISR(object):
    """
    Fetch ISR data from file.
    """

    def __init__(
        self,
        fname,
        headers=[
            "NE",
            "TI",
            "TR",
            "VO",
            "PH+",
            "PM",
            "CO",
            "GDALT",
            "RANGE",
            "AZ1",
            "AZ2",
            "EL1",
            "EL2",
            "GDLAT",
            "GLON",
        ],
    ):
        """
        Parameters:
        -----------
        fname: File name of the file
        headers: File headers
        """
        self.fname = fname
        self.headers = headers
        self.data = self.__read_dataset__()
        return

    def __read_dataset__(self):
        """
        Read dataset from file
        """
        with open(self.fname, "r") as f:
            lines = f.readlines()
        isr = []
        hdrs = list(filter(None, lines[0].replace("\n", "").split(" ")))
        hdrs_idx = [hdrs.index(h) for h in self.headers]
        for j, l in enumerate(lines[1:]):
            l = list(filter(None, l.replace("\n", "").split(" ")))
            date = dt.datetime(
                int(l[0]), int(l[1]), int(l[2]), int(l[3]), int(l[4]), int(l[5])
            )
            o = {"DATE": date}
            for idx, k in zip(hdrs_idx, rtained_hdrs):
                o[k] = float(l[idx])
            isr.append(o)
        isr = pd.DataFrame.from_records(isr)
        isr.head()
        return isr

    def ___get_nearest_times__(self, date):
        """
        Fetch nearest time index and date
        """
        dates = self.data["DATE"].astype("datetime64[s]").unique()
        dx = dates[np.argmin([abs(x - du) for x in dates])]
        return pd.to_datetime(dx).to_pydatetime()

    def fetch_ISR_hprofile(self, date, ylim=[100, 400]):
        """
        Fetch height profile data for a date
        """
        dx = self.___get_nearest_times__(date)
        o = self.data[(self.data.DATE == dx) & (self.data.GDALT <= ylim[1])]
        o = o.dropna()
        return o

    def fetch_ISR_TS(self, hs=[150, 240]):
        """
        Fetch TS of the ISR data for fixed height
        """
        TS = {}
        for h in hs:
            o = self.data[(self.data.GDALT >= h - 10) & (self.data.GDALT <= h + 10)]
            TS[h] = o
        return TS


class FetchModel(object):
    """
    Fetching dataset from NETCDF files.
    """

    def __init__(self, fname, params, event=dt.datetime(2017, 8, 21)):
        """
        Parameters:
        -----------
        fname: .nc file name
        params: List of parameter details for fetching
        event: Event date
        -----------
        Eaxmple:
        params = [
            {
                "name": "e",
                "unit_multiplier": 1.,
                "unit": "cc",
                "mol_to_cc": "T",
                "interpolate": {
                    "scale": "log",
                    "type": "linear"
                }
            }
        ]
        """
        logger.info(f"Run extraction for {fname}")
        self.fname = fname
        self.params = params
        self.event = event
        self.load_file()
        self.extract_dimension()
        self.run_extractions()
        return

    def transform_density(self, T, var, unit="cc"):
        """
        Transform electron density from mol/mol to /cc or /cm
        Formula
        -------
        P = nKT, K = Boltzmann constant in SI
        lev: pressure in hPa (1hPa=100Pa)
        T: neutral temperature
        e/O2/...: density in mol/mol

        T in K, P in Pa, n is electron density in /cubic meter or /cc
        """
        logger.info(f"Transform density scale/unit {self.fname}")
        P = self.lev * 1e2
        u = 1 if unit == "cm" else 1e-6
        den = np.zeros_like(var)
        for i, p in enumerate(P):
            na = p * pconst["avo"] / (pconst["R"] * T[:, i, :, :])
            n = u * na * var[:, i, :, :]
            den[:, i, :, :] = n
        return den

    def load_file(self):
        """
        Load files from local dataset.
        """
        logger.info(f"Load {self.fname}")
        self.ds = nc.Dataset(self.fname)
        return

    def extract_dimension(self):
        """
        Extract data dimensions
        """
        logger.info(f"Run extraction dimension: {self.fname}")
        self.time = [
            self.event + dt.timedelta(seconds=float(d))
            for d in self.ds.variables["time"][:]
        ]
        self.latx = self.ds.variables["lat"][:]
        self.lonx = self.ds.variables["lon"][:]
        self.lonx = np.where(self.lonx > 180, self.lonx - 360, self.lonx)
        self.lev = self.ds.variables["lev"][:]
        self.Z3 = self.ds.variables["Z3GM"][:] / 1e3
        return

    def run_extractions(self):
        """
        Extract informations from
        the file.
        """
        self.dataset = {}
        for i, p in enumerate(self.params):
            logger.info(f"Load {self.fname}, {p['name']}")
            var = self.ds.variables[p["name"]][:]
            if "mol_to_cc" in p.keys():
                var = (
                    self.transform_density(
                        self.ds.variables[p["mol_to_cc"]][:], var, p["unit"]
                    )
                    * p["unit_multiplier"]
                )
            self.dataset[i] = {}
            self.dataset[i]["value"] = var
        return

    def __fetch_latlon_by_station__(self, stn):
        """
        Local RAM store latlons
        """
        locations = {
            "lusk": {"lat": 42.7625, "lon": -104.4522},
            "boulder": {"lat": 40.015, "lon": -105.2705},
            "mit": {
                "lat": 42.61,
                "lon": -71.49,
            },
        }
        return locations[stn]

    def __get_latlon_index__(self, lat, lon):
        """
        Get latitude, longitude index
        """
        _ij_ = (np.argmin(np.abs(self.latx - lat)), np.argmin(np.abs(self.lonx - lon)))
        return _ij_

    def __intp_heights__(self, h, param, scale="log", kind="linear"):
        """
        Height interpolation
        """
        if scale == "linear":
            pnew = utils.extrap1d(h, param, kind=kind)(self.intp_height)
        if scale == "log":
            pnew = 10 ** utils.extrap1d(h, np.log10(param), kind=kind)(self.intp_height)
        return pnew

    def run_height_interpolate(self, stn="lusk", loc=None, t_start=None, t_end=None):
        """
        Run fitting algorithm to fit the data
        into equal height bin for a lat-lon.
        Parameters:
        -----------
        stn: Station name
        loc: Lat/Lon position
        """
        logger.info(f"Run ineterpolation {self.fname} {stn}")
        loc = self.__fetch_latlon_by_station__(stn) if loc is None else loc
        lat, lon = loc["lat"], loc["lon"]
        it_start = self.time.index(t_start) if t_start is not None else 0
        it_end = self.time.index(t_end) + 1 if t_end is not None else -1
        idx, jdx = self.__get_latlon_index__(lat, lon)
        self.intp_height = np.arange(50, 600, 1)
        Z3 = self.Z3[it_start:it_end, :, idx, jdx]
        for i, p in enumerate(self.params):
            scale, typ = (
                p["interpolate"]["scale"],
                p["interpolate"]["type"],
            )
            logger.info(f"Run ineterpolation {self.fname} {stn} {p['name']}")
            val = self.dataset[i]["value"][it_start:it_end, :, idx, jdx]
            px = np.zeros((val.shape[0], len(self.intp_height)))
            for tj in range(val.shape[0]):
                px[tj, :] = self.__intp_heights__(Z3[tj, :], val[tj, :], scale, typ)
            self.dataset[i]["interpol_value"] = {
                "stn": stn,
                "value": px,
                "time": self.time[it_start:it_end],
            }
        return

    def fetch_data(self, kind="TS", hs=[150, 240]):
        """
        Get data in TS/Altitude profiles
        Parameters:
        -----------
        kind: TS/A
        """
        dct = {}
        if kind == "TS":
            for i, p in enumerate(self.params):
                name = p["name"]
                param = self.dataset[i]["interpol_value"]
                time = param["time"]
                for h in hs:
                    o = pd.DataFrame()
                    iH = self.intp_height.tolist().index(h)
                    o["time"], o[name] = (
                        time,
                        param["value"][:, iH],
                    )
                    dct[name + "_" + str(h)] = o
        return dct


class DiffWACCMX(object):
    """
    Class is dedicated to
    differential parameters.
    """

    def __init__(
        self,
        eclipse_file,
        bgc_file,
        params,
        event=dt.datetime(2017, 8, 21),
        stn="lusk",
        loc=None,
        t_start=None,
        t_end=None,
    ):
        """
        Parameters:
        -----------
        eclipse_file: .nc file name containg eclipse run
        bgc_file: .nc file name containg background run
        params: List of parameter details for fetching
        event: Event date
        stn: Station locations
        loc: Lat-lon dictionary
        t_strat: Start time
        t_end: End time
        -----------
        Eaxmple:
        params = [
            {
                "name": "e",
                "unit_multiplier": 1.,
                "unit": "cc",
                "mol_to_cc": "T",
                "interpolate": {
                    "scale": "log",
                    "type": "linear"
                }
            }
        ]
        """
        self.params = params
        self.stn = stn
        self.loc = loc
        self.t_start = t_start
        self.t_end = t_end
        self.eclipse = FetchModel(eclipse_file, copy.copy(params), event)
        self.bgc = FetchModel(bgc_file, copy.copy(params), event)
        self.run_extractions()
        return

    def run_extractions(self):
        """
        Run data extraction
        """
        self.eclipse.run_height_interpolate(
            self.stn, self.loc, self.t_start, self.t_end
        )
        self.bgc.run_height_interpolate(self.stn, self.loc, self.t_start, self.t_end)
        return

    def fetch_data(self, kind="TS", hs=[150, 240]):
        """
        Get data in TS/Altitude profiles
        Parameters:
        -----------
        kind: TS/A
        """
        dct = {}
        if kind == "TS":
            for i, p in enumerate(self.params):
                name = p["name"]
                ecl = self.eclipse.dataset[i]["interpol_value"]
                bgc = self.bgc.dataset[i]["interpol_value"]
                time = ecl["time"]
                for h in hs:
                    o = pd.DataFrame()
                    iH = self.eclipse.intp_height.tolist().index(h)
                    o["time"], o["ecl"], o["bgc"] = (
                        time,
                        ecl["value"][:, iH],
                        bgc["value"][:, iH],
                    )
                    o["dif"] = o["ecl"] - o["bgc"]
                    dct[name + "_" + str(h)] = o
        return dct


if __name__ == "__main__":
    pnames = [
        {
            "name": "e",
            "unit_multiplier": 1.0,
            "unit": "cc",
            "mol_to_cc": "T",
            "interpolate": {"scale": "log", "type": "linear"},
        }
    ]
    dw = DiffWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=pnames,
        t_start=dt.datetime(2017, 8, 21, 16),
        t_end=dt.datetime(2017, 8, 21, 17),
    )
    dct = dw.fetch_data()
    for k in dct.keys():
        logger.info(f"\n {dct[k]}")
