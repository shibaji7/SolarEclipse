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

import copy
import datetime as dt
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd
import plots
from loguru import logger

sys.path.extend(["py/"])

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


class Fetch2DModel(object):
    """
    Fetching dataset from NETCDF files.
    """

    def __init__(
        self, fname, params, event=dt.datetime(2017, 8, 21), lat_res=1, lon_res=1, h=150
    ):
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
        self.h = h
        self.lat_res = lat_res
        self.lon_res = lon_res
        self.dataframe = pd.DataFrame()
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
        # Change these values
        x = np.arange(0, 91, self.lat_res)
        y = np.arange(-180, 0, self.lon_res)
        X, Y = np.meshgrid(x, y)
        self.dataframe["lat"], self.dataframe["lon"] = X.ravel(), Y.ravel()
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
            "mil": {"lat": 42.61, "lon": -71.49},
        }
        return locations[stn]

    def __get_latlon_index__(self, lat, lon):
        """
        Get latitude, longitude index
        """
        _ij_ = (np.argmin(np.abs(self.latx - lat)), np.argmin(np.abs(self.lonx - lon)))
        return _ij_

    def __get_h_index__(self, Z, h=None):
        """
        Get h index
        """
        h = h if h else self.h
        _k_ = np.argmin(np.abs(h - Z))
        return _k_

    def run_2D_interpolate(self, row):
        """
        Run fitting algorithm to fit the data
        into equal height bin for a lat-lon.
        Parameters:
        -----------
        stn: Station name
        loc: Lat/Lon position
        """
        tdx = self.time.index(self.event)
        tdx0 = tdx - 1
        idx, jdx = self.__get_latlon_index__(row["lat"], row["lon"])
        zdx = self.__get_h_index__(self.Z3[tdx, :, idx, jdx])
        zdx0 = self.__get_h_index__(self.Z3[tdx0, :, idx, jdx])
        for i, p in enumerate(self.params):
            row[p["name"] + ".dt"] = (
                self.dataset[i]["value"][tdx, zdx, idx, jdx]
                - self.dataset[i]["value"][tdx0, zdx0, idx, jdx]
            )
            row[p["name"]] = self.dataset[i]["value"][tdx, zdx, idx, jdx]
        return row

    def run_2D_interpolation(self, row):
        """
        Run fitting algorithm to fit the data
        into equal height bin for a lat-lon.
        Parameters:
        -----------
        stn: Station name
        loc: Lat/Lon position
        """
        tdx = self.time.index(self.event)
        idx, jdx = self.__get_latlon_index__(row["lat"], row["lon"])
        zdx = self.__get_h_index__(self.Z3[tdx, :, idx, jdx], row["h"])
        for i, p in enumerate(self.params):
            row[p["name"]] = self.dataset[i]["value"][tdx, zdx, idx, jdx]
        return row

    def fetch_data(self, time):
        logger.info(f"Run ineterpolation {self.fname}")
        self.tm = time
        self.dataframe = self.dataframe.apply(self.run_2D_interpolate, axis=1)
        return


class Diff2DWACCMX(object):
    """
    Class is dedicated to
    differential parameters.
    """

    def __init__(
        self,
        eclipse_file,
        bgc_file,
        params,
        time,
        event=dt.datetime(2017, 8, 21),
        lat_res=1,
        lon_res=1,
        h=150,
    ):
        """
        Parameters:
        -----------
        eclipse_file: .nc file name containg eclipse run
        bgc_file: .nc file name containg background run
        params: List of parameter details for fetching
        event: Event date
        h_interpolate: Kind of height interpolation [nearest, fitted]
        Hs: Heights in km, only required for 'nearest'
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
        self.h = h
        self.lat_res = lat_res
        self.lon_res = lon_res
        self.time = time
        self.eclipse = Fetch2DModel(
            eclipse_file, copy.copy(params), event, lat_res, lon_res, h
        )
        self.bgc = Fetch2DModel(bgc_file, copy.copy(params), event, lat_res, lon_res, h)
        self.run_extractions()
        return

    def run_extractions(self):
        """
        Run data extraction
        """
        self.eclipse.fetch_data(self.time)
        self.bgc.fetch_data(self.time)
        return

    def diffential_difference_data(self):
        """
        Get data in TS/Altitude profiles
        Parameters:
        -----------
        kind: TS/A
        """
        self.dataframe = self.eclipse.dataframe.copy()
        for i, p in enumerate(self.params):
            self.dataframe[p["name"] + ".ddif"] = (
                self.eclipse.dataframe[p["name"] + ".dt"]
                - self.bgc.dataframe[p["name"] + ".dt"]
            )
        self.dataframe["Op_CHM.ddif"] = (
            self.dataframe["Op_CHMP.ddif"] - self.dataframe["Op_CHML.ddif"]
        )
        return self.dataframe


def summary_plots_2D(
    params,
    time,
    h=150,
    lat_res=1,
    lon_res=1,
    p_names=[],
    colorbar_labels=[],
    fname=None,
):
    dw = Diff2DWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        time=time,
        h=h,
        lat_res=lat_res,
        lon_res=lon_res,
    )
    dw.diffential_difference_data()
    df = dw.diffential_difference_data()
    o = plots.OverlayDataset()
    for p, cl in zip(p_names, colorbar_labels):
        ax = o.plot_data(df, time, h, p, cl)
    del df
    del dw
    if fname:
        o.save(fname)
    o.close()
    return


if __name__ == "__main__":
    params = [
        {
            "name": "e",
            "unit_multiplier": 1.0,
            "unit": "cc",
            "mol_to_cc": "T",
            "interpolate": {"scale": "log", "type": "linear"},
        }
    ]

    dw = Diff2DWACCMX(
        "dataset/29jan2022/ECL_1100_0100_all.nc",
        "dataset/29jan2022/BGC_1100_0100_all.nc",
        params=params,
        time=dt.datetime(2017, 8, 21, 16),
    )
    dw.diffential_difference_data()
