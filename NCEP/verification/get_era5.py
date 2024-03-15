import os
from datetime import datetime, timedelta
import pathlib
import argparse

from era5 import GetEra5Data
from model import GraphCastModel

if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument('startdate', type=datetime.fromisoformat, help="The date format 'YYYY-MM-DD' or 'YYYY-MM-DD HH:MM:SS'")
    argparser.add_argument('enddate', type=datetime.fromisoformat, help="The date format 'YYYY-MM-DD' or 'YYYY-MM-DD HH:MM:SS'")

    args = argparser.parse_args()
    startdate = args.startdate
    enddate = args.enddate
    print(startdate)
    print(enddate)

    sfc = f"./era5_surface_{startdate.strftime('%Y%m%d')}-{enddate.strftime('%Y%m%d')}_6hr.nc"
    pl = f"./era5_pl_{startdate.strftime('%Y%m%d')}-{enddate.strftime('%Y%m%d')}_6hr.nc"
    print(sfc)
    print(pl)

    era5 = GetEra5Data(startdate, enddate, sfc, pl)
