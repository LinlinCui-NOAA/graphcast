import os
from datetime import datetime, timedelta
import glob
import pathlib

import numpy as np
import xarray as xr
import pygrib
import matplotlib.pyplot as plt
import boto3
from botocore import UNSIGNED
from botocore.config import Config

from era5 import get_era5_data

s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

def download_from_aws(date, s3_bucket, proudct, bucket, outdir):
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if product == 'gcgdas':
        for it in np.arange(6, 242, 6):
            file_name = f'graphcastgfs.t{date.hour:02d}z.pgrb2.0p25.f{it:03d}'
            local_file_name = f'{outdir}/graphcastgfs.t{date.hour:02d}z.pgrb2.0p25.f{it:03d}'
            key = f'graphcastgfs.{date.strftime("%Y%m%d")}/{date.hour:02d}/forecasts_13_levels/{file_name}'
            print(key)
            with open(local_file_name, 'wb') as f:
                s3.download_fileobj(bucket, key, f)

    elif product == 'gfs':
        for it in np.arange(6, 242, 6):
            file_name = f'gfs.t{date.hour:02d}z.pgrb2.0p25.f{it:03d}'
            local_file_name = f'{outdir}/gfs.t{date.hour:02d}z.pgrb2.0p25.f{it:03d}'
            key = f'gfs.{date.strftime("%Y%m%d")}/{date.hour:02d}/atmos/{file_name}'
            print(key)
            with open(local_file_name, 'wb') as f:
                s3.download_fileobj(bucket, key, f)

    elif product == 'gdas':
        for it in [0, 6]:
            file_name = f'gdas.t{date.hour:02d}z.pgrb2.0p25.f{it:03d}'
            local_file_name = f'{outdir}/gdas.t{date.hour:02d}z.pgrb2.0p25.f{it:03d}'
            key = f'gdas.{date.strftime("%Y%m%d")}/{date.hour:02d}/atmos/{file_name}'
            print(key)
            with open(local_file_name, 'wb') as f:
                s3.download_fileobj(bucket, key, f)

def get_forecast_grib_data(date, product, var, level_type, desired_level):
    print(f'product is {product}')

    if product == 'gcgdas':
        path = f'forecast/forecast-gdas/graphcastgfs.{date.strftime("%Y%m%d")}/{date.hour:02d}'
        if not pathlib.Path(path).is_dir():
            download_from_aws(s3_bucket, 'gcgdas', 'noaa-nws-graphcastgfs-pds', path)
        nc_fname = f'{path}/forecast_{product}_{var}_date_{date.strftime("%Y%m%d%H")}_levels-13.nc'

    elif product == 'gfs':
        path = f'forecast/forecast-gfs/gfs.{date.strftime("%Y%m%d")}/{date.hour:02d}'
        if not pathlib.Path(path).is_dir():
            download_from_aws(s3_bucket, 'gfs', 'noaa-gfs-bdp-pds', path)
        nc_fname = f'{path}/forecast_{product}_{var}_date_{date.strftime("%Y%m%d%H")}_levels-13.nc'

    elif product == 'gdas':
        path = f'truth/gdas.{date.strftime("%Y%m%d")}/{date.hour:02d}'
        if not pathlib.Path(path).is_dir():
            download_from_aws(s3_bucket, 'gdas', 'noaa-gfs-bdp-pds', path)
        nc_fname = f'{path}/truth_{product}_{var}_date_{date.strftime("%Y%m%d%H")}_levels-13.nc'

    else:
        raise ValueError(f"{product} is not in ['gcgdas', 'gfs', 'gdas']")

    if os.path.isfile(nc_fname):
        print(f'Reading from nc file {nc_fname}:')
        ds = xr.open_dataset(nc_fname)
        return ds[var]
    
    if product in ['gcgdas', 'gfs']:
        files = glob.glob(f'{path}/*pgrb2*')
        files.sort()
    elif product == 'gdas':
        files = []
        dates = np.arange(date + timedelta(hours=6), date + timedelta(days=10.1), timedelta(hours=6)).astype(datetime)
        for date2 in dates:
            if var != 'tp':
                fname = f'truth/gdas.{date2.strftime("%Y%m%d")}/{date2.hour:02d}/gdas.t{date2.hour:02d}z.pgrb2.0p25.f000'
            elif var == 'tp':
                date2_tp = date2 - timedelta(hours=6)
                fname = f'truth/gdas.{date2_tp.strftime("%Y%m%d")}/{date2_tp.hour:02d}/gdas.t{date2_tp.hour:02d}z.pgrb2.0p25.f006'
            files.append(fname)

    daMerged = []
    for fname in files:
        grbfile = pygrib.open(fname)

        # Find the matching grib message
        #variable_message = grbfile.select(shortName=var, typeOfLevel=level_type, level=desired_level)
        if product == 'gcgdas' and var == 'tp': 
            variable_message = grbfile.select(shortName='unknown', typeOfLevel=level_type, level=desired_level)
        else:
            variable_message = grbfile.select(shortName=var, typeOfLevel=level_type, level=desired_level)
    
        # create a netcdf dataset using the matching grib message
        lats, lons = variable_message[0].latlons()
        lats = lats[:,0]
        lons = lons[0,:]
    
        #check latitude range
        reverse_lat = False
        if lats[0] > 0:
            reverse_lat = True
            lats = lats[::-1]
    
        steps = variable_message[0].validDate
        if var=='tp':
            steps = steps + timedelta(hours=6)
        #precipitation rate has two stepType ('instant', 'avg'), use 'instant')
        if len(variable_message) > 2:
            data = []
            for message in variable_message:
                data.append(message.values)
            data = np.array(data)
            if reverse_lat:
                data = data[:, ::-1, :]
        else:
            data = variable_message[0].values
            if reverse_lat:
                data = data[::-1, :]
    
        if len(data.shape) == 2:
            da = xr.Dataset(
                data_vars={
                    var: (['lat', 'lon'], data.astype('float32'))
                },
                coords={
                    'lon': lons.astype('float32'),
                    'lat': lats.astype('float32'),
                    'time': steps,  
                }
            )
        elif len(data.shape) == 3:
            da = xr.Dataset(
                data_vars={
                    var: (['level', 'lat', 'lon'], data.astype('float32'))
                },
                coords={
                    'lon': lons.astype('float32'),
                    'lat': lats.astype('float32'),
                    'level': np.array(desired_level).astype('int32'),
                    'time': steps,  
                }
            )
    
        daMerged.append(da)

    ds = xr.concat(daMerged, dim='time')

    ds.to_netcdf(nc_fname)
    return ds[var]

def get_gdas_data(times, var, level):
    ds = xr.open_dataset('truth/source-gdas_date-2024020606_res-0.25_levels-13_steps-121.nc')
    if var == 'geopotential':
        ds[var] = ds[var] / 9.80665

    if var == 'total_precipitation_6hr':
        ds[var] = ds[var] * 1000

    ds = ds.assign_coords(time=ds.datetime.values[0,:])

    if level > 50:
        ds = ds.sel(time=times).sel(level=level)
    else:
        ds = ds.sel(time=times)
   
    return ds[var].squeeze('batch') 

def get_gcera5_data(date, times, var, level):
    ds = xr.open_dataset(f'forecast/gc-era5/forecasts_era5_{date.strftime("%Y%m%d%H")}_levels-13_steps-40.nc')
    if var == 'geopotential':
        ds[var] = ds[var] / 9.80665

    if var == 'total_precipitation_6hr':
        ds[var] = ds[var] * 1000
    #ds = ds.assign_coords(time=ds.datetime.values[0,:])
    ds = ds.assign_coords(time=[np.datetime64(date) + dt for dt in ds.time.values])
    if level > 50:
        ds = ds.sel(time=times).sel(level=level)
    else:
        ds = ds.sel(time=times)
   
    return ds[var].squeeze('batch') 

if __name__ == '__main__':

    variables = {
        '2t': {'name': '2m_temperature', 'levelType': 'heightAboveGround', 'era_levelType': 'surface', 'level': 2, 'units': 'K'},
        '10u': {'name': '10m_u_component_of_wind', 'levelType': 'heightAboveGround', 'era_levelType': 'surface', 'level': 10, 'units': 'm/s'},
        '10v': {'name': '10m_v_component_of_wind', 'levelType': 'heightAboveGround', 'era_levelType': 'surface', 'level': 10, 'units': 'm/s'},
        'tp': {'name': 'total_precipitation_6hr', 'levelType': 'surface', 'era_levelType': 'surface', 'level': 0, 'units': 'kg/m2'},
        't': {'name': 'temperature', 'levelType': 'isobaricInhPa', 'era_levelType': 'pl', 'level': 850, 'units': 'K'},
        'gh': {'name': 'geopotential', 'levelType': 'isobaricInhPa', 'era_levelType': 'pl', 'level': 500, 'units': 'gpm'},
        'u': {'name': 'u_component_of_wind', 'levelType': 'isobaricInhPa', 'era_levelType': 'pl', 'level': 200, 'units': 'm/s'},
        'v': {'name': 'v_component_of_wind', 'levelType': 'isobaricInhPa', 'era_levelType': 'pl', 'level': 200, 'units': 'm/s'},
    }
    
    startdate = datetime(2024, 2, 6)
    enddate = datetime(2024, 2, 19, 12)
    datevectors = np.arange(startdate, enddate, timedelta(hours=6)).astype(datetime)

    for short_name in variables.keys():
        #short_name = '2t' #'tp' #'2t'
        name = variables[short_name]['name'] #'10m_v_component_of_wind' #'total_precipitation_6hr'
        levelType = variables[short_name]['levelType'] #'heightAboveGround' #'isobaricInhPa' #'heightAboveGround'
        era_levelType = variables[short_name]['era_levelType'] #'surface'
        level = variables[short_name]['level'] #10
        units = variables[short_name]['units'] #'m * s-1' #'kg * m-2'

        gcgdas_gdas_mse_all = []
        gfs_gdas_mse_all = []
        gcgdas_era5_mse_all = []
        gcera5_era5_mse_all = []

        for date in datevectors:
            print(f'Plotting {short_name} for date {date.strftime("%Y%m%d-%H")}')
            #get forecast-gcgdas
            gcgdas = get_forecast_grib_data(date,'gcgdas', short_name, levelType, level)
            gcgdas_mean = gcgdas.mean(dim=["lat", "lon"])

            #get forecast-gfs
            gfs = get_forecast_grib_data(date, 'gfs', short_name, levelType, level)
            gfs_mean = gfs.mean(dim=["lat", "lon"])

            times = gcgdas.time

            #get_forecast-gcera5
            gcera5 = get_gcera5_data(date, times, name, level)
            gcera5_mean = gcera5.mean(dim=["lat", "lon"])

            #get truth-era5
            truth_era5 = get_era5_data(times, name, era_levelType, level)
            truth_era5_mean = truth_era5.mean(dim=["lat", "lon"])

            #get truth-gdas
            truth_gdas = get_forecast_grib_data(date,'gdas', short_name, levelType, level)
            truth_gdas_mean = truth_gdas.mean(dim=["lat", "lon"])
            #truth_gdas = get_gdas_data(times, name, level)
            #truth_gdas_mean = truth_gdas.mean(dim=["lat", "lon"])

            #get MSE
            gcgdas_gdas_mse = ((gcgdas - truth_gdas) ** 2).mean(dim=['lat', 'lon'])
            gcgdas_gdas_mse_all.append(gcgdas_gdas_mse.values)
            gfs_gdas_mse = ((gfs - truth_gdas) ** 2).mean(dim=['lat', 'lon'])
            gfs_gdas_mse_all.append(gfs_gdas_mse.values)
            gcgdas_era5_mse = ((gcgdas - truth_era5) ** 2).mean(dim=['lat', 'lon'])
            gcgdas_era5_mse_all.append(gcgdas_era5_mse.values)
            gcera5_era5_mse = ((gcera5 - truth_era5) ** 2).mean(dim=['lat', 'lon'])
            gcera5_era5_mse_all.append(gcera5_era5_mse.values)


            #plot
            fig, axs = plt.subplots(2, 1, figsize=(6, 8))
            axs[0].plot(times, gcgdas_mean, label='gcgdas')
            axs[0].plot(times, gfs_mean, label='gfs')
            axs[0].plot(times, gcera5_mean, label='gcera5')
            axs[0].plot(times, truth_era5_mean, label='era5')
            axs[0].plot(times, truth_gdas_mean, label='gdas')
            #axs[0].set_xticklabels(rotation=45)
            axs[0].set_xticks(axs[0].get_xticks(), [])
            axs[0].set_ylabel(f'{name}')
            axs[0].legend()
            axs[0].set_title(f'{name} ({units}) for {date.strftime("%Y%m%d-%H:00:00")}')

            axs[1].plot(times, gcgdas_gdas_mse, label='gcgdas_vs_gdas')
            axs[1].plot(times, gfs_gdas_mse, label='gfs_vs_gdas')
            axs[1].plot(times, gcgdas_era5_mse, label='gcgdas_vs_era5')
            axs[1].plot(times, gcera5_era5_mse, label='gcera5_vs_era5')
            axs[1].legend()
            axs[1].set_xticks(axs[1].get_xticks(), axs[1].get_xticklabels(), rotation=45, ha='right')
            axs[1].set_ylabel(f'MSE of {short_name}')
            plt.savefig(f'MSE_{short_name}_{level}m_{date.strftime("%Y%m%d%H")}.png')
            plt.close()

        x = np.arange(6, 246, 6)
        plt.plot(x, np.array(gcgdas_gdas_mse_all).mean(axis=0), label='gcgdas_vs_gdas')
        plt.plot(x, np.array(gfs_gdas_mse_all).mean(axis=0), label='gfs_vs_gdas')
        plt.plot(x, np.array(gcgdas_era5_mse_all).mean(axis=0), label='gcgdas_vs_era5')
        plt.plot(x, np.array(gcera5_era5_mse_all).mean(axis=0), label='gcera5_vs_era5')
        plt.legend()
        plt.title(f'{short_name} ({units}) from {startdate.strftime("%Y%m%d%H")} to {enddate.strftime("%Y%m%d%H")}')
        #plt.xticks(plt.xticks(), x, rotation=45, ha='right')
        plt.ylabel(f'MSE of {short_name}')
        plt.savefig(f'MSE_{short_name}_{level}m_mean_over_{startdate.strftime("%Y%m%d%H")}_{enddate.strftime("%Y%m%d%H")}.png')
        plt.close()
