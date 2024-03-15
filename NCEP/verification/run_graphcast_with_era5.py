import os
from datetime import datetime, timedelta
import pathlib
import glob

import numpy as np
import xarray as xr

from model import GraphCastModel

def gen_graphcast_ic_from_era5( model_startdate, sfc_fname, pl_fname):
    #get two time steps for IC
    time_6hr_ago = model_startdate - timedelta(hours=6)
    end_date = model_startdate + timedelta(days=0.1) 
    dates = np.arange(
        np.datetime64(time_6hr_ago, "ns"),
        np.datetime64(end_date, "ns"),
        np.timedelta64(6, "h")
    ) 

    ds_surface = xr.open_dataset(sfc_fname)
    #get 2D lsm/z
    lsm = ds_surface.lsm.sel(time=np.datetime64(model_startdate, "ns"))
    orog = ds_surface.z.sel(time=np.datetime64(model_startdate, "ns"))

    ds_surface_slice = ds_surface.sel(time=dates)
    ds_surface_slice = ds_surface_slice.drop_vars(['lsm', 'z', 'tisr'])

    #get vars on pl
    ds_pl = xr.open_dataset(pl_fname)
    ds_pl_slice = ds_pl.sel(time=dates)
    ds_pl_slice = ds_pl_slice.rename({
        'z': 'geopotential',
    })

    ds_merged = xr.merge([ds_surface_slice, ds_pl_slice, lsm, orog])
    ds_merged = ds_merged.rename({
        'latitude': 'lat',
        'longitude': 'lon',
        'msl': 'mean_sea_level_pressure',
        'u10': '10m_u_component_of_wind',
        'v10': '10m_v_component_of_wind',
        't2m': '2m_temperature',
        'tp': 'total_precipitation_6hr',
        'u': 'u_component_of_wind',
        'v': 'v_component_of_wind',
        'w': 'vertical_velocity',
        #'z': 'geopotential',
        't': 'temperature',
        'q': 'specific_humidity', 
        'z': 'geopotential_at_surface',
        'lsm': 'land_sea_mask',
    })
    ds_merged = ds_merged.assign_coords(datetime=ds_merged.time)

    # Convert data types
    ds_merged['lat'] = ds_merged['lat'].astype('float32')
    ds_merged['lon'] = ds_merged['lon'].astype('float32')
    ds_merged['level'] = ds_merged['level'].astype('int32')
    
    var_exclude = ['time', 'datetime']
    for var in ds_merged.variables:
        if var not in var_exclude:
            ds_merged[var] = ds_merged[var].astype('float32')

    # Adjust time values relative to the first time step
    ds_merged['time'] = ds_merged['time'] - ds_merged.time[0]

    # Expand dimensions
    ds_merged = ds_merged.expand_dims(dim='batch')
    ds_merged['datetime'] = ds_merged['datetime'].expand_dims(dim='batch')

    # Squeeze dimensions
    ds_merged['geopotential_at_surface'] = ds_merged['geopotential_at_surface'].squeeze('batch')
    ds_merged['land_sea_mask'] = ds_merged['land_sea_mask'].squeeze('batch')

    ds_merged_reversed = ds_merged.reindex(lat=list(reversed(ds_merged.lat)))

    steps = ds_merged_reversed.time.shape[0]
    outfile_name = f'era5_date_{model_startdate.strftime("%Y%m%d%H")}_res-0.25_levels-13_steps-{steps}.nc'
    ds_merged_reversed.to_netcdf(outfile_name)
  
    return outfile_name

if __name__ == "__main__":

    startdate = datetime(2024, 2, 6)
    rundir = '/scratch1/NCEPDEV/nems/Linlin.Cui/Tests/graphcast/GCERA5'

    sfc = glob.glob(f"{rundir}/era5_surface_*_6hr.nc")[0]
    pl = glob.glob(f"{rundir}/era5_pl_*_6hr.nc")[0]
    stats = '/scratch1/NCEPDEV/nems/AIML/gc_weights'

    input_filename = gen_graphcast_ic_from_era5(startdate, sfc, pl)
        
    runner = GraphCastModel(stats, input_filename)
    
    runner.load_pretrained_model()
    runner.load_gdas_data()
    runner.extract_inputs_targets_forcings()
    runner.load_normalization_stats()
    runner.get_predictions()
