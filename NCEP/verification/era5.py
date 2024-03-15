import os
from datetime import datetime, timedelta

import numpy as np
import xarray as xr
import cdsapi

def calculate_6hr_accumulation(series):
    nt, ny, nx = series.shape

    nt_6hr = int(np.floor(nt / 6))
    accumulated = np.full((nt_6hr, ny, nx), fill_value=np.nan)

    for i in range(nt_6hr):
        accumulated[i, :, :] = np.sum(series[i*6:(i+1)*6, :, :], axis=0)

    return accumulated

class GetEra5Data:
    def __init__(self, start_date, end_date, sfc_file, pl_file): 
        self.start_date = start_date
        self.end_date = end_date
        self.c = cdsapi.Client()

        self.sfc_fname = sfc_file

        self.pl_fname = pl_file

        if (not os.path.isfile(self.sfc_fname)) | (not os.path.isfile(self.pl_fname)):
            self.get_all_fields_6hr()
        else:
            print(f'Both sfc and pl files exist!')


    def gen_graphcast_ic_from_era5(self, model_startdate):
        #get two time steps for IC
        time_6hr_ago = model_startdate - timedelta(hours=6)
        end_date = model_startdate + timedelta(days=0.1) 
        dates = np.arange(
            np.datetime64(time_6hr_ago, "ns"),
            np.datetime64(end_date, "ns"),
            np.timedelta64(6, "h")
        ) 

        ds_surface = xr.open_dataset(self.sfc_fname)
        #get 2D lsm/z
        lsm = ds_surface.lsm.sel(time=np.datetime64(model_startdate, "ns"))
        orog = ds_surface.z.sel(time=np.datetime64(model_startdate, "ns"))
 
        ds_surface_slice = ds_surface.sel(time=dates)
        ds_surface_slice = ds_surface_slice.drop_vars(['lsm', 'z', 'tisr'])
  
        #get vars on pl
        ds_pl = xr.open_dataset(self.pl_fname)
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

    def get_all_fields_6hr(self):

        sfc_params = [165, 166, 167, 151, 212, 129, 172, 228]
        self.sfc_fields(sfc_params)
        
        pl_params = [130, 131, 132, 135, 133, 129]
        self.pl_fields(pl_params)
   

    def sfc_fields(self, params):

        #'param': [165, 166, 167, 151, 212, 129, 172, 228],

        sfc_startdate = self.start_date - timedelta(days=1)
        
        r = self.c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'param': params,
                'date': f"{self.start_date.strftime('%Y-%m-%d')}/{self.end_date.strftime('%Y-%m-%d')}",
                'time': [
                     '00:00','01:00','02:00','03:00','04:00','05:00',
                     '06:00','07:00','08:00','09:00','10:00','11:00',
                     '12:00','13:00','14:00','15:00','16:00','17:00',
                     '18:00','19:00','20:00','21:00','22:00','23:00'
                ],
                'format': 'netcdf'
            },
        )

        fname = f"era5_surface_{self.start_date.strftime('%Y%m%d')}-{self.end_date.strftime('%Y%m%d')}_hourly.nc"
        r.download(fname)

        # down-sample to 6-hr
        ds = xr.open_dataset(fname)
        tp_6hr = calculate_6hr_accumulation(ds.tp)

        dates = np.arange(
            ds.time[0].values,
            ds.time[-1].values,
            np.timedelta64(6, 'h')
        )

        ds_slice = ds.sel(time=dates)
        ds_slice.tp.values[:,:,:] = tp_6hr

        #change dtype to np.float32 to avoid encoding to type "short"
        for var in ds_slice.variables:
            if ds_slice[var].values.dtype == 'float32':
                ds_slcie[var] = ds_slice[var].astype(np.float32)

        ds_slice.to_netcdf(self.sfc_fname)
        ds.close()
        ds_slice.close()


    def pl_fields(self, params):
        r = self.c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'param': params,
                'pressure_level': [
                    '50', '100', '150',
                    '200', '250', '300',
                    '400', '500', '600',
                    '700', '850', '925',
                    '1000',
                ],
                'date': f"{self.start_date.strftime('%Y-%m-%d')}/{self.end_date.strftime('%Y-%m-%d')}",
                'time': [
                    '00:00',
                    '06:00',
                    '12:00',
                    '18:00'
                ],
                'format': 'netcdf'
            },
        )
        r.download(self.pl_fname)
