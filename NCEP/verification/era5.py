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
