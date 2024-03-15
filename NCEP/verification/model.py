'''
Description: Script to call the graphcast model using gdas products
Author: Sadegh Sadeghi Tabas (sadegh.tabas@noaa.gov)
Revision history:
    -20231218: Sadegh Tabas, initial code
    -20240118: Sadegh Tabas, S3 bucket module to upload data, adding forecast length, Updating batch dataset to account for forecast length
    -20240125: Linlin Cui, added a capability to save output as grib2 format
    -20240205: Sadegh Tabas, made the code clearer, added 37 pressure level option, updated upload to s3
'''
import os
import argparse
from datetime import datetime, timedelta
import pathlib
import glob
import dataclasses
import functools
import math
import re

import jax
import numpy as np
import xarray
#import boto3
import pandas as pd
import haiku as hk

from graphcast import autoregressive
from graphcast import casting
from graphcast import checkpoint
from graphcast import data_utils
from graphcast import graphcast
from graphcast import normalization
from graphcast import rollout

class GraphCastModel:
    def __init__(self, pretrained_model_path, gdas_data_path, output_dir=None, num_pressure_levels=13, forecast_length=40):
        self.pretrained_model_path = pretrained_model_path
        self.gdas_data_path = gdas_data_path
        self.forecast_length = forecast_length
        self.num_pressure_levels = num_pressure_levels
        
        if output_dir is None:
            self.output_dir = os.path.join(os.getcwd(), f"forecasts_{str(self.num_pressure_levels)}_levels")  # Use current directory if not specified
        else:
            self.output_dir = os.path.join(output_dir, f"forecasts_{str(self.num_pressure_levels)}_levels")
        os.makedirs(self.output_dir, exist_ok=True)
        
        self.params = None
        self.state = {}
        self.model_config = None
        self.task_config = None
        self.diffs_stddev_by_level = None
        self.mean_by_level = None
        self.stddev_by_level = None
        self.current_batch = None
        self.inputs = None
        self.targets = None
        self.forcings = None
        self.dates = None
        

    def load_pretrained_model(self):
        """Load pre-trained GraphCast model."""
        if self.num_pressure_levels==13:
            model_weights_path = f"{self.pretrained_model_path}/params/GraphCast_operational - ERA5-HRES 1979-2021 - resolution 0.25 - pressure levels 13 - mesh 2to6 - precipitation output only.npz"
        else:
            model_weights_path = f"{self.pretrained_model_path}/params/GraphCast - ERA5 1979-2017 - resolution 0.25 - pressure levels 37 - mesh 2to6 - precipitation input and output.npz"

        with open(model_weights_path, "rb") as f:
            ckpt = checkpoint.load(f, graphcast.CheckPoint)
            self.params = ckpt.params
            self.state = {}
            self.model_config = ckpt.model_config
            self.task_config = ckpt.task_config

    def load_gdas_data(self):
        """Load GDAS data."""
        #with open(gdas_data_path, "rb") as f:
        #    self.current_batch = xarray.load_dataset(f).compute()
        self.current_batch = xarray.load_dataset(self.gdas_data_path).compute()
        self.dates =  pd.to_datetime(self.current_batch.datetime.values)
        
        if (self.forecast_length + 2) > len(self.current_batch['time']):
            print('Updating batch dataset to account for forecast length')
            
            diff = int(self.forecast_length + 2 - len(self.current_batch['time']))
            ds = self.current_batch

            # time and datetime update
            curr_time_range = ds['time'].values.astype('timedelta64[ns]')
            new_time_range = (np.arange(len(curr_time_range) + diff) * np.timedelta64(6, 'h')).astype('timedelta64[ns]')
            ds = ds.reindex(time = new_time_range)
            curr_datetime_range = ds['datetime'][0].values.astype('datetime64[ns]')
            new_datetime_range = curr_datetime_range[0] + np.arange(len(curr_time_range) + diff) * np.timedelta64(6, 'h')
            ds['datetime'][0]= new_datetime_range

            self.current_batch = ds
            print('batch dataset updated')
            
        
    def extract_inputs_targets_forcings(self):
        """Extract inputs, targets, and forcings from the loaded data."""
        self.inputs, self.targets, self.forcings = data_utils.extract_inputs_targets_forcings(
            self.current_batch, target_lead_times=slice("6h", f"{self.forecast_length*6}h"), **dataclasses.asdict(self.task_config)
        )

    def load_normalization_stats(self):
        """Load normalization stats."""
        
        diffs_stddev_path = f"{self.pretrained_model_path}/stats/diffs_stddev_by_level.nc"
        mean_path = f"{self.pretrained_model_path}/stats/mean_by_level.nc"
        stddev_path = f"{self.pretrained_model_path}/stats/stddev_by_level.nc"
        
        with open(diffs_stddev_path, "rb") as f:
            self.diffs_stddev_by_level = xarray.load_dataset(f).compute()
        with open(mean_path, "rb") as f:
            self.mean_by_level = xarray.load_dataset(f).compute()
        with open(stddev_path, "rb") as f:
            self.stddev_by_level = xarray.load_dataset(f).compute()
    
    # Jax doesn't seem to like passing configs as args through the jit. Passing it
    # in via partial (instead of capture by closure) forces jax to invalidate the
    # jit cache if you change configs.
    def _with_configs(self, fn):
        return functools.partial(fn, model_config=self.model_config, task_config=self.task_config,)

    # Always pass params and state, so the usage below are simpler
    def _with_params(self, fn):
        return functools.partial(fn, params=self.params, state=self.state)

    # Deepmind models aren't stateful, so the state is always empty, so just return the
    # predictions. This is requiredy by the rollout code, and generally simpler.
    @staticmethod
    def _drop_state(fn):
        return lambda **kw: fn(**kw)[0]

    def load_model(self):
        def construct_wrapped_graphcast(model_config, task_config):
            """Constructs and wraps the GraphCast Predictor."""
            # Deeper one-step predictor.
            predictor = graphcast.GraphCast(model_config, task_config)

            # Modify inputs/outputs to `graphcast.GraphCast` to handle conversion to
            # from/to float32 to/from BFloat16.
            predictor = casting.Bfloat16Cast(predictor)

            # Modify inputs/outputs to `casting.Bfloat16Cast` so the casting to/from
            # BFloat16 happens after applying normalization to the inputs/targets.
            predictor = normalization.InputsAndResiduals(predictor, diffs_stddev_by_level=self.diffs_stddev_by_level, mean_by_level=self.mean_by_level, stddev_by_level=self.stddev_by_level,)

            # Wraps everything so the one-step model can produce trajectories.
            predictor = autoregressive.Predictor(predictor, gradient_checkpointing=True,)
            return predictor

        @hk.transform_with_state
        def run_forward(model_config, task_config, inputs, targets_template, forcings,):
            predictor = construct_wrapped_graphcast(model_config, task_config)
            return predictor(inputs, targets_template=targets_template, forcings=forcings,)
        
        jax.jit(self._with_configs(run_forward.init))
        self.model = self._drop_state(self._with_params(jax.jit(self._with_configs(run_forward.apply))))
    
 
    def get_predictions(self):
        """Run GraphCast and save forecasts to a NetCDF file."""

        print (f"start running GraphCast for {self.forecast_length} steps --> {self.forecast_length*6} hours.")
        self.load_model()
           
        # output = self.model(self.model ,rng=jax.random.PRNGKey(0), inputs=self.inputs, targets_template=self.targets * np.nan, forcings=self.forcings,)
        forecasts = rollout.chunked_prediction(self.model, rng=jax.random.PRNGKey(0), inputs=self.inputs, targets_template=self.targets * np.nan, forcings=self.forcings,)
        filename = f"forecasts_era5_{self.gdas_data_path.split('_')[2]}_levels-{self.num_pressure_levels}_steps-{self.forecast_length}.nc"
        output_netcdf = os.path.join(self.output_dir, filename)
        
        # save forecasts
        forecasts.to_netcdf(output_netcdf)
        # print (f"GraphCast run completed successfully, you can find the GraphCast forecasts in the following directory:\n {output_netcdf}")

        #self.save_grib2(forecasts)
        forecasts.close()
