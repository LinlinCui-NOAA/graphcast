# GraphCast verification

In order to do GraphCastGFS verification, we'll need the following data:

- ERA5
- GraphCast-ERA5
- GDAS
- GraphCast-GFS
- GFS

### Get ERA5 data with

    python get_era5.py {startdate, format: yyyy-mm-dd} {enddate,format: yyyy-mm-dd)

### Run graphcast with ERA5 
Change the startdate in the code and submit the job:

    sbatch run_fcst.sh

To run multiple cycles, after submitting the job, change the `last_cycle` and other input variables 
in `auto_run_gc.py` and run it in the background (detached terminal):

    python auto_run_gc.py    

It will check if the model is running. If not, get the last forecast cycle from output folder and re-submit
the job for the next cycle. 
