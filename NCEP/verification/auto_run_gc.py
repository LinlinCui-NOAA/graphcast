import os
import time
import subprocess
import glob

rundir = '/scratch1/NCEPDEV/nems/Linlin.Cui/Tests/graphcast/GCERA5'
queue_query_str = "squeue -u Linlin.Cui"
run_job_name = 'graphcas'
last_cycle = '2024022918'
last_file = f'{rundir}/forecasts_13_levels/forecasts_era5_{last_cycle}_levels-13_steps-40.nc'
my_print_suffix = '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'

def get_startdate(outdir):
    if '~' in outdir:
        outdir = outdir.replace('~', os.path.expanduser('~'))
    outdir = os.path.abspath(outdir)

    forecasts_files = glob.glob(f"{outdir}/forecasts_*.nc")
    forecasts_files.sort()
    if len(forecasts_files) == 0:
        raise Exception(f"No forecast results before run stopped")

    print(forecasts_files[-1])
    date = forecasts_files[-1].split('_')[4] 
    print(date)
    return date


while (not os.path.exists(last_file)):
    if run_job_name in subprocess.getoutput(queue_query_str):
        print(f'{my_print_suffix}{run_job_name} running/queueing, wait ... {my_print_suffix}', flush=True)
        time.sleep(120)
    else:
        last_date = get_startdate(f'{rundir}/forecasts_13_levels')
        year, month, day, hour = int(last_date[0:4]), int(last_date[4:6]), int(last_date[6:8]), int(last_date[8:])
        print(f'Last date is {last_date} before time out!')
        os.system(f"sed -i '/startdate = datetime/c\\    startdate = datetime({year}, {month}, {day}, {hour}) + timedelta(hours=6)\' run_graphcast_with_era5.py")
        #os.system(f'sbatch {rundir}/run_fcst.sh')
