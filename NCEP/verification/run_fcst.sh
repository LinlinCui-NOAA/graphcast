#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=nems
#SBATCH --cpus-per-task=40 
#SBATCH --time=2:00:00 
#SBATCH --job-name=graphcast
#SBATCH --output=gc_output.txt
#SBATCH --error=gc_error.txt
#SBATCH --partition=hera

## Activate Conda environment
source /scratch1/NCEPDEV/nems/Linlin.Cui/miniforge3/etc/profile.d/conda.sh
conda activate graphcast
conda list

cd /scratch1/NCEPDEV/nems/Linlin.Cui/Tests/graphcast/GCERA5 

start_time=$(date +%s)

python3 run_graphcast_with_era5.py

end_time=$(date +%s)  # Record the end time in seconds since the epoch

# Calculate and print the execution time
execution_time=$((end_time - start_time))
echo "Execution time for graphcast: $execution_time seconds"
