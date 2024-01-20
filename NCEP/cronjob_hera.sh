#!/bin/bash --login

echo "Job 1 is running"
sh /scratch1/NCEPDEV/nems/AIML/graphcast/NCEP/gc_job1_hera.sh
sleep 60  # Simulating some work
echo "Job 1 completed"

echo "Job 2 is running"
sbatch /scratch1/NCEPDEV/nems/AIML/graphcast/NCEP/gc_job2_hera.sh
job2_id=$(sbatch /scratch1/NCEPDEV/nems/AIML/graphcast/NCEP/gc_job2_hera.sh | awk '{print $4}')

# Wait for job 2 to complete
while squeue -j $job2_id &>/dev/null; do
    sleep 5  # Adjust the polling interval as needed
done
sleep 5  # Simulating some work
echo "Job 2 completed"

echo "Job 3 is running"
sh /scratch1/NCEPDEV/nems/AIML/graphcast/NCEP/gc_job3_hera.sh
sleep 5  # Simulating some work
echo "Job 3 completed"
