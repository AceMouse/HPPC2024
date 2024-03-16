#!/bin/bash

make
job_id=$(sbatch job.sh)
job_out="slurm-${job_id##* }.out"
while [ ! -e ${job_out} ]
do
  sleep 1
done

cat "$job_out"
rm "$job_out"
