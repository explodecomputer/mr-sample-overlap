#!/bin/bash

#SBATCH --job-name=scenario3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1212
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --partition=mrcieu

echo "Running on ${HOSTNAME}"
module add languages/r/3.6.0

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript sim_overlap.r 20000 ${i}

