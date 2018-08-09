#!/bin/bash

#SBATCH --job-name=scenario1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-83
#SBATCH --output=job_reports/slurm-%A_%a.out
#SBATCH --partition=mrcieu

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}

Rscript sim1.r 10000 ${i}

