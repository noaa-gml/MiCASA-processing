#!/bin/bash
#SBATCH --account co2
#SBATCH --time 0:15:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --job-name profile_ingest_day

set -eu
cd "${SLURM_SUBMIT_DIR}"
mkdir -p jobs
exec Rscript lib/profile_ingest_day.r
