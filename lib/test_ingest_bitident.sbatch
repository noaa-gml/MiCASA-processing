#!/bin/bash
#SBATCH --account co2
#SBATCH --time 0:10:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 8g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --job-name test_ingest_bitident

set -eu
cd "${SLURM_SUBMIT_DIR}"
mkdir -p jobs
exec Rscript lib/test_ingest_bitident.r
