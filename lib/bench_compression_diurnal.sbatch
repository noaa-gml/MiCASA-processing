#!/bin/bash
#SBATCH --account co2
#SBATCH --time 0:30:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16g
#SBATCH --output jobs/%x.o%j
#SBATCH --partition orion
#SBATCH --job-name bench_compression_diurnal

set -eu
cd "${SLURM_SUBMIT_DIR}"
mkdir -p jobs
exec Rscript lib/bench_compression_diurnal.r
