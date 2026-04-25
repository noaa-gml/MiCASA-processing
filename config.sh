#!/bin/sh
# MiCASA pipeline configuration — sourced by every shell script in this repo.
#
# Override any value on the command line, e.g.:
#     MICASA_YEAR=2026 MICASA_VERSION=vNRT ./run_year.sh
#
# After sourcing, $WORK_DIR is the repo/working directory. All other paths
# (DAILY_1X1_DIR etc.) are relative names; combine with $WORK_DIR as needed.

# ---- Per-invocation knobs ---------------------------------------------------

# Year being processed.
: "${MICASA_YEAR:=2025}"

# Product version: "v1" (final, archived) or "vNRT" (near-real-time).
: "${MICASA_VERSION:=v1}"

# Year range for bulk operations (climatologies, multi-year ingest, etc.).
: "${MICASA_YEAR_START:=2001}"
: "${MICASA_YEAR_END:=${MICASA_YEAR}}"

# ---- Site config ------------------------------------------------------------

# SLURM mail contact.
: "${MAIL_USER:=ashley.pera@noaa.gov}"

# Top of the GFED-CASA tree on /work2. Each year lives at $BASE_DIR/$YEAR/MiCASA_v1.
: "${BASE_DIR:=/work2/noaa/co2/GFED-CASA}"

# Working directory — defaults to this year's tree.
: "${WORK_DIR:=${BASE_DIR}/${MICASA_YEAR}/MiCASA_v1}"

# Upstream NCCS portal (no trailing slash).
: "${PORTAL_URL_BASE:=https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA}"

# ---- Layout (paths relative to WORK_DIR) ------------------------------------

DAILY_1X1_DIR="daily_1x1"
MONTHLY_1X1_DIR="monthly_1x1"
ERA5_DIR="ERA5"
RAW_SRC_DIR="portal.nccs.nasa.gov"
JOBS_DIR="jobs"

export MICASA_YEAR MICASA_VERSION MICASA_YEAR_START MICASA_YEAR_END
export MAIL_USER BASE_DIR WORK_DIR PORTAL_URL_BASE
export DAILY_1X1_DIR MONTHLY_1X1_DIR ERA5_DIR RAW_SRC_DIR JOBS_DIR
