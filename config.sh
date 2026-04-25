#!/bin/sh
# MiCASA pipeline configuration — sourced by every shell script in this repo.
#
# Override any value on the command line, e.g.:
#     MICASA_YEAR=2026 MICASA_VERSION=vNRT ./run_year.sh
#
# WORK_DIR is the working/checkout directory (the location of these scripts
# and where intermediate outputs land). MICASA_YEAR is the *data year* being
# processed — these are independent: a single checkout can process any year.

# ---- Per-invocation knobs ---------------------------------------------------

# Year being processed.
: "${MICASA_YEAR:=2025}"

# Product version: "v1" (final, archived) or "vNRT" (near-real-time).
: "${MICASA_VERSION:=v1}"

# Year range for bulk operations (climatologies, multi-year ingest, etc.).
: "${MICASA_YEAR_START:=2001}"
: "${MICASA_YEAR_END:=${MICASA_YEAR}}"

# Month range for partial-year operations (e.g. mid-year NRT updates).
# Defaults to all of $MICASA_YEAR; override for FastTrack-style runs.
: "${MICASA_MONTH_START:=1}"
: "${MICASA_MONTH_END:=12}"

# ---- Site config ------------------------------------------------------------

# SLURM mail contact.
: "${MAIL_USER:=ashley.pera@noaa.gov}"

# Top of the GFED-CASA tree on /work2. Each year may live at $BASE_DIR/$YEAR/MiCASA_v1
# (used by link_old_micasa_*.sh and check_unchanged.sh).
: "${BASE_DIR:=/work2/noaa/co2/GFED-CASA}"

# Working directory — defaults to the directory of the invoking script.
# A single checkout can process any year by exporting MICASA_YEAR; WORK_DIR
# need not change. Set explicitly to point at a different checkout.
if [ -z "${WORK_DIR:-}" ]; then
    if [ -n "${BASH_SOURCE:-}" ]; then
        WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    else
        WORK_DIR="$(pwd)"
    fi
fi

# Upstream NCCS portal (no trailing slash).
: "${PORTAL_URL_BASE:=https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA}"

# ---- Layout (paths relative to WORK_DIR) ------------------------------------

DAILY_1X1_DIR="daily_1x1"
MONTHLY_1X1_DIR="monthly_1x1"
ERA5_DIR="ERA5"
RAW_SRC_DIR="portal.nccs.nasa.gov"
JOBS_DIR="jobs"

export MICASA_YEAR MICASA_VERSION MICASA_YEAR_START MICASA_YEAR_END
export MICASA_MONTH_START MICASA_MONTH_END
export MAIL_USER BASE_DIR WORK_DIR PORTAL_URL_BASE
export DAILY_1X1_DIR MONTHLY_1X1_DIR ERA5_DIR RAW_SRC_DIR JOBS_DIR
