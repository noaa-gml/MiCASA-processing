#!/usr/bin/env python
"""Verify SHA-256 checksums for downloaded MiCASA daily + monthly files.

Year range comes from $MICASA_YEAR_START / $MICASA_YEAR_END (set by config.sh
or run_year.sh). Falls back to $MICASA_YEAR (single year) if those aren't set.
Run standalone (no config sourced) defaults to 2001..current MICASA_YEAR.
"""

import os
import sys
from glob import glob
from os.path import join
from os import popen


def year_range_from_env():
    y_start = os.environ.get('MICASA_YEAR_START')
    y_end   = os.environ.get('MICASA_YEAR_END')
    y_one   = os.environ.get('MICASA_YEAR')

    if y_start and y_end:
        return list(range(int(y_start), int(y_end) + 1))
    if y_one:
        return [int(y_one)]
    # Standalone fallback — verify everything we might have on disk.
    return list(range(2001, 2026))


years = year_range_from_env()
year_glob = '{' + ','.join(f'{y:04d}' for y in years) + '}'  # brace expansion via glob

print(f'CHECKING HASHES (years: {years[0]}..{years[-1]})')

netcdf_dir = './portal.nccs.nasa.gov'

dailies_month_dirs = []
montlies_year_dirs = []
for y in years:
    dailies_month_dirs += sorted(glob(join(netcdf_dir, f'daily/{y:04d}/??/')))
    montlies_year_dirs += sorted(glob(join(netcdf_dir, f'monthly/{y:04d}/')))

dailies_month_dirs.sort()
montlies_year_dirs.sort()

incorrect_daily_hashes = 0
correct_daily_hashes = 0
for month_dir in dailies_month_dirs:
    daily_filenames_full = sorted(glob(join(month_dir, 'MiCASA_v1_flux_x3600_y1800_daily_????????.nc4')))
    daily_filenames_NRT = sorted(glob(join(month_dir, 'MiCASA_vNRT_flux_x3600_y1800_daily_????????.nc4')))
    daily_filenames = daily_filenames_full + daily_filenames_NRT

    checksum_filename_full = glob(join(month_dir, 'MiCASA_v1_flux_x3600_y1800_daily_??????_sha256.txt')) #should be in same order
    checksum_filename_NRT = glob(join(month_dir, 'MiCASA_vNRT_flux_x3600_y1800_daily_??????_sha256.txt')) #should be in same order
    checksum_candidates = checksum_filename_full + checksum_filename_NRT
    if not checksum_candidates:
        print(f'WARNING: no checksum file found in {month_dir}, skipping')
        continue
    checksum_filename = checksum_candidates[0]

    with open(checksum_filename, 'r') as checksum_file:
        checksum_lines = checksum_file.readlines()

    month_id = "/".join(month_dir.split("/")[-3:-1])
    print(f'{len(daily_filenames)} in month {month_id}[', end='')
    for daily_filename, saved_checksum in zip(daily_filenames, checksum_lines):
        saved_checksum = saved_checksum.strip().split(' ')[0] # only get actual checksum
        file_checksum = popen(f'sha256sum {daily_filename}').read().strip().split(' ')[0] # only get actual checksum
        print('*', end='', flush=True)
        if not file_checksum == saved_checksum:
            print(f'{daily_filename} does not match')
            incorrect_daily_hashes += 1
        else:
            correct_daily_hashes += 1
    print(']')
if incorrect_daily_hashes == 0:
    print(f'All {correct_daily_hashes} daily hashes correct')
else:
    print(f'{incorrect_daily_hashes} failed daily hashes ({correct_daily_hashes} correct)')


incorrect_monthly_hashes = 0
correct_monthly_hashes = 0
for year_dir in montlies_year_dirs:
    monthly_filenames = sorted(glob(join(year_dir, 'MiCASA_v1_flux_x3600_y1800_monthly_??????.nc4')))
    checksum_filenames= sorted(glob(join(year_dir, 'MiCASA_v1_flux_x3600_y1800_monthly_??????_sha256.txt')))
    
    for monthly_filename, checksum_filename in zip(monthly_filenames, checksum_filenames):
        file_checksum = popen(f'sha256sum {monthly_filename}').read().strip().split(' ')[0] # only get actual checksum
        with open(checksum_filename, 'r') as checksum_file:
            saved_checksum = checksum_file.read().strip().split(' ')[0] # only get actual checksum
        if not file_checksum == saved_checksum:
            print(f'{monthly_filename} does not match')
            incorrect_monthly_hashes += 1
        else:
            correct_monthly_hashes += 1

if incorrect_monthly_hashes == 0:
    print(f'All {correct_monthly_hashes} monthly hashes correct')
else:
    print(f'{incorrect_monthly_hashes} failed montly hashes ({correct_monthly_hashes} correct)')
