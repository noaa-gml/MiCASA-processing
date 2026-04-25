#!/usr/bin/env python

from glob import glob
from os.path import join
from os import popen

print('CHECKING HASHES (only 2024 & 2025)')

netcdf_dir = './portal.nccs.nasa.gov'
dailies_month_dirs = sorted(glob(join(netcdf_dir, 'daily/202[4-5]/??/')))
montlies_year_dirs = sorted(glob(join(netcdf_dir, 'monthly/202[4-5]/')))

incorrect_daily_hashes = 0
correct_daily_hashes = 0
for month_dir in dailies_month_dirs:
    daily_filenames_full = sorted(glob(join(month_dir, 'MiCASA_v1_flux_x3600_y1800_daily_????????.nc4')))
    daily_filenames_NRT = sorted(glob(join(month_dir, 'MiCASA_vNRT_flux_x3600_y1800_daily_????????.nc4')))
    daily_filenames = daily_filenames_full + daily_filenames_NRT

    checksum_filename_full = glob(join(month_dir, 'MiCASA_v1_flux_x3600_y1800_daily_??????_sha256.txt')) #should be in same order
    checksum_filename_NRT = glob(join(month_dir, 'MiCASA_vNRT_flux_x3600_y1800_daily_??????_sha256.txt')) #should be in same order
    checksum_filename = (checksum_filename_full + checksum_filename_NRT)[0]

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



