#!/bin/sh


sh download.sh
Rscript check_daily_downloads.r
python check_hashes.py
sh check_unchanged.sh

