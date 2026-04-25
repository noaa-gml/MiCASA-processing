#!/bin/bash

month_names=("" "Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")

for mon in `seq 1 12`; do
    monf=`printf '%02d' $mon`

    # choose 2004 to intentionally include 29 Feb
    ndays=`cal $mon 2004 | awk 'NF {DAYS = $NF}; END {print DAYS}'`

    echo -n "${month_names[$mon]} ["
    for day in `seq 1 $ndays`;do
        dayf=`printf "%02d" $day`

	# Bug found 17 July 2025. If _0000 (clim year) files already
	# exist in the daily_1x1 directory, they would be included in
        # the climatological average.
#	fls=`ls daily_1x1/MiCASA_v1_flux_x360_y180_daily_????${monf}${dayf}.nc`
        fls=`ls daily_1x1/MiCASA_v1_flux_x360_y180_daily_2???${monf}${dayf}.nc`
	
        echo ${fls[*]}| ncea -O -o daily_1x1/MiCASA_v1_flux_x360_y180_daily_0000${monf}${dayf}.nc
        echo -n "*"
    done
    echo "]"
    
done
