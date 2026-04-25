#!/bin/bash

cd daily_1x1
for yr in 2000 2024; do
    for mon in `seq 1 12`; do
	monf=`printf '%02d' $mon`

	ndays=`cal $mon ${yr} | awk 'NF {DAYS = $NF}; END {print DAYS}'`

	for day in `seq 1 $ndays`;do
	    dayf=`printf "%02d" $day`

	    targetfile=MiCASA_v1_flux_x360_y180_daily_${yr}${monf}${dayf}.nc
	    if [ ! -e $targetfile ]; then
		ln -s MiCASA_v1_flux_x360_y180_daily_0000${monf}${dayf}.nc $targetfile
	    else
		echo File exists $targetfile
	    fi
	done
    done
done


