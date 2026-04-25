#/bin/bash
# Time-stamp: <orion-login-3.hpc.msstate.edu:/work2/noaa/co2/GFED-CASA/2025/MiCASA_v1/link_vNRT_to_v1_for_CT2025.sh: 10 Jul 2025 (Thu) 15:22:32 UTC>

cd daily_1x1
year=2025
for mon in {1..12}; do
    monf=`printf '%02d' $mon`
    ndays=`cal $mon ${year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
    for day in `seq 1 $ndays`; do
	dayf=`printf '%02d' $day`
	target=MiCASA_v1_flux_x360_y180_daily_${year}${monf}${dayf}.nc
	if [ -e ${target} ]; then
	    echo "Skipping existing $target"
	else
	    src=MiCASA_vNRT_flux_x360_y180_daily_${year}${monf}${dayf}.nc
	    if [ -e ${src} ]; then
		ln -s $src $target
		echo "Linked $target"
	    else
		echo "WARNING source does not exist $src"
	    fi
	fi

    done

done
