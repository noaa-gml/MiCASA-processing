#/bin/bash

for year in {2000..2023}; do
    for mon in {1..12}; do
        monf=`printf '%02d' $mon`
        ndays=`cal $mon ${year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
        for day in `seq 1 $ndays`; do
	    dayf=`printf '%02d' $day`
            ln -s "../../../2024/MiCASA_v1/ERA5/MiCASA_v1.nee.${year}${monf}${dayf}.nc" "ERA5/MiCASA_v1.nee.${year}${monf}${dayf}.nc"
	done
	ln -s "../../../2024/MiCASA_v1/ERA5/fluxes_${year}${monf}.nc" "ERA5/fluxes_${year}${monf}.nc"
    done
done
