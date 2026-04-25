#!/bin/sh

calc() {
    awk " function ceiling(x) {print int(x+0.9999999)} function round(x) {print int(x+0.4999999)} BEGIN{OFMT = \"%.12g\"; print $* }"
}

#for yr in `seq 2000 2024`;do
#for yr in `seq 2024 2024`;do
for yr in 2025;do

    for mon in `seq 4 12`; do

	monf=`printf '%02d' $mon`
	srcfile=ERA5/fluxes_${yr}${monf}.nc
	nslots=`ncdump -h ${srcfile} | grep UNLIMITED | sed -e 's/.*(//' -e 's/currently.*//'`
	ndays=`calc $nslots/24`
	echo -n "$srcfile: $ndays days ["

	for day in `seq 1 $ndays`; do
	    dayf=`printf '%02d' $day`
	    daym1=`calc $day - 1`
	    min=`calc ${daym1}*24`
	    max=`calc $min + 23`
	    targetfile=ERA5/MiCASA_v1.nee.${yr}${monf}${dayf}.nc
	    ncks -v NEE -O -d time,$min,$max $srcfile $targetfile
	    echo -n '.'
	done
	echo ']'

    done

done
