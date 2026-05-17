#!/bin/sh
# Split each ERA5/fluxes_YYYYMM.nc into per-day MiCASA_<ver>.nee.YYYYMMDD.nc
# files keeping only NEE. NEE = GPP + RESP, positive = source to atm.
#
# Range: years [$MICASA_YEAR_START, $MICASA_YEAR_END],
#        months [$MICASA_MONTH_START, $MICASA_MONTH_END].

set -e

. "$(dirname "$0")/config.sh"
. "$(dirname "$0")/lib/manifest.sh"

calc() {
    awk " function ceiling(x) {print int(x+0.9999999)} \
          function round(x)   {print int(x+0.4999999)} \
          BEGIN{OFMT = \"%.12g\"; print $* }"
}

# Run manifest (lib/manifest.sh): record start now, ok/fail on exit.
_ds_t0=$(date +%s)
_ds_detail="years ${MICASA_YEAR_START}-${MICASA_YEAR_END} months ${MICASA_MONTH_START}-${MICASA_MONTH_END}"
manifest_record daysplitter.sh start - "$_ds_detail"
trap '_ds_rc=$?; [ "$_ds_rc" -ne 0 ] && manifest_record daysplitter.sh fail "$(($(date +%s)-_ds_t0))" "exit $_ds_rc"' EXIT

for yr in $(seq "${MICASA_YEAR_START}" "${MICASA_YEAR_END}"); do
    for mon in $(seq "${MICASA_MONTH_START}" "${MICASA_MONTH_END}"); do
        monf=$(printf '%02d' "$mon")
        srcfile="${ERA5_DIR}/fluxes_${yr}${monf}.nc"
        if [ ! -e "$srcfile" ]; then
            echo "Skipping missing $srcfile"
            continue
        fi
        nslots=$(ncdump -h "${srcfile}" | grep UNLIMITED | sed -e 's/.*(//' -e 's/currently.*//')
        ndays=$(calc "$nslots/24")
        echo -n "$srcfile: $ndays days ["

        for day in $(seq 1 "$ndays"); do
            dayf=$(printf '%02d' "$day")
            daym1=$(calc "$day - 1")
            min=$(calc "${daym1}*24")
            max=$(calc "$min + 23")
            targetfile="${ERA5_DIR}/MiCASA_${MICASA_VERSION}.nee.${yr}${monf}${dayf}.nc"
            # ncks copies the source file's global attributes (so the daily
            # file inherits diurnalize-ERA5.r's CF/ACDD provenance) and
            # appends its own command to `history`. The --gaa markers below
            # record that this file is a single-day NEE subset.
            ncks --gaa "daily_split_from=fluxes_${yr}${monf}.nc" \
                 --gaa "daily_split_tool=daysplitter.sh" \
                 -v NEE -O -d time,"$min","$max" "$srcfile" "$targetfile"
            echo -n '.'
        done
        echo ']'
    done
done

manifest_record daysplitter.sh ok "$(($(date +%s)-_ds_t0))" "$_ds_detail"
trap - EXIT
