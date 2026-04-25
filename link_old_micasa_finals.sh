#!/bin/bash
# Symlink the previous year's diurnalized (final) ERA5 outputs into this year's
# tree, so we don't have to re-run diurnalize-ERA5.r for the unchanged history.
# Range: [$MICASA_YEAR_START, $MICASA_YEAR-1] (the year before the one we're
# now processing). Source: $BASE_DIR/$((MICASA_YEAR-1))/MiCASA_v1/ERA5/.

set -e

. "$(dirname "$0")/config.sh"

prev_year=$((MICASA_YEAR - 1))
src_rel="../../../${prev_year}/MiCASA_v1/ERA5"

mkdir -p "${ERA5_DIR}"

for year in $(seq "${MICASA_YEAR_START}" "${prev_year}"); do
    for mon in {1..12}; do
        monf=$(printf '%02d' "$mon")
        ndays=$(cal "$mon" "$year" | awk 'NF {DAYS = $NF}; END {print DAYS}')
        for day in $(seq 1 "$ndays"); do
            dayf=$(printf '%02d' "$day")
            ln -sf "${src_rel}/MiCASA_${MICASA_VERSION}.nee.${year}${monf}${dayf}.nc" \
                   "${ERA5_DIR}/MiCASA_${MICASA_VERSION}.nee.${year}${monf}${dayf}.nc"
        done
        ln -sf "${src_rel}/fluxes_${year}${monf}.nc" \
               "${ERA5_DIR}/fluxes_${year}${monf}.nc"
    done
done
