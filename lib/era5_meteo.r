## lib/era5_meteo.r -- ERA5 hourly-meteo path-resolution helpers.
##
## Pure functions: no ct.setup(), no CT helper library, no script globals,
## so they can be unit-tested standalone -- see tests/test_era5_meteo.r.
##
## Used by diurnalize-ERA5.r to resolve each day's meteo from the primary
## ERA5 tree, falling back to the FastTrack (ea_0005) tree for the NRT
## trailing window. See docs/PROPOSALS.md item (12).

## Relative path of one (yr, mon, day, var) ERA5 file, built from a
## filename template containing YYYY / MM / DD / VVV placeholders.
era5.relpath <- function(template, yr, mon, day, varnm) {
  e5nm <- gsub("YYYY", sprintf("%d",   yr),  template)
  e5nm <- gsub("MM",   sprintf("%02d", mon), e5nm)
  e5nm <- gsub("DD",   sprintf("%02d", day), e5nm)
  gsub("VVV", varnm, e5nm)
}

## Resolve a day to the first meteo tree (a named entry of `era5dirs`)
## that holds ALL `varnms` for it. Returns that era5dirs name, or NA if
## no tree has the complete set.
resolve.era5.source <- function(era5dirs, template, yr, mon, day, varnms) {
  for (src in names(era5dirs)) {
    paths <- file.path(era5dirs[[src]],
                       vapply(varnms,
                              function(v) era5.relpath(template, yr, mon, day, v),
                              character(1)))
    if (all(file.exists(paths))) return(src)
  }
  NA_character_
}

## Compact run-length encoding of a per-day source vector, e.g.
## "primary:1-30 fasttrack:31". `days` are the day numbers to encode;
## `srcvec` is indexed by day number.
encode.day.runs <- function(days, srcvec) {
  if (length(days) == 0) return("")
  parts <- character(0)
  for (s in unique(srcvec[days])) {
    ds <- sort(days[srcvec[days] == s])
    runs <- character(0); i <- 1
    while (i <= length(ds)) {
      j <- i
      while (j < length(ds) && ds[j + 1] == ds[j] + 1) j <- j + 1
      runs <- c(runs, if (j > i) sprintf("%d-%d", ds[i], ds[j])
                      else        sprintf("%d", ds[i]))
      i <- j + 1
    }
    parts <- c(parts, sprintf("%s:%s", s, paste(runs, collapse = ",")))
  }
  paste(parts, collapse = " ")
}
