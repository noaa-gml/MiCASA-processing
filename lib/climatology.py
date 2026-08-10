"""lib/climatology.py -- modulo-month climatology helpers.

The pure-arithmetic core behind two pipeline stages:

  compute_clim.py             the NPP/Rh day-of-year climatology that
                              diurnalize-ERA5.r falls back on for months with
                              no published monthly file.

  make_climatology_series.py  the CLIMATOLOGY PRIOR: a synthetic monthly series
                              in which every calendar month carries its
                              multi-year mean, so the diurnalized product has
                              climatological magnitude but real-meteorology
                              sub-daily structure.

Both need "the mean of each calendar month", and the second also needs to
restrict that mean to a baseline year window. Kept free of netCDF/xarray so it
is unit-testable standalone -- tests/test_climatology.py.

Standard library + numpy only.
"""
import calendar
import datetime as _dt
import warnings

import numpy as np


def _select(months, years=None, year_start=None, year_end=None):
    """Boolean mask over the time axis for the requested baseline window."""
    months = np.asarray(months)
    sel = np.ones(months.shape[0], dtype=bool)
    if years is not None:
        years = np.asarray(years)
        if years.shape[0] != months.shape[0]:
            raise ValueError("years length must match months length")
        if year_start is not None:
            sel &= years >= year_start
        if year_end is not None:
            sel &= years <= year_end
    elif year_start is not None or year_end is not None:
        raise ValueError("year_start/year_end need `years`")
    return sel


def modulo_month_mean(values, months, years=None,
                      year_start=None, year_end=None):
    """Mean of each calendar month across the time axis.

    values : ndarray, axis 0 is time.
    months : 1-D int array (calendar month 1..12), length == values.shape[0].
    years  : optional 1-D int array parallel to `months`; required to use
             year_start / year_end.
    year_start, year_end : inclusive baseline window. Samples outside it are
             excluded from the mean. Omit both to average the whole record
             (the historical behaviour, unchanged).

    returns: ndarray, axis 0 length 12 (Jan..Dec), other axes preserved.
             NaN / fill values are skipped per cell; a month with no data
             (or an all-missing cell) yields NaN.

    Equivalent to PyFerret's <var>[GT=MONTH_IRREG@MOD] and to xarray's
    `da.groupby("time.month").mean("time")`, optionally restricted in year.
    """
    values = np.asarray(values, dtype="float64")
    months = np.asarray(months)
    if months.shape[0] != values.shape[0]:
        raise ValueError("months length must match values' time axis")
    base = _select(months, years, year_start, year_end)
    out = np.full((12,) + values.shape[1:], np.nan, dtype="float64")
    for m in range(1, 13):
        sel = base & (months == m)
        if np.any(sel):
            # An all-NaN cell (e.g. ocean) yields NaN via nanmean -- that is
            # the intended result, so silence the "Mean of empty slice"
            # RuntimeWarning rather than letting it spam the log.
            with warnings.catch_warnings(), np.errstate(invalid="ignore"):
                warnings.simplefilter("ignore", RuntimeWarning)
                out[m - 1] = np.nanmean(values[sel], axis=0)
    return out


def modulo_month_counts(months, years=None, year_start=None, year_end=None):
    """Samples contributing to each calendar month, as a length-12 int array.

    A climatology built from an uneven number of years per calendar month is
    not wrong, but it is rarely intended -- callers should report this.
    """
    months = np.asarray(months)
    base = _select(months, years, year_start, year_end)
    return np.array([int((base & (months == m)).sum()) for m in range(1, 13)])


def parse_yyyymm(s):
    """'2021-07' or '202107' -> (2021, 7)."""
    s = str(s).strip()
    if "-" in s:
        y, m = s.split("-", 1)
    else:
        y, m = s[:4], s[4:]
    y, m = int(y), int(m)
    if not 1 <= m <= 12:
        raise ValueError("month out of range in %r" % (s,))
    return y, m


def month_span(start, end):
    """Inclusive list of (year, month) from `start` to `end`.

    Accepts (y, m) tuples or 'YYYY-MM' strings.
    """
    if not isinstance(start, tuple):
        start = parse_yyyymm(start)
    if not isinstance(end, tuple):
        end = parse_yyyymm(end)
    if end < start:
        raise ValueError("end %r precedes start %r" % (end, start))
    out = []
    y, m = start
    while (y, m) <= end:
        out.append((y, m))
        m += 1
        if m == 13:
            y, m = y + 1, 1
    return out


def month_mid_epoch(year, month):
    """Seconds since 1970-01-01 at the exact midpoint of the month.

    Matches the convention the MiCASA 1-degree monthly files already use
    (2001-01 -> 979646400.0 = 2001-01-16 12:00:00 UTC), so a synthetic series
    is indistinguishable from a real one to any time-axis reader.
    """
    start = _dt.datetime(year, month, 1, tzinfo=_dt.timezone.utc)
    ndays = calendar.monthrange(year, month)[1]
    return start.timestamp() + ndays * 86400.0 / 2.0


def seconds_in_month(year, month):
    return calendar.monthrange(year, month)[1] * 86400.0
