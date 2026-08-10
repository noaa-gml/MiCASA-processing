#!/usr/bin/env python3
"""Unit tests for lib/climatology.py (numpy-only, CI-runnable).

Run:  python3 tests/test_climatology.py
Exits non-zero on any failure.

The baseline-window behaviour is the point of this module: the climatology
prior exists to REMOVE a trend, so the tests assert both that a restricted
window is honoured and that ignoring it gives a materially different answer.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
from climatology import (modulo_month_mean, modulo_month_counts, month_span,
                         month_mid_epoch, parse_yyyymm, seconds_in_month)

_failures = []


def check(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        _failures.append(name)


# --- 1. Unwindowed behaviour is unchanged (compute_clim.py's contract) ------
months = np.tile(np.arange(1, 13), 3)
values = np.empty((36, 2, 5))
for i, m in enumerate(months):
    values[i] = m
clim = modulo_month_mean(values, months)
check("shape is (12, 2, 5)", clim.shape == (12, 2, 5))
check("month m averages to m", np.allclose(clim, np.arange(1, 13)[:, None, None]))

months2 = np.array([1, 1, 2])
values2 = np.array([10.0, 20.0, 7.0]).reshape(3, 1, 1)
clim2 = modulo_month_mean(values2, months2)
check("Jan = mean(10,20) = 15", clim2[0, 0, 0] == 15.0)
check("absent months are NaN", np.all(np.isnan(clim2[2:, 0, 0])))

v = np.array([[[1.0, np.nan]], [[3.0, np.nan]]])
clim3 = modulo_month_mean(v, np.array([1, 1]))
check("NaN-skipping cell = mean(1,3) = 2", clim3[0, 0, 0] == 2.0)
check("all-NaN cell -> NaN", np.isnan(clim3[0, 0, 1]))

try:
    modulo_month_mean(np.zeros((4, 1, 1)), np.array([1, 2, 3]))
    check("mismatched months raises ValueError", False)
except ValueError:
    check("mismatched months raises ValueError", True)

# --- 2. Baseline window ----------------------------------------------------
# 25 years x 12 months; every cell of year y holds the value y.
yrs = np.repeat(np.arange(2001, 2026), 12)
mons = np.tile(np.arange(1, 13), 25)
vals = yrs.astype("float64").reshape(-1, 1, 1) * np.ones((1, 2, 2))

c_base = modulo_month_mean(vals, mons, yrs, 2001, 2020)
check("baseline mean == mean(2001..2020)",
      np.allclose(c_base, np.arange(2001, 2021).mean()))

# The negative control: without the window the answer is different. If the
# window were silently ignored the climatology would carry the very years it
# exists to exclude.
c_all = modulo_month_mean(vals, mons, yrs)
check("unwindowed differs from windowed (window is actually applied)",
      not np.allclose(c_all, c_base))
check("unwindowed mean == mean(2001..2025)",
      np.allclose(c_all, np.arange(2001, 2026).mean()))

# One-sided windows
c_from = modulo_month_mean(vals, mons, yrs, 2021, None)
check("open-ended window (2021-) == mean(2021..2025)",
      np.allclose(c_from, np.arange(2021, 2026).mean()))

# --- 3. Trend removal ------------------------------------------------------
# A pure linear ramp in time: after modulo-month averaging only the within-year
# offset survives, i.e. the 12 monthly values span exactly 11 units.
ramp = np.arange(len(yrs), dtype="float64").reshape(-1, 1, 1) * np.ones((1, 2, 2))
cr = modulo_month_mean(ramp, mons, yrs, 2001, 2020)
spread = cr[:, 0, 0].max() - cr[:, 0, 0].min()
check("linear trend collapses to the seasonal offset only", abs(spread - 11.0) < 1e-9)
check("every calendar month is a single value across years",
      np.allclose(cr[:, 0, 0], cr[:, 0, 1]))

# --- 4. Counts -------------------------------------------------------------
cnt = modulo_month_counts(mons, yrs, 2001, 2020)
check("counts are 20 per calendar month", np.all(cnt == 20))
cnt2 = modulo_month_counts(np.array([1, 1, 2]), np.array([2001, 2002, 2001]),
                           2001, 2001)
check("counts honour the window", cnt2[0] == 1 and cnt2[1] == 1)

try:
    modulo_month_counts(np.array([1, 2]), None, 2001, 2020)
    check("year window without years raises ValueError", False)
except ValueError:
    check("year window without years raises ValueError", True)

# --- 5. Time helpers -------------------------------------------------------
# The archive's own value for 2001-01 -- a synthetic file must be
# indistinguishable from a real one on the time axis.
check("month_mid_epoch(2001,1) matches the archive",
      month_mid_epoch(2001, 1) == 979646400.0)
check("Feb midpoint differs between leap and common years",
      month_mid_epoch(2024, 2) != month_mid_epoch(2023, 2))
check("seconds_in_month leap Feb == 29 days",
      seconds_in_month(2024, 2) == 29 * 86400.0)

check("month_span crosses a year boundary",
      month_span("2020-11", "2021-02") ==
      [(2020, 11), (2020, 12), (2021, 1), (2021, 2)])
check("month_span single month", month_span("2021-07", "2021-07") == [(2021, 7)])
check("month_span length 2020-01..2026-12 is 84",
      len(month_span("2020-01", "2026-12")) == 84)
try:
    month_span("2021-05", "2021-04")
    check("reversed span raises ValueError", False)
except ValueError:
    check("reversed span raises ValueError", True)

check("parse_yyyymm accepts both forms",
      parse_yyyymm("2021-07") == (2021, 7) and parse_yyyymm("202107") == (2021, 7))
try:
    parse_yyyymm("2021-13")
    check("month 13 raises ValueError", False)
except ValueError:
    check("month 13 raises ValueError", True)

if _failures:
    print(f"\n{len(_failures)} FAILED: {', '.join(_failures)}")
    sys.exit(1)
print("\nall climatology tests passed")
