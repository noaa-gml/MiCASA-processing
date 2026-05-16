#!/usr/bin/env python3
"""Unit tests for compute_clim.modulo_month_mean (numpy-only, CI-runnable).

Run:  python3 tests/test_compute_clim.py
Exits non-zero on any failure.
"""
import os
import sys

import numpy as np

# Import modulo_month_mean from compute_clim.py one level up. compute_clim's
# `main()` imports xarray lazily, so this import needs only numpy.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from compute_clim import modulo_month_mean

_failures = []


def check(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        _failures.append(name)


# --- 1. Each calendar month recovers its own value -------------------------
# 3 years x 12 months; every cell of month m holds the value m.
months = np.tile(np.arange(1, 13), 3)
values = np.empty((36, 2, 5))
for i, m in enumerate(months):
    values[i] = m
clim = modulo_month_mean(values, months)
check("shape is (12, 2, 5)", clim.shape == (12, 2, 5))
check("month m averages to m", np.allclose(clim, np.arange(1, 13)[:, None, None]))

# --- 2. Genuine averaging across years -------------------------------------
# Two Januaries holding 10 and 20 -> clim Jan = 15.
months2 = np.array([1, 1, 2])
values2 = np.array([10.0, 20.0, 7.0]).reshape(3, 1, 1)
clim2 = modulo_month_mean(values2, months2)
check("Jan = mean(10,20) = 15", clim2[0, 0, 0] == 15.0)
check("Feb = 7", clim2[1, 0, 0] == 7.0)
check("absent months (Mar..Dec) are NaN", np.all(np.isnan(clim2[2:, 0, 0])))

# --- 3. NaN handling -------------------------------------------------------
# One Jan cell is NaN; nanmean must skip it. An all-NaN cell -> NaN.
v = np.array([[[1.0, np.nan]],
              [[3.0, np.nan]]])          # 2 Januaries, cells [ok] and [all-NaN]
clim3 = modulo_month_mean(v, np.array([1, 1]))
check("NaN-skipping cell = mean(1,3) = 2", clim3[0, 0, 0] == 2.0)
check("all-NaN cell -> NaN", np.isnan(clim3[0, 0, 1]))

# --- 4. Uneven month counts (the 2026-partial case) ------------------------
# 3 Januaries, 2 Februaries -- a trailing partial year. Means still correct.
months4 = np.array([1, 1, 1, 2, 2])
values4 = np.array([2.0, 4.0, 6.0, 100.0, 200.0]).reshape(5, 1, 1)
clim4 = modulo_month_mean(values4, months4)
check("uneven: Jan = mean(2,4,6) = 4", clim4[0, 0, 0] == 4.0)
check("uneven: Feb = mean(100,200) = 150", clim4[1, 0, 0] == 150.0)

# --- 5. Length-mismatch is rejected ----------------------------------------
try:
    modulo_month_mean(np.zeros((4, 1, 1)), np.array([1, 2, 3]))
    check("mismatched months raises ValueError", False)
except ValueError:
    check("mismatched months raises ValueError", True)

if _failures:
    print(f"\n{len(_failures)} FAILED: {', '.join(_failures)}")
    sys.exit(1)
print("\nall compute_clim tests passed")
