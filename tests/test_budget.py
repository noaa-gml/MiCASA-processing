#!/usr/bin/env python3
"""Unit tests for lib/budget.py (numpy-only, CI-runnable).

Run:  python3 tests/test_budget.py
Exits non-zero on any failure.

Includes the negative control that matters most: an unmasked -1e34 sentinel
must produce an absurd answer, so that the masking in clean() is demonstrably
doing something rather than being decorative.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
from budget import (grid_area, clean, pgc, nee_from_npp_rh, mol_to_gc,
                    seconds_in_month, days_in_year, R_EARTH, GC_PER_MOL)

_failures = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}{('  -- ' + detail) if detail else ''}")
    if not ok:
        _failures.append(name)


# --- Area ------------------------------------------------------------------
area = grid_area()
check("grid_area shape", area.shape == (180, 360))
check("total area == 4 pi R^2",
      abs(area.sum() - 4 * np.pi * R_EARTH ** 2) < 1.0,
      f"{area.sum():.6e}")
check("area is symmetric about the equator",
      np.allclose(area[:90], area[90:][::-1]))
check("polar cells are smaller than equatorial",
      area[0, 0] < area[90, 0])
check("coarse grid also closes",
      abs(grid_area(18, 36).sum() - 4 * np.pi * R_EARTH ** 2) < 1.0)

# --- clean() ---------------------------------------------------------------
a = np.array([1.0, np.nan, -1e34, 2.0, np.inf])
c = clean(a)
check("clean zeroes NaN, inf and the -1e34 sentinel",
      np.array_equal(c, np.array([1.0, 0.0, 0.0, 2.0, 0.0])))
m = np.ma.masked_array([1.0, 5.0], mask=[False, True])
check("clean handles masked arrays", np.array_equal(clean(m), np.array([1.0, 0.0])))

# --- pgc -------------------------------------------------------------------
uni = np.ones((180, 360))
want = area.sum() * 31 * 86400.0 / 1e15
check("uniform 1 gC m-2 s-1 over a 31-day month",
      abs(pgc(uni, area, seconds_in_month(2021, 1)) - want) < 1e-9 * want,
      f"{pgc(uni, area, seconds_in_month(2021,1)):.6e} PgC")

# NEGATIVE CONTROL: without masking, one sentinel cell destroys the total.
contaminated = np.ones((180, 360))
contaminated[90, 180] = -1.0e34
masked_total = pgc(clean(contaminated), area, seconds_in_month(2021, 1))
raw_total = pgc(contaminated, area, seconds_in_month(2021, 1))
one_cell = area[90, 180] * 31 * 86400.0 / 1e15
check("masked total drops exactly the sentinel cell",
      abs(masked_total - (want - one_cell)) < 1e-6)
check("negative control: unmasked sentinel gives an absurd total",
      abs(raw_total) > 1e6, f"{raw_total:.3e} PgC")

# --- Sign convention -------------------------------------------------------
npp = np.zeros((4, 4)); npp[1, 1] = 2.0
rh = np.zeros((4, 4)); rh[1, 1] = 0.5
check("NEE = Rh - NPP is negative where uptake dominates",
      nee_from_npp_rh(npp, rh)[1, 1] == -1.5)
npp2 = np.zeros((2, 2)); rh2 = np.zeros((2, 2)); rh2[0, 0] = 3.0
check("NEE positive where respiration dominates",
      nee_from_npp_rh(npp2, rh2)[0, 0] == 3.0)

# --- Units -----------------------------------------------------------------
check("mol -> gC uses 12.011", mol_to_gc(np.array([1.0]))[0] == GC_PER_MOL)
check("mol_to_gc also cleans sentinels", mol_to_gc(np.array([-1e34]))[0] == 0.0)

# --- Calendar --------------------------------------------------------------
check("Feb 2024 has 29 days", seconds_in_month(2024, 2) == 29 * 86400.0)
check("Feb 2023 has 28 days", seconds_in_month(2023, 2) == 28 * 86400.0)
check("2024 is a leap year", days_in_year(2024) == 366)
check("2100 is not a leap year", days_in_year(2100) == 365)

if _failures:
    print(f"\n{len(_failures)} FAILED: {', '.join(_failures)}")
    sys.exit(1)
print("\nall budget tests passed")
