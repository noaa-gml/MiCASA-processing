#!/usr/bin/env python3
"""Unit tests for check_bounds.flux_to_tgc_per_year (numpy-only, CI-runnable).

check_bounds.py imports xarray lazily inside main(), so importing the module
for these tests needs only numpy.

Run:  python3 tests/test_check_bounds.py
Exits non-zero on any failure.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from check_bounds import flux_to_tgc_per_year, EARTH_AREA, SECONDS_IN_YEAR

_failures = []


def check(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        _failures.append(name)


# --- zero / sign / linearity ----------------------------------------------
check("zero flux -> 0 TgC/yr", flux_to_tgc_per_year(0.0) == 0.0)
check("negative flux -> negative TgC/yr", flux_to_tgc_per_year(-1.0) < 0.0)
check("doubling the flux doubles the result",
      abs(flux_to_tgc_per_year(2.0) - 2.0 * flux_to_tgc_per_year(1.0)) < 1e-6)

# --- conversion is mean_flux * EARTH_AREA * SECONDS_IN_YEAR / 1e15 ---------
expected = 1.0 * EARTH_AREA * SECONDS_IN_YEAR / 1e15
check("flux_to_tgc_per_year(1.0) == EARTH_AREA * SECONDS_IN_YEAR / 1e15",
      abs(flux_to_tgc_per_year(1.0) - expected) <= 1e-9 * expected)

# --- constants are physically sane ----------------------------------------
# Earth surface area ~ 5.1e14 m^2; a (Julian) year ~ 3.156e7 s.
check("EARTH_AREA is ~5.1e14 m^2", 5.0e14 < EARTH_AREA < 5.2e14)
check("SECONDS_IN_YEAR is ~3.156e7 s", 3.10e7 < SECONDS_IN_YEAR < 3.20e7)

# --- order of magnitude: 1 gC/m2/s globally for a year ~ 1.6e7 TgC --------
# Guards against a wrong scale constant (e.g. /1e12 instead of /1e15).
v = flux_to_tgc_per_year(1.0)
check("1 gC/m2/s over Earth for a year is ~1.6e7 TgC/yr", 1.5e7 < v < 1.7e7)

if _failures:
    print(f"\n{len(_failures)} FAILED: {', '.join(_failures)}")
    sys.exit(1)
print("\nall check_bounds tests passed")
