"""lib/budget.py -- global carbon budget arithmetic for MiCASA products.

Area weighting, unit conversion and the NEE sign convention, in one place, so
that a climatology and the product it replaces are always measured the same
way. A difference between two numbers computed by two different routines is
unattributable; this module exists so there is only one routine.

Sign convention throughout: POSITIVE = source to the atmosphere. The MiCASA
1-degree files store NPP positive-as-uptake, so NEE = Rh - NPP. FIRE and FUEL
are separate channels and are never included here.

Pure numpy (no netCDF), so it is unit-testable standalone --
tests/test_budget.py.
"""
import calendar

import numpy as np

R_EARTH = 6371000.0          # m; 4*pi*R^2 = 5.100645e14 m^2
GC_PER_MOL = 12.011          # gC per mol C
FILL_ABS = 1e30              # anything larger in magnitude is a sentinel


def grid_area(nlat=180, nlon=360):
    """Cell area (m^2) on a regular lat/lon grid, shape (nlat, nlon).

    area = R^2 * dlon_rad * (sin(phi_upper) - sin(phi_lower))
    """
    dlon_rad = np.deg2rad(360.0 / nlon)
    edges = np.deg2rad(np.linspace(-90.0, 90.0, nlat + 1))
    band = R_EARTH ** 2 * dlon_rad * (np.sin(edges[1:]) - np.sin(edges[:-1]))
    return np.repeat(band[:, None], nlon, axis=1)


def clean(a):
    """Mask fill/NaN to zero and return a plain float64 ndarray.

    The house -1e34 sentinel must never reach a sum: a single unmasked cell
    swamps a global total by ~20 orders of magnitude, and the result still
    looks like a number.
    """
    if isinstance(a, np.ma.MaskedArray):
        a = a.filled(np.nan)
    a = np.asarray(a, dtype="float64")
    bad = ~np.isfinite(a) | (np.abs(a) > FILL_ABS)
    return np.where(bad, 0.0, a)


def seconds_in_month(year, month):
    return calendar.monthrange(year, month)[1] * 86400.0


def days_in_year(year):
    return 366 if calendar.isleap(year) else 365


def pgc(flux_gc_m2_s, area, seconds):
    """Integrate a gC m-2 s-1 field over `area` and `seconds` -> PgC."""
    return float(np.sum(flux_gc_m2_s * area) * seconds / 1e15)


def nee_from_npp_rh(npp, rh):
    """NEE = Rh - NPP, in the source files' gC m-2 s-1."""
    return clean(rh) - clean(npp)


def mol_to_gc(a):
    """mol m-2 s-1 -> gC m-2 s-1."""
    return clean(a) * GC_PER_MOL
