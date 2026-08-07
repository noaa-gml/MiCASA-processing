#!/usr/bin/env python3
"""Fail if our local raw mirror is STALE relative to the NCCS portal.

Upstream silently *republishes* raw files (NASA reprocesses; the filename does
not change). Our mirror is pinned by `wget --no-clobber`, so a republished file
is never re-fetched -- and because the per-directory `_sha256.txt` manifest is
skipped by --no-clobber too, the stale manifest validates the stale data and
`check_hashes.py` goes green. The result is a silently outdated product.

    2025-06-17: NASA republished 47 daily vNRT files (2025-01..2025-05),
    fixing a 35-day spurious FIRE smear (2025-01-07..02-10, ~20-29x the
    global fire flux). Our mirror, downloaded 2025-05-08, kept the broken
    vintage; the June-2026 rebuild re-ingested it and shipped it.

This check compares every local raw file against the portal directory listing
(exact filename, size, publication date) and exits non-zero if any local file
differs or is missing. It needs outbound HTTPS -- on Orion that means a LOGIN
or DTN node, not a compute node. With no network it exits 0 with a warning
(it must not break offline reprocessing), so run it where the network is.

Usage:
    ./check_upstream_fresh.py [--version v1|vNRT|both] [--quiet]

Environment (as exported by config.sh):
    MICASA_YEAR_START / MICASA_YEAR_END   years to audit (default: MICASA_YEAR)
    MICASA_VERSION                        stream to audit (default: v1)
    RAW_SRC_DIR                           local mirror root
    PORTAL_URL_BASE                       portal root

Exit codes:
    0  mirror matches upstream (or upstream unreachable / not yet published)
    1  at least one local file is STALE or MISSING  <-- refresh before building
    2  usage / configuration error
"""
import os
import re
import sys
import urllib.request
import urllib.error

ROW = re.compile(
    r'href="(?P<f>MiCASA_[^"]+\.nc4)".*?'
    r'class="size">(?P<s>[^<]+)</td><td class="date">(?P<d>[^<]+)</td>',
    re.S)

PORTAL = os.environ.get(
    "PORTAL_URL_BASE",
    "https://portal.nccs.nasa.gov/datashare/gmao/geos_carb/MiCASA")
RAW = os.environ.get("RAW_SRC_DIR", "portal.nccs.nasa.gov")


def mib(nbytes):
    """Render a byte count the way the portal listing does ('17.0 MiB')."""
    return f"{nbytes / 1048576:.1f} MiB"


def listing(url):
    """Return [(filename, size_str, date_str)] or None if unreachable/absent."""
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            html = r.read().decode("utf8", "ignore")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None          # this stream has not published this period
        raise
    except Exception:
        return "UNREACHABLE"
    return [(m.group("f"), m.group("s").strip(), m.group("d").strip())
            for m in ROW.finditer(html)]


def audit(version, kind, year, month=None):
    """Compare one upstream directory against the local mirror."""
    sub = f"{year}/{month}" if month else f"{year}"
    rows = listing(f"{PORTAL}/{version}/netcdf/{kind}/{sub}/")
    if rows is None:
        return None
    if rows == "UNREACHABLE":
        return "UNREACHABLE"
    local_dir = os.path.join(RAW, kind, sub)
    stale, missing, ok = [], [], 0
    for fn, size, date in rows:
        p = os.path.join(local_dir, fn)
        if not os.path.exists(p):
            missing.append((fn, date))
        elif mib(os.path.getsize(p)) != size:
            stale.append((fn, mib(os.path.getsize(p)), size, date))
        else:
            ok += 1
    return ok, stale, missing


def main(argv):
    version = os.environ.get("MICASA_VERSION", "v1")
    quiet = "--quiet" in argv
    if "--version" in argv:
        try:
            version = argv[argv.index("--version") + 1]
        except IndexError:
            print("--version needs a value (v1|vNRT|both)", file=sys.stderr)
            return 2
    versions = ["v1", "vNRT"] if version == "both" else [version]
    if any(v not in ("v1", "vNRT") for v in versions):
        print(f"bad --version '{version}' (v1|vNRT|both)", file=sys.stderr)
        return 2

    y0 = int(os.environ.get("MICASA_YEAR_START",
                            os.environ.get("MICASA_YEAR", "2001")))
    y1 = int(os.environ.get("MICASA_YEAR_END",
                            os.environ.get("MICASA_YEAR", y0)))

    n_stale = n_missing = n_ok = 0
    unreachable = False
    for v in versions:
        for year in range(y0, y1 + 1):
            targets = [("daily", f"{m:02d}") for m in range(1, 13)]
            targets.append(("monthly", None))
            for kind, month in targets:
                r = audit(v, kind, year, month)
                if r is None:
                    continue
                if r == "UNREACHABLE":
                    unreachable = True
                    continue
                ok, stale, missing = r
                n_ok += ok
                n_stale += len(stale)
                n_missing += len(missing)
                for fn, have, want, date in stale:
                    print(f"STALE   {fn}  local={have} upstream={want} "
                          f"(republished {date})")
                if missing and not quiet:
                    for fn, date in missing[:5]:
                        print(f"MISSING {fn}  (published {date})")
                    if len(missing) > 5:
                        print(f"MISSING ... and {len(missing) - 5} more in "
                              f"{v}/{kind}/{year}{'/' + month if month else ''}")

    if unreachable:
        print("WARN: portal unreachable (no outbound HTTPS?) -- freshness "
              "NOT verified. Run this on a login/DTN node.")
        return 0

    print(f"\nupstream freshness: {n_ok} match, {n_stale} STALE, "
          f"{n_missing} missing  [{','.join(versions)} {y0}..{y1}]")
    if n_stale or n_missing:
        print("\n*** The local raw mirror does not match upstream. Ingesting it "
              "will\n*** reproduce a superseded vintage. Refresh first:\n"
              "***     MICASA_REFRESH=1 ./download.sh      # re-fetches "
              "republished files\n*** then re-ingest the affected years.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
