#!/usr/bin/env python3
"""Verify SHA-256 checksums of MiCASA raw .nc4 files against the
per-directory aggregate manifests NCCS publishes alongside.

NCCS layout (as of 2026):
  daily/<YYYY>/<MM>/MiCASA_<v|vNRT>_flux_x3600_y1800_daily_<YYYYMM>_sha256.txt
  monthly/<YYYY>/MiCASA_<v|vNRT>_flux_x3600_y1800_monthly_<YYYYMM>_sha256.txt

Each manifest is a plain `<sha256>  <filename>` per line, listing the
.nc4 files NCCS produced for that directory. v1 and vNRT have separate
manifests (and separate .nc4 files) that co-exist in the same dir.

Year range comes from $MICASA_YEAR_START / $MICASA_YEAR_END, falling
back to $MICASA_YEAR (single year) if range isn't set, or to
2001..(current calendar year) if neither is set.

Exits non-zero if any .nc4 hash mismatched its manifest entry. Files
with no manifest entry are reported as warnings but don't fail the
run -- a manifest can lag a recently-published file, and it's better
for download.sh to keep moving.

Replaces an earlier version that paired v1 dailies against v1
manifests but also paired vNRT dailies against the *same v1 manifest*
via a single-list zip(), which produced spurious mismatches whenever
both streams co-existed. Also adds vNRT monthly support (the prior
version only handled v1 monthlies).
"""

from __future__ import annotations

import datetime
import hashlib
import os
import re
import sys
from glob import glob
from os.path import join


def year_range_from_env() -> list[int]:
    y_start = os.environ.get("MICASA_YEAR_START")
    y_end   = os.environ.get("MICASA_YEAR_END")
    y_one   = os.environ.get("MICASA_YEAR")
    if y_start and y_end:
        return list(range(int(y_start), int(y_end) + 1))
    if y_one:
        return [int(y_one)]
    return list(range(2001, datetime.date.today().year + 1))


def parse_manifest(path: str) -> dict[str, str]:
    """Return {filename: sha256_hex} from one _sha256.txt manifest."""
    out: dict[str, str] = {}
    with open(path) as fp:
        for line in fp:
            parts = line.strip().split()
            if len(parts) >= 2:
                # parts[0] = hash, parts[1] = filename (may have leading ./ etc.)
                out[os.path.basename(parts[1])] = parts[0].lower()
    return out


def merge_manifests(directory: str) -> dict[str, str]:
    """Merge every _sha256.txt manifest found in *directory* into a single
    {filename: hash} map. v1 and vNRT manifests merge cleanly because
    the filenames they reference are disjoint."""
    merged: dict[str, str] = {}
    for mp in sorted(glob(join(directory, "MiCASA_*_sha256.txt"))):
        for fn, hx in parse_manifest(mp).items():
            merged[fn] = hx
    return merged


def sha256_file(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fp:
        for block in iter(lambda: fp.read(chunk), b""):
            h.update(block)
    return h.hexdigest().lower()


def verify_dir(directory: str, glob_pat: str) -> tuple[int, int, int, list[str]]:
    """Verify every .nc4 in *directory* matching *glob_pat*. Returns
    (n_verified, n_no_record, n_mismatch, mismatch_lines)."""
    manifest = merge_manifests(directory)
    files = sorted(glob(join(directory, glob_pat)))
    n_verified = n_no_record = n_mismatch = 0
    mismatches: list[str] = []
    for f in files:
        bn = os.path.basename(f)
        expected = manifest.get(bn)
        if expected is None:
            n_no_record += 1
            continue
        actual = sha256_file(f)
        if actual != expected:
            n_mismatch += 1
            mismatches.append(f"{f}: expected {expected[:12]}.. got {actual[:12]}..")
        else:
            n_verified += 1
    return n_verified, n_no_record, n_mismatch, mismatches


def main() -> int:
    years = year_range_from_env()
    print(f"CHECKING SHA-256 (years: {years[0]}..{years[-1]})")
    portal = "./portal.nccs.nasa.gov"
    if not os.path.isdir(portal):
        print(f"ERROR: {portal} not present (run from MiCASA working dir)")
        return 2

    total_verified = total_no_record = total_mismatch = 0
    all_mismatches: list[str] = []

    # Dailies: <portal>/daily/YYYY/MM/
    for y in years:
        month_dirs = sorted(glob(join(portal, f"daily/{y:04d}/??/")))
        for d in month_dirs:
            short = "/".join(d.rstrip("/").split("/")[-3:])
            v, nr, m, mm = verify_dir(d, "MiCASA_*_flux_x3600_y1800_daily_*.nc4")
            total_verified += v
            total_no_record += nr
            total_mismatch += m
            all_mismatches.extend(mm)
            tag = " ".join(filter(None, [
                f"verified={v}",
                f"no_record={nr}" if nr else "",
                f"MISMATCH={m}" if m else "",
            ]))
            print(f"  daily {short}: {tag}")

    # Monthlies: <portal>/monthly/YYYY/
    for y in years:
        d = join(portal, f"monthly/{y:04d}/")
        if not os.path.isdir(d):
            continue
        v, nr, m, mm = verify_dir(d, "MiCASA_*_flux_x3600_y1800_monthly_*.nc4")
        total_verified += v
        total_no_record += nr
        total_mismatch += m
        all_mismatches.extend(mm)
        tag = " ".join(filter(None, [
            f"verified={v}",
            f"no_record={nr}" if nr else "",
            f"MISMATCH={m}" if m else "",
        ]))
        if v + nr + m > 0:
            print(f"  monthly {y}: {tag}")

    print()
    print(f"SUMMARY: {total_verified} verified, "
          f"{total_no_record} files lacking a manifest entry, "
          f"{total_mismatch} MISMATCH")
    if all_mismatches:
        print("\nMismatch detail:")
        for line in all_mismatches[:20]:
            print(f"  {line}")
        if len(all_mismatches) > 20:
            print(f"  ... and {len(all_mismatches) - 20} more")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
