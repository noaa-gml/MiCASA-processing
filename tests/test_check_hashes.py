#!/usr/bin/env python3
"""Unit tests for check_hashes.py helpers (standard library only, CI-runnable).

check_hashes.py imports only the standard library, so importing it for these
tests pulls in nothing extra.

Run:  python3 tests/test_check_hashes.py
Exits non-zero on any failure.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import check_hashes

_failures = []


def check(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        _failures.append(name)


_ABC = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
_EMPTY = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

with tempfile.TemporaryDirectory() as tmp:
    # --- sha256_file -------------------------------------------------------
    f_abc = os.path.join(tmp, "abc")
    with open(f_abc, "wb") as fh:
        fh.write(b"abc")
    check("sha256_file of 'abc' is the known digest",
          check_hashes.sha256_file(f_abc) == _ABC)
    f_empty = os.path.join(tmp, "empty")
    open(f_empty, "wb").close()
    check("sha256_file of an empty file is the known digest",
          check_hashes.sha256_file(f_empty) == _EMPTY)

    # --- parse_manifest ----------------------------------------------------
    man = os.path.join(tmp, "MiCASA_v1_flux_x3600_y1800_daily_202001_sha256.txt")
    with open(man, "w") as fh:
        fh.write("aaa111  MiCASA_v1_flux_x3600_y1800_daily_20200101.nc4\n"
                 "bbb222  ./MiCASA_v1_flux_x3600_y1800_daily_20200102.nc4\n"
                 "shortline_no_filename\n")
    pm = check_hashes.parse_manifest(man)
    check("parse_manifest reads <hash> <file> lines",
          pm.get("MiCASA_v1_flux_x3600_y1800_daily_20200101.nc4") == "aaa111")
    check("parse_manifest basenames a ./-prefixed filename",
          pm.get("MiCASA_v1_flux_x3600_y1800_daily_20200102.nc4") == "bbb222")
    check("parse_manifest skips malformed lines (2 entries kept)", len(pm) == 2)

    # --- merge_manifests: v1 + vNRT manifests in one directory -------------
    man2 = os.path.join(tmp, "MiCASA_vNRT_flux_x3600_y1800_daily_202001_sha256.txt")
    with open(man2, "w") as fh:
        fh.write("ccc333  MiCASA_vNRT_flux_x3600_y1800_daily_20200103.nc4\n")
    merged = check_hashes.merge_manifests(tmp)
    check("merge_manifests merges the v1 and vNRT manifests",
          len(merged) == 3 and
          merged.get("MiCASA_vNRT_flux_x3600_y1800_daily_20200103.nc4") == "ccc333")

    # --- year_range_from_env ----------------------------------------------
    saved = {k: os.environ.get(k)
             for k in ("MICASA_YEAR_START", "MICASA_YEAR_END", "MICASA_YEAR")}

    def _setenv(**kw):
        for k, v in kw.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    _setenv(MICASA_YEAR_START="2010", MICASA_YEAR_END="2012", MICASA_YEAR=None)
    check("year_range from START/END is the inclusive range",
          check_hashes.year_range_from_env() == [2010, 2011, 2012])
    _setenv(MICASA_YEAR_START=None, MICASA_YEAR_END=None, MICASA_YEAR="2020")
    check("year_range from a single MICASA_YEAR",
          check_hashes.year_range_from_env() == [2020])
    _setenv(**saved)

    # --- verify_dir: matching / mismatching / no-manifest .nc4 -------------
    vd = os.path.join(tmp, "vdir")
    os.makedirs(vd)
    fa = os.path.join(vd, "MiCASA_v1_flux_x3600_y1800_daily_20200101.nc4")
    fb = os.path.join(vd, "MiCASA_v1_flux_x3600_y1800_daily_20200102.nc4")
    fc = os.path.join(vd, "MiCASA_v1_flux_x3600_y1800_daily_20200103.nc4")
    for f, content in ((fa, b"AAA"), (fb, b"BBB"), (fc, b"CCC")):
        with open(f, "wb") as fh:
            fh.write(content)
    # manifest lists fa with its CORRECT hash, fb with a WRONG hash; fc absent
    with open(os.path.join(vd, "MiCASA_v1_x_sha256.txt"), "w") as fh:
        fh.write(f"{check_hashes.sha256_file(fa)}  {os.path.basename(fa)}\n")
        fh.write(f"0000000000000000  {os.path.basename(fb)}\n")
    v, nr, m, mm = check_hashes.verify_dir(vd, "MiCASA_*.nc4")
    check("verify_dir: the matching file is verified", v == 1)
    check("verify_dir: the file absent from the manifest is no_record", nr == 1)
    check("verify_dir: the wrong-hash file is a mismatch", m == 1)
    check("verify_dir: the mismatch detail names the bad file",
          len(mm) == 1 and "20200102" in mm[0])

if _failures:
    print(f"\n{len(_failures)} FAILED: {', '.join(_failures)}")
    sys.exit(1)
print("\nall check_hashes tests passed")
