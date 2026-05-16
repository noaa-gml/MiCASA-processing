#!/usr/bin/env python3
"""stamp_provenance.py -- write CF/ACDD provenance global attributes onto
existing netCDF files, using the helper in lib/provenance.py.

Two uses:

  * In-pipeline -- stamp a file the current run just produced. cat_monthly.sh
    calls this on the concatenated monthly file. The git commit recorded is
    the commit that produced the file, which is correct since it was just
    made:

        stamp_provenance.py FILE --step cat_monthly.sh --title "..." \\
            [--summary "..."]

  * Retrofit -- add provenance to files generated *before* the pipeline
    stamped its own outputs. Such files cannot recover their true generating
    commit or input checksums, so --retrofit writes only the static
    citation/pipeline attributes plus an explicit `provenance_note` saying
    so; it never asserts a generating commit it does not know:

        stamp_provenance.py --retrofit FILE [FILE ...]

Requires the netCDF4 package (e.g. the `sci` conda env on Orion). Idempotent:
re-stamping a file overwrites the provenance attributes and appends one line
to the CF `history` attribute.
"""
import argparse
import os
import sys

# lib/provenance.py lives next to this script, under lib/.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "lib"))
from provenance import provenance_attrs, git_commit, git_version, timestamp

# In --retrofit mode only these static citation/pipeline keys are asserted --
# everything run-specific (commit, timestamp, inputs, step) is unknowable for
# a file whose true provenance was never recorded. `source` is set
# separately so it does not carry a misleading processing step.
_RETROFIT_STATIC_KEYS = (
    "Conventions", "institution", "references", "license",
    "creator_name", "creator_url", "doi",
    "processing_pipeline", "processing_pipeline_url",
)


def _append_history(ds, line):
    """Append one line to the CF `history` attribute (create if absent)."""
    existing = getattr(ds, "history", "")
    ds.history = (existing + "\n" + line) if existing else line


def stamp(path, step, work_dir, title=None, summary=None, retrofit=False):
    """Write provenance global attributes onto one netCDF file in place."""
    import netCDF4
    ds = netCDF4.Dataset(path, "r+")
    try:
        attrs = provenance_attrs(step=step, work_dir=work_dir,
                                 title=title, summary=summary)
        if retrofit:
            for k in _RETROFIT_STATIC_KEYS:
                if attrs.get(k):
                    ds.setncattr(k, attrs[k])
            ds.setncattr("source", "%s pipeline" % attrs["processing_pipeline"])
            if title:
                ds.setncattr("title", title)
            ds.setncattr(
                "provenance_note",
                "Static citation/pipeline attributes added retroactively by "
                "stamp_provenance.py; this file predates in-pipeline "
                "provenance, so its generating git commit and input "
                "checksums are not recorded. Regenerate via MiCASA-processing "
                "for complete provenance.")
            ds.setncattr("provenance_retrofit_commit", git_commit(work_dir))
            ds.setncattr("provenance_retrofit_date", timestamp())
            _append_history(
                ds, "%s: provenance metadata added retroactively by "
                    "stamp_provenance.py [MiCASA-processing %s]"
                    % (timestamp(), git_version(work_dir)))
        else:
            hist = attrs.pop("history", None)
            for k, v in attrs.items():
                ds.setncattr(k, v)
            if hist:
                _append_history(ds, hist)
    finally:
        ds.close()


def main(argv=None):
    p = argparse.ArgumentParser(
        description="stamp CF/ACDD provenance global attributes onto netCDF files")
    p.add_argument("files", nargs="+", help="netCDF file(s) to stamp")
    p.add_argument("--step", default="stamp_provenance.py",
                   help="producing pipeline step (recorded as processing_step)")
    p.add_argument("--title", default=None, help="ACDD title attribute")
    p.add_argument("--summary", default=None, help="ACDD summary attribute")
    p.add_argument("--retrofit", action="store_true",
                   help="add only static attributes (file predates "
                        "in-pipeline provenance)")
    p.add_argument("--work-dir", default=os.environ.get("WORK_DIR", _HERE),
                   help="pipeline checkout dir (git repo root; holds lib/)")
    args = p.parse_args(argv)

    fail = 0
    for path in args.files:
        if not os.path.isfile(path):
            print("stamp_provenance: not a file: %s" % path, file=sys.stderr)
            fail += 1
            continue
        try:
            stamp(path, step=args.step, work_dir=args.work_dir,
                  title=args.title, summary=args.summary,
                  retrofit=args.retrofit)
            print("stamped %s%s" % ("(retrofit) " if args.retrofit else "", path))
        except Exception as exc:                            # noqa: BLE001
            print("stamp_provenance: FAILED %s: %s" % (path, exc),
                  file=sys.stderr)
            fail += 1
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
