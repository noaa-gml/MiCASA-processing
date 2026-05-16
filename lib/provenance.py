"""lib/provenance.py -- pipeline provenance for netCDF global attributes.

Builds the CF/ACDD-style global-attribute dict that every netCDF the MiCASA
pipeline writes should carry: producing software (git commit + version),
timestamp, host, input files (path + SHA-256), and citation metadata.
compute_clim.py calls provenance_attrs() and sets it on the output Dataset.

Citation constants (institution, DOI, ...) come from lib/provenance.conf --
a KEY="VALUE" file shared with lib/provenance.r and the shell stampers, so
the DOI lives in exactly one place.

Standard library only (no numpy / xarray), so provenance_attrs() is
unit-tested standalone -- tests/test_provenance.py.
"""
import datetime
import hashlib
import os
import socket
import subprocess


def load_conf(work_dir):
    """Parse lib/provenance.conf (KEY="VALUE" lines, '#' comments) -> dict.

    Missing file -> empty dict (callers fall back to defaults).
    """
    path = os.path.join(work_dir, "lib", "provenance.conf")
    conf = {}
    if not os.path.exists(path):
        return conf
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, val = line.split("=", 1)
            key, val = key.strip(), val.strip()
            if len(val) >= 2 and val[0] == val[-1] and val[0] in "\"'":
                val = val[1:-1]
            conf[key] = val
    return conf


def _git(repo_dir, *args, default="unknown"):
    """First stdout line of `git -C repo_dir <args>`, or `default`."""
    try:
        out = subprocess.run(["git", "-C", repo_dir, *args],
                             capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.SubprocessError):
        return default
    lines = out.stdout.strip().splitlines()
    return lines[0] if out.returncode == 0 and lines else default


def git_commit(repo_dir):
    return _git(repo_dir, "rev-parse", "HEAD")


def git_version(repo_dir):
    return _git(repo_dir, "describe", "--tags", "--always", "--dirty")


def file_sha256(path):
    """SHA-256 hex digest of a file, or None if it cannot be read."""
    if not path or not os.path.isfile(path):
        return None
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


def timestamp():
    """ISO-8601 UTC timestamp, e.g. '2026-05-16T05:21:09Z'."""
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")


def provenance_attrs(step, work_dir, title=None, summary=None,
                     inputs=None, extra=None):
    """Ordered dict of netCDF global attributes for a file written by `step`.

    step     producing script, e.g. "compute_clim.py"
    work_dir pipeline checkout dir (git repo root; holds lib/)
    title    ACDD title    (optional)
    summary  ACDD summary  (optional)
    inputs   mapping name -> path; each emits input_<name> and
             input_<name>_sha256
    extra    mapping merged last (wins on key collision)
    """
    conf = load_conf(work_dir)
    doi      = conf.get("MICASA_DOI", "PENDING")
    landing  = conf.get("MICASA_LANDING_PAGE", "PENDING")
    pipe     = conf.get("MICASA_PROV_PIPELINE", "MiCASA-processing")
    pipe_url = conf.get("MICASA_PROV_PIPELINE_URL", "")
    commit   = git_commit(work_dir)
    version  = git_version(work_dir)
    has_doi     = bool(doi) and doi != "PENDING"
    has_landing = bool(landing) and landing != "PENDING"

    refs = ["%s pipeline: %s" % (pipe, pipe_url)]
    if has_doi:
        refs.append("dataset DOI: https://doi.org/%s" % doi)
    if has_landing:
        refs.append("dataset landing page: %s" % landing)

    ts = timestamp()
    a = {}
    a["Conventions"] = conf.get("MICASA_PROV_CONVENTIONS", "CF-1.10, ACDD-1.3")
    if title:
        a["title"] = title
    if summary:
        a["summary"] = summary
    a["institution"] = conf.get("MICASA_PROV_INSTITUTION", "")
    a["source"] = "%s pipeline, step %s" % (pipe, step)
    a["references"] = " ; ".join(refs)
    a["license"] = conf.get("MICASA_PROV_LICENSE", "")
    a["creator_name"] = conf.get("MICASA_PROV_CREATOR_NAME", "")
    a["creator_url"] = conf.get("MICASA_PROV_CREATOR_URL", "")
    a["date_created"] = ts
    a["processing_pipeline"] = pipe
    a["processing_pipeline_url"] = pipe_url
    a["processing_pipeline_commit"] = commit
    a["processing_pipeline_version"] = version
    a["processing_step"] = step
    a["processing_host"] = socket.gethostname()
    if has_doi:
        a["doi"] = doi

    if inputs:
        for name, path in inputs.items():
            a["input_%s" % name] = str(path)
            digest = file_sha256(path)
            a["input_%s_sha256" % name] = digest or "unavailable"

    # history: one CF audit line; downstream NCO tools append their own.
    a["history"] = "%s: created by %s [%s %s, commit %s]" % (
        ts, step, pipe, version, commit[:12])

    if extra:
        for k, v in extra.items():
            a[k] = v if isinstance(v, str) else str(v)

    # drop empty values so we never write blank attributes
    return {k: v for k, v in a.items() if v not in (None, "")}
