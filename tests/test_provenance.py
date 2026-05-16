#!/usr/bin/env python3
"""Unit tests for lib/provenance.py (standard library only, CI-runnable).

Run:  python3 tests/test_provenance.py
Exits non-zero on any failure.
"""
import os
import sys
import tempfile

# Import provenance from lib/ one level up.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, "lib"))
import provenance

_failures = []


def check(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        _failures.append(name)


with tempfile.TemporaryDirectory() as tmp:
    # --- load_conf ---------------------------------------------------------
    os.makedirs(os.path.join(tmp, "lib"))
    with open(os.path.join(tmp, "lib", "provenance.conf"), "w") as fh:
        fh.write('# a comment\n\n'
                 'MICASA_DOI="10.0/ABCD-1234"\n'
                 'MICASA_PROV_INSTITUTION="Test Institution"\n'
                 'KEY_WITH_EQ="https://x/?id=9"\n'
                 "SINGLE='quoted'\n")
    conf = provenance.load_conf(tmp)
    check("load_conf reads KEY=VALUE", conf.get("MICASA_DOI") == "10.0/ABCD-1234")
    check("load_conf strips double quotes",
          conf.get("MICASA_PROV_INSTITUTION") == "Test Institution")
    check("load_conf keeps '=' inside values",
          conf.get("KEY_WITH_EQ") == "https://x/?id=9")
    check("load_conf strips single quotes", conf.get("SINGLE") == "quoted")
    check("load_conf skips comments and blank lines (4 keys)", len(conf) == 4)
    check("load_conf on a missing file -> empty dict",
          provenance.load_conf(os.path.join(tmp, "nope")) == {})

    # --- file_sha256 -------------------------------------------------------
    f_abc = os.path.join(tmp, "abc.txt")
    with open(f_abc, "wb") as fh:
        fh.write(b"abc")
    check("sha256 of 'abc' matches the known digest",
          provenance.file_sha256(f_abc) ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
    f_empty = os.path.join(tmp, "empty.txt")
    open(f_empty, "wb").close()
    check("sha256 of an empty file matches the known digest",
          provenance.file_sha256(f_empty) ==
          "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
    check("sha256 of a missing file -> None",
          provenance.file_sha256(os.path.join(tmp, "nope")) is None)

    # --- provenance_attrs --------------------------------------------------
    a = provenance.provenance_attrs("test-step.py", _REPO)
    need = ["Conventions", "institution", "source", "references", "license",
            "date_created", "processing_pipeline", "processing_pipeline_commit",
            "processing_pipeline_version", "processing_step",
            "processing_host", "history"]
    check("attrs has every required key", all(k in a for k in need))
    check("attrs processing_step is the step",
          a["processing_step"] == "test-step.py")
    check("attrs history line mentions the step", "test-step.py" in a["history"])
    check("attrs Conventions names CF and ACDD",
          "CF-" in a["Conventions"] and "ACDD" in a["Conventions"])
    check("attrs has no title when none is passed", "title" not in a)

    a2 = provenance.provenance_attrs("s.py", _REPO, title="My Title",
                                     summary="My summary")
    check("attrs includes title when passed", a2.get("title") == "My Title")
    check("attrs includes summary when passed", a2.get("summary") == "My summary")

    a3 = provenance.provenance_attrs("s.py", _REPO, inputs={"thing": f_abc})
    check("attrs inputs -> input_<name>", a3.get("input_thing") == f_abc)
    check("attrs inputs -> input_<name>_sha256",
          a3.get("input_thing_sha256") ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")

    a4 = provenance.provenance_attrs("s.py", _REPO,
                                     extra={"micasa_version": "v9",
                                            "processing_step": "overridden"})
    check("attrs extra adds new keys", a4.get("micasa_version") == "v9")
    check("attrs extra wins on key collision",
          a4.get("processing_step") == "overridden")

    check("references omit doi.org when the DOI is PENDING",
          "doi.org" not in a["references"])
    a5 = provenance.provenance_attrs("s.py", tmp)   # tmp conf carries a real DOI
    check("references include doi.org when the DOI is registered",
          "doi.org/10.0/ABCD-1234" in a5["references"])
    check("attrs expose a `doi` key when the DOI is registered",
          a5.get("doi") == "10.0/ABCD-1234")

    # --- git helpers -------------------------------------------------------
    check("git_commit returns a non-empty string",
          isinstance(provenance.git_commit(_REPO), str)
          and bool(provenance.git_commit(_REPO)))
    check("git_version returns a non-empty string",
          isinstance(provenance.git_version(_REPO), str)
          and bool(provenance.git_version(_REPO)))
    check("git on a non-repository -> 'unknown'",
          provenance.git_commit(tmp) == "unknown")

if _failures:
    print(f"\n{len(_failures)} FAILED: {', '.join(_failures)}")
    sys.exit(1)
print("\nall provenance tests passed")
