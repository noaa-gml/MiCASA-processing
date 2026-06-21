# Runbook — producing a PIQS-fitted product update (Path 1)

**When to use this.** You need a MiCASA product update that uses the **legacy PIQS**
sub-monthly fitter, but produced on the **current V2 pipeline** (with its
robustness fixes — corrected 0.1°→1° aggregation, FastTrack meteo, per-month
climatology auto-detect, provenance, manifest). This is *not* a bit-faithful
reproduction of the archived V1 product — for that, use the `legacy` branch
(PIQS-only). It *is* a PIQS-shaped product on the maintained pipeline.

## One-line setup

The fitter is a flag on `run_year.sh` (no hand-edit):

```sh
./run_year.sh <YEAR> <VERSION> --fitter piqs
```

- `<VERSION>` = `v1` (final/archived line) or `vNRT` (near-real-time).
- `--fitter piqs` swaps stage 4 from `write_pchip.r` to `write_piqs.r` and
  auto-sets `MICASA_PIQS_PAD_RIGHT=2` (the NRT trailing-edge padding, proposal #1).
  Stages 1–3, 5, 6 are the unchanged V2 pipeline.
- `--fitter` also accepts `pchip` (default), `ppm`, `linmm`, `mss`, `atpk`; all
  write `fit.piqs.rda` (recording the fitter in `piqsfit.meta$fitter`), which the
  diurnalize stage reads transparently.

Prerequisites (env, per `config.sh`): `MAIL_USER`, `BASE_DIR`. Runs on Orion; the
ingest and diurnalize stages submit SBATCH jobs with `--wait`.

## The one thing to design around: PIQS is non-local

PIQS is a **single global solve over the whole record**, so adding or revising
months **re-fits everything and rewrites every historical month** (`verify_v2`
NRT footprint: all 303; this is the exact property the V2 PCHIP default was chosen
to avoid — see [`V1_TO_V2_JUSTIFICATION.md`](V1_TO_V2_JUSTIFICATION.md) §1).
Practical consequences:

1. **Stage 4 refits the full record**, not just `<YEAR>` — expect it to read the
   whole concatenated monthly file.
2. **Re-diurnalize the affected tail**, not only the new month: the global refit
   shifts prior coefficients across the record (`PAD_RIGHT=2` mutes but does not
   remove the trailing-edge effect). Safest: re-diurnalize the whole window you
   republish, or at minimum the last ~12 months + the new months.
3. **The published past will shift slightly each cycle** — budget for that
   downstream. (PCHIP would not; that is the V2 trade-off.)

## Update procedure

1. **Full run for a year** (download → ingest → aggregate → PIQS fit → diurnalize
   → daysplit):
   ```sh
   ./run_year.sh 2026 v1 --fitter piqs
   ```
2. **Incremental** (data already on disk) — go straight to fit + diurnalize +
   daysplit:
   ```sh
   ./run_year.sh 2026 v1 --fitter piqs --skip-download --skip-ingest --skip-aggregate
   ```
   Because the PIQS fit is global, then **re-diurnalize the changed tail** (rerun
   stage 5 for the affected months, or invoke `diurnalize-ERA5.r` with the
   relevant `MICASA_YEAR_START/END`).
3. **Verify**: run `verify_v2`. The PIQS-relevant checks are 2.1 (integral
   preservation), 3.1 (sub-monthly sign-flips — expect the higher PIQS rates, not
   PCHIP's), 11.2 (tail-coefficient stability), and 6.1 (v2-vs-v1 spatial sanity).

## What you get vs the archived V1

- ✓ PIQS sub-monthly shapes (the V1 fitter), recorded `piqsfit.meta$fitter == "piqs"`
  and in the output provenance attributes.
- ✗ **Not bit-identical to the original archived V1**: this runs on the V2 pipeline,
  which fixed the 0.1°→1° aggregation latitude-weight bug and added FastTrack
  meteo + per-month climatology. For a bug-for-bug V1 reproduction, check out the
  `legacy` branch instead.

## Dry run first

```sh
./run_year.sh 2026 v1 --fitter piqs --dry-run
```
prints the banner (`FITTER piqs (stage 4: write_piqs.r)`) and every stage command
without executing — use it to confirm the wiring before committing compute.
