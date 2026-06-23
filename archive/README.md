# archive/

One-off and year-pinned scripts kept for reproducibility, not active use. They
are still syntax-checked by CI (the `bash -n` / `parse()` jobs `find` the whole
tree), but the maintained entry points live in the repository root.

| Archived script | Superseded by | Why kept |
|---|---|---|
| `produce_2025_2026.sh` | [`../run_record.sh`](../run_record.sh) | Year-pinned 2025–2026 NRT-tail builder. Kept for the NRT trailing-edge completion it special-cases (day-of-year clim-fill of the partial final month + a synthetic monthly file with a patched time coordinate), which `run_record.sh` omits. Run from the repo root; it `cd`s to `$(dirname "$0")/..` so its relative calls resolve. |
| `diag_v1_vNRT_handoff.r` | — | One-off splice-continuity diagnostic for the v1→vNRT (2024-12 / 2025-01) handoff. Paths are cwd-relative; run `Rscript archive/diag_v1_vNRT_handoff.r` from the working directory. |

For a clean multi-year (re)build of the whole record, use
`../run_record.sh YEAR_START YEAR_END [--v1-through YEAR]` instead.
