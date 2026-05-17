# lib/manifest.sh -- append structured run records to jobs/run_manifest.tsv.
#
# Source this file, then call:
#     manifest_record <step> <status> <elapsed_s> <detail>
#
#   step       producing script / stage, e.g. "diurnalize-ERA5.r"
#   status     start | ok | fail | info
#   elapsed_s  integer seconds, or "-" when not applicable
#   detail     free text (tabs / newlines are squashed to spaces)
#
# The manifest is the pipeline's structured run record -- verify_v2 reads it
# instead of globbing job logs. The file is tab-separated; the columns are
#     timestamp  step  status  host  commit  elapsed_s  detail
#
# manifest_record never fails its caller: a logging call must not abort the
# pipeline, so the whole body runs in a subshell guarded by `|| true` and is
# safe to call under `set -e`.

manifest_record() {
    (
        _mf_step=${1:-unknown}
        _mf_status=${2:-info}
        _mf_elapsed=${3:--}
        _mf_detail=${4:-}
        _mf_root=${WORK_DIR:-$(pwd)}
        _mf_dir=$_mf_root/jobs
        _mf_file=$_mf_dir/run_manifest.tsv
        mkdir -p "$_mf_dir" 2>/dev/null || exit 0
        if [ ! -f "$_mf_file" ]; then
            printf '# timestamp\tstep\tstatus\thost\tcommit\telapsed_s\tdetail\n' \
                > "$_mf_file" 2>/dev/null || exit 0
        fi
        _mf_ts=$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo unknown)
        _mf_host=$(hostname 2>/dev/null || echo unknown)
        _mf_commit=$(git -C "$_mf_root" rev-parse --short HEAD 2>/dev/null || echo unknown)
        _mf_step=$(printf '%s' "$_mf_step"   | tr '\t\n' '  ')
        _mf_detail=$(printf '%s' "$_mf_detail" | tr '\t\n' '  ')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$_mf_ts" "$_mf_step" "$_mf_status" "$_mf_host" "$_mf_commit" \
            "$_mf_elapsed" "$_mf_detail" >> "$_mf_file" 2>/dev/null || exit 0
    ) || true
    return 0
}
