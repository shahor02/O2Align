#!/usr/bin/env bash
# Merge equally named ROOT files from consecutive directories of each run.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: merge_runs.sh [-f] DIRECTORY [DIRECTORY ...]

Merge all root files supplied */LHC*/<RUN>/*/ directories into */LHC*/<RUN>/merged/.
Directories belonging to each run must be consecutive; their order is preserved.
Only files directly inside the supplied directories are included.
Use -f to overwrite existing outputs. Run in an O2 environment.
Set HADDV to override the merger executable (default: $O2Align/shared/macro/haddv.py).
Note that haddv.py script works as the roor hadd utility but also concatenates the std::vectors.
EOF
}

fail() { printf 'merge_runs: %s\n' "$*" >&2; exit 1; }

options=()
while (($#)); do
    case $1 in
        -f) options=(-f); shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; break ;;
        -*) fail "Unknown option: $1" ;;
        *) break ;;
    esac
done
(($#)) || { usage >&2; exit 1; }

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
merger=${HADDV:-}
if [[ -z $merger ]]; then
    if command -v haddv >/dev/null 2>&1; then
        merger=haddv
    elif command -v haddv.py >/dev/null 2>&1; then
        merger=haddv.py
    else
        merger=$script_dir/haddv.py
    fi
fi
command -v "$merger" >/dev/null 2>&1 || fail "Merger executable not found: $merger"

work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
shopt -s nullglob
declare -A lists=()
names=()
run_dir=

merge_run() {
    [[ -n $run_dir ]] || return 0
    ((${#names[@]})) || fail "No ROOT files found for $run_dir"
    mkdir -p -- "$run_dir/merged"
    local name
    for name in "${names[@]}"; do
        printf 'Merging %s -> %s/merged/%s\n' "${run_dir##*/}" "$run_dir" "$name"
        "$merger" "${options[@]}" "$run_dir/merged/$name" "@${lists[$name]}"
    done
}

for directory in "$@"; do
    [[ -d $directory ]] || fail "Not a directory: $directory"
    # Make paths absolute without resolving symlinks in the LHC/run prefix.
    [[ $directory = /* ]] || directory=$PWD/$directory
    directory=${directory%/}
    if [[ $directory =~ ^(.*?/LHC[^/]+/([0-9]+))(/|$) ]]; then
        next_run=${BASH_REMATCH[1]}
    else
        fail "Cannot extract LHC*/RUN from: $directory"
    fi
    if [[ $next_run != "$run_dir" ]]; then
        merge_run
        run_dir=$next_run
        lists=()
        names=()
    fi
    for file in "$directory"/*.root; do
        [[ -f $file ]] || continue
        name=${file##*/}
        if [[ ! ${lists[$name]+present} ]]; then
            lists[$name]=$(mktemp "$work/inputs.XXXXXXXX")
            names+=("$name")
        fi
        printf '%s\n' "$file" >> "${lists[$name]}"
    done
done
merge_run
