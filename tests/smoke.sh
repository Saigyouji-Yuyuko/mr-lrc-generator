#!/usr/bin/env bash
set -euo pipefail

build_dir="${1:-build}"
log_file="$(mktemp)"
trap 'rm -f "$log_file"' EXIT

cmake -S . -B "$build_dir"
cmake --build "$build_dir"

run_case() {
    local name="$1"
    shift

    if "$build_dir/mr-lrc-generator" "$@" >"$log_file" 2>&1; then
        printf 'ok %s\n' "$name"
        return
    fi

    printf 'not ok %s\n' "$name" >&2
    sed -n '1,120p' "$log_file" >&2
    return 1
}

run_case "small_cauchy" -k 4 -g 2 -r 1 -p 1 -s 7 --construction false -m cauchy
run_case "three_group_cauchy" -k 6 -g 3 -r 1 -p 2 -s 123 --construction false --local-method cauchy --global-method cauchy
run_case "small_vandermonde" -k 4 -g 2 -r 1 -p 1 -s 7 --construction false --local-method vandermonde --global-method vandermonde
run_case "local_a2_random_global" -k 6 -g 2 -r 2 -p 1 -s 99 --construction false --local-method cauchy --global-method random
run_case "random_random_limit" -k 6 -g 2 -r 1 -p 2 -s 1 --construction false --local-method random --global-method random --random-limit 10
run_case "vandermonde_random_threads" -k 6 -g 2 -r 1 -p 2 -s 1 --construction false --local-method vandermonde --global-method random --thread-count 2
run_case "random_cauchy_threads" -k 6 -g 2 -r 1 -p 2 -s 1 --construction false --local-method random --global-method cauchy --random-limit 10 --thread-count 4
run_case "construction_default_6_2_2" -k 6 -g 2 -r 1 -p 2 -s 1 --random-limit 1
run_case "construction_disabled_6_2_2" -k 6 -g 2 -r 1 -p 2 -s 1 --construction false -m cauchy --random-limit 1

# Industrial data-local LRC(k,l,h) style candidates. These use one local
# parity per local group, so l maps to --groups and h maps to --global-parity.
run_case "azure_lrc_6_2_2" -k 6 -g 2 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "azure_lrc_8_2_2" -k 8 -g 2 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "azure_lrc_10_2_2" -k 10 -g 2 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "azure_lrc_12_2_2" -k 12 -g 2 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "azure_lrc_14_2_2" -k 14 -g 2 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "azure_lrc_12_4_2" -k 12 -g 4 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "azure_lrc_8_2_3" -k 8 -g 2 -r 1 -p 3 -s 1 --construction false -m cauchy
run_case "azure_lrc_10_2_3" -k 10 -g 2 -r 1 -p 3 -s 1 --construction false -m cauchy
run_case "storage_spaces_8_2_1" -k 8 -g 2 -r 1 -p 1 -s 1 --construction false -m cauchy
run_case "storage_spaces_12_2_1" -k 12 -g 2 -r 1 -p 1 -s 1 --construction false -m cauchy

# Representative all-symbol MR-LRC (g,r,a,h) scales, mapped to this generator as
# data = g * (r - a) - h, groups = g, local parity = a, global parity = h.
run_case "mr_3_5_1_2" -k 10 -g 3 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "mr_4_4_1_2" -k 10 -g 4 -r 1 -p 2 -s 1 --construction false -m cauchy
run_case "mr_2_7_1_3" -k 9 -g 2 -r 1 -p 3 -s 1 --construction false -m cauchy
run_case "mr_3_5_1_3" -k 9 -g 3 -r 1 -p 3 -s 1 --construction false -m cauchy
run_case "mr_2_6_2_2" -k 6 -g 2 -r 2 -p 2 -s 1 --construction false -m cauchy
run_case "mr_3_6_2_2" -k 10 -g 3 -r 2 -p 2 -s 1 --construction false -m cauchy
