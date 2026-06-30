#!/usr/bin/env bash
set -euo pipefail

build_dir="${1:-build}"
log_file="$(mktemp)"
json_file="$(mktemp)"
bad_json_file="$(mktemp)"
bad_hex_json_file="$(mktemp)"
trap 'rm -f "$log_file" "$json_file" "$bad_json_file" "$bad_hex_json_file"' EXIT

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

run_fail_case() {
    local name="$1"
    shift

    if "$build_dir/mr-lrc-generator" "$@" >"$log_file" 2>&1; then
        printf 'not ok %s\n' "$name" >&2
        sed -n '1,120p' "$log_file" >&2
        return 1
    fi

    printf 'ok %s\n' "$name"
}

expect_log_text() {
    local needle="$1"

    if grep -Fq "$needle" "$log_file"; then
        return
    fi

    printf 'not ok log_contains_%s\n' "$needle" >&2
    sed -n '1,120p' "$log_file" >&2
    return 1
}

run_case "small_cauchy" -k 4 -g 2 -r 1 -p 1 -s 7 --construction 0 -m cauchy
run_case "seed_optional" -k 4 -g 2 -r 1 -p 1 --construction 0 -m cauchy --random-limit 1
run_case "json_output" -k 4 -g 2 -r 1 -p 1 -s 7 --construction 0 -m cauchy --json "$json_file"
run_case "check_json_output" --check-json "$json_file"

expect_json_text() {
    local needle="$1"

    if grep -Fq "$needle" "$json_file"; then
        return
    fi

    printf 'not ok json_output_contains_%s\n' "$needle" >&2
    sed -n '1,160p' "$json_file" >&2
    return 1
}

expect_json_text '"candidate_source":'
expect_json_text '"gf256_backend":'
expect_json_text '"groups": ['
expect_json_text '    { "data": [0, 1], "local": [4] },'
expect_json_text '"matrix_hex": ['

cat >"$bad_json_file" <<'JSON'
{
  "data_cnt": 4,
  "group_cnt": 2,
  "local_parity_cnt": 1,
  "global_parity_cnt": 1,
  "matrix": [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]
  ]
}
JSON
run_fail_case "check_json_rejects_non_mr" --check-json "$bad_json_file"

cat >"$bad_hex_json_file" <<'JSON'
{
  "data_cnt": 4,
  "group_cnt": 2,
  "local_parity_cnt": 1,
  "global_parity_cnt": 1,
  "matrix_hex": [
    ["01", "00", "00", "00"],
    ["00", "01", "00", "00"],
    ["00", "00", "01", "00"],
    ["00", "00", "00", "01"],
    ["00", "00", "00", "00"],
    ["00", "00", "00", "00"],
    ["00", "00", "00", "00"]
  ]
}
JSON
run_fail_case "check_json_hex_rejects_non_mr" --check-json "$bad_hex_json_file"

run_case "cauchy_dedup" -k 4 -g 2 -r 1 -p 1 -s 7 --construction 0 -m cauchy --cauchy-dedup
run_case "three_group_cauchy" -k 6 -g 3 -r 1 -p 2 -s 123 --construction 0 --local-method cauchy --global-method cauchy
run_case "global_column_multiplier_cauchy" -k 6 -g 2 -r 1 -p 2 -s 17 --construction 0 --local-method cauchy --global-method column_multiplier_cauchy
run_fail_case "local_column_multiplier_cauchy_rejected" -k 4 -g 2 -r 1 -p 1 -s 7 --construction 0 --local-method column_multiplier_cauchy --global-method cauchy
run_case "small_vandermonde" -k 4 -g 2 -r 1 -p 1 -s 7 --construction 0 --local-method vandermonde --global-method vandermonde
run_case "local_a2_random_global" -k 6 -g 2 -r 2 -p 1 -s 99 --construction 0 --local-method cauchy --global-method random
run_case "random_random_limit" -k 6 -g 2 -r 1 -p 2 -s 1 --construction 0 --local-method random --global-method random --random-limit 10
run_case "random_stress_prefilter" -k 6 -g 2 -r 1 -p 2 -s 1 --construction 0 -m cauchy --prefilter-count 8
run_fail_case "targeted_prefilter_10_2_4_rejects" -k 10 -g 2 -r 1 -p 4 -s 1 --construction 0 --local-method cauchy --global-method column_multiplier_cauchy --random-limit 1 --prefilter-count 1
run_case "vandermonde_random_threads" -k 6 -g 2 -r 1 -p 2 -s 1 --construction 0 --local-method vandermonde --global-method random --thread-count 2
run_case "random_cauchy_threads" -k 6 -g 2 -r 1 -p 2 -s 1 --construction 0 --local-method random --global-method cauchy --random-limit 10 --thread-count 4
run_case "construction_default_off_6_2_2" -k 6 -g 2 -r 1 -p 2 -s 1 --random-limit 1
run_case "construction_disabled_6_2_2" -k 6 -g 2 -r 1 -p 2 -s 1 --construction 0 -m cauchy --random-limit 1
run_case "difference_pack_h2_subspace_28_2_2" -k 28 -g 2 -r 1 -p 2 -s 1 --construction 1 --random-limit 1 --step-time 0
expect_log_text "candidate_source=difference_pack_h2"
expect_log_text "attempt=1"
run_case "difference_pack_h2_spread_32_4_2" -k 32 -g 4 -r 1 -p 2 -s 1 --construction 1 --random-limit 1 --step-time 0
expect_log_text "candidate_source=difference_pack_h2"
expect_log_text "attempt=1"
run_case "difference_pack_h2_spread_48_4_2" -k 48 -g 4 -r 1 -p 2 -s 1 --construction 1 --random-limit 1 --step-time 0
expect_log_text "candidate_source=difference_pack_h2"
expect_log_text "attempt=1"
run_case "difference_pack_h2_spread_64_8_2" -k 64 -g 8 -r 1 -p 2 -s 1 --construction 1 --random-limit 1 --step-time 0
expect_log_text "candidate_source=difference_pack_h2"
expect_log_text "attempt=1"
run_case "feasible_points_h3_12_4_3" -k 12 -g 4 -r 1 -p 3 -s 100 --construction 1 --random-limit 0 --step-time 0 --local-method random --global-method random
expect_log_text "candidate_source=feasible_points_h3"
expect_log_text "attempt=1"

# Industrial data-local LRC(k,l,h) style candidates. These use one local
# parity per local group, so l maps to --groups and h maps to --global-parity.
run_case "azure_lrc_6_2_2" -k 6 -g 2 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "azure_lrc_8_2_2" -k 8 -g 2 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "azure_lrc_10_2_2" -k 10 -g 2 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "azure_lrc_12_2_2" -k 12 -g 2 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "azure_lrc_14_2_2" -k 14 -g 2 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "azure_lrc_12_4_2" -k 12 -g 4 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "azure_lrc_8_2_3" -k 8 -g 2 -r 1 -p 3 -s 1 --construction 0 -m cauchy
run_case "azure_lrc_10_2_3" -k 10 -g 2 -r 1 -p 3 -s 1 --construction 0 -m cauchy
run_case "storage_spaces_8_2_1" -k 8 -g 2 -r 1 -p 1 -s 1 --construction 0 -m cauchy
run_case "storage_spaces_12_2_1" -k 12 -g 2 -r 1 -p 1 -s 1 --construction 0 -m cauchy

# Representative all-symbol MR-LRC (g,r,a,h) scales, mapped to this generator as
# data = g * (r - a) - h, groups = g, local parity = a, global parity = h.
run_case "mr_3_5_1_2" -k 10 -g 3 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "mr_4_4_1_2" -k 10 -g 4 -r 1 -p 2 -s 1 --construction 0 -m cauchy
run_case "mr_2_7_1_3" -k 9 -g 2 -r 1 -p 3 -s 1 --construction 0 -m cauchy
run_case "mr_3_5_1_3" -k 9 -g 3 -r 1 -p 3 -s 1 --construction 0 -m cauchy
run_case "mr_2_6_2_2" -k 6 -g 2 -r 2 -p 2 -s 1 --construction 0 -m cauchy
run_case "mr_3_6_2_2" -k 10 -g 3 -r 2 -p 2 -s 1 --construction 0 -m cauchy
