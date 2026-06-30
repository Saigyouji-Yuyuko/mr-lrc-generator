# MR-LRC Generator

[中文说明](README_zh.md)

MR-LRC Generator is a C++17 command-line tool for generating and exactly
checking GF(256) data-local MR-LRC generator matrices. It uses ISA-L for
GF(256) arithmetic and matrix inversion.

The tool builds a systematic generator matrix:

```text
[ I_k ]
[ local parity rows ]
[ global parity rows ]
```

Data symbols are split as evenly as possible across local groups. Local parity
rows only cover data symbols in their own group, while global parity rows cover
the full stripe. Candidate matrices are generated with independent local and
global construction methods, then the checker enumerates maximum erasure cases
and verifies recoverability over GF(256).

## Features

- Data-local systematic LRC layout.
- Independent local and global construction methods. Local rows support
  `cauchy`, `vandermonde`, or `random`; global rows also support
  `column_multiplier_cauchy`.
- Exact MR check over all maximum erasure cases.
- Optional `--seed` for deterministic random search; omitted seeds are
  generated and printed.
- Search bound via `--random-limit`.
- Parallel search via `--thread-count`.
- Dedicated `h=2` difference-pack construction for large data-local LRCs, plus
  an experimental `h=3` feasible-points construction for small GF(256) cases.
- Optional stress-test prefilter before exact verification.
- Exact verification of a matrix supplied in JSON.
- ISA-L-backed GF(256) arithmetic.

## Build

```bash
cmake -S . -B build -DMR_LRC_REQUIRE_ISAL=ON
cmake --build build
```

If ISA-L is installed in a non-default location:

```bash
cmake -S . -B build \
  -DMR_LRC_REQUIRE_ISAL=ON \
  -DISAL_INCLUDE_DIR=/usr/include/isa-l \
  -DISAL_LIBRARY=/usr/lib/x86_64-linux-gnu/libisal.so
```

## Usage

Generate and check a `6 data + 2 local parity + 2 global parity` code. This
example uses 2 local groups, 1 local parity row per group, and 2 global parity
rows:

```bash
./build/mr-lrc-generator \
  --data 6 \
  --groups 2 \
  --local-parity 1 \
  --global-parity 2 \
  --seed 1 \
  --construction 0 \
  --local-method cauchy \
  --global-method cauchy
```

You can choose different construction methods for local and global parity rows:

```bash
./build/mr-lrc-generator \
  --data 6 \
  --groups 2 \
  --local-parity 1 \
  --global-parity 2 \
  --seed 1 \
  --construction 0 \
  --local-method random \
  --global-method cauchy \
  --random-limit 100 \
  --thread-count 4
```

The shorter `--method` alias sets both local and global methods:

```bash
./build/mr-lrc-generator --data 6 --groups 2 --local-parity 1 \
  --global-parity 2 --seed 1 --method cauchy
```

Because `column_multiplier_cauchy` is global-only, select it with
`--global-method` rather than `--method`.

For large data-local `LRC(k,l,2)` parameters, use the dedicated h=2 construction
described in [H=2 Difference-Pack Construction](#h2-difference-pack-construction).

## Options

Required parameters:

| Option | Meaning |
| --- | --- |
| `-k`, `--data K` | Number of systematic data symbols. |
| `-g`, `--groups G` | Number of local groups. Data is split nearly evenly. |
| `-r`, `--local-parity R` | Local parity symbols per group. |
| `-p`, `--global-parity M` | Global parity symbols. |

Optional parameters:

| Option | Meaning |
| --- | --- |
| `-s`, `--seed S` | Random seed. If omitted, one is generated and printed. |
| `--local-method M` | Local parity construction: `cauchy`, `vandermonde`, or `random`. |
| `--global-method M` | Global parity construction: `cauchy`, `column_multiplier_cauchy`, `vandermonde`, or `random`. |
| `-m`, `--method M` | Alias that sets both local and global methods. |
| `--construction N` | Registered construction attempt budget before random search. Default: `0` disables construction. Current paths cover h=2 difference-pack and experimental h=3 feasible-points candidates. |
| `--random-limit N` | Maximum candidate attempts. Default: unbounded `uint64` max. |
| `--prefilter-count N`, `--random-prefilter N`, `--stress-prefilter N` | Stress-test patterns to run before exact verification for each candidate. Default: `0` disables the prefilter. |
| `-t`, `--thread-count N` | Parallel search worker count. Default: `1`, max: `256`. |
| `--step_time N`, `--step-time N` | Print timestamped searched and strict-complete counts to stderr every N seconds. Default: `30`; `0` disables it. |
| `--json FILE`, `--matrix-json FILE` | Write the found result as pretty JSON to `FILE`, including stdout metadata, local groups, and decimal/hex matrices. Standard output keeps the normal text report. |
| `--check-json FILE`, `--verify-json FILE`, `--input-json FILE` | Read a matrix JSON file and exactly verify whether it is MR-LRC. Search parameters are not required in this mode. |
| `--cauchy-dedup` | Skip duplicate all-Cauchy candidates using canonical Cauchy parameter keys. Default: disabled. |
| `-h`, `--help` | Print CLI help. |

## H=2 Difference-Pack Construction

`--construction N` tries up to `N` registered constructive candidates before
the generic random search path. For h=2 data-local LRCs, use
`--construction 10000`; deterministic cases usually finish on attempt 1, while
the larger value gives the randomized fallback useful room.

h=2 construction successes are accepted by the difference-pack construction
proof and do not run the exact erasure verifier in `generate()`. Their output
therefore reports `strict_complete=false` and `patterns_checked=0`. Use
`--check-json` on the emitted matrix when you want an explicit exact verifier
run for a small enough h=2 instance.

Value semantics:

| Value | Meaning |
| --- | --- |
| `0` | Disable construction. This is the default. |
| `N > 0` | Try up to `N` registered construction attempts before generic random search. |
| `true`, `on`, `yes` | Compatibility aliases for `1`. |
| `false`, `off`, `no` | Compatibility aliases for `0`. |

Basic requirements for `difference_pack_h2`:

| Requirement | Value |
| --- | --- |
| Field | GF(256) |
| `--global-parity` | Exactly `2` |
| `--local-parity` | Exactly `1` per group |
| Per-group data symbols | `group_data + 1 <= 256` |
| Matrix shape | Local row is all `1`; global rows are `(t, t^2)` |
| Success marker | Output contains `candidate_source=difference_pack_h2` |

Construction paths:

| Path | Applies when | Guarantee |
| --- | --- | --- |
| Complementary subspaces | `sum ceil(log2(group_data + 1)) <= 8` | Deterministic packing. |
| 4-dimensional spread | `groups <= 17` and `group_data <= 15` | Deterministic packing. If `r` means local-group size including the one local parity symbol, this is `r <= 16`. |
| Randomized difference packing | Basic requirements hold, but deterministic paths do not apply | Consumes the `--construction` attempt budget; no deterministic size guarantee. |

At the spread limit, the largest evenly split data-local case is:

```text
groups = 17
group_data = 15
data = 255
symbols = 255 + 17 local parity + 2 global parity = 274
data / symbols = 255 / 274 ~= 93.07%
```

Example:

```bash
./build/mr-lrc-generator \
  --data 64 \
  --groups 8 \
  --local-parity 1 \
  --global-parity 2 \
  --construction 10000 \
  --random-limit 1 \
  --step-time 0
```

For the `data=255` edge, the default `--global-method cauchy` is rejected by
the generic method validator because `data + global_parity > 256`; use
`--global-method random` with the h=2 construction:

```bash
./build/mr-lrc-generator \
  --data 255 \
  --groups 17 \
  --local-parity 1 \
  --global-parity 2 \
  --global-method random \
  --construction 10000 \
  --random-limit 1 \
  --step-time 0
```

The experimental all-symbol skew-polynomial prototype is isolated in
`src/all_symbol_skew_lrc.cpp` and is not part of the build.

## H=3 Feasible-Points Construction

For data-local `local_parity=1`, `global_parity=3`, `--construction N` also
enables an experimental GF(256)^3 feasible-points construction before generic
random search. It builds each local group as a point set with the local parity
at the origin, filters residual vector/plane compatibility across groups, and
then runs the exact MR checker. This is a heuristic candidate generator, not a
proof shortcut; the exact checker remains the authority.

Basic requirements for `feasible_points_h3`:

| Requirement | Value |
| --- | --- |
| Field | GF(256) |
| `--global-parity` | Exactly `3` |
| `--local-parity` | Exactly `1` per group |
| Matrix shape | Local row is all `1`; global rows are 3D point coordinates |
| Success marker | Output contains `candidate_source=feasible_points_h3` |

Observed status from the current experiments:

| Status | Parameters |
| --- | --- |
| Exact-verified by the construction/search experiments | `LRC(12,2,3)`, `LRC(14,2,3)`, `LRC(16,2,3)`, `LRC(12,4,3)`, `LRC(16,3,3)`, `LRC(16,4,3)` |
| Not found by the current construction/search experiments | `LRC(18,2,3)`, `LRC(18,3,3)`, `LRC(20,4,3)` |

Example:

```bash
./build/mr-lrc-generator \
  --data 16 \
  --groups 4 \
  --local-parity 1 \
  --global-parity 3 \
  --construction 100 \
  --random-limit 0 \
  --local-method random \
  --global-method random \
  --step-time 0
```

## H=4 Experimental Status

No h=4 construction path is registered in the generator yet.
`global_parity=4` should be treated as experimental. Small cases such as
`LRC(8,2,4)` and `LRC(10,2,4)` are useful baselines, but the current GF(256)^4
point-set experiments, h=3 lifting attempts, and GF(16)/Moore-style probes have
not found `LRC(12,2,4)`, `LRC(12,3,4)`, `LRC(12,4,4)`, or `LRC(14,2,4)`.

The observed bottleneck is cross-group residual-subspace packing, especially
fourth-group triple compatibility, rather than constructing valid individual
local groups. Treat any larger h=4 run as research work until a dedicated
CSP/hypergraph construction is added and exact-verified.

## Construction Methods

`cauchy` and `vandermonde` use the usual GF(256) matrix families as coefficient
generators. The implementation picks distinct nonzero field parameters so the
base construction parameters do not collide.

`column_multiplier_cauchy` is a global-only method. It builds global rows from a
Cauchy matrix, then multiplies every data column by an independently selected
nonzero GF(256) scalar.

With `--cauchy-dedup`, when both local and global methods are `cauchy`, the
search keeps a canonical Cauchy-parameter key for each candidate. The key
removes the common Cauchy xor translation and ignores permutations of global
parity rows and of local parity rows within each local group, so equivalent
candidates are not verified twice.

`random` fills selected parity rows with random GF(256) coefficients. Local
random rows use nonzero coefficients to preserve the local pivot invariant used
by residual elimination; global random rows may still include `00`
coefficients. The exact MR checker is still the source of truth for whether the
candidate is accepted.

Local and global methods are independent. For example, `--local-method random
--global-method column_multiplier_cauchy` creates nonzero random local parity
rows and column-scaled Cauchy global parity rows.

## Checking a JSON Matrix

Use `--check-json FILE` to verify an existing systematic data-local generator
matrix without running the random search:

```bash
./build/mr-lrc-generator --check-json matrix.json
```

The command exits with status `0` when the matrix is MR-LRC and `2` when it is
not. The JSON written by `--json` can be fed back directly. A minimal input can
also use the compatibility field names below:

```json
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
    [1, 1, 0, 0],
    [0, 0, 1, 1],
    [1, 2, 3, 4]
  ]
}
```

`matrix` entries are decimal GF(256) bytes in `0..255`. Alternatively,
`matrix_hex` may be used with hex strings such as `"0a"`. If a `groups` array
is present, its `data` layout is honored; otherwise data symbols are split
evenly across groups, matching the generator. The verifier requires identity
data rows and local parity rows that only touch data symbols in their own
group.

## Scale Notes

- For `local_parity=1` and `global_parity=2`, prefer `--construction 10000`.
  This path has deterministic coverage for complementary binary subspaces and
  for the 4-dimensional spread case with up to 17 groups and local-group size
  `r <= 16` (`group_data <= 15`), so wide data-local `LRC(k,l,2)` layouts such
  as `LRC(64,8,2)` are practical.
- For `local_parity=1` and `global_parity=3`, `--construction` enables the
  experimental feasible-points path. Current experiments support small GF(256)
  layouts up to the tested `k=16` cases, with the exact checker still required
  for every accepted matrix.
- The generic search/check path still becomes expensive as `global_parity`
  grows. For unconstructed or larger `global_parity >= 3` parameters, use small
  parameters first and consider `--prefilter-count`.
- No h=4 construction is currently registered. `global_parity >= 4` remains
  experimental and is not recommended for large layouts.

## Verification Model

The total parity count is computed as:

```text
total_parity = groups * local_parity + global_parity
```

A maximum erasure pattern has size `total_parity` and must erase at least
`local_parity` symbols in every local group, then spend the remaining erasure
budget globally. The checker enumerates these maximum erasure patterns and
accepts a candidate only when every enumerated case is recoverable.

`--prefilter-count` builds that many stress-probe recipes before the search
starts, then runs those recipes for each candidate before the exact checker. For
the common `local_parity=1`, `groups=2..6`, `global_parity=2..4` range, the
probe schedule prioritizes residual shapes such as line collisions,
vector-in-plane cases, plane intersections, hyperplane-vector cases, and small
multi-vector rank defects. If the targeted schedule is exhausted, remaining
budget falls back to residual random maximum-erasure samples. A candidate
rejected by the prefilter is genuinely non-MR because each witness is confirmed
with the same local-elimination residual rank criterion used by the exact
verifier; candidates that pass the prefilter still go through the exact checker
and full gate.

`--random-limit` bounds the total number of candidate matrices. `--thread-count`
only controls how many workers draw from that shared attempt budget; it does not
multiply the total number of attempts. With multiple threads, the first worker
that finds a passing candidate wins, so the reported attempt is the first
completed success, not necessarily the lowest passing attempt number.

## Parameter Examples

Industrial `LRC(k,l,h)` style parameters map directly to:

```text
--data k --groups l --local-parity 1 --global-parity h
```

Representative smoke-test cases:

| Name | Notation | CLI parameters |
| --- | --- | --- |
| `azure_lrc_6_2_2` | `LRC(6,2,2)` | `--data 6 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_8_2_2` | `LRC(8,2,2)` | `--data 8 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_12_4_2` | `LRC(12,4,2)` | `--data 12 --groups 4 --local-parity 1 --global-parity 2` |
| `h2_subspace_28_2_2` | `LRC(28,2,2)` | `--data 28 --groups 2 --local-parity 1 --global-parity 2 --construction 10000` |
| `h2_spread_64_8_2` | `LRC(64,8,2)` | `--data 64 --groups 8 --local-parity 1 --global-parity 2 --construction 10000` |
| `h3_points_12_4_3` | `LRC(12,4,3)` | `--data 12 --groups 4 --local-parity 1 --global-parity 3 --construction 1 --random-limit 0 --local-method random --global-method random` |
| `h3_points_16_4_3` | `LRC(16,4,3)` | `--data 16 --groups 4 --local-parity 1 --global-parity 3 --construction 100 --random-limit 0 --local-method random --global-method random` |
| `storage_spaces_8_2_1` | `LRC(8,2,1)` | `--data 8 --groups 2 --local-parity 1 --global-parity 1` |
| `mr_2_6_2_2` | `(g,r,a,h) = (2,6,2,2)` | `--data 6 --groups 2 --local-parity 2 --global-parity 2` |
| `mr_3_6_2_2` | `(g,r,a,h) = (3,6,2,2)` | `--data 10 --groups 3 --local-parity 2 --global-parity 2` |

Larger h=2 data-local LRCs should use `--construction 10000`; larger generic
searches may require a higher `--random-limit`, and the number of enumerated
checks can still grow quickly with larger `global_parity` values. For h=3 and
h=4, keep the status notes above as the current boundary of what has been
observed.

## Test

```bash
./tests/smoke.sh
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
