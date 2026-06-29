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
- Optional stress-test prefilter before exact verification.
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
  --construction false \
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
  --construction false \
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
| `--construction BOOL` | Enable registered data-local constructions. Default: `true`. |
| `--random-limit N` | Maximum candidate attempts. Default: unbounded `uint64` max. |
| `--prefilter-count N`, `--random-prefilter N`, `--stress-prefilter N` | Stress-test patterns to run before exact verification for each candidate. Default: `0` disables the prefilter. |
| `-t`, `--thread-count N` | Parallel search worker count. Default: `1`, max: `256`. |
| `--step_time N`, `--step-time N` | Print timestamped searched and strict-complete counts to stderr every N seconds. Default: `30`; `0` disables it. |
| `--json FILE`, `--matrix-json FILE` | Write the found result as pretty JSON to `FILE`, including stdout metadata, local groups, and decimal/hex matrices. Standard output keeps the normal text report. |
| `--cauchy-dedup` | Skip duplicate all-Cauchy candidates using canonical Cauchy parameter keys. Default: disabled. |
| `-h`, `--help` | Print CLI help. |

`--construction` accepts `true/false`, `on/off`, `1/0`, and `yes/no`. The
current data-local flow does not register a closed-form construction, so the
main path falls back to random candidate search. The experimental all-symbol
skew-polynomial prototype is isolated in `src/all_symbol_skew_lrc.cpp` and is
not part of the build.

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

## Scale Notes

- For `data_cnt > 12`, using `global_parity >= 3` is not recommended.
- Using `global_parity >= 4` is not recommended.
- Wide LRC layouts are not supported; `data_cnt > 24` has not been tested.

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
| `storage_spaces_8_2_1` | `LRC(8,2,1)` | `--data 8 --groups 2 --local-parity 1 --global-parity 1` |
| `mr_2_6_2_2` | `(g,r,a,h) = (2,6,2,2)` | `--data 6 --groups 2 --local-parity 2 --global-parity 2` |
| `mr_3_6_2_2` | `(g,r,a,h) = (3,6,2,2)` | `--data 10 --groups 3 --local-parity 2 --global-parity 2` |

Larger codes may require a higher `--random-limit`; the number of enumerated
checks can still grow quickly with larger `global_parity` values.

## Test

```bash
./tests/smoke.sh
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
