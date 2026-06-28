# MR-LRC Generator

[中文说明](README_zh.md)

MR-LRC Generator is a C++17 command-line tool for generating and exhaustively
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
global construction methods, then every maximum erasure pattern is enumerated.
If every remaining `k x k` decoding matrix is nonsingular over GF(256), the
candidate is reported as an MR-LRC matrix.

## Features

- Data-local systematic LRC layout.
- Independent local and global construction methods:
  `cauchy`, `vandermonde`, or `random`.
- Exhaustive MR check over all maximum erasure cases.
- Deterministic random search from `--seed`.
- Search bound via `--random-limit`.
- Parallel search via `--thread-count`.
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

## Options

Required parameters:

| Option | Meaning |
| --- | --- |
| `-k`, `--data K` | Number of systematic data symbols. |
| `-g`, `--groups G` | Number of local groups. Data is split nearly evenly. |
| `-r`, `--local-parity R` | Local parity symbols per group. |
| `-p`, `--global-parity M` | Global parity symbols. |
| `-s`, `--seed S` | Random seed. |

Optional parameters:

| Option | Meaning |
| --- | --- |
| `--local-method M` | Local parity construction: `cauchy`, `vandermonde`, or `random`. |
| `--global-method M` | Global parity construction: `cauchy`, `vandermonde`, or `random`. |
| `-m`, `--method M` | Alias that sets both local and global methods. |
| `--construction BOOL` | Enable registered data-local constructions. Default: `true`. |
| `--random-limit N` | Maximum candidate attempts. Default: `1000`. |
| `-t`, `--thread-count N` | Parallel search worker count. Default: `1`, max: `256`. |
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

`random` fills the selected parity rows with fully random GF(256) coefficients,
including possible `00` coefficients. The exhaustive MR checker is still the
source of truth for whether the candidate is accepted.

Local and global methods are independent. For example, `--local-method random
--global-method cauchy` creates fully random local parity rows and Cauchy global
parity rows.

## Search Model

The total parity count is computed as:

```text
total_parity = groups * local_parity + global_parity
```

A maximum erasure pattern has size `total_parity` and must erase at least
`local_parity` symbols in every local group, then spend the remaining erasure
budget globally. The checker enumerates every such erasure pattern and verifies
that the surviving `k` rows form a nonsingular GF(256) matrix.

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

Larger codes may require a higher `--random-limit` and can become expensive
because erasure pattern enumeration grows quickly.

## Test

```bash
./tests/smoke.sh
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
