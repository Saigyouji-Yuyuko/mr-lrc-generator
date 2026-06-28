# MR-LRC Generator

[English](README.md)

MR-LRC Generator 是一个 C++17 命令行工具，用于生成并精确校验 GF(256)
data-local MR-LRC 生成矩阵。工具使用 ISA-L 做 GF(256) 运算和矩阵求逆。

工具生成系统码形式的编码矩阵：

```text
[ I_k ]
[ local parity rows ]
[ global parity rows ]
```

数据符号会尽量均匀地分到各个 local group。local parity row 只覆盖本组
data symbol，global parity row 覆盖整个 stripe。候选矩阵由相互独立的
local/global 构造方法生成；之后，checker 会枚举最大擦除情形，并在 GF(256) 上
校验这些擦除是否可恢复。

## 功能

- Data-local 系统码 LRC 布局。
- local 和 global 构造方法可独立选择。local row 支持 `cauchy`、`vandermonde`
  或 `random`；global row 额外支持 `column_multiplier_cauchy`。
- 对所有最大擦除情形做精确 MR 校验。
- 可选 `--seed` 用于确定性随机搜索；不传时会自动生成并打印 seed。
- 通过 `--random-limit` 限制搜索次数。
- 通过 `--thread-count` 并行搜索。
- 使用 ISA-L 做 GF(256) 运算。

## 构建

```bash
cmake -S . -B build -DMR_LRC_REQUIRE_ISAL=ON
cmake --build build
```

如果 ISA-L 安装在非默认路径：

```bash
cmake -S . -B build \
  -DMR_LRC_REQUIRE_ISAL=ON \
  -DISAL_INCLUDE_DIR=/usr/include/isa-l \
  -DISAL_LIBRARY=/usr/lib/x86_64-linux-gnu/libisal.so
```

## 使用

生成并校验一个 `6 data + 2 local parity + 2 global parity` 的码。这个例子使用
2 个 local group，每组 1 个 local parity row，并使用 2 个 global parity row：

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

local 和 global parity row 可以使用不同构造方法：

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

简写 `--method` 会同时设置 local 和 global 方法：

```bash
./build/mr-lrc-generator --data 6 --groups 2 --local-parity 1 \
  --global-parity 2 --seed 1 --method cauchy
```

因为 `column_multiplier_cauchy` 只支持 global，请使用 `--global-method` 来选择它，
不要用 `--method`。

## 参数

必填参数：

| 参数 | 含义 |
| --- | --- |
| `-k`, `--data K` | 系统 data symbol 数。 |
| `-g`, `--groups G` | local group 数；data 会尽量均匀分组。 |
| `-r`, `--local-parity R` | 每组 local parity symbol 数。 |
| `-p`, `--global-parity M` | global parity symbol 数。 |

可选参数：

| 参数 | 含义 |
| --- | --- |
| `-s`, `--seed S` | 随机种子；不传时会自动生成并打印。 |
| `--local-method M` | local parity 构造方法：`cauchy`、`vandermonde` 或 `random`。 |
| `--global-method M` | global parity 构造方法：`cauchy`、`column_multiplier_cauchy`、`vandermonde` 或 `random`。 |
| `-m`, `--method M` | 同时设置 local 和 global 方法的别名。 |
| `--construction BOOL` | 启用已注册的 data-local 构造。默认：`true`。 |
| `--random-limit N` | 最大候选尝试次数。默认：无界 `uint64` 最大值。 |
| `-t`, `--thread-count N` | 并行搜索 worker 数。默认：`1`，最大：`256`。 |
| `--step_time N`, `--step-time N` | 每 N 秒向 stderr 打印带时间戳的搜索进度。默认：`30`；`0` 表示关闭。 |
| `--json FILE`, `--matrix-json FILE` | 将找到的矩阵以 pretty JSON 写入 `FILE`。标准输出仍保持普通文本报告。 |
| `--cauchy-dedup` | 对全 Cauchy 候选使用规范化 Cauchy 参数 key 跳过重复验证。默认关闭。 |
| `-h`, `--help` | 打印 CLI 帮助。 |

`--construction` 支持 `true/false`、`on/off`、`1/0` 和 `yes/no`。当前
data-local 主流程没有注册闭式构造，因此会回退到随机候选搜索。实验性的
all-symbol skew-polynomial 原型隔离在 `src/all_symbol_skew_lrc.cpp`，没有加入构建。

## 构造方法

`cauchy` 和 `vandermonde` 使用常见的 GF(256) 矩阵族生成系数。实现会选择不同的
非零域参数，避免基础构造参数碰撞。

`column_multiplier_cauchy` 只支持 global。它会先生成 global Cauchy 矩阵，再为每个
data column 乘一个独立选择的非零 GF(256) 系数。

开启 `--cauchy-dedup` 后，当 local 和 global 方法都为 `cauchy` 时，搜索会为每个
候选保存规范化后的 Cauchy 参数 key。这个 key 会消去 Cauchy 的整体 xor 平移，并
忽略 global parity row 之间、以及同一个 local group 内 local parity row 之间的
排列，因此等价候选不会重复验证。

`random` 会用完全随机的 GF(256) 系数填充对应 parity row，允许出现 `00` 系数。
精确 MR checker 仍然是候选是否被接受的唯一依据。

local 和 global 方法彼此独立。例如，`--local-method random --global-method
column_multiplier_cauchy` 会生成随机 local parity row 和 column-scaled Cauchy
global parity row。

## 规模说明

- 当 `data_cnt > 12` 时，不建议选择 `global_parity >= 3`。
- 当前不支持 wide LRC。
- `data_cnt >= 20` 的情形没有测试过。

## 校验模型

总校验块数由工具自动计算：

```text
total_parity = groups * local_parity + global_parity
```

最大擦除模式大小为 `total_parity`，且必须在每个 local group 至少擦除
`local_parity` 个符号，然后把剩余擦除预算用于全局。checker 会枚举这些最大擦除
模式；只有所有被枚举的情形都可恢复时，候选才会被接受。

`--random-limit` 限制候选矩阵的总尝试次数。`--thread-count` 只控制从同一尝试预算
中领取任务的 worker 数，不会放大总尝试次数。多线程时，第一个完成并通过校验的
worker 获胜，因此报告的 attempt 是第一个完成的成功尝试，不一定是编号最低的成功尝试。

## 参数示例

工业里常见的 `LRC(k,l,h)` 风格参数可以直接映射为：

```text
--data k --groups l --local-parity 1 --global-parity h
```

代表性 smoke-test 参数：

| 名称 | 记号 | CLI 参数 |
| --- | --- | --- |
| `azure_lrc_6_2_2` | `LRC(6,2,2)` | `--data 6 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_8_2_2` | `LRC(8,2,2)` | `--data 8 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_12_4_2` | `LRC(12,4,2)` | `--data 12 --groups 4 --local-parity 1 --global-parity 2` |
| `storage_spaces_8_2_1` | `LRC(8,2,1)` | `--data 8 --groups 2 --local-parity 1 --global-parity 1` |
| `mr_2_6_2_2` | `(g,r,a,h) = (2,6,2,2)` | `--data 6 --groups 2 --local-parity 2 --global-parity 2` |
| `mr_3_6_2_2` | `(g,r,a,h) = (3,6,2,2)` | `--data 10 --groups 3 --local-parity 2 --global-parity 2` |

更大的码可能需要更高的 `--random-limit`；被枚举的校验数量仍会随着更大的
`global_parity` 快速增长。

## 测试

```bash
./tests/smoke.sh
```

## 许可证

本项目使用 MIT License。见 [LICENSE](LICENSE)。
