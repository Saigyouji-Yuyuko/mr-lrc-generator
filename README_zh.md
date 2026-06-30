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
- 为大规模 data-local LRC 提供专门的 `h=2` difference-pack 构造。
- 可选 stress-test 预筛，再进入精确校验。
- 支持从 JSON 读取已有矩阵并做精确 MR 校验。
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
  --construction 0 \
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
  --construction 0 \
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

对常见的 data-local `LRC(k,l,2)` 参数，`--construction 1` 会启用专门的 `h=2`
difference-pack 构造。它可以比通用随机搜索路径更容易构造大得多的 LRC：

```bash
./build/mr-lrc-generator \
  --data 64 \
  --groups 8 \
  --local-parity 1 \
  --global-parity 2 \
  --construction 1 \
  --random-limit 1 \
  --step-time 0
```

成功输出里会出现 `candidate_source=difference_pack_h2`。

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
| `--construction N` | 随机搜索前最多尝试 N 次 `h=2` difference-pack 构造。默认：`0`，表示关闭构造。 |
| `--random-limit N` | 最大候选尝试次数。默认：无界 `uint64` 最大值。 |
| `--prefilter-count N`, `--random-prefilter N`, `--stress-prefilter N` | 每个候选进入精确校验前先跑的 stress-test 模式数。默认：`0`，表示关闭预筛。 |
| `-t`, `--thread-count N` | 并行搜索 worker 数。默认：`1`，最大：`256`。 |
| `--step_time N`, `--step-time N` | 每 N 秒向 stderr 打印带时间戳的搜索数和严格完成检查数。默认：`30`；`0` 表示关闭。 |
| `--json FILE`, `--matrix-json FILE` | 将找到的结果以 pretty JSON 写入 `FILE`，包含 stdout 元数据、local group 布局，以及十进制/十六进制矩阵。标准输出仍保持普通文本报告。 |
| `--check-json FILE`, `--verify-json FILE`, `--input-json FILE` | 从 JSON 读取矩阵，并精确校验它是否为 MR-LRC。该模式不需要搜索参数。 |
| `--cauchy-dedup` | 对全 Cauchy 候选使用规范化 Cauchy 参数 key 跳过重复验证。默认关闭。 |
| `-h`, `--help` | 打印 CLI 帮助。 |

`--construction N` 会在随机候选搜索前最多尝试 N 次已注册的 data-local 构造；`0`
表示关闭构造。当前注册的构造是 `difference_pack_h2`，适用于 GF(256) 上
`local_parity=1`、`global_parity=2` 的 data-local 参数：它为每个 group 选择标签，
要求不同 group 的两两 xor 差集互不相交，local row 使用全 1，global row 使用
`(t, t^2)`。当这些 group 能放进互补二进制子空间时
（`sum ceil(log2(group_data + 1)) <= 8`），会先使用确定性的子空间打包；随后会尝试
确定性的 4 维 spread 打包，它覆盖最多 17 个 group、且每组最多 15 个 data symbol。
等价地，如果 `r` 表示包含 1 个 local parity symbol 的本地组总大小，那么这条确定性
spread 路径要求 `r <= 16`。如果两个确定性路径都不适用，再回退到 GF(256) 上的随机
差集打包；实现层面只要求 `group_data + 1 <= 256`，但这条随机回退不提供同等的确定性
规模保证。在确定性的 spread 情形下，例如 `LRC(64,8,2)` 这样的参数可以被直接构造出来，
而不是靠随机搜索碰出来。旧布尔写法仍作为兼容别名保留（`true` = `1`，`false` = `0`）。
实验性的
all-symbol skew-polynomial 原型隔离在
`src/all_symbol_skew_lrc.cpp`，没有加入构建。

## 构造方法

`cauchy` 和 `vandermonde` 使用常见的 GF(256) 矩阵族生成系数。实现会选择不同的
非零域参数，避免基础构造参数碰撞。

`column_multiplier_cauchy` 只支持 global。它会先生成 global Cauchy 矩阵，再为每个
data column 乘一个独立选择的非零 GF(256) 系数。

开启 `--cauchy-dedup` 后，当 local 和 global 方法都为 `cauchy` 时，搜索会为每个
候选保存规范化后的 Cauchy 参数 key。这个 key 会消去 Cauchy 的整体 xor 平移，并
忽略 global parity row 之间、以及同一个 local group 内 local parity row 之间的
排列，因此等价候选不会重复验证。

`random` 会用随机 GF(256) 系数填充对应 parity row。local random row 使用非零
系数，以保持 residual 消元使用的 local pivot 不变量；global random row 仍允许
出现 `00` 系数。精确 MR checker 仍然是候选是否被接受的唯一依据。

local 和 global 方法彼此独立。例如，`--local-method random --global-method
column_multiplier_cauchy` 会生成非零随机 local parity row 和 column-scaled Cauchy
global parity row。

## 校验 JSON 矩阵

使用 `--check-json FILE` 可以直接校验已有的系统 data-local 生成矩阵，不启动随机
搜索：

```bash
./build/mr-lrc-generator --check-json matrix.json
```

如果矩阵是 MR-LRC，命令返回 `0`；如果不是，返回 `2`。`--json` 输出的文件可以
直接作为输入。也可以使用下面这种最小格式：

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

`matrix` 使用十进制 GF(256) 字节，范围为 `0..255`。也可以提供
`matrix_hex`，其中元素是类似 `"0a"` 的十六进制字符串。如果 JSON 里有 `groups`
数组，会按其中的 `data` 布局校验；否则 data symbol 会和生成器一样均匀分到各组。
校验器要求 data row 是 identity，且 local parity row 只能覆盖本组 data symbol。

## 规模说明

- 当 `local_parity=1` 且 `global_parity=2` 时，优先使用 `--construction 1`。
  这条路径对互补二进制子空间，以及“最多 17 个 group、本地组总大小 `r <= 16`
  （即 `group_data <= 15`）”的 4 维 spread 都有确定性覆盖，因此像 `LRC(64,8,2)`
  这样的 wide data-local LRC 已经是实际可构造的参数。
- 通用搜索/校验路径仍会随着 `global_parity` 增大而迅速变贵。对于
  `global_parity >= 3`，建议先从小参数开始，并考虑打开 `--prefilter-count`。
- `global_parity >= 4` 仍偏实验性质，不建议直接用于大参数。

## 校验模型

总校验块数由工具自动计算：

```text
total_parity = groups * local_parity + global_parity
```

最大擦除模式大小为 `total_parity`，且必须在每个 local group 至少擦除
`local_parity` 个符号，然后把剩余擦除预算用于全局。checker 会枚举这些最大擦除
模式；只有所有被枚举的情形都可恢复时，候选才会被接受。

`--prefilter-count` 会在搜索开始前生成对应数量的 stress-probe recipe，然后每个
候选进入精确 checker 前都按这些 recipe 预筛。对常见的 `local_parity=1`、
`groups=2..6`、`global_parity=2..4` 区间，probe schedule 会优先尝试 residual
line 冲突、vector-in-plane、plane 相交、hyperplane-vector，以及少量多 vector
rank defect。targeted schedule 用完后，剩余 budget 会回退到 residual random 最大
擦除抽样。预筛拒绝的候选一定不是 MR，因为每个 witness 都用精确 verifier 相同的
local 消元 residual rank 标准确认；通过预筛的候选仍会继续进入精确 checker 和
full gate。

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
| `h2_subspace_28_2_2` | `LRC(28,2,2)` | `--data 28 --groups 2 --local-parity 1 --global-parity 2 --construction 1` |
| `h2_spread_64_8_2` | `LRC(64,8,2)` | `--data 64 --groups 8 --local-parity 1 --global-parity 2 --construction 1` |
| `storage_spaces_8_2_1` | `LRC(8,2,1)` | `--data 8 --groups 2 --local-parity 1 --global-parity 1` |
| `mr_2_6_2_2` | `(g,r,a,h) = (2,6,2,2)` | `--data 6 --groups 2 --local-parity 2 --global-parity 2` |
| `mr_3_6_2_2` | `(g,r,a,h) = (3,6,2,2)` | `--data 10 --groups 3 --local-parity 2 --global-parity 2` |

更大的 `h=2` data-local LRC 应该使用 `--construction 1`；更大的通用搜索可能需要
更高的 `--random-limit`，且被枚举的校验数量仍会随着更大的 `global_parity` 快速增长。

## 测试

```bash
./tests/smoke.sh
```

## 许可证

本项目使用 MIT License。见 [LICENSE](LICENSE)。
