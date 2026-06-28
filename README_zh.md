# MR-LRC Generator

[English](README.md)

这是一个 C++17 写的 GF(256) MR-LRC 矩阵生成和穷举校验工具，依赖 ISA-L 做 GF(256) 运算和矩阵求逆。

工具暴露一组基础参数，以及局部/全局两种构造方法：

- 数据块数 `--data`
- local 组数 `--groups`
- 每组 local 校验块数 `--local-parity`
- 全局校验块数 `--global-parity`
- 随机种子 `--seed`
- 随机候选次数上限 `--random-limit`
- 搜索线程数 `--thread-count`
- 局部块随机构造方法 `--local-method`
- 全局块随机构造方法 `--global-method`
- 候选生成模式 `--construction`

编码矩阵是系统码形式：

```text
[ I_k ]
[ local parity rows ]
[ global parity rows ]
```

数据块会按组数尽量均分到各个 local group。当前主流程保持 data-local 输出语义：local parity row 只依赖本组 data。搜索流程里，局部校验行和全局校验行可以分别选择 `cauchy`、`vandermonde` 或 `random` 构造，默认都是 `cauchy`。候选矩阵生成后，工具会枚举所有最大损失模式，并检查每种模式剩余的 `k x k` 行矩阵在 GF(256) 上是否非奇异。全部通过时，认为这组矩阵通过 MR-LRC 校验。

## 构建

```bash
cmake -S . -B build -DMR_LRC_REQUIRE_ISAL=ON
cmake --build build
```

如果需要显式指定系统 ISA-L：

```bash
cmake -S . -B build \
  -DMR_LRC_REQUIRE_ISAL=ON \
  -DISAL_INCLUDE_DIR=/usr/include/isa-l \
  -DISAL_LIBRARY=/usr/lib/x86_64-linux-gnu/libisal.so
```

## 使用

`6+2+2`，即 `6 data + 2 local parity + 2 global parity`。这里有 2 个 local group，每组 1 个 local parity，所以总校验块数是 `2 + 2 = 4`：

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

若要分开指定 local 和 global 构造方法：

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

`--local-method` 和 `--global-method` 都支持：

- `cauchy`
- `vandermonde`
- `random`

其中 `random` 表示对应校验行的所有系数都在 GF(256) 上全随机生成，允许系数为 `00`。工具会继续用完整 erasure pattern 枚举来判断这次随机候选是否通过 MR-LRC 校验。

`--construction` 是布尔开关，支持 `true/false/on/off/1/0/yes/no`，默认是 `true`。它用于启用已注册的 data-local 显式构造；当前还没有启用任何 data-local 构造，因此主流程会直接进入搜索。之前实验过的 skew-polynomial 论文构造是 all-symbol LRC 拓扑，global parity 符号也属于 local group；相关原型已经隔离在 `src/all_symbol_skew_lrc.cpp`，没有加入构建，也没有任何当前调用路径。

默认搜索示例：

```bash
./build/mr-lrc-generator \
  --data 6 \
  --groups 2 \
  --local-parity 1 \
  --global-parity 2 \
  --seed 1
```

也可以显式关闭构造注册点：

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

`--random-limit` 是全局候选矩阵尝试次数上限。`--thread-count` 只控制并行搜索线程数，不会放大总尝试次数；例如 `--random-limit 100 --thread-count 4` 表示最多尝试 100 个候选，由 4 个线程并行领取任务。

总校验块数由工具自动计算：

```text
total_parity = groups * local_parity + global_parity
```

最大损失模式定义为：每个 local group 至少丢失 `local_parity` 个符号，再额外消耗 `global_parity` 个全局预算。工具会枚举所有这类大小为 `groups * local_parity + global_parity` 的 erasure pattern。

## 候选参数

当前 CLI 使用系统码 / data-local 布局：

```text
data symbols + per-group local parity symbols + global parity symbols
```

工业和论文里常见的 `LRC(k,l,h)` 记号可以直接映射到：

```text
--data k --groups l --local-parity 1 --global-parity h
```

其中 `l` 是 local parity 总数。对于当前测试里的 Azure/WAS/YTsaurus 风格参数，每个 local group 使用 1 个 local parity，因此 `l` 同时就是 `--groups`。

`./tests/smoke.sh` 已加入一组短码候选，全部使用 `--seed 1 --local-method cauchy --global-method cauchy` 做确定性回归：

| 候选名 | 工业记号 | CLI 参数 |
| --- | --- | --- |
| `azure_lrc_6_2_2` | `LRC(6,2,2)` | `--data 6 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_8_2_2` | `LRC(8,2,2)` | `--data 8 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_10_2_2` | `LRC(10,2,2)` | `--data 10 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_12_2_2` | `LRC(12,2,2)` | `--data 12 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_14_2_2` | `LRC(14,2,2)` | `--data 14 --groups 2 --local-parity 1 --global-parity 2` |
| `azure_lrc_12_4_2` | `LRC(12,4,2)` | `--data 12 --groups 4 --local-parity 1 --global-parity 2` |
| `azure_lrc_8_2_3` | `LRC(8,2,3)` | `--data 8 --groups 2 --local-parity 1 --global-parity 3` |
| `azure_lrc_10_2_3` | `LRC(10,2,3)` | `--data 10 --groups 2 --local-parity 1 --global-parity 3` |
| `storage_spaces_8_2_1` | `LRC(8,2,1)` | `--data 8 --groups 2 --local-parity 1 --global-parity 1` |
| `storage_spaces_12_2_1` | `LRC(12,2,1)` | `--data 12 --groups 2 --local-parity 1 --global-parity 1` |

标准 all-symbol MR-LRC 常写成 `(g,r,a,h)`，其中维度是：

```text
k = g * (r - a) - h
n = g * r
```

这个工具不是 all-symbol parity-check 布局，但可以用相同的 `n/k` 规模做系统码搜索测试：

```text
--data "g * (r - a) - h" --groups g --local-parity a --global-parity h
```

默认 smoke 里还加入了这些代表性规模：

| 候选名 | `(g,r,a,h)` | CLI 参数 |
| --- | --- | --- |
| `mr_3_5_1_2` | `(3,5,1,2)` | `--data 10 --groups 3 --local-parity 1 --global-parity 2` |
| `mr_4_4_1_2` | `(4,4,1,2)` | `--data 10 --groups 4 --local-parity 1 --global-parity 2` |
| `mr_2_7_1_3` | `(2,7,1,3)` | `--data 9 --groups 2 --local-parity 1 --global-parity 3` |
| `mr_3_5_1_3` | `(3,5,1,3)` | `--data 9 --groups 3 --local-parity 1 --global-parity 3` |
| `mr_2_6_2_2` | `(2,6,2,2)` | `--data 6 --groups 2 --local-parity 2 --global-parity 2` |
| `mr_3_6_2_2` | `(3,6,2,2)` | `--data 10 --groups 3 --local-parity 2 --global-parity 2` |

下面这些参数更适合作为手动压力搜索或系统仿真，不放进默认 smoke：

```bash
# 当前构造下不一定能在短时间或默认 random-limit 内找到矩阵
./build/mr-lrc-generator --data 12 --groups 2 --local-parity 1 --global-parity 3 \
  --seed 1 --local-method cauchy --global-method cauchy
./build/mr-lrc-generator --data 10 --groups 2 --local-parity 1 --global-parity 4 \
  --seed 1 --local-method cauchy --global-method cauchy

# wide LRC 数值规模很大，完整 erasure pattern 枚举会非常重
./build/mr-lrc-generator --data 48 --groups 4 --local-parity 1 --global-parity 3 \
  --seed 1 --local-method cauchy --global-method cauchy
./build/mr-lrc-generator --data 120 --groups 5 --local-parity 1 --global-parity 3 \
  --seed 1 --local-method cauchy --global-method cauchy
```

## 快速测试

```bash
./tests/smoke.sh
```

## 许可证

本项目使用 MIT License，见 [LICENSE](LICENSE)。
