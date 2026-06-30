# 项目优化记录

本文记录 MR-LRC Generator 当前已经完成的主要优化，覆盖校验器、搜索流程、构造路径、CLI/JSON 可用性、测试和文档。记录依据是当前 `main` 分支历史，以及 2026-07-01 工作区内关于 `strict_complete` 输出语义的最新改动。

## 总览

项目从最初的“生成候选矩阵后直接枚举最大擦除模式”逐步优化为以下结构：

1. 先通过专门构造路径直接生成一部分参数的候选，例如 h=2 difference-pack 和 h=3 feasible-points。
2. 对通用随机搜索，先用可选 stress prefilter 快速淘汰明显失败的候选。
3. 对进入精确校验的候选，使用 residual verifier、dual verifier、h=2/h=4/双组快路径减少枚举量。
4. 最后保留 full maximal erasure gate，确保被报告为严格完成的候选确实通过完整最大擦除校验。

## 校验器优化

### Residual verifier

核心校验从直接枚举完整擦除矩阵，改为先对每个 local group 做 local 消元，把问题压缩到 global parity 维度的 residual rank 检查。

- 每个 group 按 `extra` 擦除数预计算 residual block。
- 对 local parity 已有 MDS 保证的 Cauchy 或专门构造，跳过重复 local rank 检查。
- 对相同 column space 的 residual block 去重，减少后续组合数量。
- 复用 scratch buffer，避免在大量 rank 检查中反复分配内存。
- 失败时仍记录第一组失败擦除符号，方便调试。

相关实现：`ResidualVerifier`、`build_residual_erasure_block_into()`、`precompute_group()`。

### Full maximal erasure gate

Residual verifier 通过后，仍会运行 full gate 对最大擦除模式做最终确认。

- 这保证 `strict_complete=true` 的结果不是只通过了优化条件，而是通过了完整最大擦除语义。
- 如果 residual checker 通过但 full gate 失败，会把 `strict_complete` 标为 `false` 并返回失败擦除模式。

相关实现：`FullErasureGate`、`check_mr_with_residual_verifier()`。

### h=2 residual line 快路径

当 `global_parity=2` 时，residual block 可以压缩为 projective line key。

- 单个 residual column 必须在保留的 global row 上非零。
- 不同 group 不能产生相同 projective line，否则存在不可恢复最大擦除模式。
- 用 packed key 和 dense hash table 做碰撞检测，避免展开大量组合。

相关实现：`verify_h2_fast_path()`、`projective_column_key_packed()`。

### Streaming dual verifier

对 h=3/h=4 等场景，加入 dual-space 检查，把“是否存在坏 residual 组合”转换为 projective dual vector 的可达性问题。

- 对每个 residual block 枚举左零空间的 projective vector。
- 对每个 dual vector 用动态规划判断各 group 的 `extra` 选择是否能组成坏模式。
- 先 streaming 检查，不必一开始构建完整 dual index。
- h=4 设置 streaming 上限，达到上限后回退到 indexed dual，避免内存和时间失控。

相关实现：`verify_by_streaming_dual()`、`build_dual_index()`、`verify_by_dual_index()`、`check_dual_options()`。

### 双 group residual projection 快路径

对两组参数，直接枚举保留的 global rows，并投影每组 residual block 后做拼接 rank 检查。

- 避免通用多组递归状态。
- 对 `groups=2` 的工业常见小参数更直接。
- 仍使用同一套 residual rank 标准，失败语义与通用 checker 一致。

相关实现：`verify_two_group_projection_fast_path()`、`check_two_group_projection_case()`。

### 小 h / local parity 1 快路径

对 `local_parity=1` 且 `global_parity<=4` 的高频参数，预筛和 residual 检查使用小矩阵数组路径。

- `build_residual_erasure_block_a1_small()` 避免动态 vector 热路径开销。
- `gf256_rank_small()` 用于 4x4 以内矩阵的快速 rank 判断。
- 这条路径直接服务 h=2、h=3、h=4 的搜索和预筛。

## 搜索流程优化

### 独立 local/global 构造方法

local 和 global parity row 可以独立选择构造族。

- local 支持 `cauchy`、`vandermonde`、`random`。
- global 支持 `cauchy`、`column_multiplier_cauchy`、`vandermonde`、`random`。
- `--method` 保留为同时设置 local/global 的简写。
- random local row 使用非零系数，保持 residual 消元依赖的 local pivot 不变量。

这让搜索可以组合不同结构，例如 random local + Cauchy global，或 random local + column-multiplier Cauchy global。

### Column-multiplier Cauchy global

新增 `column_multiplier_cauchy` global 方法。

- 先生成 global Cauchy 矩阵。
- 再为每个 data column 乘独立非零 GF(256) 系数。
- 只允许用于 global row，CLI 会拒绝把它用于 local row。

该方法为 global row 增加了列缩放自由度，扩大 Cauchy 系候选空间。

### Cauchy 候选去重

`--cauchy-dedup` 针对 local/global 都为 Cauchy 的搜索跳过等价候选。

- 为每个候选生成规范化 Cauchy 参数 key。
- key 消去 Cauchy 参数整体 xor 平移。
- key 忽略 global parity row 之间、以及同一 local group 内 local parity row 之间的排列。
- 使用固定大小 key 和 `ankerl::unordered_dense`，降低 hash set 开销。
- 输出 `duplicate_candidates_skipped` 和 `cauchy_dedup_key_bytes`，方便判断去重效果。

相关实现：`CauchyKeyBytes`、`CauchyDeduper`、`append_cauchy_canonical_key()`。

### Stress prefilter

`--prefilter-count` 在精确 verifier 前增加一批 stress probes。

- h=2 优先检查 line collision 和 single rank。
- h=3 优先检查 vector-in-plane、three-vectors、projected line。
- h=4 优先检查 single rank、plane-plane、hyperplane-vector、plane-vectors、four-vectors 等结构。
- targeted probe 用完后，剩余预算回退到 residual random 最大擦除抽样。
- prefilter 拒绝候选时使用与精确 checker 相同的 residual rank 标准，所以拒绝结果是可靠 witness。
- 通过 prefilter 的候选仍会进入精确 residual verifier 和 full gate。

输出字段包括 `prefilter_candidates_checked`、`prefilter_candidates_rejected`、`prefilter_patterns_checked`。

相关实现：`PrefilterPlan`、`build_prefilter_plan()`、`RandomStressPrefilter`。

### 并行搜索

`--thread-count` 支持多个 worker 并行搜索。

- 所有 worker 共享同一个 `--random-limit` 尝试预算，不会把预算乘以线程数。
- 第 `worker` 个线程检查 `worker+1, worker+1+thread_count, ...` 这些 attempt。
- 每个 worker 复用自己的 `Code`、`ResidualVerifier` 和可选 prefilter，减少线程间共享状态。
- 成功候选用原子 `stop` 通知其他 worker 尽快退出。
- 最终合并每个 worker 的 attempts、checked、strict、failed、duplicate、prefilter 统计。

相关实现：`WorkerContext`、`WorkerStats`、`PublishedWorkerStats`。

### 进度报告

`--step-time` 会周期性输出搜索进度。

- 报告 searched attempt 数和 strict-complete 候选数。
- 使用 worker 端定期发布的原子统计，避免进度线程频繁加锁。
- `--step-time 0` 可关闭进度输出，适合 smoke test 或脚本调用。

## 专门构造路径

### h=2 difference-pack 构造

`--construction N` 对 `local_parity=1`、`global_parity=2` 的 data-local LRC 优先尝试 difference-pack 构造。

矩阵形态：

- local row 全为 `1`。
- global rows 使用 `(t, t^2)`。
- 每个 local group 的标签差集不能与其他 group 冲突。

已实现路径：

| 路径 | 适用条件 | 作用 |
| --- | --- | --- |
| 互补子空间 | `sum ceil(log2(group_data + 1)) <= 8` | 确定性打包，通常 attempt 1 成功。 |
| 4 维 spread | `groups <= 17` 且 `group_data <= 15` | 覆盖宽 data-local h=2 参数，例如 `LRC(64,8,2)`。 |
| 随机差集打包 | 满足基本 h=2 形态但不在确定性覆盖内 | 使用 `--construction` 预算做随机 fallback。 |

最新语义调整：

- h=2 difference-pack 成功后由构造证明接受。
- `generate()` 不再额外跑完整擦除枚举 verifier。
- 输出 `strict_complete=false`、`patterns_checked=0`。
- 如果需要对小参数显式完整校验，可把矩阵写成 JSON 后用 `--check-json` 再跑。

相关实现：`build_difference_pack_h2_candidate()`、`check_difference_pack_h2_proof()`。

### h=3 feasible-points 构造

`--construction N` 对 `local_parity=1`、`global_parity=3` 启用实验性的 GF(256)^3 feasible-points 构造。

- 把每个 local group 看作以 local parity 为原点的点集。
- 先过滤跨组 residual vector/plane 兼容性。
- 生成候选后仍交给精确 checker，h=3 路径不是证明捷径。
- 成功输出 `candidate_source=feasible_points_h3`，且通过 exact checker 时 `strict_complete=true`。

当前已观察并精确验证的小参数包括：

- `LRC(12,2,3)`
- `LRC(14,2,3)`
- `LRC(16,2,3)`
- `LRC(12,4,3)`
- `LRC(16,3,3)`
- `LRC(16,4,3)`

当前构造/搜索实验尚未找到：

- `LRC(18,2,3)`
- `LRC(18,3,3)`
- `LRC(20,4,3)`

相关实现：`build_feasible_points_h3_candidate()`。

### h=4 实验边界

当前没有注册 h=4 构造路径。

- 小参数如 `LRC(8,2,4)`、`LRC(10,2,4)` 可作为基线。
- 当前瓶颈主要在跨 group residual-subspace packing，尤其是更多 group 下的 triple compatibility。
- 后续若要推进 h=4，建议优先设计专门 CSP/hypergraph 构造，再由现有 exact checker 验证。

## JSON 和可复现性优化

### JSON 输出

`--json FILE` / `--matrix-json FILE` 可以把找到的矩阵写成 pretty JSON。

包含内容：

- 参数元数据：data、groups、local/global parity、seed、attempt、method 等。
- 搜索统计：patterns、strict、prefilter、dedup 等。
- group 布局和 row label。
- 十进制矩阵 `matrix`。
- 十六进制矩阵 `matrix_hex`。

最新工作区改动已把 `strict_complete` 加入普通 stdout 和 JSON 输出，避免 h=2 proof shortcut 与 exact verifier 结果混淆。

### JSON 校验

`--check-json FILE` / `--verify-json FILE` / `--input-json FILE` 可以直接读取矩阵并精确校验。

- 不需要搜索参数。
- 支持 `matrix` 十进制字节。
- 支持 `matrix_hex` 十六进制字符串。
- 如果 JSON 中有 `groups`，按显式 group 布局校验；否则按生成器均匀分组。
- 要求 data rows 是 identity，local parity row 只能覆盖本组 data symbol。
- MR 返回 `0`，非 MR 返回 `2`。

这让外部系统矩阵、构造路径输出、实验样例都能用同一套 exact checker 复验。

### Seed 和 attempt 语义

- 未传 `--seed` 时会生成并打印 seed。
- 每个 attempt 使用 `seed_for_attempt(seed, attempt)` 派生独立随机流。
- 多线程下报告的 attempt 是第一个完成并通过的成功 attempt，不保证是编号最小的成功 attempt。
- `--random-limit` 是总 attempt 预算，和线程数无乘法关系。

## CLI 和文档优化

CLI 已增加或统一以下入口：

- `--local-method` / `--global-method` 独立选择构造族。
- `--method` 作为兼容简写。
- `--construction` 接受数值，也兼容 `true/on/yes` 和 `false/off/no`。
- `--random-limit` 控制通用随机搜索预算。
- `--thread-count` 控制 worker 数。
- `--prefilter-count`、`--random-prefilter`、`--stress-prefilter` 为同一预筛预算的别名。
- `--cauchy-dedup` / `--no-cauchy-dedup` 控制 Cauchy 去重。
- `--step-time` / `--step_time` 控制进度输出。
- `--json` / `--matrix-json` 输出矩阵。
- `--check-json` / `--verify-json` / `--input-json` 校验矩阵。

README 和 README_zh 已补充：

- h=2 difference-pack 的适用条件、构造路径和规模边界。
- h=3 feasible-points 的实验状态。
- h=4 尚未注册构造的边界说明。
- prefilter、thread-count、random-limit 的语义。
- representative smoke-test 参数。

## 测试覆盖

`tests/smoke.sh` 已覆盖以下优化入口：

- Cauchy、Vandermonde、Random 方法。
- local/global 方法独立选择。
- global-only `column_multiplier_cauchy`，以及错误使用时拒绝。
- JSON 输出和 `--check-json` 校验。
- 十进制/十六进制 JSON 非 MR 拒绝。
- `--cauchy-dedup`。
- `--prefilter-count` 预筛。
- 多线程搜索。
- h=2 difference-pack 构造，包括 subspace 和 spread 参数。
- h=3 feasible-points 构造。
- 工业常见 data-local LRC 参数和 `(g,r,a,h)` 参数。
- 最新 smoke 断言包含 `strict_complete=false` 的 h=2 proof shortcut，以及 `strict_complete=true` 的 h=3 exact-verified 构造。

## 推荐使用方式

| 场景 | 推荐参数 |
| --- | --- |
| 小参数基线验证 | `--construction 0 --method cauchy` |
| h=2 data-local 大参数 | `--construction 10000 --random-limit 1` |
| h=2 且 `data + h > 256` 的边界参数 | 加 `--global-method random` |
| h=3 小参数实验 | `--construction 100 --random-limit 0 --local-method random --global-method random` |
| 通用随机搜索 | 设置 `--random-limit`，按机器核数设置 `--thread-count` |
| 随机搜索失败很多 | 增加 `--prefilter-count` 先淘汰明显非 MR 候选 |
| 全 Cauchy 随机搜索 | 可试 `--cauchy-dedup` |
| 外部矩阵复验 | `--check-json matrix.json` |

## 当前限制和后续方向

- h=2 difference-pack 当前按构造证明接受，不在 `generate()` 内跑 exact verifier；这是有意的性能和规模取舍，输出用 `strict_complete=false` 明确标记。
- h=3 feasible-points 仍是启发式候选生成器，最终依赖 exact checker。
- h=4 没有注册构造，通用搜索会很快遇到 residual-subspace packing 瓶颈。
- prefilter 是拒绝优化，不是接受优化；通过 prefilter 不代表 MR，只代表需要继续 exact checker。
- 多线程提高吞吐，但不会减少单个候选的 exact checker 成本。
- 后续最值得投入的是 h=4 专门构造、更多参数族的 deterministic construction，以及对 residual/dual checker 的参数化 benchmark。

## 变更索引

| 提交/状态 | 优化主题 |
| --- | --- |
| `1618247` | residual verifier、scratch 复用、校验器结构优化。 |
| `b87e0cc` | streaming dual checks before indexing。 |
| `1afd379` | two-group residual projection fast path。 |
| `f0da6dd` | column-multiplier Cauchy global 方法。 |
| `a1359ce` | 进度报告和规模说明调整。 |
| `000e2b6` | random stress prefilter、并行统计、搜索流程优化。 |
| `91fd9f8` | JSON verifier、h=2 difference-pack 构造。 |
| `4d3820a` | h=2 构造文档整理。 |
| `7d9ec73` | h=3 feasible-points 构造。 |
| 当前工作区 | h=2 proof shortcut 的 `strict_complete=false` / `patterns_checked=0` 输出语义和 smoke 断言。 |
