# IBEX BF 解码详解（C++ 参考实现）

目标：完整梳理 `src/ldpc_codec.cpp` / `src/ldpc_codec.h` 中 IBEX 位翻转解码的流程、关键参数与与 RTL 对应关系，便于对照硬件验证与调参。

## 模块与调用链
- `ldpc_dec_bf_ibex`：核心解码循环（列/比特扫描、似然更新、syndrome 更新、早停/后处理）。
- `ldpc_ibex_input`：把信道输出映射到 `s_ldpc_decoder_input`（硬判决 + 软位）。
- `ldpc_ibex_parameters`：装载后处理开关、初始似然系数、软位映射表、早停阈值等。
- `ldpc_config`（512-bit 分支）：生成 `s_h_matrix`（占用/fade、移位量、mask、行列权），供编码/解码。
- 辅助：`f_likelihood_levels`（初始似然）、`f_update_vn_post`（翻转/后处理更新）、`f_check_nodes`（syndrome）、`f_check_node_weight`（syndrome weight）、`f_256_bit_lfsr`/`f_512_bit_lfsr`（后处理随机扰动）。

## 矩阵与数据结构
- `s_h_matrix`：`rows(≤13)`、`cols(=80)`、`bits(=256/512)`；`occupied/fade` 标记循环移位子矩阵；`element[i][j]` 为移位量；`mask[j][k]` 处理尾部非整字节 parity；`parity_column` 标记校验列；`row_weight/col_weight` 用于调度。
- `s_variable_nodes`：每比特 `bit_hard`、`likelihood`（宽度由 `VN_BITS`）、`flipped`。
- `s_check_nodes`：syndrome 展开到 bit 级（与 `element` 对齐）。
- `s_likelihood_levels`：`level[4]`、`flip_thr`、`weak/strong/min/max`。
- `s_ldpc_decoder_input`：`iteration_limit`、`post_iteration`、`nand_strobes`(0/1/2，对应 rd_num 1/3/≥5)、`soft_bits`(0/1/2)、`syndrome_cal_only` 以及 `corrupted_codeword`。
- `s_ldpc_decoder_parameters`：后处理开关/阈值 (`post_process_en`, `syndrome_weight_thr_qc/post`, `post_ratio`)、似然初值系数 `likelihood_init_coef_all`（按行数选择）、分数插值 `likelihood_init_fraction`、软位映射 `likelihood_map`、翻转强度阈值 `likelihood_thr`、早停阈值表 `early_terminate_thr[2][8]`（按行数与 `nand_strobes` 选）。

## 输入与初始似然
1. 硬判决：`rx_blk>=0 ? 0 : 1` 写入 `bit_hard`；尾部未用位填 0。
2. 软位：根据 `rd_num` 选择 1 或 2 bit（`bit_questionable/bit_questionable2`），未配置的复制最弱等级（`likelihood_map[3]`）。
3. 初始似然 `f_likelihood_levels(strobes=soft_bits, rows, syndrome_weight)`：
   该函数负责在解码开始前，根据当前的错误严重程度（Syndrome Weight）和信道软信息，动态计算并初始化所有 VN 的似然值等级。它相当于解码器的“初始状态配置器”。

   - **基础初始化**: 
     * 设定 `max` (7) 和 `min` (1)。
     * 设定 `flip_thr = (1<<VN_BITS)-1-3/7`（3bit → max-3，即 4）。这意味着初始确信的比特需要至少 3 个校验报错才会翻转。
     * 初始 `weak` / `strong` 设为 `flip_thr - 4`。

   - **动态调整 (`strobes>0`)**:
     * 根据 `syndrome_weight`（错误总数）微调初始置信度。错误越多，信道越差，初始置信度应设得越低（越容易翻转）。
     * 计算逻辑：`address = syndrome_weight >> 5`，查表 `likelihood_init_coef_all` 获取系数，计算衰减量 `delta_total`，进而调整 `strong` 等级。

   - **软信息插值 (`strobes>1`)**:
     * 若 `weak` 和 `strong` 差值足够大（`>=8`），利用 `likelihood_init_fraction` 对中间等级 `level[1]/[2]` 进行插值，拉开梯度。

   - **输出**: 返回包含 `flip_thr`, `min/max` 和 `level[4]` 查找表的结构体。VN 初始化时根据软信息查表：`likelihood = level[soft_bit_val]`。

**代码摘录：初始似然生成（`src/ldpc_codec.cpp`）**
```cpp
likelihood_levels.max = (1 << VN_BITS) - 1;
likelihood_levels.flip_thr = (VN_BITS == 3) ? likelihood_levels.max - 3 : likelihood_levels.max - 7;
likelihood_levels.weak = likelihood_levels.flip_thr - 4;
...
address = syndrome_weight >> 5;
...
delta_total = (VN_BITS == 8) ? ((delta_sum >> 1) + (delta_sum >> 2))
                             : ((delta_sum >> 2) + (delta_sum >> 3));
likelihood_levels.strong = likelihood_levels.weak - delta_total;
...
likelihood_levels.min = likelihood_levels.strong;
return likelihood_levels;
```

## 前置检查
- `f_check_nodes` 计算初始 syndrome；`f_check_node_weight` 得 `syndrome_weight`。
- 早停：若 `early_terminate_dis==0` 且 `syndrome_weight >= early_terminate_thr[strobes][rows_bin]`（rows_bin 由行数分档 6..13），立即返回失败（`early_termination=1`）。
- `syndrome_cal_only`：只算 syndrome，不做翻转，直接返回。
- `syndrome_weight==0`：直接成功。

## 迭代主循环（按列、按比特）
循环条件：`iteration < iteration_limit` 且未完成/未放弃。

### A. PRNG / 后处理调度
- 在 `post_iteration` 首列加载种子 `prng_init[32]` 到 256/512 位 LFSR；之后每列 `f_256_bit_lfsr` / `f_512_bit_lfsr` 递推。
- `syndrome_weight_delayed`：5 级移位寄存器模拟流水线，后处理判断使用延迟值。
- 触发判定与生效：
  - `post_trigger`: `(syndrome_weight_delayed < syndrome_weight_thr_qc) && (iteration%16 < post_ratio) && post_process_en`（不含 PRNG）
  - `post_trigger2`: `(syndrome_weight_delayed < syndrome_weight_thr_post) && (iteration%16 >= post_ratio) && post_process_en`（不含 PRNG）
  - 实际生效标志：`prng_post_process = post_trigger && PRNG_bit`，`prng_post_process2 = post_trigger2 && PRNG_bit`，后续似然更新使用这两个带随机性的标志。

**代码摘录：后处理触发与 PRNG 掩码**
```cpp
post_trigger = hamming_weight_lt_circ_thr &&
               ((iteration % 16) < ldpc_decoder_parameters.post_ratio) &&
               ldpc_decoder_parameters.post_process_en;
post_trigger2 = hamming_weight_lt_post_thr &&
                ((iteration % 16) >= ldpc_decoder_parameters.post_ratio) &&
                ldpc_decoder_parameters.post_process_en;
...
prng_post_process = post_trigger &&
                    (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_qc) &&
                    ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
prng_post_process2 = post_trigger2 &&
                     (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_post) &&
                     ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
```

### B. 列/比特扫描
- 跳过无效尾部：当存在 `extra_bits_of_parity` 且列为最后两列时，超出有效位不参与更新。
- 权重计算 `weight`：遍历所有行，使用 `(k + bits - element[i][j]) % bits` 取对应 syndrome 位；按 `occupied/fade/mask` 规则计数未满足的校验。
- 侵略模式 `be_aggressive`：`soft_bits>0 && likelihood_levels.min < likelihood_thr && !post_trigger && !post_trigger2`。仅对“未翻转过”的比特放大权重：0→0, 1→1, 2→3, 3→5, 4→7（代码中 1→1 也走专门分支，效果相同）。

**代码摘录：侵略模式权重放大**
```cpp
bool aggr = (ldpc_decoder_input.soft_bits > 0) &&
            (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev);
int w = weight;
if (aggr && (weight == 0))
    w = weight + 0;
if (aggr && (weight == 1))
    w = weight + 0; // =1
if (aggr && (weight == 2))
    w = weight + 1; // =3
if (aggr && (weight == 3))
    w = weight + 2; // =5
if (aggr && (weight == 4))
    w = weight + 3; // =7
```

### C. 似然更新（`f_update_vn_post`）

该函数是 IBEX BF 算法中计算下一个状态似然值（Likelihood Update）的核心逻辑单元。它的输入包括当前的似然值、计算出的权重、以及各种后处理和激进模式的标志。它的输出是新的似然值。

#### 1. 基础更新逻辑 (Basic Update)
这是算法的主干，决定了似然值的基本走向。
*   **如果当前已翻转 (`flipped == true`)**：
    *   说明之前认为该比特是错的。现在的 `weight` 代表依然有多少个校验方程报错。
    *   **逻辑**：`likelihood - weight`。报错越多，对“翻转”这个决定的信心越低（似然值下降）。
*   **如果当前未翻转 (`flipped == false`)**：
    *   说明之前认为该比特是对的。现在的 `weight` 代表有多少个校验方程指控该比特是错的。
    *   **逻辑**：`likelihood + weight - 1`。报错越多，该比特是错的可能性越大（似然值上升）。此处 `-1` 为阻尼策略，防止增长过快。

#### 2. 后处理强制干预 (Post-Processing)
如果触发了随机扰动（`post_process`），算法会强行修改似然值：
*   **巩固翻转 (`do_post_flipped`)**: 若 `new_likelihood == flip_thr`，强制设为 `flip_thr + 1`。既然已到边缘，就往里推一步，防止轻易退回。
*   **诱导翻转 (`do_post_unflipped`)**: 若 `new_likelihood < flip_thr`，强制拉到 `flip_thr - 1`。强行拉到悬崖边，下一轮只要有微小报错（Weight>=1）就会翻转。

#### 3. 激进后处理 (Post-Process 2 / Kick)
针对 `Weight=1` 的顽固错误进行的特殊打击。
*   若未翻转、`Weight=1` 且似然值在阈值边缘 (`flip_thr - 1`)，强制 `+1` 触发翻转。这用于消除单线连接的死锁错误。

#### 4. 饱和截断
最后保证数值不溢出 `[min, max]` 范围。

**代码摘录：似然更新与后处理**
```cpp
likelihood_new = flipped ? likelihood - weight : likelihood + weight - 1;
do_post_flipped = post_process && (likelihood_new == flip_threshold);
do_post_unflipped = post_process && (likelihood_new < flip_threshold);
if (do_post_unflipped)
    likelihood_new = flip_threshold - 1;
else if (do_post_flipped)
    likelihood_new = flip_threshold + 1;

if (post_process2 && !flipped && (weight == 1) && (likelihood_new == flip_threshold - 1))
    likelihood_new++;
...
return likelihood_new;
```

### D. 2-bit 位宽下的逻辑影响分析
当 `VN_BITS` 缩减为 2 时（范围 0~3，`flip_thr=3`），`f_update_vn_post` **函数本身的代码逻辑不需要修改**，因为它依赖的是参数传入的 `flip_threshold`、`min`、`max`，具有自适应性。但其实际行为会发生以下变化：

1.  **数值空间极度压缩**：
    *   `flip_thr - 1` 变为 2 (Weak)。
    *   `flip_thr + 1` 变为 4 (超出 Max=3，会被最后的饱和逻辑截断回 3)。
    *   这意味着后处理中的“巩固翻转”操作 (`flip_thr+1`) 实际上变成了维持在最大值 3，不再有额外的缓冲空间。

2.  **后处理行为**：
    *   `do_post_unflipped` (拉向边缘): 会将似然值强制设为 2。下一轮只要 `Weight >= 2` (未翻转时 `L_new = 2 + 2 - 1 = 3`) 即可触发翻转。
    *   `post_process2` (临门一脚): 当 `Likelihood=2` 且 `Weight=1` 时，强制 +1 变为 3，触发翻转。这在 2-bit 下依然有效且关键。

3.  **饱和截断**：
    *   由于动态范围小，`likelihood_new` 很容易触碰 `min` (1) 或 `max` (3)。代码末尾的 `if (likelihood_new <= min)` 和 `if (likelihood_new >= max)` 将频繁生效，保证数值安全。

**结论**：该函数在 2-bit 模式下逻辑完备，无需改动代码，只需确保传入正确的阈值参数（`min=1`, `max=3`, `flip_thr=3`）。

### E. 2-bit 自适应更新策略 (Adaptive Update Strategy)
针对 2-bit 位宽下数值空间极小（0~3）导致的震荡问题，IBEX 引入了一种基于全局收敛趋势的自适应更新策略。该策略通过动态调整更新步长（Delta），在“维持翻转”和“撤销翻转”之间取得平衡。

#### 1. 核心逻辑
利用全局信号 `pushing`（Syndrome Weight 是否未下降）来决定奇数权重的取整方向。

**Delta 计算公式**:
```cpp
int delta;
if (VN_BITS <= 2)
    // 2-bit 模式: 基础步长减半 (阻尼)。
    // 若趋势不好 (pushing=true) 且权重为奇数，则向上取整 (加大力度)；否则向下取整 (保守)。
    delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
else
    delta = weight; // 常规模式
```

#### 2. 状态更新
```cpp
likelihood_new = flipped ? (likelihood - delta) : (likelihood + delta - 1);
```
*   **已翻转 (`flipped`)**: 执行减法。
    *   **趋势好 (`!pushing`)**: `delta` 较小（向下取整）。少减一点，**保护当前的翻转状态**（惯性）。
    *   **趋势坏 (`pushing`)**: `delta` 较大（向上取整）。多减一点，**果断撤销翻转**（纠错）。
*   **未翻转 (`!flipped`)**: 执行加法。趋势不好时加速翻转，打破僵局。

#### 3. 效果示例 (以 Weight=3, Likelihood=3 为例)
*   **场景 A (收敛中)**: `delta = 1`。`New = 3 - 1 = 2` (Weak)。状态温和回退，避免直接变回 Strong。
*   **场景 B (震荡中)**: `delta = 2`。`New = 3 - 2 = 1` (Strong)。状态大幅回退，强制撤销翻转。

### F. 翻转与 syndrome 更新
- 若翻转状态改变：沿该比特关联的所有 1 边把 `cn` 对应位取反（同 RTL mask/toggle 逻辑，含 `fade/mask` 分支）。
- 每列末重算 `syndrome_weight`；若 0 则提前成功。

### E. 迭代推进
- 列完成递增 `clock_cycles`；外层循环结束后再次检查 syndrome 作为最终状态。

## 输出与计数
- `dec_do_blk = corrupted ^ flipped`（早停或 syndrome-only 不改写）。
- 统计翻转数：`errors_in_codeword`（全局），`errors_in_userdata`（仅信息位）；上限 4095。
- `iterations` 为实际循环次数，`clock_cycles` 粗略计数（2*cols + 迭代列步）。
- `cw_fail` 由残余 syndrome 是否为 0 决定；早停直接标失败，syndrome-only 标成功但不纠错。

## 与 RTL 文档关键对照
- Step1/2：`ldpc_ibex_input` + `f_check_nodes` + `f_likelihood_levels` 对应初始化、syndrome 计算与似然分级。
- Step3：列/比特扫描、侵略模式权重放大、后处理扰动、翻转更新 syndrome；`syndrome_weight_delayed` 对应 RTL 管线延迟。
- Step4：`dec_do_blk`、翻转计数与失败标志输出；翻转上限保护与早停逻辑与 RTL 行为一致。

## Mermaid 流程图
```mermaid
flowchart TD
    A[输入映射\n初始 syndrome] --> B{syndrome==0?}
    B -- 是 --> Z[成功输出]
    B -- 否 --> C[早停阈值检查]
    C -- 触发 --> Y[早停失败返回]
    C -- 继续 --> D[似然初值 f_likelihood_levels]
    D --> E[迭代 iteration<limit]
    E --> F[列扫描 j]
    F --> G[比特扫描 k\n跳过无效尾部]
    G --> H[统计 weight\n(occupied/fade/mask)]
    H --> I[侵略模式放大?]
    I --> J[似然更新 f_update_vn_post\n含 post_trigger/post_trigger2+PRNG]
    J --> K{likelihood>=flip_thr?}
    K -- 状态变化 --> L[翻转比特\n更新 syndrome]
    K -- 无变化 --> M[syndrome 保持]
    L --> N[重算 syndrome_weight]
    M --> N
    N --> O{syndrome==0?}
    O -- 是 --> Z
    O -- 否 --> P[列/迭代继续]
    P --> E
    Z --> Q[纠错输出\n翻转计数/迭代数]
```
