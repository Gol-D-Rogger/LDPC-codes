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
   - `flip_thr = (1<<VN_BITS)-1-3/7`（3bit → max-3，8bit → max-7），`weak=flip_thr-4`，`strong=weak` 起始。
   - 若 `strobes>0`：`address = syndrome_weight>>5`（饱和 63），按行数索引 `likelihood_init_coef_all` 得系数，分段累加 `delta_sum`，按 `VN_BITS` 进行衰减得到 `strong`（8bit 场景约 0.75 系数），四级 `level` 以 `strong/weak` 填充。
   - 若 `strobes>1`：用 `likelihood_init_fraction` 对 `level[1]/[2]` 做分数插值，拉开 4 个等级。
   - `likelihood_map[8]` 把软位组合 (00/01/10/11/...) 映射到 `level[x]`。

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
- 基本更新：`flipped ? (likelihood - weight) : (likelihood + weight - 1)`。
- 后处理：  
   - 若 `post_trigger`：  
     * `likelihood_new == flip_thr` → 设为 `flip_thr+1`（巩固已翻转状态，增加“反转回去”的难度）  
     * `likelihood_new < flip_thr` → 拉到 `flip_thr-1`（推向阈值边缘，降低翻转门槛）
   - 若 `post_trigger2`：未翻转且 `weight==1` 且 `likelihood_new==flip_thr-1` → +1（轻推向阈值，加速收敛）。
- 饱和到 `[min_likelihood, max_likelihood]`。
- 判决翻转：`likelihood_new >= flip_thr`。

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

### D. 翻转与 syndrome 更新
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
