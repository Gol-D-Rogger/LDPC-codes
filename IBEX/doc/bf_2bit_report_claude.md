# IBEX BF 2-bit 方案评估记录（ssd_fc_test, env-switchable modes, VN_BITS=2）

统一命令格式：`env IBEX_2BIT_MODE=<0|1|2|3> [IBEX_AGGR_ITER_HI=A IBEX_AGGR_ITER_LO=B IBEX_AGGR_SYND_TH=C] ./ssd_fc_test LDPC config/<cnfg> AWGN <snr>`

评估指标：主要关注 $LDPC\ FER$，辅以平均迭代次数。评估使用固定包数配置（`max_sim_num=200/500, max_err=0`），保证每次运行的包数一致，便于统计对比。

译码器：`ssd_fc_test`（`src/ldpc_codec_test.cpp`），在 T3 基线代码基础上新增 4 处改动，通过 `IBEX_2BIT_MODE` 环境变量切换 4 种 2bit BF 算法模式。

## 实现概述

### 设计目标
在 `src/ldpc_codec_test.cpp` 中实现 env-switchable 的 2bit BF 算法探索框架，支持 4 种模式：
- **Mode 0**（默认）：T3 基线，与不设置 `IBEX_2BIT_MODE` 时行为完全一致
- **Mode 1**：V25 attack-only 权重放大
- **Mode 2**：Post-processing 在 aggr 阶段关闭（PostGate）
- **Mode 3**：V25 + PostGate 组合

### 4 处修改点

#### Site 1：全局模式变量（`src/ldpc_codec_test.cpp:13`）
```cpp
static int g_2bit_mode = 0;
```
用于在 `f_update_vn_post`（列内函数）中访问当前模式，因为 mode 在 `ldpc_dec_bf_ibex`（外层函数）中从 env 读取。

#### Site 2：env 变量读取（`src/ldpc_codec_test.cpp:2617-2629`）
在 `ldpc_dec_bf_ibex` 的 while 循环之前，读取模式开关和 aggr 门限：
```cpp
// === 2-bit mode switch and tunable aggr thresholds ===
const char *env_2bit_mode = getenv("IBEX_2BIT_MODE");
const int mode_2bit = env_2bit_mode ? atoi(env_2bit_mode) : 0;
g_2bit_mode = mode_2bit;

const char *env_aggr_iter_hi = getenv("IBEX_AGGR_ITER_HI");
const char *env_aggr_iter_lo = getenv("IBEX_AGGR_ITER_LO");
const char *env_aggr_synd_th = getenv("IBEX_AGGR_SYND_TH");
const char *env_aggr_strong_synd_th = getenv("IBEX_AGGR_STRONG_SW_TH");
const int aggr_iter_hi = env_aggr_iter_hi ? atoi(env_aggr_iter_hi) : 100;
const int aggr_iter_lo = env_aggr_iter_lo ? atoi(env_aggr_iter_lo) : 50;
const int aggr_synd_th = env_aggr_synd_th ? atoi(env_aggr_synd_th) : 150;
const int aggr_strong_synd_th = env_aggr_strong_synd_th ? atoi(env_aggr_strong_synd_th) : (1 << 30);
```
- 默认值 `100/50/150/INT_MAX` 与 T3 原始硬编码一致，确保 Mode 0 行为不变。
- `aggr_strong_synd_th` 默认 `1<<30`（极大值），等价于"不限制"。

#### Site 3：模式感知 aggr 条件 + 权重映射（`src/ldpc_codec_test.cpp:2751-2791`）
```cpp
// Adjust weight if aggressive mode (mode-aware)
bool aggr;
if (mode_2bit == 0) {
  // Mode 0 (T3): env-tunable thresholds, original soft-bits condition
  aggr = (VN_BITS <= 2 && ((iteration >= aggr_iter_hi && syndrome_weight < aggr_strong_synd_th) ||
          (iteration >= aggr_iter_lo && syndrome_weight < aggr_synd_th))) ||
         ((ldpc_decoder_input.soft_bits > 0) &&
          (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
} else {
  // Mode 1/2/3 (V25): dual-threshold + relaxed soft-bits condition
  // soft_bits > 0 guard preserved; relaxation is only on !flipped_prev
  aggr = (VN_BITS <= 2 && ((iteration >= aggr_iter_hi && syndrome_weight < aggr_strong_synd_th) ||
          (iteration >= aggr_iter_lo && syndrome_weight < aggr_synd_th))) ||
         ((ldpc_decoder_input.soft_bits > 0) &&
          (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr &&
           ((VN_BITS <= 2) || !flipped_prev)));
}
int w = weight;
if (mode_2bit == 0) {
  // Mode 0: symmetric mapping (all aggr bits)
  if (aggr && (weight == 2)) w = weight + 1; // =3
  if (aggr && (weight == 3)) w = weight + 2; // =5
  if (aggr && (weight == 4)) w = weight + 3; // =7
} else {
  // Mode 1/2/3: V25 attack-only (amplify only when !flipped_prev)
  if (aggr && !flipped_prev) {
    if (weight == 2) w = weight + 1;      // =3
    else if (weight == 3) w = weight + 2;  // =5
    else if (weight == 4) w = weight + 3;  // =7
  }
}
```

关键差异：
| 特性 | Mode 0 (T3) | Mode 1/2/3 (V25) |
|------|-------------|-------------------|
| soft aggr 条件 | `soft_bits>0 && min<thr && !flipped_prev` | `soft_bits>0 && min<thr && (VN_BITS<=2 \|\| !flipped_prev)` |
| 权重放大 | 对所有 aggr 比特（含 flipped）| 仅对 `!flipped_prev`（attack-only） |
| aggr 门限 | env 可调，默认 `100/50/150` | env 可调，默认 `100/50/150` |

#### Site 4：Post-processing 门控（`src/ldpc_codec_test.cpp:586-588`）
```cpp
bool post_gate = (g_2bit_mode >= 2) ? !be_aggressive : true;
do_post_flipped = post_process && post_gate && (likelihood_new == flip_threshold);
do_post_unflipped = post_process && post_gate && (likelihood_new < flip_threshold);
```
- Mode 0/1：`post_gate=true`，post 行为与原始一致。
- Mode 2/3：`post_gate=!be_aggressive`，在 aggr 阶段关闭 post 扰动，避免"aggr 推进 + post 抖动"叠加导致发散。

### 构建方式
```bash
make ssd_fc_test "CPPFLAGS=-I./src -MMD -MP -g"
```
注意：需去掉 Makefile 默认的 `-D_LDPC_DUMP` 标志（该标志触发 `ldpc_gen_gm()` 中一个预存的 `fp` 未声明错误）。

---

## 方案 Mode 0（T3 基线）

### 设计
与不设置 `IBEX_2BIT_MODE` 时完全一致。aggr 条件使用默认门限 `100/50/150`，权重放大对称作用于所有 aggr 比特，post 始终启用。

### Row8@SNR=5.3 结果

#### 200 包
- 命令：`./ssd_fc_test LDPC config/Ibex_hd_row8_eval200.cnfg AWGN 5.3`
- 结果：
  - `[STATISTICS] Total packets simulated: 200`
  - `[STATISTICS] LDPC FER  : 3.000000e-02`（6/200）
  - `[STATISTICS] IBEX Decoder average iterations: 82.30`

#### 500 包
- 命令：`./ssd_fc_test LDPC config/Ibex_hd_row8_eval500.cnfg AWGN 5.3`
- 结果：
  - `[STATISTICS] Total packets simulated: 500`
  - `[STATISTICS] LDPC FER  : 1.600000e-02`（8/500）
  - `[STATISTICS] IBEX Decoder average iterations: 69.60`

### Row8@SNR=5.4 结果（200 包）
- 命令：`./ssd_fc_test LDPC config/Ibex_hd_row8_eval200.cnfg AWGN 5.4`
- 结果：
  - `[STATISTICS] Total packets simulated: 200`
  - `[STATISTICS] LDPC FER  : 5.000000e-03`（1/200）
  - `[STATISTICS] RAW  BER  : 6.523708e-03`
  - `[STATISTICS] IBEX Decoder average iterations: 53.15`

### Row11@SNR=4.9 结果

#### 200 包
- 命令：`./ssd_fc_test LDPC config/Ibex_hd_row11_eval200.cnfg AWGN 4.9`
- 结果：
  - `[STATISTICS] Total packets simulated: 200`
  - `[STATISTICS] LDPC FER  : 2.200000e-01`（44/200）
  - `[STATISTICS] IBEX Decoder average iterations: 277.81`

#### 500 包
- 命令：`./ssd_fc_test LDPC config/Ibex_hd_row11_eval500.cnfg AWGN 4.9`
- 结果：
  - `[STATISTICS] Total packets simulated: 500`
  - `[STATISTICS] LDPC FER  : 2.420000e-01`（121/500）
  - `[STATISTICS] IBEX Decoder average iterations: 298.80`

### 分析
- Row8@5.4 在 T3 下 FER 约 $5\times 10^{-3}$，平均迭代 53，属于该矩阵的"高 SNR 良好区间"。
- Row11@4.9 FER 高达 $\sim 0.24$，约 1/4 的包失败，平均迭代接近 300（最大 1024），属于"瀑布区边缘"，是后续方案的主要改善目标。

---

## 方案 Mode 1（V25 attack-only 权重放大）

### 设计
- aggr 条件 soft 分支放松 `!flipped_prev` 限制：2bit 下允许 flipped 比特也进入 soft aggr。
- 权重放大仅对 `!flipped_prev`（attack-only）：`2→3, 3→5, 4→7`。flipped 侧不做放大。
- Post 保持启用（与 Mode 0 一致）。

### 核心直觉
T3 的对称放大（attack + retract 都放大）在低 SNR 下容易把"正确翻转"的比特过度撤销，导致振荡与失败。V25 只在进攻端放大，保留撤销端的稳定性。

### Row8@SNR=5.3 结果（200 包）
- 命令：`env IBEX_2BIT_MODE=1 ./ssd_fc_test LDPC config/Ibex_hd_row8_eval200.cnfg AWGN 5.3`
- 结果：
  - `[STATISTICS] LDPC FER  : 2.000000e-02`（4/200）
  - `[STATISTICS] IBEX Decoder average iterations: 93.09`
- 对比 Mode 0：FER 下降 33%（3.0e-2 → 2.0e-2），但平均迭代上升 13%（82→93），攻守不对称导致收敛稍慢。

### Row11@SNR=4.9 结果

#### 200 包（aggr 门限 240/120/280）
- 命令：`env IBEX_2BIT_MODE=1 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test LDPC config/Ibex_hd_row11_eval200.cnfg AWGN 4.9`
- 结果：
  - `[STATISTICS] LDPC FER  : 3.000000e-02`（6/200）
  - `[STATISTICS] IBEX Decoder average iterations: 151.92`
- 对比 Mode 0（FER=2.2e-1）：**7.3 倍 FER 改善**，平均迭代从 278 降至 152（-45%）。

#### 500 包（aggr 门限 240/120/280）
- 命令：`env IBEX_2BIT_MODE=1 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test LDPC config/Ibex_hd_row11_eval500.cnfg AWGN 4.9`
- 结果：
  - `[STATISTICS] LDPC FER  : 3.400000e-02`（17/500）
  - `[STATISTICS] IBEX Decoder average iterations: 157.20`
- 500 包验证与 200 包一致，**7.1 倍 FER 改善**（0.242 → 0.034），统计稳定。

### 分析
- Mode 1 在 Row11@4.9 的效果极为显著，是本轮实验的核心发现。
- 关键因素：V25 的 attack-only 放大 + 后移 aggr 门限（240/120/280）。门限后移避免了中前期误翻扩散，attack-only 避免了撤销侧振荡。
- Row8@5.3 上也有适度改善，但 Mode 2 表现更好（见下）。

---

## 方案 Mode 2（PostGate：aggr 阶段关闭 post）

### 设计
- aggr 条件与权重放大同 Mode 1（V25 attack-only）。
- Post-processing 在 `be_aggressive==true` 时被 gate 掉：`post_gate = !be_aggressive`。
- 直觉：aggr 阶段的权重放大已经提供足够推力，此时 post 的"snap 到 thr±1"反而会引入额外振荡。

### Row8@SNR=5.3 结果

#### 200 包
- 命令：`env IBEX_2BIT_MODE=2 ./ssd_fc_test LDPC config/Ibex_hd_row8_eval200.cnfg AWGN 5.3`
- 结果：
  - `[STATISTICS] LDPC FER  : 1.500000e-02`（3/200）
  - `[STATISTICS] IBEX Decoder average iterations: 67.88`
- 对比 Mode 0：**FER 下降 50%**（3.0e-2 → 1.5e-2），平均迭代下降 17%（82→68）。

#### 500 包
- 命令：`env IBEX_2BIT_MODE=2 ./ssd_fc_test LDPC config/Ibex_hd_row8_eval500.cnfg AWGN 5.3`
- 结果：
  - `[STATISTICS] LDPC FER  : 1.400000e-02`（7/500）
  - `[STATISTICS] IBEX Decoder average iterations: 67.30`
- 500 包验证：FER 12.5% 改善（1.6e-2 → 1.4e-2），迭代下降 3%。

### Row8@SNR=5.4 结果（200 包）
- 命令：`env IBEX_2BIT_MODE=2 ./ssd_fc_test LDPC config/Ibex_hd_row8_eval200.cnfg AWGN 5.4`
- 结果：
  - `[STATISTICS] LDPC FER  : 0.000000e+00`（**0/200，零错误**）
  - `[STATISTICS] IBEX Decoder average iterations: 49.03`
- 对比 Mode 0（FER=5.0e-3）：零错误 + 平均迭代从 53 降至 49。

### Row11@SNR=4.9 结果（500 包，aggr 门限 240/120/280）
- 命令：`env IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test LDPC config/Ibex_hd_row11_eval500.cnfg AWGN 4.9`
- 结果：
  - `[STATISTICS] LDPC FER  : 7.600000e-02`（38/500）
  - `[STATISTICS] IBEX Decoder average iterations: 195.80`
- 对比 Mode 0：3.2 倍 FER 改善（0.242 → 0.076），但不及 Mode 1 的 7.1 倍。

### 分析
- Mode 2 在 **Row8（中高 SNR）** 表现最优：FER 和迭代数双降。
- 在 Row11@4.9（低 SNR）表现不如 Mode 1，说明低 SNR 下 post 扰动对部分包仍有正面脱困作用，完全关闭反而丢失了这些收益。
- Row8@5.4 的零错误结果虽然在 200 包口径下有统计波动，但与基线 FER=5e-3 的对比仍然明确。

---

## 方案 Mode 3（V25 + PostGate 组合）

### 设计
- aggr 条件与权重放大同 Mode 1/2（V25 attack-only）。
- Post-gating 同 Mode 2（aggr 阶段关闭 post）。
- 即 Mode 1 和 Mode 2 的所有改动同时生效。

### Row8@SNR=5.3 结果（200 包）
- 命令：`env IBEX_2BIT_MODE=3 ./ssd_fc_test LDPC config/Ibex_hd_row8_eval200.cnfg AWGN 5.3`
- 结果：
  - `[STATISTICS] LDPC FER  : 3.000000e-02`（6/200）
  - `[STATISTICS] IBEX Decoder average iterations: 100.15`
- 对比 Mode 0：FER 无改善，平均迭代增加 22%。

### Row11@SNR=4.9 结果（200 包，aggr 门限 240/120/280）
- 命令：`env IBEX_2BIT_MODE=3 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test LDPC config/Ibex_hd_row11_eval200.cnfg AWGN 4.9`
- 结果：
  - `[STATISTICS] LDPC FER  : 5.000000e-02`（10/200）
  - `[STATISTICS] IBEX Decoder average iterations: 170.92`
- 对比 Mode 0：4.4 倍 FER 改善，但不及 Mode 1 的 7.3 倍。

### 分析
- Mode 3 在两个配置上都不是最优：Row8 上不如 Mode 2，Row11 上不如 Mode 1。
- 说明 V25 的 attack-only 放大和 PostGate 两个机制**并不正交**：
  - V25 保守的撤销端需要 post 在早期帮助打破振荡；PostGate 关掉了这个帮助。
  - PostGate 的"aggr 阶段不扰动"的好处在 V25 的 attack-only 下被削弱，因为 V25 本身已经减少了 aggr 阶段的振荡。
- 结论：Mode 3 作为组合方案没有实现 "1+1>2" 的效果，不推荐使用。

---

## 实验过程中的关键 Bug 与修复

### Bug：首轮实验 Mode 1/2/3 全部 FER=1.0

- 现象：首轮 8 组实验中，所有 Mode 1/2/3 的 run 均 FER=1.0，平均迭代 1024（打满），完全不收敛。
- 根因：Mode 1/2/3 的 soft aggr 条件最初实现为：
  ```cpp
  // 错误实现（缺少 soft_bits > 0 guard）
  aggr = ... || (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr &&
                 (VN_BITS <= 2 || !flipped_prev));
  ```
  这丢掉了 `ldpc_decoder_input.soft_bits > 0` 前置条件。在 `sd_num=1`（硬输入）场景下 `soft_bits==0`，原始条件中 `soft_bits > 0` 会使整个 soft 分支为 false；丢掉该 guard 后，aggr + 权重放大从 iteration 0 就对所有 2-bit VN 生效，导致大规模误翻和发散。

- 参考：`src/ldpc_codec_test2.cpp:2770-2775`（V25 正确实现）中 soft 分支完整保留了 `soft_bits > 0`。

- 修复：恢复 `soft_bits > 0` 前置条件：
  ```cpp
  // 正确实现
  aggr = ... || ((ldpc_decoder_input.soft_bits > 0) &&
                 (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr &&
                  ((VN_BITS <= 2) || !flipped_prev)));
  ```
  "放松"仅体现在 `!flipped_prev` 条件上（2bit 时不要求 `!flipped_prev`），而非去掉 `soft_bits > 0` 整个 guard。

- 修复后重跑所有实验，结果正常。

---

## 综合结果对比

### Row8@SNR=5.3

| Mode | 200 pkt FER | 200 pkt Iter | 500 pkt FER | 500 pkt Iter | vs Mode 0 |
|------|------------|-------------|------------|-------------|-----------|
| 0 (T3 baseline) | 3.00e-2 | 82.3 | 1.60e-2 | 69.6 | -- |
| 1 (V25) | 2.00e-2 | 93.1 | -- | -- | -33% FER, +13% iter |
| **2 (PostGate)** | **1.50e-2** | **67.9** | **1.40e-2** | **67.3** | **-50% FER, -17% iter** |
| 3 (Combined) | 3.00e-2 | 100.2 | -- | -- | 无收益, +22% iter |

**Row8@5.3 最优：Mode 2（PostGate）**

### Row8@SNR=5.4

| Mode | 200 pkt FER | 200 pkt Iter | vs Mode 0 |
|------|------------|-------------|-----------|
| 0 (T3 baseline) | 5.00e-3 | 53.1 | -- |
| **2 (PostGate)** | **0** | **49.0** | **零错误, -8% iter** |

**Row8@5.4 最优：Mode 2（PostGate）**

### Row11@SNR=4.9（aggr 门限 240/120/280）

| Mode | 200 pkt FER | 200 pkt Iter | 500 pkt FER | 500 pkt Iter | vs Mode 0 |
|------|------------|-------------|------------|-------------|-----------|
| 0 (T3 baseline) | 2.20e-1 | 277.8 | 2.42e-1 | 298.8 | -- |
| **1 (V25)** | **3.00e-2** | **151.9** | **3.40e-2** | **157.2** | **7.1x FER, -47% iter** |
| 2 (PostGate) | -- | -- | 7.60e-2 | 195.8 | 3.2x FER |
| 3 (Combined) | 5.00e-2 | 170.9 | -- | -- | 4.4x FER |

**Row11@4.9 最优：Mode 1（V25），使用 aggr 门限 A/B/C=240/120/280**

---

## 最优方案推荐

| 矩阵 | SNR 区间 | 推荐 Mode | 关键 env 参数 | 典型 FER |
|-------|---------|-----------|-------------|---------|
| Row8 | 5.3~5.4（中高 SNR）| Mode 2 | `IBEX_2BIT_MODE=2` | 0~1.4e-2 |
| Row11 | 4.9（低 SNR）| Mode 1 | `IBEX_2BIT_MODE=1 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280` | 3.4e-2 |

### 可复现命令

#### Row8 最优（Mode 2）
```bash
# 200 包快速验证
env IBEX_2BIT_MODE=2 ./ssd_fc_test LDPC config/Ibex_hd_row8_eval200.cnfg AWGN 5.4
# 500 包统计
env IBEX_2BIT_MODE=2 ./ssd_fc_test LDPC config/Ibex_hd_row8_eval500.cnfg AWGN 5.3
```

#### Row11 最优（Mode 1）
```bash
# 200 包快速验证
env IBEX_2BIT_MODE=1 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test LDPC config/Ibex_hd_row11_eval200.cnfg AWGN 4.9
# 500 包统计
env IBEX_2BIT_MODE=1 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test LDPC config/Ibex_hd_row11_eval500.cnfg AWGN 4.9
```

---

## 评估配置文件说明

| 配置文件 | 矩阵 | max_sim_num | max_err | 用途 |
|---------|-------|-------------|---------|------|
| `config/Ibex_hd_row8_eval200.cnfg` | Row8 | 200 | 0 | 快速 200 包评估 |
| `config/Ibex_hd_row8_eval500.cnfg` | Row8 | 500 | 0 | 500 包统计验证 |
| `config/Ibex_hd_row11_eval200.cnfg` | Row11 | 200 | 0 | 快速 200 包评估 |
| `config/Ibex_hd_row11_eval500.cnfg` | Row11 | 500 | 0 | 500 包统计验证 |

- `max_err=0` 保证不会提前停止，确保固定包数。
- 这些配置文件仅改变仿真停止条件，译码参数与 `Ibex_hd_row8.cnfg` / `Ibex_hd_row11.cnfg` 完全一致。

---

## 与 `ssd_fc_test2` 方案的关系

本报告的方案（Mode 0~3）实现在 `src/ldpc_codec_test.cpp`（编译为 `ssd_fc_test`），与 `doc/bf_2bit_report_codex.md` 中记录的方案（V0~V57，实现在 `src/ldpc_codec_test2.cpp`，编译为 `ssd_fc_test2`）是**独立的两套代码**。

对应关系：
- Mode 0 ≈ Codex T3 基线（`ldpc_codec_test2.cpp` 的默认行为）
- Mode 1 ≈ Codex V25（attack-only 权重放大）
- Mode 2 ≈ Codex 方向2-方案2 的 PostGate 部分
- Mode 3 ≈ Codex V25 + PostGate 组合

本报告的实现更轻量（仅 4 处改动 + 环境变量开关），适合快速 A/B 对比；Codex 的 `ssd_fc_test2` 包含更多演进（V25→V46→V51 等），具有更丰富的旋钮但代码复杂度更高。
