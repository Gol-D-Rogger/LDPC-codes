# IBEX BF 解码器 2-bit 优化方案 (VN_BITS = 2)

本文档描述了将 IBEX BF (Bit-Flipping) 解码器的变量节点位宽 (`VN_BITS`) 从默认的 3-bit (或 8-bit) 缩减为 **2-bit** 的修改方案。此优化旨在减少硬件存储资源（Memory Footprint），同时尽可能保持解码性能。

## 1. 背景与目标

*   **当前状态**: 默认 `VN_BITS = 3`，变量节点似然值范围为 0~7 (3 bits)。每个 VN 节点在硬件中至少需要 5 bits 存储 (1 Hard + 3 Likelihood + 1 Flip)。
*   **目标状态**: `VN_BITS = 2`，变量节点似然值范围为 0~3 (2 bits)。
*   **预期收益**: 每个变量节点 (VN) 的核心存储需求从 5 bits 降低为 **4 bits** (1 Hard + 2 Likelihood + 1 Flip)，存储节省约 **20%**。对于大码长（如 N=32768），这将显著减少片上 SRAM 的面积。

## 2. 代码修改方案详解

### 2.1 VN_BITS修改
传入参数的config中进行修改

### 2.2 似然值分级初始化逻辑适配
**文件**: `src/ldpc_codec.cpp`
**函数**: `ldpc_packet::f_likelihood_levels`

**修改原因**:
原逻辑使用固定偏移量（如 `max - 3`）计算翻转阈值。在 2-bit 模式下，`max` 仅为 3。
*   如果沿用 `max - 3`，则 `flip_thr` = 0。这意味着所有比特在初始化时就会处于翻转边缘，导致解码器立即开始全网翻转，产生剧烈震荡且无法收敛。
*   我们需要手动重新规划 0~3 的数值空间，使其既能区分强/弱置信度，又能保留合理的翻转门槛。

**修改代码段**:

```cpp
likelihood_levels.max = (1 << VN_BITS) - 1;
likelihood_levels.min = 1; // 保持最小值 min=1，有效范围 {1, 2, 3}

// [新增] 2-bit 模式专用初始化
if (VN_BITS <= 2) {
    // 2-bit 动态范围极小配置: 
    // 1: Strong (强置信度) - 初始化值，离翻转阈值有 2 的距离 (需要 Weight>=3 才能翻转)
    // 2: Weak   (弱置信度) - 离翻转阈值有 1 的距离 (需要 Weight>=2 才能翻转)
    // 3: Flip Threshold (翻转阈值) - 达到此值即翻转
    likelihood_levels.flip_thr = 3; 
    likelihood_levels.weak     = 2;
    likelihood_levels.strong   = 1;
} 
// [原有逻辑保持不变]
else if (VN_BITS == 3) {
    // ...
```

### 2.3 软信息分级插值阈值调整
**文件**: `src/ldpc_codec.cpp`
**函数**: `ldpc_packet::f_likelihood_levels`

**修改原因**:
原代码中包含一个逻辑：`if (weak_minus_strong >= 8)`，用于判断是否有足够的动态范围来对软信息进行更细粒度的插值（Fractional Interpolation）。
*   在 2-bit 模式下，`weak (2) - strong (1) = 1`。
*   由于 `1 < 8`，原逻辑会直接跳过软信息插值步骤，导致软信息（Soft Bits）无法被充分利用，所有比特都将被初始化为相同的似然值，损失了纠错性能。
*   必须将此门限下调为 `1`，以允许在极其有限的 2-bit 空间内也能区分不同的软信息等级。

**修改代码段**:

```cpp
    // ...
    weak_minus_strong = likelihood_levels.level[3] - likelihood_levels.level[0];
    
    // [修改] 引入自适应阈值
    // 对于 2-bit，只要有 1 的差值 (2-1=1) 就允许进行分级映射
    int split_threshold = (VN_BITS <= 2) ? 1 : 8; 
    
    if (weak_minus_strong >= split_threshold) {
        // ... 内部插值计算逻辑保持不变 ...
    }
```

### 2.4 禁用激进模式 (Aggressive Mode)
**文件**: `src/ldpc_codec.cpp`
**函数**: `ldpc_packet::ldpc_dec_bf_ibex`

**修改原因**:
激进模式的核心思想是人为放大权重（例如 Weight=2 视为 3，Weight=3 视为 5）。
*   在 3-bit+ 模式下，似然值范围大（0~7），+2 或 +3 的增量只是推动似然值快速靠近阈值。
*   在 2-bit 模式下，似然值范围极小（1~3）。
    *   当前似然值通常为 `1`。
    *   如果 Weight=2 被激进地放大为 `w=3`。
    *   新似然值 = `1 + 3 - 1 = 3` (直接翻转)。
*   这会导致只有 2 个校验方程报错的比特也被立即翻转。在低信噪比下，这种“过度反应”会引入大量新错误，导致解码发散。
*   因此，2-bit 模式下必须保持“保守”，严格遵守物理校验结果。

**修改代码段**:

```cpp
bool aggr = (ldpc_decoder_input.soft_bits > 0) &&
            (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev);

// [新增] 2-bit 模式下强制关闭激进增益，防止数值饱和导致的震荡
if (VN_BITS <= 2) {
    aggr = false;
}

int w = weight;
// ... 后续 if (aggr) w = ... 逻辑保持不变
```

## 3. 验证与预期行为

修改完成后，解码器的预期行为如下：

1.  **数值范围**: 所有 VN 的 `likelihood` 值将被严格限制在 `[1, 3]` 之间。
2.  **翻转逻辑**:
    *   初始状态下，VN 似然值为 `1` (Strong)。
    *   **Case 1 (Weight=1)**: `New_L = 1 + 1 - 1 = 1` (不变)。
    *   **Case 2 (Weight=2)**: `New_L = 1 + 2 - 1 = 2` (变 Weak，但不翻转)。这意味着 2-bit 方案比原方案（原方案可能在激进模式下翻转）更**保守**。
    *   **Case 3 (Weight=3)**: `New_L = 1 + 3 - 1 = 3` (>= Flip_Thr，**翻转**)。
    *   **结论**: 只有当至少 3 个校验方程同时报错时，才会触发比特翻转。这种高门槛有助于消除 Error Floor，但可能会略微增加所需的迭代次数。

## 4. 风险与对策

*   **收敛速度变慢**: 由于不再使用激进模式，翻转动作变少。建议在仿真中适当增加 `max_iteration` 上限（例如从 50 增加到 60）来补偿。
*   **Trap Sets 敏感**: 状态空间的缩减可能导致更容易落入某些 Trap Sets。IBEX 原有的 `Post-Processing` 机制（基于 PRNG 的随机扰动）在 2-bit 模式下至关重要，**务必不要禁用 Post-Processing 功能**。