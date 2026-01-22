# IBEX BF 2-bit 调优候选方案（用于 ssd_fc_test 评估）

## 背景问题
当前 2-bit IBEX BF 解码器存在**纠错性能下降过缓**的问题：
- 配置：`VN_BITS=2`，`flip_thr=2`，`max=3`，`min=0`
- 当前半步逻辑：$delta = \lfloor weight/2 \rfloor + (weight \mod 2) \land pushing$
- 问题：当 `weight=3` 时，$delta \in [1,2]$，未翻转分支 $likelihood\_new = likelihood + delta - 1$，最多只能增加 0~1，**推力严重不足**

## 限制条件
1. **不改动现有源码文件**，只复制 `src/ldpc_codec.cpp` → `src/ldpc_codec_test.cpp`
2. Makefile 新增目标 `ssd_fc_test` 链接该副本
3. 运行命令：`./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
4. 先跑 200 包观测 FER，目标：FER 接近 0

---

## 方案清单

### 方案 A：存2算3（扩展域）
**思路**：存储 2-bit，在 `f_update_vn_post` 内部将 likelihood 左移 1 位运算（变成 3-bit 域），运算完再右移回写。

**修改位置**：`f_update_vn_post()` 函数（约 561-594 行）

**伪代码**：
```cpp
int f_update_vn_post(...) {
  int scale = (VN_BITS <= 2) ? 1 : 0;
  int like = (likelihood << scale) + (scale ? 1 : 0);  // 0,1,2,3 → 1,3,5,7
  int thr = scale ? 4 : flip_threshold;                 // 3-bit 阈值
  int mn = min_likelihood << scale;
  int mx = (max_likelihood << scale) + (scale ? 1 : 0);
  
  bool flipped = (like >= thr);
  int delta = weight;  // 恢复全步长
  int like_new = flipped ? (like - delta) : (like + delta - 1);
  
  // clamp & post process...
  
  return like_new >> scale;  // 回落到 2-bit
}
```

**预期**：恢复梯度和步长精度，推力增强。

---

### 方案 B：攻守分离（全权重进攻，半权重撤销）
**思路**：保持 2-bit 存储，不扩展域。
- **未翻转分支**：全权重 `likelihood + weight - 1`（保证 weight=3 能推过阈值）
- **已翻转分支**：半权重 $\lfloor weight/2 \rfloor$ 撤销（防止快速乒乓）

**修改位置**：`f_update_vn_post()` 函数

**伪代码**：
```cpp
int delta;
if (VN_BITS <= 2) {
  if (flipped)
    delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);  // 半步撤销
  else
    delta = weight;  // 全步进攻
} else {
  delta = weight;
}
likelihood_new = flipped ? (likelihood - delta) : (likelihood + delta - 1);
```

**预期**：推得动，又不乒乓。

---

### 方案 C：半步 + 强制下沉后处理
**思路**：保持半步增量 $\lfloor(weight-1)/2\rfloor$（奇数受趋势控制）。
在 2-bit 下，若 `post_process && like_new < thr`，直接将 `like_new` 置为 `min_likelihood`（强制按死），避免 Weak/Flip 震荡。

**修改位置**：`f_update_vn_post()` 函数，`do_post_unflipped` 分支

**伪代码**：
```cpp
if (do_post_unflipped) {
  if (VN_BITS <= 2)
    likelihood_new = min_likelihood;  // 强制按死
  else
    likelihood_new = flip_threshold - 1;
}
```

**预期**：打破 error floor 死锁。

---

### 方案 D：分时后处理/扰动
**思路**：不改更新公式，按迭代分时控制后处理：
- 前半轮（iteration < post_iteration/2）：关闭 post_process
- 后半轮：开启 post_process 或减小 `post_ratio`

**修改位置**：`ldpc_dec_bf_ibex()` 函数中 post_trigger 计算逻辑（约 2628-2640 行）

**伪代码**：
```cpp
// 原逻辑
post_trigger = hamming_weight_lt_circ_thr && ((iteration % 16) < ldpc_decoder_parameters.post_ratio) && ...;

// 新逻辑：前半轮关闭
bool early_phase = (iteration < ldpc_decoder_input.post_iteration + 100);
post_trigger = !early_phase && hamming_weight_lt_circ_thr && ...;
```

**预期**：早期稳定收敛，后期解锁死锁。

---

### 方案 E：直接全权重（恢复原始逻辑）
**思路**：2-bit 下不做半步，直接用 `delta = weight`，与 3-bit 及以上完全一致。

**修改位置**：`f_update_vn_post()` 函数

**伪代码**：
```cpp
// 原代码
if (VN_BITS <= 2)
  delta = (weight>>1) + ((weight&1) && pushing ? 1 : 0);
else
  delta = weight;

// 改为
delta = weight;  // 所有位宽统一
```

**预期**：推力完全恢复，但可能出现乒乓震荡。可作为基准对比。

---

### 方案 F：动态步长（趋势自适应）
**思路**：根据 syndrome_weight 连续下降趋势动态调整步长。
- 若 `syndrome_weight` 连续下降（`trend_good`）：使用全步长
- 否则：使用半步长

**修改位置**：
1. `ldpc_dec_bf_ibex()` 中计算 `pushing` 的逻辑（约 2488 行）
2. `f_update_vn_post()` 中根据新增参数选择步长

**伪代码**：
```cpp
// 在 ldpc_dec_bf_ibex() 中计算趋势
static int last_sw[3] = {9999, 9999, 9999};
bool trend_good = (syndrome_weight < last_sw[0]) && (last_sw[0] < last_sw[1]);
last_sw[2] = last_sw[1]; last_sw[1] = last_sw[0]; last_sw[0] = syndrome_weight;

// 在 f_update_vn_post() 中
int delta;
if (VN_BITS <= 2 && !trend_good)
  delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
else
  delta = weight;
```

**预期**：收敛好时加速，停滞时保守，自适应平衡。

---

## 激进加权处理说明
代码中存在激进加权机制（约 2714-2728 行）：
```cpp
bool aggr = (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min < likelihood_thr && !flipped_prev);
int w = weight;
if (aggr && (weight == 0)) w = weight + 0;
if (aggr && (weight == 1)) w = weight + 0;
if (aggr && (weight == 2)) w = weight + 1;  // 2 → 3
if (aggr && (weight == 3)) w = weight + 2;  // 3 → 5
if (aggr && (weight == 4)) w = weight + 3;  // 4 → 7
```

**建议**：在 2-bit 方案中**关闭激进加权**，避免与新步长逻辑冲突。
```cpp
bool aggr = (VN_BITS > 2) && (soft_bits > 0) && ...;  // 仅 3-bit 及以上启用
```

---

## 趋势判据说明
变量 `pushing` 定义（约 2488-2489 行）：
```cpp
int prev_sw = (iteration == 0) ? syndrome_weight_delayed : syndrome_weight_r[3];
bool pushing = (syndrome_weight_delayed >= prev_sw);
```
含义：若当前 syndrome_weight **未下降**（停滞或上升），则 `pushing = true`。

---

## 评估方式
1. 复制 `src/ldpc_codec.cpp` → `src/ldpc_codec_test.cpp`
2. Makefile 新增 target `ssd_fc_test`
3. 每次只实现**一种方案**，运行 500 包
4. 记录 `[STATISTICS] LDPC FER`
5. 若 FER ≈ 0，增加包数（5000+）进一步验证

---

## 实现优先级建议
| 优先级 | 方案 | 理由 |
|-------|------|------|
| 1 | E | 最简单，作为基准 |
| 2 | B | 攻守分离，平衡推力与稳定 |
| 3 | A | 扩展域，理论最优但改动大 |
| 4 | F | 自适应，可能效果好 |
| 5 | C | 针对 error floor |
| 6 | D | 不改核心公式 |

