# IBEX BF “存2算3”混合精度方案（VN_BITS=2，计算域扩展）

目标：在保持存储位宽 2-bit 的前提下，内部运算域临时扩展到 3-bit（×2 缩放），以恢复梯度和步长，提升收敛速度，同时尽量不改动上层调用接口。

## 核心思路
- 存储域保持 2bit：取值 {0,1,2,3}。
- 计算域左移一位（×2）：取值 {0,2,4,6}，可容纳“半步”变化（1、3、5）。
- 翻转阈值/上下限/权重统一在计算域运算，最后右移回写。
- 激进模式在 2bit 扩展时默认关闭，防止一步顶格。

## 代码修改点

### 1) `f_update_vn_post` 支持计算域扩展
**文件**: `src/ldpc_codec.cpp`  
**声明**: `src/ldpc_codec.h`

改造签名，增加 `pushing`（半步取整判据）和 `scale2x`（是否开启×2扩展）：
```cpp
int f_update_vn_post(int likelihood, int weight, int min_likelihood, int max_likelihood,
                     bool post_process, bool post_process2, bool be_aggressive,
                     int flip_threshold, bool pushing, bool scale2x);
```

实现要点：
```cpp
const int scale = (scale2x && VN_BITS <= 2) ? 1 : 0;
const int unit  = 1 << scale;

int like = likelihood << scale;
int thr  = flip_threshold << scale;
int mn   = min_likelihood << scale;
int mx   = max_likelihood << scale;

// 半步增量：基于 (weight-1)/2，奇数权重由 pushing 决定是否进一档
int delta_half = (weight - 1) >> 1;
if ((weight & 1) && pushing) delta_half++;
int delta_scaled = delta_half * unit;

bool flipped = (like >= thr);
int like_new = flipped ? (like - delta_scaled) : (like + delta_scaled - unit);

// 后处理拉动/推高（保持原逻辑，单位换成 unit）
if (post_process && (like_new < thr))       like_new = thr - unit;
else if (post_process && (like_new == thr)) like_new = thr + unit;
if (post_process2 && !flipped && (weight == 1) && (like_new == thr - unit))
    like_new += unit;

// 饱和并回写
like_new = std::min(std::max(like_new, mn), mx);
return like_new >> scale;
```

### 2) 调用处传递 `pushing` 和 `scale2x`
**文件**: `src/ldpc_codec.cpp` (`ldpc_dec_bf_ibex` 主循环)

- 计算 `pushing`：可用宏观趋势（当前迭代 syndrome_weight 与上一迭代结束值比较），或简单列前后对比；示例：
  ```cpp
  int prev_sw = (iteration == 0) ? syndrome_weight_delayed : syndrome_weight_r[3];
  bool pushing = (syndrome_weight_delayed >= prev_sw);
  ```
- 调用：
  ```cpp
  vn.c[j].b[k].likelihood =
      f_update_vn_post(..., pushing, scale2x);
  ```
- 2bit 入口 `ldpc_dec_bf_ibex_2bit` 传 `scale2x=true`，其他位宽传 `false`。

### 3) 激进模式在扩展域默认关闭
在主循环计算 `aggr` 时加门控：
```cpp
bool aggr = (ldpc_decoder_input.soft_bits > 0) &&
            (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr) &&
            !flipped_prev && !scale2x; // 2bit 扩展时禁用激进加权
```

### 4) 可选：初值与映射的安全配置
- `f_likelihood_levels`：VN_BITS<=2 时保留梯度，推荐 `level[0]=Strong(1)`, `level[1]=Weak(2)`, `level[2]=Weak(2)`, `level[3]=Weak/边缘`，避免最强软位映射到翻转阈值。
- `likelihood_map`：确认最强软位组合映射到 `level[0]`（强），最弱映射到弱，避免直接落在阈值档。

## 行为预期
- 存储仍为 2bit，但运算域允许半步调整，推力/阻尼更细腻。
- 激进模式在 2bit 扩展下默认关闭，降低误翻风险；可按需要基于错误量或迭代阶段再打开。
- 后处理逻辑沿用原判据，但单位已换算到扩展域，随机扰动仍可打破死锁。

## 接线与验证
- VN_BITS 配置为 2；调用 2bit 入口 `ldpc_dec_bf_ibex_2bit` 时传 `scale2x=true`。
- 其他位宽走原入口，`scale2x=false`，不受影响。
- 回归验证：观察 Strong/Weak 比例、syndrome_weight 下降曲线；如仍保守，可在未翻转分支对 weight=3 做小幅特判（如 delta_half++）。
