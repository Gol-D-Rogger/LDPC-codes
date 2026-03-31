# IBEX Matrix LUT Semantics

本文档统一描述 IBEX QC-LDPC 矩阵 RTL 参数的生成语义。

当前脚本体系采用两阶段管线，完全对齐 `IBEX/src/ldpc_codec.cpp` 的 RTL 语义。所有公式均已在 `origin_matrix.csv` 上完成 13 行 × 4 档 K 的数值验证。

---

## 1. 目标与范围

当前目标是基于 `IBEX/ibex_matrix/{M}x{N}` 家族矩阵生成参数表：

- `M = 5..17`（parity 行数 / 码率配置数）
- `K = 64..67`（payload 列数，即 Extra User Data Column 的变化范围）
- `N = M + K`（总列数）
- `Z = 512`（circulant size）

生成的三个 RTL 参数：

| 参数 | 位宽 | 索引 | 语义 |
|------|------|------|------|
| `WrapBase[r]` | 9-bit | 每行一个 | Delta 链虚拟起点（外部输入或自动搜索） |
| `WrapBaseDelta[r][(M,K)]` | 6-bit | 每行每(M,K)配置一个 | 从 last_element 环绕到 0 所需的 delta 步数 |
| `FirstMaskShift[M][K]` | 9-bit | 每个 M 的最后一行，每档 K 一个 | 非整列 parity 掩码窗口偏移 |

---

## 2. 源矩阵与 (M, K) 裁剪规则

### 2.1 源矩阵

对每个 `M`，权威输入是 `K=67` 的 family 矩阵：

```text
IBEX/ibex_matrix/{M}x{M+67}/matrix/LDPC_{M}x{M+67}ex512_w4_dense5_QC_H_1_1.txt
```

Stage 1 脚本（`generate_family_artifacts.py`）将所有 M 的矩阵归一化为统一的 `ldpc_matrix_assignments_maxk67.svh`——以 packed column bitmaps 形式存储每列每行的 occupied/fade 状态，Stage 2 脚本从此文件读取。

### 2.2 裁剪规则

目标 `(M, K)` 矩阵由 `M × (M+67)` 矩阵裁剪得到：

- 保留 payload 列 `0..K-1`
- 保留最右侧 `M` 个 parity 列
- 删除 payload 尾部列 `K..66`

等价的 C model 列映射（`ldpc.h:573-580`）：

```c
if (j < (h_matrix.cols - h_matrix.rows))
    src_col = j;
else
    src_col = N - h_matrix.cols + j;
```

---

## 3. 基础不变量

### 3.1 Delta：行级固定步长

每行的所有 active CPM 的 shift 值构成一条模 Z 的等差数列，步长为 `delta[r]`。在 shift 链上，从左到右每命中一个 active CPM，shift 值增加一个 delta；反方向则减少一个 delta。

Delta 值（13 行 / 17 行）：

```text
delta_13 = {0, 13, 19, 29, 41, 67, 73, 79, 91, 97, 103, 111, 119}
delta_17 = {0, 13, 19, 29, 41, 67, 73, 79, 91, 97, 103, 111, 119, 127, 131, 137, 149}
```

Delta 是矩阵的结构常量，不随 `(M, K)` 配置改变。

### 3.2 First / Last Element

对任意 `(M, K)` 矩阵、任意行 `r`：

- `first_element` = 该行最左侧 active CPM 的 shift 值
- `last_element` = 该行最右侧 active CPM 的 shift 值
- `T` = 该行 active CPM 总数（occupied + fade）

三者满足：

```text
first_element = (last_element - (T - 1) × delta) mod Z
last_element  = (first_element + (T - 1) × delta) mod Z
```

### 3.3 Last Element 的固定常量性质

`last_element` 是每行最右侧 active CPM（位于 parity 区域）的 shift 值。由于 parity 区在裁剪过程中始终保留，`last_element` 不随 K 变化，仅由 parity 子矩阵的结构决定。IBEX 矩阵的 parity 结构由 G 子矩阵（rows 0..4）和 E 子矩阵（rows 5+）两部分组成，这两部分的 shift 赋值规则是完全确定性的，可以从 delta 常量直接推导。

**G 子矩阵（rows 0..g-1，g=5）** 占据最右侧 5 个 parity 列。矩阵生成代码（`GenLDPC_T4/src/gen_ldpc.c`）按如下规则为 G 赋值：对于行 i，在 5 个 parity 列中，满足 `(k % 5) != i`（其中 k 是从右端数的列偏移）的位置被标记为 occupied，且 row 4 的最后一列额外排除以保证矩阵可逆。每个 occupied 位置的 shift 按 `g_ind * delta[i]` 递增赋值，其中 g_ind 从 1 开始随列从左到右递增。因此，每行最右侧 G 条目的 shift 值等于该行 occupied 条目总数乘以 delta。

**E 子矩阵（rows g..M-1）** 是对角线结构，每行恰好有一个 parity 条目，其 shift 固定为 0（identity circulant）。

由此得到 `last_element` 的封闭公式：

```text
last_element[r] =
  (g - 1) × delta[r] mod Z    当 r < g - 1 时  (rows 0..3, 4个active parity条目)
  (g - 2) × delta[r] mod Z    当 r = g - 1 时  (row 4, 3个条目, 一个被排除)
  0                            当 r ≥ g 时      (E 对角线, identity shift)
```

其中 g = 5。展开为具体数值：

| row | delta | G/E | active parity 数 | last_element |
|-----|-------|-----|-------------------|-------------|
| 0 | 0 | G | 4 | 0 |
| 1 | 13 | G | 4 | 52 |
| 2 | 19 | G | 4 | 76 |
| 3 | 29 | G | 4 | 116 |
| 4 | 41 | G | 3 | 123 |
| 5 | 67 | E | 1 | 0 |
| 6 | 73 | E | 1 | 0 |
| ... | ... | E | 1 | 0 |
| 16 | 149 | E | 1 | 0 |

### 3.4 行平移不变性——为什么 WrapBase 可以自由选择

QC-LDPC 码的每个非零位置是一个 Z×Z 循环置换矩阵（CPM），其"shift"值决定了对角线的旋转量。一行中的所有 CPM 如果同时旋转相同的量——比如每个 shift 加上常数 c——效果等价于将该行对应的 Z 个校验方程集体做循环移位。由于校验矩阵的零空间（码空间）不受行排列影响，这种整行平移不改变码的纠错能力。

用一个具体例子来说明：假设某行有 3 个 active CPM，shift 分别为 `{100, 141, 182}`（delta=41）。如果把整行加上 50，变成 `{150, 191, 232}`，新矩阵与原矩阵的码空间完全相同。相邻 CPM 的间距始终是 41——这才是影响性能的关键。

```text
原始:      ···  100  ···  141  ···  182  ···
                 ↓ +41      ↓ +41
平移+50:   ···  150  ···  191  ···  232  ···
                 ↓ +41      ↓ +41

两行产生相同的码——仅 delta 间距不变即可
```

这一性质是 RTL 参数设计的理论基础：WrapBase（delta 链的起点）可以作为外部输入自由选择，或通过搜索优化以减小 WrapBaseDelta 的位宽。

---

## 4. Delta 链公式（核心恒等式）

### 4.1 链生成公式（RTL 语义：从左到右）

RTL 中每行的 shift 链按**从左到右**的方向生成，起点由 WrapBase 决定：

```text
first_active_shift = (WrapBase[r] + delta[r]) mod Z
next_active_shift  = (current_shift + delta[r]) mod Z
```

对应 `ldpc_codec.cpp` RTL 路径以及 `generate_family_artifacts_from_occ_fade_rtl.py` 中的 `rebuild_shifts_from_wrap_base()` 函数。

等价地，可以用位置索引 n 表示整条 shift 链：

```text
S_r(n) = (WrapBase[r] + (n + 1) × delta[r]) mod Z,   n = 0, 1, 2, ..., T-1
```

- `S_r(0)` = `first_element`（最左侧 active CPM 的 shift）
- `S_r(T-1)` = `last_element`（最右侧 active CPM 的 shift）

### 4.2 WrapBase 的语义

`WrapBase[r]` 是 delta 链上 `first_element` **前一步**的虚拟位置，即 `S_r(-1)`。它不对应矩阵中的任何实际列，而是作为 shift 链的"种子"——RTL 从 `(WrapBase + delta)` 开始分配 shift。

与 `first_element` 的关系：

```text
WrapBase[r] = (first_element_r - delta[r] + Z) mod Z
first_element_r = (WrapBase[r] + delta[r]) mod Z
```

**关键区分**：WrapBase 在当前的脚本体系中有两种来源——

- **固定输入模式**（`generate_family_artifacts_from_occ_fade_rtl.py`）：调用者直接提供 WrapBase 表，脚本不做搜索。默认值为 13-row RTL 参考值 `{0, 270, 415, 337, 182, 301, 445, 70, 490, 65, 85, 75, 216}`。
- **自动搜索模式**（`generate_family_artifacts_from_occ_fade_rtl_autobase.py`）：对每行在 [0, Z) 范围内暴力搜索，选择使 WrapBaseDelta 位宽最小的 WrapBase。

两种模式的下游参数计算逻辑完全相同——差异仅在于 WrapBase 的来源。

### 4.3 从左到右 vs 从右到左

Shift 链可以从两个方向重建，结果完全等价：

| 方向 | 起点 | 步进 | RTL/脚本 |
|------|------|------|---------|
| **从左到右** | `WrapBase + delta` | `+delta` | `rebuild_shifts_from_wrap_base()` — RTL 对齐 |
| **从右到左** | `last_element` | `-delta` | `regenerate_shifts()` — 早期常量路径 |

当 WrapBase 取 `(first_element - delta) % Z` 时，两个方向从同一条等差数列的两端出发，产生完全相同的 shift 赋值。RTL 脚本选择从左到右以直接对齐硬件行为。

### 4.4 全行验证

`S_r(0) = first_element` 在 `origin_matrix.csv` 全部 13 行验证通过：

| row | delta | wrap_base | S_r(0) | origin first | status |
|-----|-------|-----------|--------|--------------|--------|
| 0 | 0 | 0 | 0 | 0 | OK |
| 1 | 13 | 270 | 283 | 283 | OK |
| 2 | 19 | 415 | 434 | 434 | OK |
| 3 | 29 | 337 | 366 | 366 | OK |
| 4 | 41 | 182 | 223 | 223 | OK |
| 5 | 67 | 301 | 368 | 368 | OK |
| 6 | 73 | 445 | 6 | 6 | OK |
| 7 | 79 | 70 | 149 | 149 | OK |
| 8 | 91 | 490 | 69 | 69 | OK |
| 9 | 97 | 65 | 162 | 162 | OK |
| 10 | 103 | 85 | 188 | 188 | OK |
| 11 | 111 | 75 | 186 | 186 | OK |
| 12 | 119 | 216 | 335 | 335 | OK |

---

## 5. WrapBaseDelta（Definition A）

### 5.1 定义

WrapBaseDelta 回答一个问题：从 `last_element` 出发，沿着 delta 链继续步进，需要多少步才能环绕回 0？

```text
WrapBaseDelta[r][(M,K)] = min k ≥ 0  such that  (last_element + k × delta) mod Z == 0
```

这是 RTL 使用的**定义 A**（数学定义），直接对应 `compute_wrap_base_delta_from_last_element()` 函数。

### 5.2 直觉理解——环形跑道上的"余程"

想象一条周长为 Z=512 的环形跑道。选手从 `(WrapBase + delta)` 处起跑，每到一个 active CPM 标记处就前进 delta 步。经过 T 个标记后到达 `last_element`——这是比赛的"正式终点"。

WrapBaseDelta 回答的问题是：从终点继续跑，还差多少步才能回到跑道的 0 标记？

```text
          0                       Z = 512 的环形跑道
         ╱ ╲
        ╱   ╲
       │     │   ← last_element (终点)
       │     │
       │     │
        ╲   ╱   ← first_element (起跑线)
         ╲ ╱
          ↓
   展开为线性视图:

   WrapBase  →  FE  → ... →  LE  → (+k·delta) →  0(mod Z)
   (虚拟起点) S(0)   T个CPM  S(T-1)   k 步余程     回到零
                                    = WrapBaseDelta
```

为什么 K 变小时 WrapBaseDelta 变大？K 减小意味着 payload 列被裁掉，该行的 active CPM 数 T 减少。Shift 链提前结束，last_element 离 0 更远，需要更多步才能回到 0。

### 5.3 Wraparound 恒等式

环绕距离 `wraparound` 定义为 shift 链从 `first_element` 越过 `last_element` 后、继续步进直到回到 `first_element` 的距离：

```text
wraparound = (Z + first_element - last_element) mod Z
```

Decoder 的旋转量等于 `WrapBase + WrapBaseDelta × delta`，与 wraparound 之间的恒等关系：

```text
decoder_rotation = WrapBase + WrapBaseDelta × delta
                 = wraparound - delta
```

在 `origin_matrix.csv`（M=13, K=67）上以 Definition A 验证：

| row | delta | last_element | WrapBaseDelta | `(LE + wbd·δ) % Z` | status |
|-----|-------|-------------|---------------|---------------------|--------|
| 0 | 0 | 0 | 0 (delta=0 特殊情况) | 0 | OK |
| 1 | 13 | 52 | 32 | 0 (52+32×13=468≡0) | OK |
| 2 | 19 | 76 | 33 | 0 | OK |
| 3 | 29 | 116 | 33 | 0 | OK |
| 4 | 41 | 123 | 32 | 0 | OK |
| 5 | 67 | 0 | 0 | 0 | OK |
| 6 | 73 | 0 | 0 | 0 | OK |
| 7 | 79 | 0 | 0 | 0 | OK |
| 8 | 91 | 0 | 0 | 0 | OK |
| 9 | 97 | 0 | 0 | 0 | OK |
| 10 | 103 | 0 | 0 | 0 | OK |
| 11 | 111 | 0 | 0 | 0 | OK |
| 12 | 119 | 0 | 0 | OK |

> **注**：rows 5-12 的 last_element = 0（E 子矩阵 identity shift），Definition A 直接得 k=0。row 0 的 delta=0，Definition A 也直接得 k=0。RTL DV Verilog 中对应位置的 `wrap_num_deltas` 值（如 row 0 的 31）来自早期 DV 脚本的**行权重差分定义**，两者数值不同，但在 RTL 对齐语义下以 Definition A 为准。

### 5.4 与行权重差分定义的关系

早期 DV 脚本使用另一种等价定义：`wbd = ref_row_weight - target_row_weight`。当所有 (M,K) 配置共享同一个 WrapBase 和 last_element 时，减少一个 active CPM 等价于 delta 链缩短一步，使得到达 0 需要多一步。两个定义在特定条件下数值相同，但 RTL 脚本直接使用 Definition A（数学定义）以避免依赖参考矩阵的选择。

---

## 6. FirstMaskShift

### 6.1 定义

```text
FirstMaskShift[M][K] = (Z + FE_last_row - first_parity_CPM_last_row) mod Z
```

其中：
- `FE_last_row` = M 矩阵最后一行（row=M-1）的 first_element（整行最左侧 active shift）
- `first_parity_CPM_last_row` = 该行在 parity 区域（col ≥ K）的第一个 active CPM 的 shift

### 6.2 直觉理解——Mask 窗口的锚点

当 parity 长度不是 Z 的整数倍时（`extra_bits_of_parity > 0`），最后一行的 syndrome 向量只有前 `extra_bits_of_parity` 个 lane 有效，其余 lane 是 padding。RTL 需要知道这个有效窗口从哪里开始——这就是 FirstMaskShift。

为什么不是简单地"从 0 开始"？因为每个 CPM 都带有旋转（shift），有效位的起始位置取决于该 CPM 的 shift 值。RTL 把整行看作一条 shift 链，FirstMaskShift 告诉它：从 shift 链的起点（first_element）到 parity 区的入口（first_parity_CPM）之间隔了多远。

```text
  payload 区                    parity 区
  ┌───────────────────────┐   ┌──────────────────┐
  │ FE ··· (active CPMs)  │   │ parity_1st ··· LE│
  └───────────────────────┘   └──────────────────┘
       ←  FirstMaskShift  →

  FirstMaskShift = (Z + FE - parity_1st) mod Z
                 = payload 区 active CPM 贡献的 delta 步数 × delta
```

如果某个 (M,K) 配置的最后一行在 payload 区没有 active CPM，那么 FE 就是 parity 区的第一个 CPM，FirstMaskShift = 0。

### 6.3 K 变化时的 Delta 步进

当 K 从 67 逐步减小到 64 时，payload 区最右侧的列依次被裁掉。如果被裁列在最后一行有 active CPM，则 T 减 1，`first_element` 沿 delta 链向右移动一步（增加 delta），FirstMaskShift 随之增加 delta。如果被裁列无 active CPM，则 FirstMaskShift 不变。

以 M=5（row=4, delta=41）为例：

```text
K=64: FMS=510    K=65: FMS=469    K=66: FMS=428    K=67: FMS=387
           -41           -41           -41
```

每档差值 = delta[4] = 41，说明 col 64, 65, 66 在 row 4 均有 active CPM。

以 M=7（row=6, delta=73）为例：

```text
K=64: FMS=225    K=65: FMS=225    K=66: FMS=152    K=67: FMS=79
            0            -73           -73
```

K=64 与 K=65 的 FMS 相同（差值=0），说明 col 64 在 row 6 无 active CPM。

### 6.4 DV Verilog 参考值

```text
mask_shift_05 = {510, 469, 428, 387}   (K=64..67)
mask_shift_06 = { 57, 502, 435, 435}
mask_shift_07 = {225, 225, 152,  79}
mask_shift_08 = {307, 228, 228, 228}
mask_shift_09 = {251, 251, 160, 160}
mask_shift_10 = {453, 356, 259, 259}
mask_shift_11 = { 88, 497, 394, 291}
mask_shift_12 = {408, 408, 408, 297}
```

---

## 7. 两阶段管线与完整推理链

### 7.1 架构总览

```text
Stage 1: generate_family_artifacts.py
┌──────────────────────────────────────────────────────┐
│  输入: ibex_matrix/{M}x{N}/ (matrix/occupied/fade)   │
│  输出: ldpc_matrix_assignments_maxk67.svh            │
│        (统一的 occupied/fade bitmaps)                 │
└──────────────────────┬───────────────────────────────┘
                       │
                       ▼
Stage 2: generate_family_artifacts_from_occ_fade_rtl.py  (固定 WrapBase)
   或者: generate_family_artifacts_from_occ_fade_rtl_autobase.py (自动搜索)
┌──────────────────────────────────────────────────────┐
│  输入: assignments SVH + WrapBase 表                  │
│                                                      │
│  处理: 对每个 (M, K):                                 │
│    1. 裁剪 occupied/fade → 目标视图                   │
│    2. 从 WrapBase 左到右重建 shift 链                  │
│    3. 提取 first_element, last_element, row_weight    │
│                                                      │
│  输出: WrapBase / WrapBaseDelta / FirstMaskShift      │
│        (SVH + JSON + CSV)                            │
└──────────────────────────────────────────────────────┘
```

### 7.2 RTL 对齐的 Shift 重建（核心算法）

Stage 2 的核心是 `rebuild_shifts_from_wrap_base()` 函数，精确镜像 `ldpc_codec.cpp::ldpc_config()` 的 RTL 语义：

```text
对每行 r:
  cursor = (WrapBase[r] + delta[r]) mod Z          // 第一个 shift
  从 col = 0 到 col = N-1（从左到右）:
    如果 occupied[r][col] 或 fade[r][col]:
      shift[r][col] = cursor
      如果是第一个 active → 记录 first_element
      如果在 parity 区且是该区第一个 → 记录 first_parity_shift
      更新 last_element = cursor
      cursor = (cursor + delta[r]) mod Z            // 下一个 shift
```

这条路径完全不依赖 shift matrix 文件——仅需 occupied/fade 二值模式 + WrapBase + delta 常量。

### 7.3 三个参数的计算

从重建的 shift 视图中提取：

```text
┌─────────────────┐   ┌──────────────────────────┐
│  WrapBase[r]    │   │  Assignments SVH          │
│  (固定输入或搜索) │   │  (occupied/fade bitmaps)  │
└────────┬────────┘   └────────────┬──────────────┘
         │                         │
         └────────┬────────────────┘
                  ▼
     ┌─────────────────────────────┐
     │  rebuild_shifts_from_       │
     │  wrap_base()                │
     │  从左到右: cursor += delta  │
     │  → FE, LE, T, row_weight,  │
     │    first_parity_shift       │
     └──────────────┬──────────────┘
                    │
      ┌─────────────┼────────────────┐
      ▼             ▼                ▼
 ┌──────────┐  ┌───────────────┐  ┌────────────────┐
 │WrapBase  │  │WrapBaseDelta  │  │FirstMaskShift  │
 │(透传)    │  │= min k:       │  │= (Z + FE       │
 │          │  │  LE+k·δ≡0    │  │  - parity_1st) │
 │          │  │  (mod Z)      │  │  mod Z         │
 └──────────┘  └───────────────┘  └────────────────┘
```

### 7.4 RTL shift_dist 搜索环路

`ldpc_codec.cpp::ldpc_config()` 中的 RTL 搜索并不直接计算 Definition A，而是用一个双层循环同时搜索 `shift_dist`（即整行的额外偏移量 j）和 `delta_num`（即环绕步数 k）：

```c
init_base = (first_element - delta + bits) % bits;   // = WrapBase
for (j = 0; j < bits; j++) {                          // 外层: 尝试 shift_dist
    shift_base = init_base + j;
    for (k = 0; k < 50; k++) {                        // 内层: 尝试 delta_num
        if ((last_element + shift_base + k*delta + j) % bits == shift_base) {
            // 找到了: shift_dist = j, wrap_num_deltas = k
            break;
        }
    }
}
```

当 j=0 时，条件简化为 `(last_element + k*delta) % Z == 0`，这恰好就是 Definition A。因此对于大多数行，RTL 搜索在 j=0 就找到结果，`shift_dist = 0`。

`shift_dist` 的物理含义是：如果 delta 链的自然起点无法让环绕恰好对齐，需要给整行施加多少额外偏移。它是 autobase 脚本（`generate_family_artifacts_from_occ_fade_rtl_autobase.py`）中 `compute_rtl_shift_dist_from_ldpc_config()` 函数的输出之一，用于 debug 验证。

### 7.5 AutoBase 搜索策略

`generate_family_artifacts_from_occ_fade_rtl_autobase.py` 对每行在 [0, Z) 暴力搜索 WrapBase，评分标准（按优先级）：

1. 所有 (M,K) 的 WrapBaseDelta ≤ width_limit（默认 63）
2. 族内最大 WrapBaseDelta 最小
3. 族内 WrapBaseDelta 总和最小
4. WrapBase 值本身最小（tie-break）

搜索利用了 Section 3.4 的行平移不变性：改变 WrapBase 等价于对整行施加统一偏移，不影响编解码正确性，但会改变 last_element 的值，从而影响 WrapBaseDelta 的大小。

---

## 8. IBEX Family 的实现细节

### 8.1 多 M 矩阵的 WrapBase 不一致问题

IBEX family 的 13 个（或 17 个）矩阵由外部工具独立搜索生成，同一行的 delta 链起始点在不同 M 矩阵中可能不同。这意味着对于 13-rate family，不存在一个统一的 WrapBase 能让所有 M 矩阵的 shift 链同时对齐。

数值验证结果（`(first_element - delta) mod Z` 跨 M 的一致性）：

| row | 不同值数 | 说明 |
|-----|---------|------|
| 0 | 1 | 一致（delta=0） |
| 1 | 12 | M=5: 322, M=6: 322, M=7: 439, ..., M=17: 291 |
| 2-4 | 13 | 全部不同 |
| 5-14 | 2-12 | 部分不同 |
| 15-16 | 1 | 一致（仅 1-2 个 M 包含该行） |

### 8.2 解决方案

对于 origin_matrix.csv（单一 13×80 矩阵），WrapBase 是单一确定的值。对于多 M family，有两种策略：

- **固定 WrapBase**：选定一个参考矩阵（如 M=13, K=67），从其 first_element 推导 WrapBase。其他 M 矩阵的 shift 链会相应偏移，但 WrapBaseDelta（Definition A）仍然正确。
- **自动搜索 WrapBase**：per-row 搜索使 WrapBaseDelta 位宽最小的值。

### 8.3 参数存储

| 参数 | 维度 | 位宽 | 索引方式 |
|------|------|------|---------|
| `delta[r]` | 13 或 17 | 9-bit | row |
| `WrapBase[r]` | 13 或 17 | 9-bit | row |
| `WrapBaseDelta[r][(M,K)]` | 三角结构 | 6-bit | (row, M, K) |
| `FirstMaskShift[M][K]` | M×4 | 9-bit | (M, K) |

### 8.4 生成脚本

**Stage 1**: `scripts/ibex_matrix_family/generate_family_artifacts.py`

- 输入：`ibex_matrix/{M}x{N}/` 下的 matrix/occupied/fade 文件
- 输出：`ldpc_matrix_assignments_maxk67.svh` + ENS/RDEC 调度表

**Stage 2a**: `scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py`

- 输入：`assignments_maxk67.svh` + 固定 WrapBase 表（`--wrap-base`）
- 输出：SVH / JSON / CSV / 重建矩阵

**Stage 2b**: `scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl_autobase.py`

- 输入：`assignments_maxk67.svh`
- 处理：per-row 暴力搜索最优 WrapBase
- 额外输出：`rtl_shift_dist` 表（对齐 `ldpc_codec.cpp` 的 shift_dist 参数）

---

## 9. Origin Matrix 的角色

`IBEX/ibex_matrix/origin_matrix.csv` 是旧 13×80 矩阵的 shift 表。

它的角色：

- **已完成**：用于验证 delta 链公式、wraparound 恒等式、RTL Verilog 参数一致性
- **不适用于**：直接复用到 M=5..17 family 的参数生成

对新 family，参数必须从 `IBEX/ibex_matrix/{M}x{N}` 矩阵重新提取。但公式体系（delta 链、wraparound 恒等式、Definition A）在所有 family 上通用。

---

## 附录 A. 速查

| 参数 | 定义 | 来源 |
|------|------|------|
| `delta[r]` | 行级固定步长 | 固定常量 `(0,13,19,...,149)` |
| `last_element[r]` | 最右 active parity CPM 的 shift | 固定常量：rows 0..3 = `4×δ`, row 4 = `3×δ`, rows 5+ = `0` |
| `WrapBase[r]` | delta 链虚拟起点 | 固定输入或自动搜索 |
| `WrapBaseDelta[r][(M,K)]` | min k: `(LE + k·δ) ≡ 0 (mod Z)` | Definition A，从 last_element 反推 |
| `FirstMaskShift[M][K]` | `(Z + FE - first_parity_CPM) mod Z` | 最后一行 FE 与 parity 首 CPM 之差 |

核心恒等式：

```text
S_r(n) = (WrapBase[r] + (n+1) × delta[r]) mod Z      // RTL: 从左到右

S_r(0)   = first_element  = (WrapBase + delta) mod Z
S_r(T-1) = last_element

WrapBaseDelta = min k ≥ 0: (last_element + k·delta) ≡ 0 (mod Z)

decoder_rotation = WrapBase + WrapBaseDelta × delta
```

---

## 10. CN 更新中 WrapBase / WrapBaseDelta 的作用

本节从 syndrome 数学定义出发，推导 RTL 逐 cycle 处理的移位机制，解释 WrapBase 和 WrapBaseDelta 如何让 cn 寄存器在迭代间正确对齐。

### 10.1 Syndrome 的矩阵定义与 delta 等差结构

**符号约定**：`r` = 行索引，`c` = 列索引，`Z` = 512（circulant size），`>>` = 右循环移位，`<<` = 左循环移位。

对行 r，syndrome 向量 Sr 是所有列的 CPM 贡献之和：

```text
Sr = Hr_0 * C0 + Hr_1 * C1 + ... + Hr_(N-1) * C(N-1)
```

其中 `Hr_c` 是第 c 列的 CPM 矩阵，`Cc` 是第 c 列的 codeword 向量。非零 CPM 的乘法等价于对 `Cc` 做循环移位（右移 shift 值）。

取出行 r 所有 active（非零 circulant）列的 codeword，按从左到右顺序记为 `C' = {C'0, C'1, ..., C'(T-1)}`（T = 该行 active CPM 总数）。由 Section 3.1 的 delta 等差性质，这些 CPM 的 shift 值构成模 Z 等差数列：

```text
shift 序列（从左到右）:
  {j·δ, (j-1)·δ, ..., δ, 0, -δ, -2δ, ..., -x_r·δ}   (mod Z)
```

其中：
- `j` = parity 侧 active CPM 数（shift > 0 的部分，位于链头方向）
- `x_r` = payload 侧 active CPM 数（shift < 0 的部分，位于链尾方向）
- `T = j + x_r + 1`（包含 shift=0 的位置）
- `δ = delta[r]`

> 注：这里的 "shift=0 位置" 是 delta 链上的一个参考锚点，不一定是矩阵中的物理零。实际 shift 值由 WrapBase 平移后确定（见 Section 3.4 行平移不变性）。此处为简化推导，以 shift=0 为锚点展开。

Syndrome 展开为循环移位之和：

```text
Sr = C'0 >> (j·δ) + C'1 >> ((j-1)·δ) + ... + C'(j) >> 0 + ... + C'(T-1) >> (-x_r·δ)
```

### 10.2 逐 Cycle 处理：cn 的右移调度

RTL 每个 cycle 处理一列 active CPM。cn 寄存器存储的是 syndrome 的一个旋转版本——每 cycle 右移 delta，使当前列的 CPM 对齐到固定位置，从而只需一个固定位置的 XOR 操作。

**Cycle 0**（处理最左侧 active 列 C'0）：

cn 初始化为 Sr 左移 `j·δ`，使 C'0 的贡献对齐到位置 0：

```text
S'r0 = Sr << (j·δ)
     = C'0 + C'1 << δ + ... + C'(x_r) << (j+x_r)·δ
```

**Cycle 1**（处理 C'1）：

cn 右移 delta，使 C'1 对齐到位置 0：

```text
S'r1 = S'r0 >> δ = Sr << ((j-1)·δ)
     = C'0 >> δ + C'1 + C'2 << δ + ... + C'(x_r) << (j+x_r-1)·δ
```

**Cycle k**（处理 C'k）：

```text
S'rk = S'r0 >> (k·δ) = Sr << ((j-k)·δ)
```

**Cycle (j+x_r)**（处理最后一列 C'(T-1)）：

```text
S'r(j+x_r) = Sr << (-x_r·δ)
```

**关键恒等式**：最后一个 cycle 的 cn 状态与第 0 个 cycle 的关系：

```text
S'r(j+x_r) << (j+x_r)·δ = S'r0
```

即：经过 T-1 = j+x_r 个 cycle 的右移后，cn 需要**左移 `(j+x_r)·δ`** 才能回到 cycle 0 的状态。

### 10.3 迭代间回绕：基本方程的代数推导

一轮迭代处理完所有 active 列后，cn 停在 `S'r(j+x_r)` 状态。下一轮迭代需要回到 `S'r0`。

由 10.2 的恒等式，回到 cycle 0 需要左移 `(T-1)·δ`（其中 T = j + x_r + 1 为 active CPM 总数）。等价地，需要右移 `h·Z - (T-1)·δ`（h 为使结果为正的最小正整数）。

RTL 将这个右移量表示为 `(y_r + 1)·δ + base_r`：

```text
(y_r + 1) · δ + base_r = h · Z - (T - 1) · δ
```

移项得 **基本方程**：

```text
(T + y_r) · δ + base_r = h · Z      (h 为正整数)
```

其中 `y_r = WrapBaseDelta`，`base_r = WrapBase`。

**代数证明**：基本方程实际上是 Definition A 与 last_element 公式的直接推论——

1. 由 Section 4.1 的 shift 链公式：`last_element = (WrapBase + T · δ) mod Z`
2. 由 Definition A（Section 5.1）：`(last_element + WBD · δ) mod Z = 0`
3. 将 1 代入 2：`(WrapBase + T · δ + WBD · δ) mod Z = 0`
4. 整理：`(T + WBD) · δ + WrapBase ≡ 0 (mod Z)`

这是一个恒等式——对任何满足 Definition A 的 (WrapBase, WBD, T, delta) 组合都成立。

直觉理解：delta 链上共 T 个 active 位置消耗 T 步 delta，从 last_element 走回 0 需要 WBD 步 delta，最后加上 WrapBase 的余量，总计恰好走完 h 整圈。

```text
环形跑道展开（周长 Z=512）：

处理阶段 (T 步 delta)      回绕阶段 (WBD 步 delta + WrapBase)
├────────────────────────┤├──────────────────────────────────┤
S_r(0) → ... → S_r(T-1)    (+δ)×WBD 步          (+WrapBase)
  FE              LE        LE→0 (Def A)         0→WrapBase

总计: (T + WBD)·δ + WrapBase ≡ 0  (mod Z)
      ─────────────────────────────────────
      即恰好走完 h 整圈
```

### 10.4 基本方程的数值验证

用 origin_matrix（M=13, K=67, Z=512）的 13-row RTL 参考 WrapBase 验证。

首先确定每行的 T（active CPM 总数）。由 `last_element = (WrapBase + T·δ) mod Z` 反解 T：

```text
T = [(last_element - WrapBase) mod Z] × δ⁻¹ mod Z
```

（`δ⁻¹` 为 delta 在模 Z 下的乘法逆元；delta=0 时 T 无关紧要。）

**完整验证表**：

| row | δ | WrapBase | LE | T | WBD(DefA) | (T+WBD)·δ+WB | h | status |
|-----|-----|---------|-----|-----|-----------|---------------|---|--------|
| 0 | 0 | 0 | 0 | - | 0 | 0 | 0 | OK (trivial) |
| 1 | 13 | 270 | 52 | 62 | 508 | 7680 | 15 | OK |
| 2 | 19 | 415 | 76 | 63 | 508 | 11264 | 22 | OK |
| 3 | 29 | 337 | 116 | 63 | 508 | 16896 | 33 | OK |
| 4 | 41 | 182 | 123 | 61 | 509 | 23552 | 46 | OK |
| 5 | 67 | 301 | 0 | 49 | 0 | 3584 | 7 | OK |
| 6 | 73 | 445 | 0 | 43 | 0 | 3584 | 7 | OK |
| 7 | 79 | 70 | 0 | 38 | 0 | 3072 | 6 | OK |
| 8 | 91 | 490 | 0 | 34 | 0 | 3584 | 7 | OK |
| 9 | 97 | 65 | 0 | 31 | 0 | 3072 | 6 | OK |
| 10 | 103 | 85 | 0 | 29 | 0 | 3072 | 6 | OK |
| 11 | 111 | 75 | 0 | 27 | 0 | 3072 | 6 | OK |
| 12 | 119 | 216 | 0 | 24 | 0 | 3072 | 6 | OK |

逐行抽检：

- **Row 4**（G 子矩阵）：`(61 + 509) × 41 + 182 = 570 × 41 + 182 = 23370 + 182 = 23552 = 46 × 512` ✓
- **Row 5**（E 子矩阵）：`(49 + 0) × 67 + 301 = 3283 + 301 = 3584 = 7 × 512` ✓
- **Row 1**（G 子矩阵）：`(62 + 508) × 13 + 270 = 570 × 13 + 270 = 7410 + 270 = 7680 = 15 × 512` ✓

> **关于 WBD 的数量级**：注意到 G 行（rows 1-4）的 Definition A 值为 508/509，远超 6-bit 范围（max 63）。这是因为 `(last_element + k·δ) mod Z = 0` 在 Z=512 且 gcd(δ, Z) = 1 时，解 k 在 [0, Z) 范围内是唯一的——如果 last_element 不"靠近" 0（mod δ），k 就会很大。E 行因 last_element=0 而 WBD=0。RTL 通过 **shift_dist 机制**（见 10.5）解决 G 行的位宽溢出问题。

### 10.5 RTL shift_dist 机制：调整 WrapBase 以压缩 WBD

`ldpc_codec.cpp::ldpc_config()` 中的 RTL 搜索（Section 7.4）并不直接使用 Definition A 的 WBD，而是引入 `shift_dist`（记为 j）允许对 WrapBase 施加额外偏移：

```text
WrapBase_actual = init_base + j     (init_base = FE - δ = 原始 WrapBase)
```

搜索条件变为：

```text
(last_element + j + k · δ) mod Z = 0,    j ≥ 0, k < 50
```

当 j > 0 时，等价于将整行平移 j（Section 3.4 行平移不变性），使 last_element 变为 `last_element + j`，从而找到更小的 k。

**Origin matrix RTL 搜索结果**：

| row | δ | shift_dist (j) | wrap_num_deltas (k) | WB_actual | 验证 |
|-----|-----|---------------|--------------------|-----------| -----|
| 0 | 0 | 0 | 0 | 0 | trivial |
| 1 | 13 | 5 | 35 | 275 | (52+5+35×13)%512=0 ✓ |
| 2 | 19 | 17 | 49 | 432 | (76+17+49×19)%512=0 ✓ |
| 3 | 29 | 9 | 31 | 346 | (116+9+31×29)%512=0 ✓ |
| 4 | 41 | 19 | 34 | 201 | (123+19+34×41)%512=0 ✓ |
| 5-12 | >0 | 0 | 0 | =init_base | LE=0, trivial ✓ |

G 行（1-4）的 shift_dist > 0，通过平移 WrapBase 将 WBD 压缩到 k < 50。E 行的 LE=0，j=0 即可满足。

基本方程对 shift_dist 调整后的值同样成立：

```text
Row 1: (T + k)·δ + WB_actual = (62 + 35)×13 + 275 = 1261 + 275 = 1536 = 3×512 ✓
Row 4: (T + k)·δ + WB_actual = (61 + 34)×41 + 201 = 3895 + 201 = 4096 = 8×512 ✓
```

### 10.6 DV 推导的核心恒等式

DV 团队的推导从 syndrome 逐 cycle 展开出发，得到的关键结论是：

```text
S'r(j+x_r) << (j+x_r)·δ = S'r0
```

即最后一个 cycle 的 cn 左移 `(T-1)·δ` 等于第一个 cycle 的 cn。

RTL 将回绕实现为右移 `(y_r + 1)·δ + base_r`，这个右移与上述左移互补：

```text
h · Z - [(y_r + 1) · δ + base_r] = (T - 1) · δ
⟹ (T + y_r) · δ + base_r = h · Z
```

这与 10.3 的基本方程完全一致，其中 DV 的记号对应关系为 `j + x_r + 1 = T`，`y_r = WBD`，`base_r = WrapBase`。

### 10.7 两种调度方式的对比

| 方面 | Version 1（逐列旋转） | Version 2（迭代间一次性旋转） |
|------|----------------------|------------------------------|
| cn 旋转时机 | 每列右移 delta | 一轮结束后一次性右移 |
| 处理阶段 cn 状态 | 每 cycle 更新 | 不变（toggle 写入 flip_ram） |
| 回绕机制 | 隐式（逐步累积已包含回绕） | 显式：`(WBD+1)·δ + WrapBase` |
| 需要 WrapBase/WBD | 不需要 | 需要 |
| RTL 对应 | 概念简单但移位器利用率高 | 流水线友好，移位器在迭代间使用 |

Version 2 对应流水线化的 RDEC 调度——内层循环只做 VN 决策和 flip_ram 写入，cn 的旋转延迟到迭代边界，用一个乘加运算完成。

### 10.8 为什么拆分为 WrapBase + WrapBaseDelta × delta

RTL 不直接存储完整旋转量（它随 M/K 配置变化），而是拆分为：

| 分量 | 位宽 | 随 (M,K) 变化？ | 说明 |
|------|------|-----------------|------|
| `WrapBase[r]` | 9-bit | 否（每行固定） | delta 链虚拟起点（含 shift_dist 调整），存一次 |
| `WrapBaseDelta[r][(M,K)]` | 6-bit | 是（每配置一个） | 环绕步数（RTL 搜索保证 < 50），存多个 |

乘法 `WrapBaseDelta × delta` 由 RTL 硬件在运行时完成。这样只需存 6-bit 的 WrapBaseDelta 表而不是 9-bit 的完整旋转量，节省 LUT 面积。

K 变小 → active CPM 数 T 减少 → 处理阶段消耗的 delta 步数减少 → 回绕需要更多步 → WrapBaseDelta 增大。这与 Section 5.2 的"余程"直觉一致。

> **关于 Section 5.3 验证表的勘误**：Section 5.3 表中 rows 1-4 的 WrapBaseDelta 值（32, 33, 33, 32）与 Definition A 的真实计算结果（508, 508, 508, 509）不符，也与 RTL 搜索结果（35, 49, 31, 34）不同。这些值可能来自早期 DV 脚本的行权重差分定义，已在 Section 5.4 中说明。本节的验证表（10.4）使用经过脚本验证的 Definition A 正确值。

---

## 11. WrapBase 变换对 Mask/Fade 机制的影响分析

行平移不变性（Section 3.4）说：行内所有 shift 同加常量 j 产生等价码。但当 `extra_bits_of_parity > 0` 时，最后一行只有 `L = extra_bits_of_parity` 个 CN 索引对应真实 parity bit，行平移不变性不完全成立。本节分析 WrapBase 变换（包括 RTL shift_dist 调整）对 mask/fade 机制的影响。

### 11.1 Mask/Fade 的两个独立约束域

系统通过两个独立机制分别约束 CN 域和 VN 域：

**CN 域 — mask/fade 分裂**（参见 `ibex_mask_mechanism.md` Section 2）：

```text
mask[col][k] = 1  iff  m = (k - element[last_row][col]) % Z < L
```

始终选中 CN 索引 {0, 1, ..., L-1} 为 occupied-last-row 有效区，{L, ..., Z-1} 为 fade 有效区。mask 与 !mask 互为补集，覆盖全部 Z 个位置。

**VN 域 — do_not_use_this_bit**（参见 `ibex_mask_mechanism.md` Section 6.1）：

```text
第一列 parity (col = n-m): k >= L → padding, 物理不存在, 永远为 0
```

这是绝对物理约束，不随 WrapBase 变化。

### 11.2 WrapBase 变化导致 Mask 窗口偏移

设 last row 在首列 parity 的 shift 原为 s，WrapBase 变化 +j 后变为 s+j。

```text
First parity column 的 k 域示意 (Z=16, L=6):

变化前 (s=0):
k:     0 1 2 3 4 5 | 6 7 8 9 A B C D E F
real:  R R R R R R | . . . . . . . . . .   ← 物理约束, 不变
mask:  M M M M M M | f f f f f f f f f f   ← M=occupied-lastrow, f=fade
       ^^^^^^^^^^^^   完美对齐

变化后 (s+j=3):
k:     0 1 2 3 4 5 | 6 7 8 9 A B C D E F
real:  R R R R R R | . . . . . . . . . .   ← 不变
mask:  f f f M M M | M M M f f f f f f f   ← 窗口右移 3
```

错位区域：
- k=0,1,2：真实 bit，但走了 fade 路径（不经过 occupied-last-row CN）
- k=6,7,8：padding bit，但走了 mask=1 路径（参与 occupied-last-row CN syndrome）

### 11.3 影响分析：什么会变，什么不变

#### 不变的性质

| 性质 | 原因 |
|------|------|
| **mask/fade 互补性** | `mask + !mask = Z` 始终成立，每个 VN bit 不遗漏不重复 |
| **编码/解码自洽性** | 编码器基于变化后的 H 矩阵计算 parity，syndrome 编码后为 0；解码器使用同一套 mask/fade 规则 |
| **padding bit 安全性** | do_not_use_this_bit 保证 padding VN bit 永远为 0，即使落入 mask=1 区域也不被翻转 |
| **非 last-row 的行** | 不涉及 mask，纯等价码变换 |

#### 变化的性质

**有效 H 矩阵改变**：mask 窗口偏移导致 last row 的 occupied/fade 连接模式变化：

```text
WrapBase +j
  → last row 所有 shift +j
  → mask 窗口在 k 域整体移动 j 位
  → first parity 列: 原先 mask=1 的某些真实 bit 变成 fade
                      原先 fade 的某些 padding bit 变成 mask=1
  → last row 的有效连接模式改变
  → 产生不同的码（不是等价码）
```

虽然 padding bit 不会被翻转（do_not_use_this_bit 兜底），但 last row CN 的 syndrome 中多了一些永远为 0 的常数贡献、少了一些原本应参与的真实连接。这改变了码的 girth、最小距离等特性，可能影响 FER 性能。

### 11.4 各场景的风险等级

| 场景 | 影响 | 风险 |
|------|------|------|
| 非 last-row 的 WrapBase 变化 | 纯等价码变换 | 无风险 |
| last row, `extra_bits_of_parity = 0` | mask 不参与, 全部 Z 个 CN 都 active | 无风险 |
| **last row, `extra_bits_of_parity > 0`, shift_dist = 0** | mask 窗口不移动, 保持原始对齐 | **无风险** |
| **last row, `extra_bits_of_parity > 0`, shift_dist > 0** | mask 窗口偏移, 产生不同的码 | **需 FER 验证** |

### 11.5 实践建议

1. **当 `extra_bits_of_parity = 0` 时**：所有行的 WrapBase/shift_dist 调整均安全，无需额外验证。

2. **当 `extra_bits_of_parity > 0` 时**：
   - **推荐**：约束 last row 的 `shift_dist = 0`（保持原始 WrapBase），仅对 rows 0..(M-2) 做 shift_dist 搜索。
   - **若 last row 必须使用 shift_dist > 0**：需通过 FER 仿真验证变化后的码性能，不能仅依赖行平移不变性假设。

3. **autobase 搜索的安全策略**：autobase 脚本搜索 WrapBase 时，可以将 last row 独立处理——先固定 last row 的 WrapBase 使 mask 对齐（shift_dist=0），再对其余行搜索最优 WrapBase。

---

## 12. 对齐版 C++ RTL 实现

本节说明 `DVCtrans/IBEXsrc/MP_Framework/ibex_rtl_engine.cpp` 中，对齐版 `BF_IBEX_RTL_CN` 是如何把 C model 的解码过程拉到 RTL 坐标系上的。

先给结论：

- `ldpc_config()` 仍然只负责从矩阵文本文件生成基准 `h_matrix`
- LUT 不是在 `ldpc_config()` 中加载，而是在 `ibex_rtl_engine.cpp` 中作为 sidecar 读取
- RTL 对齐版的关键，不是改 encoder，而是让 **syndrome 初始化、列更新、迭代回卷** 全部使用同一份 RTL 视图
- 这份 RTL 视图由 `build_rtl_view()` 构造，核心内容是 `shift / wrap_base / wrap_base_delta / mask`

换句话说，对齐版 C++ 做的不是“替换整个矩阵生成器”，而是：

```text
基准 h_matrix (来自 ldpc_config)
  + LUT sidecar (WrapBase / WrapBaseDelta / FirstMaskShift)
  -> build_rtl_view()
  -> 局部覆盖 h_matrix.element/mask
  -> BF_IBEX_RTL_CN 在 RTL 坐标系中完成整轮解码
```

### 12.1 片段 1：LUT 侧参数不是在 `ldpc_config()` 里读，而是在 RTL 引擎里读

`ldpc_config()` 只读文本矩阵文件：

```cpp
snprintf(matrix_path, sizeof(matrix_path),
         "%s/matrix/LDPC_%dx%dex%d_w4_dense5_QC_H_1.txt", ...);
snprintf(occupied_path, sizeof(occupied_path),
         "%s/occupied_matrix/LDPC_%dx%dex%d_w4_dense5_occupied_1.txt", ...);
snprintf(fade_path, sizeof(fade_path),
         "%s/fade_matrix/LDPC_%dx%dex%d_w4_dense5_fade_1.txt", ...);
```

而 RTL sidecar LUT 的接入点在 `get_ibex_family_lut()`：

```cpp
const char *default_paths[] = {
    "IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
    "../IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
    "../../../../IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
    "../../../../../IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
};
```

该函数读取的不是完整 H 矩阵，而是三类 RTL 参数：

- `wrap_base[row]`
- `wrap_base_deltas[row][(M,K)]`
- `first_mask_shift[M][K]`

因此当前 C++/RTL 对齐链路的结构是：

```text
文本矩阵文件 = 基准拓扑
LUT sidecar  = RTL 运行期参数
```

两者在 decoder 中合成，而不是由 `ldpc_config()` 一步生成。

### 12.2 片段 2：`build_rtl_view()` 用 LUT 重建 RTL 视图

`build_rtl_view()` 是对齐版的核心。它做三件事：

1. 决定 `wrap_base`
2. 从 `wrap_base + delta` 开始按 active CPM 从左到右重建整条 `shift` 链
3. 为当前 `(M,K)` 计算 `wrap_base_delta` 和 `mask`

对应代码片段：

```cpp
const bool use_lut_shifts = lut.loaded;

const int wb = (use_lut_shifts && (row < LDPC_MAX_ROWS))
                   ? lut.wrap_base[row]
                   : h_matrix.wrap_base[row];
view->wrap_base[row] = wb;

int next_shift = (wb + h_matrix.delta[row]) % h_matrix.bits;
for (int col = 0; col < h_matrix.cols; col++) {
  if (h_matrix.occupied[row][col] || h_matrix.fade[row][col]) {
    const int rebuilt_shift = next_shift;
    view->shift[row][col] = use_lut_shifts
                                ? rebuilt_shift
                                : h_matrix.element[row][col];
    next_shift = (next_shift + h_matrix.delta[row]) % h_matrix.bits;
  }
}
```

这段代码正对应本文前面 Section 4 的公式：

```text
S_r(0) = WrapBase + delta
S_r(n+1) = S_r(n) + delta
```

然后，`wrap_base_delta` 优先从 LUT 查；查不到时再搜索或回退：

```cpp
const int lut_wbd = use_lut_shifts
                        ? lookup_wrap_base_delta(lut, row, m, k_payload)
                        : -1;
if (lut_wbd >= 0) {
  view->wrap_base_delta[row] = lut_wbd;
} else {
  const int searched = search_wrap_base_delta(le, wb, h_matrix.delta[row],
                                              h_matrix.bits);
  view->wrap_base_delta[row] =
      (searched >= 0) ? searched : h_matrix.wrap_num_deltas[row];
}
```

最后一行如果存在 fractional parity，则按 RTL 语义重建 `mask`，并与 LUT 中的 `FirstMaskShift` 做一致性检查：

```cpp
if (offset < h_matrix.extra_bits_of_parity)
  view->mask[col][bit] = true;
```

这里的核心思想是：

- `h_matrix.occupied/fade` 仍定义“哪些列 active”
- `view->shift` 重新定义“这些 active 列在 RTL 坐标系中的 shift 是多少”
- `view->mask` 重新定义“last row 的哪个窗口是 occupied，哪个窗口是 fade”

### 12.3 片段 3：`ldpc_dec_bf_ibex_rtl_cn()` 在 syndrome 初始化前切换坐标系

第一版问题的根因，是 syndrome 初始化仍然在基准矩阵坐标系中做，后续 RTL 旋转却在另一套坐标系中进行，导致“初始 syndrome”与“列内/迭代内旋转”不在同一参考系里。

对齐版修复点就是在 `f_check_nodes()` 前，先把局部 `h_matrix` 切到 RTL 视图：

```cpp
const ibex_family_lut &family_lut = get_ibex_family_lut();
ibex_rtl_view rtl_view;
build_rtl_view(h_matrix, family_lut, &rtl_view);

if (family_lut.loaded) {
  for (int i = 0; i < h_matrix.rows; i++)
    for (int j = 0; j < h_matrix.cols; j++)
      if (rtl_view.shift[i][j] >= 0)
        h_matrix.element[i][j] = rtl_view.shift[i][j];
  if (h_matrix.extra_bits_of_parity > 0) {
    for (int j = 0; j < h_matrix.cols; j++)
      for (int k = 0; k < h_matrix.bits; k++)
        h_matrix.mask[j][k] = rtl_view.mask[j][k];
  }
}

cn = f_check_nodes(h_matrix, hard_codeword);
```

这里有两个重要点：

1. **覆盖的是值传递进来的局部 `h_matrix`**  
   不会污染 `ldpc_packet` 中保存的原始基准矩阵。

2. **`f_check_nodes()` 之后的初始 syndrome 已经在 RTL 坐标系里**  
   这样后续旋转、weight 统计、回卷就不会再和初始 syndrome 打架。

这一步就是交接文档里说的“全链路 RTL 坐标系对齐”。

### 12.4 片段 4：进入列循环前，先把每行 CN 旋到该行 first active shift

RTL 版本不是直接在“基准 element 坐标”上操作 `cn.r[i]`，而是每行先旋到 RTL 的 first active shift：

```cpp
for (int i = 0; i < h_matrix.rows; i++) {
  const int first_shift =
      first_active_shift_for_row(h_matrix, rtl_view, i);
  rotate_cn_row(&cn.r[i], h_matrix.bits, first_shift);
  accum_rot[i] = first_shift;
}
```

这里的含义是：

- `cn.r[i]` 从基准 syndrome 变成“RTL 当前参考坐标下的 syndrome”
- `accum_rot[i]` 记录当前行已经累计旋了多少
- 后续每列处理时，都要求 `accum_rot[i] == rtl_view.shift[i][j]`

因此 `accum_rot` 本质上是一个“当前 CN 参考系光标”。

### 12.5 片段 5：列更新时，weight 统计与 toggle 都按 `rtl_view.shift/mask` 执行

对齐版的列更新不再从 `h_matrix.element[i][j]` 反推索引，而是直接在“已经旋好的 CN 行”上，用当前列的 `rtl_view` 语义统计 weight：

```cpp
for (int i = 0; i < h_matrix.rows; i++) {
  const int shift = rtl_view.shift[i][j];
  if (shift < 0)
    continue;
  if (h_matrix.extra_bytes_of_parity == 0) {
    if (h_matrix.occupied[i][j] && cn.r[i].b[k])
      weight++;
  } else {
    if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)) &&
        cn.r[i].b[k])
      weight++;
    if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
        cn.r[i].b[k] && rtl_view.mask[j][k])
      weight++;
    if (h_matrix.fade[i][j] && cn.r[i].b[k] && !rtl_view.mask[j][k])
      weight++;
  }
}
```

然后 VN likelihood 更新仍然和基准 `BF_IBEX` 一样，是对 `vn.c[j].b[k].likelihood` 做更新：

```cpp
vn.c[j].b[k].likelihood = f_update_vn_post(
    vn.c[j].b[k].likelihood, adjusted_weight, likelihood_levels.min,
    likelihood_levels.max, prng_post_process, prng_post_process2,
    aggressive, likelihood_levels.flip_thr, pushing);
vn.c[j].b[k].flipped =
    (vn.c[j].b[k].likelihood >= likelihood_levels.flip_thr);
```

这意味着：

- **VN 决策语义与基准 BF_IBEX 保持一致**
- **变化的是 CN 权重统计所在的坐标系**

如果翻转状态变化，则 toggle 也按 `rtl_view.shift/mask` 对应的 active CPM 更新到 `cn.r[i].b[k]` 上。

### 12.6 片段 6：列尾只转 `delta`，迭代尾再做 `WrapBase + WrapBaseDelta × delta` 回卷

这是 RTL 行为最关键的结构特征。

每完成一列，当前活动行只前进一个 `delta`：

```cpp
for (int i = 0; i < h_matrix.rows; i++) {
  if (rtl_view.shift[i][j] >= 0) {
    rotate_cn_row(&cn.r[i], h_matrix.bits, h_matrix.delta[i]);
    accum_rot[i] = (accum_rot[i] + h_matrix.delta[i]) % h_matrix.bits;
  }
}
```

到一轮结束时，再做一次显式回卷：

```cpp
const int wrap_shift =
    (rtl_view.wrap_base[i] +
     (rtl_view.wrap_base_delta[i] * h_matrix.delta[i])) %
    h_matrix.bits;
rotate_cn_row(&cn.r[i], h_matrix.bits, wrap_shift);
accum_rot[i] = (accum_rot[i] + wrap_shift) % h_matrix.bits;
```

这正对应本文 Section 10 中讨论的 Version 2 语义：

- 列内只消费 `delta`
- 迭代边界显式做 `WrapBase + WrapBaseDelta × delta`

因此，对齐版 C++ 的 RTL_CN 实际上是在软件里显式复刻了 RTL 的“列内推进 + 迭代边界回卷”机制。

### 12.7 片段 7：诊断逻辑验证“坐标系真的闭合”

对齐版还保留了一组非常有用的诊断钩子：

- 初始旋转后，对比 `cn` 与 `cn_ref` 是否等价
- 每列开始前，检查 `accum_rot[i] == rtl_view.shift[i][j]`
- 每列前 8 bit，对比 RTL weight 与基准 ref weight
- 迭代回卷后，检查 `accum_rot[i]` 是否回到 `first_active_shift`

这组检查说明，当前实现不是“看起来像 RTL”，而是在软件中持续检查：

```text
初始坐标系
  -> 列内 delta 推进
  -> 迭代末 wrap 回卷
  -> 回到下一轮起始 first_active_shift
```

也就是要求旋转参考系在每轮结束时闭合。

### 12.8 这一版对齐究竟修复了什么

可以把修复前后的差异概括成一句话：

```text
修复前：Wrap 参数来自 LUT，但 syndrome 初始化和部分 shift 语义仍来自基准矩阵
修复后：shift / mask / syndrome_init / column_update / wrap 都在同一 RTL 坐标系
```

更具体地说，修复的是三处“坐标系混用”：

1. `build_rtl_view()` 中 `shift` 与 `wrap_base / wrap_base_delta` 来源不一致
2. `f_check_nodes()` 仍使用基准 `h_matrix.element`
3. last-row 的 `mask` 与 RTL first mask shift 不在同一条链上

对齐版通过 `build_rtl_view()` + 局部覆盖 `h_matrix.element/mask`，把这三处收成了一套统一语义。

### 12.9 一个容易混淆的点：encoder 没有切到 LUT 视图

当前实现中，encoder 入口仍是 `ldpc_ibex_encoder() -> f_ldpc_encode(..., h_matrix)`，它吃的是 `ldpc_config()` 产出的基准 `h_matrix`。

也就是说，本节讨论的“RTL 对齐版”是：

- **解码器内部视图对齐**
- 不是“整个工程把矩阵源切成 LUT 生成矩阵”

因此当前的设计意图是：

```text
编码端：沿用基准矩阵
RTL_CN：在 decoder 内部构造与 RTL 行为一致的局部视图
```

只要这份局部视图与基准矩阵处于同一 code space，这两条链就可以互相兼容。
