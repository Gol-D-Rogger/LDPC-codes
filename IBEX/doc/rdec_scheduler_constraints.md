# RDEC 调度中的 `rdec_hdmem_cont_thrshd` / `rdec_cmem_cont_thrshd` 约束与含义

本文解释 `src/ldpc_codec.cpp` 的 `print_hm()` 在导出 `rdec_sched` 时做的两类“硬件约束检查”，并用一个小基矩阵例子说明：

- `rdec_hdmem_cont_thrshd`：同一 **列（col）** 的访问间隔约束（HD‑MEM/列相关存储）
- `rdec_cmem_cont_thrshd`：同一 **层/行（row）相关 C‑MEM** 的“上一行读回”间隔约束

> 重要定位：这两个阈值在当前代码里用于 **dump/验证调度**，不参与译码算法本身。它们在 `#ifdef _LDPC_DUMP` 下被赋值并触发 `print_hm()`（见 `src/ldpc_codec.cpp:1712`）。

---

## 1. 背景：为什么 layer‑decoder 需要调度

在 layer（按行/层）译码里，每处理一个非零 circulant（一个基矩阵非零位置 `(row=i, col=c)`）都会发生典型的两类存储访问：

1) **按列的存储（这里统称 HD‑MEM）**：用 `col=c` 索引，读/改/写该列的 Q/APP/硬判决相关状态。
2) **按行的存储（这里统称 C‑MEM）**：用 `row` 索引，读上一层（`pre_row`）的 C‑msg，再写回当前层 `i` 的 C‑msg。

如果硬件中这些存储是单端口/有流水写回延迟，就会出现“过近重复访问同一资源”的 hazard。因此需要：

- 对 **同一列** 的访问在时间上隔开足够多的调度条目；
- 对 **紧邻上一行**（`i-1`）的 C‑MEM 读回尽量延后，给上一行写回留出空档。

---

## 2. 调度 entry 是怎么生成的（对应代码路径）

`print_hm()` 中生成 `rdec_sched` 的核心结构是三层循环（见 `src/ldpc_codec.cpp:309` 起）：

1) 外层：按行 `i=0..m-1`（`bm_m`）
2) 中层：`j=1..m-1`
3) 内层：遍历行 `i` 上的每个非零 entry `e`

对每个 entry `e`，代码先找到它在同列中的“前驱 entry”：

- `e_pre = mod2sparse_prev_in_col(e)`，若无则回绕到 `last_in_col`（见 `src/ldpc_codec.cpp:325`）
- 因此 `pre_row = e_pre->row` 表示：**同一列里，上一层（上一轮写过该列的层）是哪一行**

然后只在满足下面条件时输出该条目：

- `pre_row == (i + j) % m`（见 `src/ldpc_codec.cpp:329`）

这意味着：对固定的行 `i`，它会把本行的 circulant 按 `pre_row` 的“相对距离”分组输出，顺序是：

`pre_row = i+1, i+2, ..., i+(m-2), i+(m-1)=i-1`（均按模 $m$）

也就是说：**依赖“紧邻上一行 $i-1$”的条目总是最后输出**，这是 C‑MEM 约束的关键设计点。

> 时间间隔的定义：调度文件里的每一条 entry 按顺序对应硬件的一个“调度时隙/周期”（至少是一个原子步骤）。因此“间隔”就是两条 entry 在输出序列里的索引差。

---

## 3. C‑MEM 约束：`rdec_cmem_cont_thrshd` 在检查什么

### 3.1 `row_dis` 的定义来自哪里

在每行 `i` 的中层循环里，当 `j == m-2` 时（见 `src/ldpc_codec.cpp:411`）：

- 这时刚好已经输出完 `pre_row = i+1 .. i+(m-2)` 这些组；
- **下一轮 `j=m-1` 就要开始输出 `pre_row=i-1`（紧邻上一行）这一组**。

代码把此刻累计输出的条目数 `tmp` 记作 `row_dis`：

`row_dis = tmp`

含义：本行 `i` 的调度里，在开始访问 `pre_row=i-1` 这组条目之前，已经插入了 `row_dis` 条“其它 pre_row 的条目”。

你也可以把它理解为：

- `row_dis = row_wt - count(pre_row==i-1)`

因为到 `j=m-2` 时，本行剩下没输出的只可能是 `pre_row=i-1` 这一组。

### 3.2 `rdec_cmem_cont_thrshd` 的约束含义

代码判定（见 `src/ldpc_codec.cpp:415`）：

- 若 `row_dis <= rdec_cmem_cont_thrshd`：报 “C‑MEM constraint violation”

因此 `rdec_cmem_cont_thrshd` 可以直观理解为：

> 硬件要求：在进入“读取紧邻上一行 $i-1$ 的 C‑MEM（作为 pre_row）”之前，至少要先执行超过 `rdec_cmem_cont_thrshd` 条其它条目，给上一行写回留足时间。

单位是“调度条目数/时隙数”（通常近似等于周期数，取决于你硬件每条 entry 实际耗时）。

---

## 4. HD‑MEM 约束：`rdec_hdmem_cont_thrshd` 在检查什么

`print_hm()` 会把每条输出 entry 的 `col` 依次记到 `sch_col[]`，然后做环形回看检查（见 `src/ldpc_codec.cpp:423`）：

对每个调度位置 `i`：

- 回看 `j=1..T`（其中 $T=rdec\_hdmem\_cont\_thrshd$）
- 若发现 `sch_col[i] == sch_col[i-j]`（用 `%cir_cnt` 回绕），则报错

等价语义：

> 同一列 `col` 在调度序列中两次出现之间的距离必须满足 $distance > T$，并且序列是环形的（需要满足尾→头回绕处的距离）。

这通常对应“按列存储”的流水/端口约束：同一列的读‑改‑写链路需要至少 $T$ 个时隙才能安全再次访问。

---

## 5. 小基矩阵例子：调度顺序、`row_dis`、HD‑MEM 间隔怎么计算

取一个小基矩阵（只关心非零位置，shift 值不影响“顺序”）：

基矩阵大小：$m=4, n=6$

|      | c0 | c1 | c2 | c3 | c4 | c5 |
|------|----|----|----|----|----|----|
| r0   | X  | X  | X  | .  | .  | .  |
| r1   | X  | .  | .  | X  | X  | .  |
| r2   | .  | X  | .  | X  | .  | X  |
| r3   | .  | .  | X  | .  | X  | X  |

同列“前驱行（pre_row）”由 `prev_in_col` 决定（行号升序 + 回绕）。例如：

- col0 在行 {0,1}：`pre_row(r0,c0)=1`，`pre_row(r1,c0)=0`
- col3 在行 {1,2}：`pre_row(r1,c3)=2`，`pre_row(r2,c3)=1`
- col5 在行 {2,3}：`pre_row(r2,c5)=3`，`pre_row(r3,c5)=2`

### 5.1 按代码规则生成每行的输出顺序

对每行 `i`，按 `pre_row=(i+j)%m`（$j=1,2,3$）分组输出：

- 行 r0（$i=0$）：
  - $j=1$ 目标 `pre_row=1`：输出 c0
  - $j=2$ 目标 `pre_row=2`：输出 c1
  - $j=3$ 目标 `pre_row=3=i-1`：输出 c2
  - 本行列序：`[0,1,2]`
  - `row_dis`（在 $j=m-2=2$ 时统计）：已输出 2 条 ⇒ `row_dis=2`

  这里“$j=1$ 输出 c0”的含义是：在本行 r0 里，遍历到列 c0/c1/c2 三个非零 entry 后，先计算它们的 `pre_row`：
  - c0：同列只在 r0、r1 有 X，因此对 (r0,c0) 来说 `prev_in_col` 会回绕到 r1 ⇒ `pre_row=1`
  - c1：同列只在 r0、r2 有 X，因此对 (r0,c1) 来说 `pre_row=2`
  - c2：同列只在 r0、r3 有 X，因此对 (r0,c2) 来说 `pre_row=3`

  然后在 `j=1` 这一轮，代码的筛选条件是 `pre_row == (i+j)%m == 1`，所以只会命中 c0 并输出该 entry（注意是输出 (r0,c0) 这一个 circulant 的调度条目，不是“输出整列”）。

- 行 r1（$i=1$）：
  - $j=1$ 目标 `pre_row=2`：输出 c3
  - $j=2$ 目标 `pre_row=3`：输出 c4
  - $j=3$ 目标 `pre_row=0=i-1`：输出 c0
  - 本行列序：`[3,4,0]`
  - `row_dis=2`

- 行 r2（$i=2$）：列序 `[\!5,1,3]`，`row_dis=2`
- 行 r3（$i=3$）：列序 `[\!2,4,5]`，`row_dis=2`

将 4 行串接起来得到整个调度列序列（`sch_col[]`）：

`[0,1,2, 3,4,0, 5,1,3, 2,4,5]`

### 5.2 这个序列如何用于两类约束

**(A) C‑MEM：看每行的 `row_dis`**

在这个例子里每行 `row_dis=2`：

- 若 `rdec_cmem_cont_thrshd=1`：因为 $2>1$，通过
- 若 `rdec_cmem_cont_thrshd=2`：因为 $2\le2$，会报 violation

**(B) HD‑MEM：看同列在序列中的最小间隔**

例如 col2 出现于位置 2 和 9，间隔为 7；col5 出现于位置 6 和 11，间隔为 5。

- 若 `rdec_hdmem_cont_thrshd=2`：所有列的最小间隔都 $>2$，通过
- 若 `rdec_hdmem_cont_thrshd=6`：col5 的最小间隔 5 不满足 $>6$，会报 violation  
  （并且代码是环形检查，还要同时满足尾→头回绕处的距离）

---

## 6. 使用建议（面向 RTL）

1) 这两个阈值本质上是**硬件时序/端口约束参数**，应从 RTL 的“资源复用周期/流水写回延迟/冲突消隐需求”反推，而不是从译码性能推。
2) 若你的 RTL 约束更严格（需要更大间隔），但 `print_hm()` 报 violation，说明：
   - 需要更强的调度算法（不仅是按 `pre_row` 分组），或
   - 需要硬件侧增加 bank/双端口/旁路，降低对间隔的要求。
3) 由于检查是“按条目数”计数，若硬件实现中“一条 entry 实际耗多个 cycle”，需要把阈值按比例换算到“条目尺度”或改检查口径。
