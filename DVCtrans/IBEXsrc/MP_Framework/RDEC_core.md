# RDEC（LAYER_G2 / `ldpc_dec_layer2`）Dump 文件说明（DV 视角）

本文档面向 DV 使用与问题定位，目的在于解释 `ldpc_codec.cpp` 在打开 `_LDPC_DEBUG_DUMP` 时生成的 RDEC dump 文件**各自包含什么数据**、**输出格式如何解析**、以及这些数据与译码器内部变量（Q/R/APP/C-memory/HD 等）的对应关系。

> 约束：本文档中涉及数学符号时仅使用 `$...$` 形式，不使用双美元定界符。

## 0. 适用范围与前提

- 代码位置：`DVCtrans/IBEXsrc/MP_Framework/ldpc_codec.cpp`。
- 译码器：`ldpc_packet::ldpc_dec_layer2()`（`dec_mode == LAYER_G2`）。
- Dump 触发条件：仅在编译时定义 `_LDPC_DEBUG_DUMP` 时才会生成这些文件。
- Dump 目录：`./output/c_code/`（若不存在会在运行时创建）。
- 打开方式：`fopen(..., "w")`，即每次进入一次 `ldpc_dec_layer2()` 会**覆盖**旧文件内容。

## 1. 通用名词与维度约定（用于理解四类 dump）

RDEC（layered min-sum）在 QC-LDPC 下主要的维度与符号如下：

- $Z$：circulant size，对应代码里的 `cir_sz`。IBEX 常见 $Z=512$。
- `itr`：迭代编号（从 0 开始）。
- `layer`：base-matrix 的行编号（从 0 到 `bm_m-1`）。
- `col`：base-matrix 的列编号（从 0 到 `bm_n-1`）。
- 对某个 `col`，其 QC 码字视图在内存里是连续的 $Z$ 个 bit：`dec_do_blk[col*Z + i]`（$i\in[0,Z-1]$）。

### 1.1 软信息与量化打印约定

`rdec_log_dump.txt` 中对浮点软信息（`float`）的打印遵循以下规则：

- 先将数值按 `finite_f_num` 缩放：`vtmp = int(value * 2^{finite_f_num})`（代码中用 `pow(2, finite_f_num)`）。
- 再以符号 `+/-` 表示正负，幅度以十六进制显示：例如 `  7-2A` 表示索引 7 的值为负，量化幅度约为 `0x2A / 2^{finite_f_num}`。

### 1.2 512bit 向量的十六进制打包规则（HDMEM/STOT 通用）

`rdec_hdmem_dump.txt` 与 `rdec_stot_dump.txt` 都采用“每 4bit 输出 1 个 hex 字符”的打包方式：

- 外层按 nibble 逆序输出：从 `i = Z/4-1` 到 `0`。
- 每个 nibble 内的 bit 也按 MSB->LSB 累积：`j = 3..0`，累积式 `stmp = stmp*2 + bit`。

因此，对 $Z=512$：

- 每行 hex 字符数为 $512/4 = 128$。
- **最左侧**的第 1 个 hex 字符对应 bit 索引 `[511:508]`（其中 `vn_dec_hd[511]` 是该 nibble 的 MSB）。
- **最右侧**的最后 1 个 hex 字符对应 bit 索引 `[3:0]`。

> 若要与 RTL dump 对齐，必须先确认双方对 bit 索引与 nibble 输出顺序的约定是否一致；本文档描述的是 C-model 的实际输出规则。

## 2. dump 文件总览

在 `_LDPC_DEBUG_DUMP` 打开时，`ldpc_dec_layer2()` 会创建以下文件：

- `rdec_cmem_dump.txt`：每次调度到一个 `(itr, layer, col)` entry 后，输出该 circulant 喂给 C updater 的本地输入，也就是 `sign_tmp / val_tmp`。
- `rdec_stot_dump.txt`：每个 `layer` 的 `sign_tot` 位图（从 C-memory 派生，压成 128 个 hex 字符）。
- `rdec_hdmem_dump.txt`：每次调度到一个 `(itr, layer, col)` entry 时，对应 `col` 的 512bit 硬判决向量（HD）。
- `rdec_log_dump.txt`：每次调度到一个 `(itr, layer, col)` entry 时，输出该点位的 Q/R/APP 等软信息（以及可能的 bit 翻转日志）。
- `rdec_tblog_dump.txt`：面向 testbench 的“Q NEW”专用 dump（见 §6.4），用固定 token 形式输出符号位与幅度，便于脚本/仿真环境解析。

> 注：如果你的本地版本在 `rdec_log_dump.txt` 中额外插入了以 `[DVC]` 开头的对齐/探针日志，可忽略这些行；不影响本文对 “Q PRE / R NEW / APP-*” 的解析。

## 3. `rdec_cmem_dump.txt`（本地 C-update 输入 dump）

### 3.1 文件内容是什么

该文件记录每次处理完一个 `(itr, layer, col)` entry 时，当前 circulant 喂给 C updater 的本地输入（$i\in[0,Z-1]$）：

- `sign_tmp`：`sign(Q_new)`，映射成 0/1 打印（正号为 0，负号为 1）。
- `val_tmp`：`abs(Q_new)`；若 shortening/mask 逻辑判定该 bit 无效，则会被置成哨兵值 `100000`。

它不是累计后的 `cn_c_mem[layer]`，而是更新 `cn_c_updt_cur` 前的本地输入流。

### 3.2 输出触发时刻

在 per-circulant（per entry）循环内，每处理完一个当前 `col` 的 Q/R/C 更新后，立刻输出一次该 circulant 的 `sign_tmp / val_tmp`。

### 3.3 输出格式（逐行）

每行对应一个 circulant 内位置 `i`；连续的 $Z$ 行组成一个 block，对应“某个 circulant 的本地 C-update 输入”。不同 circulant block 之间插入一个空行分隔。

```
ITR<itr>/L<layer>/C<i>: sign_tmp <0|1>, val_tmp <hex>
```

字段说明：

- `ITR<itr>`：迭代号（0-based）。
- `L<layer>`：layer（base-matrix 行号）。
- `C<i>`：circulant 内索引 $i$（0..$Z-1$），这里的 `C` 不是 base-matrix column。
- `sign_tmp`：把内部 $\pm 1$ 映射为 0/1：正号（+1）打印为 0，负号（-1）打印为 1。
- `val_tmp`：以 `int(val_tmp * 2^{finite_f_num})` 打印的十六进制近似值；对被 mask/fade 移除的 bit，会看到较大的哨兵值。

> 由于现在是 per-circulant 快照，文件体积会显著增大，约为“原先 per-layer 版本乘以每层实际 circulant 数”。

## 4. `rdec_stot_dump.txt`（sign_tot bit-map）

### 4.1 文件内容是什么

对每个 `layer`，把 `cn_c_mem[layer][i].sign_tot`（$i\in[0,Z-1]$）压成 512bit 的位图后输出一行。

### 4.2 输出触发时刻

每个 `layer` 更新结束后输出一行；与 `rdec_cmem_dump.txt` 的 per-circulant 频率不同。

### 4.3 输出格式（逐行）

- 每行仅由十六进制字符组成（无 `ITR/LAYER` 头）。
- 对 $Z=512$，每行应为 128 个 hex 字符。
- bit 映射：`bit = (1 - sign_tot)/2`，即 `sign_tot=+1 -> 0`，`sign_tot=-1 -> 1`。
- nibble/bit 顺序：见 §1.2。

> 实务建议：`rdec_stot_dump.txt` 没有行头，解析时需用“第几行”对应到 `(itr,layer)`。其顺序严格跟随 `for (itr) for (layer)` 的循环顺序。

## 5. `rdec_hdmem_dump.txt`（hard-decision memory）

### 5.1 文件内容是什么

每次调度到一个 `(itr, layer, col)` entry（即遍历 `qc_bm` 的某个 `mod2entry e`），都会输出该 `col` 对应的 512bit 硬判决向量 `vn_dec_hd[0..Z-1]`：

- `vn_dec_hd` 的生成：由 `cn_app_pre` 的符号得到（`>=0 -> 0`，`<0 -> 1`）。
- 随后 `vn_dec_hd` 会被写入 `dec_do_blk[col*Z + i]`，因此它也可视为该调度点上该列的 HD 写回内容。

### 5.2 输出触发时刻

在 per-circulant（per entry）循环内，写回 `dec_do_blk` 后立刻输出一行。

### 5.3 输出格式（逐行）

```
ITRxx/LAYERyy/COLzz: <hex_512bits>
```

- `xx/yy/zz` 分别是 `itr / layer / col`（十进制，宽度为 2）。
- `<hex_512bits>` 为 128 个 hex 字符，对应 512bit，打包规则见 §1.2。

> 对齐提示：`LAYERyy` 是“当前遍历到的 base-matrix 行号”。若 RTL 的 schedule 把 `fade` entry 也作为一步输出，而 C-model 的 `qc_bm` 没有插入该 `fade` entry，则会出现“同一 `col` 的图样在 C-model dump 中落到下一条连接行（layer）才出现”的现象。该现象属于调度集合差异，而非 bit 打包错误。

## 6. `rdec_log_dump.txt`（Q/R/APP 等软信息日志）

该文件用于解释 `ldpc_dec_layer2()` 在某个 `(itr,layer,col)` 点位上的软信息计算过程。

### 6.1 主要输出块（每个 entry 一组）

每个 entry 先输出一行头：

```
ITRxx/LAYERyy/COLzz:
```

随后按组输出（每组会覆盖 0..$Z-1$ 的索引，按 8 个为一行分组）：

- `Q PRE MSG:`：`cn_q_sel_pre[idx]`，即该 `col` 的 Q-memory（变量节点到校验节点消息）。
- `R NEW MSG:`：`cn_r_new_pre[idx]`，由上一连接行（`e_pre->row`）的 C-memory 与 Q 的符号组合得到的 check-to-variable 消息。
- `APP-C MSG:`：`cn_app_pre[idx]`，近似为 `Q + Rnew`（在“上一连接行口径/CN 顺序”下的 APP）。
- `APP-S MSG:`：`cn_app_cur[idx]`，将 `APP-C` 做循环移位对齐到当前 `(layer,col)` 边的口径（S=shifted）。

### 6.2 每个数值的打印格式

每行包含多个 `idx±XX` token：

- `idx`：0..$Z-1$。
- `+/-`：该值的符号。
- `XX`：`abs(int(value * 2^{finite_f_num}))` 的十六进制表示。

### 6.3 “flip bit” 行（可选）

若启用了该段日志，则可能输出：

```
ITRxx/LAYERyy/COLzz: flip bit <i> (<old>--> <new>)
```

表示该 `col` 的某个 bit 在该调度点硬判决发生翻转（`dec_do_blk[col*Z+i]` 与新 `vn_dec_hd[i]` 不同）。

### 6.4 `rdec_tblog_dump.txt`（testbench 友好的 Q_NEW 日志）

一些版本会额外生成 `rdec_tblog_dump.txt`，其定位是“给 testbench/脚本喂数据或做逐 token 对齐”：

- 仅输出 `Q NEW`（对应 `cn_q_updt_cur`），不输出 `Q PRE/R NEW/APP-*`。
- 仍以 `ITRxx/LAYERyy/COLzz:` 作为 entry 的行头，便于与 `rdec_hdmem_dump.txt` 对齐同一个调度点。
- 每 8 个样本为一行，每个样本用固定 token 表示：
  - `SIGN:0/Q_NEW=..`：表示该位置 `Q_NEW` 为非负；
  - `SIGN:1/Q_NEW=..`：表示该位置 `Q_NEW` 为负；
  - `Q_NEW` 的数值为 `abs(int(cn_q_updt_cur[idx] * 2^{finite_f_num}))` 的十六进制打印（与 `rdec_log_dump.txt` 的幅度口径一致）。

与 `rdec_log_dump.txt` 的差异：

- `rdec_log_dump.txt`：面向人读，带 `idx±XX`，含 Q/R/APP/flip 等全链路信息；
- `rdec_tblog_dump.txt`：面向 TB 解析，去掉 `idx`、只保留 `SIGN/Q_NEW` 的稳定格式，便于 SystemVerilog/脚本逐 token 扫描。

## 7. `ldpc_dec_layer2()` 算法流程（逐段解释，DV 视角）

`ldpc_dec_layer2()` 实现的是 layered min-sum（下称 RDEC）译码。其主结构是：

- 外层迭代 `itr`：最多 `ldec_max_itr` 次；
- 每次迭代按 `layer=0..bm_m-1` 顺序处理；
- 每个 `layer` 内按 `qc_bm` 中该行的 entry（`col` 与 `shift`）依次更新：
  - 读该列的 Q（`cn_q_mem[col][*]`）；
  - 用该列“上一连接行”的 C-memory（`cn_c_mem[e_pre->row][*]`）生成 `R NEW`；
  - 合成 `APP` 并形成硬判决 `vn_dec_hd`，写回 `dec_do_blk[col*Z + *]`；
  - 更新该列的 Q（写回 `cn_q_mem[col][*]`）；
  - 在本层内积累新的 C-memory（`cn_c_updt_cur[*]`），层结束后写回 `cn_c_mem[layer][*]`。

下面按代码分段解释关键变量与数据链路，并给出与 dump 的对应关系。

### 7.1 输入/输出与关键成员（DV 需要知道的边界）

#### 7.1.1 主要输入

- `dec_di_blk[0..hm_n-1]`：译码输入，**LLR bin 的索引数组**（不是 float LLR）。该数组由 `ldpc_packet::ldpc_decoder(LAYER_G2)` 事先从 `det_blk` 拼装而来，并按 IBEX 的 shortening/尾部裁剪语义填入必要的“未发送 bit”的占位值（通常用 `max_llr_bin`）。
- `llr_tbl[0..bin_num-1]`：bin→LLR 的映射表（`float`）。`cn_q_mem` 初始化时会查表获得初始软信息。
- `qc_bm`：稀疏 base-matrix（`mod2sparse`），决定每个 `layer` 的调度 entry 集合（从 DV 角度可理解为 RDEC 的 schedule 基础）。
- `h_matrix`：包含 `occupied/fade/mask/extra_bits_of_parity` 等，用于 shortening 时按 bit gating。
- 配置参数：`alpha`（min-sum 缩放）、`finite_mode` 与 `finite_*`（量化参数）、`ldec_early_term_en`（早停开关）、`ldec_max_itr`（最大迭代）。

#### 7.1.2 主要输出

- `dec_do_blk[0..hm_n-1]`：译码的 full-QC 视图硬判决（0/1）。随后由 `ldpc_decoder()` 去 pad / 裁剪 parity tail，得到 DV 可见的 `dec_blk`。
- `cw_fail`：是否失败（0 表示满足收敛判据，1 表示未收敛）。
- `init_synd_wt` / `fina_synd_wt`：初始/最终 syndrome weight（以 0/1 硬比特计算）。
- `cnvg_itr` / `cnvg_lyr`：收敛点（早停启用时为实际收敛点；早停禁用时会被写成“跑满”的固定语义）。

### 7.2 初始化阶段（建立 Q-memory 与基线 syndrome）

#### 7.2.1 `dec_do_blk` 初值与 `dec_init` 的作用

函数开头把 `dec_di_blk` 直接复制到 `dec_do_blk` 作为占位：

- `vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);`

同时分配 `dec_init[bm_n]` 并全部置 1。它的语义是：**该列的 Q-memory 是否仍处于“初始 VN-order”**。第一次访问到某列时，需要用不同的 shift 规则把 APP/Q 变换到当前 layer 的坐标系；访问过后该列被标记为 0，后续用“相对 shift 差”做变换（见 §7.4.4）。

> 重要：由于 `dec_do_blk` 初始不是严格的硬判决，代码用 `hd_init=(vec_sum(dec_init)!=0)` 来禁止在“列尚未全部访问一遍”时触发早停计数。

#### 7.2.2 SDLite LLR override（如启用）

若 `reg_sdlite_llr_config` 置位且 `bin_num>=4`，代码会覆盖 `llr_tbl[0..3]`：

- `llr_tbl[k] = reg_sdlite_llrk * pow(2, -finite_f_num)`（即按 $2^{-finite\_f\_num}$ 缩放）

这会直接影响后续 Q 初始化与所有迭代。

#### 7.2.3 `init_synd_wt` 的计算

`init_synd_wt` 不直接用 `dec_di_blk`（bin）计算，而是先得到 hard 判决：

- `hard_init[t] = (llr_tbl[dec_di_blk[t]] >= 0) ? 0 : 1`

再计算：

- `init_synd_wt = dvc_ibex_syndrome_weight(this, hard_init);`

它表示“信道初判决码字”的 syndrome weight。

#### 7.2.4 Q-memory 初始化（bin→LLR）

`cn_q_mem[bm_n][Z]` 初始化为：

- `cn_q_mem[col][i] = llr_tbl[ dec_di_blk[col*Z + i] ]`

若 `finite_mode==1`，对其做饱和量化（`Sat_Quan`）。

### 7.3 主循环结构（`itr → layer → entry`）

可用以下伪代码概括：

```
for itr = 0..ldec_max_itr while (!early_term || cw_fail):
  cir_cnt = 0
  for layer = 0..bm_m-1 while (!early_term || cw_fail):
    hd_init = any(dec_init[col]==1)
    init cn_c_updt_cur[*] (min=INF, sign_tot=+1)
    hd_updated = 0
    layer_synd[*] = 0
    for e in qc_bm.row(layer):  // e=(row=layer, col, shift)
      ... per-entry update ...
      cir_cnt++
    write back cn_c_mem[layer][*] from cn_c_updt_cur[*]
    early-termination bookkeeping
```

### 7.4 每个 entry（`e=(layer,col,shift)`）的核心计算

以下解释对应 per-entry 循环体中几个关键步骤；同时也说明 `rdec_log_dump` 与 `rdec_hdmem_dump` 的对应变量。

#### 7.4.1 `e_pre`：同一列的上一连接行（循环意义）

对同一 `col`，代码会寻找该列在 `qc_bm` 中“上一条 entry”：

- `e_pre = mod2sparse_prev_in_col(e)`，若到头则 wrap 到 `mod2sparse_last_in_col(qc_bm, col)`。

语义：layered decoder 需要使用“该列在上一连接行刚更新过的 check 信息”，以实现逐层更新的快速收敛。

随后读取“上一连接行”的 C-memory：

- `cn_c_sel_pre = cn_c_mem[e_pre->row]`

#### 7.4.2 `R NEW`：由 C-memory 生成的 min-sum 校验消息

对每个 $i\in[0,Z-1]$：

- `sign_Q = sign(cn_q_sel_pre[i])`（`>=0 -> +1`，`<0 -> -1`）
- `Rnew = (min1 or min2) * sign_tot * sign_Q`
  - 若 `min1_pos == col` 则取 `min2`，否则取 `min1`

这等价于标准 min-sum 形式：

- $R_{c\to v} = \left(\prod_{v'\in N(c)\setminus v}\text{sign}(Q_{v'\to c})\right)\cdot \min_{v'\in N(c)\setminus v}|Q_{v'\to c}|$

其中 `sign_tot * sign_Q` 用来实现“排除当前列后的符号”（因为 $\text{sign}(Q)^2=1$）。

#### 7.4.3 `APP-C`：$APP = Q + Rnew$（shortening 时按 bit gating）

full-parity（`h_matrix.extra_bytes_of_parity == 0`）：

- `cn_app_pre[i] = cn_q_sel_pre[i] + cn_r_new_pre[i]`

shortening（`extra_bytes_of_parity != 0`）：

- 代码会按 `occupied/fade/mask` 对 `Rnew` 是否参与 APP 做 gating。其效果可以理解为：
  - **有效 bit**：`APP = Q + Rnew`
  - **无效 bit**：`APP = Q`（禁用该 check 反馈）

若 `finite_mode==1`，对 `cn_app_pre` 做饱和量化。

#### 7.4.4 `APP-S` 与硬判决：shift 对齐与写回 `dec_do_blk`

代码用两个 shift 值把 `cn_app_pre` 从“上一连接行坐标系”变换到：

- 当前连接（当前 layer）的坐标系（得到 `cn_app_cur`，即 `APP-S`）
- 列的 VN 顺序坐标系（得到 `vn_dec_hd`，即硬判决）

计算规则：

- 若该列首次访问（`dec_init[col]==1`）：
  - `shift_val1 = e->shift`
  - `shift_val2 = 0`
  - `dec_init[col]=0`
- 否则：
  - `shift_val1 = -e_pre->shift + e->shift`
  - `shift_val2 = -e_pre->shift`

随后：

- `cn_app_cur[i] = cn_app_pre[(i + shift_val1) mod Z]`
- `vn_dec_hd[i] = (cn_app_pre[(i + shift_val2) mod Z] >= 0) ? 0 : 1`

并写回：

- `dec_do_blk[col*Z + i] = vn_dec_hd[i]`

因此：

- `rdec_hdmem_dump.txt` 一行对应一次 `(itr,layer,col)` 写回的 `vn_dec_hd`（按 §1.2 打包）。
- `rdec_log_dump.txt` 的 `Q PRE/R NEW/APP-C/APP-S` 就对应 `cn_q_sel_pre/cn_r_new_pre/cn_app_pre/cn_app_cur`。

#### 7.4.5 layer syndrome：`layer_synd` 的累积与含义

对当前 entry，先把 `vn_dec_hd` 按 `e->shift` 转到 CN 口径：

- `vec_shift(vn_dec_hd, cn_dec_hd, Z, -e->shift)`

结合 `vec_shift` 的实现，可等价写作：

- `cn_dec_hd[t] = vn_dec_hd[(t + e->shift) mod Z]`

然后将其 XOR 累积到本层 syndrome：

- `layer_synd ^= cn_dec_hd`

shortening 时，代码会在 XOR 前按 `mask` 将某些位置清零，避免无效 bit 影响 syndrome 统计与早停判据。

#### 7.4.6 `R OLD`、`Q NEW` 与 `cn_q_sign[cir_cnt][i]`

layered 更新需要从 `APP-S` 中“去掉旧的 R，再加入新的 R”。因此要计算本层对该列的旧消息：

- `cn_c_sel_cur = cn_c_mem[layer]`（该层旧 C-memory）
- `cn_q_sign[cir_cnt][i]`：该 entry 上一次更新时的 `sign(Q_new)`（每个 i 一个符号）

构造：

- `Rold = (min1 or min2) * sign_tot * sign(Q)`（与 Rnew 同形，但来源是“当前层旧 C-memory”）

再更新：

- $Q_{\text{new}} = APP\_S - R_{\text{old}}$

在代码中即：

- 先生成 `cn_r_old_cur[i]`（对应 $R_{\text{old}}$），再计算 `cn_q_updt_cur[i]`（对应 $Q_{\text{new}}$）：

```
if (cn_c_sel_cur[i].min1_pos == e->col)
  cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
else
  cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

// Q -= Rold
cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
```

shortening（`h_matrix.extra_bytes_of_parity != 0`）时，代码会在某些“mask 无效位”场景下将 `cn_r_old_cur[i]` 强制置 0，再做减法（等价于该 bit 不存在旧消息、避免把无效位的旧反馈减回去）：

```
if (h_matrix.occupied[e_pre->row][e_pre->col] && (e_pre->row == (h_matrix.rows - 1)) &&
    (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
  cn_r_old_cur[i] = 0;
}
if (h_matrix.fade[e_pre->row][e_pre->col] &&
    (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
  cn_r_old_cur[i] = 0;
}
cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
```

并可选量化。

> 为什么必须存 `cn_q_sign[cir_cnt][i]`：C-memory 里只存了整层的 `sign_tot` 与 `min1/min2`，要恢复“排除当前列后的符号”必须知道当前列的符号；该符号随迭代变化，不能从 C-memory 反推，因此必须显式保存。

> 初始化语义：`cn_q_sign` 由 `calloc` 分配，初值为 0。第一次迭代/第一次访问某个 entry 时，`cn_r_old_cur[i]` 会自然变为 0（相当于“旧消息不存在”），随后在本 entry 末尾 `cn_q_sign[cir_cnt][i]` 会被更新为 `+1/-1`，供后续迭代使用。

#### 7.4.7 写回 Q-memory，并累积更新本层 C-memory

对该 entry：

1) 写回该列 Q-memory：
   - `cn_q_mem[col][i] = cn_q_updt_cur[i]`

2) 用 `cn_q_updt_cur[i]` 更新本层累积 C-memory（`cn_c_updt_cur[i]`）：
   - `sign_tmp = sign(Q_new)`
   - `val_tmp = |Q_new|`
   - 更新 `sign_tot` 与 `min1/min2/min1_pos`
   - 保存 `cn_q_sign[cir_cnt][i] = sign_tmp`

shortening 时，对无效 bit 位置会采用：

- `sign_tmp = +1`
- `val_tmp = 100000`

其效果是“该 bit 在该边上等价被移除”：既不参与 min1/min2 竞争，也不影响 sign_tot。

当该层所有 entry 处理完后，写回该层 C-memory：

- `cn_c_mem[layer].min* = cn_c_updt_cur.min* * alpha`
- 必要时用 `finite_c_*` 进行量化
- `min1_pos/sign_tot` 直接拷贝

因此：

- `rdec_cmem_dump.txt` 记录的是每个 `(itr,layer,col)` 处理完成后的本地 C-update 输入，也就是 `sign_tmp / val_tmp`；
- `rdec_stot_dump.txt` 是其 `sign_tot` 的 512bit 压缩表示。

#### 7.4.8 早停判据（`hd_init/hd_updated/layer_synd` 的组合）

每层结束后计算：

- `layer_synd_wt = sum(layer_synd)`（0/1 求和）

并维护两个计数器：

- `synd_pass_cnt`：连续满足“`layer_synd_wt==0` 且 `hd_updated==0`”的层数；
- `hd_stable_cnt`：同上，但阈值不同。

同时用 `hd_init = (vec_sum(dec_init)!=0)` 保护：只要还有列未在本次迭代被访问过，就认为硬判决尚在初始化阶段，不累计早停条件。

当满足：

- `synd_pass_cnt >= bm_m` 且 `hd_stable_cnt >= bm_m-1`

则宣布收敛：

- `cw_fail = 0; cnvg_itr = itr; cnvg_lyr = layer;`

若 `ldec_early_term_en==1`，外层循环会在 `cw_fail==0` 后停止；若 `ldec_early_term_en==0`，函数末尾会把 `cnvg_itr/cnvg_lyr` 固定写为“跑满”语义（但 `cw_fail` 仍可能为 0）。

### 7.5 shortening（`extra_bits_of_parity > 0`）下 `occupied/fade/mask` 的逐段说明

shortening 的核心问题是：有些 parity 列（尤其 first parity column）只有部分 bit 是有效码字位，其余 bit 在 DV 视角是“未发送/被裁剪”的。为避免这些 bit 干扰译码，代码通过 `mask[col][k]` 做按 bit gating。

`mask` 的典型使用形态是 `mask[col][(i+shift)%Z]`：因为 CN 侧索引 `i` 对应的 VN 侧 bit 位置是 `(i+shift) mod Z`。

在 `ldpc_dec_layer2()` 内部，mask 主要影响 3 处：

1) **APP 计算阶段：决定 `APP-C` 是否包含 `R NEW`**  
   代码在 `extra_bytes_of_parity != 0` 分支中使用如下判定（逐条对应 if/else）：
   - 若 `occupied[e_pre->row][e_pre->col]` 且 `e_pre->row < rows-1`：`APP = Q + Rnew`
   - 否则若 `occupied[e_pre->row][e_pre->col]` 且 `e_pre->row == rows-1` 且 `mask[e_pre->col][(i + e_pre->shift) % Z]`：`APP = Q + Rnew`
   - 否则若 `fade[e_pre->row][e_pre->col]` 且 `!mask[e_pre->col][(i + e_pre->shift) % Z]`：`APP = Q + Rnew`
   - 否则：`APP = Q`

2) **syndrome 累积阶段：对无效 bit 清零后再 XOR**  
   在 `extra_bytes_of_parity != 0` 时，代码会在某些条件下对 `cn_dec_hd` 执行：
   - 若 `!mask[e->col][(i + e->shift) % Z]` 则 `cn_dec_hd[i]=0`
   以避免无效 bit 影响 `layer_synd`。

3) **`R OLD / Q NEW` 与 C-memory 更新阶段：对无效 bit “边移除”**  
   在 `extra_bytes_of_parity != 0` 时，代码会对无效 bit：
   - 将 `cn_r_old_cur[i]=0`（避免减去不存在/不应存在的旧消息）
   - 将 `val_tmp=100000` 且 `sign_tmp=+1`（避免该 bit 参与 min1/min2 与 sign_tot）

> 重要提示：shortening 分支里部分 gating 条件同时引用 `e` 与 `e_pre`（上一连接行）。从 DV 视角可把它理解成“在不同连接语义/不同坐标系下做 enable 判定”，其最终效果是：无效 bit 不参与译码约束。若要与 RTL 逐拍对齐，应以 RTL 对 mask/fade 的定义为准逐条对照。

### 7.6 `qc_bm` / RDEC 调度与 `fade=1` 的关系（影响 dump 对齐）

`ldpc_dec_layer2()` 的 per-entry 遍历完全由 `qc_bm` 决定：

- `for (e = mod2sparse_first_in_row(qc_bm, layer); ...)`

而 `qc_bm` 在 `ldpc_ibex_phck()` 中的插入条件依赖 `extra_bits_of_parity`：

- full-parity（`extra_bits_of_parity==0`）：只插入 `occupied==1` 的位置（`fade==1` 不进入 `qc_bm`）。
- shortening（`extra_bits_of_parity>0`）：插入 `element>=0` 的位置（通常会包含 `fade` 边，因为 fade 边的 `element` 也被赋值）。

因此，当 RTL trace 把 `fade` entry 也作为 schedule 步骤输出（例如 `layer=0,col=42` 恰为 `fade=1`）而 C-model 当前处于 full-parity 分支时：

- `rdec_hdmem_dump.txt` 在该 `layer` 不会出现该 `col` 的 entry；
- 你会看到“同一列的 data pattern 在 C-model dump 中落到其它 layer 才出现”的现象。其根因是**调度集合差异**，而非硬判决打包或数值计算错误。
