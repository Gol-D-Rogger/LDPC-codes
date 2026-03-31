# MP_Framework：IBEX DV C-Model（DPI）对齐说明

本文档面向 DV（SystemVerilog）环境的使用者，描述 `DVCtrans/IBEXsrc/MP_Framework/ldpc_c_model.c` 这一套 DPI-C 接口如何与 DV 对齐，重点解释 `ldpc_config(...)` 的参数语义、单位/缩放、与内部状态的对应关系，并对 `ldpc_c_model.c` 中每个函数的代码段逐段说明其作用与必要的对拍前提。

> 说明（必要的“怀疑态度”）：  
> 1) 本目录下 C-Model 采用 **MSB-first** 的 bit 打包/解包约定，并在 $blk\_len \bmod 32 = 16$ 时做 **2Byte-align** 修正（`<<16`）。若 DV 侧对 word 内 bit 的解释不是 MSB-first，则必须在 SV 侧做等价变换，否则会在“最后半 word”处出现系统性差异。  
> 2) `BF_IBEX` 模式依赖内部的 `rx_blk[]`（模拟接收值）与 `det_blk[]`（bin index）保持一致；在当前 DPI API 设计下，**推荐 DV 使用 `sd_err_inj(...)` 生成 soft-bin 输入**，不要绕开注错函数直接喂 `ldpc_dec(...)`，否则 `rx_blk` 可能不一致。

---

## 1. 目录与文件角色

- `ldpc_c_model.c`：DPI-C 入口层（DV 调用的 C API）。负责：
  - 解释 DV 传入的配置/数据格式（32-bit word 或逐 bit bin index）；
  - 计算长度域（Gen4 口径：DV“物理码字长度”与内部 padding 的关系）；
  - 调用 `ldpc_packet`（`ldpc_codec.cpp/.h`）与 `ch_packet`（IBEX `transceiver`）完成编解码与信道/软读建模；
  - 将 C/C++ 内部 bit/bin 向量重新打包回 DV 侧结构，并输出指标（收敛迭代、syndrome weight 等）。
- `ldpc_codec.cpp/.h`：IBEX 版 LDPC 内核（包含 `LAYER`、`BF_IBEX` 等译码路径，以及 `_LDPC_DBG_DUMP` 打印迁移）。

---

## 2. DV 对齐关键点总览（强制对齐项）

### 2.1 Bit 打包/解包：MSB-first

`ldpc_c_model.c` 的所有 “word ↔ bitstream” 都采用 MSB-first：

- 解包：第 $j$ 个 bit（$j\in[0,31]$）来自 `(word >> (31-j)) & 1`
- 打包：`tmp = (tmp << 1) | bit`（代码里写成 `tmp = tmp * 2 + bit`）

因此 DV 侧如果采用 LSB-first 视角做日志/比对，需要特别处理“尾 word”。

### 2.2 2Byte-align（仅当 $blk\_len \bmod 32 = 16$）

`ldpc_enc(...)` 与 `ldpc_dec(...)` 输出后，对最后一个 word 做：

- `last_word <<= 16`

其目标是与既有 Gen4 DVC 对拍策略一致（参见 `DVCtrans/doc/2Byte_Align_DVC.md` 的分析背景）。

### 2.3 长度域：Gen4 口径（DV 物理码字不包含内部 padding 0）

该 C-Model 在 DPI 层维护如下长度语义（用来与 DV 对齐）：

- DV 输入输出的码字长度：`blk_len = info_len + hm_m_phys`（**不包含** `pad_len` 个内部补零）
- LDPC 内核内部的系统信息段长度：`hm_k = (bm_n - bm_m) * h_sc`
- 内部 padding（shortening）：`pad_len = hm_k - info_len`

### 2.3.1 支持范围：仅 `h_sc==512`

MP_Framework 当前仅支持 `h_sc==512`，但矩阵来源已经切换为**外部矩阵文件加载**：

- `h_sc==512`：支持（从 `IBEX/ibex_matrix_flat_13rate` 或 `DVC_IBEX_MATRIX_ROOT` 指向的目录读取）
- `h_sc!=512`：不支持（直接报错返回）

外部矩阵加载规则如下：

- 优先读取与当前 `(bm_m, bm_n)` 完全同名的矩阵文件；
- 若目录中只有 full-width 版本（payload 固定 67 列），则按 **payload 取左侧、parity 取最右侧 `bm_m` 列** 的规则裁剪；
- 读取文件包括 `matrix`、`occupied_matrix`、`fade_matrix` 三类；
- 运行目录不固定时，建议显式设置环境变量 `DVC_IBEX_MATRIX_ROOT`。

当前实现对裁剪后的 `shift` 采用 **“先裁形状，再重生 shift”** 的口径，而不是直接沿用 full-width 文件中的 `shift`：

- 先从外部文件读取 full-width 的 `occupied_matrix` / `fade_matrix` / `matrix`；
- 裁剪时，真正决定目标矩阵形状的是裁剪后的 `occupied/fade` 分布；
- `matrix` 文件只用于给每一行提供 **最右侧 surviving 非零 CPM 的相位参考**；
- 之后在目标矩阵上按固定 `delta[row]` 从右往左重新生成 `h_matrix.element`；
- 因此，裁剪后同一行新的相邻非零 CPM（`occupied` 或 `fade`）仍保持单步 `delta[row]` 关系。

这点与“直接把 full-width 的 shift 抄到窄矩阵”不同。后者虽然简单，但在裁掉 payload 尾部列后，新的相邻非零 CPM 之间可能变成 `2*delta`、`3*delta` 等，不再符合 `IBEX/src` 的矩阵语义。

### 2.4 Soft decision 数据形态（`sd_num>=2`）

当 `sd_num>=2`：

- `sd_err_inj(...)` 输出的 `rx_data_sv` 是 **逐 bit 的 bin index**（`int` 数组，长度 `blk_len`），不是打包 word。
- `ldpc_dec(...)` 的 `det_data_sv` 在该模式下也期望 **逐 bit bin index**，并会对输入做夹紧到 `[0, bin_num-1]`。

### 2.5 `BF_IBEX` 与 `sd_num` 的强耦合

- `dec_mode==2` 且 `sd_num>=2` 时，`ldpc_dec(...)` 会选择 `BF_IBEX`。
- `BF_IBEX` 必须 `sd_num>=2`（否则直接报错返回）。
- `BF_IBEX` 在 decode 前会调用 `ldpc_ibex_input(...)`，其实现依赖 `rx_blk[]` 与 `det_blk[]`；因此推荐 DV 先走 `sd_err_inj(...)` 更新内部信道状态。

---

## 3. DPI API 一览（DV 调用序列建议）

### 3.1 推荐调用序列（硬判）

1. `ldpc_config(...)`
2. `ldpc_enc(usr_data, enc_data, ...)`
3. `ch_err_inj(enc_data, rx_data, ...)`（输出硬判 bitstream，打包 word）
4. `ldpc_dec(rx_data, dec_data, ..., dec_mode, sd_num=1, ...)`
5. `ldpc_cleanup()`

### 3.2 推荐调用序列（软读 + `BF_IBEX`）

1. `ldpc_config(..., sd_num>=2, ..., nand_strobes=sd_num, ...)`
2. `ldpc_enc(usr_data, enc_data, ...)`
3. `sd_err_inj(enc_data, bin_data, rd_data, llr_tbl, split_bin, bin_id, sd_num, v_ref, ...)`
4. `ldpc_dec(bin_data, dec_data, ..., dec_mode=2 或 BF_IBEX, sd_num>=2, ...)`
5. `ldpc_cleanup()`

> 反例（不推荐）：DV 绕过 `sd_err_inj`，直接构造 bin_data 喂 `ldpc_dec`。  
> 原因：`ldpc_ibex_input` 仍会读 `rx_blk[]`，此时 `rx_blk[]` 可能不是与 bin_data 同一组样本，导致不可解释的偏差。

### 3.3 DV 缓冲区维度（必须满足，否则容易“误判为算法错误”）

本实现大量使用 `svGetArrayPtr(...)` 直接取得连续内存指针，因此 DV 侧建议传入 **一维、连续** 的 `int` open-array（或等价 packed array，经 DPI 映射后为连续 `int`）。

记：

- $N_{info}=info\_num$（bit）
- $N_{blk}=blk\_len$（bit，见 §4.1 计算）
- $N_{cw}=h\_n\cdot h\_{sc}$（bit，`ldpc_dec` 会把输出补零到该长度）
- 32-bit word 数：$W(L)=\lceil L/32\rceil$

各 DPI 函数对数组维度的要求如下：

- `ldpc_enc(usr_data_sv, enc_data_sv, ...)`
  - `usr_data_sv`：`int` 数组，长度 $\ge W(N_{info})$
  - `enc_data_sv`：`int` 数组，长度 $\ge W(N_{blk})$
- `ch_err_inj(tx_data_sv, rx_data_sv, ...)`
  - `tx_data_sv`：`int` 数组，长度 $\ge W(N_{blk})$
  - `rx_data_sv`：`int` 数组，长度 $\ge W(N_{blk})$
- `sd_err_inj(tx_data_sv, rx_data_sv, rd_data_sv, llr_tbl_sv, split_bin_sv, bin_id_sv, sd_num, v_ref_sv, ...)`
  - `tx_data_sv`：`int` 数组，长度 $\ge W(N_{blk})$
  - `rx_data_sv`：`int` 数组（逐 bit bin index），长度 $\ge N_{blk}$
  - `rd_data_sv`：`int` 数组（打包后的 read bits），长度 $\ge W(N_{blk}\cdot sd\_num)$
  - `llr_tbl_sv`：`int` 数组，长度 $\ge 128$（输出固定写满 128）
  - `split_bin_sv`：`int` 数组，长度 $\ge 128$（输出固定写满 128）
  - `bin_id_sv`：`int` 数组，长度 $\ge 128$（输出固定写满 128）
  - `v_ref_sv`：`int` 数组，长度 $\ge sd\_num$
- `ldpc_dec(det_data_sv, dec_data_sv, dec_unc_sv, init_synd_wt_sv, dec_cnvg_itr_sv, dec_cnvg_col_sv, fina_synd_wt_sv, dec_mode, sd_num, ..., h_n, h_sc)`
  - `det_data_sv`：
    - 若 `sd_num<2`：`int` 数组（打包 bitstream），长度 $\ge W(N_{blk})$
    - 若 `sd_num>=2`：`int` 数组（逐 bit bin index），长度 $\ge N_{blk}$
  - `dec_data_sv`：`int` 数组（打包输出），长度 $\ge W(N_{cw})$
  - `dec_unc_sv/init_synd_wt_sv/dec_cnvg_itr_sv/dec_cnvg_col_sv/fina_synd_wt_sv`：指向 DV 侧 `int` 标量的指针（DPI `output int`/`inout int` 等形式）

---

## 4. `ldpc_config(...)` 参数含义（逐项解释 + 单位/缩放）

`ldpc_config(...)` 的参数数量较多，可按 6 组理解：矩阵/长度、信道、LLR/定点、译码器迭代、SDLite 覆盖、IBEX BF 寄存器组。

> 口径提醒：本实现采用 Gen4 文档中的“口径 A”（logical 维度）：`bm_m=h_m+1`、`bm_n=h_n+1`。

### 4.1 矩阵与长度域参数

| 参数 | 类型 | 含义（DV 侧） | C-Model 用途 |
|---|---:|---|---|
| `h_m` | `int` | base-matrix 行数的 logical 维度减 1（即 `bm_m-1`） | `bm_m=h_m+1`，用于长度计算与内部矩阵生成 |
| `h_n` | `int` | base-matrix 列数的 logical 维度减 1（即 `bm_n-1`） | `bm_n=h_n+1` |
| `h_sc` | `int` | circulant size（本实现仅支持 512） | 决定 `hm_*` 维度与 codec 路径（`cir_sz`） |
| `h_st` | `int` | H 矩阵生成/结构相关参数（IBEX codec 使用） | 透传至 `sim_pckt->ldpc_config(...)` |
| `h_wt` | `int` | H 矩阵权重/连接度相关参数 | 透传至 `sim_pckt->ldpc_config(...)` |
| `info_num` | `int` | 用户信息 bit 数 | `g_info_len=info_num` |
| `pad_bit` | `int` | “最后 circulant 中有效 bit 数”（Gen4 pad_bit 语义） | `mask_len=h_sc-pad_bit`，进入 parity/码长物理域计算；建议 `pad_bit` 为 8 的倍数以保持 byte 对齐 |

由此得到（DPI 层）物理长度域：

- $bm\_m=h\_m+1,\; bm\_n=h\_n+1$
- $mask\_len=h\_{sc}-pad\_bit$
- $hm\_m^{phys}=bm\_m\cdot h\_{sc}-mask\_len$
- $hm\_n^{phys}=bm\_n\cdot h\_{sc}-mask\_len$
- $hm\_k=(bm\_n-bm\_m)\cdot h\_{sc}$
- $pad\_len=hm\_k-info\_len$
- $blk\_len=info\_len+hm\_m^{phys}$

### 4.2 信道参数

| 参数 | 类型 | 含义 | 缩放/取值 |
|---|---:|---|---|
| `ch_mode` | `int` | 信道模式选择 | `0=CLEAN, 1=AWGN, 2=BSC, 3=ERR_INJ, 4=MAX_ERR` |
| `ch_para` | `int` | 信道参数 | 进入 C 后：`m_ch_para = ch_para / 100.0` |
| `sd_num` | `int` | 软读 read 次数（或硬判时=1） | `rd_num=sd_num`，并用于 LLR/bin 结构 |
| `v_ref_sv` | `svOpenArrayHandle` | vref 数组（长度 ≥ `sd_num`） | 每个元素按 `vref[i]=v_ref[i]/1000.0` 转成 `float` |

注意：本实现对 `sd_num` 做了夹紧 `1..127`，并对 `sd_num>=2` 才会读取 `v_ref_sv`。

### 4.3 LLR/定点与译码器共用参数

| 参数 | 类型 | 含义 | 缩放/说明 |
|---|---:|---|---|
| `llr0` | `int` | 硬判 0 的 LLR 基准（DV 定点） | `m_llr0 = llr0 / 16.0` |
| `llr1` | `int` | 硬判 1 的 LLR 基准（DV 定点） | `m_llr1 = llr1 / 16.0` |
| `alpha` | `int` | layer decoder 的 `alpha`（DV 定点） | `m_alpha = alpha / 100.0` |
| `finite_q_num` | `int` | LLR 定点总 bit（或量化等级数相关） | 传入 `ch_llr_gen(..., finite_q_num-1, 4)` 与 `ldpc_dec_config(..., finite_q_num, finite_r_num, 4, ...)` |
| `finite_r_num` | `int` | decoder 内部 R-message 定点位宽相关 | 透传至 `ldpc_dec_config` |

本实现固定 `finite_f_num=4`（分数位），因此：

- LLR 表内部 `float` 与 DV fixed-point 的换算，采用 $\times 2^{4}=16$（见 `sd_err_inj` 输出 `llr_tbl[i]=(int)(llr_tbl_float*16)`）。

### 4.4 译码器迭代参数（非 IBEX BF 寄存器组）

| 参数 | 类型 | 含义 | 用途 |
|---|---:|---|---|
| `fdec_max_itr` | `int` | BF 类 decoder 的最大迭代数 | 传入 `ldpc_dec_config(max_fdec_itr, ...)` |
| `ldec_max_itr` | `int` | `LAYER` decoder 的最大迭代数 | 传入 `ldpc_dec_config(..., max_ldec_itr, ...)`；若为负数则关闭 early-term（见 `ldpc_codec.cpp`） |

### 4.5 SDLite LLR override（透传给 Layer decoder）

| 参数 | 类型 | 含义 | 影响范围 |
|---|---:|---|---|
| `sdlite_llr_config` | `int` | SDLite override 使能 | `LAYER`/`LAYER2` 初始化时生效 |
| `sdlite_llr0..3` | `int` | 覆盖 `llr_tbl[0..3]` 的寄存器值 | 在 `bin_num>=4` 时生效；缩放为 `reg * 2^{-finite_f_num}` |

对应代码位置：

- `ldpc_dec_layer()`：`ldpc_codec.cpp:3214`
- `ldpc_dec_layer2()`：`ldpc_codec.cpp:3490`

### 4.6 IBEX BF（`BF_IBEX`）寄存器组（对齐 `ldpc_decoder_inv`）

这组参数用于复用 IBEX/DVsrc 的 `ldpc_decoder_inv` 风格寄存器解包逻辑：`ldpc_config` 阶段调用 `ldpc_ibex_parameters(...)` 下发到 `ldpc_decoder_parameters`；`ldpc_dec` 阶段在 `BF_IBEX` 模式调用 `ldpc_ibex_input(...)` 下发 `iteration_limit/post_iteration/nand_strobes` 并构造 `corrupted_codeword`。

| 参数 | 类型 | 含义 | 备注 |
|---|---:|---|---|
| `dv_user_data_bytes` | `int` | DV 侧认为的 userdata 字节数 | 仅用于与本地计算值比对并打印告警 |
| `dv_parity_bytes` | `int` | DV 侧认为的 parity 字节数 | 同上 |
| `post_iter` | `int` | post-processing 起始迭代 | 进入 `ldpc_ibex_input(..., post_iter, ...)` |
| `max_iter` | `int` | 最大迭代数 | 进入 `ldpc_ibex_input(..., max_iter, ...)` |
| `nand_strobes` | `int` | 软读路数（影响 soft_bits） | 建议与 `sd_num` 一致；不一致会打印 WARN |
| `codeword_4k_8k` | `int` | DV 侧区分 4k/8k codeword 的标志 | 当前仅打印保留（不影响算法） |
| `syndrome_weight_thr_qc` | `int` | syndrome 阈值（QC 阶段） | 下发至 `ldpc_decoder_parameters` |
| `syndrome_weight_thr_post` | `int` | syndrome 阈值（post 阶段） | 下发至 `ldpc_decoder_parameters` |
| `early_termination_dis` | `int` | 禁用 early termination | 下发至 `ldpc_decoder_parameters` |
| `post_process_en` | `int` | 启用 post-process | 下发至 `ldpc_decoder_parameters` |
| `ldpc_decoder_control_likelihood_0..3` | `int` | likelihood init 系数寄存器 | 在 `ldpc_ibex_parameters` 中按 nibble 解包 |
| `ldpc_decoder_control_post` | `int` | post 配置寄存器 | 解出 `likelihood_thr/post_ratio` |
| `ldpc_early_term_0..6` | `unsigned int` | early-term/likelihood map 寄存器 | `0..4` 用于 fraction/map；`5/6` 目前保留 |

---

## 5. `ldpc_c_model.c` 逐函数逐段说明（代码段级别）

本节按函数出现顺序解释每一段代码“做什么、为什么需要、对拍前提是什么”。

### 5.1 `ldpc_config(...)`

1) **重复配置保护**：若已配置（`g_configured!=0`），先调用 `ldpc_cleanup()`，避免内存泄漏与旧配置残留。  
2) **绑定全局句柄**：`sim_pckt=&ldpc_pckt`，后续所有 DPI API 均通过该对象访问 LDPC/CH 状态。  
3) **DV 定点 → 浮点换算**：将 `ch_para/alpha/llr0/llr1` 从 DV 定点转为 C 侧 `float`（缩放见 §4.2/4.3）。  
4) **长度域计算（Gen4 口径）**：根据 `h_m/h_n/h_sc/pad_bit/info_num` 计算 `g_info_len/g_blk_len/g_hm_k/g_hm_m/g_hm_n/g_pad_len`。  
5) **保存 BF_IBEX 运行期参数**：将 `post_iter/max_iter/nand_strobes` 保存到静态全局（供 `ldpc_dec` 的 `BF_IBEX` 路径使用）。  
6) **cir_sz 支持范围检查**：MP_Framework 仅支持 `h_sc==512`，若 `h_sc!=512` 则直接报错返回；`h_sc==512` 时从外部矩阵文件读取，并在需要时执行 payload-left/parity-right 裁剪。  
   - 裁剪后不会直接复用 full-width `shift`；  
   - 当前实现会依据裁剪后的 `occupied/fade` 形状、固定 `delta[row]` 和右侧 surviving 相位参考，重新生成目标矩阵的 `shift`，以保持一行内相邻非零 CPM 的 `delta` 连续性。  
7) **关键尺寸打印与一致性告警**：打印本地计算的 bit/byte 长度，并与 `dv_user_data_bytes/dv_parity_bytes` 做一致性检查（不一致只 WARN，不改变内部逻辑）。  
8) **vref 读取与 rd_num 夹紧**：当 `sd_num>=2` 从 `v_ref_sv` 读出 vref（单位转换），否则使用 `vref[0]=0`。  
9) **信道配置**：将 `ch_mode` 映射到 `enum ch_model` 并调用 `sim_pckt->ch_config(g_info_len, g_blk_len, ...)`。  
10) **LDPC 配置与缓存分配**：依次调用 `sim_pckt->ldpc_config(...)`（外部矩阵读取 + payload/parity 裁剪 + 裁后 shift 重建）、`ldpc_pckt_alloc()`。  
11) **LLR/soft-bin 结构初始化**：`ch_llr_alloc(MANUAL, rd_num, vref)` 使 `bin_num=rd_num+1`，随后 `ch_llr_gen(m_llr0, m_llr1, finite_q_num-1, 4)` 生成 `llr_tbl/bin_split/bin_id`。  
12) **译码器配置**：`ldpc_dec_config(..., sdlite_llr_*)` 透传 layer/BF 的迭代与 SDLite 覆盖寄存器。  
13) **IBEX BF 参数下发**：调用 `ldpc_ibex_parameters(...)`，将 `ldpc_decoder_inv` 风格寄存器解包到 `ldpc_decoder_parameters`。  
14) **完成标记**：`g_configured=1`。

### 5.2 `ldpc_cleanup()`

按“内核→缓存→LLR”的顺序释放并清空指针：

- `ldpc_clean()`：释放 H/矩阵相关
- `ldpc_pckt_clean()`：释放 LDPC 缓冲
- `ch_llr_clean()`：释放 vref/llr/bin 结构

最后 `g_configured=0`。

### 5.3 `ldpc_enc(usr_data_sv, enc_data_sv, debug)`

1) **配置检查**：未调用 `ldpc_config` 则报错返回。  
2) **解包 user bits（MSB-first）**：`usr_data` 为 32-bit word 流，解包写入 `sim_pckt->usr_blk[0..info_len-1]`。  
3) **调用编码器**：  
   - MP_Framework 仅支持 `cir_sz==512`，走 `ldpc_ibex_encoder()`（支持 parity 非整列/unused bytes 处理）。  
4) **打包 TX（MSB-first）**：从 `tx_blk[0..blk_len-1]` 打包回 `enc_data`。  
5) **2Byte-align 修正**：当 $blk\_len\bmod 32=16$，对 `enc_data` 最后一个 word `<<16`。

### 5.4 `ldpc_dec(det_data_sv, dec_data_sv, ..., dec_mode, sd_num, ..., h_n, h_sc)`

1) **配置检查**：未配置则报错返回。  
2) **输入格式分流**：  
   - `sd_num<2`：把 `det_data_sv` 当作打包 word（0/1 bitstream）解包到 `det_blk[]`。  
   - `sd_num>=2`：把 `det_data_sv` 当作逐 bit bin index（`int` 数组），并夹紧到 `[0, bin_num-1]` 写入 `det_blk[]`。  
3) **译码模式映射**：兼容 Gen4 旧接口与 IBEX 原生枚举值：  
   - `dec_mode==0 → BF_P3`  
   - `dec_mode==1 → LAYER`  
   - `dec_mode==2 → (sd_num>=2 ? BF_IBEX : BF_G2)`  
   - 或直接传 `BF_P3/BF_G2/BF_IBEX/LAYER/SKIP` 的枚举值。  
4) **soft 输入 + hard decoder 的兜底**：当 `sd_num>=2` 且选择了非 `LAYER/BF_IBEX/SKIP` 的 hard-decoder，会根据 `llr_tbl[bin]` 的符号把 bin 转为 hard bit（并打印 WARN）。  
5) **`BF_IBEX` 输入结构例化（关键）**：若选择 `BF_IBEX`：  
   - 强制 `sd_num>=2`；  
   - 若 `nand_strobes!=sd_num` 打印 WARN；  
   - 调用 `ldpc_ibex_input(0, g_max_iter, g_post_iter, g_nand_strobes)`：设置 `iteration_limit/post_iteration/nand_strobes/soft_bits`，并构造 `corrupted_codeword`。  
6) **运行译码器**：`sim_pckt->ldpc_decoder(model)`。  
7) **输出指标**：写回 `cw_fail/init_synd_wt/cnvg_itr/cnvg_lyr/fina_synd_wt`。其中 `BF_IBEX` 下的 `cnvg_lyr` 等价于 DVsrc 的 `col_cnt`（0-based 列号；未收敛或未进入列处理则为 `-1`）。  
8) **打包解码码字（MSB-first）**：从 `dec_blk[0..blk_len-1]` 打包到 `dec_data`。  
9) **2Byte-align 修正**：当 $blk\_len\bmod 32=16$，对最后一个 word `<<16`。  
10) **补零到固定宽度**：继续把 `dec_data` 填 0，直到覆盖 `h_n*h_sc` bit（以适配 DV 侧固定宽度 buffer 的对拍需求）。

### 5.5 `ch_update(ch_mode, ch_para, debug)`

仅更新信道参数（不重建矩阵、不重分配缓冲）：

- 将 `ch_mode/ch_para` 映射并调用 `sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->blk_len, ...)`。

### 5.6 `ch_err_inj(tx_data_sv, rx_data_sv, debug)`

硬判注错链路：

1) 解包 TX（MSB-first）到 `tx_blk[]`。  
2) `ch_transmit()` 生成 `rx_blk[]`。  
3) 以 0 为阈值生成 `det_blk[]`：`rx>=0 → 0`，`rx<0 → 1`。  
4) **输出打包策略（对齐尾 word）**：以输入的 `tx_data` word 为基线，仅覆盖前 `blk_len` 个 bit；超出 `blk_len` 的尾部 bit 保持 TX 原值，避免 `blk_len%32!=0` 时尾字被“清零”导致 DV 误判（该策略与 Gen4 适配一致）。

### 5.7 `sd_err_inj(tx_data_sv, rx_data_sv, rd_data_sv, llr_tbl_sv, split_bin_sv, bin_id_sv, sd_num, v_ref_sv, debug)`

软读注错链路（要求 `sd_num>=2`）：

1) `sd_num` 合法性检查与告警：`rd_num` 夹紧到 `<=127`；并对 `sd_num/vref` 与配置阶段的 `sim_pckt->rd_num/vref` 做一致性检查（不一致 WARN）。  
2) 解包 TX 到 `tx_blk[]`。  
3) `ch_transmit()` 生成 `rx_blk[]`。  
4) **Gen4/DVC 风格 soft detector（与 `DVCtrans/src` 对齐）**：  
   - 初始化 `det_blk[i]=0`；  
   - 对每次 read `rd_indx`：`rd_bit = (rx_blk[i] >= vref[rd_indx]) ? 0 : 1`；  
   - 若当前 `det_blk[i] == bin_split[rd_indx]` 且 `rd_bit==1`，则更新 `det_blk[i]=rd_indx+1`；  
   - 同时把每次 read 的 `rd_bit` 按“rd_indx 优先、bit_index 次之”的顺序打包输出到 `rd_data_sv`。  
5) 输出 `rx_data_sv`：逐 bit bin index（`int`），来自 `det_blk[]`。  
6) 输出 `llr_tbl_sv`：长度固定 128，按 `llr_tbl_float * 16` 转成定点（不足 `bin_num` 的填 0）。  
7) 输出 `split_bin_sv/bin_id_sv`：长度固定 128，多余位置填 0；其中 `bin_id` 输出 `rd_num+1` 个元素（含第 0 bin）。

---

## 6. `_LDPC_DBG_DUMP` 打印（已迁移）

当编译时定义 `_LDPC_DBG_DUMP`，`BF_IBEX` 路径会：

- 创建 `./output/` 目录
- 输出：
  - `output/ldpc_dbg_sw.txt`（syndrome weight 轨迹）
  - `output/ldpc_dbg_trace_iter*.txt`（每轮迭代 trace）

对应代码在 `ldpc_codec.cpp` 的 `ldpc_dec_bf_ibex(...)` 内（多处 `#ifdef _LDPC_DBG_DUMP` 块）。

---

## 7. DV 侧自检建议（避免“看起来像算法错，其实是接口错”）

1) **长度一致性**：确保 `info_num/ pad_bit/ h_m/ h_n/ h_sc` 与 DV 内部长度域一致；观察 `ldpc_config` 打印的 `user/parity/blk` 与 DV 传入的 `dv_user_data_bytes/dv_parity_bytes` 是否一致。  
2) **`sd_num` 一致性**：soft 模式时建议 `nand_strobes == sd_num`；不一致会改变 `soft_bits` 的意义。  
3) **vref 一致性**：`sd_err_inj` 会检测 call 侧 `vref` 与 config 侧 `vref` 是否一致（允许浮点微小误差）。  
4) **尾 word 对齐**：若出现“最后 16bit 翻转/挪位”，优先复核 MSB-first 与 2Byte-align 的假设（参见 `DVCtrans/doc/2Byte_Align_DVC.md`）。  
5) **`BF_IBEX` 输入来源**：若不调用 `sd_err_inj`，`rx_blk` 可能不是当前码字的模拟样本；建议先用推荐调用序列验证对拍，再考虑优化链路。

---

## 8. 问题记录：`LAYER_G2` 在 full-parity 场景的 parity 列错位（已修复）

> 本节只记录“数据链路/拼装逻辑”类问题；不记录与 `alpha` 缩放相关的性能问题。

### 8.1 现象（症状）

在 `dec_mode=LAYER_G2`（`enum dec_model=9`）且 **非 shortening/full-parity** 的配置下（典型特征：$Z=512$、`pad_len=0`、`unused_bytes_of_parity=0`，因此 $parity\_bytes = M\cdot Z/8$），出现：

- 全 0 `user_data` 时看起来“可以解码”，但全 1 `user_data` 时即便注错很少也会解码失败；
- `init_synd`/`fina_synd` 呈现“异常偏大且不随注错数量线性变化”的特征；
- 解码输出的 parity 段发生以 $Z=512$ 为单位的列级错位，并伴随“第一列 parity 变成全 0/最后一列被丢弃”的表象。

该类现象优先判定为 **码字视角的列拼装错误**，而非译码器算法本身性能不足。

### 8.2 证据（如何用最小打印快速判定）

在 `DVCtrans/IBEXsrc/MP_Framework/ldpc_c_model.c` 里，`debug>=2 && model==LAYER_G2` 时会打印：

- `parity_col_sig`：把 parity 段按列（每列 $Z=512$ bit）计算 FNV-1a 签名，分别对 `tx/rx/dec` 输出 10 个 32-bit 十六进制数；
  - 对 $Z=512$：全 0 列的签名恒为 `0x4d7705c5`；全 1 列恒为 `0x89c627c5`。

若观察到满足以下模式（示意）：

- `dec[0] == 0x4d7705c5`（第一列 parity 被置为全 0）
- 对 $c=1..9$：`dec[c] == rx[c-1]`（整体右移一列，最后一列丢失）

则可直接锁定为“parity 列起始偏移了 $Z$ bit”这一类拼装 bug。

### 8.3 根因（为何会出现“首列全 0 + 整体移位”）

根因位于 IBEX 原始实现的 `ldpc_packet::ldpc_decoder()` 中 `LAYER_G2` 的 **decoder 输入拼装**：

- `LAYER_G2` 的设计初衷是支持 shortening（第一列 parity 只有 `extra_bits_of_parity` 个有效 bit）；
- 但在 full-parity（`extra_bits_of_parity==0`）场景下，原逻辑仍按“fractional first parity column”来拼装，导致：
  - 第 0 列 parity 没有从 `det_blk[]` 拷贝进 `dec_di_blk[]`，而 `dec_di_blk` 由 `calloc` 初始化，因此保持全 0；
  - 后续 parity 从第 1 列开始拷贝，形成以 $Z$ 为单位的整体移位，并等价于“丢掉最后一列”。

> 注：该问题在全 0 `user_data` 上容易被掩盖（许多 parity 列本身可能接近全 0），而在全 1 `user_data` 上更容易暴露。

### 8.4 修复（只在 full-parity case 改拼装分支）

修复方式是：在 `LAYER_G2` 分支中显式区分 shortening 与 full-parity：

- 若 `extra_bits_of_parity > 0`：保留原本的 fractional-first-column 拼装；
- 若 `extra_bits_of_parity == 0`：按连续 parity 拷贝 `hm_m` bit（从 `det_blk[info_len]` 到 `dec_di_blk[hm_k]`），禁止跳过首列。

对应补丁位置：`IBEX/src/ldpc_codec.cpp` 的 `ldpc_packet::ldpc_decoder(enum dec_model)`（`dec_mode==LAYER_G2` 的输入拼装段）。

### 8.5 回归建议（最小集合）

建议用以下用例做“接口级回归”，以避免同类拼装问题再次引入：

1) **全 0 / 全 1 `user_data` 各一组**，固定注错 bit 数（例如 20）：
   - 观察 `parity_col_sig` 是否存在 `dec[0]=all0` 且 `dec[c]=rx[c-1]` 的模式；
   - 观察 `init_synd` 的量级是否与注错数量同阶（不应出现“注错很少但 `init_synd` 上千”的情况）。
2) **C-Model 内部一致性**：若 `tx_synd_wt_ibex`（`ch_err_inj` 打印）非 0，应先排查码字拼装/矩阵裁剪是否一致，再谈译码性能。
