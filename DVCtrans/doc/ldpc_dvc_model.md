# DVCtrans LDPC DV C-Model 链路说明（DPI 接口视角）

## 1. 目的与适用范围

`DVCtrans` 的核心职责是作为 **DV（SystemVerilog）环境中的 DPI-C 参考模型**：将 DV 侧的 32-bit word 流（打包比特/软读数据）转换为 C/C++ 内部的 bit/LLR 向量，调用 LDPC 编解码与信道/软读模型，然后将结果重新打包回 DV 可直接对拍的数据结构，并输出收敛与 syndrome 等可观测指标。

本说明以 `DVCtrans/src/ldpc_c_model.c`（Gen3）与 `DVCtrans/src/ldpc_c_model_gen4.c`（Gen4 适配版）为入口，串联其对 `ldpc_packet`/`ch_packet` 的调用链路与长度语义。

## 2. 代码组织与关键模块

- DPI 入口（SV 调用的 C API）  
  - `DVCtrans/src/ldpc_c_model.c`：Gen3 口径  
  - `DVCtrans/src/ldpc_c_model_gen4.c`：Gen4 口径（引入 `pad_bit`、2Byte 对齐等）
- LDPC 内核（编解码 + 矩阵处理）  
  - `DVCtrans/src/ldpc_codec.cpp`, `DVCtrans/src/ldpc_codec.h`
- 信道/软读模型（注错、软读、LLR 表）  
  - `DVCtrans/src/transceiver.cpp`, `DVCtrans/src/transceiver.h`
- 工具与基础库  
  - `mod2sparse/mod2dense/mod2convert/vec_op/finite_lib/rand` 等

## 3. 核心对象：`sim_pckt` 与内存/长度域

DPI 入口文件中存在全局指针 `sim_pckt`，指向一个 `ldpc_packet` 实例（该实例继承自 `ch_packet`，因此同时持有信道与 LDPC 状态）：

- 信道域（`ch_packet`，见 `DVCtrans/src/transceiver.h`）  
  - `info_len`：信息比特数  
  - `blk_len`：DV 侧“物理码字长度”（通常是 $info\_len + parity\_len$，不含内部 padding）  
  - `tx_blk[]/rx_blk[]/det_blk[]/rd_blk[]`：发送比特、接收模拟值、检测结果、软读结果
- LDPC 域（`ldpc_packet`，见 `DVCtrans/src/ldpc_codec.h`）  
  - $hm_m,hm_n,hm_k$：展开后的 QC-H 维度（bit 级）  
  - `pad_len`：系统信息段内部 padding（shortening）长度  
  - `usr_blk/enc_di_blk/enc_do_blk/dec_di_blk/dec_do_blk/dec_blk`：编解码内部缓冲

### 3.1 Gen3（`ldpc_c_model.c`）的长度关系

配置输入为 base matrix 维度 `h_m,h_n` 与循环块大小 `h_sc`：

- $hm_m = h_m \cdot h_{sc}$  
- $hm_n = h_n \cdot h_{sc}$  
- $hm_k = hm_n - hm_m = (h_n-h_m)\cdot h_{sc}$  
- $pad\_len = hm_k - info\_len$

DV 对拍的“物理码字”（DVC 接口上传输/返回）采用 **shortened 码字**：

- $blk\_len = info\_len + hm_m$  
- 内部 padding 的 $pad\_len$ **在编解码内部存在，但 DV 不传输/不对拍**。

`ldpc_c_model.c` 中等价实现（见 `DVCtrans/src/ldpc_c_model.c:49-52`）：

- $h_k=h_n-h_m$  
- $pad\_num = h_k\cdot h_{sc}-info\_num$（即 $pad\_len$）  
- $blk\_num = h_n\cdot h_{sc}-pad\_num = info\_num + h_m\cdot h_{sc}$（即 $blk\_len$）

### 3.2 Gen4（`ldpc_c_model_gen4.c`）的长度关系（接口域）

Gen4 适配版 `ldpc_config(...)` 额外引入 `pad_bit`，并在 DPI 层维护“DV 物理长度域”的一套量（见 `DVCtrans/src/ldpc_c_model_gen4.c:62-71`）：

- **维度口径有两种（必须与 SV 侧一致）：**  
  - **口径 A（SV 传 logical）：**`h_m/h_n` 为逻辑 base-matrix 维度（不含最后的 padding&mask 行/列），则 DPI 内部取：  
    - $bm_m = h_m+1,\; bm_n=h_n+1$（物理 base matrix 扩展维度）  
  - **口径 B（SV 传 physical）：**`h_m/h_n` 已经是物理维度（含最后的 padding&mask 行/列），则 DPI 内部应取：  
    - $bm_m = h_m,\; bm_n=h_n$（不再 +1）  
  - 注：两种口径的区别直接影响 `pchk_file` 命名与 `blk_len` 计算；当前仓库内 `DVCtrans/src/ldpc_c_model_gen4.c` 采用的是 **口径 A**（见 `DVCtrans/src/ldpc_c_model_gen4.c:63-64`）。
- $mask\_len = h_{sc}-pad\_bit$  
- $parity\_len^{phys} = bm_m\cdot h_{sc} - mask\_len$  
  - 在口径 A 下：$parity\_len^{phys} = (h_m+1)\cdot h_{sc} - (h_{sc}-pad\_bit)= h_m\cdot h_{sc} + pad\_bit$  
  - 在口径 B 下：$parity\_len^{phys} = bm_m\cdot h_{sc} - (h_{sc}-pad\_bit)= (bm_m-1)\cdot h_{sc} + pad\_bit$
- $blk\_len^{phys} = info\_len + parity\_len^{phys}$

该长度域用于 **DV 输入输出与信道注错**；而 LDPC 内核实际采用的 `hm_*` 维度取决于链接的 codec 版本（在部门 Gen4 DV 环境下通常链接到支持 `mask_len/pad_bit` 的实现）。

## 4. 32-bit 打包约定（DV 已确认：MSB-first）

### 4.1 解包（word → bitstream）

以 MSB-first 方式从 `int word` 中取 bit：

- 第 $j$ 个 bit（$j\in[0,31]$）取自 $(word >> (31-j)) \& 1$  
- 代码模式见 `DVCtrans/src/ldpc_c_model.c:118-127`、`DVCtrans/src/ldpc_c_model.c:273-283` 等。

### 4.2 打包（bitstream → word）

多数接口用累乘方式构造 word：

- `tmp = tmp*2 + bit;` 重复 $j$ 次后写出 `tmp`。

注意：当 $blk\_len \bmod 32 \neq 0$ 时，最后一个 word 若直接写出 `tmp`，其有效比特会落在低位；而 DV 的 MSB-first 解包期望有效比特左对齐到高位。因此 Gen4 适配版额外做了 **2Byte 对齐修正**（示例：$blk\_len \bmod 32 = 16$ 时左移 16），见：

- `DVCtrans/src/ldpc_c_model_gen4.c:184-189`（ENC 输出）  
- `DVCtrans/src/ldpc_c_model_gen4.c:302-307`（DEC 输出）

### 4.3 DV 侧 buffer 维度建议（按“bit 长度”到“word/元素数”映射）

为避免 DPI 越界，建议 DV 侧显式按以下规则分配数组：

- 打包 word 数：$N_{word}(L)=\lceil L/32 \rceil$（承载 $L$ bit 的 `int` 个数）
- `usr_data_sv`：$N_{word}(info\_len)$  
- `enc_data_sv` / `tx_data_sv` / `rx_data_sv`（硬判码字）：$N_{word}(blk\_len)$  
- `dec_data_sv`：优先按 $N_{word}(blk\_len)$；若使用 Gen3 版 `ldpc_dec` 末尾补零逻辑，则按 $N_{word}(hm_n)$ 更安全  
- `sd_err_inj` 的 `rd_data_sv`：$N_{word}(blk\_len\cdot sd\_num)$  
- `sd_err_inj` 的 `rx_data_sv`（bin index 输出）：**逐 bit 一个元素**，应分配 `int rx_data[blk_len]`（非 word 打包）

### 4.4 DV 接口契约（建议作为 SV 侧约束）

本节将 DPI 层的 **位序**、**长度口径**、**数组维度**写成可直接用于 DV 对拍的“契约”。若 SV 端与本契约不一致，最常见表现是尾字对齐错误、`ch_err_inj` 后尾部被清零、或 `ldpc_dec` 输入/输出长度理解不一致。

#### 4.4.1 位序契约：MSB-first

- bitstream 序号 $i$（从 0 开始）映射到 32-bit word 的规则为：  
  - `word_index = i / 32`  
  - `bit_in_word = 31 - (i % 32)`  
  - 解包：`bit = (word[word_index] >> bit_in_word) & 1`（见 `DVCtrans/src/ldpc_c_model.c:124`、`DVCtrans/src/ldpc_c_model_gen4.c:147`）。
- 打包（`tmp = tmp*2 + bit`）对满 32bit 的 word 与上述解包规则严格互逆；但当 $L \bmod 32 \neq 0$ 时，**最后一个 word** 的有效比特会先落在低位，需要做“左对齐”修正（见 4.4.3）。

#### 4.4.2 长度口径契约：DV 对拍统一为 shortened `[Info|Parity]`

DV 输入输出（`ldpc_enc/ch_err_inj/ldpc_dec`）传输与对拍的码字统一理解为：

- `blk_len = info_len + parity_len_phys`
- 码字结构为 `[Info | Parity]`（不包含内部 Pad 段）

其中：

- Gen3：  
  - $hm_m = h_m\cdot h_{sc}$  
  - $hm_k = (h_n-h_m)\cdot h_{sc}$  
  - $pad\_len = hm_k - info\_len$（仅存在于编解码内部）  
  - $blk\_len = info\_len + hm_m$（DV 对拍长度，见 `DVCtrans/src/ldpc_codec.cpp:579-580`）
- Gen4（DVC wrapper 的“接口域”）：  
  - $mask\_len = h_{sc}-pad\_bit$  
  - $parity\_len^{phys} = bm_m\cdot h_{sc} - mask\_len$（其中 $bm_m/bm_n$ 的取值见 3.2 的两种口径）  
  - $blk\_len^{phys} = info\_len + parity\_len^{phys}$  
  - 注：Gen4 内核的 $hm\_*$ 维度与 `pad_len`/`mask_len` 的联动属于 codec 内核逻辑；但对 DV 侧而言，**对拍口径仍然是** `[Info|Parity]` 的 `blk_len_phys`。

#### 4.4.3 尾字对齐契约：$L \bmod 32$ 的处理

当 $L=blk\_len$（或 $L=info\_len$）不是 32 的整数倍时：

- MSB-first 解包期望：最后一个 word 的有效比特落在高位（bit 31 向下填充），低位为“无效填充位”。
- 因此 DPI 的打包实现应做左对齐：对最后一个 word 执行 `<< (32 - (L % 32))`。

当前 Gen4 wrapper 已实现最常见情况 $L \bmod 32 = 16$ 的修正：`last <<= 16`（见 `DVCtrans/src/ldpc_c_model_gen4.c:184-189`、`DVCtrans/src/ldpc_c_model_gen4.c:302-307`）。若未来出现其他余数（例如 8/24），建议推广为通用移位量。

**最小示例（$L \bmod 32 = 16$）：**  
设最后 16bit 的 bit 序列在数值上等于 `0xE70C`（按 bitstream 顺序从高位到低位）。若直接写出 `tmp`，最后一个 word 会变成 `0x0000E70C`，此时 MSB-first 解包读取 `bit[32]` 会从 bit31 开始读，得到的是 0（错位）。做 `<<16` 后得到 `0xE70C0000`，有效 16bit 被左对齐到 bit31..16，与 MSB-first 解包一致。

#### 4.4.4 DPI 函数的数组“形状”契约（硬判/软读）

以 `int` 为单位描述 DV 侧数组长度（元素个数）：

- `ldpc_enc(usr_data_sv, enc_data_sv, ...)`  
  - `usr_data_sv`：打包 bits，长度 $\lceil info\_len/32 \rceil$  
  - `enc_data_sv`：打包 bits，长度 $\lceil blk\_len/32 \rceil$  
  - 输出内容为 `[Info|Parity]`，不是 `[Info|Pad|Parity]`。
- `ch_err_inj(tx_data_sv, rx_data_sv, ...)`（硬判注错）  
  - `tx_data_sv`/`rx_data_sv`：打包 bits，长度 $\lceil blk\_len/32 \rceil$  
  - Gen4 wrapper 输出侧以 TX word 为基线，仅覆盖前 `blk_len` 个 bit（见 `DVCtrans/src/ldpc_c_model_gen4.c:398-418`），这是为避免尾字不足 32bit 时“末尾被清零”。
- `sd_err_inj(tx_data_sv, rx_data_sv, rd_data_sv, ..., sd_num, ...)`（软读/多阈值）  
  - `tx_data_sv`：打包 bits，长度 $\lceil blk\_len/32 \rceil$  
  - `rx_data_sv`：**逐 bit 一个 int**，长度 `blk_len`（内容为 bin index：$0..sd\_num$）  
  - `rd_data_sv`：打包 bits，长度 $\lceil blk\_len\cdot sd\_num/32 \rceil$
- `ldpc_dec(det_data_sv, dec_data_sv, ..., sd_num, ...)`  
  - 若 `sd_num < 2`（硬判）：`det_data_sv` 为打包 bits，长度 $\lceil blk\_len/32 \rceil$  
  - 若 `sd_num >= 2`（软读）：`det_data_sv` 为**逐 bit 的 bin index int**，长度 `blk_len`（见 `DVCtrans/src/ldpc_c_model.c:186-190`、`DVCtrans/src/ldpc_c_model_gen4.c:229-233`）  
  - `dec_data_sv`：打包 bits，长度 $\lceil blk\_len/32 \rceil$；输出码字结构为 `[Info|Parity]`。

## 5. Gen3 链路（`ldpc_c_model.c`）：DV 调用序列与数据语义

### 5.1 推荐调用序列（硬判 + 注错）

1. `ldpc_config(...)`：配置长度/信道/译码器并分配内存  
2. `ldpc_enc(usr_data, enc_data, ...)`：编码，输出 TX 码字  
3. `ch_err_inj(enc_data, rx_data, ...)`：信道传输/注错，输出硬判 RX 码字  
4. `ldpc_dec(rx_data, dec_data, ..., dec_mode, ...)`：译码，输出 `[Info|Parity]` 与指标  
5. `ldpc_cleanup()`：释放矩阵与缓存

### 5.2 `ldpc_config(...)`（Gen3）

职责：

- 计算 $pad\_len$ 与 $blk\_len$（shortening 口径）  
- `sim_pckt->ch_config(info_num, blk_num, ch_mode, ...)`：配置并确定 `info_len/blk_len`（DV 物理长度域）  
- `sim_pckt->ldpc_config(h_m, h_n, h_sc, h_st, h_wt, pchk_file)`：配置 QC-H 并生成编码所需矩阵  
- 分配缓冲：`ldpc_pckt_alloc()`、`ch_llr_alloc()`  
- LLR 表与译码器参数：`ch_llr_gen()`、`ldpc_dec_config()`

其中 `pchk_file` 命名规则见 `DVCtrans/src/ldpc_c_model.c:53`，要求仿真工作目录可访问该矩阵文件。

补充：Gen3 的 `ch_mode` 映射（见 `DVCtrans/src/ldpc_c_model.c:70-80` 与 `DVCtrans/src/transceiver.h`）：

- `0=CLEAN`, `1=AWGN`, `2=BSC`, `3=ERR_INJ`, `4=MAX_ERR`

### 5.3 `ldpc_enc(...)`（Gen3）

输入：`usr_data_sv`（32-bit word 数组，承载 $info\_len$ bit）  
输出：`enc_data_sv`（32-bit word 数组，承载 $blk\_len$ bit 的 TX 码字）

内部语义：

- 解包 `usr_data` → `usr_blk[0..info_len-1]`  
- 编码器内部会形成 `enc_di_blk[0..hm_k-1] = [Info | Pad(0)]`  
- 但对 DV 只输出 `tx_blk = [Info | Parity]`，长度 $blk\_len=info\_len+hm_m$（见 `DVCtrans/src/ldpc_codec.cpp:579-580`）

### 5.4 `ch_err_inj(...)`（Gen3，硬判注错）

输入：TX 码字 `tx_data`（$blk\_len$ bit）  
处理：

- 解包写入 `sim_pckt->tx_blk`（`DVCtrans/src/ldpc_c_model.c:273-283`）
- `ch_transmit()` 生成 `rx_blk[]`（模拟值，见 `DVCtrans/src/transceiver.cpp:104`）  
- `ch_detector(1,&vref)` 生成 `det_blk[]`（硬判，见 `DVCtrans/src/transceiver.cpp:253`）

输出：RX 码字 `rx_data`（打包自 `det_blk[0..blk_len-1]`）

### 5.5 `sd_err_inj(...)`（Gen3，软读/多阈值）

当 `sd_num>=2` 时，`ch_detector(sd_num, vref)` 的 `det_blk[i]` 不再是 0/1，而是 **bin index**（$0..sd\_num$），并同时生成 `rd_blk[rd_indx*blk_len + i]`（每次 read 的 0/1 结果）。

接口输出：

- `rx_data_sv`：逐 bit 的 bin index（长度 $blk\_len$）  
- `rd_data_sv`：打包后的 `rd_blk`（长度 $blk\_len\cdot sd\_num$ bit）  
- `llr_tbl_sv`：LLR 表（代码中做了 `*16` 的定点缩放，见 `DVCtrans/src/ldpc_c_model.c:401-403`）  
- `split_bin_sv/bin_id_sv`：bin 切分与 bin 编号信息

### 5.6 `ldpc_dec(...)`（Gen3）

输入：RX 码字 `det_data`（$blk\_len$ bit，硬判时为 0/1）  
内部会先在 codec 中重构 `dec_di_blk = [Info | Pad(0) | Parity]`（见 `DVCtrans/src/ldpc_codec.cpp:595-599`），再运行 BF/LAYER 等译码器，最终形成：

- `dec_blk = [Info | Parity]`（对 DV 输出的 shortened 码字口径）

输出：

- `dec_data_sv`：打包后的 `dec_blk`（$blk\_len$ bit）  
- 指标：`cw_fail/init_synd_wt/cnvg_itr/cnvg_lyr/fina_synd_wt`

注：`ldpc_c_model.c` 末尾还会把 `dec_data` 补零到 $hm_n=h_n\cdot h_{sc}$ 的 32-bit 边界（见 `DVCtrans/src/ldpc_c_model.c:224-231`），以适配某些 DV 端固定宽度 buffer。

补充：Gen3 的 `dec_mode` 映射（见 `DVCtrans/src/ldpc_c_model.c:192-199`）：

- `0 → BF_P3`（bit-flipping）  
- `1 → LAYER`（layered min-sum）  
- `2 → TBFDEC`（带 soft-bit 的 BF 变体）

### 5.7 `ch_update(...)`（Gen3，可选）

`DVCtrans/src/ldpc_c_model.c:235-255` 提供 `ch_update(ch_mode, ch_para, ...)`，用于在不重建矩阵的情况下动态切换信道参数（例如扫 BER/SNR），其内部重新调用 `ch_config(info_len, blk_len, ...)` 并保留已分配的 `tx_blk/det_blk` 等缓冲。

## 6. Gen4 链路（`ldpc_c_model_gen4.c`）：相对 Gen3 的关键差异

Gen4 入口与 Gen3 保持同名函数，主要差异集中在：

1) `ldpc_config` 增加 `pad_bit`，并在 DPI 层显式计算“物理 parity 长度”与 `blk_len`（见 `DVCtrans/src/ldpc_c_model_gen4.c:62-71`）。  
2) `ldpc_dec` 在进入 codec 之前，显式把输入从 `[Info|Parity]` 重构为解码器常用的 `[Info|Pad|Parity]`（见 `DVCtrans/src/ldpc_c_model_gen4.c:235-249`）。  
3) 2Byte 对齐（典型 $blk\_len\bmod32=16$）的 MSB-left 对齐修正：ENC/DEC 输出都对最后一个 word 做左移（`<<16`）。  
4) `ch_err_inj` 输出侧不再“从 0 重新打包”，而是以 TX 的 32-bit word 为基线，仅覆盖前 $blk\_len$ 个 bit，超出部分保持 TX 原值（见 `DVCtrans/src/ldpc_c_model_gen4.c:398-418`），以避免尾字不足 32bit 时出现“尾部被清零”的 DV 对齐问题。

补充：Gen4 的 `pchk_file` 命名使用扩展维度（见 `DVCtrans/src/ldpc_c_model_gen4.c:72-74`），因此矩阵文件名应与 $(bm_m,bm_n,h_{sc})$ 一致。

## 7. 常见踩坑清单（DV 对拍优先）

1) **长度域混用**：区分 $hm_n$（含 Pad 的数学码长）与 $blk\_len$（DV shortened 口径）。DV 对拍通常只看 `[Info|Parity]` 的 $blk\_len$。  
2) **MSB-first + 尾字不足 32bit**：若不做左对齐，DV 解包会把尾部读成 0。Gen4 已对 $blk\_len\bmod32=16$ 做修正；若出现其他余数，应推广为 `<< (32 - (blk_len % 32))` 的通用策略。  
3) **软读语义**：`sd_err_inj` 的 `rx_data_sv` 是 bin index（0..sd_num），不是 0/1 硬判。  
4) **矩阵文件路径**：`pchk_file` 在运行目录必须可见，否则 `ldpc_rd_phck` 会直接退出。
