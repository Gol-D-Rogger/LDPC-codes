# Gen4 原生 `ldpc_codec.cpp` 编码器矩阵结构与 `_LDPC_DEBUG_DUMP` 说明

本文针对 `gen4_ldpc_sim/src/ldpc_codec.cpp` 中的 `ldpc_packet::ldpc_encoder()`（原生 Gen4 encoder，而非 KY 版），详细说明其矩阵结构、编码流程以及 `_LDPC_DEBUG_DUMP` 下各 dump 文件的数据长度和含义。与 KY 版说明（`gen4_encoder.md`）分开记录，便于对比两套实现。

---

## 1. 矩阵维度与结构概览

### 1.1 维度定义

和 KY 版类似，原生 `ldpc_codec.cpp` 也使用 QC 结构，关键维度如下：

- $Z = cir\_{sz}$：循环块大小（代码中 `cir_sz`）；
- 基矩阵维度：
  - $bm_m$：行数；
  - $bm_n$：列数；
  - $bm_k = bm_n - bm_m$；
- 展开的 H 矩阵维度（bit 域）：
  - $hm_m = bm_m \cdot Z$；
  - $hm_n = bm_n \cdot Z$；
  - $hm_k = hm_n - hm_m$；
- 其他相关：
  - `info_len`：真实用户信息长度；
  - `pad_len = hm_k - info_len`：信息端内部 padding 长度；
  - `blk_len`：物理码长（`ch_packet::blk_len`）；
  - `pad_bit`：最后一个循环块中有效位数，`cir_sz - pad_bit` 对应被“截掉”的部分。

在 `ldpc_pckt_alloc()` 中，encoder/decoder 内部使用的缓冲长度为：

- `usr_blk`：信息输入，长度 `info_len`；
- `enc_di_blk`：信息段 $Z_1$，长度 `hm_k`；
- `enc_do_blk`：内部完整码字，长度 `hm_n`；
- `dec_di_blk` / `dec_do_blk`：解码内部缓冲，长度 `hm_n`；
- `dec_blk`：最终解码输出，与物理码长 `blk_len` 对齐。

### 1.2 矩阵分块：A/B/C/D/E/F/G 结构

在读取 H 矩阵且完成系统分块（通过 `ldpc_gen_gm` 等函数）后，原生 encoder 使用如下子矩阵：

- `qc_a`：上部对系统段 $Z_1$ 的连接（A）；
- `qc_b`：上部对校验段 $Z_2$ 的连接（B）；
- `qc_c`：下部对 $Z_1$ 的连接（C）；
- `qc_d`：下部对 $Z_2$ 的连接（D）；
- `qc_e`：下部对 $Z_3$ 的连接（E）；
- `qc_f1` / `qc_f2`：与 parity $P_2$（记为 `gz`）相关的变换矩阵（对应论文中 F1/F2）；
- `qc_g`：从 $Z_1$ 生成 $P_2$ 的矩阵；
- `qc_fi`：$F^{-1}$，其中 $F = E B + D$。

相应地，编码时的向量为：

- $Z_1$：`enc_di_blk`，长度 `hm_k`；
- $Z_2$：`z2`，长度 `(bm_m - tm_sz - 1) * Z`；
- $Z_3$：`z3`，长度 `tm_sz * Z`；
- $P_2$：`gz`，长度 `Z`；
- 上述中间量 `az1/eaz1/cz1/sumz1/f1z/f2z` 等分别表示 $A Z_1$、$E(AZ_1)$、$C Z_1$、它们的组合以及 F1/F2 与 $P_2$ 的作用结果。

注意：原生 encoder 未直接参与 DQ/KY 的 pad/mask 设计，因此这里的 `Z3` 长度为 `tm_sz * Z`，而 KY 版中则是 `tm_sz*Z - mask_len`。

---

## 2. 编码流程与中间变量长度

以下基于被注释的原始 encoder（`ldpc_codec.cpp` 中的旧版本）及带 dump 的版本，说明关键步骤及每个向量的长度。

### 2.1 信息端 padding：构造 $Z_1$

```c++
vec_copy(usr_blk, enc_di_blk, 0, 0, info_len);
for (int i=0; i<pad_len; i++)
    enc_di_blk[hm_k-pad_len+i] = 0;
```

- `usr_blk`：长度 `info_len`；
- `enc_di_blk`（$Z_1$）：长度 `hm_k`，前 `info_len` 位为用户数据，后 `pad_len` 位为内部补零。

### 2.2 A 分支：`az1 = A Z_1` + 增量验证

```c++
mod2sparse_mulvec(qc_a, enc_di_blk, az1);
// az1: 长度 tm_sz * Z
```

同时，encoder 用一个 per‑edge 的增量实现对 `az1` 做 sanity check（`cn_synd_mem`）并在 `_LDPC_DEBUG_DUMP` 下写入 `din_spare_mul_dump.txt`：

- `vn_flp_sel` / `cn_flp_sel`：长度 `cir_sz`；
- `cn_synd_mem`：长度 `hm_m`；
- dump 的每一行（per edge）包括：
  - `enc_din data = ...`：`vn_flp_sel`，长度 `cir_sz`；
  - `shft_vec = ...`：`cn_flp_sel`，长度 `cir_sz`；
  - `spare mulvec = ...`：`cn_synd_mem[e->row*cir_sz ..]`，长度 `cir_sz`。

这一块主要是验证 $A Z_1$ 的构造是否与 `qc_bm + barrel shift` 的实现等价。

### 2.3 P2 分支：`gz = G Z_1`、`f1z = F1 P2`、`f2z = F2 P2`

```c++
mod2sparse_mulvec(qc_g, enc_di_blk, gz);       // gz: 长度 Z
vec_clr(gz + pad_bit, cir_sz - pad_bit);      // 清零 P2 的 padding 部分

mod2sparse_mulvec(qc_f1, gz, f1z);            // f1z: 长度 tm_sz * Z
mod2sparse_mulvec(qc_f2, gz, f2z);            // f2z: 长度 (bm_m - tm_sz - 1) * Z
```

这里 $P_2$ 是一个长度为 `Z` 的 parity 段，后 `cir_sz - pad_bit` 位被清 0，保留前 `pad_bit` 位有效。

### 2.4 上部合成：$A Z_1 + F_1 P_2$、$E(AZ_1)$

```c++
vec_mod2_add(az1, f1z, az1, tm_sz*cir_sz);   // az1 ← A Z1 + F1 P2
mod2sparse_mulvec(qc_e, az1, eaz1);          // eaz1: 长度 (bm_m - tm_sz - 1) * Z
```

### 2.5 下部合成：$C Z_1 + F_2 P_2$、$E(AZ_1) + C Z_1$

```c++
mod2sparse_mulvec(qc_c, enc_di_blk, cz1);          // cz1: 长度 (bm_m - tm_sz - 1) * Z
vec_mod2_add(cz1, f2z, cz1, (bm_m-tm_sz-1)*cir_sz);// cz1 ← C Z1 + F2 P2
vec_mod2_add(eaz1, cz1, sumz1, (bm_m-tm_sz-1)*cir_sz); // sumz1 = E(AZ1)+C Z1
```

此时 `sumz1` 表示 $F Z_2 = E A Z_1 + C Z_1$ 的右端项（F 的定义被折叠进 `qc_fi` 中）。

### 2.6 求解 $Z_2$：`z2 = F^{-1} sumz1`

```c++
mod2sparse_mulvec(qc_fi, sumz1, z2); // z2: 长度 (bm_m - tm_sz - 1) * Z
```

### 2.7 计算 $B Z_2$、$Z_3 = B Z_2 + A Z_1$

```c++
mod2sparse_mulvec(qc_b, z2, bz2);             // bz2: 长度 tm_sz * Z
vec_mod2_add(az1, bz2, z3, tm_sz*cir_sz);     // z3: 长度 tm_sz * Z
```

### 2.8 组装内部码字与输出

```c++
// encoded data
vec_copy(enc_di_blk, enc_do_blk, 0, 0, hm_k);
vec_copy(z2,         enc_do_blk, 0, hm_k, (bm_m-tm_sz)*cir_sz);
vec_copy(z3,         enc_do_blk, 0, (bm_m-tm_sz-1)*cir_sz, tm_sz*cir_sz);
vec_copy(gz,         enc_do_blk, 0, (bm_n-1)*cir_sz, pad_bit);

// enc_do_blk: 长度 hm_n

// removing 0 padding
vec_copy(enc_do_blk, tx_blk, 0, 0, info_len);
vec_copy(enc_do_blk, tx_blk, hm_k, info_len, hm_m-cir_sz+pad_bit);
// tx_blk: 长度 blk_len，与物理码长一致
```

在原生 encoder 中，`Z_3` 长度为 `tm_sz * Z`，末尾 ISOLATED 的 $P_2$ 段（`gz`）通过最后一次 `vec_copy` 放入 `enc_do_blk` 的最高位区域，长度仅为 `pad_bit`。

---

## 3. 原生 `_LDPC_DEBUG_DUMP`：文件与数据长度

在 `_LDPC_DEBUG_DUMP` 条件下，原生 `ldpc_encoder()` 打开多个文件用于 dump 中间量。以下给出各文件中主要条目及其长度。

### 3.1 `din_spare_mul_dump.txt`（`spfp`）—— per‑edge A·Z1 增量路径

每个非零基矩阵 entry $(row, col, shift)$ 会输出：

- `shft_val=<shift>, row=<row>, col=<col>`；
- `enc_din data = ...`：长度 `cir_sz`；
- `shft_vec = ...`：长度 `cir_sz`；
- `spare mulvec = ...`：长度 `cir_sz`。

上述向量均以 circulant 为单位，一行对应 `cir_sz` 比特，一般配置下（`cir_sz=256`）即为 32Byte（64 个 hex）。

### 3.2 `au_dump.txt`（`aufp`）—— A 分支

两类输出：

- `au calcued:`：`az1 = A Z_1`：
  - 外层 `j = 0 .. tm_sz-1`，每个 `j` 一行；
  - 行长：`cir_sz` 比特 → `cir_sz/4` 个 hex；
  - 对应上部 A·Z1 的每个 circulant。
- `au + f1z calcued:`：`az1 ← A Z_1 + F_1 P_2`：
  - 行数同上；
  - 每行长度仍为 `cir_sz` 比特。

### 3.3 `gu_dump.txt`（`gufp`）—— G/F1/F2 分支

多次输出：

- `g*u calcued:`：`gz = G * Z_1`：
  - 通常只打印若干 block；
  - 行长 `cir_sz` 比特。
- `f1z:`：`f1z = F_1 P_2`，长度 `tm_sz * Z`，按每个 circulant 一行输出；
- `f2z:`：`f2z = F_2 P_2`，长度 `(bm_m - tm_sz - 1)*Z`，同样每个 circulant 一行。

### 3.4 `cu_dump.txt`（`cufp`）—— C 分支

两波输出：

- `cu calcued:`：`cz1 = C Z_1`：
  - 行数：`bm_m - tm_sz - 1`（下半部 circulant 数）；
  - 每行长度：`cir_sz` 比特；
  - 总比特数：`(bm_m - tm_sz - 1) * Z`。
- `cz1 + f2z calcued:`：`cz1 ← C Z_1 + F_2 P_2`：
  - 行数、行长同上。

### 3.5 `eau_dump.txt`（`eaufp`）—— E/B/Z3 路径

按顺序输出：

- `eaz1:`：`eaz1 = E(A Z_1)`：
  - 行数：`bm_m - tm_sz - 1`；
  - 行长：`cir_sz` 比特。
- `sumz1:`：`sumz1 = eaz1 + cz1`：
  - 同上；
  - 总比特数：`(bm_m - tm_sz - 1)*Z`。
- `bz2 calcued:`：`bz2 = B Z_2`：
  - 行数：`tm_sz`；
  - 行长：`cir_sz` 比特；
  - 总比特数：`tm_sz * Z`。
- `z3 calcued:`：`z3 = B Z_2 + A Z_1`：
  - 行数：`tm_sz`；
  - 行长：`cir_sz` 比特。

### 3.6 `fi_in_dump.txt`（`fi_infp`）—— Fi 输入向量

内容：

- `sumz1 calcued:` 接着 `sumz1` 的各行值；
- 行数、行长与 `eaufp` 中的 `sumz1` 相同：
  - 行数：`bm_m - tm_sz - 1`；
  - 行长：`cir_sz` 比特。

### 3.7 `fi_out_dump.txt`（`fi_outfp`）—— Fi 输出向量

- 文本头：`qc_fi calcued:`；
- 输出 `z2 = qc_fi * sumz1`：
  - 行数：`bm_m - tm_sz - 1`；
  - 行长：`cir_sz` 比特；
  - 总比特数：`(bm_m - tm_sz - 1)*Z`。

### 3.8 `fi_val.txt` / `lenc_fi_mem_*.txt` / `lenc_fi_mem_*_hex.txt`—— F⁻¹ 矩阵本身

- `fi_val.txt`：
  - 对若干 row block（例如前 6 个），构造一个 0/1 向量 `ve` 表示每行的非零列；
  - 行长：`mod2sparse_rows(qc_fi)` 比特；
  - 每 4bit→1 hex，按高位到低位输出，行数为选定 block 数；
- `lenc_fi_mem_*`（`_LDPC_DUMP_FI` / `_LDPC_DUMP_FI_HEX`）：
  - 将 `qc_fi` 的某些 row block 以纯 0/1 或 Verilog `32'h...` 的形式输出，方便 RTL 直接装入 BRAM/ROM。

---

## 4. 与 KY 版的对比与互补

- 原生 `ldpc_codec.cpp` encoder：
  - 使用更完整的矩阵分支（含 G/F1/F2/F⁻¹）；
  - `_LDPC_DEBUG_DUMP` 覆盖了从 per‑edge A·Z1 构造、P2 分支、上下部合成到 Fi 运算的全部中间向量；
  - 每个 circulant block 一行，行长 `cir_sz` 比特（常见为 32 Byte、64 hex）。

- KY 版 `ldpc_codec_ky.cpp` encoder：
  - 采用直接矩阵乘法路径，不展开 G/F1/F2 的细节；
  - 目前 `_LDPC_DEBUG_DUMP` 在 KY 版中只覆盖 A/C/E/B/Z2/Z3 及编码前/后码字（`AZ1/CZ1/EAZ1/SUMZ1/Z2/BZ2/Z3/ENC_DO_BLK/TX_BLK`），并拆写到 `au_dump.txt` / `cu_dump.txt` / `eau_dump.txt` 三个文件中；
  - 打印格式同样是每行固定 32Byte（64 hex），但 vector 是 flatten 之后的整体，而不是按 circulant 分行。

在 DV 场景中，可以：

- 用原生 `ldpc_codec.cpp` 的 dump 作为设计/算法开发阶段的基准，逐个模块（A/B/C/E/F/G/F⁻¹）对 RTL 实现做单元验证；
- 在集成 KY 内核时，用 KY 版 `ldpc_encoder` 的 dump 验证“整体路径”在新的 pad/mask 结构下的正确性，必要时再根据原生 encoder 的细粒度 dump 做二次排查。 
