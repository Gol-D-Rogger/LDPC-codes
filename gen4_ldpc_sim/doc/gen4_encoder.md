# Gen4 KY LDPC 编码器原理说明（基于 Padding & Masking 机制）

本文针对 `gen4_ldpc_sim/src/ldpc_codec_ky.cpp` 中的

- `ldpc_packet::ldpc_encoder()`

进行结构化说明，并与 `padding_masking_mechanism.md` 中的 Padding / Masking 机制建立一一对应，帮助理解 Gen4 KY 架构下编码端如何在“逻辑 QC 维度”与“物理比特流长度”之间做映射。

---

## 1. 长度体系与符号约定

为避免混淆，先统一几个关键长度符号（均可在 `ldpc_codec_ky.cpp` / `ldpc_codec.h` 中找到对应成员）：

### 1.1 逻辑维度（QC / H 矩阵内部）

- $bm_m$：基矩阵行数（物理实现中为 $h_m + 1$，多出的 1 行专门用于承载 `pad_bit` 对应的“部分循环块”）。
- $bm_n$：基矩阵列数（为 $h_n + 1$）。
- $cir_{sz}$：循环块大小（代码中的 `cir_sz`，配置中的 `h_sc`）。
- $hm_m = bm_m \cdot cir_{sz}$：H 矩阵逻辑行数，对应**逻辑校验位长度**。
- $hm_n = bm_n \cdot cir_{sz}$：H 矩阵逻辑列数，对应**逻辑码长**。
- $hm_k = hm_n - hm_m$：逻辑信息位长度（含内部 padding）。

这些维度是从 QC 结构出发的“完整维度”，其中最后一个循环块会配合 `pad_bit` / `mask_len` 做部分有效。

### 1.2 物理维度（外部比特流 / 通道）

- $info\_len$：用户真实信息比特数（`ldpc_packet::info_len`）。
- $pad\_len = hm_k - info\_len$：信息侧内部补零长度（逻辑空间内的 padding）。
- $blk\_len$：外部实际传输的码长（`ch_packet::blk_len`，在 Gen4 DQ/KY 配置中为 $H_N = h_n \cdot cir_{sz} + pad\_bit$）。
- $phys\_parity\_len = blk\_len - info\_len$：物理上真正发出的校验位长度。

与 `padding_masking_mechanism.md` 一致，还需要：

- `pad_bit`：最后一个循环块中**有效位数**；
- $drop\_len = pad\_bit$：保留长度；
- $mask\_len = cir_{sz} - pad\_bit$：无效 / 被 Mask 的长度。

因此，有如下关系：

- 逻辑校验长度 $hm_m$ 与物理校验长度 $phys\_parity\_len$ 的差恰为 $mask\_len$ 的某种组合；
- 逻辑码长 $hm_n$ 与物理码长 $blk\_len$ 之间通过 $mask\_len$ 关联：
  - $hm_n = blk\_len + mask\_len$（在典型配置中成立）。

---

## 2. H 矩阵分块与 A/B/C/D/E/F 结构

在 `ldpc_codec_ky.cpp` 中，读取 QC H（`ldpc_rd_phck`）并经 `ldpc_gen_gm()` 进行分块后，构造出若干子矩阵：

- `qc_a`, `qc_b`, `qc_c`, `qc_d`, `qc_e`：对应 H 在行、列方向上的分块；
- `qc_fi`：经由
  - $F = E \cdot B + D$，
  - $F^{-1}$
 计算得到，用于求解校验部分。

在抽象层面，可以认为 H 被划分为：

- 列方向：`[Z_1 | Z_2 | Z_3]`
  - $Z_1$：系统部分（信息 + 内部 padding），长度 $hm_k$；
  - $Z_2$：一部分校验，长度 `(bm_m - tm_sz) \cdot cir_{sz}`；
  - $Z_3$：剩余校验，长度 $tm\_sz \cdot cir_{sz} - mask\_len$。
- 行方向：`[上部 tm\_sz \cdot cir_{sz} 行 | 下部 (bm_m - tm_sz) \cdot cir_{sz} 行]`。

其中最后一个循环块的 `mask_len` 部分被“截断”（Partial Circulant），即只对 $pad\_bit$ 个位置有约束，其余 $mask\_len$ 位被 Mask 掉。这一点在 `padding_masking_mechanism.md` 中已有详细说明。

上述分块保证了：

- 在逻辑维度 $hm_n$ 上，H 仍然是一个完整的校验矩阵；
- 但在实际编码输出中，最后 $mask\_len$ 个“校验位”不会被物理发送，而是通过 mask/pad 机制进行处理。

---

## 3. `ldpc_encoder()` 的总体目标

`ldpc_packet::ldpc_encoder()` 的核心目标可以概括为：

1. 在逻辑空间上，根据 QC H 的结构，为给定的信息比特 $Z_1$ 构造一个满足 $H \cdot cw^T = 0$ 的完整码字 `$cw`，其内部表示为：
   - $cw = [Z_1, Z_2, Z_3]$，长度 $hm_n$；
2. 再将该逻辑码字映射到外部物理比特流 `tx_blk`：
   - 前 $info\_len$ 位：直接取 $Z_1$ 的前 $info\_len$ 位；
   - 后 $phys\_parity\_len$ 位：从 `[Z_2, Z_3]` 中截取物理存在的校验位；
   - 舍弃最后 $mask\_len$ 个仅存在于逻辑空间的“掩码校验位”，不在物理链路上传输。

这与解码端的做法形成镜像：解码时再根据 `pad_len`、`mask_len` 把缺失的逻辑位置补成“强 0”或最大 LLR，保证 syndrome 计算正确。

---

## 4. 编码流程逐步解析

以下基于 `ldpc_packet::ldpc_encoder()` 源码，从上到下分步说明。

### 4.1 临时向量的分配与维度

```c++
char *az1, *eaz1, *cz1, *sumz1, *z2, *bz2, *z3;

az1 = (char *)calloc(tm_sz*cir_sz - mask_len, sizeof(*az1));
eaz1 = (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*eaz1));
cz1  = (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*cz1));
sumz1= (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*sumz1));
z2   = (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*z2));
bz2  = (char *)calloc(tm_sz*cir_sz - mask_len, sizeof(*bz2));
z3   = (char *)calloc(tm_sz*cir_sz - mask_len, sizeof(*z3));
```

- $Z_1$：内部编码输入 `enc_di_blk`，长度 $hm_k$；
- $A Z_1$ 存放在 `az1`，长度 $tm\_sz \cdot cir\_sz - mask\_len$，对应 H 上部与 $Z_1$ 的连接；
- $E A Z_1$ 存放在 `eaz1`，长度 $(bm_m - tm_sz) \cdot cir\_sz$；
- $C Z_1$ 预留在 `cz1`；
- `sumz1` 存放 $E A Z_1 + C Z_1$；
- $Z_2$ 存在 `z2`，长度 $(bm_m - tm_sz) \cdot cir\_sz$；
- $B Z_2$ 在 `bz2` 中；
- $Z_3$ 在 `z3` 中，长度 $tm\_sz \cdot cir\_sz - mask\_len$。

这些维度与 `ldpc_gen_gm()` 中对 H 的分块完全对齐。

### 4.2 构造系统部分 $Z_1$：信息 + 内部 padding

```c++
// padding 0s
vec_copy(enc_di_blk, usr_blk, 0, 0, info_len);
for (int i=0; i<pad_len; i++)
    enc_di_blk[hm_k-pad_len+i] = 0;
```

- `enc_di_blk` 即内部逻辑空间中的 $Z_1$；
- 前 $info\_len$ 比特直接拷贝用户数据 `usr_blk`；
- 尾部 $pad\_len = hm_k - info\_len$ 位填 0，用于“补足逻辑维度”；
- 这样得到的 $Z_1$ 长度为 $hm_k$，与 H 的列分块匹配。

结合 `padding_masking_mechanism.md` 来看，这一步是在**逻辑空间**中把信息扩展到完整的 QC 信息维度 $hm_k$，物理端仍然只会发送前 $info\_len$ 位。

### 4.3 计算 $A Z_1$ 与 $E A Z_1$

```c++
// A*Z1
mod2sparse_mulvec(qc_a, enc_di_blk, az1);
// E*(A*Z1)
mod2sparse_mulvec(qc_e, az1, eaz1);
```

含义：

- 上半部分校验方程中，子矩阵 $A$ 作用在 $Z_1$ 上，得到 `az1`；
- 再通过 $E$ 映射，为构造 $F = E B + D$ 的右端项提前计算 $E A Z_1$。

### 4.4 理论上的 $C Z_1$ 与 $Z_2$ 求解

```c++
// C*Z1
// E*(A*Z1)+C*Z1
vec_mod2_add(eaz1, cz1, sumz1, (bm_m - tm_sz) * cir_sz);
// Z2 = F_inv*[E*(A*Z1)+C*Z1]
mod2sparse_mulvec(qc_fi, sumz1, z2);
```

按 QC-LDPC 系统编码的标准推导，下半部分方程大致为：

$$
C Z_1 + D Z_2 = 0
$$

再配合上半部分及 $E$ 的变换，可以写成：

$$
F Z_2 = E A Z_1 + C Z_1,\quad F = E B + D
$$

因此：

$$
Z_2 = F^{-1} (E A Z_1 + C Z_1)
$$

代码中的：

- `sumz1 = eaz1 XOR cz1` 对应 $E A Z_1 + C Z_1$；
- `z2 = qc_fi * sumz1` 即 $Z_2 = F^{-1} (E A Z_1 + C Z_1)$。

需要注意的是：当前源码中并未显式调用 `mod2sparse_mulvec(qc_c, enc_di_blk, cz1)`，因此 `cz1` 初始为全 0。对当前使用的矩阵，这可能是结构上允许的简化（例如 `C` 部分被吸收到 $E$ 相关结构中），但从理论上理解时应当记住：**概念上存在 $C Z_1$ 这一项**。

### 4.5 计算 $B Z_2$，得到 $Z_3$

```c++
// B*Z2
mod2sparse_mulvec(qc_b, z2, bz2);
// Z3 = BZ2 +AZ1
vec_mod2_add(az1, bz2, z3, tm_sz * cir_sz - mask_len);
```

上半部分校验方程可理解为：

$$
A Z_1 + B Z_2 + T Z_3 = 0
$$

在 Gen4 结构中，`T` 接近单位矩阵，但最后一个循环块由于 `mask_len` 的存在是“缺了一截”的单位块，这就是所谓的 Partial Circulant。

在适当消元后，可以推得：

$$
Z_3 = B Z_2 + A Z_1
$$

因为运算在 GF(2) 上，`+` 即为异或，多数实现中直接写成 `vec_mod2_add(az1, bz2, z3, ...)`。

这里的 `z3` 长度为 $tm\_sz \cdot cir_{sz} - mask\_len$，恰好体现了 `padding_masking_mechanism.md` 中“最后一个循环块只保留 $pad\_bit$ 有效位、抛弃 $mask\_len$ 无效位”的思想。

### 4.6 组装内部完整码字 enc\_do\_blk

```c++
// encoded data
vec_copy(enc_do_blk, enc_di_blk, 0, 0, hm_k);
vec_copy(enc_do_blk, z2, hm_k, 0, (bm_m - tm_sz) * cir_sz);
// place Z3 at the tail after Z2, length without masked part
vec_copy(enc_do_blk, z3, hm_k + (bm_m - tm_sz) * cir_sz, 0, tm_sz * cir_sz - mask_len);
```

执行之后，`enc_do_blk` 的结构为：

- `[0 .. hm_k - 1]`：$Z_1$，即系统部分（信息 + 内部 padding）；
- `[hm_k .. hm_k + (bm_m - tm_sz) \cdot cir_{sz} - 1]`：$Z_2$；
- `[hm_k + (bm_m - tm_sz) \cdot cir_{sz} .. hm_n - mask\_len - 1]`：$Z_3$。

在逻辑空间 $[0 .. hm_n-1]$ 内，最后 $mask\_len$ 个位置在编码输出中**没有被填入任何有效校验比特**，它们对应的是被 mask 掉的部分 circulant 行列。这与解码端的做法呼应：解码时会对这些位置统一赋“强 0”或“最大 LLR”，不参与实际的约束。

### 4.7 从 enc\_do\_blk 到 tx\_blk：去除 padding 与 mask

```c++
// removing 0 padding
vec_copy(tx_blk, enc_do_blk, 0, 0, info_len);
vec_copy(tx_blk, enc_do_blk, hm_k, info_len, blk_len - info_len);
```

这一步是编码器对外的关键接口，决定了物理链路上的比特流构成：

1. 前 $info\_len$ 位：
   - `tx_blk[0 .. info_len-1] = enc_do_blk[0 .. info_len-1]`；
   - 即直接将 $Z_1$ 的信息部分映射到输出；
   - $pad\_len$ 位（`enc_do_blk[hm_k - pad_len .. hm_k - 1]`）完全不出现在 `tx_blk` 中；
2. 后 $phys\_parity\_len = blk\_len - info\_len$ 位：
   - 从 `enc_do_blk[hm_k]` 开始，连续拷贝 $phys\_parity\_len$ 个比特；
   - 剩余的 $hm_m - phys\_parity\_len = mask\_len$ 位（尤其是最后一块中被 mask 的部分）全部被丢弃；
   - 这对应 `padding_masking_mechanism.md` 中的结论：**被 mask 的尾部校验位不存在于物理码字中，只在逻辑空间存在，用于保证 H 的结构与求解的完备性。**

最终，`tx_blk` 的布局为：

- $[0 .. info\_len-1]$：真实用户信息；
- $[info\_len .. blk\_len-1]$：物理发送的校验位；
- 整个码字长度为 $blk\_len$，适配 SSD / DV 所需的非 $cir_{sz}$ 整数倍长度。

---

## 5. 与解码端 Padding/Masking 的一致性

解码时（`ldpc_packet::ldpc_decoder()`）会做与编码端“相反方向”的映射：

1. 从 `det_blk` 中取出前 $info\_len$ 位，作为 $Z_1$ 的信息部分；
2. 对 $Z_1$ 中的内部 padding 区（长度 $pad\_len$）填“强 0”或最大 LLR；
3. 从 `det_blk` 中取出 $phys\_parity\_len$ 个校验比特，填入 $Z_2/Z_3$ 对应的物理位置；
4. 对最后 $mask\_len$ 个“仅逻辑存在”的校验位置，统一填“强 0”或最大 LLR；
5. 使用与编码相同的 H、A/B/C/D/E/F 结构和 mask 规则计算 syndrome 与迭代更新。

这样，编码与解码在“逻辑空间”上对 $Z_1/Z_2/Z_3$ 的处理完全一致，通过 `pad_len` / `mask_len` 等参数，把非整数循环块长度的码字嵌入 QC 结构中，从而在不改变 QC-LDPC 数学本质的前提下，支持精细的 `pad_bit` 粒度。

---

## 6. 小结

- Gen4 KY 编码器 `ldpc_encoder()` 的本质，是在扩展后的 QC 逻辑维度 $hm_n$ 上，为 $Z_1$ 构造满足 $H \cdot cw^T = 0$ 的完整码字 $cw = [Z_1, Z_2, Z_3]$；
- 通过 `pad_len`，在信息端实现内部 padding，使 $Z_1$ 适配 $hm_k$ 这一 QC 结构长度；
- 通过 `mask_len`（以及 `mask_matrix`），在校验端实现部分 circulant 的裁剪：$Z_3$ 长度为 $tm\_sz \cdot cir_{sz} - mask\_len$，最后 $mask\_len$ 个校验位置仅存在于逻辑空间；
- 通过 `tx_blk` 的映射，将逻辑码字投影到物理码长 $blk\_len$ 上，只发送 $info\_len$ 信息位和 $phys\_parity\_len$ 校验位；
- 解码端再利用相同的 Padding & Masking 机制，将物理接收序列嵌回逻辑 QC 空间，保证 syndrome 与迭代算法正确运行。

从工程角度看，这一机制实现了：

- **QC 结构与任意 `pad_bit` 粒度码长之间的解耦**：设计者可以在 QC 维度上使用标准工具生成 H 矩阵，再通过 `pad_bit` / `mask_len` 调整物理码长；
- **对 DV / 控制器友好的接口**：外部只需关心 $info\_len$ 和 $blk\_len$，内部的 $hm_m/hm_n/hm_k/mask\_len$ 完全由 C-Model 和 KY 内核封装处理。

---

## 7. `_LDPC_DEBUG_DUMP` 编码中间数据导出（DV 对接）

Gen4 KY 版本的 `ldpc_encoder()` 可以通过 `_LDPC_DEBUG_DUMP` 宏，将编码流程中的关键中间量按 QC block 的 16 进制形式 dump 到文件。这些文件主要用于 DV/RTL 验证对照，而不是面向普通仿真用户的日志。

### 7.1 文件概览与路径

在 `_LDPC_DEBUG_DUMP` 条件下，编码器在当前工作目录下打开以下文件（路径固定为 `./output`）：

- `au_dump.txt`：上部路径 $A Z_1$ 及 $A Z_1 + F_1 P_2$；
- `cu_dump.txt`：下部路径 $C Z_1$ 及 $C Z_1 + F_2 P_2$；
- `gu_dump.txt`：$G Z_1$、$F_1 P_2$、$F_2 P_2$ 等与 parity $P_2$ 相关路径；
- `eau_dump.txt`：$E(A Z_1)$、$E(A Z_1) + C Z_1$、$B Z_2$、$Z_3$；
- `fi_in_dump.txt`：传入 $F^{-1}$ 的向量 $sumz1 = E(A Z_1) + C Z_1$；
- `fi_out_dump.txt`：$Z_2 = F^{-1} sumz1$；
- `fi_val.txt`：$F^{-1}$ 矩阵 `qc_fi` 的若干行的 0/1 pattern；
- `din_spare_mul_dump.txt`：基于 `qc_bm` + barrel shifter 的 $A Z_1$ 增量实现相关向量。

在 `_LDPC_DUMP_FI` / `_LDPC_DUMP_FI_HEX` 条件下，还会在 `./fi_mtx` 中输出：

- `lenc_fi_mem_%dx%dex%d_w7_dense6_QC_H.txt`：`qc_fi` 的 bit‑pattern（纯 0/1 串）；
- `lenc_fi_mem_%dx%dex%d_w7_dense6_QC_H_hex.txt`：`qc_fi` 的 bit‑pattern，以 Verilog `32'h...` 形式输出，方便直接作为 RTL ROM/BRAM 初始值。

所有这些文件都采用“每个 circulant 一行，每 4bit→1 个 hex 字，从高 bit 到低 bit 组合”的编码方式，便于 RTL 仿真侧按 QC block 对照。

### 7.2 各文件具体内容

**1）`din_spare_mul_dump.txt`（`spfp`）—— A·Z1 增量实现路径**

在用 `qc_bm` + barrel shifter 增量计算 $A Z_1$ 时，针对每个非零基矩阵元素 `e`（`e->row, e->col, e->shift`）输出：

- 行首：`shft_val=<shift>, row=<row>, col=<col>`；
- `enc_din data = ...`：当前列对应的 `vn_flp_sel`（VN 域的 Z1 子块，长度为 `cir_sz`），按每 4bit→1 hex 打印；
- `shft_vec = ...`：`vec_shift(vn_flp_sel, cn_flp_sel, ...)` 的结果，即 CN 域视角的翻转模式；
- `spare mulvec = ...`：`cn_synd_mem[e->row*cir_sz ..]` 更新后的行 syndrome 片段，同样按 4bit→hex 输出。

用途：验证“`qc_bm` + barrel shift + 累加”这一实现是否与直接 `qc_a * Z_1` 等价，是 encoder 中 A·Z1 生成路径的细节 dump。

**2）`au_dump.txt`（`aufp`）—— 上部 A 分支**

两次输出：

- `au calcued:`：`az1 = A * Z_1` 的结果；
- `au + f1z calcued:`：`az1 ← az1 + f1z`，即 $A Z_1 + F_1 P_2$。

每行对应一个 circulant 行块（上部 `tm_sz` 行），一行包含 `cir_sz/4` 个 hex 字，按 bit 从高到低组合。

**3）`gu_dump.txt`（`gufp`）—— G/F1/F2 分支**

多次复用同一个文件句柄：

- `g*u calcued:`：`gz = G * Z_1`；
- `f1z:`：`f1z = F_1 * P_2`；
- `f2z:`：`f2z = F_2 * P_2`。

用途：覆盖所有与 parity $P_2$ 相关的路径，RTL 中 G/F1/F2 模块可以以此为黄金参考。

**4）`cu_dump.txt`（`cufp`）—— 下部 C 分支**

两波输出：

- `cu calcued:`：`cz1 = C * Z_1`；
- `cz1 + f2z calcued:`：`cz1 ← cz1 + f2z`，即 $C Z_1 + F_2 P_2$。

与 `au_dump` 对应，但针对 H 的下半部分（`bm_m - tm_sz - 1` 行）。

**5）`eau_dump.txt`（`eaufp`）—— E/B/Z3 路径**

按顺序输出：

- `eaz1:`：`eaz1 = E * (A Z_1)`；
- `sumz1:`：`sumz1 = eaz1 + cz1`（即 $E A Z_1 + C Z_1$）；
- `bz2 calcued:`：`bz2 = B * Z_2`；
- `z3 calcued:`：`z3 = B Z_2 + A Z_1`（即 $Z_3$）。

用途：覆盖 encoder 中所有通过 E、B 构造的上/下部等效方程路径。

**6）`fi_in_dump.txt`（`fi_infp`）—— Fi 输入向量**

- 文本头：`sumz1 calcued:`；
- 输出 `sumz1 = E A Z_1 + C Z_1`，内容与 `eaufp` 中的 `sumz1` 相同，但单独成文件。

用途：明确标记 $F^{-1}$ 运算的输入向量，DV 可以直接将此向量喂给 RTL 中的 Fi 模块。

**7）`fi_out_dump.txt`（`fi_outfp`）—— Fi 输出向量**

- 文本头：`qc_fi calcued:`；
- 输出 `z2 = qc_fi * sumz1` 的结果，即 $Z_2$。

用途：Fi 模块输出的黄金结果，用于 RTL bit‑exact 对照。

**8）`fi_val.txt`（`fi_val_fp`）以及 `_LDPC_DUMP_FI*` 相关文件—— $F^{-1}$ 矩阵本身的 bit‑pattern**

`fi_val.txt` 中：

- 对若干选定的 row block（例如前 6 个 block）：
  - 构造一个 0/1 向量 `ve`，表示 `qc_fi` 对应行所有非零列的位置；
  - 然后按 `(bit3..bit0)`→1 hex 的方式输出。

`_LDPC_DUMP_FI` / `_LDPC_DUMP_FI_HEX` 变体则将同样信息：

- 以纯 0/1 流形式输出到 `lenc_fi_mem_*.txt`；
- 或以 Verilog `32'h...` 格式输出到 `lenc_fi_mem_*_hex.txt`。

用途：用于生成 RTL 中 $F^{-1}$ 的 ROM/BRAM 初始化内容，或检查矩阵装载是否正确。

### 7.3 对 DV/RTL 验证的意义

综合上述 `_LDPC_DEBUG_DUMP` 系列输出：

- 提供了 `Z_1 → A Z_1 → E A Z_1 → sumz1 → Z_2 → B Z_2 → Z_3` 的全路径向量级黄金数据；
- 暴露了 encoder 内所有关键矩阵乘法（G、F1、F2、C、E、B、F_i）的中间结果；
- 给出了 $F^{-1}$ 矩阵本身的 bit‑pattern 及其 ROM/BRAM 友好格式；
- 通过 `din_spare_mul_dump` 对 QC barrel‑shifter + pad/mask 相关逻辑提供了细粒度 bit trace。

因此，Gen4 KY 编码器的 `_LDPC_DEBUG_DUMP` 宏可以看作是：

- 为 DV/RTL 场景专门设计的一组“内部观测点导出接口”，
- 使硬件实现可以逐级（每个 QC block、每条数据路径）与 C‑Model 做 bit‑exact 对照，
- 尤其在存在 `pad_bit` / `mask_len` / `drop_col_bit_map` 等复杂裁剪逻辑时，有助于快速定位矩阵装载、部分 circulant、padding & masking 处理是否与软件一致。

### 7.4 KY 版 dump 的实现方式与当前覆盖范围

对于 KY 版 `ldpc_encoder()`（`ldpc_codec_ky.cpp`），当前实现的 `_LDPC_DEBUG_DUMP` 是在保持算法简洁前提下，对原 `ldpc_codec.cpp` dump 体系的一个子集适配：

- KY 版 encoder 内部使用的是“整体矩阵乘法”（`mod2sparse_mulvec`）实现 $A Z_1$ / $E(A Z_1)$ / $F^{-1}(\cdot)$ 等，而非原版中 per‑edge 的增量构造；
- 因此，KY 版更自然的选择是对每个关键向量（而不是每个 edge）做 block 级别 dump。

当前在 `_LDPC_DEBUG_DUMP` 下，KY 版 encoder 做了如下处理：

1. 在函数开头按需打开三个文件指针：

   ```c
   FILE *aufp  = fopen("./output/au_dump.txt",  "w");
   FILE *cufp  = fopen("./output/cu_dump.txt",  "w");
   FILE *eaufp = fopen("./output/eau_dump.txt", "w");
   ```

2. 定义一个通用宏 `KY_PRINT_VEC_HEX`，用于将任意长度的 `char` 向量以 4bit→1 hex 的形式打印到指定文件：

   ```c
   #define KY_PRINT_VEC_HEX(fp, label, vec, len)                          \
       do {                                                               \
           if ((fp) != NULL) {                                            \
               int _len = (len);                                          \
               fprintf((fp), "%s: ", (label));                            \
               int _nib = (_len + 3) / 4;                                 \
               for (int _i = _nib - 1; _i >= 0; --_i) {                   \
                   int _tmp = 0;                                          \
                   for (int _k = 3; _k >= 0; --_k) {                      \
                       int _idx = _i * 4 + _k;                            \
                       _tmp <<= 1;                                        \
                       if (_idx < _len)                                   \
                           _tmp += ((vec)[_idx] & 1);                     \
                   }                                                      \
                   fprintf((fp), "%x", _tmp);                             \
               }                                                          \
               fprintf((fp), "\n");                                       \
           }                                                              \
       } while (0)
   ```

   - 其中 `len` 以 bit 为单位，`vec[i]` 必须是 0/1；
   - 输出格式为：`LABEL: <hex串>\n`，hex 串最左侧对应最高索引的 4 个 bit。

3. 在 KY 版编码流程的关键节点调用该宏，将中间向量分配到对应文件：

   - 在完成 `A Z_1` 后：

     ```c
     mod2sparse_mulvec(qc_a, enc_di_blk, az1);
     KY_PRINT_VEC_HEX(aufp, "AZ1", az1, tm_sz*cir_sz - mask_len);
     ```

   - 在完成 `E(A Z_1)` 后：

     ```c
     mod2sparse_mulvec(qc_e, az1, eaz1);
     KY_PRINT_VEC_HEX(eaufp, "EAZ1", eaz1, (bm_m - tm_sz)*cir_sz);
     ```

   - 为了完整覆盖 C 分支，KY 版编码器在原有基础上补算了：

     ```c
     mod2sparse_mulvec(qc_c, enc_di_blk, cz1);
     KY_PRINT_VEC_HEX(cufp, "CZ1", cz1, (bm_m - tm_sz)*cir_sz);
     ```

   - 然后按顺序 dump：
     - `SUMZ1 = EAZ1 + CZ1` → `eau_dump.txt` 中的 `"SUMZ1"`；
     - `Z2 = F^{-1} SUMZ1` → `"Z2"`；
     - `BZ2 = B Z2` → `"BZ2"`；
     - `Z3 = B Z2 + A Z1` → `"Z3"`；
   - 最后对组装好的码字和物理输出码字也做 dump：

     ```c
     KY_PRINT_VEC_HEX(aufp, "ENC_DO_BLK",
                      enc_do_blk, hm_k + (bm_m - tm_sz)*cir_sz + (tm_sz*cir_sz - mask_len));
     KY_PRINT_VEC_HEX(aufp, "TX_BLK", tx_blk, blk_len);
     ```

4. 在函数结尾关闭文件并释放宏定义：

   ```c
   if (aufp)  fclose(aufp);
   if (cufp)  fclose(cufp);
   if (eaufp) fclose(eaufp);
   #undef KY_PRINT_VEC_HEX
   ```

### 7.5 KY 版 dump 与原 `ldpc_codec.cpp` dump 的对应关系

与原 `ldpc_codec.cpp` 中更细粒度的 dump 相比，KY 版目前覆盖的是其子集：

- 已对齐的部分：
  - `au_dump.txt`：KY 版 dump 了 `AZ1`、`ENC_DO_BLK`、`TX_BLK`，对应原版的 `A*Z1`、编码输出、去 padding 输出；
  - `cu_dump.txt`：KY 版补算并 dump 了 `CZ1 = C*Z1`，对应原版 `cu calcued`；
  - `eau_dump.txt`：KY 版 dump 了 `EAZ1`、`SUMZ1`、`Z2`、`BZ2`、`Z3`，与原版 `eaz1`、`sumz1`、`Z2`、`BZ2`、`Z3` 的含义一致。

- 当前未在 KY 版实现的部分：
  - `gu_dump.txt` 中的 `G*Z1`、`F1*P2`、`F2*P2`；
  - `fi_in_dump.txt` / `fi_out_dump.txt` / `fi_val.txt` 中针对 $F^{-1}$ 的更细致矩阵/向量 dump；
  - `din_spare_mul_dump.txt` 中 per‑edge 的 barrel shift + syndrome 增量更新轨迹。

原因在于：KY 版 encoder 使用的是直接矩阵乘法路径，暂未在 C 侧展开 per‑edge 或 G/F1/F2 的实现。若需要完全对齐原版的 dump 粒度，可以进一步在 KY 版本中增加相应的中间向量和矩阵乘法路径，再复用 `KY_PRINT_VEC_HEX` 宏将其拆写到与原版相同的文件名中。

### 7.6 KY 版各 dump 文件数据长度分析

为方便 DV 在解析 KY 版 dump 时做边界检查，本小节给出三个文件中每个标签对应向量的 bit 长度。记：

- $Z = cir\_{sz}$；
- 上部行块数：$tm\_sz$；
- 基矩阵行数：$bm_m$；
- 逻辑信息长度：$hm_k$；
- 逻辑校验长度：$hm_m = bm_m \cdot Z$；
- 物理码长：$blk\_len$；
- Mask 长度：$mask\_len = Z - pad\_bit$。

#### 7.6.1 `au_dump.txt`（`aufp`）—— 上部 A 分支与最终码字

KY 版在 `au_dump.txt` 中写入：

- `USR_BLK`：
  - 来源：`usr_blk`；
  - 长度：`info_len` bit。
- `ENC_DI_BLK`：
  - 来源：`enc_di_blk`（padding 后的 `$Z_1$`）；
  - 长度：`hm_k` bit。
- `AZ1`：
  - 来源：`az1 = A * Z_1`；
  - 长度：`tm_sz \cdot Z - mask\_len` bit；
  - 对应上部 `tm_sz` 个 circulant 中，最后一块去掉 `mask\_len` 的 partial circulant 部分。
- `ENC_DO_BLK`：
  - 来源：`enc_do_blk`（内部完整码字 `[Z_1, Z_2, Z_3]`，不含尾部被 mask 的逻辑校验位）；
  - dump 长度：
    $$
    hm_k + (bm_m - tm\_sz)\,Z + (tm\_sz\,Z - mask\_len)
    = hm_k + hm_m - mask\_len
    $$
  - 即逻辑码长 `hm_n = hm_k + hm_m` 减去 `mask\_len` 个被 mask 的尾部校验位。
- `TX_BLK`：
  - 来源：`tx_blk`（物理发送码字）；
  - 长度：`blk_len` bit（`info_len + phys\_parity\_len`）。

每个标签下的向量由宏按 32 Byte（64 个 hex）一行输出：

- 每行 bit 数：`32 * 8 = 256`；
- 每行 hex 数：`64`；
- 不足 256bit 的最后一行自动截断。

#### 7.6.2 `cu_dump.txt`（`cufp`）—— 下部 C 分支

KY 版在 `cu_dump.txt` 中写入：

- `CZ1`：
  - 来源：`cz1 = C * Z_1`；
  - 长度：`(bm_m - tm\_sz) \cdot Z` bit；
  - 对应 H 下半部分对 `$Z_1$` 的全部连接。

同样按 32 Byte 为一行输出。

#### 7.6.3 `eau_dump.txt`（`eaufp`）—— E/B/Z2/Z3 路径

KY 版在 `eau_dump.txt` 中写入：

- `EAZ1`：
  - 来源：`eaz1 = E * (A Z_1)`；
  - 长度：`(bm_m - tm\_sz) \cdot Z` bit。
- `SUMZ1`：
  - 来源：`sumz1 = eaz1 + cz1 = E A Z_1 + C Z_1`；
  - 长度：`(bm_m - tm\_sz) \cdot Z` bit；
  - 即传入 `$F^{-1}$` 的右端向量。
- `Z2`：
  - 来源：`z2 = qc_fi * sumz1`；
  - 长度：`(bm_m - tm\_sz) \cdot Z` bit。
- `BZ2`：
  - 来源：`bz2 = B Z_2`；
  - 长度：`tm\_sz \cdot Z - mask\_len` bit；
  - 对应上部对 `$Z_2$` 的贡献，最后一个 circulant 去掉 `mask\_len` 位。
- `Z3`：
  - 来源：`z3 = B Z_2 + A Z_1`；
  - 长度：`tm\_sz \cdot Z - mask\_len` bit；
  - 即 `$Z_3$`，与 `AZ1` / `BZ2` 长度一致。

这些向量同样按“32 Byte 一行（64 个 hex）”的方式输出，对 DV 来说，只要事先知道各标签的 bit 长度，就可以精确切分每一行对应的 bit 区间。  
