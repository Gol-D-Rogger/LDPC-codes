# DVC C-Model 中 2Byte 对齐与 LSB/MSB 处理问题分析报告

## 1. 背景与现象

在 Gen4 DVC C-Model (`DVCtrans/src/ldpc_c_model_gen4.c`) 适配过程中，出现如下现象：

- 配置中使用 `pad_bit = 16`、`48` 等时（即非 32bit 整数倍的补位），编码输出和经 `ch_err_inj` 后的硬判输出在 **最后一个 32bit word** 上出现：
  - 期望值（按照 SV 侧 LSB-first 视角）：`0000E70C`
  - 实际值（C-Model 打包结果）：`E70C0000` 或相反
- 当 `pad_bit = 32` 时，问题消失（`blk_len` 为 32 的整数倍）。

同时，在 `ch_err_inj` 中，`tx_data` 和 `rx_data` 对比时也出现类似“最后 16bit 清零”的现象：当 `blk_len % 32 = 16` 时，最后一个 32bit word 的一部分位为 0。

这些现象集中在“最后一个非整 32bit 的 word”，提示问题与：

- 32bit 对齐（2Byte/4Byte 对齐）
- C-Model 与 SV 之间的 bit 顺序约定（MSB-first vs LSB-first）

有关，而非 LDPC 编解码数学逻辑本身。

---

## 2. 现有实现的 bit 打包/解包约定

### 2.1 C-Model 侧：统一的 MSB-first 约定

在 `ldpc_c_model_gen4.c` 中，所有 DPI 接口的打包/解包逻辑都使用类似模式：

```c
// 解包例子（ldpc_enc, ldpc_dec, ch_err_inj 等）
if (j == 0)
    tmp = (unsigned int)(data[i/32]);
bit = (tmp >> (31 - j)) & 1;

// 打包例子
tmp = tmp * 2 + bit;       // 等价于 tmp = (tmp << 1) | bit
if (j == 0 || i == blk_len - 1)
    word = tmp;
```

这等价于：

- **解包**：假定 SV 侧把第一个 bit 放在 32bit word 的 bit31（MSB），然后 bit30, bit29 ...直到 bit0。
- **打包**：从 bit 流头开始连续左移，把第 0 个 bit 放在 word 的最高位，第 31 个 bit 放在最低位。

因此，C-Model 内部是自洽的 MSB-first：所有入口/出口都按“MSB-first”打包/解包。

### 2.2 SV 侧：LSB-first 视角

当前工程中，SV 侧按照 **LSB-first** 的方式解析 32bit word，即从 bit0 开始数“第一个 bit”。因此：

- 对完整 32bit word，如果只做码字级比较（相同打包方式），差异不显见。
- 对最后一个 **只包含 16 个有效 bit** 的 word，MSB/LSB 差异会直接表现为：
  - C 侧看：有效 bit 在高 16 位，高位或低位填 0 无所谓。
  - SV 侧看：有效 bit 在低 16 位才是“正常”，出现在高 16 位就被认为“位置错了”。

---

## 3. ch_err_inj 中的 32bit 打包问题

`ch_err_inj` 的输出打包逻辑如下：

```c
// data out
j = 0;
k = 0;
tmp = 0;
for (i = 0; i < sim_pckt->blk_len; i++)
{
    tmp = tmp*2 + sim_pckt->det_blk[i];
    j++;
    j = j%32;

    if ((j==0) || (i==sim_pckt->blk_len-1))
    {
        det_data[k] = tmp;
        tmp = 0;
        k++;
    }
}
```

这里 `sim_pckt->blk_len` 是码字实际有效 bit 数（Info + Parity）。对于两种情况：

- 若 `blk_len % 32 == 0`：每次 `tmp` 都被 32 个 bit 填满，没有“空位”，`det_data` 的每一 word 都是“满字”。
- 若 `blk_len % 32 == 16`：最后一次循环只写入 16 个 bit，`tmp` 的其余 16 位保持初始值 0。因此最后一个 word 的一半是有效 bit，一半是填充 0。

对 C-Model（MSB-first）来说，这是合乎逻辑的“尾部对齐”：**码字之外的 bit 是无定义/填充 0**。

但在 SV（LSB-first + 32bit word 比对）视角下，如果简单地对比整 32bit word，会误以为“最后 16bit 被清零”，特别是在 pad16/pad48 下：

- pad32 使得 `blk_len` 恰好是 32 的倍数 → 不出现半 word，提高了问题“隐身率”。
- pad16/48 导致 `blk_len % 32 = 16` → 尾 word 半数是 0，问题暴露。

严格来说，从“有效 bit 数 = blk_len”的角度看，这些额外 16bit 填 0 不属于码字，只是对齐填充，但如果 SV 端按整 word 比较，确实会观察到：

> `tx_data` 最后 word 有随机数据，`err_data`/`dec_data` 最后 16bit 为 0。

---

## 4. ldpc_enc / ldpc_dec 中的最后半 word 问题

类似的 32bit 打包逻辑也出现在：

- `ldpc_enc` 输出编码码字到 `enc_data_sv`
- `ldpc_dec` 输出解码结果到 `dec_data_sv`

都使用相同的 MSB-first 模式：

```c
for (i = 0; i < sim_pckt->blk_len; i++)
{
    tmp = tmp*2 + sim_pckt->tx_blk[i];
    ...
    if ((j==0) || (i==sim_pckt->blk_len-1))
        enc_data[k] = tmp;
}
```

因此，当 `blk_len` 不是 32 的整数倍时，这里也会出现和 `ch_err_inj` 一样的尾 word 对齐问题。  
在 pad16 这种 2Byte 对齐但非 4Byte 对齐的场景下，最后 16bit 的位置（高 16 或低 16）与 SV 的 LSB-first 期望不一致，就直接表现为：

- 例如你期望最后是 `0000E70C`（LSB-first 看低 16bit 是 E70C），但 C-Model 打包的 32bit 值在 SV 视角下显示为 `E70C0000`。

这不是 LDPC pad 算法的问题，而是 **word 内 bit 顺序 + 最后半 word 对齐策略** 与 SV 不一致。

---

## 5. 问题本质小结

1. C-Model 当前是 **MSB-first 打包/解包**，SV 是 **LSB-first 视角**；
2. 对齐 32bit 时，针对最后一个非整 32bit 的 word，C-Model 将有效 bit 放在一侧，另一侧填 0；
3. pad32bit 时 `blk_len` 往往恰好整除 32，尾 word 是“满字”，问题不显；
4. pad16/48 使得 `blk_len % 32 = 16`，尾 word 的一半是填充 0，LSB-first 视角下就觉得“最后 16bit 被清零”或“有效 bit 在错误一边”。

---

## 6. 实际修正的思路与折中方案

从理论上讲，最干净的做法有两种：

1. **统一改成 LSB-first**：重写所有 DPI 边界的打包/解包（enc/dec/ch_err_inj/sd_err_inj），用 LSB-first；
2. **接受 MSB-first**：在 SV 侧也按 MSB-first 理解 word 内 bit 排布，不再要求“最后 16bit 在低位”之类的视觉惯性。

考虑工程成本和兼容性，目前采用的是一个折中方案：

- 保留 C-Model 内部的 MSB-first 实现（不动 LDPC 内部逻辑）；
- 在 DPI 出口（`ldpc_enc` 和 `ldpc_dec`），对最后一个只含 16bit 的 word 做一个 **2Byte-align 对齐修正**，让 SV 的 LSB-first 视角看到的最后 16bit 落在期待的位置。

### 6.1 在 ldpc_enc 中的修正示意

原打包代码：

```c
j = 0;
k = 0;
tmp = 0;
for (i = 0; i < sim_pckt->blk_len; i++)
{
    tmp = tmp*2 + sim_pckt->tx_blk[i];
    j++;
    j=j%32;

    if ((j==0) || (i==sim_pckt->blk_len-1))
    {
        enc_data[k] = tmp;
        k++;
        tmp = 0;
    }
}
```

为兼容 2Byte-align（例如 pad16/pad48），在循环后增加对最后一个 word 的 special handling，例如：

```c
// 仅当最后只剩 16bit 有效数据时做 2Byte-align 调整
if ((sim_pckt->blk_len % 32) == 16) {
    unsigned int last = (unsigned int)enc_data[k-1];

    // 根据实际情况调整方向：
    // 如当前 last = 0x0000E70C，而 SV LSB-first 希望看到最后 16bit 为 E70C，可用：
    last = last << 16;  // 0x0000E70C -> 0xE70C0000

    enc_data[k-1] = last;
}
```

具体左移还是右移，需要根据 SV/RTL 日志中实际的 bit 落位来决定。实务中观察到：

- 原始编码后最后一 word 的日志中 bit 顺序为 E70C（从 0000 方向读起），  
  而 C-Model 数值为 0x0000E70C，  
  因此采用 `last <<= 16` 将有效 16bit 移到 RTL 期望的半区。

### 6.2 在 ldpc_dec 中的对称修正

`ldpc_dec` 在打包 `dec_blk` 时使用相同的 MSB-first 逻辑，因此需要对最后 half-word 做同样位移，使解码结果在 SV 日志中也对齐：

```c
// 原打包循环
for (i = 0; i < sim_pckt->blk_len; i++)
{
    tmp = tmp*2 + sim_pckt->dec_blk[i];
    ...
    if ((j==0) || (i==sim_pckt->blk_len-1))
    {
        dec_data[k] = tmp;
        tmp = 0;
        k++;
    }
}

// 尾 word 调整（2Byte-align）
if ((sim_pckt->blk_len % 32) == 16) {
    unsigned int last = (unsigned int)dec_data[k-1];
    last = last << 16;   // 与 enc 同步
    dec_data[k-1] = last;
}
```

这样，TX / ERR / DEC 在 SV 的 LSB-first 视角下，最后 16bit 的位置保持一致。

### 6.3 ch_err_inj 的处理建议

`ch_err_inj` 的输出打包目前仍是原始 MSB-first 实现，对于 CLEAN 信道，期望 `rx_data` 与 `tx_data` 一致。要保证在 pad16/48 时最后一 word 也对齐，可采取类似策略：

```c
// 原打包循环后
if ((sim_pckt->blk_len % 32) == 16) {
    unsigned int last = (unsigned int)det_data[k-1];
    last = last << 16;   // 或按实际需要调整方向
    det_data[k-1] = last;
}
```

或者更严格一些：以 `tx_data[w]` 为基线，仅覆盖前 `blk_len` 个 bit，让对齐位继承 TX 的原值，而不是填 0。此方案更精确，但改动略大。

---

## 7. 风险与注意事项

1. **这是接口层的格式修正**  
   - LDPC 内部 `[Info | Pad | Parity]` 的数学逻辑没有变；
   - 仅改变了“最后半个 32bit word 在 C 侧数值表示中的位移方式”，以便匹配 SV 的 LSB-first + 2Byte-align 习惯。

2. **需保持前后对称**  
   - 编码输出（enc）、解码输出（dec）、以及信道硬判输出（err/ch_err_inj）在最后 half-word 上的位移方向必须一致，才能保证 TX/ERR/DEC 对齐；
   - 否则，“编码看起来对了，解码/错误注入仍然错位”。

3. **后续扩展建议**  
   - 若工程后期需要支持多种 bit 排布方式，建议在文档和代码中显式标记“打包约定”（MSB-first/LSB-first），并在 SV+C 两侧统一；
   - 当前这套方案可以视为一个实用的 workaround，解决 Gen4 DVC 在 pad16/48 时 2Byte/4Byte 对齐导致的可见差异。

---

## 8. 结论

- 问题根源在于 **C-Model 内部 MSB-first 打包** 与 **SV 侧 LSB-first + 2Byte-align 期望** 不一致，在 `blk_len % 32 != 0` 时的最后 half-word 上被放大；
- ch_err_inj 与 ldpc_enc/dec 中的 32bit 打包逻辑在数学上是自洽的，但从“观察十六进制 word”的角度可能与 SV 的视觉习惯不符；
- 通过在 `ldpc_enc`、`ldpc_dec`（以及必要时 `ch_err_inj`）打包结束后，对最后只含 16bit 有效数据的 word 做一次统一位移修正，可以在不改变 LDPC 内部算法的前提下，解决 pad16/48 场景下的“最后 16bit 方向错位”和“尾 16bit 清零”的问题。***

## 9. Gen4 KY pad\_bit ≠ 256 时 ldpc\_cleanup 崩溃问题分析与解决流程

> 本节针对 Gen4 KY 版本（`gen4_ldpc_sim/src/ldpc_codec_ky.cpp`）与 DVC C-Model（`DVCtrans/src/ldpc_c_model_gen4.c`）组合时，在 `pad_bit` 非 256 的配置下出现的 `ldpc_cleanup` 崩溃问题，记录完整的长度链路分析与修正步骤。

### 9.1 现象与配置场景

- 配置场景：
  - 使用 Gen4 KY 版本 LDPC 内核；
  - 调用 DVC C-Model：`ldpc_c_model_gen4.c`；
  - 配置参数示例：`h_m = 20`，`h_n = 149`，`h_sc = 256`，`pad_bit = 16`（`pad_bit ≠ 256`）；
  - KY 内核经 `ldpc_config_dq` 计算得到：
    - $bm_m = h_m + 1 = 21$，$bm_n = h_n + 1 = 150$（逻辑维度，含额外一行一列）；
    - $hm_m = bm_m \cdot h_{sc} = 21 \cdot 256 = 5376$；
    - $hm_n = bm_n \cdot h_{sc} = 150 \cdot 256 = 38400$；
    - $hm_k = hm_n - hm_m = 33024$；
    - $mask\_len = h_{sc} - pad\_bit = 240$。
- 通道与 C-Model 外部链路：
  - 实际物理码长（外部见到的块长）为 $blk\_len = H_N = h_n \cdot h_{sc} + pad\_bit = 149 \cdot 256 + 16 = 38160$；
  - 满足：
    - $blk\_len = info\_len + phys\_parity\_len$，
    - $phys\_parity\_len = blk\_len - info\_len = 5136$，
    - $hm_m - phys\_parity\_len = mask\_len = 240$。
- 实际观察到的现象：
  1. 在 `ldpc_pckt_alloc` 和 `ldpc_config` 打印中，`blk_len = 38160`；
  2. 在 `ldpc_decoder` 打印中，开始时 `blk_len = 38160`，中途有代码把 `blk_len` 改为 `38400`，结束时又改回 `38160`；
  3. 在 `ldpc_cleanup` 中打印，`blk_len = 38160`，执行到 `ldpc_clean` 中 `mod2sparse_free(qc_e)` 即崩溃，报 `corrupted size vs. prev_size`；
  4. 当 `pad_bit = 256`（即无“部分 circulant”块）时，不出现崩溃。

初步结论：问题只在 `pad_bit ≠ 256`、存在最后“部分 circulant 块”（多出 $pad\_bit$）时出现，并与 KY 内部 $hm_m$ 与外部 $blk\_len$ 长度不一致高度相关。

### 9.2 长度链路梳理：逻辑维度 vs 物理维度

在 KY+DQ 结构中，可以区分两类长度：

1. **逻辑维度（H 矩阵内部维度）**
   - $bm_m, bm_n$：QC 矩阵块维度（行/列个数），包含额外的一行一列承载 `pad_bit` 对应的“部分 circulant 块”；
   - $hm_m = bm_m \cdot h_{sc}$：逻辑上完整的校验位长度（包括最后一块的 `h_{sc}` 位，其中有 $mask\_len$ 实际被 mask）；
   - $hm_n = bm_n \cdot h_{sc}$：逻辑上完整的码长；
   - $hm_k = hm_n - hm_m$：逻辑上的信息位长度（包括 padding）。

2. **物理维度（外部链路实际传输/缓存长度）**
   - 通道与 C-Model 只处理长度为 $blk\_len = H_N$ 的 bit 流，其中：
     - $H_N = h_n \cdot h_{sc} + pad\_bit$；
     - $H_M = h_m \cdot h_{sc} + pad\_bit$；
   - `sim_pckt->blk_len` 在 DVC C-Model 中被设置为 $blk\_len$，并用于：
     - `tx_blk` / `det_blk` / `dec_blk` 等外部缓存的访问上界；
     - DPI 接口（SV ↔ C）数组长度。

二者之间的关系为：

- $hm_n = blk\_len + mask\_len$；
- $hm_m = H_M + mask\_len$；
- $hm_m - phys\_parity\_len = mask\_len$。

换言之：最后一个 circulant 块在逻辑上有 $h_{sc}$ 位，但物理上只传输 $pad\_bit$ 位，其余 $mask\_len$ 位由 mask/pad 机制处理。

### 9.3 根本问题：以 hm\_m 为界访问 tx\_blk/det\_blk/dec\_blk 导致越界

在原始实现中，存在几个关键问题点：

1. **C-Model 解码入口临时修改 blk\_len**
   - 在 `ldpc_c_model_gen4.c::ldpc_dec` 中曾有代码将：
     - `sim_pckt->blk_len` 临时改为 $hm_k + hm_m = 38400$，以为方便 KY 内核使用；
   - 这违反了“外部物理长度统一使用 $blk\_len$”的原则，导致：
     - KY 内部认为 `tx_blk/det_blk/dec_blk` 长度为 `38400`；
     - 实际 `tx_blk/det_blk/dec_blk` 分配长度仍为 `38160`；
     - 所有使用 `hm_m` 作为长度的拷贝都可能越界。

2. **KY 解码器中以 hm\_m 为界的 vec\_copy**
   - 在 `ldpc_codec_ky.cpp::ldpc_decoder` 中，存在类似逻辑（伪代码）：

     ```c++
     // det_blk → dec_di_blk
     // 期望结构：[Info (info_len) | Pad (pad_len) | Parity (hm_m)]
     vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);

     // dec_do_blk → dec_blk
     vec_copy(dec_do_blk, dec_blk, hm_k, info_len, hm_m);
     ```

   - 然而外部 `det_blk/dec_blk` 实际长度只有 $blk\_len = info\_len + phys\_parity\_len$，而：
     - $info\_len + hm_m = info\_len + phys\_parity\_len + mask\_len = blk\_len + mask\_len > blk\_len$，
   - 导致：
     - 读取 `det_blk` 或写回 `dec_blk` 时，最后 $mask\_len$ 位越界访问（刚好对应最后 circulant 的 mask 区）。

3. **KY 编码器中去 padding 时同样使用 hm\_m**
   - 在 `ldpc_codec_ky.cpp::ldpc_encoder` 中，去 padding 时也存在：

     ```c++
     // removing 0 padding
     vec_copy(tx_blk, enc_do_blk, 0,       0,        info_len);
     vec_copy(tx_blk, enc_do_blk, hm_k,    info_len, hm_m);  // 以 hm_m 拷贝 parity
     ```

   - `tx_blk` 对外同样只有 $blk\_len$ 位有效，以 $hm_m$ 为长度拷贝 parity 时，也会从 `tx_blk` 尾部越界读取 $mask\_len$ 位。

4. **heap corruption 在 ldpc\_clean 中暴露**
   - 上述越界写/读会破坏堆上的结构（特别是 `mod2sparse` 相关结构体及链表），不一定立即崩溃；
   - 实际运行中，问题在 `ldpc_clean` 调用 `mod2sparse_free(qc_e)` 时触发 `malloc` / `free` 堆一致性检查，报：

     > corrupted size vs. prev\_size

   - 当 `pad_bit = 256` 时，`mask_len = 0`，`phys\_parity\_len = hm_m`，上述越界区间长度为 0，因此看起来“一切正常”，实则只是特殊配置掩盖了问题。

### 9.4 初始尝试：只在 C-Model 端修正 blk\_len 的失败

问题暴露后，首先做了两个尝试：

1. **去掉 C-Model 中修改 blk\_len 的逻辑**
   - 将 `ldpc_c_model_gen4.c::ldpc_dec` 中修改 `sim_pckt->blk_len = hm_k + hm_m` 的代码删除，使：
     - `blk_len` 在 C-Model 全流程中始终保持为物理长度 `38160`；
   - 删除后，崩溃时机略有变化，但 KY 内核在 `ldpc_decoder` 开头仍打印出 warning：

     > info\_len + hm\_m > blk\_len

   - 说明 **根本问题在于 KY 内核自身使用了 `hm_m` 作为外部缓冲访问长度**。

2. **在 C-Model 解码入口使用 hm\_k/hm\_m 重构 parity**
   - C-Model 在重构 `[Info | Pad | Parity]` 时，曾按 `hm_m` 长度从 `det_blk` 拷贝 parity：

     ```c
     for (i = 0; i < hm_m; i++)
         tmp_blk[hm_k + i] = det_blk[info_len + i];
     ```

   - 这同样越界读取 `det_blk` 尾部 $mask\_len$ 位，因而需一并修正。

仅修 C-Model 而不动 KY 内核时，`ldpc_clean` 崩溃仍然存在，说明 KY 编码器/解码器内部的长度使用也必须调整。

### 9.5 正确修正思路：引入 phys\_parity\_len 作为物理 parity 长度

通过对长度关系的梳理，得到关键约束：

- 外部物理链路（包含 `tx_blk/det_blk/dec_blk`，以及 SV ↔ C DPI 接口）只能使用：

  - $blk\_len$ 作为 **bit 流总长度上界**；
  - $info\_len$ 与 $phys\_parity\_len = blk\_len - info\_len$ 作为 **信息段/校验段的物理长度**。

- KY 内部逻辑仍可以使用：

  - $hm_m, hm_n, hm_k, mask\_len$ 处理“部分 circulant 块”的 pad & mask；
  - 但凡涉及 `tx_blk/det_blk/dec_blk` 的访问，都必须以 $blk\_len$ 或 $phys\_parity\_len$ 为边界。

据此，修正应遵循以下原则：

1. 在 `ldpc_c_model_gen4.c` 中，引入：

   ```c
   int info_len        = sim_pckt->info_len;
   int hm_k            = sim_pckt->hm_k;          // 逻辑信息位长度
   int blk_len         = sim_pckt->blk_len;       // 物理码长 = H_N
   int phys_parity_len = blk_len - info_len;      // 物理 parity 长度
   int hm_m            = sim_pckt->hm_m;          // 逻辑 parity 长度
   int mask_len        = hm_m - phys_parity_len;  // = mask_len
   ```

2. 在 KY 编/解码器中，同样计算并使用 `phys_parity_len`，所有与 `tx_blk/det_blk/dec_blk` 相关的拷贝不得再使用 `hm_m` 作为长度。

### 9.6 具体修正步骤（C-Model + KY 内核）

#### 9.6.1 C-Model 解码入口 ldpc\_dec 的修正

1. **不再修改 `sim_pckt->blk_len`**
   - 删除 `ldpc_dec` 中任何赋值 `sim_pckt->blk_len = hm_k + hm_m` 的语句，确保 `blk_len` 始终为物理长度。

2. **重构 `[Info | Pad | Parity]` 时使用 phys\_parity\_len**

   重构逻辑示意：

   ```c
   int info_len        = sim_pckt->info_len;
   int hm_k            = sim_pckt->hm_k;
   int blk_len         = sim_pckt->blk_len;
   int phys_parity_len = blk_len - info_len;
   int hm_m            = sim_pckt->hm_m;

   // 先拷贝 Info
   for (int i = 0; i < info_len; i++)
       tmp_blk[i] = sim_pckt->det_blk[i];

   // 填充 Pad 区（逻辑上是 hm_k - info_len 位）
   for (int i = info_len; i < hm_k; i++)
       tmp_blk[i] = 0;

   // 拷贝物理存在的 parity
   for (int i = 0; i < phys_parity_len; i++)
       tmp_blk[hm_k + i] = sim_pckt->det_blk[info_len + i];

   // 剩余 mask 区用 0 或最大 LLR 填充（视软判/硬判实现）
   for (int i = phys_parity_len; i < hm_m; i++)
       tmp_blk[hm_k + i] = 0;
   ```

   这里的关键是：**从 `det_blk` 读出的 parity 长度只能是 `phys_parity_len`，绝不能用 `hm_m`。**

#### 9.6.2 KY 解码器 ldpc\_decoder 的修正

1. 在解码入口计算：

   ```c++
   int info_len        = this->info_len;
   int blk_len         = this->blk_len;
   int phys_parity_len = blk_len - info_len;
   int hm_k            = this->hm_k;
   int hm_m            = this->hm_m;
   ```

2. **det\_blk → dec\_di\_blk 时只读 phys\_parity\_len**

   ```c++
   // Info 段
   vec_copy(det_blk, dec_di_blk, 0,      0,      info_len);

   // Parity 段：只拷贝物理存在的部分
   vec_copy(det_blk,     dec_di_blk,
            info_len,    hm_k,
            phys_parity_len);

   // 剩余 mask 区：用最大 LLR 或 0 填充
   for (int i = phys_parity_len; i < hm_m; ++i)
       dec_di_blk[hm_k + i] = max_llr_bin;   // 软判；硬判场景可用 0
   ```

3. **dec\_do\_blk → dec\_blk 时同样只写 phys\_parity\_len**

   ```c++
   // Info 段
   vec_copy(dec_do_blk, dec_blk, 0,      0,      info_len);

   // Parity 段：只写回物理存在的部分
   vec_copy(dec_do_blk, dec_blk,
            hm_k,       info_len,
            phys_parity_len);
   ```

4. 这样就保证：
   - 所有对 `det_blk/dec_blk` 的访问都在 `[0, blk_len)` 范围内；
   - KY 内部仍可使用 `hm_m` 处理 mask，但不再向外部缓冲溢出。

#### 9.6.3 KY 编码器 ldpc\_encoder 的修正

编码器在去 padding 时也需用 `phys_parity_len` 限制对 `tx_blk` 的访问：

```c++
int info_len        = this->info_len;
int blk_len         = this->blk_len;
int phys_parity_len = blk_len - info_len;
int hm_k            = this->hm_k;
int hm_m            = this->hm_m;

// removing 0 padding
vec_copy(tx_blk, enc_do_blk, 0,       0,        info_len);

// 只拷贝物理存在的 parity
vec_copy(tx_blk, enc_do_blk,
         hm_k,   info_len,
         phys_parity_len);

// 对 hm_m - phys_parity_len 区域可按需要清零或忽略
for (int i = phys_parity_len; i < hm_m; ++i)
    enc_do_blk[info_len + i] = 0;
```

在实际排查过程中，验证结果为：

- **仅修解码（ldpc\_dec + ldpc\_decoder）仍会出现 `ldpc_clean` 崩溃**；
- **同时修正编码器（ldpc\_encoder）后，`pad_bit = 16/48` 等所有配置下崩溃完全消失**。

说明编码阶段对 `tx_blk` 的越界读同样会破坏堆结构，只是崩溃延后到 `ldpc_clean` 才暴露。

### 9.7 问题解决流程小结

1. **问题暴露**
   - 在 `pad_bit ≠ 256` 配置下，仿真在 `ldpc_cleanup` 中崩溃；
   - 打印发现 `ldpc_decoder` 内部曾把 `blk_len` 改为 `hm_k + hm_m`，并在 `ldpc_clean` 中 `mod2sparse_free(qc_e)` 崩溃。

2. **初步分析**
   - 打印 `bm_m, bm_n, hm_m, hm_n, hm_k, pad_len, pad_bit, mask_len, blk_len`，发现：
     - `hm_m = 5376`，`blk_len - info_len = 5136`，差值恰为 `mask_len = 240`；
   - 判定为“逻辑 parity 长度 `hm_m` 与物理 parity 长度 `phys_parity_len` 不一致导致的越界访问”。

3. **排除 C-Model 内部 blk\_len 修改**
   - 删除 `ldpc_dec` 中修改 `blk_len` 的代码，使 `blk_len` 在 C-Model 内部保持为 `38160`；
   - 崩溃仍存在，且 KY 解码器打印了 `info_len + hm_m > blk_len` 的 warning。

4. **引入 phys\_parity\_len 概念**
   - 在 C-Model 与 KY 内核中显式定义：
     - $phys\_parity\_len = blk\_len - info\_len$，
   - 用来约束所有对外部缓冲（`tx_blk/det_blk/dec_blk`）的访问。

5. **解码路径修正**
   - 在 `ldpc_c_model_gen4.c::ldpc_dec` 中：
     - 重构 `[Info | Pad | Parity]` 时，以 `phys_parity_len` 为界读取 `det_blk`；
   - 在 `ldpc_codec_ky.cpp::ldpc_decoder` 中：
     - `det_blk → dec_di_blk`、`dec_do_blk → dec_blk` 两处 vec\_copy 的长度从 `hm_m` 改为 `phys_parity_len`，并对剩余 `hm_m - phys_parity_len` 区域填充默认值。

6. **编码路径修正**
   - 在 `ldpc_codec_ky.cpp::ldpc_encoder` 中：
     - `tx_blk → enc_do_blk` 去 padding 时的 parity 部分长度从 `hm_m` 改为 `phys_parity_len`；
   - 修正后再测试，`pad_bit = 16/48` 等配置均不再触发 `ldpc_clean` 崩溃。

7. **最终结论**
   - 问题根因：**Gen4 KY 内核在 pad\_bit ≠ 256 时仍按逻辑 parity 长度 `hm_m` 访问外部物理缓冲 `tx_blk/det_blk/dec_blk`，导致越界写/读，堆结构被破坏，在 `ldpc_clean` 中释放 `mod2sparse` 结构时崩溃。**
   - 解决思路：在 C-Model 与 KY 内核中统一引入物理 parity 长度 `phys_parity_len = blk_len - info_len`，所有面向外部 bit 流的访问均以此为边界，内部“部分 circulant”处理仍由 `hm_m/hm_n/mask_len` 完成。


---

## 9. Gen4 KY/DQ 长度失配导致 `ldpc_clean` 崩溃的 8D 分析

本节记录在 Gen4 KY/DQ 适配过程中，因 **逻辑码长与物理码长失配** 导致 `ldpc_clean` 阶段崩溃（`corrupted size vs. prev_size`）的完整 8D 分析。  
该问题发生在使用 KY/DQ 内核（`gen4_ldpc_sim/src/ldpc_codec_ky.cpp`）与 Gen4 C-Model (`ldpc_c_model_gen4.c`) 组合时，特别是在 `pad_bit != 256`（例如 16/48）配置下。

### D1：问题描述

- 现象：
  - 在 `pad_bit = 16` / `48` 等非整 circulant 配置下，仿真结束调用 `ldpc_cleanup()` 时，进到 `sim_pckt->ldpc_clean()` 立刻报错：
    - glibc 抛出 `corrupted size vs. prev_size`，进程以 `SIGABRT` 终止；
  - 在 `ldpc_clean` 中逐段加打印后发现：
    - 打印到 `start mod2sparse_free(qc_e)` 后立即 abort；
    - `qc_bm/qc_hm/qc_a/...` 等前面的 `mod2sparse_free` 可以正常执行。
- 约束条件：
  - Gen4 C-Model 使用的是 KY/DQ 内核（`ldpc_codec_ky.cpp`），并非 Gen3 的 `ldpc_codec.cpp`；
  - C-Model 前端 `ldpc_c_model_gen4.c` 中引入了 `pad_bit`，试图模仿 Gen4 DQ 的行为。

### D2：组建问题解决小组

- 代码层参与组件：
  - C-Model 外壳：`ldpc_c_model_gen4.c`
  - Gen4 KY/DQ 内核：`gen4_ldpc_sim/src/ldpc_codec_ky.cpp`
  - 通道与打包：`gen4_ldpc_sim/src/transceiver.*`、`vec_op.c`
- 分工：
  - C-Model 层负责确认 `info_len` / `blk_len` 配置及 DPI 打包逻辑；
  - KY 内核负责 H 矩阵维度、编解码逻辑及资源释放；
  - 通道层负责 `tx_blk` / `det_blk` 的物理长度。

### D3：问题遏制措施

- 为避免进一步误判，首先：
  - 在 `ldpc_cleanup()` 增加分段打印，确认崩溃点位于 `sim_pckt->ldpc_clean()` 内；
  - 在 `ldpc_packet::ldpc_clean()` 中增加入口打印，记录 `bm_m, bm_n, hm_m, hm_n, hm_k, pad_len, pad_bit_num, mask_len, blk_len` 以及指针地址；
  - 在 `ldpc_decoder()` 开头增加打印，记录同样的长度信息。
- 通过这些打印确认：
  - `ldpc_clean` 入口时：
    - `bm_m=21, bm_n=150, hm_m=5376, hm_n=38400, hm_k=33024, pad_len=0, pad_bit=16, mask_len=240, blk_len=38160`；
  - `ldpc_decoder` 入口时（C-Model 曾经修改 blk_len 时）：
    - `blk_len=38400`；
  - 最终 fix 后，`ldpc_decoder` 入口也看到 `blk_len=38160`。

### D4：根本原因分析

**4.1 两套码长定义：逻辑 vs 物理**

- KY 内核 (`ldpc_codec_ky.cpp::ldpc_config_dq`) 中：
  - `bm_m = h_m + 1, bm_n = h_n + 1, cir_sz = sc`
  - 逻辑 H 矩阵尺寸：
    - `hm_m = bm_m * sc`，例如 `21 * 256 = 5376`
    - `hm_n = bm_n * sc`，例如 `150 * 256 = 38400`
    - `hm_k = hm_n - hm_m = 33024`
  - DQ 特有参数：
    - `pad_bit_num = pad_bit`
    - `drop_len = pad_bit`
    - `mask_len = sc - drop_len`，例如 `256 - 16 = 240`
  - 因此：
    - 逻辑 parity 长：`hm_m = H_M + mask_len`；
    - 逻辑码长：`hm_n = H_N + mask_len`；
    - 其中 `H_M = h_m * sc + pad_bit`，`H_N = h_n * sc + pad_bit`。

- 通道与 C-Model 视角：
  - C-Model 中根据 Gen4 DQ 设计计算：
    - `H_M = h_m * sc + pad_bit`，`H_N = h_n * sc + pad_bit`；
    - `blk_len = H_N`（信道存储码字长度）。
  - 通道分配（`transceiver.cpp::ch_pckt_alloc`）：
    - `tx_blk` / `det_blk` 长度均为 `blk_len = H_N = 38160`。

因此存在 **逻辑维度与物理维度的天然差异**：

- 内核认为码长是 `hm_n = 38400`；
- 信道与 C-Model 认为码长是 `blk_len = H_N = 38160`；
- 差异量恰好是 `mask_len = sc - pad_bit = 240`。

**4.2 C-Model 对 blk_len 的临时改写**

- 初始版本在 `ldpc_c_model_gen4.c::ldpc_dec` 中，为了让 KY 内核看到完整 `[Info|Pad|Parity]`，曾加入：

  ```c
  char *old_det    = sim_pckt->det_blk;
  int   old_blk_len = sim_pckt->blk_len;

  sim_pckt->det_blk = tmp_blk;
  sim_pckt->blk_len = hm_k + hm_m;   // 38400

  sim_pckt->ldpc_decoder(...);

  sim_pckt->det_blk = old_det;
  sim_pckt->blk_len = old_blk_len;   // 38160
  ```

- 结果：
  - `ldpc_decoder` 内看到的 `blk_len = 38400`；
  - 但 `tx_blk` / `det_blk` / `dec_blk` 在 `ldpc_pckt_alloc` 中是按 `blk_len = 38160` 分配的；
  - 即：**数组真实容量是 38160，但逻辑上被当作 38400 使用**。

**4.3 KY 解码器中的越界访问**

在 `ldpc_codec_ky.cpp::ldpc_decoder` 中（原始实现）：

```c++
// add 0 padding
vec_copy(det_blk, dec_di_blk, 0,      0,        info_len);
...
vec_copy(det_blk, dec_di_blk, info_len, hm_k,   hm_m);
...
//remove padding
vec_copy(dec_do_blk, dec_blk, 0,      0,        info_len);
vec_copy(dec_do_blk, dec_blk, hm_k,   info_len, hm_m);
```

- 对 `det_blk` 的访问：
  - 第二行读取 `det_blk[info_len .. info_len + hm_m - 1]`；
  - 具体为 `det_blk[33024 .. 38399]`；
  - 但 `det_blk` 长度是 `blk_len = H_N = 38160`，最大合法索引是 `38159`；
  - 因此读越界 `38400 - 38160 = 240 = mask_len` 个 bit。

- 对 `dec_blk` 的访问：
  - 最后一行写 `dec_blk[info_len .. info_len + hm_m - 1]`；
  - 目标索引同样是 `33024 .. 38399`，而 `dec_blk` 当初按 `blk_len=38160` 分配；
  - 导致写越界 240 个 bit，破坏堆。

这些越界不会立刻报错，而是在后续 `ldpc_clean` 中 `mod2sparse_free(qc_e)` 首次触碰被破坏 chunk 时，glibc 检测到 heap 元数据不一致，引发 `corrupted size vs. prev_size`。

**4.4 C-Model 重构阶段对 det_blk 的越界访问**

- C-Model 的 `ldpc_dec` 在调用内核前，会重构 `[Info | Pad | Parity]`：

  ```c
  for (i = 0; i < info_len; i++)
      tmp_blk[i] = det_blk[i];
  for (i = 0; i < pad_len; i++)
      tmp_blk[info_len + i] = 0;
  for (i = 0; i < hm_m; i++)
      tmp_blk[hm_k + i] = det_blk[info_len + i];
  ```

- 这里第三个 `for` 同样用 `hm_m` 作为长度，从 `det_blk[info_len..info_len+hm_m-1]` 读取 parity；
- 因为 `info_len + hm_m = hm_n = 38400 > blk_len = 38160`，也会发生读越界 240 bit；
- 即便后来不再修改 `blk_len`，这一写仍旧早早破坏 `det_blk` 附近的 heap metadata。

**4.5 KY 编码器中的越界访问**

- `ldpc_encoder()` 末尾的 “removing 0 padding” 部分：

  ```c++
  vec_copy(tx_blk, enc_do_blk, 0, 0, info_len);
  vec_copy(tx_blk, enc_do_blk, hm_k, info_len, hm_m);
  ```

- 按 `vec_copy(vs,vd,vs_si,vd_si,d)` 的定义：
  - 源 `vs = tx_blk`，长度 `blk_len`；
  - 第二行读 `tx_blk[hm_k .. hm_k+hm_m-1] = tx_blk[33024..38399]`，仍然超出 `38159`；
  - 这里是读越界（不会直接损坏堆，但在工具下会报 invalid read）。

综上，**根本原因是：KY 内核在逻辑维度上使用 `hm_m/hm_n`，而 C-Model 与通道只分配了 `blk_len = H_N` 的物理 buffer，且在编解码入口/出口对 `hm_m` 与 `blk_len` 的差异未做映射/裁剪，导致 det_blk/dec_blk/tx_blk 被按 `hm_n` 维度访问，越界长度正好是 `mask_len`。**

### D5：永久性纠正措施

整体修复策略：

1. 在 C-Model 层不再改写 `blk_len`，让 KY 内核看到的 `blk_len` 始终是物理码长 `H_N`；
2. 在 C-Model 重构阶段，引入“物理 parity 长度” `phys_parity_len = blk_len - info_len = H_M`，只从 `det_blk` 读取 `phys_parity_len` 个 parity；对于逻辑上多出来的 `hm_m - phys_parity_len = mask_len` 部分，在临时 buffer 中补 0；
3. 在 KY 解码入口（`ldpc_decoder`）中，同样以 `phys_parity_len` 为界，从 `det_blk` 读 parity、向 `dec_blk` 写 parity，剩余部分在内部 `dec_di_blk/dec_do_blk` 中用 `max_llr_bin` 表示“虚拟强 0”；
4. 在 KY 编码出口（`ldpc_encoder`）中，以 `phys_parity_len` 为界，从内部 `enc_do_blk` 只向 `tx_blk` 写 `info_len + phys_parity_len` 个 bit，其余填 0 或忽略。

**具体措施（示意）：**

1. C-Model `ldpc_dec` 重构：

   ```c
   int phys_parity_len = sim_pckt->blk_len - info_len;          // H_M
   if (phys_parity_len > hm_m) phys_parity_len = hm_m;

   // Info
   for (i = 0; i < info_len; i++)
       tmp_blk[i] = det_blk[i];
   // Pad
   for (i = 0; i < pad_len; i++)
       tmp_blk[info_len + i] = 0;
   // Parity 实际收到部分
   for (i = 0; i < phys_parity_len; i++)
       tmp_blk[hm_k + i] = det_blk[info_len + i];
   // Parity 剩余 mask_len 部分补 0
   for (i = phys_parity_len; i < hm_m; i++)
       tmp_blk[hm_k + i] = 0;
   ```

2. KY 解码入口：

   ```c++
   int phys_parity_len = blk_len - info_len;  // H_M
   if (phys_parity_len > hm_m) phys_parity_len = hm_m;

   vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
   for (int i = 0; i < pad_len; i++)
       dec_di_blk[info_len + i] = max_llr_bin;

   // Parity 的物理部分
   vec_copy(det_blk, dec_di_blk,
            info_len, hm_k,
            phys_parity_len);
   // Parity 剩余 mask_len 填 max_llr_bin
   for (int i = phys_parity_len; i < hm_m; i++)
       dec_di_blk[hm_k + i] = max_llr_bin;

   for (int i = 0; i < mask_len; i++)
       dec_di_blk[hm_n + i] = max_llr_bin;

   // ... 解码过程 ...

   vec_copy(dec_do_blk, dec_blk, 0, 0, info_len);
   vec_copy(dec_do_blk, dec_blk,
            hm_k, info_len,
            phys_parity_len);
   ```

3. KY 编码出口：

   ```c++
   int phys_parity_len = blk_len - info_len;  // H_M
   if (phys_parity_len > hm_m) phys_parity_len = hm_m;

   // 从 enc_do_blk → tx_blk
   vec_copy(enc_do_blk, tx_blk, 0,      0,        info_len);
   vec_copy(enc_do_blk, tx_blk,
            hm_k, info_len,
            phys_parity_len);

   // tx_blk 剩余尾部清零（如有）
   for (int i = info_len + phys_parity_len; i < blk_len; i++)
       tx_blk[i] = 0;
   ```

通过上述调整，所有涉及 `det_blk` / `dec_blk` / `tx_blk` 的访问都被限制在 `0..blk_len-1` 范围内，越界访问被彻底消除。

### D6：验证与效果

- 修改后重新运行：
  - `ldpc_decoder` 起始打印中 `info_len + hm_m > blk_len` 的 warning 消失（转为检查 `info_len + phys_parity_len == blk_len`）；
  - C-Model 不再修改 `blk_len`，KY 内核始终看到 `blk_len = H_N`；
  - `ldpc_clean` 中调用 `mod2sparse_free(qc_e)` 不再触发 heap corruption，仿真正常结束；
  - pad_bit 为 256、16、48 等配置均通过基本功能和 DV 烟测。

- 进一步结合 2Byte-align LSB/MSB 修正（见本报告前文），编码/解码/错误注入路径在 SV 日志下的最后 16bit 也表现一致。

### D7：防止再次发生的措施

1. 在 C-Model 文档中明确区分：
   - 内核逻辑维度：`hm_m/hm_n`；
   - 物理传输维度：`H_M/H_N`；
   - 并强调：**任何对 `tx_blk/det_blk/dec_blk` 的访问，都必须以 `blk_len=H_N` 为上界，仅在内部 buffer 中使用 `hm_n`。**
2. 在 KY 内核关键入口处保留边界检查打印（编译开关控制）：
   - `KY_DEC` 打印 `info_len + phys_parity_len` 与 `blk_len` 的关系；
   - 一旦发现 `info_len + phys_parity_len != blk_len`，直接报错退出。
3. 在未来适配新的矩阵/码长时，优先遵循 “**blk_len = H_N**，内部使用 `hm_n`” 的模式，避免再次出现“先按 A 分配，后按 B 访问”的情况。

### D8：经验教训

- 对于带有“逻辑扩展维度”（如 KY/DQ 额外一行/一列）的 LDPC 方案，**必须从一开始就把逻辑码长与物理码长分离清楚**，并在接口层做显式映射，而不是混用一个长度变量。
- C-Model 尝试“补救”逻辑维度时（例如临时修改 `blk_len`），如果不同步修改所有 buffer 的实际分配，反而更容易掩盖真正的越界点，导致问题只在 cleanup 时才暴露。
- 对于此类问题，系统性地画出“长度链路”（从顶层配置参数 → C-Model → 通道 → 内核 → cleanup），再结合局部的 `index vs capacity` 打印，是比盲目 valgrind 更高效的定位方式。
