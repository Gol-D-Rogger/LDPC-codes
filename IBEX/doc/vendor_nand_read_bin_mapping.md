# Vendor NAND Read 与 Bin 编号映射说明

本文整理 `src/transceiver.cpp` 中 `VENDOR0` / `VENDOR1` 模式下 NAND read pattern、soft bit 生成和 bin 编号的关系。

这里采用 DV/硬件更常用的视角：**每一刀 Vref 跨所有电压区间形成一条 read 向量**。

## 代码位置

相关逻辑在 `ch_packet::ch_llr_gen()`：

```cpp
hard_bit = nand_read[i][sd_num - 1];

for (int j = 0; j < sd_num - 1; j++) {
    if (nand_read[i][j] == nand_read[i][rd_num - 1 - j]) {
        if (sd_type == VENDOR1)
            soft_bits += 1;
    } else {
        if (sd_type == VENDOR0)
            soft_bits += 1;
    }
}

sd_asc_ord[i] = soft_bits * 2 + hard_bit;
bin_asc_ord[i] = sd_asc_ord[i];
```

当前代码的 bin bit packing 是：

```text
bin = {soft_outer, soft_next, ..., soft_inner, hard}
```

也就是说，**hard bit 放在 bin 编号最低位**。

## 基本约定

对于 `sd_num = n`：

```text
rd_num = 2 * sd_num - 1
```

因此：

```text
Vref 数量 = 2n - 1
电压区间数量 = 2n
```

每一刀 Vref 的 read 规则是：

```text
左边读 1，右边读 0
```

例如一刀位于 A/B 中间，则它跨所有区间的 read vector 是：

```text
1000...
```

如果一刀位于倒数第二个区间边界，例如 E/F 中间，则 read vector 是：

```text
111110...
```

## 代码中的矩阵视角

DV 常看的视角是：

```text
read_index -> 所有电压区间
```

例如：

```text
r0 = 100000
r1 = 110000
r2 = 111000
```

但 C 代码里存的是转置视角：

```cpp
nand_read[region][read_index]
```

也就是：

```text
region -> 所有 read_index
```

这两个视角是同一张 read matrix 的行列互换。

## Vendor1 Soft Bit 规则

`VENDOR1` 使用对称 read pair 做 NXOR：

```text
soft_outer = left_outer NXOR right_outer
soft_next  = left_next  NXOR right_next
...
```

`VENDOR0` 使用同样的对称 pair，但做 XOR。

hard bit 固定取中间 read：

```text
H = r[sd_num - 1]
```

## 1H1S 推导

`sd_num = 2` 时：

```text
rd_num = 3
```

3 个 Vref 切 4 个区间：

```text
      r0    r1    r2
A  |  B  |  C  |  D
```

每一刀跨所有区间的 read vector 是：

```text
r0 = 1000
r1 = 1100
r2 = 1110
```

hard bit：

```text
H = r1 = 1100
```

Vendor1 soft bit：

```text
S = r0 NXOR r2
  = 1000 NXOR 1110
  = 1001
```

当前 bin packing 是 `{S,H}`：

```text
region: A  B  C  D
S:      1  0  0  1
H:      1  1  0  0
bin:    3  1  0  2
```

所以 1H1S 的 Vendor1 bin 顺序是：

```text
3, 1, 0, 2
```

## 1H2S 推导

`sd_num = 3` 时：

```text
rd_num = 5
```

5 个 Vref 切 6 个区间：

```text
      r0    r1    r2    r3    r4
A  |  B  |  C  |  D  |  E  |  F
```

每一刀跨所有区间的 read vector 是：

```text
r0 = 100000
r1 = 110000
r2 = 111000
r3 = 111100
r4 = 111110
```

hard bit：

```text
H = r2 = 111000
```

Vendor1 soft bits：

```text
S1 = r0 NXOR r4 = 100001
S0 = r1 NXOR r3 = 110011
```

当前 bin packing 是 `{S1,S0,H}`：

```text
region: A  B  C  D  E  F
S1:     1  0  0  0  0  1
S0:     1  1  0  0  1  1
H:      1  1  1  0  0  0
bin:    7  3  1  0  2  6
```

所以 1H2S 的 Vendor1 bin 顺序是：

```text
7, 3, 1, 0, 2, 6
```

## 1H3S 推导

`sd_num = 4` 时：

```text
rd_num = 7
```

7 个 Vref 切 8 个区间：

```text
      r0    r1    r2    r3    r4    r5    r6
A  |  B  |  C  |  D  |  E  |  F  |  G  |  H
```

每一刀跨所有区间的 read vector 是：

```text
r0 = 10000000
r1 = 11000000
r2 = 11100000
r3 = 11110000
r4 = 11111000
r5 = 11111100
r6 = 11111110
```

hard bit：

```text
H = r3 = 11110000
```

Vendor1 soft bits：

```text
S2 = r0 NXOR r6 = 10000001
S1 = r1 NXOR r5 = 11000011
S0 = r2 NXOR r4 = 11100111
```

当前 bin packing 是 `{S2,S1,S0,H}`：

```text
region: A   B  C  D  E  F  G   H
S2:     1   0  0  0  0  0  0   1
S1:     1   1  0  0  0  0  1   1
S0:     1   1  1  0  0  1  1   1
H:      1   1  1  1  0  0  0   0
bin:    15  7  3  1  0  2  6   14
```

所以 1H3S 的 Vendor1 bin 顺序是：

```text
15, 7, 3, 1, 0, 2, 6, 14
```

这也解释了 `sd_num = 4` 时，最大 LLR 所在的最右侧 region 从旧 packing 的 `bin7` 变成当前 packing 的 `bin14`。

## sd_num = 7 的结果

`sd_num = 7` 时：

```text
rd_num = 13
```

13 个 Vref 切 14 个区间。当前 packing 为：

```text
bin = {S5,S4,S3,S2,S1,S0,H}
```

Vendor1 从左到右的 bin 顺序是：

```text
127, 63, 31, 15, 7, 3, 1, 0, 2, 6, 14, 30, 62, 126
```

二进制形式：

```text
1111111 = 127
0111111 = 63
0011111 = 31
0001111 = 15
0000111 = 7
0000011 = 3
0000001 = 1
0000000 = 0
0000010 = 2
0000110 = 6
0001110 = 14
0011110 = 30
0111110 = 62
1111110 = 126
```

## 通用公式

对 `sd_num = n`：

```text
rd_num = 2n - 1
region_num = 2n
```

第 `k` 刀 read vector 的形状是：

```text
rk = 111...1100...00
     k+1 个 1，后面全 0
```

hard bit：

```text
H = r[n - 1]
```

Vendor1 soft bits：

```text
S[n-2] = r0       NXOR r[2n-2]
S[n-3] = r1       NXOR r[2n-3]
...
S0     = r[n - 2] NXOR r[n]
```

当前代码打包：

```text
bin = {S[n-2], S[n-3], ..., S0, H}
```

Vendor0 只需要把上面的 NXOR 换成 XOR。

## 和 `nand_read[i][j]` 的对应关系

在 C 代码中：

```cpp
nand_read[i][j] = (vref_asc_ord[i] <= vref_asc_ord[j]) ? 1 : 0;
```

其中：

```text
i = region index
j = read / Vref index
```

因此代码里一行是：

```text
某一个 region 被所有 read 读出来的结果
```

而 DV 中常看的 `r0 = 100000...` 是：

```text
某一个 read 跨所有 region 的结果
```

两者是同一张矩阵的转置。
