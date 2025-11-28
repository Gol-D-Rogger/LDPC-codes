# Gen3 LDPC & DVCtrans 比特流与数据流详解

本文档详细描述了 Gen3 LDPC 仿真器 (`gen3_ldpc_sim`) 及其对应的 DPI C-model (`DVCtrans`) 中的比特流结构、长度计算及数据全生命周期流转细节。

## 1. 核心参数与长度公式 (Gen3 Standard)

Gen3 采用标准的 QC-LDPC 结构，没有 Gen4 中的 Padding/Masking 复杂机制。所有维度均基于完整的循环块大小 `h_sc`。

假设配置文件 (`.cnfg` 或 DPI 输入) 参数如下：
*   **`h_m`**: 基矩阵校验行数
*   **`h_n`**: 基矩阵总列数
*   **`h_sc`**: 循环块大小 (扩展因子 Z)
*   **`info_num`** (或 `dsp_info_len`): 用户有效信息长度

### 长度公式表

| 参数名 | 计算公式 | 说明 |
| :--- | :--- | :--- |
| **bm_m** | `h_m` | 基矩阵行数 (无扩展) |
| **bm_n** | `h_n` | 基矩阵列数 (无扩展) |
| **hm_m** | `h_m * h_sc` | **物理校验位总长** (Parity Length) |
| **hm_n** | `h_n * h_sc` | **物理码字总长** (Codeword Length) |
| **hm_k** | `(h_n - h_m) * h_sc` | **系统信息位容量** (Systematic Bits Capacity) |
| **pad_num** | `hm_k - info_num` | **补零数量** (Shortening Zeros) |
| **blk_num** | `hm_n - pad_num`<br>= `info_num + hm_m` | **有效传输块长** (Effective Block Length) |

> **关键点**: 
> *   `hm_n` 是编码器/译码器内部处理的完整长度。
> *   `blk_num` 是信道传输、DPI 接口交互的实际有效长度（不含 Padding）。
> *   `pad_num` 个 0 仅在编解码器内部存在，属于 Shortening 操作。

---

## 2. DVCtrans (DPI C-Model) 数据流详解

DVCtrans 是 Gen3 仿真器的 DPI 封装版本，用于对接 SystemVerilog (SV) 验证环境。

### 2.1 接口数据流向

1.  **配置 (`ldpc_config`)**:
    *   SV 传入 `h_m`, `h_n`, `h_sc`, `info_num` 等参数。
    *   C-Model 内部计算 `hm_k`, `pad_num`, `hm_m`。
    *   分配内存 `usr_blk` (info_num), `enc_di` (hm_k), `enc_do` (hm_n) 等。

2.  **编码 (`ldpc_enc`)**:
    *   **输入**: `usr_data_sv` (32-bit 打包数组, 长度 `info_num`).
    *   **解包**: 将 SV 数据解包至 `usr_blk`。
    *   **补零 (Padding)**: 将 `usr_blk` 拷贝至 `enc_di_blk`，并在尾部填充 `pad_num` 个 0。
        ```
        enc_di_blk = [ User_Data (info_num) | 00...00 (pad_num) ]
        ```
    *   **核心编码**: `ldpc_encoder` 计算校验位 `Parity` (长度 `hm_m`)。
    *   **输出构造**:
        ```
        enc_do_blk = [ User_Data | 00...00 | Parity ] (总长 hm_n)
        ```
    *   **去零与打包**: 去除中间的 `0`，将有效数据打包回 `enc_data_sv`。
        ```
        enc_data_sv = [ User_Data | Parity ] (总长 blk_num)
        ```

3.  **信道/错误注入 (`ch_err_inj` / `sd_err_inj`)**:
    *   **输入**: `tx_data_sv` (即编码输出的 `blk_num` 长度数据)。
    *   **信道处理**:
        *   内部恢复为 `tx_blk`。
        *   添加噪声/翻转比特。
        *   生成 `rx_blk` (软判决 LLR 或 硬判决 bit)。
    *   **输出**: `rx_data_sv` (软判决 bin index 或 硬判决 bit)，长度仍为 `blk_num`。

4.  **译码 (`ldpc_dec`)**:
    *   **输入**: `det_data_sv` (来自信道的 `blk_num` 长度数据)。
    *   **恢复填充 (Re-Padding)**:
        *   将 `det_data_sv` 解包至 `dec_di_blk`。
        *   在信息位 (`info_num`) 之后，校验位之前，**插入** `pad_num` 个强判决 LLR (如最大值 0)。
        ```
        dec_di_blk = [ Rx_User | Max_LLR(0)... | Rx_Parity ] (总长 hm_n)
        ```
    *   **核心译码**: `ldpc_dec_bf` / `ldpc_dec_layer` 对 `hm_n` 长度的数据进行迭代。
    *   **输出处理**:
        *   得到 `dec_do_blk` (长度 `hm_n`, 含 0)。
        *   去除 Padding 0，提取有效用户数据 `info_num`。
    *   **返回**: `dec_data_sv` (32-bit 打包的 `info_num` 长度用户数据)。

### 2.2 DVCtrans 特有数据结构
*   **`ch_packet`**: 包含 `rd_blk` (原始读数据), `split_bin`, `bin_distr` 等，用于支持软判决 (Soft Decision) 的细粒度建模。
*   **`ldpc_packet`**: 继承自 `ch_packet`，增加了 QC 矩阵及编解码 buffer。

---

## 3. Gen3 仿真器 (`gen3_ldpc_sim`) 全流程

与 DVCtrans 相比，Gen3 仿真器包含更多外围模块（加扰、MCRC）。

### 阶段 1: 源数据生成
*   **组成**: `Meta` + `LBA Data`。
*   **长度**: `dsp_src_len = dsp_meta_size + dsp_lba_len`。

### 阶段 2: 预处理 (MCRC & Scramble)
*   **MCRC**: `dsp_src_len` + 32-bit CRC -> `info_len`。
*   **Scramble**: 对 `info_len` 数据进行加扰。
*   **变量**: `usr_blk` (长度 `info_len`)。

### 阶段 3: LDPC 适配 (Shortening)
*   **Padding**: `pad_len = hm_k - info_len`。
*   **Input**: `enc_di_blk = [ usr_blk | 0...0 ]`。

### 阶段 4: 编码与传输
*   **Encoding**: 产生 `hm_m` 长度校验位。
*   **Tx Block**: `tx_blk = [ usr_blk | Parity ]` (去除 0)。
*   **长度**: `dsp_blk_len = info_len + hm_m`。

### 阶段 5: 译码与校验
*   **Rx**: 接收 `dsp_blk_len` 长度数据。
*   **Re-Pad**: 补回 `pad_len` 个 0。
*   **Decoding**: 译码得到 `dec_blk`。
*   **Post-process**: 去 0 -> 解扰 -> MCRC 校验 -> 数据比对。

---

## 4. Gen3 vs Gen4 (DQ) 数据流对比

| 特性 | Gen3 (Standard) | Gen4 (DQ / Partial Circulant) |
| :--- | :--- | :--- |
| **基矩阵维度** | `h_m` x `h_n` | `(h_m+1)` x `(h_n+1)` |
| **校验位长度** | `h_m * h_sc` | `h_m * h_sc + pad_bit` |
| **码字总长** | `h_n * h_sc` | `h_n * h_sc + pad_bit` |
| **Padding 机制** | 纯逻辑补 0 (Shortening) | 逻辑补 0 + **物理结构填充 (Partial Block)** |
| **Masking** | 无 | 需要 `mask_matrix` 和 `vec_mask` 处理边缘块 |
| **DPI 接口长度** | `info + hm_m` | `info + hm_m` (但 `hm_m` 计算含 `pad_bit`) |

```mermaid
graph LR
    subgraph Gen3_Flow
    A[User Info] -->|补0到整块倍数| B(Encoding Input)
    B --> C[Standard QC-LDPC]
    C --> D[Codeword = Info + 00 + Parity]
    D -->|去0| E[Tx Block = Info + Parity]
    end

    subgraph Gen4_DQ_Flow
    A2[User Info] -->|补0到 hm_k| B2(Encoding Input)
    B2 --> C2[Extended QC-LDPC (M+1, N+1)]
    C2 -->|Masking Edge Blocks| D2[Codeword]
    D2 -->|去0| E2[Tx Block = Info + Parity]
    end
    
    style C fill:#e1f5fe
    style C2 fill:#fff9c4
```
