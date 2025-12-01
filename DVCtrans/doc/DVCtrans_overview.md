# DVCtrans 框架与 `ldpc_c_model.c` 说明

## 1. 框架定位与组成
DVCtrans 作为部门 DV 使用的 DPI C-model，提供可直接对接 SystemVerilog 的 LDPC 编码、解码与信道注入能力。相较 `gen3_ldpc_sim` 的自驱动仿真平台，DVCtrans 保留数学核心与信道建模，去除数据格式/扰码/MCRC/统计输出，新增 DPI 友好的接口层。

主要目录与组件：
- `src/ldpc_c_model.c`：DPI 封装层，提供 `extern "C"` 接口给 SV 调用。
- `src/transceiver.{h,cpp}`：信道与软判决生成，支持 CLEAN/AWGN/BSC/ERR_INJ/MAX_ERR/ALL_ZERO，含软读 bin 划分与 LLR 表生成。
- `src/ldpc_codec.{h,cpp}`：QC-LDPC 读入、编码、BF/Layer/TBF 解码，寄存器化的解码配置（含 SDLite LLR 配置）。
- 公共库：`mod2*`（稀疏/稠密矩阵与变换）、`vec_op`、`rand`、`finite_lib`、`alloc`、`intio`，承担 GF(2) 运算、向量操作、随机与定点量化。

整体调用链：SV 侧通过 DPI 入口 → `ldpc_config` 完成信道与 LDPC 初始化 → `ldpc_enc`/`ldpc_dec` 完成编解码 → `ch_update`/`ch_err_inj`/`sd_err_inj` 支持不同的信道场景与错误注入 → `ldpc_cleanup` 释放资源。

## 2. 数据结构与差异要点
- `ldpc_packet` 继承 `ch_packet`，内部持有 QC 矩阵、G 矩阵、解码阈值、量化上下界，以及数据块缓冲（`usr_blk/tx_blk/det_blk/dec_blk` 等）。  
- `ch_packet` 在 DVCtrans 中新增 `rd_blk/split_bin/bin_id/bin_distr` 等软判决缓存，`ch_llr_alloc/ch_llr_gen` 固定 128 大小的表，用于可重复的软读建模。  
- 解码模式集：`SKIP/BF_P0/BF_P3/LAYER/TBFDEC`，其中 `TBFDEC` 是面向软/三值 BF 的变体。  
- 量化参数：`ldpc_dec_config` 接收 `finite_q_num/finite_r_num/finite_f_num` 及 `reg_sdlite_llr*`，用于对齐硬件寄存器配置；饱和量化通过 $Sat\_Quan(x)$，总比特宽 `tot_num`、小数位 `frac_num`，饱和值 $\pm (2^{tot\_num-1}-1)/2^{frac\_num}$。

## 3. `src/ldpc_c_model.c` 接口详解
### 3.1 典型调用顺序
1. `ldpc_config(...)`：一次性配置信道与 LDPC，并分配缓冲与 LLR 表。  
2. 编码：`ldpc_enc(usr_data_sv, enc_data_sv, debug)`。  
3. 信道与错误注入（可选）：`ch_update` 动态更新信道；`ch_err_inj` 硬判注入；`sd_err_inj` 软判注入并返回 bin/LLR 信息。  
4. 解码：`ldpc_dec(det_data_sv, dec_data_sv, ... , dec_mode, sd_num, debug, h_n, h_sc)`。  
5. 资源释放：`ldpc_cleanup()`。

### 3.2 `ldpc_config`（初始化）
函数签名（参数顺序与缩放）：
- `h_m, h_n, h_sc, h_st, h_wt`：QC 基矩阵行、列、扩展因子、T 子矩阵行数、列重。块信息长度为 `info_num`，码字长度为 `h_n*h_sc - pad_num`，其中 $pad\_num = (h_n - h_m)h\_sc - info\_num$。  
- 信道：`ch_mode`（0..4 映射 CLEAN/AWGN/BSC/ERR_INJ/MAX_ERR），`ch_para` 以百分比输入，内部换算 $m\_ch\_para = ch\_para/100$（AWGN 为 dB，BSC/ERR 为概率/错误数）。  
- 解码参数：`fdec_max_itr, ldec_max_itr, alpha`（百分比输入，内部 `m_alpha = alpha/100`），`finite_q_num/finite_r_num` 定点宽度，`sd_num` 软读次数。`llr0/llr1` 以 $1/16$ 为步长输入。  
- SDLite 配置：`sdlite_llr_config`、`sdlite_llr0..3` 直接写入寄存器字段。  
- `v_ref_sv`：当 `sd_num>=2` 时读取 SV 数组，按 $1/1000$ 转为电压阈值；否则默认硬判。  
- 内部步骤：构造校验矩阵文件名 `ldpc_h_m_n_sc_wt_wt.txt` → `ch_config` → `ldpc_config` → `ldpc_pckt_alloc` → `ch_llr_alloc` → `ch_llr_gen(sd_num, vref, m_llr0, m_llr1, finite_q_num-1, 4)` → `ldpc_dec_config(...)`。

### 3.3 `ldpc_cleanup`
释放 LDPC 结构、数据块与 LLR/信道缓冲，对应 `ldpc_clean`、`ldpc_pckt_clean`、`ch_llr_clean`。

### 3.4 `ldpc_enc`
- 输入 `usr_data_sv`（SV 数组，32b 打包，MSB 先行），解包到 `usr_blk`。  
- 调用 `ldpc_encoder` 生成码字；输出 `enc_data_sv`（同样 32b 打包）。填充顺序：信息比特在前，校验比特追加。

### 3.5 `ldpc_dec`
- 输入判决：若 `sd_num<2`，按 32b 硬判打包解出；否则直接读取软判数组。  
- `dec_mode`：0→`BF_P3`，1→`LAYER`，2→`TBFDEC`。  
- 解码后返回：`dec_unc_sv`（是否失败）、`init_synd_wt_sv`、`dec_cnvg_itr_sv`（收敛迭代）、`dec_cnvg_col_sv`（收敛列/层）、`fina_synd_wt_sv`，并输出解码结果 `dec_data_sv`（32b 打包）。多余位填 0 直至 $h_n \times h\_sc$ 对齐。

### 3.6 `ch_update`
在不重新分配缓冲的情况下，更新信道模式/参数并重算噪声方差；参数缩放同 `ldpc_config`。

### 3.7 `ch_err_inj`
- 输入 `tx_data_sv` 硬判比特（32b 打包），写入 `tx_blk`。  
- 调用 `ch_transmit` 执行信道/错误注入，再用 `ch_detector(1, &vref)` 硬判。  
- 输出 `rx_data_sv` 为判决比特（32b 打包）。适合硬判 DV 场景。

### 3.8 `sd_err_inj`
- 用软读数 `sd_num>=2` 与 `v_ref_sv`（$1/1000$ 缩放）生成软判。  
- `ch_detector(sd_num, vref)` 返回：  
  - `rx_data_sv`（bin 序列，长度等于 CW 比特），  
  - `rd_data_sv`（各次读的原始 0/1，按 32b 打包，长度 `blk_len*sd_num`），  
  - `llr_tbl_sv`（128 元素，LLR 放大 16 倍存储），  
  - `split_bin_sv`（被分裂的 bin ID），  
  - `bin_id_sv`（按 Vref 升序的 bin 编号）。  
- 便于 DV 侧重现软读分布、LLR 量化与 bin 划分。

## 4. 典型集成建议
- 初始化与资源管理需成对调用：`ldpc_config` 后在仿真结束调用 `ldpc_cleanup`。  
- 软判量化：`tot_num = finite_q_num-1`，`frac_num=4`；若 SV 侧需要验证量化，可对照输出的 `llr_tbl_sv` 与输入 Vref。  
- 校验矩阵文件需放置在仿真工作目录或通过路径配置，命名需匹配 `ldpc_h_<m>_<n>_<sc>_<wt>_<wt>.txt`。  
- 若只做硬判，`sd_num` 设 1，可节省 Vref 传递与软读分布计算。  
- 解码模式选择：BF_P3（硬判 BF）、LAYER（层迭代最小和）、TBFDEC（软/三值 BF），可根据 RTL 顶层模式寄存器对应选择。


## gen4 to DVC
Gen4 DQ (Partial Circulant) 模式的核心在于：物理矩阵扩展了一行一列，但最后一个 Circulant 是部分有效的（由 pad_bit
  决定有效长度，mask_len 决定无效长度）。

  ---

  1. 参数计算与设置 (ldpc_config)

  在 ldpc_c_model_gen4.c 的 ldpc_config 函数中，你需要计算以下全局参数 (g_ 开头变量)。

  输入参数：
   * h_m, h_n: 逻辑基矩阵行列数 (e.g., 20, 149)
   * h_sc: 循环块大小 (Z, e.g., 256)
   * info_num: 用户实际数据长度 (e.g., 4KB + Meta)
   * pad_bit: 最后一个 Partial Block 的有效长度 (e.g., 128)

  计算公式：

```
  // 1. 物理基矩阵维度 (Gen4 扩展)
  int bm_m = h_m + 1;
  int bm_n = h_n + 1;
  
  // 2. Mask 参数
  int mask_len = h_sc - pad_bit; // 无效尾部长度
  
  // 3. 系统信息位总容量 (Systematic Bits Capacity)
  // 注意: hm_k 不受 mask 影响，它总是完整的 (N-M)*Z
  // Gen4 结构保证了前 bm_k 列是完整的 Circulant
  g_hm_k = (bm_n - bm_m) * h_sc; 
  
  // 4. 物理校验位长度 (Parity Length)
  // 校验位对应矩阵的行。最后一行被 Mask 截断，所以总长减去 mask_len
  g_hm_m = bm_m * h_sc - mask_len;
  
  // 5. 物理码字总长 (Codeword Length)
  // 码字 = 系统位 + 校验位
  g_hm_n = g_hm_k + g_hm_m;
  
  // 6. 补零长度 (Padding Zeros)
  // 为了填满 hm_k 系统位容量，需要在用户数据后补的 0
  g_pad_num = g_hm_k - info_num;
  
  // 7. 传输块长度 (Interface Block Length)
  // 实际在 DPI 接口上传输的有效数据 = 用户数据 + 校验位 (不传 Padding)
  g_blk_len = info_num + g_hm_m; 
  // 或者等价于: g_hm_n - g_pad_num
```
  ---

  2. 比特数据流详解

  以下展示数据在 编码 -> 传输 -> 解码 全过程中的形态变化。

  阶段 A: 编码 (Encoding)

   1. 输入 (User Data):
       * 来源: ldpc_enc 的 usr_data_sv
       * 变量: sim_pckt->usr_blk
       * 内容: [ User_Data ]
       * 长度: info_num

   2. 补零 (Padding / Shortening):
       * 操作: 在 usr_blk 尾部填充 0
       * 变量: sim_pckt->enc_di_blk (编码器输入)
       * 内容: [ User_Data (info_num) | 00...00 (pad_num) ]
       * 总长: g_hm_k

   3. 核心编码 (LDPC Encoder):
       * 操作: 计算校验位
       * 变量: sim_pckt->enc_do_blk
       * 内容: [ User_Data | 00...00 | Parity (hm_m) ]
       * 总长: g_hm_n

   4. 打包输出 (Packing):
       * 操作: 去除 Padding，仅输出有效数据
       * 目标: enc_data_sv
       * 内容: [ User_Data (info_num) | Parity (hm_m) ]
       * 总长: g_blk_len

  阶段 B: 传输 (Channel)

   * 接口: ch_err_inj 或 sd_err_inj
   * 数据: 保持 [ User | Parity ] 结构，长度 g_blk_len。
   * 注入: 错误注入在此长度范围内进行。

  阶段 C: 解码 (Decoding)

  这是最容易出错的地方，必须显式重构数据结构。

   1. 输入 (Detection Data):
       * 来源: ldpc_dec 的 det_data_sv
       * 内容: [ Rx_User (info_num) | Rx_Parity (hm_m) ]
       * 长度: g_blk_len

   2. 重构 (Reconstruction) - 关键步骤:
       * 操作: 将 Rx_Parity 向后搬移，中间插入强判决 0 (或最大LLR)
       * 变量: sim_pckt->det_blk (或拷贝到 dec_di_blk)
       * 内存布局变化:

```
          输入 Buffer: [ Rx_User ... | Rx_Parity ... ]
                                     ^ 
                                     |
          (搬移 Parity) --------------+
          |
          v
          目标 Buffer: [ Rx_User ... | 00...00 | Rx_Parity ... ]
                       <--info_num--> <-pad_num-> <---hm_m----->
```
       * 总长: g_hm_n
   3. 核心译码 (LDPC Decoder):
       * 输入: sim_pckt->dec_di_blk (长度 g_hm_n)
       * 算法: 使用 (bm_m, bm_n) 矩阵进行迭代，期间利用 mask_matrix 屏蔽无效节点。
       * 输出: sim_pckt->dec_do_blk (长度 g_hm_n)

   4. 提取输出 (Unpacking):
       * 操作: 仅提取前面的用户数据
       * 内容: dec_do_blk[0 ... info_num-1]
       * 目标: dec_data_sv

  ---

  3. 关键代码实现片段 (参考)

  在 ldpc_c_model_gen4.c 中：

  `ldpc_config`:

``` 
  // 必须计算并保存这些全局变量
  g_info_len = info_num;
  g_hm_k     = (h_n + 1 - (h_m + 1)) * h_sc; 
  g_hm_m     = (h_m + 1) * h_sc - (h_sc - pad_bit);
  g_pad_num  = g_hm_k - info_num;
  g_blk_len  = info_num + g_hm_m; 
```
  `ldpc_dec` (重构逻辑):
```
  // 假设 det_blk 已经读入了 g_blk_len 长度的数据 [User | Parity]
  
  // 1. 从后往前搬移 Parity，防止覆盖
  // src_start = info_num
  // dst_start = hm_k (即 info_num + pad_num)
  // len = hm_m
  for (int i = g_hm_m - 1; i >= 0; i--) {
      sim_pckt->det_blk[g_hm_k + i] = sim_pckt->det_blk[g_info_len + i];
  }
  
  // 2. 填充中间的 Padding 区域
  // start = info_num
  // len = pad_num
  for (int i = 0; i < g_pad_num; i++) {
      sim_pckt->det_blk[g_info_len + i] = 0; // 硬判决填0，软判决填最大确信度
  }
```