# HREModel 仿真链路与 HRE 逻辑

## HRE 概念回顾
- **High-Reliability Error (HRE)** 指本质为错误的比特却被信道以极高置信度报告成某值的情况，对依赖软信息的 LDPC 解码器是最具破坏性的错误类型。
- NAND 退化会产生两类 HRE：
  - *固定 HRE*：由硬件缺陷或复制错误导致，位置稳定。
  - *概率性 HRE*：由阈值漂移、保持力下降等造成，位置随读次变化。
- NAND 原始错位数 $e$ 常建模为二项式分布，尾部概率
  $P(e > T) = 1 - \sum_{e=0}^{T} \binom{n}{e} p^e (1-p)^{n-e}$
  可通过实验测得，以反推 RBER 并驱动仿真。

## 整体仿真流程（`src/SSD_FC.cpp`）
1. `read_arg` / `read_config_file` 读取矩阵、信道、HRE、LLR 量化等配置。`hre_bit`、`hre_mode`、`hre_dec` 在此解析。
2. `ch_config` 初始化信道（SNR/BER/误注入参数）与 HRE 控制字段，并分配 `ch_packet` 缓冲。
3. `ch_llr_alloc → ch_llr_gen` 根据 Vref 或直接模式生成 LLR 表。
4. `dfmt_config`、`ldpc_config`、`ldpc_dec_config`、`rand_config`、`mcrc_config` 初始化数据路径。
5. 每个仿真循环依次执行：数据生成 → 随机化 → `ecc_encoder` → `ch_transmit`（含 HRE 注入）→ `ch_detector` → `ecc_decoder`。

## HRE 配置项（`src/SSD_FC.h`）
- `hre_bit`：每个码字强制注入的 HRE 数，本工程中视为“固定平均数量”。
- `hre_mode`：信道端 HRE 模型，传递给 `hre_model`。
- `hre_dec`：译码端是否对 HRE LLR 做降权处理（`hre_llr`）。

## 信道侧 HRE 注入（`src/transceiver.cpp`）
1. `ch_transmit` 每次调用前执行 `randomBinError(hre_vec, hre_pos, blk_len, hre_num)`：
   - `hre_vec` 设为长度 `blk_len` 的 0/1 向量；
   - `hre_pos` 保存抽中的索引（目前仅用于调试追踪，未被下游使用）。
2. 仅在 AWGN 通道中，若 `hre_vec[i]==1` 则覆盖噪声后的 `rx_blk[i]`：
   - `hre_model==0`：固定写入 `+1.0`，即高置信度“0”；
   - `hre_model==1`：固定写入 `-1.0`，即高置信度“1”；
   - 其他值：随机选择 `±1.0`，实现概率性方向。
3. 其它通道（BSC/ERR_INJ/MAX_ERR）路径目前不检查 `hre_vec`，因此不会注入 HRE。
4. `ch_detector` 直接基于修改后的 `rx_blk` 输出软判决或硬判决，HRE 因此完整地传递到 `dec_di_blk`。

## 译码侧 HRE 处理（`src/ldpc_codec.cpp`）
- 层译码器 `ldpc_dec_layer` 在初始化 `cn_q_mem` 时读取 `llr_tbl[dec_di_blk[*]]` 得到原始 LLR。
- 若配置 `hre_llr==1`（即 `hre_dec` 置位），则对被标记的 HRE 比特将 `cn_q_mem` 强制设为 0，从而在 APP 层面移除其“高置信度”属性。
  - 对信息位与校验位分别检查：索引 `< info_len` 或 `>= hm_k` 的区间各自匹配 `hre_vec`。
  - 译码器其他阶段无额外 HRE 特殊逻辑，因此未启用 `hre_llr` 时，HRE 将以强信赖值扰乱迭代直至失败。

## 现状与局限
- `hre_num` 恒定，尚无“平均10.88个”等混合分布或 RBER 统计回灌逻辑；若需尾部分布，可在配置层外部调度或扩展 `randomBinError`。
- 信道注入仅覆盖 AWGN 路径，若需在 BSC 或错误注入模式下评估 HRE，需要补充对应分支。
- `hre_pos` 目前未回传给译码器或统计模块，可用来做 targeted logging。
- RBER 估计和 HRE 混合建模尚停留在文档层，需要在仿真启动前由用户自行计算并配置 `ch_para`、`hre_bit` 等参数。

### HRE 注入方式（AWGN 分支两种模式）
- 模式 0：随机挑选 `hre_bit` 个位置，翻转无噪声符号且不叠加 AWGN，形成必错且高置信的固定型错误。
- 模式 1：随机挑选 `hre_bit` 个位置，翻转无噪声符号后继续叠加 AWGN，错误方向保持但置信度受噪声波动，极小概率被噪声拉回。


## HRE 注入为何不影响理论 RBER
- **注入原理**：在 AWGN 通道完成调制与噪声叠加后，程序会按照 `hre_num` 数目随机挑出比特，并把这些位置的波形强制设置为固定的 ±1（或随机 ±1）。这一步与噪声过程独立，也不关心该比特在 AWGN 作用下本来是否出错，因此可以视为“在 AWGN 结果上叠加一个定量的、高置信度的干扰”。
- **理论 RBER 定义**：`transceiver.cpp:30-36` 里计算的 `rber = 0.5*(1 + erf(-1/(awgn_sigma*sqrt(2))))` 只跟 AWGN 噪声有关，等价于经典 BPSK/AWGN 的误码概率 $Q(1/\sigma)$。它描述的是“如果没有 HRE，仅凭 AWGN 会导致多少错误”，与后续注入的 HRE 无因果关系。
- **为什么不改 RBER**：理论 RBER 作为基线，用于衡量 AWGN 失真程度；HRE 注入则是额外的、可控的故障模型。即使有 HRE，我们仍希望知道：在同样 SNR 下，纯 AWGN 应该达到的误码率是多少。这样才能通过比较“统计 RAW BER”与“理论 RBER”来量化 HRE 带来的额外劣化。
- **实际误码率的组成**：最终的 Raw BER = AWGN 错误 + HRE 强制产生的错误（部分比特可能刚好与原始符号一致，因此不一定全部变成错误，但平均而言会额外增加 `hre_num/blk_len` 量级的错误概率）。统计结果必然高于理论值，这种差距正是评估 HRE 影响的依据。

## `config/` 目录下 6 份示例配置
文件名遵循 `21x150_hreX_{en,dis}.cnfg`。`X∈{0,1,01}` 对应 `hre_mode` 的设置，`en/dis` 对应 `hre_dec`（译码时是否将 HRE LLR 置零）。六份配置除 HRE 相关字段外，其余参数完全一致：均采用 21×150 QC-LDPC、256 倍扩展、`FC_RDEC` 解码、MANUAL LLR 量化等。以下解释默认 `hre_bit` 为 0；若要真正注入 HRE，需要在相应配置中把 `HRE bit count` 设为所需数量。

- `config/21x150_hre0_dis.cnfg` / `..._en.cnfg`  
  `hre_mode=0` 表示将被标记的 HRE 比特强制为 `rx_blk=+1.0`（译码器看来是“极可信的0”）；`_dis` 关闭译码端 LLR 降权（`hre_dec=0`），`_en` 则开启（`hre_dec=1`）。适合评估“仅有 0 型 HRE”对译码的影响以及译码端补丁的收益。

- `config/21x150_hre1_dis.cnfg` / `..._en.cnfg`  
  `hre_mode=1`，被标记的 HRE 比特强制为 `rx_blk=-1.0`（极可信的1）。同样 `_dis` 表示译码端不降权，`_en` 表示降权。用于模拟 NAND 失效导致“全部朝 1 偏移”的场景。

- `config/21x150_hre01_dis.cnfg` / `..._en.cnfg`  
  `hre_mode=2`，每个 HRE 位独立随机选择 ±1.0，等价于概率性地生成 0 型或 1 型 HRE。`_en/_dis` 的含义与前述一致，便于对比“概率混合 HRE + 译码缓解”组合。

> 注意：上述配置里的 `HRE bit count` 默认为 0，即便启动 `_en` 的译码处理也不会触发任何 HRE 注入。实际仿真时请根据实验计划在配置中写入期望的 HRE 数量。

## LLR 量化与 bin 映射（MANUAL 模式示例）
- 生成阶段（`ch_llr_gen`）：输入参考电压 `vref`（当前示例为 0, 0.15, -0.15, 0.3, -0.3, 0.5, -0.5），先按大小排序为 `vref_asc_ord`，形成 `rd_num+1` 个区间 $(-\infty, v_0), (v_0,v_1), \ldots, (v_{rd\_num-1},+\infty)$，并为每个区间分配 bin ID，当前排序结果为 `[7,5,3,1,2,4,6,0]`。每个 bin 的 LLR 取区间中点的对数似然估计，并按量化位宽饱和。
- 判决阶段（`ch_detector`）：对每个采样 `rx_blk[i]`，从最小 `vref` 起查找第一个满足 `rx < vref_asc_ord[rd_indx]` 的区间，输出对应 bin ID；若全部不满足则落入最右区间（bin 0）。DIRECT 模式则直接用公式 $LLR=\text{SatQuan}(2\,rx/\sigma^2)$。
- 当前 Vref 对应的区间→bin 映射（`sd_type=MANUAL`）：
  - $(-\infty, -0.5)$ → bin 7
  - $[-0.5, -0.3)$ → bin 5
  - $[-0.3, -0.15)$ → bin 3
  - $[-0.15, 0)$ → bin 1
  - $[0, 0.15)$ → bin 2
  - $[0.15, 0.3)$ → bin 4
  - $[0.3, 0.5)$ → bin 6
  - $[0.5, +\infty)$ → bin 0  
  左侧 bin（7/5/3/1）产出负 LLR，对应“高置信度 1”；右侧 bin（2/4/6/0）产出正 LLR，对应“高置信度 0”。更换 Vref 时，bin 排序与阈值会随之调整。

## new logic
新逻辑的目的是让 hre_bit 不再代表“高可靠 LLR 的注入”，而是“确定数量的硬性 bit-flip”，这些 flip 位不再受 AWGN 影响，而是直接把调制输出反转后送往译码器，相当于人为注入 hre_bit 个绝对错误。                     
  - 为了实现这一点，仍然需要生成整帧 AWGN 噪声向量，确保非 HRE 位的行为与原来一致；但对 HRE 位，必须在符号层面把 tx 的 ±1 映射直接翻转（例如原本 -1*(tx*2-1) 给出的是 ±1 的无噪声基值，我们对选中的位直接乘以 -1 或
  强制为固定值），并跳过噪声加成。                                                                                                                                                                                 
  - 这样能让后续统计中，FBC（Raw error count）约等于“AWGN 产生的错 + hre_bit”，因而在做 FER vs FBC/HRE 的表格时，只需从 CW_len × real_rber 中减掉 hre_bit 即可得到 AWGN 部分，再加上确定的 hre_bit 个硬错误。
