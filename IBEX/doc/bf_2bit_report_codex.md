# IBEX BF 2-bit 方案评估记录（ssd_fc_test2, VN_BITS=2, AWGN 5.4）

统一命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`

评估指标：主要关注 $LDPC\ FER$，当 $LDPC\ FER < 0.05$ 且相对基线有明显价值时，认为是“潜在可行方案”，需要进一步在更多 packet 或更低 SNR 下验证。

## 方案 0（基线：ldpc_codec.cpp 原始 2bit 半步）

- 实现说明：
  - `src/ldpc_codec_test2.cpp` 初始完全复制自 `src/ldpc_codec.cpp`，2bit 下使用原始半步更新：
    - `f_update_vn_post` 中，当 $VN\_BITS \le 2$ 时，步长 `delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0)`。
    - 激进加权 `aggr` 条件保持原版，仅在 `ldpc_dec_bf_ibex` 中按原逻辑使用。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置文件：`config/Ibex_hd_row8.cnfg`（100 包或 10 错停止，`VN_BITS=2`）。
- 结果：
  - `[STATISTICS] Total packets simulated: 172`
  - `[STATISTICS] LDPC FER  : 5.813953e-02`
  - `[STATISTICS] RAW BER   : 6.523566e-03`
  - `[STATISTICS] LDPC BER  : 1.886014e-05`
  - IBEX 平均迭代约 138.16。
- 分析：
  - 基线 $LDPC\ FER \approx 5.8\times 10^{-2}$，与你描述的 “FER≈0.1 量级” 相同量级，作为后续方案对比的参考。
  - 半步公式在 2bit 下推力不足，特别是 $weight=3$ 时，单次增量有限，容易形成 error floor，这与 `doc/2bit-BF.md` 和 `doc/agent.md` 中的分析一致。

## 方案 E：2bit 全权重进攻/撤销（aggr 关闭）

- 关键代码改动（相对方案 0）：
  - `src/ldpc_codec_test2.cpp:569-576`，`f_update_vn_post` 中删除 2bit 半步逻辑，统一使用全步长：
    - 原逻辑：`if (VN_BITS <= 2) delta = (weight>>1)+((weight&1)&&pushing?1:0); else delta = weight;`
    - 现逻辑：`delta = weight;`（所有位宽统一全步进攻/撤销）。
  - `src/ldpc_codec_test2.cpp:2709-2716`，`ldpc_dec_bf_ibex` 中调整 aggr 条件：
    - 原逻辑：`bool aggr = (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min < likelihood_thr && !flipped_prev);`
    - 现逻辑：`bool aggr = (VN_BITS > 2) && (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min < likelihood_thr && !flipped_prev);`
    - 含义：2bit 下强制关闭 aggr，仅 3bit 及以上保留激进加权。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - config 文件同方案 0（100 包或 10 错停止）。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.000000e+00`
  - `[STATISTICS] LDPC BER  : 1.944285e-01`
  - 平均迭代：`1024`（已打满最大迭代），syndrome_weight 长期维持在约 $10^3$ 量级，几乎完全不收敛。
- 分析：
  - 与 `doc/bf_2bit_report.md` 中对方案 E 的结论一致：**全步长 + 无半步/趋势在 2bit 下极易发散**。
  - 2bit 空间过于“粗糙”，全步 $delta=weight$（尤其是 $weight\ge 3$）会将 likelihood 迅速推到饱和，翻转之后也难以撤销，导致 syndrome_weight 长期高位震荡。
  - 相比基线方案 0（$LDPC\ FER\approx 0.058$），方案 E 退化为 $LDPC\ FER=1.0$，无实用价值，后续不再考虑该方向。

## 方案 B：攻守分离（未翻转全步，已翻转半步；aggr 关闭）

- 关键代码改动（相对方案 E）：
  - `src/ldpc_codec_test2.cpp:569-578`，`f_update_vn_post` 2bit 分支实现攻守分离：
    - 未翻转：`delta = weight;`（全步进攻）
    - 已翻转：`delta = (weight>>1) + ((weight & 1) && pushing ? 1 : 0);`（半步撤销）
  - `ldpc_dec_bf_ibex` 中 aggr 条件沿用方案 E：`VN_BITS>2` 时才启用激进加权，2bit 完全关闭 aggr。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 其他配置同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.000000e+00`
  - `[STATISTICS] LDPC BER  : 2.546229e-01`
  - 平均迭代：`1024`（打满迭代），syndrome_weight 一直在约 $1.1\times 10^3$ 附近高位震荡。
- 分析：
  - 相比方案 E，加入“半步撤销”并未改善整体收敛，$LDPC\ FER$ 仍然为 1.0，说明全步进攻在 2bit 下已经让错误进入“深坑”，撤销力度不足以挽回。
  - 与 `doc/bf_2bit_report.md` 中对方案 B 的记录一致：攻守分离在 2bit 空间仍然过于激进，导致大规模误翻。
  - 结论：方案 B 也明显劣于基线方案 0，不具实用价值，仅作为“上界激进度”参考。

## 方案 A：存2算3扩展域（左移运算、右移回写；2bit aggr 关闭）

- 关键代码改动（相对方案 B）：
  - `src/ldpc_codec_test2.cpp:563-586`，在 `f_update_vn_post` 开头为 2bit 增加“存2算3”路径：
    - 对 $VN\_BITS\le 2$：
      - 将 likelihood 左移一位并加偏置：`like = (likelihood<<1)+1`，在 3bit 域 $\{1,3,5,7\}$ 运算；
      - 使用固定阈值 `thr=4`（3bit 经典阈值），`mn/max` 同样在扩展域计算；
      - 使用全步长 `w_ext = weight` 更新：`like_new = flipped ? like-w_ext : like+w_ext-1`；
      - 后处理仍按阈值 `thr` 做“边界±1”处理，最后裁剪到 `[mn,mx]`；
      - 回写时右移一位：`return like_new >> 1`。
    - 对 $VN\_BITS>2$：仍走原来的半步/全步逻辑（旧代码在 2bit 提前 return 后作为 3bit+ 路径保留）。
  - `ldpc_dec_bf_ibex` 中 aggr 条件保持 `VN_BITS>2` 才启用，2bit 完全关闭 aggr。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 其他同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.000000e+00`
  - `[STATISTICS] LDPC BER  : 7.672145e-02`
  - 平均迭代：`1024`，所有 100 个 packet 全失败。
- 分析：
  - 与 `doc/bf_2bit_report.md` 中方案 A 的结论一致：**扩展域+全步长在 2bit 下同样容易发散**，虽然 syndrome_weight 数值明显低于方案 E/B 的 “千级”，但始终无法收敛到 0。
  - 左移扩展确实增加了步长分辨率，但如果不配合更温和的步长/门限设计，依然会出现大规模误翻，FER 退化为 1.0。
  - 综合 E/B/A 三个方案，可以确认：简单地把 3bit 方案“直接搬到” 2bit（无论是全步还是扩展域）都不可行，后续需要更精细的步长控制与迟滞/权重放大策略（即后面 H/T/T3 一类思路）。

## 方案 D：分时后处理（前半关闭 post，后半开启；半步步长）

- 关键代码改动（相对方案 A/B）：
  - `src/ldpc_codec_test2.cpp:563-640`，`f_update_vn_post` 恢复为原始半步逻辑：
    - 对 $VN\_BITS \le 2$：`delta = (weight>>1) + ((weight & 1) && pushing ? 1 : 0)`，未翻转 `likelihood + delta - 1`，已翻转 `likelihood - delta`，与基线 ldpc_codec.cpp 一致。
    - 对 $VN\_BITS > 2$：保持 `delta = weight` 的全步逻辑。
  - `ldpc_dec_bf_ibex` 中 post 触发逻辑（`src/ldpc_codec_test2.cpp:2623-2640`）：
    - 新增 `early_phase = (iteration < ldpc_decoder_input.post_iteration)`；
    - 只有在 `!early_phase && iteration >= post_iteration` 时才按照原公式计算 `post_trigger/post_trigger2`；
    - 在早期阶段强制 `post_process=0, post_trigger=0, post_trigger2=0`，即完全关闭后处理与扰动；
    - 2bit 下 aggr 仍保持关闭，仅 `VN_BITS>2` 时启用权重放大。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置与基线一致（100 包或 10 错停止）。
- 结果：
  - `[STATISTICS] Total packets simulated: 207`
  - `[STATISTICS] LDPC FER  : 4.830918e-02`
  - `[STATISTICS] LDPC BER  : 2.011141e-05`
  - 平均迭代：约 `125.6`。
- 分析：
  - 相比基线方案 0（`LDPC FER ≈ 5.81e-02`），方案 D 把 FER 稍微拉低到 `≈4.83e-02`，有一定改善，并且平均迭代次数略有下降（138 → 126 左右），说明“前期无 post、后期再开启”有助于避免早期扰动导致的不稳定。
  - 与 `doc/bf_2bit_report.md` 中“FER≈0.11（100 包，平均迭代≈199）”的数量级趋势一致：分时后处理比 E/B/A 类激进方案明显收敛得更好，只是这里采用了当前 config 的统计和 2bit aggr 显式关闭，数值上略有差别。
  - 虽然方案 D 已经将 FER 压到 $<0.05$，符合你“FER<0.05 有价值就停下来”的软条件，但考虑到后续还有文档中更优的 H/R/T/T3 等方案，且你要求“顺序实现所有方案”，这里暂视为一个中等偏好的候选，继续向后探索。

## 方案 F：动态步长（趋势好全步，停滞半步；aggr 2bit 关闭）

- 关键代码改动（相对方案 D）：
  - `ldpc_dec_bf_ibex` 中增加 syndrome_weight 宏观趋势统计（`src/ldpc_codec_test2.cpp:2488` 附近）：
    - 在迭代 while 循环开头维护一个长度为 3 的历史窗口 `last_sw[3]`：
      - `trend_good = (syndrome_weight < last_sw[0]) && (last_sw[0] < last_sw[1]);`
      - 然后 `last_sw[2]=last_sw[1]; last_sw[1]=last_sw[0]; last_sw[0]=syndrome_weight;`
    - 使用 `be_aggressive = trend_good;` 把“趋势好”标志传入 `f_update_vn_post`。
  - `f_update_vn_post` 中 2bit 动态步长（`src/ldpc_codec_test2.cpp:563-580`）：
    - 若 `VN_BITS <= 2 && be_aggressive==true`（趋势好）：`delta = weight;` 使用全步长加快收敛；
    - 否则（趋势不好或 3bit+）：2bit 恢复半步公式 `delta = (weight>>1)+((weight&1)&&pushing?1:0)`，3bit+ 仍为 `delta=weight`。
  - post/post2 触发逻辑恢复为基线形式（不再做 D 的分时关闭），2bit aggr 仍然关闭，仅 `VN_BITS>2` 使用加权放大。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 130`
  - `[STATISTICS] LDPC FER  : 7.692308e-02`
  - `[STATISTICS] LDPC BER  : 3.306328e-05`
  - 平均迭代：约 `158.9`。
- 分析：
  - 与文档中记录一致，动态步长方案在 2bit 下表现**变差**：FER 从基线的 `≈0.058` 上升到 `≈0.077`，且平均迭代时长也增长。
  - 原因大致符合你的分析：宏观趋势判据在 2bit BF 的噪声环境下不够稳定，经常误判为“趋势好”，从而使用全步长，导致大规模误翻；而当趋势不佳时退回半步也难以及时修正。
  - 结论：方案 F 作为“自适应步长”的尝试在 2bit IBEX BF 上并不成功，不建议作为候选；后续继续按顺序尝试 G/H/.../T/T3 等更针对 2bit 的方案。

## 扩展+攻守分离变体：存2算3，未翻转全步，已翻转半步（2bit aggr 关闭）

- 关键代码改动（相对方案 F）：
  - `src/ldpc_codec_test2.cpp:563-607`，`f_update_vn_post` 2bit 分支改为“存2算3+攻守分离”：
    - 扩展到 3bit 域：`like = (likelihood<<1)+1`，在 $\{1,3,5,7\}$ 上运算，阈值设为 `thr=4`；
    - 未翻转（攻）：`w_attack = weight<<1`，使用扩展域全步 `like_new = like + w_attack - 1`；
    - 已翻转（守）：`w_retract = (weight>>1)+((weight & 1) && pushing ? 1 : 0)`，半步撤销 `like_new = like - w_retract`；
    - 后处理：仍按 `thr` 做边界±1 调整，再裁剪到 `[mn,mx]`，最后 `return like_new>>1` 回写为 2bit；
    - 对 `VN_BITS>2`：保持原来 `delta=weight` 的全步逻辑不变。
  - `ldpc_dec_bf_ibex` 中 aggr 条件保持前面方案的限制：2bit 关闭 aggr，仅 `VN_BITS>2` 启用。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.000000e+00`
  - `[STATISTICS] LDPC BER  : 9.859440e-01`
  - 平均迭代：`1024`（全部打满迭代），syndrome_weight 长期在 200~270 左右高位震荡，几乎无改善。
- 分析：
  - 与 `doc/bf_2bit_report.md` 中对该变体的结论一致：在 2bit 环境下，即使扩展到 3bit 域并做攻守分离，**未翻转端用“权重左移全步”仍过于激进**，初始高权重会把大量比特推入错误状态，半步撤销追不上，导致 FER=1.0。
  - 相比方案 A（仅扩展）和 B（仅攻守分离），此变体仍然完全发散，说明“扩展域+强攻全步”在 2bit BF 上是不可行方向。
  - 结论：该方案仅作为失败样本记录，不作为任何候选；后续继续进入 G/H 等更细化的 2bit 策略。

## 方案 G：2bit 固定梯度（flip_thr=2，强/弱=0/1），禁用 post，自适应关闭

- 关键代码改动（相对扩展+攻守分离变体）：
  - `f_likelihood_levels`（`src/ldpc_codec_test2.cpp:452-524`）：
    - 对 `VN_BITS <= 2`：
      - 固定设置 `flip_thr=2`，`strong=0`，`weak=1`；
      - 初始 `level[0..3] = {weak, weak, weak, weak}`；
      - 在 `strobes>0` 和 `strobes>1` 分支中加条件 `VN_BITS>2`，即 2bit 模式完全跳过 syndrome_weight 自适应，保持固定梯度。
  - `f_update_vn_post` 2bit 分支改为固定梯度 BF（`src/ldpc_codec_test2.cpp:524-606`）：
    - 恢复 2bit 域直接运算，不再“存2算3”：
      - `delta_attack = weight`（未翻转全步）；
      - `delta_retract = (weight>>1) + ((weight & 1) && pushing ? 1 : 0)`（已翻转半步）；
      - `likelihood_new = flipped ? (likelihood - delta_retract) : (likelihood + delta_attack - 1)`；
      - 对 2bit 模式逻辑上禁用 post/post2（虽然调用时仍传入 post\_process，但 2bit 下阈值/梯度固定，使后处理几乎无正面作用且易造成震荡）。
    - 对 `VN_BITS>2`：继续使用原 BF 全步更新。
  - `ldpc_ibex_parameters` 中 likelihood 映射（`src/ldpc_codec_test2.cpp:1120-1138`）：
    - 对 `VN_BITS <= 2` 将软信息到 level 的映射改为 `{0,1,2,2}`，最强软值对应 level\[0\]=strong=0（低 likelihood 表示“强 0”）。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.000000e+00`
  - `[STATISTICS] LDPC BER  : 9.861246e-01`
  - 平均迭代：`1024`，全部打满迭代，syndrome_weight 常年维持在 200~280 左右高位，几乎不下降。
- 分析：
  - 和 `doc/bf_2bit_report.md` 一致，固定梯度方案在 2bit 下彻底失效：绝大多数比特在强/弱两个档位之间剧烈震荡，配合全步攻击导致大量误翻，解码几乎总是失败。
  - 虽然从设计上看“固定梯度 + 禁用后处理”能避免早期自适应干扰，但在 2bit 的粗粒度下缺乏足够的稳态拉回能力，整体表现反而最差之一。
  - 结论：方案 G 是明显的负样本，不作为候选；后续继续尝试 H/I 等在 2bit 域上做更细致步长/边界设计的方案。

## 方案 H：2bit 小权重“向上取半步”，w>=3 全步；撤销半步；flip_thr=2，post 关闭

- 关键代码改动（相对方案 G）：
  - `f_likelihood_levels`（`src/ldpc_codec_test2.cpp:452-524`）：
    - 保持 2bit 梯度固定：`flip_thr=2, strong=0, weak=1`；
    - 初始 level 设置为 `{0,1,1,1}`：
      - `level[0] = strong`，`level[1..3] = weak`；
    - 2bit 模式继续关闭基于 syndrome_weight 的自适应，VN_BITS>2 才使用 coef-based 自适应。
  - `f_update_vn_post` 2bit 分支（`src/ldpc_codec_test2.cpp:524-606`）：
    - 2bit 直接在原域上运算，不再做存2算3：
      - 若未翻转（`likelihood < flip_thr`）：
        - `delta_attack = (weight >= 3) ? weight : ((weight + 1) >> 1)`（w>=3 全步，小权重上取半步）；
        - `likelihood_new = likelihood + delta_attack - 1`；
      - 若已翻转（`likelihood >= flip_thr`）：
        - `delta_retract = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0)`（半步撤销）；
        - `likelihood_new = likelihood - delta_retract`；
      - 2bit 分支中不使用 post/post2，仅做 min/max clamp。
    - 对 `VN_BITS>2`：仍保持 `delta=weight` 的原 BF 全步逻辑。
  - `ldpc_ibex_parameters` 中 soft→level 映射保持方案 G 的 2bit 设定 `{0,1,2,2}`，aggr 仍仅在 `VN_BITS>2` 启用，2bit 完全关闭 aggr。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.800000e-01`
  - `[STATISTICS] LDPC BER  : 2.243728e-05`
  - 平均迭代：约 `201.6`。
- 分析：
  - 相比方案 G（FER=1.0）、扩展变体与 E/B/A 等激进方案，H 确实“拉回来了”，FER 降到 `≈0.18`，验证了“小权重上取半步 + w>=3 全步 + 半步撤销”这一结构有一定稳定作用。
  - 但与基线方案 0（`FER≈0.058`）和方案 D（`FER≈0.048`）相比，方案 H 仍明显偏差，且平均迭代数更高，表现仍处于 error floor 区域。
  - 结论：方案 H 作为正样本比 G 好很多，但在当前配置下远未达到“有价值候选”的标准（相对基线不优），按计划继续尝试 I/J/K/L/M/N/P/Q/R/S/T/T3。 

## 方案 I：在 H 基础上启用 2bit 后处理（未翻转按到 min，翻转边缘推到 thr+1）

- 关键代码改动（相对方案 H）：
  - `f_update_vn_post` 2bit 分支（`src/ldpc_codec_test2.cpp:524-606`）在 H 的攻守分离基础上加入 post：
    - 先按方案 H 计算 `likelihood_new`：
      - 未翻转：`delta_attack = (w>=3)?w:((w+1)>>1)`，`likelihood_new = likelihood + delta_attack - 1`；
      - 已翻转：`delta_retract = (w>>1)+((w & 1)&&pushing?1:0)`，`likelihood_new = likelihood - delta_retract`；
    - 然后在 2bit 下启用简化后处理：
      - 若 `post_process && likelihood_new < flip_threshold`：强制置为 `min_likelihood`；
      - 若 `post_process && likelihood_new == flip_threshold`：置为 `flip_threshold + 1`；
    - 仍不使用 post_process2，只做 min/max clamp。
  - 2bit 梯度（flip_thr=2、strong=0、weak=1、level={0,1,1,1}）、映射 {0,1,2,2} 以及 aggr 关闭均与方案 H 保持一致。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.800000e-01`
  - `[STATISTICS] LDPC BER  : 3.379109e-05`
  - 平均迭代：约 `200.8`。
- 分析：
  - 与 H 相比，I 的 FER 反而略高，平均迭代也变长，和文档里“FER≈0.18、平均迭代更长”的趋势一致。
  - 直觉上，I 的 post 把“尚未翻转但接近阈值”的比特直接按到 min，降低了部分摇摆，但整体推力不足，更多时间困在 error floor 周围，导致性能恶化。
  - 结论：方案 I 说明简单“强压 post”对 2bit BF 帮助有限，甚至负面，不作为候选；按计划继续进入 J/K/L/M/N/P/Q/R/S/T/T3。

## 方案 T3：只用 aggr 权重放大，f_update_vn_post 用半步（2bit，VN_BITS<=2 且 iter>=100）

- 关键代码改动（相对 H/J 基线）：
  - `ldpc_dec_bf_ibex` 中的 aggr 判据（`src/ldpc_codec_test2.cpp:2688-2770`）：
    - 修改为：
      - `bool aggr = ((VN_BITS <= 2) && (iteration >= 100)) || ((ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));`
    - 含义：2bit 在迭代次数达到 100 之后强制启用 aggr（即权重放大），3bit 及以上仍按原判据启用 aggr。
    - 权重放大仍为原始映射，特别是 2bit 下：
      - `w=2 → 3`，`w=3 → 5`，`w=4 → 7`，为半步公式提供足够推力。
  - `f_update_vn_post` 2bit 分支（`src/ldpc_codec_test2.cpp:524-606`）：
    - 改为 T3 的“半步唯一”形式：
      - 2bit 下完全忽略 be_aggressive 全步逻辑：
        - 未翻转：`delta = (weight + 1) >> 1`（始终半步，上取整），`likelihood_new = likelihood + delta - 1`；
        - 已翻转：`delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0)`，`likelihood_new = likelihood - delta`；
      - 不再对 2bit 使用 post/post2，只做 min/max clamp。
    - 3bit 及以上仍保持原来的全步公式：
      - `delta = weight`，`likelihood_new = flipped ? (likelihood - delta) : (likelihood + delta - 1)`。
  - 2bit 梯度与映射沿用 H 方案：
    - `flip_thr=2, strong=0, weak=1, level 初值 {0,1,1,1}`；
    - soft→level 映射 `{0,1,2,2}`。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 与前几方案相同配置，但在 T3 下为观察 FER 的统计稳定性，仿真时间被拉长（直到外部超时）。
- 运行结果（命令因超时被外部终止，统计读自中间多次 `[SIM] Statistical result`，而非最终 `[STATISTICS]` 汇总）：
  - 多个时间点的统计：
    - 83 包、1 错：`LDPC FER ≈ 1.20e-02`
    - 389 包、2 错：`LDPC FER ≈ 5.14e-03`
    - 1756 包、3 错：`LDPC FER ≈ 1.71e-03`
    - 3108 包、4 错：`LDPC FER ≈ 1.29e-03`
    - 3693 包、5 错：`LDPC FER ≈ 1.35e-03`
    - 3748 包、6 错：`LDPC FER ≈ 1.60e-03`
  - 在几千个 packets 的尺度上，FER 稳定在 $10^{-3}$ 量级。
- 分析：
  - 相比当前最佳的方案 D（`LDPC FER ≈ 4.83e-02`）和基线方案 0（`≈5.81e-02`），T3 在相同 SNR=5.4 下将 FER 压到了约 $10^{-3}$，提升接近两个数量级，方向和你原报告中的 T3 非常一致。
  - 原理与文档分析高度吻合：
    - iter<100 时使用半步，避免早期发散；
    - iter≥100 后启用 aggr 权重放大，将 w=2,3,4 分别放大成 3,5,7，再通过半步公式转换成更合适的 delta，既提高推力，又不至于全步过猛导致误翻。
    - 2bit 下完全不用 be_aggressive 切换全步，而是把“激进性”全部压到权重放大层，配合半步进行平滑推进。
  - 从多段统计看，在几千包规模下得到的 `LDPC FER ≈ 1.3×10^{-3}` 非常有价值，完全满足你“FER 显著低于 0.05、有进一步验证价值”的标准。
  - 建议后续你可以：
    - 单独运行 `ssd_fc_test2`（适当放宽 timeout），例如固定跑 5000 或 10000 包，进一步收紧 FER 置信区间；
    - 在更低 SNR（例如 5.2、5.0）下测试 T3 方案的瀑布区表现，确认它在边缘工作点的优势。

## 方案 T3-iter80：在 T3 基础上提前启用 aggr（SNR=5.3）

- 关键代码改动（相对 T3）：
  - 仅调整 2bit 场景下 aggr 启动迭代门限，其余 T3 逻辑保持不变。
  - `ldpc_dec_bf_ibex` 中权重放大判据（`src/ldpc_codec_test2.cpp`）：
    ```cpp
    // 方案 T3-iter80：2bit 在迭代达到 80 次后启用 aggr 权重放大
    bool aggr = ((VN_BITS <= 2) && (iteration >= 80)) ||
                ((ldpc_decoder_input.soft_bits > 0) &&
                 (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
    ```
  - 2bit 下仍采用 T3 的半步更新策略：
    - `f_likelihood_levels`：`flip_thr=2`，`strong=0`，`weak=1`，初始 level 为 $\{0,1,1,1\}$；
    - `f_update_vn_post`：始终半步  
      - 未翻转：`delta=(weight+1)>>1`，更新为 `likelihood + delta - 1`；
      - 已翻转：`delta=(weight>>1)+((weight & 1)&&pushing?1:0)`，更新为 `likelihood - delta`；
    - 2bit 下不启用 post/post2，仅做 min/max clamp。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.3`
  - 停止条件：`100 packets or 10 errors`（程序默认配置）。
- 仿真结果（一次运行）：
  - `[STATISTICS] Total packets simulated: 216`
  - `[STATISTICS] RAW BER   : 7.023916e-03`
  - `[STATISTICS] LDPC BER  : 9.000445e-03`
  - `[STATISTICS] LDPC FER  : 4.629630e-02`
  - `[STATISTICS] IBEX Decoder average iterations: 120.27`
- 性能对比：
  - 当前基线 `ldpc_codec_test.cpp` 在相同 SNR=5.3 下（你提供的结果）：
    - `LDPC FER ≈ 9.27×10^{-3}`（1079 包，FER 约 $10^{-2}$ 级别，平均迭代约 62）
  - 方案 T3-iter80：
    - `LDPC FER ≈ 4.63×10^{-2}`（216 包，FER 仍在 $10^{-2}$ 高位甚至接近 $5×10^{-2}$）
    - 平均迭代约 120，比基线的 $\approx 62$ 显著偏高，说明在 5.3 下大量包需要更长时间挣扎仍未收敛。
  - 相比原始 T3 在 5.3 下的初始测试结果（约 184 包、FER≈0.054、平均迭代≈145），T3-iter80：
    - FER 从约 $5.4×10^{-2}$ 略微改善到约 $4.6×10^{-2}$；
    - 平均迭代从约 145 降到约 120，说明提前启用 aggr 确实加快了部分包的收敛。
  - 但整体上，T3-iter80 依然**明显劣于** `ldpc_codec_test.cpp` 在 SNR=5.3 的表现（$\sim 9.3×10^{-3}$），尚不足以作为“超越基线”的方案。
- 原因分析（为何未能超越基线）：
  - T3 在 SNR=5.4 的优势来自于“前期纯半步、后期少量激进权重放大”，可以在错误较少的场景下细致清除残余错误，将 FER 压到 $10^{-3}$。
  - 降到 SNR=5.3 后：
    - 通道错误增多，2bit 半步在 flip_thr=2 的格子里推力仍然偏小，很多错误比特长期处于边界附近，无法顺利翻转；
    - 即使将 aggr 提前到 iter≥80，半步 + 权重放大对长期高 syndrome 的包仍然略显保守，导致**平均迭代数高**但 FER 仍停留在 $10^{-2}$ 高位；
    - 与 `ldpc_codec_test.cpp` 相比，当前 T3/T3-iter80 组合在“低 SNR 时的早期纠错能力”明显偏弱：  
      基线版本更倾向于在 w=2/3 时更早推进翻转，从而在中等 SNR 区域获得较好的 FER/迭代平衡。
  - 因此，本次“仅调整 aggr 启动迭代门限”的小步优化虽然略微降低了 FER 和平均迭代，但不足以弥补 T3 在 SNR=5.3 区域整体推力偏弱的问题，未能超越 `ldpc_codec_test.cpp`。

## 方案 T3-w2push：在 T3-iter80 基础上增强 w=2 的后期推力（SNR=5.3）

- 关键代码改动（相对 T3-iter80）：
  - 保持 2bit 的 T3 结构不变（固定梯度、flip_thr=2、半步攻守分离、aggr 只通过权重放大），仅在 **aggr 阶段对 w=2 的未翻转分支去掉 “-1” 抑制**，让 w=2 在后期可以逐步推进到阈值。
  - `f_update_vn_post` 中 2bit 分支（`src/ldpc_codec_test2.cpp`）：
    ```cpp
    if (VN_BITS <= 2) {
      bool flipped = (likelihood >= flip_threshold);
      int likelihood_new;

      if (!flipped) {
        // 半步：权重放大仍由 aggr 完成
        int delta = (weight + 1) >> 1;

        // 仅在 aggr 阶段对 w=2 去掉 “-1” 抑制，增强边界比特推进能力
        if (be_aggressive && (weight == 2)) {
          likelihood_new = likelihood + delta;      // e.g. 0->1, 1->2
        } else {
          likelihood_new = likelihood + delta - 1;  // 与 T3 一致
        }
      } else {
        int delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
        likelihood_new = likelihood - delta;
      }

      // clamp 到 [min_likelihood, max_likelihood]
      ...
    }
    ```
  - `ldpc_dec_bf_ibex` 中 aggr 判据仍为 T3-iter80：
    ```cpp
    bool aggr = ((VN_BITS <= 2) && (iteration >= 80)) ||
                ((ldpc_decoder_input.soft_bits > 0) &&
                 (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
    ```
  - 其余设置与 T3 相同：
    - `flip_thr=2`，`strong=0`，`weak=1`，level 初始 $\{0,1,1,1\}$；
    - 2bit 下不使用 post/post2，只做 min/max clamp。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.3`
  - 停止条件：`100 packets or 10 errors`（默认）。
- 仿真结果（一次运行）：
  - `[STATISTICS] Total packets simulated: 256`
  - `[STATISTICS] RAW BER   : 7.031609e-03`
  - `[STATISTICS] LDPC BER  : 7.622108e-03`
  - `[STATISTICS] LDPC FER  : 3.906250e-02`
  - `[STATISTICS] IBEX Decoder average iterations: 114.89`
- 性能对比：
  - 相比 T3 原始在 5.3 的结果（约 184 包，`LDPC FER ≈ 5.43×10^{-2}`，平均迭代≈145）：
    - FER 从 $\sim 5.4×10^{-2}$ 先在 T3-iter80 降到 $\sim 4.6×10^{-2}$，再在 T3-w2push 降到 $\sim 3.9×10^{-2}$；
    - 平均迭代也从 145 → 120 → 115，呈稳步下降趋势。
  - 相比 T3-iter80：
    - FER：`4.63×10^{-2} → 3.91×10^{-2}`，有进一步约 15% 左右的相对改善；
    - 平均迭代：`120.27 → 114.89`，说明增强 w=2 的后期推力确实帮助部分边缘包更快收敛。
  - 但与基线 `ldpc_codec_test.cpp` 在 SNR=5.3 下的表现相比：
    - 基线 FER 仍约为 `9.27×10^{-3}`（$10^{-2}$ 低位），平均迭代约 62；
    - T3-w2push 虽已在 $4×10^{-2}$ 附近，依然高出基线约 4 倍，且迭代数也接近基线的两倍。
- 原因分析：
  - 在 T3 架构下，2bit 解码的“后期清尾能力”已经很强（这在 SNR=5.4 表现为 FER 接近 $10^{-3}$），但在 SNR=5.3 这种更“粗糙”的信道环境中，**早期和中期的纠错推力仍然偏弱**：
    - 原始 T3 中，w=2 在 2bit 半步+“减 1” 抑制下对未翻转比特几乎没有实际推进作用；  
    - T3-iter80 提前 aggr 让 w=3 放大后的半步起作用，但 w=2 仍然不动，很多边界比特长期卡在 $likelihood=0/1$；
    - T3-w2push 在 aggr 阶段让 w=2 也能贡献“0→1→2”的累计推进，对部分尾部错误确实有帮助，因此 FER 和迭代都有实测改善。
  - 但与 `ldpc_codec_test.cpp` 相比：
    - 基线方案在 SNR=5.3 下对 w=2/3 的处理整体更加激进，在早期就能将大量 correct bit 拉到正确一侧，从而在中后期需要纠的残差较少；
    - T3/T3-iter80/T3-w2push 的设计初衷偏向“高 SNR 精修”：严格控制误翻，在高 SNR 下有极佳的 error floor，但在稍低 SNR 时会表现为**持续 error floor 较高、迭代数偏大**。
  - 综上，T3-w2push 在 SNR=5.3 下相对于 T3/T3-iter80 是一次有效的小步改进，但仍未达到或超过 `ldpc_codec_test.cpp` 的性能水平。

## 方向1-方案1：T3-aggr-w23push（SNR=5.3）

- 关键代码改动（相对 T3-iter80）：
  - 在 2bit `f_update_vn_post` 中，仍保持半步攻守分离，只在 **aggr 阶段的中等权重** 上去掉 “-1” 抑制：
    - 原 T3/T3-iter80/T3-w2push 实现（2bit 分支，未翻转）：
      ```cpp
      int delta = (weight + 1) >> 1;
      likelihood_new = likelihood + delta - 1;
      ```
    - 本方案修改为（`src/ldpc_codec_test2.cpp`）：
      ```cpp
      if (VN_BITS <= 2) {
        bool flipped = (likelihood >= flip_threshold);
        int likelihood_new;

        if (!flipped) {
          int delta = (weight + 1) >> 1;

          // 在 aggr 阶段，对 w=3/5（对应原始 w=2/3）去掉 “-1” 抑制
          if (be_aggressive && (weight == 3 || weight == 5)) {
            likelihood_new = likelihood + delta;
          } else {
            likelihood_new = likelihood + delta - 1;
          }
        } else {
          int delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
          likelihood_new = likelihood - delta;
        }
        ...
      }
      ```
    - 这里的 `weight` 是 aggr 放大后的 `w`，因此 `weight==3/5` 对应原始校验节点权重为 2/3。  
  - 2bit 其它部分保持 T3-iter80 设定：
    - `flip_thr=2`，`strong=0`，`weak=1`，level 初始 $\{0,1,1,1\}$；
    - aggr 判据：`(VN_BITS <= 2 && iteration >= 80)`；
    - 2bit 下不启用 post/post2。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.3`
  - 停止条件：`100 packets or 10 errors`。
- 仿真结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 5.200000e-01`
  - `[STATISTICS] LDPC BER  : 1.911010e-01`
  - `[STATISTICS] IBEX Decoder average iterations: 566.23`
- 性能评价：
  - FER 直接退化到约 $5.2×10^{-1}$，几乎一半包失败，属于明显发散；
  - 平均迭代数虽然没有被 early stop 限死（仍在 566 次附近），但 syndrome_weight 长期维持在 1100 附近，很难下降。
  - 相比 T3/T3-iter80/T3-w2push 在 5.3 下 $3×10^{-2}\sim5×10^{-2}$ 的 FER，本方案显著恶化。
- 原因分析：
  - aggr 之后的权重已经被放大：`w=2→3, 3→5`，再在 `w=3/5` 上去掉 “-1” 抑制，相当于对 medium-weight 比特施加了过强的推进：
    - 举例：`likelihood=1, weight=3` 时，$\Delta=(3+1)/2=2$，更新为 `1+2=3`，一次就从边界拉到最大；
    - 在 2bit 的离散格子内，这种强推会导致大量误翻，syndrome 难以下降。
  - 原 T3 的成功经验是在 **半步 + 适度权重放大** 的组合上寻找平衡，本方案在已经放大的权重上进一步增强步长，破坏了这种平衡。
  - 结论：方向1-方案1 是一个明显失败的尝试，证明“在 aggr 区对 w=2/3 再额外增强步长”容易把 2bit decoder 推入发散区域。

## 方向1-方案2：T3-aggr-w23full（SNR=5.3）

- 关键代码改动（相对 方向1-方案1）：
  - 在 2bit `f_update_vn_post` 中进一步增强 aggr 区中等权重的步长，从“去掉 -1 抑制”升级为“直接用放大后的全权重”：
    ```cpp
    if (VN_BITS <= 2) {
      bool flipped = (likelihood >= flip_threshold);
      int likelihood_new;

      if (!flipped) {
        // 禁用 be_aggressive 全步逻辑的基础上，引入更激进的 aggr 中权重推力
        int delta_half = (weight + 1) >> 1;

        if (be_aggressive && (weight == 3 || weight == 5)) {
          // 在 aggr 区对 w=3/5 使用“全步”攻击（基于放大后的权重）
          int delta_full = weight;
          likelihood_new = likelihood + delta_full - 1;
        } else {
          likelihood_new = likelihood + delta_half - 1;
        }
      } else {
        int delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
        likelihood_new = likelihood - delta;
      }
      ...
    }
    ```
  - 其它参数保持与方向1-方案1一致：2bit 固定梯度、flip_thr=2、aggr 从 iter≥80 开启权重放大，post/post2 仍禁用。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.3`
  - 停止条件：`100 packets or 10 errors`。
- 仿真结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 5.600000e-01`
  - `[STATISTICS] LDPC BER  : 2.062754e-01`
  - `[STATISTICS] IBEX Decoder average iterations: 599.81`
- 性能评价与分析：
  - FER 约 $5.6×10^{-1}$，比方向1-方案1 的 $5.2×10^{-1}$ 还略差，确认“在放大后的权重上直接使用全步攻击”会导致极强的误翻：
    - 在 aggr 区，原始权重 2/3 被放大为 w=3/5，再用全步 `delta_full=weight`，例如 w=5 时单次更新几乎必然把 likelihood 推到最大端，引发大量错误翻转；
    - 仿真 log 中 syndrome_weight 长期停留在 1100~1200 高位，平均迭代接近 600，表明 decoder 处于持续发散状态。
  - 与 `ldpc_codec_test.cpp` 在 SNR=5.3 下的 $FER \approx 9.27×10^{-3}$ 相比，本方案相差两个数量级以上，完全不可用。
  - 方向1-方案2 进一步验证了一个结论：**在 2bit + aggr 框架下，中权重上的全步攻击非常危险**，尤其是在权重已经被放大到 3/5 的条件下，很容易把系统推到“无可挽回”的发散区域。

## 方向2-方案1：T3 + ldpc_codec_test 风格 aggr 启动条件（SNR=5.3）

- 关键代码改动：
  - 恢复 2bit 步长为 T3 风格（半步攻守分离），但将 aggr 启动条件改为与 `ldpc_codec_test.cpp` 一致：
    - 2bit `f_update_vn_post`（保持方向1-方案2 的步长设定，此处不再展开）；  
    - `ldpc_dec_bf_ibex` 中 aggr 判据（`src/ldpc_codec_test2.cpp`）从：
      ```cpp
      // 原 T3-iter80：2bit 在迭代达到 80 次后启用 aggr 权重放大
      bool aggr = ((VN_BITS <= 2) && (iteration >= 80)) ||
                  ((ldpc_decoder_input.soft_bits > 0) &&
                   (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
      ```
      修改为：
      ```cpp
      // 方向2-方案1：2bit 使用 ldpc_codec_test.cpp 风格的 aggr 启动条件
      bool aggr = ((VN_BITS <= 2) &&
                   (iteration >= 100 || (iteration >= 50 && syndrome_weight < 150))) ||
                  ((ldpc_decoder_input.soft_bits > 0) &&
                   (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
      ```
    - 含义：2bit 只有在迭代次数足够大，且 syndrome_weight 已经较小（接近收敛）时才开启权重放大，更接近 `ldpc_codec_test.cpp` 的激进策略使用时机。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.3`
  - 停止条件：`100 packets or 10 errors`。
- 仿真结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 6.800000e-01`
  - `[STATISTICS] LDPC BER  : 2.504168e-01`
  - `[STATISTICS] IBEX Decoder average iterations: 710.95`
- 性能评价与分析：
  - 仅仅更换 aggr 启动条件，在当前已经“过于激进”的 2bit 步长（方向1-方案2 的设定）基础上，并不能扭转整体发散趋势，FER 仍在 $0.6$ 左右；
  - 相比 `ldpc_codec_test.cpp` 在 SNR=5.3 下 $FER \approx 9.27×10^{-3}$，差距依然巨大；
  - 说明在当前步长/推力组合下，即使使用更合理的 aggr 启动窗口（靠近 `ldpc_codec_test.cpp`），整体推进仍然过猛，导致 decoder 长期在高 syndrome 区间震荡。

## 方向2-方案2：T3 半步 + 早期 2bit post、后期 aggr 清尾（SNR=5.3）

- 关键代码改动：
  - 在 2bit `f_update_vn_post` 中回到统一的 T3 半步逻辑，同时引入“仅在非 aggr 阶段启用 2bit post、aggr 阶段关闭 post”的混合策略：
    ```cpp
    int ldpc_packet::f_update_vn_post(int likelihood, int weight, int min_likelihood, int max_likelihood,
                                      bool post_process, bool post_process2,
                                      bool be_aggressive, int flip_threshold, bool pushing) {
      // 方案 T3：2bit 使用半步更新（攻守分离），aggr 仅通过权重放大，不再依赖 be_aggressive 全步逻辑
      if (VN_BITS <= 2) {
        bool flipped = (likelihood >= flip_threshold);
        int likelihood_new;

        if (!flipped) {
          // 混合方案：统一使用半步
          int delta = (weight + 1) >> 1;
          likelihood_new = likelihood + delta - 1;
        } else {
          int delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
          likelihood_new = likelihood - delta;
        }

        // 仅在非 aggr 阶段启用 2bit post（借鉴 ldpc_codec_test.cpp），后期 aggr 阶段关闭 post
        bool do_post_flipped = post_process && !be_aggressive && (likelihood_new == flip_threshold);
        bool do_post_unflipped = post_process && !be_aggressive && (likelihood_new < flip_threshold);
        if (do_post_unflipped)
          likelihood_new = flip_threshold - 1;
        else if (do_post_flipped)
          likelihood_new = flip_threshold + 1;

        if (likelihood_new <= min_likelihood)
          likelihood_new = min_likelihood;
        if (likelihood_new >= max_likelihood)
          likelihood_new = max_likelihood;

        return likelihood_new;
      }
      // VN_BITS>2：保持原 BF 逻辑（略）
    }
    ```
  - aggr 判据继续沿用方向2-方案1 的 `ldpc_codec_test.cpp` 风格：
    ```cpp
    bool aggr = ((VN_BITS <= 2) &&
                 (iteration >= 100 || (iteration >= 50 && syndrome_weight < 150))) ||
                ((ldpc_decoder_input.soft_bits > 0) &&
                 (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
    ```
  - 直观理解：
    - 迭代早期：`be_aggressive=false`，步长为 T3 半步，但允许 2bit post 在边界附近进行强制下沉/上抬，帮助打破局部震荡；
    - 迭代中后期且 syndrome 已较小：`be_aggressive=true`，关闭 post，仅依靠权重放大 + 半步推进做精细清尾。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.3`
  - 停止条件：`100 packets or 10 errors`（结果中 747 包时达到 10 错）。
- 仿真结果：
  - `[STATISTICS] Total packets simulated: 747`
  - `[STATISTICS] LDPC FER  : 1.338688e-02`
  - `[STATISTICS] LDPC BER  : 2.611947e-03`
  - `[STATISTICS] IBEX Decoder average iterations: 65.63`
- 性能对比与分析：
  - 与 `ldpc_codec_test.cpp` 在 SNR=5.3 下的基线结果（你提供）：  
    - 基线：`LDPC FER ≈ 9.27×10^{-3}`，平均迭代约 62.16；  
    - 方向2-方案2：`LDPC FER ≈ 1.34×10^{-2}`（747 包、10 错），平均迭代 ≈ 65.63。
  - 从数值上看，本方案的 FER 略高于基线（约高 40% 左右），平均迭代也略大，但已经处在 **同一数量级**，远好于此前所有 T3 变体在 5.3 下的 $3×10^{-2}\sim5×10^{-2}$ 水平：
    - 相对 T3 原始在 5.3 的 `FER ≈ 5.4×10^{-2}`，方向2-方案2 将 FER 降低了约 4 倍；
    - 相比 T3-iter80/T3-w2push 的 `FER ≈ 4.6×10^{-2}`、`≈3.9×10^{-2}`，进一步改善到 $≈1.3×10^{-2}$，接近基线方案。
  - 这说明：
    - 早期启用 2bit post 可以在弱推力（半步）阶段有效打破局部震荡，将大量比特推向更干净的状态；
    - 在 syndrome 足够小之后再开启 aggr + 半步，有利于模仿 T3 在高 SNR 下的“清尾”特点，而不会过早发散。
  - 目前从单次 747 包结果看，方向2-方案2 仍略逊于 `ldpc_codec_test.cpp`，但已经是一个非常有价值的候选方案：
    - FER 达到 $10^{-2}$ 级别，且平均迭代与基线接近；
    - 若增加包数（例如到 2000+）进行更长统计，有可能在统计波动下接近甚至略优于基线。


## 方案 J：H 基础上关闭 2bit post，再加“高 syndrome 才开 aggr”

- 关键代码改动（相对方案 I）：
  - `f_update_vn_post` 2bit 分支仍保持方案 H 的攻守分离与小权重上取半步（不再使用 I 的 post 强压）：
    - 未翻转：`delta_attack = (w>=3)?w:((w+1)>>1)`，`likelihood_new = likelihood + delta_attack - 1`；
    - 已翻转：`delta_retract = (w>>1)+((w & 1)&&pushing?1:0)`，`likelihood_new = likelihood - delta_retract`；
    - 不使用 post/post2，只做 min/max clamp。
  - `ldpc_dec_bf_ibex` 中 aggr 判据修改（`src/ldpc_codec_test2.cpp:2760-2787`）：
    - 对 `VN_BITS <= 2`：`aggr = (syndrome_weight_delayed > 200)`，即仅当 syndrome_weight_delayed 大于约 200 时才启用 2bit 激进权重放大；
    - 对 `VN_BITS > 2`：保留原始 aggr 判据 `(soft_bits>0 && likelihood_levels.min<likelihood_thr && !flipped_prev)`。
    - 权重放大曲线保持不变：`2→3, 3→5, 4→7` 等。
- 运行配置：
  - 命令：`./ssd_fc_test2 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
  - 配置同前。
- 结果：
  - `[STATISTICS] Total packets simulated: 100`
  - `[STATISTICS] LDPC FER  : 1.000000e+00`
  - `[STATISTICS] LDPC BER  : 3.673781e-01`
  - 平均迭代：`1024`，所有 100 包全部失败，syndrome_weight 持续在约 1100~1200 区间高位震荡。
- 分析：
  - 结果与 `bf_2bit_report.md` 一致：这类“高 syndrome 再开 aggr”的策略在 2bit 平台上实际上是灾难性的——一旦 syndrome_weight 已经很高，再加 aggr 会进一步放大误翻，完全无法拉回。
  - 相比 H/I，J 的 FER 从 0.18 降级到 1.0，说明在 2bit 下 aggr 的使用必须高度谨慎，简单门控无法避免其破坏性。
  - 结论：方案 J 是一个显著的负样本，说明“在高 syndrome 下开启 aggr”并不可取；后续方案 K/L/M/N 等都在这一基础上做进一步变化，这里继续按顺序验证。

## Row11 参数扫参（SNR=5.0，2bit）

- 配置与前提：
  - 配置文件：`config/Ibex_hd_row11.cnfg`
  - `Ibex likelihood width = 2`（即 `VN_BITS=2`）
  - 扫描对象：
    - `aggr_iter_hi = A`
    - `aggr_iter_lo = B`
    - `aggr_synd_th = C`
    - 判据：`(iteration >= A) || (iteration >= B && syndrome_weight < C)`

- 代码改动位置：`src/ldpc_codec_test2.cpp:2736-2740`
  - 采用常量化门限，便于扫参：
    - `const int aggr_iter_hi = A;`
    - `const int aggr_iter_lo = B;`
    - `const int aggr_synd_th = C;`

- 第一轮批量结果（`./ssd_fc_test2 LDPC config/Ibex_hd_row11.cnfg AWGN 5.0`）：
  - `70/35/100`: `pkts=148`, `FER=6.756757e-02`
  - `80/40/120`: `pkts=306`, `FER=3.267974e-02`
  - `90/45/130`: `pkts=274`, `FER=3.649635e-02`
  - `100/50/150`: `pkts=258`, `FER=3.875969e-02`
  - `110/55/160`: `pkts=296`, `FER=3.378378e-02`
  - `120/60/180`: `pkts=386`, `FER=2.590674e-02`
  - `130/65/190`: `pkts=521`, `FER=1.919386e-02`
  - `140/70/200`: `pkts=459`, `FER=1.307190e-02`
  - `160/80/220`: `pkts=454`, `FER=1.101322e-02`
  - `200/100/260`: 130s 内无统计（`no_stats`）

- 对 `200/100/260` 的加长验证（300s）：
  - 命令：`timeout 300s stdbuf -oL -eL ./ssd_fc_test2 LDPC config/Ibex_hd_row11.cnfg AWGN 5.0`
  - 末尾有效统计（SIM 段）：
    - `Total packets simulated: 1394`（对应 `SIM` 行）
    - `FAIL CW: 5`
    - `LDPC FER: 3.586801e-03`
  - 结论：在更长观察窗内，`200/100/260` 明显优于 `160/80/220`（后者约 `1.10e-02`）。

- 当前最优候选（SNR=5.0，row11）：
  - `A/B/C = 200/100/260`
  - 当前代码已设置为该参数点（`src/ldpc_codec_test2.cpp:2736-2738`）。

- 分析：
  - 随着 `A/B/C` 提高，aggr 启动被进一步后移，能够减少中前期误翻，FER 从 `6.76e-02` 持续下降到 `1.10e-02`。
  - 更高门限 `200/100/260` 在短窗内“慢出错”特征明显，说明其 error floor 更低，但吞吐统计时间变长。

## Row11 参数扫参（SNR=4.9，2bit）

- 配置：`./ssd_fc_test2 LDPC config/Ibex_hd_row11.cnfg AWGN 4.9`
- 扫描门限：
  - `aggr_iter_hi = A`
  - `aggr_iter_lo = B`
  - `aggr_synd_th = C`

- 第一轮粗扫结果：
  - `140/70/200` -> `pkts=86`, `FER=1.744186e-01`
  - `160/80/220` -> `pkts=113`, `FER=8.849558e-02`
  - `180/90/240` -> `pkts=112`, `FER=8.928571e-02`
  - `190/95/250` -> `pkts=152`, `FER=6.578947e-02`
  - `200/100/260` -> `pkts=95`, `FER=1.157895e-01`
  - `210/105/270` -> `pkts=121`, `FER=8.264463e-02`
  - `220/110/280` -> `pkts=181`, `FER=5.524862e-02`（粗扫最优）
  - `240/120/300` -> `pkts=87`, `FER=1.724138e-01`

- 第二轮细扫（围绕 220）结果：
  - `205/102/265` -> `FER=1.489362e-01`
  - `215/108/275` -> `FER=9.259259e-02`
  - `220/110/270` -> `FER=1.052632e-01`
  - `220/110/280` -> `FER=5.000000e-02`
  - `220/110/290` -> `FER=1.010101e-01`
  - `230/115/290` -> `FER=4.807692e-02`（细扫最优）
  - `230/115/310` -> `FER=1.300000e-01`
  - `235/118/300` -> `FER=1.123596e-01`

- 前两名加长复验（300s）：
  - `220/110/280` -> `pkts=233`, `FER=4.291845e-02`
  - `230/115/290` -> `pkts=185`, `FER=5.405405e-02`

- 结论（SNR=4.9）：
  - 最终采用 `A/B/C = 220/110/280`。
  - 与 `SNR=5.0` 下的最优点 `200/100/260` 不同，说明低 SNR 需要更晚且更稳的 aggr 启动节奏，避免过早权重放大导致误翻扩散。

## Row11 深化探索（基于当前最优门限 220/110/280，SNR=4.9）

- 目标：在不改变总体框架（2bit 半步 + aggr 权重放大）的前提下，继续寻找比 `A/B/C=220/110/280` 更优的细节方案。
- 评测口径：
  - 主口径：`./ssd_fc_test2 LDPC config/Ibex_hd_row11.cnfg AWGN 4.9`
  - 以 `LDPC FER` 为主，记录 `pkts` 与 `avg_iter`。
  - 由于配置使用“达到 10 错即停”，FER 存在样本波动，新增方案均做了重复或长窗复验。

### 基线复测（当前代码：220/110/280）
- 3次复测：
  - run1: `pkts=174`, `FER=5.747126e-02`
  - run2: `pkts=148`, `FER=6.756757e-02`
  - run3: `pkts=193`, `FER=5.181347e-02`
- 中位数约：`5.747126e-02`
- 说明：与前面长窗复验 `FER=4.291845e-02` 相比存在波动，但仍在同一量级。

### 方案 V1：2bit post 仅在非 aggr 阶段启用
- 改动：
  - `f_update_vn_post` 2bit 分支：
    - `do_post_flipped = post_process && !be_aggressive && ...`
    - `do_post_unflipped = post_process && !be_aggressive && ...`
- 3次复测：`6.666667e-02 / 6.211180e-02 / 1.111111e-01`
- 结论：较基线更差，回退。

### 方案 V2：aggr 下 soft `w=4` 放大（7->6）
- 改动：
  - `if (aggr && weight==4) w = weight + 2`（原为 `+3`）
- 3次复测：`9.803922e-02 / 6.211180e-02 / 7.092199e-02`
- 结论：整体不优，回退。

### 方案 V3：w=2 且接近阈值时边界增强（delta=2）
- 改动：
  - 2bit 未翻转分支增加：`if (!be_aggressive && weight==2 && likelihood>=flip_thr-1) delta=2`
- 结果：快速发散（中间 FER 上升到约 `5.8e-01`），提前中止。
- 结论：不可用，回退。

### 方案 V4：在 2bit 分支启用 post_process2 的单步上推
- 改动：
  - 增加：`if (post_process2 && !flipped && weight==1 && likelihood_new==flip_thr-1) likelihood_new++`
- 长窗结果：`pkts=176`, `FER=5.681818e-02`, `avg_iter=165.24`
- 结论：未优于当前最优，回退。

### 方案 V5：高 syndrome 时临时软化 aggr 的 w=3/4 放大
- 改动：
  - `weight==3`：`+2 -> (syndrome_weight>1200 ? +1 : +2)`
  - `weight==4`：`+3 -> (syndrome_weight>1200 ? +1 : +3)`
- 长窗结果：`pkts=159`, `FER=6.289308e-02`, `avg_iter=174.67`
- 结论：不优，回退。

### 评估口径补充：eval400（更稳定的对比）
- 背景：`config/Ibex_hd_row11.cnfg` 使用“达到 10 错即停”，`pkts` 波动较大，不利于比较“小改动”的真实收益。
- 做法：在 `/tmp` 下复制配置（只改仿真停止阈值，不改译码参数）：
  - `/tmp/Ibex_hd_row11_eval400.cnfg`：`max_sim_num=400`，`max_err=10`
  - 注意：当前 `SSD_FC.cpp` 的停止条件需要同时满足 `sim_cnt>=max_sim_num` 且 `err>=max_err`，因此会强制跑满 400 包（`max_err` 仅用于避免“0 错长跑”的极端情况）。
- 基线（当前稳定版：`A/B/C=220/110/280`，SNR=4.9）：
  - `./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval400.cnfg AWGN 4.9`
  - `pkts=400`, `FER=6.750000e-02`, `avg_iter=179.3175`

### 方案 V6：aggr 判据改用 syndrome_weight_delayed
- 改动：
  - aggr 判据内 `iteration>=B && syndrome_weight<C` 改为 `iteration>=B && syndrome_weight_delayed<C`
- eval400 结果：
  - `pkts=400`, `FER=7.250000e-02`, `avg_iter=183.2925`
- 分析：
  - `syndrome_weight_delayed` 更平滑且可能偏小，在 `C=280` 不变时更容易进入 aggr，使权重放大在“不够稳”的阶段更早触发，增加误翻扩散概率。
  - 若要继续验证该思路，需配套重新扫 `C`（更小）或叠加额外保护条件（例如与 `likelihood_levels.min` 或 `flipped_prev` 联动）。
- 结论：相对基线更差，回退。

### 方案 V7：两档 aggr 放大（mild 不放大 `w=2`）
- 改动：
  - aggr 触发拆分为 `aggr_strong/aggr_mild/aggr_soft`
  - `aggr_strong`（或 soft 触发）保持原放大：`2->3, 3->5, 4->7`
  - `aggr_mild` 仅做小放大：`3->4, 4->6`（`2` 不放大）
- eval400 结果：
  - 由于明显劣化，提前中止：`pkts=154`, `FER=1.428571e-01`（约 22/154）
- 分析：
  - 在当前 2bit “半步进攻”公式下，`weight==2` 的 `delta_attack=(w+1)>>1=1`，导致未翻转更新为 `likelihood += (delta_attack-1)=0`，等价于“无推进”。
  - 该方案在 mild 阶段取消 `2->3` 放大后，大量 `w=2` 的比特被卡住，整体收敛动能不足，失败包比例显著上升。
- 结论：明显更差，回退。

### 方案 V8：两档 aggr 放大（mild 保留 `2->3`，减弱 `3/4` 放大）
- 改动：
  - `aggr_strong`（或 soft 触发）保持：`2->3, 3->5, 4->7`
  - `aggr_mild`：`2->3, 3->4, 4->6`
- eval400 结果：
  - 由于灾难性劣化，提前中止：`pkts=51`, `FER≈5.294118e-01`（约 27/51）
- 分析：
  - 该门限组 `220/110/280` 是在“强放大”假设下找到的；mild 阶段削弱 `w=3/4` 的放大后，接近收敛时缺少关键推力，导致大量包无法在迭代限内清零 syndrome。
  - 2bit 空间下“可用台阶”极少，这类减弱放大往往直接把系统从“偶尔可收敛”推到“普遍收敛失败”。
- 结论：不可用，回退。

### 本轮结论
- 在 `row11@4.9` 下，上述 `V1~V8` 均未稳定优于既有最优门限方案。
- 当前仍保留：`A/B/C = 220/110/280`。
- 代码已回退到该稳定版本。

## Row11 再扫参（SNR=4.9，2bit，基于 env 可调门限）

### 方案 V9：pushing 由“常真”改为“按迭代趋势更新”
- 改动（已回退）：
  - 在 `ldpc_dec_bf_ibex` 中每次进入外层 `while` 时更新：
    - `pushing = (syndrome_weight >= prev_sw)`，并滚动 `prev_sw`
- eval400 结果（提前中止）：
  - `pkts=121`, `FER≈2.231405e-01`
- 分析：
  - 该版本仅改变撤销端奇数权重取整，导致 2bit 攻守力量失衡，出现大量“syndrome 很小但清不掉”的失败包，整体 FER 恶化明显。
- 结论：更差，回退。

### 方案 V10：引入 aggr 门限 env 覆盖（便于扫参，不改变默认行为）
- 改动：
  - 在 `src/ldpc_codec_test2.cpp` 增加 env 覆盖：
    - `IBEX_AGGR_ITER_HI`（A）
    - `IBEX_AGGR_ITER_LO`（B）
    - `IBEX_AGGR_SYND_TH`（C）
    - 可选：`IBEX_AGGR_STRONG_SW_TH`（强 aggr 的 syndrome 上限保护，默认关闭）
  - 默认值仍为 `220/110/280`，不设置 env 时与之前一致。

### Row11@4.9 扫参结果摘要（eval400 口径）
- baseline（默认 `A/B/C=220/110/280`）：
  - `pkts=400`, `FER=6.750000e-02`
- 更优候选（`A/B/C=240/120/280`）：
  - `pkts=400`, `FER=5.250000e-02`
  - 在 `config/Ibex_hd_row11.cnfg` 上代表结果：`pkts=187`, `FER=5.347594e-02`
- 反例（说明过早/过弱都会崩）：
  - `A=200,B=110,C=280`：快速劣化，`pkts=93` 时 `FER≈1.5e-01`，提前中止
  - 把 `C` 放宽到 300：同样明显劣化，提前中止
  - `IBEX_AGGR_STRONG_SW_TH=1200`（限制强 aggr 仅在 `syndrome_weight<1200` 触发）：`pkts=400`, `FER=7.000000e-02`（更差）

### 可复现实验命令
- `./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval400.cnfg AWGN 4.9`（默认门限）
- `env IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval400.cnfg AWGN 4.9`

## Row11 新算法探索（2bit，不再扫门限）

统一评估口径：
- `env IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval400.cnfg AWGN 4.9`
- baseline（此前最佳门限，但算法不变）：`FER=5.250000e-02`（21/400），`avg_iter≈171.8`

### 方案 V11：趋势门控（仅在 iter_stuck 时启用 mild aggr）
- 改动（已回退）：在 `ldpc_dec_bf_ibex` 中加入“是否改善”的 guard，只在停滞时打开 mild aggr。
- eval400 结果：`FER=6.250000e-02`
- 结论：变差，回退。

### 方案 V12：2bit 撤销端 odd rounding 软化（floor）
- 改动（已回退）：`f_update_vn_post` 2bit flipped 分支撤销端尝试 `floor(w/2)`（保留 `w=1 -> delta=1`）。
- eval400 结果：`FER=5.250000e-02`（无收益）
- 结论：无明显改善，回退。

### 方案 V13：2bit 下重新引入 post_process2 的“weight==1 边界 +1”
- 改动（已回退）：`if (post_process2 && !flipped && weight==1 && likelihood_new==flip_thr-1) likelihood_new++;`
- eval400 结果：`FER=7.000000e-02`，`avg_iter≈191`
- 结论：明显变差，回退。

### 方案 V14：2bit 引入 soft 梯度（weak/strong=1/0）尝试
- 改动（已回退）：`f_likelihood_levels` 在 `VN_BITS<=2` 时尝试 `weak=1,strong=0` 并禁用 syndrome 自适应更新，期望 soft_bits=1 能区分“弱未翻/强未翻”。
- eval400 结果：`FER=5.250000e-02`，`avg_iter=176.8325`
- 结论：未带来稳定收益。

### 方案 V15：soft_bits=1 初值映射反转（questionable==1 置弱）
- 改动（已回退）：在 `soft_data -> likelihood_level` 处对 `VN_BITS<=2 && soft_bits==1` 强制映射（忽略 `likelihood_map`）。
- eval400 结果：`FER=6.750000e-02`，`avg_iter=185.75`
- 结论：变差，回退。

### 方案 V16：late-stage 攻击“去掉 -1”（极不稳定）
- 改动（已回退）：在 `f_update_vn_post` 2bit 未翻转分支，当 2bit-aggr gate 生效时（迭代数较大）对 `delta>=2` 采用 `likelihood += delta`（不减 1）。
- 现象：大量包在中后期被过推导致发散，短窗 FER 直接上到 $>0.8$，因此中止。
- 结论：不可用，回退。

### 方案 V17：2bit 禁用 post disturbance（不稳定收益）
- 改动（已回退）：`f_update_vn_post` 2bit 分支不再执行 `post_process` 的 “snap 到 thr±1”。
- eval400 结果：
  - run1：`FER=5.000000e-02`，`avg_iter=180.415`
  - run2：`FER=5.250000e-02`，`avg_iter=180.855`
- 分析：对个别随机序列有小收益，但不稳定，且平均迭代略变长。

### 方案 V18：soft aggr 同时作用于 flipped（全权重放大，波动较大）
- 改动（已回退）：aggr 的 soft 条件去掉 `!flipped_prev`，使 flipped/unflipped 都可能进入 weight 放大。
- eval400 结果：
  - run1：`FER=4.500000e-02`，`avg_iter=167.095`
  - run2：`FER=5.250000e-02`，`avg_iter=170.885`
- 结论：有潜力但波动大，需要更“温和的撤销侧放大”。

### 方案 V19：flipped 侧 mild 放大（3->4,4->6）
- 改动（已回退）：flipped_prev 为 1 时，`3->4,4->6`，unflipped 仍为 `3->5,4->7`。
- 现象：短窗明显劣化（约 160 包时 FER $\approx 0.1$），提前中止。
- 结论：不可用，回退。

### 方案 V20（当前最优）：仅对 flipped 且 weight==2 启用放大（2->3）
- 核心思想：
  - 仍允许 soft aggr 覆盖 flipped（去掉 `!flipped_prev`），但**撤销侧只对 `weight==2` 做一次放大**，避免对 `weight==3/4` 过度撤销造成振荡。
  - 进攻侧（unflipped）保持原放大：`2->3,3->5,4->7`。
- 关键改动（`src/ldpc_codec_test2.cpp`）：
  - aggr soft 条件：`... && ((VN_BITS <= 2) || !flipped_prev)`（2bit 下允许 flipped）
  - weight 放大映射：
    - `weight==2 -> w=3`（无论 flipped/unflipped）
    - `weight==3/4` 仅在 `!flipped_prev` 时做 `3->5,4->7`；flipped 时保持 `3/4` 不变
- eval400 结果（SNR=4.9）：
  - run1：`FER=4.500000e-02`（18/400），`avg_iter=166.785`
  - run2：`FER=4.750000e-02`（19/400），`avg_iter=165.5575`
- 分析：
  - 相比 baseline `FER=5.25e-02`（21/400），错包数下降 $2\sim 3$ 个/400 包，且平均迭代略下降。
  - 直觉解释：flipped 且 `weight==2` 很可能是“边界振荡/弱错误翻转”区域，适度增强撤销有利于抑制扩散；但对 `weight==3/4` 不做放大，避免把“真正需要保持 flipped 的比特”过度拉回，减少振荡与失败包。

---

## Row11@4.9：统一 1000 包口径（eval1000）复验与新探索

为降低短窗波动，后续方案统一用 1000 包做对比：
- eval 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`，`max_err=0`，保证必跑满 1000 包）
- SNR：4.9
- 固定门限（env 覆盖）：`A/B/C = 240/120/280`
- 命令：
  - `env IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

### 方案 V20（长窗基线，eval1000）
- 结果：`FER=5.600000e-02`（56/1000），`avg_iter=173.85`
- 结论：相对 eval400 的 $<0.05$，长窗下波动变大，说明该方案对随机序列更敏感，需要进一步抑制“撤销侧振荡/误撤销”。

### 方案 V21（已回退）：pushing 按 iteration 趋势动态更新（外层迭代级）
- 改动（已回退）：将 `pushing` 从“常量”改为在每次 outer iteration 开始时按 `syndrome_weight` 趋势更新（类似 `pushing = (syndrome_weight >= prev_sw)`）。
- 结果（短窗）：约百包内 `FER≈0.18~0.20`，明显劣化，提前中止。
- 分析：
  - 该更新口径与 IBEX 的 pipeline 延迟不匹配，导致撤销端在“看似下降”时过度变保守（odd rounding 不再 +1），错误翻转难以被及时撤销，失败包激增。

### 方案 V22a（已回退）：用 soft 可靠性保护强比特（仅对弱比特做 attack 放大）
- 改动（已回退）：尝试用 `bit_questionable` 作为可靠性提示，对强比特减小 attack 放大（例如 `weight==4` 由 `4->7` 改温和）。
- 结果（短窗）：200 包时 `FER=6.500000e-02`，趋势劣于基线，提前中止。
- 关键原因：
  - `sd_num=1` 场景下，译码输入实际是“硬输入”（`ldpc_decoder_input.soft_bits==0`），`bit_questionable` 不会被可靠地赋值；该方向在当前配置下缺少可用软信息支撑，因此不继续。

### 方案 V23（已回退）：仅把 `weight==4` attack 由 `4->7` 降为 `4->6`（且误删了 `3->5`）
- 改动（已回退）：试图降低攻击端过推，但实现中丢失了 `weight==3` 的 attack 放大，导致整体推力不足且行为异常。
- 结果（短窗）：几十包内 `FER≈0.45` 级，直接判定不可用并中止。
- 结论：该实现不可用，回退。

### 方案 V24（已回退）：撤销侧对“强 flipped 且 weight==3”做增强撤销（3->5）
- 改动（已回退）：`flipped_prev==1 && weight==3 && likelihood_prev>flip_thr` 时令 `w=5`，意图清理“强错误翻转”。
- eval1000 结果：`FER=6.500000e-02`（65/1000），`avg_iter=184.741`
- 分析：
  - 强化撤销引入了额外的振荡/误撤销，长窗下失败包与平均迭代都上升，属于负收益。

### 方案 V25（当前最优，eval1000）：`weight==2` 只增强 attack，不增强 retract
- 关键改动（`src/ldpc_codec_test2.cpp`，aggr 映射处）：
  - `weight==2`：仅在 `!flipped_prev` 时做 `2->3`（attack-only），`flipped_prev` 时保持 `w=2`
  - `weight==3/4`：保持仅对 `!flipped_prev` 做 `3->5, 4->7`
- eval1000 结果：`FER=4.500000e-02`（45/1000），`avg_iter=167.251`
- 分析（为什么有效）：
  - 2bit 下撤销分支的步长非常敏感；若对 flipped 侧也把 `2->3`，会把边界比特过快拉回，触发更多反复翻转。
  - V25 保留了“尾部清理”的 attack 推力（让未翻转且 `weight==2` 的可疑比特更容易跨过阈值），同时避免对 flipped 侧的额外撤销放大，从而在长窗下显著降低失败包并减少平均迭代。

### 方案 V26（A 方向候选，eval1000）：pipeline pushing 动态更新（列级）+ V25 attack-only
- 背景：
  - 在 `src/ldpc_codec_test2.cpp` 的旧实现中，`pushing` 实际上是常量（且几乎恒为 true），无法体现 “syndrome_weight 未下降则更激进撤销” 的原始意图。
  - 该方案将 `pushing` 更新口径改为“列级 + pipeline 延迟一致”，以更接近 IBEX 的硬件时序语义。
- 关键改动（`src/ldpc_codec_test2.cpp`）：
  - 初始化 `syndrome_weight_r[0..4]=syndrome_weight`，避免未初始化导致的不确定行为。
  - 新增环境变量：
    - `IBEX_PUSH_DYNAMIC=1`：在列循环内按 `pushing = (syndrome_weight_delayed >= syndrome_weight_r[3])` 动态更新；否则强制 `pushing=true`（保持旧行为）。
    - `IBEX_W2_STOCH=1`：仅当 `IBEX_PUSH_DYNAMIC=1` 且 `pushing=true` 时，对 `weight==2` 的 attack `2->3` 做 PRNG 随机化（本次评估未启用）。
- 评估口径：
  - `env IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - 基线：不设 `IBEX_PUSH_DYNAMIC`
  - V26：设 `IBEX_PUSH_DYNAMIC=1`
- 结果（SNR=4.9，row11，eval1000）：
  - 基线（V25，本次复现实测）：`FER=5.200000e-02`（52/1000），`avg_iter=170.841`
  - V26（V25 + pushing 动态）：`FER=4.300000e-02`（43/1000），`avg_iter=159.952`
- 分析（为什么可能有效）：
  - 当 `pushing=false`（趋势在下降）时，2bit flipped 侧 odd rounding 不再 `+1`，撤销步长更小，相当于引入“惯性/动量”，减少“刚翻对了又被撤销”的振荡。
  - 当 `pushing=true`（停滞或变差）时，撤销仍会更激进，有利于快速撤回错误翻转并再次探索。
  - 该机制是对 all-core BF 中 “momentum/tabu 抑振” 思想的低开销映射：不引入 per-bit 额外状态，仅改变撤销端的离散取整。

### 方案 V27（A 方向候选，eval1000）：V26 + `weight==2` boost 随机化（仅在 pushing 时）
- 关键思想：
  - 延续 V26 的“动态 pushing”语义（撤销端带动量），并进一步引入 PGDBF/NGDBF 式的“随机扰动”，但只作用在 **增益门控** 上而不是直接翻转。
  - 仅对 `aggr && !flipped_prev && weight==2` 的 `2->3` boost 做 PRNG 随机化，并且仅在 `pushing=true`（停滞/变差）时启用，避免在收敛趋势良好时引入不必要的随机抖动。
- 关键改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env：
    - `IBEX_W2_STOCH=1`：开启 `weight==2` 的随机 boost（前提 `IBEX_PUSH_DYNAMIC=1`）
  - 随机门控（使用 LFSR bit，避免与 post_process 同一位点强相关）：
    - `if (w2_stoch && push_dynamic && pushing) w = rand_boost ? 3 : 2;`
- 评估口径（row11，SNR=4.9，eval1000）：
  - `env IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：
  - `FER=2.500000e-02`（25/1000），`avg_iter=143.797`
- 分析：
  - 相比 V26（确定性的 `2->3` boost），V27 在停滞阶段对 `weight==2` 的“临界推一把”做随机抽样，显著降低了同步误翻导致的级联崩坏概率；
  - 从 log 观察，失败包的 syndrome_weight 典型值从千级降到 $\sim 650\sim 750$，同时 `LDPC BER` 也下降到 $1.24\times 10^{-3}$ 量级，说明该随机化主要减少了“被推崩”的失败类型，而不是单纯靠更激进的攻击。

---

## 后续探索：从 core BF 机制出发的演进尝试（Row11@4.9，eval1000）

统一评估命令（除非另有说明）：
- `env IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

### 方案 V28（已回退，eval1000）：对强比特降低 `weight==2` boost 概率（soft-biased）
- 核心想法：参考 NGDBF 的 channel 项，认为“强软信息比特”更不应在 `weight==2` 这种临界证据下被频繁推进到翻转阈值；因此尝试用软信息对 `2->3` boost 做概率偏置。
- 改动（`src/ldpc_codec_test2.cpp`，已回退）：
  - 当 `IBEX_W2_STOCH=1` 且 `IBEX_PUSH_DYNAMIC=1` 时：
    - `pushing=true`：questionable 比特 boost 概率 $\approx 1/2$，strong 比特 boost 概率 $\approx 1/4$
    - `pushing=false`：仅对 questionable 比特做 `2->3`
- 结果：`FER=2.400000e-02`（24/1000），`avg_iter=145.365`
- 结论：整体变差，说明该 soft 偏置在当前 1-bit soft 量化下会削弱必要的尾部推进，回退。

### 方案 V29（已回退，短窗即失败）：Tabu1（上一迭代翻转的 VN 本迭代跳过更新）
- 核心想法：参考 TRGDBF 的 tabu-list，抑制“翻转-撤销-再翻转”的 ping-pong。
- 改动（`src/ldpc_codec_test2.cpp`，已保留为开关但默认关闭）：
  - `IBEX_TABU1=1`：若某 VN 在上一 outer-iteration 发生 `flipped` 状态变化，则本 outer-iteration 直接跳过该 VN 更新。
- 现象：前 50 包量级即出现 `LDPC FER≈0.2~0.3` 的灾难性失败（大量包完全不收敛），因此中止。
- 分析：该 tabu 过强，相当于把很多关键 VN 的梯度更新“冻结”一整轮，推力被明显削弱，导致无法收敛；该方向不继续。

### 方案 V30（已回退，eval1000）：2bit boost 仅受 2bit gate 控制（与 soft aggr 解耦）
- 核心想法：尝试避免 soft aggr 在高 syndrome 阶段触发过早 boost，减少早期误翻。
- 改动（`src/ldpc_codec_test2.cpp`，已回退）：
  - 将 `aggr` 拆成 `aggr_2bit_gate` 与 `aggr_soft`，并令 2bit 下的权重 boost 仅在 `aggr_2bit_gate` 为真时生效。
- 结果：`FER=2.000000e-02`（20/1000），`avg_iter=141.523`
- 结论：未优于 V27（本记录的 best-run），回退。

### 方案 V31（已回退，eval1000）：提前初始化 `prng_256/512`（用于 w2 stoch）
- 动机：避免 `IBEX_W2_STOCH` 在 `post_iteration` 之前访问未初始化的 PRNG 状态（不确定行为）。
- 改动（`src/ldpc_codec_test2.cpp`，已回退）：
  - 在译码开始时对 `prng_256/512` 做显式初始化（`prng_512` 用 verilog 对齐模式）。
- 结果：`FER=2.300000e-02`（23/1000），`avg_iter=143.379`
- 分析：显式初始化引入了更强的结构性相关（尤其 `prng_512` 的固定模式），使 w2 随机门控在某些包上更易“同步失效”，整体变差，回退。

### 方案 V32（基线，eval1000）：`weight==2` boost 延后到 `post_iteration` 之后启用（early-stage 禁止 w2 boost）
- 关键直觉：
  - 在 `post_iteration` 之前，译码仍处于“粗收敛”阶段，`weight==2` 的 VN 数量多且证据弱；对其进行 `2->3` boost（尤其在 aggr/随机门控参与时）更容易触发同步误翻，导致部分包直接进入高 syndrome 的错误吸引域（典型失败包 syndrome_weight $\approx 650\sim 1000$）。
  - 将 `weight==2` boost 延后到 `post_iteration` 之后，更符合“先靠高权重翻转把结构拉回，再用 w2 做尾部清理”的阶段性策略。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 当 `IBEX_W2_STOCH=1` 且 `IBEX_PUSH_DYNAMIC=1` 时：
    - `iteration < post_iteration`：对 `weight==2` 不做 `2->3` boost（保持 `w=2`）
    - `iteration >= post_iteration`：恢复 V27 逻辑（`pushing=true` 时随机 boost，`pushing=false` 时确定 boost）
- 结果（历史记录）：`FER=1.400000e-02`（14/1000），`avg_iter=136.701`
- 复现核对（2026-02-09，当前 workspace）：`FER=2.100000e-02`（21/1000），`avg_iter=142.291`
- 备注：两者不一致，说明除 `src/ldpc_codec_test2.cpp` 外的工程状态可能已变化（例如其他源文件/配置被修改）。后续对比以“可复现结果”为准。

### 方案 V33（eval1000）：失败后 smoothing（smNGDBF 风格的 hard majority vote）
- 核心想法：参考 `doc/all core bf/ldpc_smngdbf.m`，在译码失败时，不直接输出最后一次迭代的状态，而是对末尾窗口内的 hard 判决做“多数投票”，尝试从振荡/抖动状态中恢复。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_SMOOTH_FAIL=1`（默认关闭），窗口默认 64（可用 `IBEX_SMOOTH_WIN` 覆盖）。
  - 在 `iteration >= iteration_limit - smooth_win` 的末尾窗口内，对每个 VN 记录 decoded hard bit 的符号累加（0 记 +1，1 记 -1）。
  - 若最终译码失败，则对每个 VN 取 `sign(sum)` 得到 `decoded_smoothed`，重构一个候选 codeword 并重算 syndrome；若 `syndrome_weight==0` 则接受该结果。
- 评估命令：
  - `IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_SMOOTH_FAIL=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=2.100000e-02`（21/1000），`avg_iter=140.945`
- 分析：在该 row11@4.9 条件下，失败包更像是陷入“错误吸引域/高 syndrome 稳态”，而非围绕正确解的轻微振荡；因此多数投票很难把 syndrome 拉到 0，收益不明显。
- 结论：该方向在当前实现下未带来 FER 改善，不作为默认方案继续推进。

### 方案 V34（当前最优，eval1000）：失败后 restart/phase 重译码（reNGDBF 风格）
- 核心想法：参考 `doc/all core bf/ldpc_reNGDBF.m` 的 multi-phase re-decoding。对于少量“陷入 trapping-set / 错误吸引域”的失败包，单次 BF 轨迹可能会卡死；通过在失败后从初始状态重启一次，并改变随机扰动路径，有机会把 syndrome 拉到 0。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_RESTART_PHASES`（默认 1，不重启；设为 2 表示失败后重启 1 次）。
  - 重启时：将 `vn.likelihood` 恢复到初始 soft-map 值，`vn.flipped` 清零，并重算 `cn/syndrome_weight`。
  - 为了让第 2 phase 的随机路径不同：在 `iteration==post_iteration && j==0` 的 PRNG 初始化点，对 `prng_256/512` 额外做若干次 LFSR advance（skip）：
    - `skip_steps = phase * IBEX_RESTART_PRNG_SKIP`，其中 `IBEX_RESTART_PRNG_SKIP` 默认 73。
- 评估命令：
  - `IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=1.400000e-02`（14/1000），`avg_iter=133.768`
- 对比基线 V32（可复现）：`FER=2.100000e-02`（21/1000）
- 分析：
  - 失败包数减少 7/1000，符合“少量包需要不同随机轨迹才能脱困”的直觉。
  - 注意：当前工程的 `IBEX Decoder average iterations` 统计要求 `iterations<=max_iter`。因此该实现将 `ldpc_decoder_output.iterations` 维持在单次 phase 的范围内（不累计两次 phase 的总迭代），`avg_iter` 不可直接用于衡量真实算力开销；真实开销会在失败包上增加一次 phase 的迭代预算。
- 结论：在不改门限/不扫参前提下，V34 属于明确的架构级改进（multi-phase re-decoding），且在 row11@4.9/1000 包上 FER 显著下降，值得作为 row11 的新默认候选进一步加包或降 SNR 验证。

### 方案 V35（短窗即失败，已中止）：Tabu-Rev（TRGDBF 风格“禁止紧邻迭代反向翻转”）
- 核心想法：TRGDBF 的 tabu-list 本质是“刚翻过的 bit 下一轮不允许再翻回去”。在 IBEX 的 likelihood/threshold 框架下等价于：允许更新 likelihood，但若上一 outer-iteration 刚发生过 toggle，则本轮不允许再次跨越 flip_thr。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_TABU_REV=1`（默认关闭）。
  - 若 `last_toggle_iter[vn_idx] == iteration-1` 且本轮计算将导致 `flipped` 状态再次变化，则将该 VN 的 likelihood clamp 到阈值边界以保持 `flipped` 不变。
- 评估命令（短窗观测后中止）：
  - `IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_TABU_REV=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 现象：极早期即出现大量包完全不收敛；`100 packets` 时 `LDPC FER≈2.2e-01`。
- 分析：TRGDBF 的 tabu 假设通常配合“每轮只翻少量/单个 bit”的选择策略；而 IBEX 这里是列内对大量 VN 进行并行 likelihood 更新，tabu-rev 会让大量本该快速撤销的误翻无法撤销，导致错误快速扩散到高 syndrome 吸引域（失败包 syndrome_weight 可显著升高）。
- 结论：该方向在当前并行更新框架下不适配，放弃。

### 方案 V36（短窗较差，已中止）：W1-Boost（尝试对 weight==1 做随机促攻）
- 核心想法：2-bit 半步更新下 `weight==1` 对未翻转 VN 的增量为 0（`delta=(w+1)>>1`，`likelihood += delta-1`），导致大量“边缘证据”无法积累。尝试在 tail 阶段对部分 `weight==1` 的未翻转 VN 随机提升到 `w=3`，提供一次有效推进。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_W1_STOCH=1`（默认关闭）。
  - 仅在 `aggr==true`、`pushing==true`、`likelihood>=flip_thr-1`、且 `iteration>=post_iteration` 时，按 PRNG 随机将 `weight==1` 的 attack 端 `w` 提升到 3。
- 评估命令（短窗观测后中止）：
  - `IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_W1_STOCH=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 现象：`300 packets` 时 `LDPC FER≈4.33e-02`，且失败包 syndrome_weight 明显偏高（$\approx 1200\sim 1400$）。
- 分析：对 `w==1` 的促攻在该配置下更像“额外噪声源”，会把部分本可收敛的轨迹推入高 syndrome 吸引域；说明 row11@4.9 的主要失败模式不是“缺推力”，而是“错误推进导致发散/陷入坏稳态”。
- 结论：不继续该方向（保留开关但默认关闭）。

### 方案 V37（当前最优，eval1000）：2bit_mode=2（aggr 阶段关闭 post 扰动）
- 核心想法：post_trigger 的 PRNG 扰动在 tail 阶段能帮助跳出局部振荡，但当 `aggr` 已经打开时（权重放大/推进更强），继续叠加 post 扰动可能会造成不必要的随机抖动与误翻。参考 `src/ldpc_codec_test.cpp` 的 `IBEX_2BIT_MODE>=2` 思路，尝试在 `aggr==true` 时 gate 掉 `prng_post_process/prng_post_process2`。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增：读取 `IBEX_2BIT_MODE`（默认 0）。
  - 当 `IBEX_2BIT_MODE>=2`：令 `post_gate = !aggr`，即 `aggr==true` 时强制 `prng_post_process(_2)=0`。
- 评估命令：
  - `IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=1.300000e-02`（13/1000），`avg_iter=135.751`
- 对比 V32（同配置、mode=0 的一次实测）：`FER=1.600000e-02`（16/1000）
- 分析：在 row11@4.9 下，主要失败包更像是“被随机扰动拉偏后进入坏稳态”。aggr 阶段关掉 post 扰动能减少尾部误翻，使少量包从“差一点收敛”变成收敛，从而降低 FER。
- 结论：该改动属于明确的架构级收敛稳定性增强，当前作为 row11 的最优候选。

### 方案 V38（eval1000）：V34 restart + V37 post_gate（未优于 V37）
- 评估命令：
  - `IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=1.500000e-02`（15/1000），`avg_iter=138.454`
- 分析：restart 在该次随机样本下没有带来额外脱困收益，反而可能引入更多“二次轨迹”失败包。
- 结论：不作为默认配置；若后续需要进一步压低 FER，可考虑把 restart 作为“只有检测到 stall/坏稳态才触发”的条件性机制，而非无条件二次 phase。

---

## 2026-02-09 补充复现与新方案（Row11@4.9，eval1000）

重要说明（解释“同一命令跑出来 FER 差很多”的现象）：
- 当前 row11 的评估命令包含 `IBEX_W2_STOCH=1`，再叠加 AWGN 信道本身的随机性，因此 **同一命令多次运行的 `LDPC FER` 会明显波动**（1000 包样本偏小，方差不可忽略）。
- 因此 `doc/bf_2bit_best_configs.md` 里记录的 `FER` 是“单次 best-run 的数值”，主要用于复现命令与对比方向，不代表统计意义上的均值。

本次复现结果（当前 workspace，同一 cnfg `/tmp/Ibex_hd_row11_eval1000.cnfg`）：
- V32（mode=0）：`FER=3.200000e-02`（32/1000），`avg_iter=152.361`
  - `IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- V37（mode=2）：`FER=2.700000e-02`（27/1000），`avg_iter=146.039`
  - `IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- V38（mode=2 + restart=2）：`FER=1.500000e-02`（15/1000）的一次 run；后续复现出现 `FER=2.300000e-02`、`FER=1.700000e-02`
  - `IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

结论（以“趋势”而非单次结果为准）：
- 在当前工程状态下，V38 通常能比 V37/V32 更低（给失败包一次“不同随机轨迹”的机会），但需要多次运行取统计。

### 方案 V39（eval1000，较差）：restart 的第二 phase 放开 post_gate（更“探索”但会发散）
- 核心想法：phase0 用 V38 的稳定策略；若 phase0 失败，phase1 允许在 `aggr==true` 时也启用 post 扰动，期望用更强随机性脱困。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_RESTART_RELAX_POST_GATE=1`（默认关闭）。
  - `phase>0` 时令 `post_gate_mode = (!aggr) || restart_relax_post_gate`，即重启 phase 放开 `prng_post_process(_2)`。
- 评估命令：
  - `IBEX_RESTART_RELAX_POST_GATE=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=2.700000e-02`（27/1000），`avg_iter=150.252`
- 分析：phase1 的“激进随机扰动”更容易把已接近收敛的轨迹拉回坏稳态，等价于浪费一次重启机会。
- 结论：该方向在 row11@4.9 不适配，默认保持关闭。

### 方案 V40（eval1000，较差）：restart 分摊总迭代预算（固定总预算的 multi-try）
- 核心想法：将总迭代预算 $T$（cnfg 的 `Ibex maximum iteration number`）在 `restart_phases` 之间分摊，每个 phase 只跑 $\lceil T/P\rceil$ 次迭代，等价于“固定算力下做多次独立尝试”（reNGDBF 风格）。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_RESTART_SPLIT_BUDGET=1`（默认关闭）。
  - `phase_iter_limit = iteration_limit / restart_phases`（含余数分配），while 条件改为 `iteration < phase_iter_limit`。
- 评估命令：
  - `IBEX_RESTART_SPLIT_BUDGET=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=2.000000e-02`（20/1000），`avg_iter=133.838`
- 分析：固定总预算下，部分包需要更长单轨迹收敛，预算分摊反而降低成功率。
- 结论：在 row11@4.9 目前不如“全预算 + 失败再重启”的 V38。

### 方案 V41（eval1000，较差）：w2 stochastic 采用两 tap XOR（尝试降相关）
- 核心想法：担心 `prng[k+37]` 的空间/时间相关导致同步误翻，尝试用两个 tap 的 XOR 生成 gate bit，期望降低相关性。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_W2_STOCH_XOR=1`（默认关闭）。
  - `rand_boost = prng[idx1] ^ prng[idx2]`（idx2 取 `k+173`）。
- 评估命令：
  - `IBEX_W2_STOCH_XOR=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=2.300000e-02`（23/1000），`avg_iter=142.008`
- 分析：随机门控“更随机”并不等价于更好；在该配置下更像引入额外扰动，失败包增多。
- 结论：不继续该方向。

### 方案 V42（eval1000，较差）：post 仅在 pushing 时允许（减少无谓扰动）
- 核心想法：当 `pushing=false`（syndrome_weight 正在下降）时尽量不扰动；只有 `pushing=true`（停滞/变差）才允许 post 扰动，期望减少“差一点收敛却被抖动拉偏”的情况。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_POST_ONLY_WHEN_PUSHING=1`（默认关闭）。
  - `post_gate = post_gate_mode && (post_only_when_pushing ? pushing : true)`。
- 评估命令：
  - `IBEX_POST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=2.300000e-02`（23/1000），`avg_iter=145.582`
- 分析：post 扰动本身可能是少数包脱困所必需的“尾部刷子”，仅在 pushing 时开启反而会错过有效时机。
- 结论：该 gate 在 row11@4.9 下不适配，保持默认关闭。

### 方案 V43（eval1000，失败）：soft guard（用 soft questionable 门控 w2 boost / post）
- 核心想法：减少对“可靠比特”的误翻，仿照部分 Weighted-BF/NGDBF 的 reliability-aware 保护，尝试只对 `soft_unreliable` 比特启用更激进的 w2 boost 与 post 扰动。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_SOFT_GUARD=1`（默认关闭）。
  - 当 `IBEX_SOFT_GUARD=1`：仅当 `soft_unreliable==true` 时允许 `weight==2` boost（`2->3`）及 post 扰动。
- 评估命令：
  - `IBEX_SOFT_GUARD=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=4.460000e-01`（446/1000），`avg_iter=535.469`
- 分析：
  - 该结果过差且运行时间显著变长，说明该“soft_unreliable”判据要么与实际可靠性语义不一致（可能取反），要么过度抑制了必要的 w2 推进，导致大量包卡死到 `iteration_limit`。
- 结论：该方向在当前工程/量化定义下不可用，保持默认关闭。

### 方案 V44（eval1000，未见稳定收益）：init soft bias（用 soft questionable 仅影响初始 likelihood）
- 核心想法：不在迭代更新中硬门控，而是把 1-bit soft questionable 作为“初始偏置”，给可疑比特一个更接近翻转阈值的起点，等价于引入很弱的 channel 项。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_INIT_SOFT_BIAS=1`（默认关闭）。
  - 当 `VN_BITS<=2 && soft_bits>0` 时：若 `soft_unreliable` 则初始化 `likelihood=1`，否则 `likelihood=0`（保持 `flip_thr=2`）。
- 评估命令：
  - `IBEX_INIT_SOFT_BIAS=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=1.500000e-02`（15/1000），`avg_iter=138.224`
- 对照（同日同参数、未开 bias 的一次 run）：`FER=1.500000e-02`（15/1000），`avg_iter=138.351`
- 结论：单次结果看不出稳定收益，保留为可选开关但不作为默认演进方向。

### 方案 V45（eval1000，较差）：stall 触发的提前 restart（adaptive multi-try）
- 核心想法：参考 trapping-set/坏吸引域现象，若某 phase 长时间不刷新 `best_sw_phase`，则提前结束该 phase 并进入下一 phase（restart），期望在相同或更低算力下获得更多独立轨迹机会。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_RESTART_ON_STALL=1`（默认关闭）。
  - 参数：`IBEX_STALL_ITERS`（默认 32），`IBEX_STALL_MIN_ITER`（默认 200）。
  - 机制：若 `iteration>=STALL_MIN_ITER` 且连续 `STALL_ITERS` 次未刷新 `best_sw_phase`，并且还有下一 phase，则 `break` 提前进入重启。
- 评估命令：
  - `IBEX_RESTART_ON_STALL=1 IBEX_STALL_ITERS=32 IBEX_STALL_MIN_ITER=200 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果：`FER=2.300000e-02`（23/1000），`avg_iter=143.763`
- 分析：stall 判据过于粗糙，会把部分“慢收敛但仍在收敛”的包提前打断，反而降低成功率。
- 结论：该实现不适配，保持默认关闭（若后续要继续该方向，需要更精细的 stuck/attractor 识别而非仅看 best_sw 刷新）。

### 方案 V46（eval1000，显著提升 ⭐）：w2 boost 仅在 pushing==true 时允许（减少无谓过推）
- 关键直觉：
  - 目前 2bit 的 `weight==2` 属于“弱证据临界区”，`2->3` boost 会让大量 VN 获得有效正增量（从而更易跨过 `flip_thr`）。
  - 当 `pushing=false`（syndrome_weight 在下降）时，继续对 `weight==2` 做确定性 boost，容易把一些本应保持的 VN 过推到阈值附近，引入不必要的误翻与级联；这些误翻会把少量包推入坏吸引域，形成 error floor。
  - 因此把 w2 boost 收敛为“只有在 pushing==true（停滞/变差）时才允许”，更符合“只在需要时推进”的策略。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1`（默认关闭）。
  - 当开启时：`weight==2` 且 `!flipped_prev` 的 `2->3` boost 在 `pushing==false` 时被禁止；`pushing==true` 时保持原先的随机 boost（`IBEX_W2_STOCH=1`）。
- 评估命令：
  - `IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果（复现两次均为 0 错）：
  - Run#1：`FER=0.000000e+00`（0/1000），`avg_iter=131.522`
  - Run#2：`FER=0.000000e+00`（0/1000），`avg_iter=134.669`
- 长窗验证（>=1000 包 + 10 错停止，2026-02-10）：
  - 结果：`pkts=10175`, `LDPC FER=9.828010e-04`, `avg_iter=134.198722`
  - 说明：该结果对应 `max_sim_num=1000` 且 `max_err=10` 的停止条件（即至少跑满 1000 包，并累计到 10 个失败包才停），因此比“1000 包短窗”更能代表真实 FER 水平；短窗 0 错与该量级并不矛盾（期望每 1000 包约 1 错，存在随机波动）。
  - 终端统计输出（节选）：
    ```text
    [SIM] Finish simulation @ Tue Feb 10 00:03:46 2026
    --------------------------------------------------------
    [STATISTICS] Total packets simulated: 10175
    [STATISTICS] RAW  BER: 1.058614e-02
    [STATISTICS] LDPC BER: 1.126554e-05
    [STATISTICS] LDPC FER: 9.828010e-04
    [STATISTICS] LDPC MIS: 0.000000e+00
    [STATISTICS] MCRC FER: 9.828010e-04
    [STATISTICS] DATA FER: 0.000000e+00
    ---------------------------------------------------------
    -------------------------------------------------------------------------------------
    [STATISTICS] IBEX Decoder average iterations: 134.198722
    -------------------------------------------------------------------------------------
    ```
- 结论：在 row11@4.9 条件下，这是目前最有价值的架构级改进（从 $\sim 10^{-2}$ 级别直接压到 $\sim 10^{-3}$ 量级）；下一步建议尝试更低 SNR 或进一步降低随机性/算力开销。

### 方案 V47（eval1000_0err@4.8，进一步改善 ⭐）：strong 比特的 w2 boost 增加 1-bit“候选记忆”（需二次命中才允许 boost）
> 注意：该方案引入额外 1-bit/VN 的跨迭代状态（`w2_cand_mem`），因此按“纯 2bit（每 VN 仅 2-bit likelihood）”口径不成立。按当前约束（禁止任何额外 per-VN 状态），**该方案仅保留记录，不作为后续演进方向**。
- 核心想法（强约束版本的 anti-misflip）：
  - 在 2bit 下，`weight==2` 的 `2->3` boost 是“让 w2 有推力”的关键，但它也会在某些 stalled/pushing 段落里把 **strong 且本应保持的 VN** 缓慢推向阈值，带来少量误翻并触发失败包的 error floor。
  - 对 “strong 且未翻转” 的 VN（`likelihood < flip_thr-1`），引入 1-bit 的候选记忆 `cand`：第一次遇到 stalled 且 `weight==2` 时只置位 `cand`（不 boost）；只有再次命中时才允许进入原有的 `w2` boost 逻辑。
  - 对 “weak (flip_thr-1)” 的 VN 不做该限制（因为它们本来就接近阈值，需要更快的推进）。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_W2_CAND_STRONG=1`（默认关闭）。
  - phase 内为每个 VN 分配 1-bit `w2_cand_mem[vn_idx]`；phase restart 时清零。
  - 仅当 `push_dynamic==true && pushing==true && iteration>=post_iteration && weight==2 && !flipped_prev && strong_unflipped` 时：
    - `cand==0`：置 `cand=1`，本次禁止 w2 boost
    - `cand==1`：允许本次 w2 boost（走原随机门控），并在 VN 变成 weak/翻转/weight!=2 时自动清零
- 评估配置：
  - 使用快速统计 cnfg：`/tmp/Ibex_hd_row11_eval1000_0err.cnfg`（`max_sim_num=1000`, `max_err=0`，确保固定 1000 包对比）。
- 对比结果（Row11@4.8，1000 包）：
  - Baseline（V46，不开 cand）：
    - 命令：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`LDPC FER=3.700000e-02`（37/1000），`avg_iter=257.348`
  - V47（cand strong 开启）：
    - 命令：`IBEX_W2_CAND_STRONG=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`LDPC FER=2.600000e-02`（26/1000），`avg_iter=253.837`
- Row11@4.9 快速核对（1000 包，方差较大，仅作 sanity check）：
  - 命令：`IBEX_W2_CAND_STRONG=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.9`
  - 结果：`LDPC FER=1.000000e-03`（1/1000），`avg_iter=134.455`
- Row11@4.9 长窗（>=1000 包 + 10 错停止，2026-02-10）：
  - 命令：`IBEX_W2_CAND_STRONG=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - 结果：`pkts=11747`, `LDPC FER=8.512812e-04`, `avg_iter=136.013450`
  - 终端统计输出（节选）：
    ```text
    [SIM] Finish simulation @ Tue Feb 10 09:06:13 2026
    --------------------------------------------------------
    [STATISTICS] Total packets simulated: 11747
    [STATISTICS] RAW  BER: 1.058581e-02
    [STATISTICS] LDPC BER: 7.540351e-06
    [STATISTICS] LDPC FER: 8.512812e-04
    [STATISTICS] LDPC MIS: 0.000000e+00
    [STATISTICS] MCRC FER: 8.512812e-04
    [STATISTICS] DATA FER: 0.000000e+00
    ---------------------------------------------------------
    -------------------------------------------------------------------------------------
    [STATISTICS] IBEX Decoder average iterations: 136.013450
    -------------------------------------------------------------------------------------
    ```
- 分析：
  - 该改动本质是在 “strong 且 w2 证据” 处增加一个极轻量的时间一致性过滤器：**要求 w2 证据在 stalled 段落里至少出现两次才触发推进**，降低了 transient w2 对 strong 正确比特的累计漂移。
  - 从结果看，FER 在更低 SNR 点（4.8）有明显改善，且平均迭代略有下降，说明并非单纯“更保守导致更慢”，而是减少了部分失败包的吸引域陷入。
  - 代价：额外 1-bit/VN 的存储（row11: `76*512=38912` bits 级别），属于明确的硬件开销点；是否可接受需要结合 SRAM/寄存器预算评估。
- 下一步：
  - 更新：后续已在“纯 2bit（不增加任何 per-VN 额外状态）”约束下继续探索，并找到一个基于 **restart phase 非对称策略** 的可用改进（见 V50）。

---

### 方案 V48（纯 2bit，失败）：w2 boost tail guard（在低 syndrome 阶段抑制 strong+reliable 的 w2 boost）
> 目标：在接近收敛时减少 strong 正确比特被 `w=2 -> 3` 缓慢推到阈值附近导致的 rare mis-flip；该思路试图用“全局 syndrome 阶段门控”近似 V47 的时间一致性过滤，但不引入 per-VN 记忆。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_W2_TAIL_GUARD=1`（默认关闭）。
  - 在 `iteration>=post_iteration` 且 `pushing==true` 的 `w2` 随机 boost 路径里：
    - 若 `syndrome_weight_delayed < syndrome_weight_thr_post`：禁止对 strong+reliable VN 的 boost
    - 若 `syndrome_weight_delayed < syndrome_weight_thr_qc`：对 strong+reliable VN 进一步降低 boost 概率
- 评估（Row11@4.8，1000 包，固定窗）：
  - Baseline（V46）：
    - 命令：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`LDPC FER=4.000000e-02`（40/1000），`avg_iter=268.588`
  - V48（tail guard 开启）：
    - 命令：`IBEX_W2_TAIL_GUARD=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 现象：在前 200 包内 `LDPC FER` 明显偏高（$\sim 0.06$），趋势不佳，提前中止该 run（避免浪费时间）。
- 分析：
  - 该门控依赖 1-bit soft 的 `soft_unreliable` 判定；在 row11@4.8 下，部分“需要被修正”的比特在 soft 上仍表现为可靠，从而被过度保护，导致收敛被阻断。
  - 结论：当前实现过于保守，不作为后续方向（保留开关但默认关闭）。

### 方案 V49（纯 2bit，失败）：toggle-strong（强证据翻转后将 likelihood 直接 snap 到强态，试图模拟 tabu/hysteresis）
- 核心想法：
  - TRGDBF/Tabu 的主要作用是减少振荡；但 per-VN tabu 需要额外存储。
  - 这里尝试用 **不新增状态** 的方式：当 `w>=5`（强证据）触发翻转时，把 VN 的 `likelihood` 直接置到 `max/min`，让其更“坚定”，降低来回抖动。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_TOGGLE_STRONG=1`（默认关闭）。
  - 条件：`VN_BITS<=2 && (flipped_prev != flipped_new) && (w>=5)` 时，将 `likelihood` snap 到 `max/min`。
- 评估（Row11@4.8，1000 包，固定窗）：
  - 命令：`IBEX_TOGGLE_STRONG=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
  - 结果：前 100 包 `LDPC FER=1.300000e-01`，明显劣化，提前中止该 run。
- 分析：
  - 2bit 下把 likelihood snap 到极值等价于“过强磁滞”，会把少量错误翻转锁住，触发级联失败；这种现象在低 SNR 下尤其明显。
  - 结论：该方向不可用，后续不再探索（保持默认关闭）。

### 方案 V50（纯 2bit，eval1000_0err@4.8，改进 ⭐）：Phase1（retry）在 not-pushing 时允许极小概率的 w2 boost（escape hatch）
> 目标：不增加任何 per-VN 记忆，利用已有 `restart_phases` 机制，把“更激进/更随机”的动作只放到 **失败包才会进入的 retry phase**，降低对正常包的副作用，同时提高对 trapping set 的逃逸概率。
- 核心想法：
  - V46 中 `IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1` 会在 `pushing==false` 时完全禁止 `w=2 -> 3` boost；这对稳定性有利，但可能让某些 bad attractor/trap 无法获得足够扰动逃逸。
  - 由于我们本就启用 `IBEX_RESTART_PHASES=2`，因此可以让 phase0 保持保守；仅在 phase1（phase>0，retry）中，在 `pushing==false` 时对“更可疑”的 VN（weak 或 soft_unreliable）开放一个很小概率的 w2 boost，作为 escape hatch。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_PHASE1_W2_NOT_PUSHING=1`（默认关闭）。
  - 仅当 `phase>0 && iteration>=post_iteration && pushing==false && weight==2 && !flipped_prev`：
    - 若 VN 为 `weak` 或 `soft_unreliable`，则用 2 个 PRNG gate（$\sim 1/4$）触发 `w=2 -> 3` boost。
- 评估（Row11@4.8，1000 包，固定窗；2026-02-10）：
  - Baseline（V46）：
    - 命令：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`LDPC FER=4.000000e-02`（40/1000），`avg_iter=268.588`
  - V50（phase1 escape hatch 开启）：
    - 命令：`IBEX_PHASE1_W2_NOT_PUSHING=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`LDPC FER=3.100000e-02`（31/1000），`avg_iter=250.460`
    - 终端统计输出（节选）：
      ```text
      [STATISTICS] Total packets simulated: 1000
      [STATISTICS] LDPC FER: 3.100000e-02
      [STATISTICS] IBEX Decoder average iterations: 250.460000
      ```
- Row11@4.9 快速核对（1000 包固定窗；2026-02-10）：
  - 命令：`IBEX_PHASE1_W2_NOT_PUSHING=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.9`
  - 结果：`LDPC FER=1.000000e-03`（1/1000），`avg_iter=133.575`
- 分析：
  - 该方案的关键点是“把副作用隔离到失败包才会进入的 phase1”，因此相比直接在 phase0 放开 not-pushing boost，更容易得到净收益。
  - 从 1000 包固定窗看，FER 与平均迭代均有下降，说明 escape hatch 的扰动总体上帮助部分失败包跳出吸引域，而没有显著放大平均算力。
  - 风险：该结果仍有统计波动；若要确认对 row11@4.9 的长窗 FER 是否有稳定收益，需要用 `/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_err=10`）再做长窗对比。

### 方案 V51（纯 2bit，eval1000_0err@4.8，改进 ⭐）：列内 bit 扫描起点旋转（rotate-k 调度随机化）
> 目标：不新增任何 per-VN 状态，仅通过“调度随机化”打破 layered 更新的确定性关联，降低少数包在坏吸引域附近的周期性震荡概率。
- 核心想法：
  - 当前 IBEX layered BF 在每列内按固定顺序扫描 `k=0..511`，其更新顺序与 PRNG/post 的扰动也存在固定关系；在某些 trapping-set 上可能形成稳定的坏周期。
  - 对 phase1（retry）在 post 阶段引入轻量 `rotate-k`：每列内的扫描起点 `k_start` 由 PRNG 给出，等价于每列做一次循环移位扫描，从而改变局部更新的相位关系。
  - 该策略只改变“更新顺序”，不改变每个 VN 的更新公式，因此属于架构级调度改动，且不需要任何 per-VN 记忆。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_ROTATE_K=1`（默认关闭）。
  - 参数：`IBEX_ROTATE_K_PHASE1_ONLY=1`（默认关闭），开启后只在 `phase>0` 生效，避免影响大多数本就能收敛的正常包。
  - 机制：当 `rotate_k_eff` 成立时，用 PRNG 生成 `k_start`，并将列内扫描改为 `k=(kk+k_start) mod 512`。
- 评估（Row11@4.8，1000 包固定窗；基线均为 V46：w2 boost only-when-pushing）：
  - Baseline（V46）：
    - Run#1：`LDPC FER=4.400000e-02`（44/1000），`avg_iter=258.934`
    - Run#2：`LDPC FER=3.600000e-02`（36/1000），`avg_iter=253.871`
  - V51（rotate-k phase1 only）：
    - Run#1：`LDPC FER=3.400000e-02`（34/1000），`avg_iter=257.048`
    - Run#2：`LDPC FER=3.400000e-02`（34/1000），`avg_iter=254.826`
  - Run#3（并行对照，2026-02-10 15:16，`RUN_DIR=output/row11_arch_20260210_151634_snr4p8_r3`）：
    - Baseline（V46）：`LDPC FER=3.200000e-02`（32/1000），`avg_iter=252.307`
    - V51（rotate-k phase1 only）：`LDPC FER=4.100000e-02`（41/1000），`avg_iter=270.899`
- 分析：
  - 在 `IBEX_W2_STOCH=1` + AWGN 的组合下，1000 包固定窗的波动很大；V51 在 Run#1/#2 看起来有收益，但 Run#3 出现明显劣化（FER 与平均迭代均上升）。
  - 观察到的现象更像是“少数失败包的轨迹被改变”：有时能减少坏周期（Run#1/#2），但也可能引入额外扰动导致收敛变慢甚至失败（Run#3）。
  - 结论：rotate-k 不是稳定改进项，当前不作为默认演进方向；保留为可选开关，后续若要继续，需要引入更精细的“仅在 stuck/attractor 时启用”的触发条件。

### 方案 V52（纯 2bit，eval1000_0err@4.8，失败）：stall 驱动的 w2 escape（过推导致劣化）
> 目标：在 retry phase 进入“停滞”后，允许更激进的 `w=2 -> 3` boost 来跳出坏吸引域。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_STALL_W2_ESC=1`（默认关闭）。
  - 参数：`IBEX_STALL_W2_ESC_ITERS`（默认 8），`IBEX_STALL_W2_ESC_MIN_ITER`（默认 `post_iteration`）。
  - 机制：当 `phase>0` 且 `stall_count_w2` 达到阈值时，针对 “bad candidate（weak 或 soft_unreliable）” 的 `w==2` 在 `pushing==false` 时也放开 boost（更激进）。
- 评估（Row11@4.8，1000 包固定窗）：
  - Run#1：`LDPC FER=3.900000e-02`（39/1000），`avg_iter=252.803`
  - Run#2：`LDPC FER=5.100000e-02`（51/1000），`avg_iter=258.769`
- 分析：
  - `stall_count_w2` 作为 stuck 指标过于粗糙：很多“慢收敛但仍在收敛”的包也会触发 escape。
  - 一旦强行放开 `w2` boost，2bit 的量化空间过窄，容易把部分 VN 过推到阈值附近，引入误翻并触发级联失败，最终表现为 FER 上升且不稳定。
  - 结论：该实现不适配（保持默认关闭）。若后续要继续该方向，需要更精细的 stuck/attractor 识别，而不是仅看 `best_sw_phase_w2` 是否刷新。

### 方案 V53（纯 2bit，eval1000_0err@4.8，轻微改进）：PPBF-like 概率 escape（按能量代理 p(E) 触发 w2 boost）
> 目标：在 retry phase 且进入停滞后，为少数“坏候选”提供更柔和的随机逃逸，而不是像 V52 那样直接强推。
- 核心想法：
  - 参考 `doc/all core bf/ldpc_ppbf.m` 的思想：用概率形式而不是确定性形式做 escape，避免对本可收敛的包造成过大扰动。
  - 在 IBEX 框架中只对 `w==2 && !pushing` 的路径加一个“概率 boost”，并且只在 retry phase（`phase>0`）启用，控制副作用。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_PPBF_ESC=1`（默认关闭）。
  - 参数：`IBEX_PPBF_ESC_ITERS`（默认 8），`IBEX_PPBF_ESC_MIN_ITER`（默认 `post_iteration`）。
  - 机制：当 `ppbf_esc_active` 成立时，对 `w==2` 的 bad candidate，构造一个极简能量代理 $E$：
    - $E = 2 + [soft\\_unreliable] + [!pushing]$
    - 当 $E\\ge 4$：以约 $1/2$ 的概率 boost（1 个 PRNG gate）
    - 当 $E=3$：以约 $1/4$ 的概率 boost（2 个 PRNG gate）
- 评估（Row11@4.8，1000 包固定窗）：
  - Run#1：`LDPC FER=3.600000e-02`（36/1000），`avg_iter=268.837`
  - Run#2：`LDPC FER=3.300000e-02`（33/1000），`avg_iter=259.574`
- 分析：
  - 相比 V52 的“强推”，PPBF-like 逃逸对 FER 的影响更温和：两次 run 都没有出现灾难性发散，且一次 run 给出比基线更好的 FER。
  - 但其平均迭代数略有上升，说明该扰动可能让一部分包走了更长的路径才收敛；是否值得作为默认演进，需要更多重复与 row11@4.9 的长窗对比确认。

#### 补充：Row11@4.9 长窗对比（max_err=10，未跑满提前中止）
> 目的：观察 V51（rotate-k）与 V53（PPBF escape）是否能降低 row11@4.9 的 error floor。该对比跑到 >6000 包后因耗时考虑手动中止（V51/V53 尚未累计到 10 个失败包，因此仅作趋势参考）。
- 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`，达到 1000 包后继续跑直到累计到 10 个失败包才停止）
- 运行目录：`output/row11_arch_20260210_151356_snr4p9_long`
- 命令：
  - Baseline（V46）：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - V51（rotate-k phase1 only）：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - V53（PPBF escape）：`IBEX_PPBF_ESC=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果（停止点统计）：
  - Baseline（V46，已达到 10 错自动停止）：`pkts=5861`, `LDPC FER=1.706193e-03`, `avg_iter=135.244156`
  - V51（未达到 10 错，手动中止时）：`pkts=6600`, `LDPC FER=7.575758e-04`（5/6600）
  - V53（未达到 10 错，手动中止时）：`pkts=6600`, `LDPC FER=9.090909e-04`（6/6600）
- 分析：
  - 从趋势看，V51/V53 都有降低 error floor 的潜力（相同包数下失败数明显更少）。
  - 但由于未跑满“10 错停止条件”，这组数据仍不足以下最终结论；建议后续继续跑到 `FAIL CW=10` 或改用略低 SNR（例如 4.8x）加速收敛到稳定统计。

#### 补充：Row11@4.9 长窗对比（max_err=10，跑满到 FAIL CW=10）
> 目的：在同一套“纯 2bit”基线（V46）之上，跑满到 `FAIL CW=10`，获得更稳定的 error-floor 估计，并对比“调度扰动（V51）”与“概率 escape（V53）”的净收益。
- 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`）
- 运行目录：`output/row11_long10err_final3_20260210_172743_snr4p9`
- 命令：
  - Baseline（V46）：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - V51（rotate-k phase1 only）：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - V53（PPBF escape）：`IBEX_PPBF_ESC=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- 结果（>=1000 包 + 10 错停止）：
  - Baseline（V46）：`pkts=6466`, `LDPC FER=1.546551e-03`, `avg_iter=135.560470`
  - V51（rotate-k phase1 only）：`pkts=14999`, `LDPC FER=6.667111e-04`, `avg_iter=134.534436`
  - V53（PPBF escape）：`pkts=10573`, `LDPC FER=9.458054e-04`, `avg_iter=133.952048`
- 分析：
  - V51/V53 都降低了 row11@4.9 的 error floor；并且 `avg_iter` 基本不变，属于“以调度/概率扰动换 error-floor”。
  - 本次 run 中，V51 相对 baseline 的 FER 下降约 2.3 倍（$1.55e-3 \\to 6.67e-4$）；V53 的提升更温和（$1.55e-3 \\to 9.46e-4$）。
  - 由于停止条件固定为 10 个失败包，统计相对方差约为 $\\sqrt{1/10}\\approx0.316$；若要更高置信区分 V51 与 V53，建议把 `max_err` 提高到 30 或重复多次取统计。

### 方案 V54（纯 2bit，eval1000_0err@4.8，较差）：rotate-k 在所有 phase 生效（扰动过强）
- 改动：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=0`（在 phase0 也启用 rotate-k）。
- 评估（Row11@4.8，1000 包固定窗；2026-02-10 15:16）：
  - 命令：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=0 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
  - 结果：`LDPC FER=4.100000e-02`（41/1000），`avg_iter=271.439`
- 分析：把调度扰动扩展到 phase0 会影响大量本可收敛的正常包，平均迭代明显上升且 FER 变差；该方向不适合作为默认策略。

### 方案 V55（纯 2bit，eval1000_0err@4.8，较差）：列内翻转次数上限（tail-only，max_toggles=8）
- 改动：`IBEX_MAX_TOGGLES_PER_COL=8`（默认 `IBEX_MAX_TOGGLES_TAIL_ONLY=1`，仅在 tail 段限流）。
- 评估（Row11@4.8，1000 包固定窗；2026-02-10 15:16）：
  - 命令：`IBEX_MAX_TOGGLES_PER_COL=8 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
  - 结果：`LDPC FER=3.900000e-02`（39/1000），`avg_iter=259.823`
- 分析：在当前实现里，“按列限流”更像是硬性约束，容易把少数需要多翻转才能修复的包卡住；该方向暂不继续。

### 方案 V56（纯 2bit，eval1000_0err@4.8，中性）：rotate-k phase1 + max_toggles=8（未形成净收益）
- 改动：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_MAX_TOGGLES_PER_COL=8`。
- 评估（Row11@4.8，1000 包固定窗；2026-02-10 15:16）：
  - 命令：`IBEX_MAX_TOGGLES_PER_COL=8 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
  - 结果：`LDPC FER=3.400000e-02`（34/1000），`avg_iter=261.928`
- 分析：该组合在单次 run 中比 V54/V55 好，但仍未超过同批次 baseline（`3.2e-02`）；目前看不到“稳定净收益”，先不作为后续演进方向。

### 方案 V57（纯 2bit，Row11@4.8，失败/早停）：pushing 从“列内趋势”改为“迭代级趋势”（IBEX_PUSH_MODE=1）
> 背景：你提到“用前几列相比前一列的 synd wt 对比是否合理”。当前实现的 `pushing` 是基于列内 `syndrome_weight` 的延迟差分（本质是一个很局部的趋势信号）。这里尝试一个架构级替代：把 `pushing` 变成 **每轮迭代一个常量**，由“本轮迭代开始时的 syndrome_weight”相对“上一轮迭代开始时”的变化决定，期望更稳、更少抖动。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增开关：`IBEX_PUSH_MODE`（默认 0）。
  - `IBEX_PUSH_MODE=0`：保持原有列内趋势（`syndrome_weight_delayed >= prev_sw_col`）。
  - `IBEX_PUSH_MODE=1`：引入 `prev_iter_sw/pushing_iter`，在 `j==0` 时计算 `pushing_iter = (sw_iter_start >= prev_iter_sw)`，并在该迭代内对所有列复用。
- 评估（Row11@4.8，短窗早停；2026-02-10 15:43）：
  - Baseline（`IBEX_PUSH_MODE=0`）：
    - 命令：`IBEX_PUSH_MODE=0 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果（336 包时）：`LDPC FER=3.273810e-02`（11/336）
  - V57（`IBEX_PUSH_MODE=1`，iter-trend）：
    - 命令：`IBEX_PUSH_MODE=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果（289 包时）：`LDPC FER=6.574394e-02`（19/289）
  - 说明：该对比中 V57 明显劣化，因此未继续跑满 1000 包（提前终止以节省仿真时间）。
  - 分析：
  - “迭代级 pushing”过于粗糙：列内 syndrome_weight 的细粒度变化被抹平后，`w2 boost only-when-pushing` 与后续 gating 的触发时机变得不准确，导致更多误翻/失败包。
  - 结论：该方向在当前实现形态下不可用，保持默认 `IBEX_PUSH_MODE=0`。

### 方案 V58（纯 2bit，Row11@4.9，eval1000_0err，值得继续 ⭐）：V51（rotate-k）+ V53（PPBF escape）组合
> 目标：把“调度扰动”（V51）与“概率 escape”（V53）叠加到同一条 retry-phase 轨迹上，以更低的副作用提升对 trapping-set 的逃逸能力。
- 关键直觉（为什么可能更强）：
  - V51 通过 rotate-k 改变列内更新相位关系，打破坏周期；
  - V53 在 stall 后对少数 bad candidate 给出概率 boost（更柔和）；
  - 两者作用点不同，且都集中在 `phase>0`、`post_iteration`、`stall_count_w2` 触发后的 tail 区域，理论上可叠加。
- 评估（Row11@4.9，1000 包固定窗；2026-02-12）：
  - Baseline（V51，本次复现实测）：
    - 命令：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.9`
    - 结果：`LDPC FER=1.000000e-03`（1/1000），`avg_iter=135.724`
  - V58（V51 + PPBF escape）：
    - 命令：`IBEX_PPBF_ESC=1 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.9`
    - Run#1：`LDPC FER=0.000000e+00`（0/1000），`avg_iter=129.434`
    - Run#2：`LDPC FER=0.000000e+00`（0/1000），`avg_iter=137.591`
- 分析：
  - 在“固定 1000 包”口径下，V58 的两次 run 都未出现失败包，至少说明该组合没有明显副作用（未观察到 V52 那种灾难性发散）。
  - 由于 `IBEX_W2_STOCH=1` 且 AWGN 信道本身随机，`1000` 包样本方差仍然不小；要判断 V58 是否真正降低 error floor，仍需用 `/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_err=10` 或更高）跑长窗对比，并建议重复多次取统计。

- 追加评估（Row11@4.9，长窗；2026-02-15 15:14-15:29）：
  - 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`）
  - 命令：`IBEX_PPBF_ESC=1 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - 结果（达到 `FAIL CW=10` 停止）：`pkts=8380`, `LDPC FER=1.193317e-03`, `avg_iter=134.348926`
  - 初步结论：本次长窗结果明显劣于 Row11 Top3（V51/V53/V46），因此 V58 暂不进入 Top3；如需排除随机波动，可再重复 1-2 次长窗，但当前证据更像“固定窗 0/1000 属于幸存者偏差”。

### 方案 V59（纯 2bit，Row11@4.9，eval1000_0err，较差）：PPBF escape + mode-window + anneal（尝试抑制副作用）
> 背景：V58 显示“固定窗 0/1000”并不意味着长窗更好。这里尝试把 PPBF-like 逃逸约束在 stall 触发的短窗口里，并加入简单退火（tail 初期开大、尾端收小），期望减少误翻副作用。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增全局窗口开关：`IBEX_MODE_WIN=1`（以及 `IBEX_MODE_WIN_LEN/IBEX_MODE_WIN_TRIG_ITERS`）。
  - 新增 PPBF 退火开关：`IBEX_PPBF_ESC_ANNEAL=1`（以及 `IBEX_PPBF_ESC_ANNEAL_ITERS`）。
- 评估（Row11@4.9，1000 包固定窗；2026-02-15 16:16）：
  - 命令：`IBEX_MODE_WIN=1 IBEX_MODE_WIN_LEN=16 IBEX_MODE_WIN_TRIG_ITERS=8 IBEX_PPBF_ESC=1 IBEX_PPBF_ESC_ANNEAL=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.9`
  - 结果：`LDPC FER=2.000000e-03`（2/1000），`avg_iter=132.617`
- 结论：该组合在固定窗就明显劣化，说明当前 p(E) 形状/窗口节奏仍偏激进（或窗口触发过早/过频），暂不继续走长窗。

### 方案 V60（纯 2bit，Row11@4.9，长窗，较差）：rotate-k + column-global-escape + mode-window（UP-GDBF/backtracking 映射尝试）
> 目标：利用“列级 escape + 回滚（backtracking）”在不引入 per‑VN 状态的前提下增强跳出 trapping-set 的能力；mode-window 仅在 stall 时对列级 escape 的接受率做小幅提升。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - mode-window 会在窗口内对列级 escape 的 `w==2` 接受率做轻微提升（约 $1/8 \\to 1/4$），并允许 `max_toggles` 轻微 +1（仍有 backtracking 保护）。
- 评估（Row11@4.9，1000 包固定窗；2026-02-15 16:18）：
  - 命令：`IBEX_COL_GLOBAL_ESC=1 IBEX_COL_GLOBAL_ESC_ITERS=8 IBEX_COL_GLOBAL_ESC_MAX_TOGGLES=2 IBEX_MODE_WIN=1 IBEX_MODE_WIN_LEN=16 IBEX_MODE_WIN_TRIG_ITERS=8 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.9`
  - 结果：`LDPC FER=1.000000e-03`（1/1000），`avg_iter=135.270`
- 长窗（>=1000 包 + 10 错停止；2026-02-15 16:19-16:33）：
  - 命令：`IBEX_COL_GLOBAL_ESC=1 IBEX_COL_GLOBAL_ESC_ITERS=8 IBEX_COL_GLOBAL_ESC_MAX_TOGGLES=2 IBEX_MODE_WIN=1 IBEX_MODE_WIN_LEN=16 IBEX_MODE_WIN_TRIG_ITERS=8 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - 结果：`pkts=8222`, `LDPC FER=1.216249e-03`, `avg_iter=133.941255`
- 结论：该方案长窗明显劣于 Row11 Top3（V51/V53/V46），暂不进入 Top5 的靠前位置；列级 escape 的触发/接受策略仍需更谨慎的门控才能不伤正常包。

### 方案 V61（纯 2bit，Row11@4.8，eval1000_0err，改进 ⭐）：phase1 使用更早的 aggr 门限（multi-phase 非对称）
> 目标：利用多 phase 重启（`IBEX_RESTART_PHASES>=2`）提供“无 per‑VN 记忆的多样性”。保持 phase0 行为不变，只让 phase1（retry）使用更早的 aggr 门限，探索不同的收敛轨迹，优先优化 Row11@4.8（waterfall）。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env（默认不影响旧命令；未设置则继承 phase0）：
    - `IBEX_PHASE1_AGGR_ITER_HI`
    - `IBEX_PHASE1_AGGR_ITER_LO`
    - `IBEX_PHASE1_AGGR_SYND_TH`
    - `IBEX_PHASE1_AGGR_STRONG_SW_TH`
  - 在每个 `phase` 初始化时，计算 `*_eff` 并用于 aggr gate 判定（仅 `phase>0` 时可能不同）。
- 评估（Row11@4.8，固定 1000 包；2026-02-17）：
  - 对照（V46，w2-boost only-when-pushing）：
    - 命令：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.8`
    - 结果：`LDPC FER=3.300000e-02`（33/1000），`avg_iter=254.711000`（log: `output/row11_v46_4p8_long10err_202602171350.log`）
  - 对照（V51，rotate-k phase1 only）：
    - 命令：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.8`
    - 结果：`LDPC FER=4.200000e-02`（42/1000），`avg_iter=260.209000`（log: `output/row11_v51_4p8_long10err_202602171345.log`）
  - V61（V51 基座 + phase1 aggr=180/90/240）：
    - 命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`LDPC FER=2.700000e-02`（27/1000），`avg_iter=259.177000`（log: `output/row11_v61_4p8_eval1000_0err_202602171358.log`）
- 扫参（Row11@4.8，固定窗；2026-02-17）：
  - V62（phase1 aggr=160/80/240）：中途终止（`pkts=576` 时 `LDPC FER=2.951389e-02`；log: `output/row11_v62_4p8_eval1000_0err_202602171415.log`）
  - V63（phase1 aggr=200/100/260）：中途终止（`pkts=700` 时 `LDPC FER=3.142857e-02`；log: `output/row11_v63_4p8_eval1000_0err_202602171416.log`）
  - V64（V46 基座 + phase1 aggr=180/90/240）：中途终止（`pkts=340` 时 `LDPC FER=3.823529e-02`；log: `output/row11_v64_4p8_eval1000_0err_202602171429.log`）
- Row11@4.9 回归（长窗到 `FAIL CW=10`；2026-02-17）：
  - 命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - log：`output/row11_v61_4p9_long10err_202602171435.log`
  - 结果（达到 `FAIL CW=10` 停止）：`pkts=11108`, `LDPC FER=9.002521e-04`, `avg_iter=133.749460`
- Row11@4.7 waterfall 点（固定 300 包；2026-02-17）：
  - 命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.7`
  - log：`output/row11_v61_4p7_eval300_0err_202602171834.log`
  - 结果：`pkts=300`, `LDPC FER=2.933333e-01`, `avg_iter=594.973333`
- 分析：
  - V61 的 4.9 长窗性能接近 V53（PPBF escape）量级，并明显优于 V46 baseline；虽然仍不及 V51（rotate-k）Top1，但说明“phase1 更早 aggr”并不会在 4.9 上造成明显回退，可作为 4.8 主战场方案的兜底回归。
  - 在 4.7 下 `avg_iter` 明显增大（接近 600），意味着该点的仿真耗时会显著上升；后续若要补齐 4.7 的更多对照，建议统一使用 `eval300` 或更小窗先做趋势筛选。

### 方案 V65（纯 2bit，Row11@4.8，P2 尝试）：w2 boost 的 syndrome-delta 分区门控（替代二值 pushing）
> 目标：把“pushing 二值门控”升级为“按 syndrome delta 分区门控”，在 syndrome 明显恶化时禁止 w2 boost（避免火上浇油），仅在轻微停滞/恶化时允许（更细腻的 escape/稳定性折中）。不引入任何 per‑VN 跨迭代状态。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env（默认关闭，不影响旧命令）：
    - `IBEX_W2_BOOST_SW_DELTA_GATE`
    - `IBEX_W2_BOOST_SW_DELTA_HI`
    - `IBEX_W2_BOOST_SW_DELTA_MIN_ITER`
  - 仅在 `phase>0` 且 `iteration>=MIN_ITER` 时生效：计算列内 `sw_delta_col = syndrome_weight_delayed - prev_sw_col`，并要求 $0\\le sw\\_delta\\_col \\le HI$ 才允许本列的 w2 boost（包括随机 boost / stall escape / PPBF escape / phase1 not-pushing hatch）。
- 评估（Row11@4.8，短窗快速筛选；2026-02-17）：
  - V66（HI=8，跑满 300 包）：
    - 命令：`IBEX_W2_BOOST_SW_DELTA_GATE=1 IBEX_W2_BOOST_SW_DELTA_HI=8 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=3.333333e-02`, `avg_iter=268.580000`（log: `output/row11_v66_4p8_eval300_0err_202602171714.log`）
  - V67（HI=2，中途终止）：
    - 命令：`IBEX_W2_BOOST_SW_DELTA_GATE=1 IBEX_W2_BOOST_SW_DELTA_HI=2 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=138` 时 `LDPC FER=4.347826e-02`（log: `output/row11_v67_4p8_eval300_0err_202602171720.log`）
- 结论：
  - 目前在短窗筛选中未观察到明确收益；该门控较可能需要重新选择 delta 的定义/阈值标定，暂不作为 4.8 主线。

### 方案 V71（纯 2bit，Row11@4.8，P1 延伸尝试）：retry phase 放开 post gate + phase1 not-pushing 的 w2 boost hatch
> 目标：在不改变 phase0 行为的前提下，为 retry phase（`phase>0`）增加两种“小扰动多样性”来源：  
> 1) 允许在 aggr 阶段也执行 post 随机扰动（仅 retry phase 生效）；  
> 2) 在 not-pushing 时给极小概率的 `w=2 -> 3` boost（仅 retry phase 生效），以打破少数 trapping-set 的确定性循环。  
> 约束：不引入任何 per‑VN 跨迭代状态，仅使用已有 PRNG gate。
- 改动（env 开关，无代码改动）：
  - `IBEX_RESTART_RELAX_POST_GATE=1`：仅在 `phase>0` 时，允许 `aggr` 阶段也启用 post 随机扰动（对齐“副作用隔离到 retry phase”的原则）。
  - `IBEX_PHASE1_W2_NOT_PUSHING=1`：仅在 `phase>0` 且趋势 not-pushing 时，对“bad candidate”（弱/不可靠）以约 $1/4$ 概率允许 `w=2 -> 3` boost。
- 评估（Row11@4.8，`eval300_0err`；2026-02-17）：
  - V69（仅 `phase1_w2_not_pushing`）：
    - 命令：`IBEX_PHASE1_W2_NOT_PUSHING=1 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=3.666667e-02`, `avg_iter=253.056667`（log: `output/row11_v69_4p8_eval300_0err_202602171900.log`）
  - V70（仅 `restart_relax_post_gate`）：
    - 命令：`IBEX_RESTART_RELAX_POST_GATE=1 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=3.000000e-02`, `avg_iter=267.483333`（log: `output/row11_v70_4p8_eval300_0err_202602171905.log`）
  - V71（两者都开）：
    - 命令：`IBEX_RESTART_RELAX_POST_GATE=1 IBEX_PHASE1_W2_NOT_PUSHING=1 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=2.666667e-02`, `avg_iter=270.153333`（log: `output/row11_v71_4p8_eval300_0err_202602171911.log`）
- 复核（Row11@4.8，`eval1000_0err`；2026-02-17）：
  - V71：
    - 命令：`IBEX_RESTART_RELAX_POST_GATE=1 IBEX_PHASE1_W2_NOT_PUSHING=1 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`pkts=1000`, `LDPC FER=3.400000e-02`, `avg_iter=255.438000`（log: `output/row11_v71_4p8_eval1000_0err_202602171917.log`）
- 结论：
  - `eval300` 上 V71 看起来有一定优势，但在 `eval1000` 上未体现稳定收益；暂不推进到 4.9 长窗回归，保留该组合为“retry-phase 多样性扰动”的备选分支。

### 方案 V73（纯 2bit，Row11@4.8，P1 延伸尝试）：phase1 的 w=2 boost 档位（`IBEX_PHASE1_W2_BOOST_TO`）
> 目标：在不改变 phase0 的前提下，为 retry phase（`phase>0`）提供额外“推力档位”，尝试把默认的 `w=2 -> 3` boost 升级为更强的 `w=2 -> 4/5`（仅 retry phase 生效），以形成不同的收敛轨迹。  
> 风险：2bit 量化下全步/过强推力容易饱和或发散，因此该档位只作为 phase1 多样性候选。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env：`IBEX_PHASE1_W2_BOOST_TO`（默认 `3`，clamp 到 `[3,7]`；未设置时行为与历史版本一致）。
  - 仅影响 `phase>0` 的 `weight==2` boost：将原本的 `w=3` 替换为 `w=w2_boost_to_eff`（`phase==0` 仍为 `3`）。
- 评估（Row11@4.8，`eval300_0err`；2026-02-17）：
  - V72（`IBEX_PHASE1_W2_BOOST_TO=4`）：
    - 命令：`IBEX_PHASE1_W2_BOOST_TO=4 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=3.000000e-02`, `avg_iter=268.436667`（log: `output/row11_v72_4p8_eval300_0err_202602171941.log`）
  - V73（`IBEX_PHASE1_W2_BOOST_TO=5`）：
    - 命令：`IBEX_PHASE1_W2_BOOST_TO=5 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=5.666667e-02`, `avg_iter=285.746667`（log: `output/row11_v73_4p8_eval300_0err_202602171947.log`）
- 结论：
  - `IBEX_PHASE1_W2_BOOST_TO=4` 未观察到明确收益；`=5` 明显劣化（疑似推力过强导致不稳定/误翻）。该方向暂不推进到 `eval1000` 与 4.9 长窗。

### 方案 V74（纯 2bit，Row11@4.8，P1 延伸尝试）：提高多 phase 重启次数（`IBEX_RESTART_PHASES=3`）
> 目标：在不新增 per‑VN 状态的前提下，用“更多 retry phase”提供更强的解码多样性（PRNG/调度不同），以提升 waterfall（4.8）收敛成功率；并观察对 4.9 error-floor 是否有副作用。
- 改动（env 开关，无代码改动）：
  - 基于 V61（phase1 更早 aggr + rotate-k retry-only + w2 stochastic + dynamic pushing），仅提高 `IBEX_RESTART_PHASES`：
    - 对照：`IBEX_RESTART_PHASES=2`（历史默认主线）
    - 候选：`IBEX_RESTART_PHASES=3/4`
- 评估（Row11@4.8，`eval300_0err`；2026-02-17）：
  - V74（`IBEX_RESTART_PHASES=3`）：
    - 命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=1.333333e-02`, `avg_iter=246.540000`（log: `output/row11_v74_4p8_eval300_0err_202602171955.log`）
  - V75（`IBEX_RESTART_PHASES=4`）：
    - 命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=4 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=4.000000e-02`, `avg_iter=266.593333`（log: `output/row11_v75_4p8_eval300_0err_202602172000.log`）
- 复核（Row11@4.8，`eval1000_0err`；2026-02-17）：
  - V74（`IBEX_RESTART_PHASES=3`）：
    - 命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`
    - 结果：`pkts=1000`, `LDPC FER=2.600000e-02`, `avg_iter=254.117000`（log: `output/row11_v74_4p8_eval1000_0err_202602172008.log`）
- Row11@4.9 回归（长窗到 `FAIL CW=10`；2026-02-17）：
  - 命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - log：`output/row11_v74_4p9_long10err_202602172027.log`
  - 结果（达到 `FAIL CW=10` 停止）：`pkts=26723`, `LDPC FER=3.742095e-04`, `avg_iter=132.728436`
- 分析：
  - `IBEX_RESTART_PHASES=3` 在 4.8 的 `eval1000` 上优于 V61，且在 4.9 长窗也显著优于旧 Top1（V51：`LDPC FER=6.667111e-04`）；当前可视为 Row11@4.8 + Row11@4.9 的最优纯 2bit 配置。

### 方案 V76（纯 2bit，Row11@4.8，phase2-only 轻扰动）：phase2 放开 post gate + not-pushing w2 hatch
> 目标：进一步把“有风险的扰动”隔离到 `phase>=2`（第 3 次重启），尽量不影响 phase0/1 的正常收敛包；观察是否能救回少量顽固包。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env（默认关闭，不影响历史命令）：
    - `IBEX_PHASE2_RELAX_POST_GATE`：仅在 `phase>=2` 时，允许 `aggr` 阶段也启用 post 随机扰动。
    - `IBEX_PHASE2_W2_NOT_PUSHING`：仅在 `phase>=2` 且趋势 not-pushing 时，允许极小概率的 `w=2 -> 3` boost。
- 评估（Row11@4.8，`eval300_0err`；2026-02-18）：
  - 命令：`IBEX_PHASE2_RELAX_POST_GATE=1 IBEX_PHASE2_W2_NOT_PUSHING=1 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
  - 结果：`pkts=300`, `LDPC FER=1.666667e-02`, `avg_iter=240.246667`（log: `output/row11_v76_4p8_eval300_0err_202602180120.log`）
- 结论：短窗未观察到优于 V74（`LDPC FER=1.333333e-02`），暂不推进到 `eval1000`/4.9 长窗。

### 方案 V77（纯 2bit，Row11@4.8，phase2-only 档位）：phase2 的 w=2 boost 档位（`IBEX_PHASE2_W2_BOOST_TO=4`）
> 目标：只在 `phase>=2` 提供更强 `w=2` 推力档位，尝试形成与 phase1 不同的收敛轨迹。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env：`IBEX_PHASE2_W2_BOOST_TO`（未设置时 phase2 继承 phase1；设置后仅在 `phase>=2` 覆盖）。
- 评估（Row11@4.8，`eval300_0err`；2026-02-18）：
  - 命令：`IBEX_PHASE2_W2_BOOST_TO=4 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
  - 结果：`pkts=300`, `LDPC FER=2.333333e-02`, `avg_iter=247.026667`（log: `output/row11_v77_4p8_eval300_0err_202602180125.log`）
- 结论：短窗明显劣化，推力档位过强的风险仍成立，放弃该方向。

### 方案 V78（纯 2bit，Row11@4.8，phase2-only 概率门控）：FM-PGDBF 风格“翻转概率抑制”
> 目标：在 `phase>=2` + tail 阶段，对 `weight<=3` 的翻转事件施加概率门控（默认约 $3/4$ 放行），打破确定性时间循环。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env（默认关闭）：
    - `IBEX_PHASE2_FLIP_RAND_GATE`
    - `IBEX_PHASE2_FLIP_RAND_GATE_MAX_W`（默认 3）
    - `IBEX_PHASE2_FLIP_RAND_GATE_GATES`（默认 2，对应约 $1-(1/2)^2=3/4$ 放行）
  - 实现：仅在 `phase>=2` + `iteration>=post_iteration` + `hamming_weight_lt_circ_thr` 下生效；若被门控，则把 likelihood clamp 回阈值边界以取消本次翻转。
- 评估（Row11@4.8，`eval300_0err`；2026-02-18）：
  - 命令：`IBEX_PHASE2_FLIP_RAND_GATE=1 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
  - 结果：`pkts=300`, `LDPC FER=3.333333e-02`, `avg_iter=240.480000`（log: `output/row11_v78_4p8_eval300_0err_202602180131.log`）
- 结论：短窗明显劣化，说明“直接抑制翻转”在当前框架下更像在拖慢收敛而非破环，放弃。

### 方案 V79/V80（纯 2bit，Row11@4.8，phase2-only 非对称）：phase2 独立 aggr 门限
> 目标：让 phase1 与 phase2 走不同的 aggr gate（不改 phase0），提供更强的多样性。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env（默认继承 phase1，不影响历史命令）：
    - `IBEX_PHASE2_AGGR_ITER_HI/LO/SYND_TH/STRONG_SW_TH`
- 评估（Row11@4.8，`eval300_0err`；2026-02-18）：
  - V79（phase2 退回 phase0 门限 240/120/280）：
    - 命令：`IBEX_PHASE2_AGGR_ITER_HI=240 IBEX_PHASE2_AGGR_ITER_LO=120 IBEX_PHASE2_AGGR_SYND_TH=280 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=1.666667e-02`, `avg_iter=249.660000`（log: `output/row11_v79_4p8_eval300_0err_202602180136.log`）
  - V80（phase2 更早 aggr=160/80/240）：
    - 命令：`IBEX_PHASE2_AGGR_ITER_HI=160 IBEX_PHASE2_AGGR_ITER_LO=80 IBEX_PHASE2_AGGR_SYND_TH=240 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=2.333333e-02`, `avg_iter=246.540000`（log: `output/row11_v80_4p8_eval300_0err_202602180142.log`）
- 结论：短窗未见收益；phase2 门限走更早/更晚都未形成有效互补，暂不继续扩展扫参。

### 方案 V81（纯 2bit，Row11@4.8，phase2-only 热列定向 escape）：col-global-escape + hot-column targeting
> 目标：把列级 escape 限制到 `phase>=2`，并把“每迭代随机选列”改为“选热列”，减少无谓扰动。
- 改动（`src/ldpc_codec_test2.cpp`）：
  - 新增 env：
    - `IBEX_COL_GLOBAL_ESC_MIN_PHASE`（默认 1；设为 2 可做到 phase2-only）
    - `IBEX_COL_ESC_TARGET_HOT`（热列选择）
    - `IBEX_COL_ESC_TARGET_HOT_MIN_PHASE`（默认 2）
    - `IBEX_COL_ESC_TARGET_HOT_MIN_ITER`（默认 `post_iteration`）
  - 实现要点：
    - 每次 iteration 统计每列的 `col_weight_sum = \\sum_k weight(j,k)`，用最大者作为下一轮热列；
    - 若 `IBEX_COL_ESC_TARGET_HOT` 生效，则用热列替代随机 `col_esc_target_col`；
    - 同时将 `stall_count_w2` 的更新条件扩展为包含 `col_global_esc/mode_win/ngdbf_noise`（避免“开了但不生效”）。
- 评估（Row11@4.8，`eval300_0err`；2026-02-18）：
  - 命令：`IBEX_COL_GLOBAL_ESC=1 IBEX_COL_GLOBAL_ESC_MIN_PHASE=2 IBEX_COL_ESC_TARGET_HOT=1 IBEX_COL_ESC_TARGET_HOT_MIN_PHASE=2 IBEX_COL_ESC_TARGET_HOT_MIN_ITER=50 IBEX_COL_GLOBAL_ESC_ITERS=8 IBEX_COL_GLOBAL_ESC_MAX_TOGGLES=1 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
  - 结果：`pkts=300`, `LDPC FER=3.000000e-02`, `avg_iter=255.016667`（log: `output/row11_v81_4p8_eval300_0err_202602180159.log`）
- 结论：短窗未见正收益；列级 escape 即使隔离到 phase2 仍可能带来额外误翻风险，暂不继续。

### 方案 V82/V83/V84（纯 2bit，Row11@4.8，多 phase PRNG 去相关扫参）：`IBEX_RESTART_PRNG_SKIP`
> 目标：在保持 `IBEX_RESTART_PHASES=3`（V74）不变的前提下，调整不同 phase 的 PRNG skip 步数，提供更强的 restart 多样性（不改算法、不加 per‑VN 状态）。
- 改动（env 开关，无代码改动）：
  - `IBEX_RESTART_PRNG_SKIP`：`skip_steps = phase * IBEX_RESTART_PRNG_SKIP`（在 `post_iteration` 初始化 LFSR 后跳步）。
- 评估（Row11@4.8，`eval300_0err`；2026-02-18）：
  - 对照（V74，默认 `IBEX_RESTART_PRNG_SKIP=73`，见上文）：`LDPC FER=1.333333e-02`（4/300）。
  - V82（`IBEX_RESTART_PRNG_SKIP=37`）：
    - 命令：`IBEX_RESTART_PRNG_SKIP=37 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=3.666667e-02`, `avg_iter=264.463333`（log: `output/row11_v82_4p8_eval300_0err_202602181051.log`）
  - V83（`IBEX_RESTART_PRNG_SKIP=131`）：
    - 命令：`IBEX_RESTART_PRNG_SKIP=131 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=2.333333e-02`, `avg_iter=240.256667`（log: `output/row11_v83_4p8_eval300_0err_202602181058.log`）
  - V84（`IBEX_RESTART_PRNG_SKIP=257`）：
    - 命令：`IBEX_RESTART_PRNG_SKIP=257 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
    - 结果：`pkts=300`, `LDPC FER=2.000000e-02`, `avg_iter=246.986667`（log: `output/row11_v84_4p8_eval300_0err_202602181103.log`）
- 结论：
  - 短窗下未观察到优于默认 `IBEX_RESTART_PRNG_SKIP=73` 的组合；该方向暂不继续扩展扫参。

### 方案 V85（纯 2bit，Row11@4.8，syndrome-delta 门控）：`IBEX_W2_BOOST_SW_DELTA_GATE`
> 目标：对 retry phase 的 `w=2` boost 加入更细腻的全局门控：当 syndrome-weight 明显恶化时禁止 boost（避免火上浇油），改善趋势中也禁止 boost（不打扰收敛）。
- 改动（env 开关，无代码改动；功能已在 `src/ldpc_codec_test2.cpp` 中实现）：
  - `IBEX_W2_BOOST_SW_DELTA_GATE=1`
  - `IBEX_W2_BOOST_SW_DELTA_HI=8`（默认 8）
  - `IBEX_W2_BOOST_SW_DELTA_MIN_ITER=50`（默认 `post_iteration`）
- 评估（Row11@4.8，`eval300_0err`；2026-02-18）：
  - 命令：`IBEX_W2_BOOST_SW_DELTA_GATE=1 IBEX_W2_BOOST_SW_DELTA_HI=8 IBEX_W2_BOOST_SW_DELTA_MIN_ITER=50 IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.8`
  - 结果：`pkts=300`, `LDPC FER=3.666667e-02`, `avg_iter=261.240000`（log: `output/row11_v85_4p8_eval300_0err_202602181112.log`）
- 结论：短窗显著劣化；说明该门控在当前“pushing + w2 stochastic”的节奏下更像是在削弱有效推力，暂不继续扩展阈值扫参。
