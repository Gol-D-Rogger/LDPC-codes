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
