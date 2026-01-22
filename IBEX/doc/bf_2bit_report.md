# IBEX BF 2-bit 方案评估记录（ssd_fc_test, VN_BITS=2, AWGN 5.4）

命令：`./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`，默认 100 包或 10 错停止，FER 目标接近 0。

## 方案 E：全权重进攻/撤销（aggr 关闭）
- 修改：`f_update_vn_post` 使用 `delta=weight`，2bit 关闭 aggr；其他保持原步长公式。
- 结果：FER=1.0（前 35 包全失败，syndrome_weight 千级）。
- 分析：全步长 + 无半步/趋势在 2bit 下过于激进，初始权重高导致发散。

## 方案 B：攻守分离（未翻转全步，已翻转半步；aggr 关闭）
- 修改：`f_update_vn_post` 2bit 下未翻转用全步，已翻转用半步；aggr 仅 VN_BITS>2。
- 结果：FER=1.0（前 30+ 包全失败，syndrome_weight 千级）。
- 分析：半步撤销不足以稳定，全步进攻仍造成大规模误翻，解码失败。

## 方案 A：存2算3扩展域（左移运算、右移回写；aggr 关闭）
- 修改：`f_update_vn_post` 对 VN_BITS<=2 左移一位运算，偏置居中，阈值=4，权重左移；回写右移；aggr 关闭。
- 结果：FER=1.0（前 30+ 包全失败，syndrome_weight 千级）。
- 分析：扩展域全步长过猛，偏置/阈值匹配不足，导致发散。

## 方案 D：分时后处理（前半关闭 post，后半开启；半步步长）
- 修改：`f_update_vn_post` 恢复半步逻辑；post_trigger/post_trigger2 在早期关闭，后期开启。
- 结果：FER≈0.11（100 包，11 失败，平均迭代≈199）。
- 分析：比前几方案略好，但仍陷入 error floor，推力不足；后处理分时未能扭转整体收敛。

## 方案 F：动态步长（趋势好全步，停滞半步；宏观趋势判据）
- 修改：趋势用上一迭代的 syndrome_weight，比当前不降则半步；步长按趋势切换；aggr 2bit 关闭。
- 结果：FER=1.0（快速发散，syndrome_weight 高）。
- 分析：宏观趋势判据未提供有效门控，权重/步长组合在 2bit 下仍过猛或失效。

## 下一步计划
- 尝试扩展+攻守分离变体（存2算3，未翻转全步，已翻转半步，aggr 关闭）。运行后补充结果。

## 方案（扩展+攻守分离变体：存2算3，未翻转全步，已翻转半步，aggr 关闭）
- 修改：VN_BITS<=2 时左移一位运算、居中偏置，阈值=4；未翻转分支全步（权重左移），已翻转分支半步（权重/2 奇数+pushing），回写右移；2bit aggr 关闭。
- 结果：FER=1.0（100 包全部失败，syndrome_weight 多在十几到千级间波动，未收敛）。
- 分析：扩展域攻守分离仍过于激进/失配，初始高权重 + 偏置导致误翻，大量包未收敛。

## 方案 G：2bit 固定梯度（flip_thr=2，强/弱=0/1），禁用 post，修正 likelihood_map
- 修改：
  - `f_likelihood_levels`：VN_BITS<=2 时设 `flip_thr=2`、`strong=0`、`weak=1`，初始 level\[0..3\]={0,1,1,1}，关闭 syndrome_weight 自适应。
  - `f_update_vn_post`：VN_BITS<=2 关闭 post/post2，未翻转用全步长 `likelihood+weight`，已翻转半步撤销 `(weight>>1)+odd&&pushing`。
  - `ldpc_ibex_parameters`：VN_BITS<=2 时 soft→level 映射改为 {0,1,2,2}，最强软值落在 level[0]。
  - `aggr`：VN_BITS<=2 强制关闭。
- 结果：FER=1.0（运行 120s 在第 35 包仍全失败，syndrome_weight ~1600，未收敛，进程被超时终止）。
- 分析：降低阈值与禁用后处理仍不足以推进收敛，全步进攻 + 半步撤销在 2bit 空间仍导致大规模误翻，syndrome 未下降。

## 方案 H：2bit 小权重“向上取半步”，w>=3 全步；撤销半步；flip_thr=2，post关闭
- 关键改动：
  - `f_update_vn_post`（src/ldpc_codec_test.cpp:575-588）：VN_BITS<=2 时 `delta_attack = (w>=3)?w:((w+1)>>1)`，撤销端 `delta_retract = (w>>1)+odd&&pushing`，未翻转更新 `likelihood + delta_attack - 1`，已翻转 `likelihood - delta_retract`，post/post2 禁用，aggr 对 2bit 关闭。
  - `f_likelihood_levels`（src/ldpc_codec_test.cpp:477-506）：维持 2bit 梯度 flip_thr=2，strong=0，weak=1，level 初值 {0,1,1,1}，不做 syndrome 自适应。
  - `likelihood_map`（src/ldpc_codec_test.cpp:1128-1143）：2bit 映射 {0,1,2,2}，最强软值→level[0]=strong。
- 结果：FER=0.12（100 包，平均迭代 141.18）。
- 分析：比方案 G 大幅好转但仍远高于目标；小权重上取半步提供有限推力，未能完全破除 error floor。

## 方案 I：方案H基础上重新启用 2bit 后处理（未翻转时按到 min，翻转边缘推到 thr+1）
- 关键改动：
  - `f_update_vn_post`（src/ldpc_codec_test.cpp:579-591）：VN_BITS<=2 保持攻守分离，新增 post_process 作用：`likelihood_new<thr` 时强制置 `min_likelihood`，`==thr` 时置 `thr+1`，post2 仍禁用。
  - 其他同方案 H：flip_thr=2，strong=0，weak=1，level 初值 {0,1,1,1}，映射 {0,1,2,2}，aggr 关闭。
- 结果：FER=0.18（100 包，平均迭代 201.11）。
- 分析：强压 post 让部分比特按到最小但整体推力仍不足，平均迭代变长且 FER 恶化。

## 方案 J：H 基础上关闭 2bit post，再加“高 syndrome 才开 aggr”
- 关键改动：
  - `f_update_vn_post`（src/ldpc_codec_test.cpp:579-586）：保持方案 H 的攻守分离与小权重上取半步，2bit 继续禁止 post/post2。
  - aggr 判据（src/ldpc_codec_test.cpp:2734-2741）：VN_BITS<=2 时仅当 `syndrome_weight_delayed>200` 才启用激进权重提升；否则沿用原 aggr。
- 结果：FER=1.0（运行 120s，仅前 44 包就全失败，syndrome_weight ~1100，完全未收敛）。
- 分析：激进门控未能提供有效推动，初期全失败表明仍然发散。

## 方案 K：在 J 基础上，未翻转去掉“-1”抑制（likelihood += delta_attack）
- 关键改动：`f_update_vn_post`（src/ldpc_codec_test.cpp:579-583）2bit 未翻转更新改为 `likelihood + delta_attack`（不减 1）；其余同方案 J（小权重上取半步，撤销半步，post 禁用，aggr 高权重才开）。
- 结果：FER=1.0（120s 内 35+ 包全部失败，syndrome_weight 约 1.5k，完全不收敛）。
- 分析：去掉抑制导致翻转更激进，结合低阈值直接发散。

## 方案 L：恢复原半步公式（delta=weight/2，未翻转含“-1”），post 仍禁用
- 关键改动：`f_update_vn_post`（src/ldpc_codec_test.cpp:579-585）2bit 分支改为 `delta=(w>>1)+odd&&pushing`，未翻转 `likelihood+delta-1`，已翻转 `likelihood-delta`；post/post2 禁用；aggr 高 synd 门控同方案 J。
- 结果：FER=1.0（120s 内前 45 包全失败，syndrome_weight ~1.0k，完全不收敛）。
- 分析：原半步推力在 flip_thr=2、strong=0 场景下仍明显不足，错误量几乎不下降。

## 方案 M：2bit 高阈值 flip_thr=3，全步进攻半步撤销（未翻转 +w，已翻转 w/2）
- 关键改动：
  - `f_likelihood_levels`：2bit flip_thr=3，strong=0，weak=1，level[2/3]=2。
  - `f_update_vn_post`：2bit 未翻转用 full weight，已翻转 half weight，post禁用，aggr 仅 VN_BITS>2。
- 结果：FER=1.0（120s 内前 35 包全失败，syndrome_weight~1600）。
- 分析：高阈值+全步导致大规模误翻，收敛更差。

## 方案 N：多数投票（enter=3 直接翻转，<2 撤销），翻转阈值=3
- 关键改动：`f_update_vn_post` 2bit 分支改为投票：未翻转 weight≥3 则设 `likelihood=flip_thr`，否则 `min`；已翻转 weight<2 则撤销为 `min`，否则保持翻转；post/agg 禁用。
- 结果：FER=1.0（约 50 包全失败，syndrome_weight ~2800，明显发散）。
- 分析：纯投票在噪声下直接发散，远劣于方案 H。

## 方案 P：双阈值磁滞+记忆（enter=2、confirm=3，连续高权才翻转；已翻转低权撤销）
- 关键改动：在 `ldpc_dec_bf_ibex` 中 2bit 分支绕过 `f_update_vn_post`，新增 `cand[j][k]` 记忆：未翻转且 `w>=2` 进入候选，下一次 `w>=3` 才翻转；已翻转且 `w<1` 撤销为 `min`，否则保持翻转。`flip_thr=3`，level 初值 0/1/2，post/aggr 禁用。
- 结果：FER=0.21（100 包，平均迭代 311.53），仍差于方案 H。
- 分析：磁滞/记忆降低误翻但推力不足，迭代数大增，未能突破 error floor。

## 方案 Q：Gallager-B 风格多数投票
- 关键改动：
  - `f_update_vn_post`（src/ldpc_codec_test.cpp）：2bit 下纯多数投票逻辑
    - 未翻转：weight>=3 直接翻转到 max，weight=2 推一步，weight<=1 保持
    - 已翻转：weight<=1 撤销到 min，weight=2 减一步但保持翻转，weight>=3 保持
  - post 禁用，aggr 禁用
- 结果：FER=1.0（100 包全失败，LDPC BER 80%，平均迭代 1024）。
- 分析：纯投票逻辑完全失效，错误大量累积，远劣于任何其他方案。

## 方案 R：H改进 - 边界微调
- 关键改动：
  - `f_update_vn_post`（src/ldpc_codec_test.cpp:566-592）：
    ```cpp
    if (!flipped) {
      if (weight >= 3)
        delta = weight;
      else if (weight == 2 && likelihood >= flip_threshold - 1)
        delta = 2;  // 边界处给足推力
      else
        delta = (weight + 1) >> 1;
    } else {
      delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
    }
    likelihood_new = flipped ? (likelihood - delta) : (likelihood + delta - 1);
    ```
  - post 禁用，使用纯 2-bit 域计算
- 结果：FER=0.11（100 包，11 失败，平均迭代 147.48）。
- 分析：边界增强（w=2 且 likelihood 接近阈值时给足 delta=2）略有帮助，但仍停留在 error floor。

## 方案 S：R + post 边界强推
- 关键改动：在方案 R 基础上启用 post 处理边界强推
- 结果：FER=0.12（100 包，略差于 R）
- 分析：post 边界强推没有明显改善

## 方案 T：迭代自适应 aggr 权重放大
- 关键改动：
  - `ldpc_dec_bf_ibex`（src/ldpc_codec_test.cpp:2735-2736）：
    ```cpp
    bool aggr = (VN_BITS <= 2 && iteration >= 100) || 
                ((ldpc_decoder_input.soft_bits > 0) &&
                 (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
    ```
  - 2-bit 在迭代次数 >= 100 后启用 aggr 权重放大（weight 0→0, 1→1, 2→3, 3→5, 4→7）
  - `f_update_vn_post` 使用 be_aggressive 参数控制全步/半步
- 结果：**FER=0.03**（332 包，10 失败，平均迭代 83.5）
- 分析：后期启用 aggr 权重放大能突破 error floor！从 FER=0.11 降到 0.03，大幅改善。
  - iter>=50：太早，FER=0.24（发散）
  - iter>=100：最佳，FER=0.03-0.04
  - iter>=150：略差，FER=0.046

## 当前最佳方案：T（iter>=100 启用 aggr）
- FER: 0.03 @ SNR=5.4
- 平均迭代: 83.5
- 关键代码修改：aggr 判据增加 `(VN_BITS <= 2 && iteration >= 100)` 条件

## 方案 T3：只用 aggr 权重放大，f_update_vn_post 用半步 ⭐最佳
- 关键改动：
  - `ldpc_dec_bf_ibex`（src/ldpc_codec_test.cpp:2735-2746）：
    - aggr 条件：`(VN_BITS <= 2 && iteration >= 100)`
    - aggr 权重放大：weight 0→0, 1→1, 2→3, 3→5, 4→7
  - `f_update_vn_post`（src/ldpc_codec_test.cpp:573）：
    - 禁用 `be_aggressive` 全步逻辑，始终用半步 `(weight + 1) >> 1`
    - 关键：权重放大后再用半步，避免双重放大导致误翻
    ```cpp
    if (0 && be_aggressive) {  // 禁用，只靠 aggr 权重放大
      delta = weight;
    } else {
      delta = (weight + 1) >> 1;  // 始终半步
    }
    ```
- 结果：**FER=0.0015**（6696 包，10 失败，平均迭代 59.1）
- 分析：**巨大突破！从 FER=0.11 降到 0.0015**
  - 原理：iter>=100 后启用 aggr 权重放大（2→3, 3→5, 4→7），但 f_update_vn_post 仍用半步
  - 这样 weight=3 放大到 5，半步 delta=(5+1)/2=3，正好能推过阈值
  - 不会双重放大导致误翻
  - LDPC BER = 2.88e-04，比方案 T 的 1.5% 好两个数量级

## 优化总结
| 方案 | FER | 平均迭代 | 关键改动 |
|------|-----|----------|----------|
| Baseline（半步） | 0.11 | ~200 | 半步公式 (weight>>1)+odd |
| H（小权重上取） | 0.12-0.18 | 141-200 | w>=3全步，否则上取半步 |
| T（aggr+全步双放大） | 0.03 | 83 | iter>=100 启用 aggr + f_update 全步 |
| **T3（aggr+半步）** | **0.0015** | 59 | iter>=100 启用 aggr，f_update 半步 |

## 方案 T3 完整代码修改

### 1. aggr 条件修改（ldpc_dec_bf_ibex 函数，约第 2735 行）
```cpp
// 原始代码
bool aggr = (ldpc_decoder_input.soft_bits > 0) &&
            (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev);

// 修改后
bool aggr = (VN_BITS <= 2 && iteration >= 100) || ((ldpc_decoder_input.soft_bits > 0) &&
            (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev));
```

### 2. f_update_vn_post 修改（约第 573 行）
```cpp
if (!flipped) {
  if (0 && be_aggressive) {  // 禁用全步，只靠 aggr 权重放大
    delta = weight;
  } else {
    delta = (weight + 1) >> 1;  // 始终半步
  }
} else {
  delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
}
```

### 原理说明
- iter < 100：使用半步，防止早期发散
- iter >= 100：启用 aggr 权重放大（2→3, 3→5, 4→7），然后在 f_update_vn_post 里再用半步
- weight=3 放大到 w=5，半步 delta=(5+1)/2=3，正好能从 likelihood=0 推到 2（flip_thr）
- 避免了双重放大导致的误翻，同时提供足够推力突破 error floor
