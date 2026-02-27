# IBEX BF 2bit（Row11）研究路线图（多 agent 同步探索）

> 目标：在“纯 2bit（每 VN 仅保留现有 2-bit likelihood + flipped）、不新增任何跨迭代 per‑VN 状态、只改 `src/ldpc_codec_test2.cpp`”的硬约束下，以 Row11@4.8（waterfall）为主战场降低 $LDPC\ FER$，并用 Row11@4.9 作为兜底回归（避免 error floor 明显回退）。

---

## 0. 硬约束（必须遵守）

1. **只允许改算法文件**：`src/ldpc_codec_test2.cpp`（Row11 路线）。  
2. **禁止新增任何跨迭代 per‑VN 状态**：tabu/momentum/候选记忆/FSM state/每 VN 计数器等一律禁止。  
   - 允许：全局状态；每列临时小缓存（仅在列级 escape/backtracking 内部使用）；每次迭代临时随机数（不跨迭代保存）。  
3. **构建**：只跑 `make ssd_fc_test2`（不要 `make clean`）。  
4. **评测口径（默认）**：
   - 筛选窗：`/tmp/Ibex_hd_row11_eval1000_0err.cnfg`（固定 1000 包，`max_err=0`，建议同命令跑 2 次）。  
   - 长窗：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000` + `max_err=10`，跑到 `FAIL CW=10`）。  
   - 指标：只关心 `[STATISTICS] LDPC FER` 与 `[STATISTICS] IBEX Decoder average iterations`。  
5. **记录**：
   - 每个新方案都写入 `doc/bf_2bit_report_codex.md`（改动点 + 一行命令 + 结果 + 分析）。  
   - `doc/bf_2bit_best_configs.md` 维护 Row8/Row11 各自 Top5（命令必须一行）。  

---

## 1. 当前 Row11@4.8 Baseline（纯 2bit、固定 1000 包）

统一口径：`./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.8`（该口径在 4.8 下几乎总是跑满 1000 包）。

- V46（w2-boost only-when-pushing）：`LDPC FER=3.300000e-02`, `avg_iter=254.711000`（log: `output/row11_v46_4p8_long10err_202602171350.log`）
- V53（PPBF escape）：`LDPC FER=3.400000e-02`, `avg_iter=253.403000`（log: `output/row11_v53_4p8_long10err_202602171355.log`）
- V51（rotate-k phase1 only）：`LDPC FER=4.200000e-02`, `avg_iter=260.209000`（log: `output/row11_v51_4p8_long10err_202602171345.log`）

当前主线（P1：phase1 aggr 非对称）：
- V61（V51 基座 + phase1 更早 aggr 门限 180/90/240）：`LDPC FER=2.700000e-02`, `avg_iter=259.177000`（log: `output/row11_v61_4p8_eval1000_0err_202602171358.log`）
  - 4.9 兜底回归（长窗到 `FAIL CW=10`）：`pkts=11108`, `LDPC FER=9.002521e-04`, `avg_iter=133.749460`（log: `output/row11_v61_4p9_long10err_202602171435.log`）
  - 4.7 waterfall 点（固定 300 包）：`pkts=300`, `LDPC FER=2.933333e-01`, `avg_iter=594.973333`（log: `output/row11_v61_4p7_eval300_0err_202602171834.log`）

当前最优（V74 = V61 + `IBEX_RESTART_PHASES=3`）：
- 4.8（`eval1000_0err`）：`pkts=1000`, `LDPC FER=2.600000e-02`, `avg_iter=254.117000`（log: `output/row11_v74_4p8_eval1000_0err_202602172008.log`）
- 4.9（长窗到 `FAIL CW=10`）：`pkts=26723`, `LDPC FER=3.742095e-04`, `avg_iter=132.728436`（log: `output/row11_v74_4p9_long10err_202602172027.log`）

---

## 1.1 当前 Row11@4.9 Top5（纯 2bit、长窗）

以 `doc/bf_2bit_best_configs.md` 为准。本路线图只作为工作台，不作为最终榜单。

> 提醒：长窗口径目前多为“10 错停止”，相对方差约 $\\sqrt{1/10}\\approx0.316$；Top5 排名仅作方向参考，必要时需重复跑或把 `max_err` 提到 30。

- Top1（V74：V61 基座 + `restart_phases=3`）：`pkts=26723`, `LDPC FER=3.742095e-04`, `avg_iter=132.728436`
- Top2（V51 rotate-k phase1 only）：`pkts=14999`, `LDPC FER=6.667111e-04`, `avg_iter=134.534436`
- Top3（V61 phase1 更早 aggr）：`pkts=11108`, `LDPC FER=9.002521e-04`, `avg_iter=133.749460`
- Top4（V53 PPBF escape）：`pkts=10573`, `LDPC FER=9.458054e-04`, `avg_iter=133.952048`
- Top5（V46 w2-boost only-when-pushing）：`pkts=10175`, `LDPC FER=9.828010e-04`, `avg_iter=134.198722`

---

## 2. 多 agent 任务分配（6 agents）

> 统一要求：每个 agent 产出“可落地映射（不加 per‑VN 状态）+ 2~3 条实验卡（完整一行命令）”。不要只给调参建议。

### Agent A：UP‑GDBF（论文）+ backtracking
- 资料：`doc/all core bf1/MinerU_markdown_Liu_等_-_2022_-_UP-GDBF_*.md`
- 关注：`t_active`（active iteration）、randomness as soft constraint、stall 触发、energy-based backtracking
- 产出：把“active iteration + backtracking”映射为**全局/列级门控**（不可引入 `l_i/hist_flag`）

### Agent B：NGDBF / re‑decoding（脚本）
- 资料：`doc/all core bf1/ldpc_mngdbf.m`, `doc/all core bf1/ldpc_reNGDBF.m`
- 关注：噪声扰动（不实现高斯，用 LFSR bit 近似）、多 phase 重启/多样性
- 产出：只在 `phase>0` + stall window 中引入“离散噪声门控”的可落地方案

### Agent C：AD‑GDBF（论文）+ diversity 框架
- 资料：`doc/all core bf1/MinerU_markdown_Brkic_et_al_2022_Adaptive_Gradient_Descent_*.md`
- 关注：多解码器多样性、失败后切换策略（不能用训练得到的 per‑VN 参数/动量）
- 产出：全局“模式切换表”（phase/iteration/syndrome_weight 分段）

### Agent D：FSM/Entropy 类（不可用点→可替代点）
- 资料：`doc/all core bf1/MinerU_markdown_entropy-27-00049_*.md`
- 关注：FSM per‑VN state 不可用；提炼 hesitation/determination 为**全局窗口化模式切换**
- 产出：stall 触发短窗口“escape 模式”，窗口结束回到保守模式

### Agent E：PPBF / PGDBF（脚本）
- 资料：`doc/all core bf1/ldpc_ppbf.m`, `doc/all core bf1/ldpc_pgdbf_with_momentum_bsc.m`
- 关注：PPBF 的 $p(E)$；动量/Tabu 不可用
- 产出：用“能量代理 + LFSR gate”实现 $p(E)$ 的分段/退火（只全局/列级）

### Agent F：TRGDBF/Tabu（脚本）+ 无记忆替代
- 资料：`doc/all core bf1/ldpc_trgdbf.m`
- 关注：tabu_list 不可用；如何用去相关调度/窗口化扰动替代“避免 ping‑pong”
- 产出：rotate‑k、列级 escape、backtracking、toggle cap 的组合边界与触发策略

---

## 2.1 Core BF 速记（只摘“可落地映射”）

> 目的：把论文/脚本里的“硬机制”翻译成 **全局/列级门控**，不引入任何 per‑VN 跨迭代状态。

- PPBF（`ldpc_ppbf.m`）：
  - 能量：$E_i=(hard_i\\oplus hard0_i)+\\sum_{j\\in N(i)} s_j$；按 $p(E_i)$ 概率翻转。
  - 映射：用 `w`（当前列内 weight）+ `flipped`（已有）+ `soft_unreliable` 做“离散能量代理”，只在 `phase>0` + stall window 做概率 boost（而不是直接概率翻转）。
- UP‑GDBF（Liu 2022）：
  - active iteration：$t>t_{active}$ 后才启用随机惩罚（避免早期随机误伤）。
  - backtracking：利用能量分类与回溯，只访问“更可能有用”的列，减少无效翻转/空转。
  - 映射：用 `*_MIN_ITER` + `IBEX_MODE_WIN` 近似 $t_{active}$；用 `IBEX_COL_GLOBAL_ESC` + `IBEX_COL_BACKTRACK_REQUIRE_IMPROVE` 做“列级 escape + 回滚保护”。
- (re)NGDBF（`ldpc_mngdbf.m`, `ldpc_reNGDBF.m`）：
  - 噪声扰动 + 多 phase 重启（re-decoding）显著影响 error floor。
  - 映射：不实现高斯，改用 LFSR gate 的离散噪声 bump（`IBEX_NGDBF_NOISE`），并依赖 `IBEX_RESTART_PHASES` 做多样性。
- FSM/momentum/Tabu（entropy paper、PGDBF w/M、TRGDBF）：
  - 依赖 per‑VN 历史（状态/动量/禁翻列表）→ **本项目硬约束下不可用**。
  - 映射：用 rotate‑k（调度去相关）+ stall window（短时加热）替代“避免 ping‑pong”。

---

## 3. Idea Cards（可落地映射）

> 规则：每条都必须写明“禁用 per‑VN 状态的映射方式”。

### Idea R01：全局窗口化模式切换（stall 触发）
- 来源：UP‑GDBF active iteration + FSM hesitation/determination 的“全局替代”
- 映射：当 `phase>0` 且进入 stall 后，打开长度为 `LEN` 的“escape 窗口”，窗口内提高 escape 概率/放宽门控；窗口结束恢复保守模式
- 风险：窗口过长会扰动正常包；窗口过短可能无效

### Idea R02：PPBF $p(E)$ 退火（tail 初期开大、尾端收小）
- 来源：PPBF（`doc/all core bf1/ldpc_ppbf.m`）+ UP‑GDBF random penalty（soft constraint）
- 映射：用极简能量代理 $E$（基于 `weight / soft_unreliable / pushing`），并让 $p(E)$ 随 iteration 分段
- 风险：热度过高会引入误翻；过低则几乎无效

### Idea R03：NGDBF‑lite 离散噪声（LFSR gate）
- 来源：mNGDBF（`doc/all core bf1/ldpc_mngdbf.m`）
- 映射：只在 `phase>0` + stall window，对少量候选施加离散噪声（例如对 $E$ +1）以增加 escape 概率；不保存任何 per‑VN 历史
- 风险：过度随机化会把长窗 FER 拉高

---

## 4. Experiment Cards（统一模板）

### 已完成（近期失败样本，避免重复踩坑）

- V59（PPBF escape + mode-window + anneal；固定窗就劣化；2026-02-15 16:16）：
  - 命令：`IBEX_MODE_WIN=1 IBEX_MODE_WIN_LEN=16 IBEX_MODE_WIN_TRIG_ITERS=8 IBEX_PPBF_ESC=1 IBEX_PPBF_ESC_ANNEAL=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.9`
  - 结果：`LDPC FER=2.000000e-03`, `avg_iter=132.617`
  - 日志：`output/row11_v59_4p9_eval1000_0err_202602151616.log`
  - 结论：窗口/退火仍偏激进或触发过早；短窗已变差，暂不继续长窗。

- V60（rotate-k + col-global-esc + mode-window；长窗劣化；2026-02-15 16:19-16:33）：
  - 命令：`IBEX_COL_GLOBAL_ESC=1 IBEX_COL_GLOBAL_ESC_ITERS=8 IBEX_COL_GLOBAL_ESC_MAX_TOGGLES=2 IBEX_MODE_WIN=1 IBEX_MODE_WIN_LEN=16 IBEX_MODE_WIN_TRIG_ITERS=8 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
  - 结果：`pkts=8222`, `LDPC FER=1.216249e-03`, `avg_iter=133.941255`
  - 日志：`output/row11_v60_4p9_long10err_202602151619.log`
  - 结论：列级 escape 触发/接受仍需更谨慎门控，否则会伤正常包；暂不作为主线。

### 日志命名约定
- `output/row11_<idea_id>_<snr>_<evaltype>_<yyyymmddHHMM>.log`
- macOS 没有 `stdbuf`：长窗建议用 `script -q <log> <cmd...>` 录制

### 统一基线（Row11@4.9）
- Baseline（V46）：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- Best1（V74）：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

### 实验模板（复制后填写）
- ID：
- 目的：
- 改动点（env 开关）：
- 命令（一行）：
- 结果（筛选窗 x2）：
- 结果（长窗到 `FAIL CW=10`）：
- 分析：
