# IBEX 2bit Best Configs（Codex）

当前各矩阵最优复现命令（必须一行，避免 env 分行失效）：
- Row8 Best1：`IBEX_2BIT_MODE=1 ./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`
- Row11 Best1：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`
- Row11@4.8 Best1（临时，单次 1000 包）：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`

说明：
- 本文件仅保留 Row8/Row11 各自 Top5；完整演进与失败方案见 `doc/bf_2bit_report_codex.md`
- `VN_BITS=2` 由 cnfg 的 `Ibex likelihood width = 2` 决定
- Row11 的命令包含 `IBEX_W2_STOCH=1`，叠加 AWGN 随机信道，单次 `1000` 包的 `LDPC FER` 会有明显波动；建议同一命令重复多次取统计
- Row11 长窗目前多为“10 错停止”，相对方差约 $\\sqrt{1/10}\\approx0.316$；Top5 排名仅供方向参考，必要时建议重复或把 `max_err` 提高到 30

## Row8：FER 对比（同一矩阵内比较）

### Best 1：Row8@5.4（2bit_mode=1）
- 译码器：`ssd_fc_test`（`src/ldpc_codec_test.cpp`）
- 配置：`config/Ibex_hd_row8.cnfg`
- 参数：`IBEX_2BIT_MODE=1`（其余默认）
- 结果：`pkts=24269`, `LDPC FER=4.120483e-04`, `avg_iter=48.782974`
- 运行命令：`IBEX_2BIT_MODE=1 ./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`

### Best 2：Row8@5.4（2bit_mode=0，默认）
- 译码器：`ssd_fc_test`（`src/ldpc_codec_test.cpp`）
- 配置：`config/Ibex_hd_row8.cnfg`
- 参数：默认（等价于 `IBEX_2BIT_MODE=0`）
- 结果：`pkts=13455`, `LDPC FER=7.432181e-04`, `avg_iter=48.805946`
- 运行命令：`./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`

### Best 3：Row8@5.4（2bit_mode=3）
- 译码器：`ssd_fc_test`（`src/ldpc_codec_test.cpp`）
- 配置：`config/Ibex_hd_row8.cnfg`
- 参数：`IBEX_2BIT_MODE=3`（其余默认）
- 结果：`pkts=6822`, `LDPC FER=1.465846e-03`, `avg_iter=49.794782`
- 运行命令：`IBEX_2BIT_MODE=3 ./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`

### Best 4：Row8@5.4（2bit_mode=2）
- 译码器：`ssd_fc_test`（`src/ldpc_codec_test.cpp`）
- 配置：`config/Ibex_hd_row8.cnfg`
- 参数：`IBEX_2BIT_MODE=2`（其余默认）
- 结果：`pkts=5547`, `LDPC FER=1.802776e-03`, `avg_iter=49.691184`
- 运行命令：`IBEX_2BIT_MODE=2 ./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`

### Best 5：Row8@5.4（2bit_mode=1 + 更保守 aggr 门限）
- 译码器：`ssd_fc_test`（`src/ldpc_codec_test.cpp`）
- 配置：`config/Ibex_hd_row8.cnfg`
- 参数：`IBEX_2BIT_MODE=1 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280`
- 结果：`pkts=5000`, `LDPC FER=2.200000e-03`, `avg_iter=62.474800`
- 运行命令：`IBEX_2BIT_MODE=1 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 ./ssd_fc_test LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`

## Row11：FER 对比（同一矩阵内比较）

> 注意：方案 V47（`IBEX_W2_CAND_STRONG=1`）引入额外 1-bit/VN 状态，按“纯 2bit”口径不纳入本 Top5，仅保留在 `doc/bf_2bit_report_codex.md` 作为记录。

### Best 1：Row11@4.9（V74：V61 基座 + `restart_phases=3` ⭐）
- 译码器：`ssd_fc_test2`（`src/ldpc_codec_test2.cpp`）
- 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`）
- 参数：V61（phase1 更早 aggr）+ `IBEX_RESTART_PHASES=3`
- 结果（长窗，>=1000 包 + 10 错停止）：`pkts=26723`, `LDPC FER=3.742095e-04`, `avg_iter=132.728436`
- 运行命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

### Best 2：Row11@4.9（V51 rotate-k phase1 only）
- 译码器：`ssd_fc_test2`（`src/ldpc_codec_test2.cpp`）
- 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`）
- 参数：在 V46 基础上新增 `IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1`（仅 retry phase 生效，改变列内 bit 扫描起点）
- 结果（长窗，>=1000 包 + 10 错停止）：`pkts=14999`, `LDPC FER=6.667111e-04`, `avg_iter=134.534436`
- 运行命令：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

### Best 3：Row11@4.9（V61 phase1 更早 aggr ⭐）
- 译码器：`ssd_fc_test2`（`src/ldpc_codec_test2.cpp`）
- 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`）
- 参数：在 V51 基础上新增 `IBEX_PHASE1_AGGR_ITER_HI/LO/SYND_TH=180/90/240`（仅 retry phase 生效，更早进入 aggr）
- 结果（长窗，>=1000 包 + 10 错停止）：`pkts=11108`, `LDPC FER=9.002521e-04`, `avg_iter=133.749460`
- 运行命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

### Best 4：Row11@4.9（V53 PPBF escape）
- 译码器：`ssd_fc_test2`（`src/ldpc_codec_test2.cpp`）
- 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`）
- 参数：在 V46 基础上新增 `IBEX_PPBF_ESC=1`（仅 retry phase + stall 后生效的概率 escape）
- 结果（长窗，>=1000 包 + 10 错停止）：`pkts=10573`, `LDPC FER=9.458054e-04`, `avg_iter=133.952048`
- 运行命令：`IBEX_PPBF_ESC=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

### Best 5：Row11@4.9（V46 w2-boost only-when-pushing）
- 译码器：`ssd_fc_test2`（`src/ldpc_codec_test2.cpp`）
- 配置：`/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_sim_num=1000`, `max_err=10`）
- 参数：`w2_boost_only_when_pushing=1`，`2bit_mode=2`，`A/B/C = 240/120/280`，`IBEX_PUSH_DYNAMIC=1`，`IBEX_W2_STOCH=1`，`IBEX_RESTART_PHASES=2`
- 结果（长窗，>=1000 包 + 10 错停止）：`pkts=10175`, `LDPC FER=9.828010e-04`, `avg_iter=134.198722`
- 运行命令：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.9`

## Row11@4.8（waterfall）：FER 对比（固定 1000 包）

说明：
- 这里的“固定 1000 包”使用 `/tmp/Ibex_hd_row11_eval1000_0err.cnfg`（`max_sim_num=1000, max_err=0`）或 `/tmp/Ibex_hd_row11_eval1000.cnfg`（`max_err=10`）；在 `SNR=4.8` 下两者通常都会跑满 1000 包，差别主要是统计输出的“错误门槛”。
- 由于 AWGN + `IBEX_W2_STOCH=1` 波动仍存在，本榜单暂按“单次 1000 包”记录，后续若需要可重复 2~3 次取均值。

### Best 1（当前）：Row11@4.8（V74：V61 基座 + `restart_phases=3`）
- 译码器：`ssd_fc_test2`（`src/ldpc_codec_test2.cpp`）
- 配置：`/tmp/Ibex_hd_row11_eval1000_0err.cnfg`
- 参数：V61（phase1 更早 aggr）+ `IBEX_RESTART_PHASES=3`
- 结果：`pkts=1000`, `LDPC FER=2.600000e-02`, `avg_iter=254.117000`
- 4.9 回归（长窗，>=1000 包 + 10 错停止）：`pkts=26723`, `LDPC FER=3.742095e-04`, `avg_iter=132.728436`（log: `output/row11_v74_4p9_long10err_202602172027.log`）
- 运行命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=3 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000_0err.cnfg AWGN 4.8`

### 补充：Row11@4.7（V61，eval300_0err）
- 配置：`/tmp/Ibex_hd_row11_eval300_0err.cnfg`
- 结果：`pkts=300`, `LDPC FER=2.933333e-01`, `avg_iter=594.973333`
- 运行命令：`IBEX_PHASE1_AGGR_ITER_HI=180 IBEX_PHASE1_AGGR_ITER_LO=90 IBEX_PHASE1_AGGR_SYND_TH=240 IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval300_0err.cnfg AWGN 4.7`

### Baseline：Row11@4.8（V46 / V53 / V51）
- V46（w2-boost only-when-pushing）：`pkts=1000`, `LDPC FER=3.300000e-02`, `avg_iter=254.711000`
  - 命令：`IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.8`
- V53（PPBF escape）：`pkts=1000`, `LDPC FER=3.400000e-02`, `avg_iter=253.403000`
  - 命令：`IBEX_PPBF_ESC=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.8`
- V51（rotate-k phase1 only）：`pkts=1000`, `LDPC FER=4.200000e-02`, `avg_iter=260.209000`
  - 命令：`IBEX_ROTATE_K=1 IBEX_ROTATE_K_PHASE1_ONLY=1 IBEX_W2_BOOST_ONLY_WHEN_PUSHING=1 IBEX_2BIT_MODE=2 IBEX_AGGR_ITER_HI=240 IBEX_AGGR_ITER_LO=120 IBEX_AGGR_SYND_TH=280 IBEX_PUSH_DYNAMIC=1 IBEX_W2_STOCH=1 IBEX_RESTART_PHASES=2 ./ssd_fc_test2 LDPC /tmp/Ibex_hd_row11_eval1000.cnfg AWGN 4.8`
