# IBEX BF 1-bit 方案探索计划

## 任务目标
探索 1-bit LDPC BF（Bit-Flipping）解码方案，评估其在 AWGN 5.4 信道下的纠错性能。

## 性能参考基准
| 配置 | SNR | FER | 备注 |
|------|-----|-----|------|
| 2-bit（baseline 半步） | 5.4 | ~0.11 | error floor，推力不足 |
| 2-bit（T3 方案 aggr+半步） | 5.4 | 0.0015 | 当前工程内较优配置 |
| 2-bit（用户提供） | 5.4 | 5e-4 | 更激进的 2-bit 调优结果 |
| 3-bit（用户提供） | 5.4 | 4.77e-6 | 高精度 BF 参考上限 |
| **1-bit 目标** | 5.4 | **< 0.05** | 初步筛选门槛（100 packets） |

## 1-bit 方案的核心挑战与初步思路

### 1. 存储空间极度受限（梯度消失）
- 1-bit 只能表示 2 个状态：0/1，对应“未翻转/已翻转”两种离散状态。
- 2-bit/3-bit 中的 $likelihood\in[0,3]$ 或更高精度空间，可以通过 $delta$ 的加减实现“从强到弱、再到翻转”的多级渐进更新。
- 在 1-bit 下，这些“强/弱”“半步推进”全部消失，只剩下“翻”与“不翻”的布尔判决，数值空间本身无法再承载细粒度信息。

### 2. 与 2-bit/3-bit 的行为差异
- 2-bit：依赖 $delta=\lfloor weight/2\rfloor+((weight\bmod 2)\land pushing)$ 和半步逻辑，在极小的 0~3 空间内仍然维持了“弱/强 + 渐进靠近阈值”的结构；配合激进加权和后处理，可以在 $5.4$dB 下做到 $FER\approx 5\times 10^{-4}$ 级别。
- 3-bit 及以上：使用 $delta=weight$ 的全步逻辑，加上 $flip\_thr$ 和多级 level 的设计，形成了“强→弱→边界→已翻转”的完整状态机，性能可以低到 $FER\approx 4.77\times 10^{-6}$。
- 1-bit 若直接照搬 2-bit 的更新公式，只是把数值范围压到 $\{0,1\}$，则所有“渐进”和“缓冲带”都会退化，容易出现两种极端：要么很难翻转（几乎不纠错），要么频繁乒乓（直接发散）。

### 3. 从 2-bit 映射到 1-bit 的问题
| 2-bit likelihood | 含义 | 1-bit 映射（直觉） |
|------------------|------|--------------------|
| 0 | 最强未翻转 | 0 |
| 1 | 弱未翻转 | 0 |
| 2 | 弱已翻转 | 1 |
| 3 | 强已翻转 | 1 |

- 如果简单采用上表映射，则 $flip\_thr$ 在 1-bit 下只能取 $1$，即 $likelihood\ge 1$ 就视为翻转。
- 这意味着无法再通过“$likelihood=1,2$ 等中间档位”实现从弱到阈值的渐进过程，只能在单次更新中直接跨越翻转与否，这会显著放大更新策略的风险。

### 4. soft 信息与后处理在 1-bit 下的尴尬境地
- 2-bit/3-bit 中，soft 信息先通过 `likelihood_map` 映射到初始 level，再通过自适应系数、插值等机制细化梯度，soft 信息对后续每一次更新都有持续影响。
- 在 1-bit 中，如果初始 $likelihood$ 全部设为 0，soft 信息完全被丢弃；如果把“软上可疑”的 bit 直接初始化为 1（已翻转），又容易一开始就把大量比特推到错误侧。
- 同理，2-bit 的后处理依赖 $flip\_thr\pm 1$ 的数值扰动；在 1-bit 空间里，这种扰动只能变成“翻/不翻”的硬切换，需要重新设计触发条件和作用方式。

### 5. 关键函数与 1-bit 改造方向（框架级而非简单压缩）
| 函数 | 2-bit/3-bit 逻辑 | 1-bit 改造方向（框架级） |
|------|------------------|--------------------------|
| `f_likelihood_levels` | 根据 $VN\_BITS$、$syndrome\_weight$、soft 生成多级 level 及 $flip\_thr$ | 在 $VN\_BITS=1$ 时仅作为“翻转初值/偏置”的配置器，输出的 level 只能是 0/1，并弱化动态自适应 |
| `f_update_vn_post` | 半步/全步增量更新，配合 $flip\_thr$、后处理形成多级状态机 | 在 $VN\_BITS=1$ 时重写为“翻转状态机”：显式使用双阈值、趋势门控、分时策略，而不是简单的整数加减 |
| `likelihood_map` | soft→$\{0,1,2,3\}$ 等多级初值 | 在 $VN\_BITS=1$ 时只提供“初始翻转偏置”，控制哪些 bit 起步更容易翻转 |
| aggr 权重放大 | $2\rightarrow 3, 3\rightarrow 5, 4\rightarrow 7$ 等，对 $delta$ 放大 | 在 $VN\_BITS=1$ 时需要重新解释为“允许翻转的条件/概率”而非简单权重放大 |

### 6. 初步总体思路（1-bit 不等于“弱化版 2-bit”）
- 不直接把 2-bit/3-bit 的 $likelihood$ 更新公式裁剪成 1bit，而是把 1-bit 看作“一个高度约束的布尔状态机”，通过**阈值、时间、全局趋势和 soft 偏置**来弥补梯度的缺失。
- 首轮探索集中在“框架级方案”（不同的状态机结构），而不是细节调参；待某个框架在 $5.4$dB 下 $FER\le 0.05$ 时，再在该框架内做参数扫描。
- 只借鉴 `doc/all core bf` 中经典 BF 方案的**机制**（如 NGDBF 的扰动、PGDBF 的概率翻转、带记忆的动量项等），但不把 IBEX 直接改造成这些算法，因为现有 IBEX 3-bit/2-bit 方案在本码下已经明显优于这些通用算法。

---

## 设计约束与已确认决策

### 1. 源码修改范围
- 不改动现有“主线”源码文件（`src/ldpc_codec.cpp` 等），避免影响既有 2-bit/3-bit 流程。
- 仅在以下位置工作：
  - 复制 `src/ldpc_codec.cpp` 为 `src/ldpc_codec_test3.cpp`，用于 1-bit 方案实验；
  - 修改 `Makefile` 新增 `ssd_fc_test3` 目标，链接 `ldpc_codec_test3.cpp`；
  - 更新 `doc/agent_1bit_BF.md` 与 `doc/bf_1bit_report_codex.md` 记录过程与结果。

### 2. 自动探索与实验流程
- 不使用自动脚本扫参数网格，不写专门的“参数扫描工具”。
- 严格采用“手动版本迭代”的节奏：
  1. 在 `ldpc_codec_test3.cpp` 中实现某一版 1-bit 框架；
  2. 编译生成 `ssd_fc_test3`；
  3. 运行一次统一命令：`./ssd_fc_test3 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`；
  4. 记录 $FER$ 和关键现象，然后再修改代码进入下一版方案。
- 每完成 5 个版本的探索，暂停并汇报阶段性结果，再决定后续方向。

### 3. 评价准则与记录策略
- 统一评价点：$SNR=5.4$dB，配置文件 `Ibex_hd_row8.cnfg` 中 $VN\_BITS=1$，默认“100 packets 或 10 错停止”。
- 初筛门槛：$FER\le 0.05$ 视为“潜在可行方案”，需要在更多 packets 下复查（由人工加包数验证）。
- 记录策略：
  - 所有版本（包括 $FER=1$ 的彻底失败方案）都在 `doc/bf_1bit_report_codex.md` 中至少保留一条记录（版本名、核心思路、FER 结果），以标记失败路径；
  - 对于 $FER\le 0.05$ 的版本，额外写详细小节：关键代码改动、算法框架、实验结果与分析。

### 4. 软信息、后处理和 aggr 的原则
- soft 信息：在 1-bit 中主要用于**初始化偏置**（如哪些比特更容易起步翻转），而非长期连续的多级梯度；优先从“是否初始设为 flipped”以及“阈值选择”两个方向使用 soft。
- post\_process：在 1-bit 空间下，更多被重解释为“停滞阶段的随机扰动/解锁机制”，而不是数值上的 $flip\_thr\pm 1$；触发条件会与迭代次数和 $syndrome\_weight$ 紧密绑定。
- aggr：不直接把 2-bit/3-bit 的权重放大公式照搬到 1-bit，而是考虑将“激进阶段”编码为：允许更低阈值翻转、更宽松的翻转条件或更强的随机扰动。

---

## 分步骤计划（当前版本）

### 阶段 1：环境准备
1. 复制 `src/ldpc_codec.cpp` → `src/ldpc_codec_test3.cpp`，保持初始版本与主线完全一致。
2. 修改 `Makefile` 添加 `ssd_fc_test3` 目标，链接 `ldpc_codec_test3.cpp`，其他目标不改。
3. 编译并用当前配置（$VN\_BITS=1$）跑一版 baseline：`./ssd_fc_test3 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`，记录“原始 1-bit C++ 实现”的 FER。

### 阶段 2：1-bit 框架级改造（首批 5 个版本）
**目标**：在仅修改 `ldpc_codec_test3.cpp` 的前提下，设计并验证 5 个“结构上差异明显”的 1-bit BF 框架，找到至少 1 个 $FER\le 0.05$ 的候选。

1. **v1：静态双阈值 Gallager-B 基线**
   - `VN_BITS=1` 时，将 likelihood 显式视为布尔状态：0 表示未翻转，1 表示已翻转，$flip\_thr=1$。
   - 使用静态双阈值：$T\_{flip}=3$、$T\_{unflip}=1$，只依据局部权重 $weight$ 决定“翻/不翻/撤销”，不使用 soft、pushing、post、aggr 等高级机制。

2. **v2：静态双阈值 + 全局趋势门控（pushing）**
   - 在 v1 的基础上利用 $pushing=(syndrome\_weight\_delayed\ge prev\_sw)$：
     - 收敛阶段（$pushing=false$）：禁止新翻转，只允许撤销；
     - 停滞/恶化阶段（$pushing=true$）：允许按 $T\_{flip},T\_{unflip}$ 判据引入新翻转。
   - 通过该门控避免在已经“往好方向走”的阶段引入多余翻转。

3. **v3：趋势门控 + soft 初始化偏置**
   - 在 v2 的基础上引入 soft 初始化（当前配置下 `soft_bits=1`）：
     - `bit_questionable==0`：初始 likelihood=0（未翻转）；
     - `bit_questionable==1`：初始 likelihood=1（已翻转）。
   - 后续迭代仍使用 v2 的门控和双阈值，把 soft 信息只用在“起步时的偏置”上。

4. **v4：趋势门控 + 分时阈值（前期保守、后期激进）**
   - 在 v2/v3 框架上进一步引入迭代分时：
     - 早期（例如 $iteration<100$）：使用更高的翻转阈值 $T\_{flip}^{early}=4$，更保守翻转；
     - 后期（$iteration\ge 100$）：降低阈值到 $T\_{flip}^{late}=3$，加大翻转力度。
   - 结合 pushing，使“何时允许翻转”同时受全局趋势和迭代阶段的双重约束。

5. **v5：趋势门控 + 分时阈值 + 停滞扰动（简化后处理）**
   - 在 v4 基础上，重新解释 `post_process/post_process2` 为“在停滞阶段触发的随机扰动”：
     - 当 $iteration\ge post\_iteration$ 且 $syndrome\_weight$ 长时间徘徊在某一区间时，使用 PRNG 对部分满足特定 $weight$ 条件（如 $w=2$）的比特进行随机翻转或撤销。
   - 目标是在 1-bit 空间中找出一套“低频扰动”机制，用少量随机翻转打破 error floor 或长 plateau。

### 阶段 3：进一步优化与机制扩展
- 在首批 5 个框架中挑选表现最好的一两个（特别是 $FER\le 0.05$ 或接近的方案），针对以下方向做精细优化：
  - 细化翻转/撤销阈值（$T\_{flip},T\_{unflip}$）和分时边界；
  - 引入轻量记忆机制（例如对多次翻转失败的比特提高撤销门槛）；
  - 试探性引入 NGDBF/PGDBF 风格的扰动或动量项，但仅借鉴机制，不整体更换为这些算法。

### 阶段 4：性能评估与报告
- 每完成 5 个版本的实现与仿真，整理一次阶段性小结（包括失败路径），与人工确认下一步探索方向。
- 对 $FER\le 0.05$ 的方案，在 `doc/bf_1bit_report_codex.md` 中写详细评估小节并标记为候选。
- 最终在 `doc/bf_1bit_report_codex.md` 中总结 1-bit 方案在当前码率和噪声条件下的“可行性能边界”。

---

## 测试命令
```bash
./ssd_fc_test3 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4
```
- 默认 100 包或 10 错停止
- 观察 `[STATISTICS] LDPC FER` 和 `avg_iteration`

---

## 预期结果与风险评估

### 乐观情况
- 某个 1-bit 方案达到 FER < 0.1，证明 1-bit BF 可用于某些场景

### 现实预期
- 1-bit 方案 FER 在 0.05~0.5 之间，需要精细调参

### 悲观情况
- 所有 1-bit 方案 FER > 0.1，1-bit BF 在当前码率/噪声下不可行

---

## 配置说明

- **VN_BITS 配置**：通过 cnfg 文件设置为 1-bit，无需修改代码中的宏定义
- **迭代上限**：使用默认 max_iteration，无需调整
- **性能目标**：FER < 0.05 为初步目标

---

## 关于借鉴经典 BF 算法的说明

- `doc/all core bf/` 目录下包含了一些经典 BF 变种（如 PGDBF/NGDBF/带动量的 GDBF 等）的 Matlab 实现，这些算法是通用框架的优秀代表，但在本 IBEX 码上的性能通常弱于当前的 3-bit/2-bit IBEX 专用方案。
- 在 1-bit 方案设计中，只计划从这些算法中**借鉴具体机制**，例如：
  - 在停滞阶段使用随机扰动打破局部极小（NGDBF 思路）；
  - 利用“动量”或“翻转历史”抑制乒乓（带动量 BF 思路）；
  - 对不同权重或局部结构使用差异化的翻转概率。
- 不会将 IBEX 直接改写为这些通用算法的完整 C++ 版本，而是以现有 IBEX 结构为主线，在 1-bit 环境下做适度机制移植和裁剪。

---

## 文件结构
```
IBEX/
├── src/
│   ├── ldpc_codec.cpp          # 原始版本
│   ├── ldpc_codec_test.cpp     # 2-bit T3 方案
│   ├── ldpc_codec_test2.cpp    # 其他测试版本
│   └── ldpc_codec_test3.cpp    # 1-bit 方案测试
├── doc/
│   ├── agent.md                # 2-bit 方案候选清单
│   ├── bf_2bit_report.md       # 2-bit 方案测试记录
│   ├── bf_2bit_report_codex.md # 2-bit 方案 Codex 记录
│   ├── agent_1bit_BF.md        # 1-bit 方案计划与过程 [本文档]
│   └── bf_1bit_report_codex.md # 1-bit 方案测试记录（结果与分析）
└── Makefile                    # 已添加 ssd_fc_test3 目标
```

---

## 首批 1-bit 框架实验小结（v0~v5）

- 实验配置：`./ssd_fc_test3 LDPC config/Ibex_hd_row8.cnfg AWGN 5.4`，$VN\_BITS=1$，100 包或 10 错停止。
- 方案集合：
  - v0：`ldpc_codec.cpp` 原始逻辑直接复制，1-bit 无专门适配；
  - v1：静态双阈值（$T\_{flip}=3,T\_{unflip}=1$），无趋势/soft/post；
  - v2：v1 + pushing 门控（仅在 pushing 为真时允许新翻转）；
  - v3：v2 + soft 初始化偏置（`bit_questionable==1` 起步设为 flipped）；
  - v4：v2 + 分时阈值（迭代早期更保守，后期稍激进）；
  - v5：v4 + 基于 `post_process2` 的低频随机扰动。
- 实验结果（详细见 `doc/bf_1bit_report_codex.md`）：
  - 所有方案在 $SNR=5.4$dB 下均为 $LDPC\ FER = 1.0$，平均迭代数 $1024$；
  - v0/v4/v5 的 $LDPC\ BER$ 与 $RAW\ BER$ 接近，几乎“不解码”；
  - v1/v2/v3 的 $LDPC\ BER$ 分别逼近 $0.5$ 和 $1$，表现为“全局发散”。
- 阶段性结论：
  - 简单的 1-bit 静态双阈值、趋势门控、分时策略、soft 初始化和无记忆扰动都不足以在当前码率/SNR 下取得 $FER \le 0.05$；
  - 下一步需要引入更强的结构化机制（翻转记忆、动量、局部结构判据等），而不仅是对 2-bit/3-bit 策略的直接压缩或轻微修补。
