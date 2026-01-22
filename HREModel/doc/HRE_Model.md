# HREModel：SSD\_FC / transceiver 多版本与“特殊仿真场景”梳理

本文用于把当前工程里**不同可执行程序（SSD\_FC 族）**与**不同信道实现（transceiver 族）**的行为对齐，避免把同名配置项（如 `hre_mode`）在不同版本中误解为同一语义。

约定符号（便于统一描述）：

- 码字长度 $N$（工程中通常为 `dsp_blk_len`，例如 38272）
- 发送比特 $b\in\{0,1\}$
- BPSK 映射 $x = -(2b-1)$，即 $b=0\to x=+1,\; b=1\to x=-1$
- AWGN：$y=x+\sigma n$，$n\sim\mathcal{N}(0,1)$，硬判决 $\hat b = \mathbf{1}\{y<0\}$
- 理论 RBER（仅 AWGN/TAWGN 下定义）：$rber = 0.5\cdot(1+\mathrm{erf}(-1/(\sigma\sqrt{2})))$

---

## 1. Makefile 里"版本"的真实含义：4 个可执行目标

当前 `Makefile` 明确生成 4 个目标：

| 可执行文件 | 入口源文件 | transceiver 实现 | 编译宏 | 是否跑编码/译码 | 主要用途 |
|---|---|---|---|---|---|
| `ssd_fc_hre` | `src/SSD_FC.cpp` | `src/transceiver.cpp` | 无 | 是 | 主链路全仿真（包含 LDPC 编/译码），使用"旧式 HRE 注入语义" |
| `ssd_fc_hre2` | `src/SSD_FC_hre2.cpp` | `src/transceiver.cpp` | `-DHRE2_MODE` | 是 | 主链路全仿真（包含 LDPC 编/译码），使用"新式 HRE2 注入语义"，并输出更多 HRE 统计 |
| `ssd_fc_hre_only` | `src/SSD_FC_hre_only.cpp` | `src/transceiver_hre_only.cpp` | 无 | 否 | **仅信道+量化+统计**，不跑 LDPC（用于快速统计 RAW/FBC/HRE-like） |
| `ssd_fc_flip_eq` | `src/SSD_FC_flip_eq.cpp` | `src/transceiver_flip_eq.cpp` | 无 | 是 | 主链路全仿真（包含 LDPC 编/译码），使用 `TRUNC_AWGN/TAWGN` fixed-FBC 条件仿真（支持 `target_fbc` 覆盖），并输出详细电压 bin 统计 |

要点（容易混淆）：

- `ssd_fc_hre2` **不是**"只跑信道"的版本：它链接了 `src/ldpc_codec.cpp`，完整跑编/译码与 FER。
- `ssd_fc_hre` 与 `ssd_fc_hre2` 使用**同一个** `src/transceiver.cpp`，但通过 `HRE2_MODE` 宏在编译期切换两套 HRE 注入/统计语义。
- `ssd_fc_hre_only` 使用一套**独立复制**的 `transceiver_hre_only.*`，其行为并不等同于 `transceiver.cpp` 的任一分支（尤其在 HRE-like 统计上）。
- `ssd_fc_flip_eq` 使用**独立实现**的 `transceiver_flip_eq.*`，专门用于 fixed-FBC（TAWGN）条件仿真与电压 bin 分布统计，并支持用 `target_fbc` 指定固定硬错误数 $k$（`target_fbc==0` 则默认 $k=k_{\text{theo}}$）。

---

## 2. 三套 transceiver 的核心差异（最重要的分歧来源）

### 2.1 `src/transceiver.cpp`（被 `ssd_fc_hre` 与 `ssd_fc_hre2` 共用）

它同时包含：

1. 多信道：`CLEAN/AWGN/BSC/ERR_INJ/MAX_ERR/TRUNC_AWGN(TAWGN)`；
2. 两套 HRE 注入语义：
   - 未定义 `HRE2_MODE`：旧式 HRE
   - 定义 `HRE2_MODE`：新式 HRE2（"确定数量 bit-flip"）

### 2.2 `src/transceiver_hre_only.cpp`（仅被 `ssd_fc_hre_only` 使用）

它是"简化版 transceiver"，特点是：

- **不支持** `TRUNC_AWGN(TAWGN)`（`enum ch_model` 里也没有该枚举值）
- AWGN 下的 HRE 注入语义固定为"翻转且不加噪声"，并且**忽略** `hre_mode/hre_model` 的多种取值
- HRE-like 统计写死为"bin0/bin7"（对不同 Vref 排列不一定泛化）

### 2.3 `src/transceiver_flip_eq.cpp`（仅被 `ssd_fc_flip_eq` 使用）

它是 `ssd_fc_flip_eq` 专用的 transceiver，特点是：

- **支持** `TRUNC_AWGN/TAWGN` 信道模型，并允许通过 `target_fbc` 指定固定硬错误数 $k$（见 §3.5.2.1）
- 在 `TRUNC_AWGN` 模式下：构造一个**纯 TAWGN** 帧，使得相对于 $x_0$ 的硬判决错误数恒为 $k$
- 提供详细的电压 bin 统计：`tx0_bin7/5/3/1`、`tx1_bin0/2/4/6`、`tx0_self/tx1_self` 等（便于做 bin 分解分析）
- 适合做“固定错误重量/固定 FBC”下的条件性能曲线与 bin-level breakdown

因此：只要你切换可执行文件（`ssd_fc_hre` vs `ssd_fc_hre2` vs `ssd_fc_hre_only` vs `ssd_fc_flip_eq`），同一份 `.cnfg` 的 `hre_mode` 可能对应完全不同的物理含义。

---

## 3. `ch_model`（信道模式）与对应的“特殊仿真场景”

### 3.1 CLEAN

- 行为：$y=x$（无噪声）
- HRE：不注入（代码只在 AWGN 分支处理 HRE）
- 适用：基线功能验证

### 3.2 AWGN（标准主用分支）

所有版本都用 BPSK+AWGN 生成 `rx_blk`，但 **HRE 的定义因版本而异**（见第 4 节）。

### 3.3 TRUNC\_AWGN / TAWGN（截断 AWGN，用于“固定 FBC”的条件仿真）

在 `ssd_fc_hre` / `ssd_fc_hre2`（`src/transceiver.cpp`）中支持该模式。  
`ssd_fc_flip_eq` 也支持 `TRUNC_AWGN/TAWGN`，并允许用 `target_fbc` 覆盖固定 FBC 的目标值（`target_fbc==0` 则默认 $k=k_{\text{theo}}$，见 §3.5.2.1）。

TAWGN 的目标是：在 AWGN 统计意义下，直接从条件分布 $\{y_i\}\mid E=k$ 采样，使得每帧硬判决错误数恒为
$k=\lfloor N\cdot rber\rfloor$，从而避免逐帧拒绝采样。

- 行为：先随机选择大小为 $k$ 的错误位置集合 $S$，再对 $i\in S$ 从“错误侧截断高斯”采样，对 $i\notin S$ 从“正确侧截断高斯”采样（工程实现为单符号拒绝采样 `sample_trunc_awgn`）
- 当前实现明确说明：**TAWGN 分支不叠加 HRE**

更完整的推导与实现说明见 `doc/truncated_AWGN.md`。

### 3.4 BSC / ERR\_INJ / MAX\_ERR（离散翻转类）

- `BSC`：以概率 `ber` 翻转符号（硬翻转，非高斯幅度）
- `ERR_INJ`：每帧固定注入 `err_num` 个随机位置错误
- `MAX_ERR`：每帧注入 $[0,err\_num]$ 均匀随机个错误

注意：这些分支目前都**不检查** `hre_vec`，因此**不会注入 HRE**（即使配置文件里 `hre_bit>0`）。

### 3.5 `ssd_fc_flip_eq` 的 TRUNC\_AWGN/TAWGN：固定 FBC + bin 分解

仅 `src/transceiver_flip_eq.cpp` 支持（因此只对 `ssd_fc_flip_eq` 有效）。

#### 3.5.1 设计动机

当你想研究“在硬错误数固定为 $k$ 的条件下，LDPC 的 FER/迭代次数/失败模式如何变化”时，标准 AWGN 的逐帧 Monte Carlo 会受到 $E$（每帧错误数）自然波动的影响，尤其在尾部很难采样。  
TAWGN（截断 AWGN）通过直接从条件分布 $\{y_i\}\mid E=k$ 采样，使得每帧 FBC 固定，更适合做 **FER vs. error weight/FBC** 与 **bin-level breakdown**。

#### 3.5.2 数学模型

**符号定义**：
- $x_0$：原始发送比特（$b\in\{0,1\}$）对应的 BPSK 符号，$x_0 = -(2b-1)$
- $\sigma$：AWGN 标准差
- $k_{\text{theo}}$（theo\_FBC）：在当前 SNR 下理论期望硬错误数，$k_{\text{theo}}=\lfloor N\cdot rber\rfloor$，其中 $N$ 为本仿真实际送入信道的比特数（`blk_len`）；在 `ssd_fc_flip_eq` 中对应 `dsp_blk_len = H_N - dsp_pad_len`
- `target_fbc`：希望固定的**总**FBC（相对于 $x_0$），仅 `ssd_fc_flip_eq` 使用（来自 `.cnfg` 的“Target total FBC”一行或环境变量 `TARGET_FBC`）
- $k$：本次 TAWGN 要固定的硬错误数（相对于 $x_0$，输出为 `raw_err_num` / `raw_err_awgn`）

##### 3.5.2.1 固定 FBC（Fix FBC）如何做到（统一策略）

`target_fbc` 只决定“本次要固定的硬错误数 $k$ 取多少”，算法本身完全一致：

- 未设置 `target_fbc` 或 `target_fbc==0`：取 $k=k_{\text{theo}}=\lfloor N\cdot rber\rfloor$
- `target_fbc>0`：取 $k=\texttt{target\_fbc}$（并夹到 $[0,N]$）

因此当 `target_fbc < k_theo` 或 `target_fbc > k_theo` 时，并不需要不同策略：只是你选择了不同的条件事件 $E=k$。  
直观上：`target_fbc > k_theo` 对应 AWGN 尾部更“差”的帧；`target_fbc < k_theo` 对应更“好”的帧；两者都保持 AWGN 的**条件分布一致性**（AWGN-consistent conditioned on $E=k$）。

**信道传输过程**（`TRUNC_AWGN` 模式，纯 TAWGN on $x_0$）：

1. **设定目标错误数 $k$**：按上面的 `target_fbc` 规则取值
2. **选择错误位置集合 $S$**：在 $\{1,\dots,N\}$ 中均匀随机选择 $|S|=k$
3. **逐符号截断高斯采样**：
   - 若 $i\in S$：强制硬判决错误，从“错误侧截断高斯”采样
   - 若 $i\notin S$：强制硬判决正确，从“正确侧截断高斯”采样
4. **得到严格固定的错误数**：由构造可得 $d_H(x_0,\hat{x}) = k$（因此 `raw_err_num == raw_err_awgn == k`）

**实现细节**（`src/transceiver_flip_eq.cpp` 的 `TRUNC_AWGN` 分支）：

```c
// 伪代码示意：构造 err_mask（集合 S，大小为 k）
randomBinError(err_mask, err_pos_trunc, blk_len, k);

for (int i = 0; i < blk_len; i++)
{
    int want_error = err_mask[i]; // i ∈ S ?
    rx_blk[i] = sample_trunc_awgn(tx_blk[i], awgn_sigma, want_error);

    int hard_dec = (rx_blk[i] >= 0.0f) ? 0 : 1;
    if (hard_dec != tx_blk[i]) err_cnt++; // 相对于 x0 的错误（恒等于 k）
}
```

其中 `sample_trunc_awgn()` 函数通过**拒绝采样**实现截断高斯分布：

```c
float sample_trunc_awgn(char tx_bit, float sigma, int want_error)
{
    float noiseless = -(tx_bit * 2.0f - 1.0f);  // BPSK映射
    while (1)
    {
        float n = sigma * rand_gaussian();
        float y = noiseless + n;
        int hard_dec = (y >= 0.0f) ? 0 : 1;
        int is_error = (hard_dec != tx_bit);
        if ((want_error && is_error) || (!want_error && !is_error))
            return y;  // 接受样本
        // 否则拒绝，重新采样
    }
}
```

注：当前实现是“逐符号拒绝采样”，因此当 SNR 很高且你强制 `target_fbc>0`（对应极小的 `rber`）时，生成“错误侧截断样本”的期望重采样次数约为 $1/rber$，仿真会明显变慢。

#### 3.5.3 统计量输出

FLIP_EQ 版本提供详细的**电压bin分布统计**（假设8-bin量化，bin编号0-7）：

| 统计量 | 含义 |
|--------|------|
| `tx0_bin7` | 发送0时落入bin7（最高置信1区域，即HRE-like错误）的次数 |
| `tx0_bin5` | 发送0时落入bin5的次数 |
| `tx0_bin3` | 发送0时落入bin3的次数 |
| `tx0_bin1` | 发送0时落入bin1的次数 |
| `tx0_self` | 发送0时落入"自身bin"（bin1/3/5/7，即正确侧）的总次数 |
| `tx1_bin0` | 发送1时落入bin0（最高置信0区域，即HRE-like错误）的次数 |
| `tx1_bin2` | 发送1时落入bin2的次数 |
| `tx1_bin4` | 发送1时落入bin4的次数 |
| `tx1_bin6` | 发送1时落入bin6的次数 |
| `tx1_self` | 发送1时落入"自身bin"（bin0/2/4/6，即正确侧）的总次数 |

**输出示例**（`src/SSD_FC_flip_eq.cpp:426-430`）：

```
[STATISTICS] BINSTAT0: 0->7(HRE) 12.34 | 0->5 23.45 | 0->3 34.56 | 0->1 45.67 | 0->self(1/3/5/7) 103.68
[STATISTICS] BINSTAT1: 1->0(HRE) 56.78 | 1->6 67.89 | 1->4 78.90 | 1->2 89.01 | 1->self(1/3/5/7) 292.58
```

这些统计量可用于：
- 分析阈值电压分布的非对称性
- 评估不同电压bin的错误贡献
- 验证量化方案的合理性

#### 3.5.4 与标准TAWGN的区别

| 特性 | `transceiver.cpp` 的 TAWGN | `transceiver_flip_eq.cpp` 的 TAWGN |
|------|---------------------------|-----------------------------------|
| 目标 $k$ | 固定为 $k_{\text{theo}}=\lfloor N\cdot rber\rfloor$ | 可选：`target_fbc==0` 用 $k_{\text{theo}}$；`target_fbc>0` 用 $k=\texttt{target\_fbc}$ |
| HRE 注入 | **不叠加 HRE**（明确说明） | **不叠加 HRE**（`TRUNC_AWGN` 分支内忽略 `hre_bit/hre_vec`） |
| 错误数统计 | $d_H(x_0, \hat{x}) = k_{\text{theo}}$（恒定） | $d_H(x_0, \hat{x}) = k$（恒定） |
| 应用场景 | 固定FBC条件下的对称信道仿真 | 固定FBC条件仿真 + 允许扫不同的固定 $k$ 值（并输出 bin 分解统计） |

---

## 4. HRE 注入语义对照（AWGN 下最关键）

下面只讨论 AWGN 分支（因为当前 HRE 仅在 AWGN 分支生效）。

### 4.1 `ssd_fc_hre`（旧式 HRE，来自 `transceiver.cpp` 的非 `HRE2_MODE` 分支）

流程（概念上）：

1. 先按 AWGN 生成 $y=x+\sigma n$；
2. 再随机选 `hre_bit` 个位置（`hre_vec=1`）对接收值做覆盖写：
   - `hre_mode==0`：强制 `rx_blk=+1.0`（倾向落入“高置信 0”区域）
   - `hre_mode==1`：强制 `rx_blk=-1.0`（倾向落入“高置信 1”区域）
   - 其他：随机写入 `±1.0`

**重要陷阱：`raw_err_num` 的计数点在覆盖写之前。**

也就是说，这个版本的 `raw_err_num` 更接近“AWGN 硬判决误错数”，而不是“叠加 HRE 后的真实误错数”。若你用 `raw_err_num` 做筛选/分桶，会与“注入后真正送入译码器的错误分布”产生偏差。

### 4.2 `ssd_fc_hre2`（新式 HRE2，来自 `transceiver.cpp` 的 `HRE2_MODE` 分支）

该版本把 `hre_bit` 明确解释为“**确定数量的符号翻转（bit-flip）**”，并提供 2 种注入模式（由 `hre_mode` 控制）：

- `hre_mode==0`：选中的 HRE 位 **翻转且不加噪声**  
  $rx=-x$，因此这些位在硬判决意义上几乎“必错”，且置信度较高（取决于量化/Vref）
- `hre_mode!=0`：选中的 HRE 位 **翻转后仍叠加 AWGN**  
  $rx=-x+\sigma n$，错误方向被固定但置信度有噪声波动，极小概率被噪声拉回

同时它提供更细的统计量：

- `raw_err_num`：按最终 `rx_blk` 硬判决得到的**总错误数**（包含 HRE 注入）
- `raw_err_awgn`：仅统计 `hre_vec==0` 的位置的硬判决错误数（意图近似“纯 AWGN 贡献”）
- `fixed_hre_cnt`：当前代码直接赋值为 `hre_bit`（**注意**：除 AWGN 外其他信道分支并未真正使用 `hre_vec`）
- `hre_like_cnt`：见第 5 节（HRE-like 定义）

### 4.3 `ssd_fc_hre_only`（HRE-only transceiver 的注入语义）

`transceiver_hre_only.cpp` 的 AWGN+HRE 注入目前固定为：

- `hre_vec==1`：$rx=-x$（翻转且不加噪声）
- `hre_vec==0`：$rx=x+\sigma n$

---

## 5. 统计量的“定义差异”与 HRE-like 的严格含义

这一节专门回答：同名统计量在不同版本里究竟代表什么、是否可比。

### 5.1 `raw_err_num` 与 `raw_err_awgn`

- 在 `ssd_fc_hre`（旧式 HRE）里：
  - `raw_err_num`：在 HRE 覆盖写之前统计，因此更接近“纯 AWGN 的硬判决误错数”
  - `raw_err_awgn`：当前实现为 0（不可用）
- 在 `ssd_fc_hre2`（新式 HRE2）里：
  - `raw_err_num`：按最终 `rx_blk` 硬判决统计（包含 HRE 注入）
  - `raw_err_awgn`：只统计 `hre_vec==0` 的位置，近似分离“纯 AWGN 贡献”
- 在 `ssd_fc_hre_only` 里（`transceiver_hre_only.cpp`）：
  - `raw_err_num`：按最终 `rx_blk` 硬判决统计（包含 HRE 注入）
  - `raw_err_awgn`：只统计 `hre_vec==0` 的位置

因此：若你要在不同版本之间比较 “FBC(no HRE)” 或 “AWGN 部分错误数”，请优先使用 `ssd_fc_hre2` 的 `raw_err_awgn`，不要直接拿 `ssd_fc_hre` 的 `raw_err_num` 当作“注入后总 FBC”。

### 5.0 `.cnfg` 里三项 HRE 配置在不同版本的含义（必须显式区分）

配置文件中连续三行通常为：

- `hre_bit`（HRE bit count）
- `hre_mode`（HRE bit configure mode）
- `hre_dec`（HRE llr process enable）

其版本相关含义如下：

| 字段 | `ssd_fc_hre`（旧式 HRE） | `ssd_fc_hre2`（新式 HRE2） | `ssd_fc_hre_only`（hre_only transceiver） |
|---|---|---|---|
| `hre_bit` | 每帧随机标记 `hre_bit` 个位置用于覆盖写 | 每帧随机标记 `hre_bit` 个位置用于 bit-flip 注入 | 每帧随机标记 `hre_bit` 个位置用于 bit-flip 注入 |
| `hre_mode` | 0→强制 `rx=+1`；1→强制 `rx=-1`；其他→随机 `±1` | 0→翻转且不加噪；非0→翻转后仍加 AWGN | **当前实现忽略该字段**（固定为“翻转且不加噪”） |
| `hre_dec` | 传入 LDPC 译码端：若启用，会对标记位做 LLR 降权/置零（见下） | 同左；但会与“固定 bit-flip”语义交织（需谨慎） | 无译码器，因此该字段在该可执行文件中不生效 |

译码端的处理（仅 `ssd_fc_hre/ssd_fc_hre2` 有意义）：

- `hre_dec==1` 时，`ldpc_codec.cpp` 的层译码初始化会把被标记位置的 LLR 置为 0（等价于“把这些位当作擦除/不可信”）。
- 这是一种**缓解策略**：它并不是“注入高置信度”，而是恰恰相反——主动移除高置信信息，避免 HRE 以高置信错误污染迭代。

### 5.2 `hre_like_cnt`（“概率型 HRE-like”并不自动成立）

`hre_like_cnt` 的设计目标是统计：

- $b=1$ 却落入“最高置信 0”的量化 bin
- $b=0$ 却落入“最高置信 1”的量化 bin

即“错且极可信”的事件（HRE-like）。

在 `src/transceiver.cpp`（`HRE2_MODE`）中：

- “最高置信 0 bin” 用 `hi_conf0_bin` 表示（通常对应最大正 LLR 的 bin）
- “最高置信 1 bin” 用 `hi_conf1_bin` 表示（通常对应最大负 LLR 的 bin）
- 判定条件（只在非 DIRECT 量化下统计）：
  - `tx==1 && det==hi_conf0_bin` 或 `tx==0 && det==hi_conf1_bin`

在 `src/transceiver_hre_only.cpp` 中：

- 直接写死为 `tx==1 && det==0` 或 `tx==0 && det==7`
- 这隐含假设“bin0/bin7 恰好是两端最极端置信 bin”，对某些 Vref/映射成立，但不保证泛化。

**关键注意：当前实现并未排除 `hre_vec==1` 的位置。**  
也就是说，只要你注入了固定翻转型 HRE，这些位很可能也会被 `hre_like_cnt` 计入，从而导致：

- `hre_like_cnt` 不是“纯概率型 HRE”计数
- 它更准确的定义是：“落入极端错误 bin 的总次数（可能包含固定注入与 AWGN 极端事件）”

如果你严格想分离：

- 固定型 HRE：由 `hre_vec` 注入（计数应为 `fixed_hre_cnt`）
- 概率型 HRE-like：由 AWGN/量化自然产生

那么实现上需要把 `hre_like_cnt` 的统计限定在 `hre_vec==0` 的位置上（当前代码未做）。

### 5.3 `rber` 的有效范围

- 在 `src/transceiver.cpp` 里，`rber` 仅在 `CLEAN/AWGN/TAWGN` 分支中赋值
- 在 `BSC/ERR_INJ/MAX_ERR` 中未显式赋值（因此上层若依赖 `rber`，需谨慎）

### 5.4 不同版本的统计量输出对比

下表系统性地对比四个可执行程序输出的统计量及其含义差异：

| 统计量 | `ssd_fc_hre` | `ssd_fc_hre2` | `ssd_fc_hre_only` | `ssd_fc_flip_eq` |
|--------|-------------|--------------|------------------|-----------------|
| **`raw_err_num`** | HRE注入**前**的AWGN硬判决误错数（在覆盖写之前统计） | 注入**后**的总硬判决误错数（包含HRE注入） | 注入**后**的总硬判决误错数（包含HRE注入） | `TRUNC_AWGN`: 相对于 $x_0$ 的总误错数，恒等于固定 $k$（`target_fbc==0` 用 $k_{\text{theo}}$，否则用 `target_fbc`） |
| **`raw_err_awgn`** | 不可用（恒为0） | 仅 `hre_vec==0` 位置的硬判决误错数（近似"纯AWGN贡献"） | 仅 `hre_vec==0` 位置的硬判决误错数 | `TRUNC_AWGN`: 恒等于 $k$（与 `raw_err_num` 相同）；`AWGN`: 不可用（恒为0） |
| **`fixed_hre_cnt`** | 不输出 | 等于 `hre_bit`（固定HRE注入数） | 不输出 | `TRUNC_AWGN`: 恒为 0（忽略 `hre_bit`）；`AWGN`: 等于 `hre_bit` |
| **`hre_like_cnt`** | 不输出 | 落入极端错误bin的总次数（**包含**固定HRE位置） | 落入bin0/bin7的总次数（写死判定条件） | 输出：落入“错误侧极端 bin”（高置信错）的总次数（同时还输出更细的 bin 统计） |
| **`dec_err_num`** | 译码后残留误码数 | 译码后残留误码数 | 不输出（无译码器） | 译码后残留误码数 |
| **`cw_fail`** | 译码失败帧数（FER统计） | 译码失败帧数（FER统计） | 不输出（无译码器） | 译码失败帧数（FER统计） |
| **电压bin统计** | 无 | 无 | 无 | **详细输出**：`tx0_bin7/5/3/1`、`tx1_bin0/2/4/6`、`tx0_self`、`tx1_self` |
| **迭代次数统计** | `fdec_itr`、`ldec_itr` | `fdec_itr`、`ldec_itr` | 不输出（无译码器） | `fdec_itr`、`ldec_itr` |

**关键注意事项**：

1. **`raw_err_num` 的语义差异最大**：
   - `ssd_fc_hre`：统计点在HRE覆盖写**之前**，因此更接近"纯AWGN误错数"
   - 其他版本：统计点在HRE注入**之后**，反映"送入译码器的真实误错数"
   - **跨版本对比时务必注意这一差异**，否则会得出错误结论

2. **`raw_err_awgn` 的可用性**：
   - 仅 `ssd_fc_hre2`、`ssd_fc_hre_only`、`ssd_fc_flip_eq` 提供有效值
   - 若要分析"纯AWGN贡献"，应使用这些版本的 `raw_err_awgn`，而非 `ssd_fc_hre` 的 `raw_err_num`

3. **`hre_like_cnt` 的混淆性**：
   - 当前实现**未排除** `hre_vec==1` 的位置，因此该统计量**不是**"纯概率型HRE-like"
   - 它更准确的定义是："落入极端错误bin的总次数（固定HRE + AWGN极端事件）"
   - 若需严格分离固定HRE与概率HRE-like，需修改代码限定统计范围

4. **电压bin统计的独特性**：
   - 仅 `ssd_fc_flip_eq` 提供详细的bin级统计
   - 可用于分析阈值电压分布、验证量化方案、评估非对称性

5. **FBC（Flipped Bit Count）的计算**：
   - `ssd_fc_hre_only` 输出 `FBC(bit)` 和 `FBC(no HRE)`，分别对应 `raw_err_num` 和 `raw_err_awgn` 的平均值
   - 其他版本需自行从 `raw_err_num` 计算平均FBC

---

## 6. SSD\_FC 主链路里的“样本筛选仿真”（容易忽略但影响巨大）

`src/SSD_FC.cpp` 与 `src/SSD_FC_hre2.cpp` 都包含如下逻辑（位于 `ch_detector()` 之后、`ecc_decoder()` 之前）：

```c
if (sim_pckt->raw_err_num < (int)(dsp_blk_len * sim_pckt->rber))
{
    sim_cnt--;
    continue;
}
```

这意味着：

- 你并没有在做“标准 AWGN 的无条件 FER”，而是在做一种**带筛选的条件仿真**：只保留“硬错误数不少于理论期望值 $\lfloor N\cdot rber\rfloor$”的帧
- 该筛选会改变误差分布（倾向保留更糟糕的帧），其结果更接近一种“偏保守”的性能评估，而不是通信意义上的真实 FER

并且它在不同版本里筛选对象不同：

- `ssd_fc_hre`：`raw_err_num` 更接近“AWGN 误错数（在 HRE 覆盖写之前统计）”，因此筛选条件主要约束 AWGN 部分
- `ssd_fc_hre2`：`raw_err_num` 是“注入后总误错数”，因此筛选条件也会把 HRE 注入纳入
- 在 `TAWGN` 模式下：由于信道端已构造 `raw_err_num = \lfloor N\cdot rber\rfloor`，该筛选基本恒通过，因此不会再丢帧

如果你的目标是：

- **真实 AWGN FER 曲线**：这段筛选应当关闭（否则得到的是偏置结果）
- **固定/近固定 FBC 条件下的性能**：建议使用 `TAWGN`（见 `doc/truncated_AWGN.md`），而不是依赖逐帧筛选

---

## 7. 已知不一致与使用建议（基于当前代码状态）

### 7.1 `ssd_fc_hre_only` 当前无法重新编译（TRUNC\_AWGN 枚举不一致）

现状：

- `src/SSD_FC_hre_only.cpp` 已加入对 `TRUNC_AWGN/TAWGN` 的解析与配置打印
- 但 `src/transceiver_hre_only.h` 的 `enum ch_model` **不包含** `TRUNC_AWGN`

因此目前重新 `make ssd_fc_hre_only` 会报 `TRUNC_AWGN was not declared`。  
结论：在未修复前，`ssd_fc_hre_only` 实际仅支持 `CLEAN/AWGN/BSC/ERR_INJ/MAX_ERR`。

### 7.2 选择哪个可执行程序更“语义清晰”

- 你要“全链路 + 旧式 HRE（强制写 ±1）”：用 `ssd_fc_hre`
- 你要“全链路 + 固定数量 bit-flip（两种翻转模式）+ 更完整统计”：用 `ssd_fc_hre2`
- 你要“只看信道/量化统计，不跑译码”：用 `ssd_fc_hre_only`（但当前需先修复 7.1 才能继续扩展到 TAWGN）
