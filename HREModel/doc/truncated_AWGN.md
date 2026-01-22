# 截断 AWGN 信道（Truncated-AWGN）用于固定 FBC 仿真

本说明记录如何在当前 C 模型中，用“截断高斯分布”构造一个新的信道分支，用于在 AWGN 场景下高效仿真“硬判决错误数固定为某值 $k$”时 LDPC 的纠错性能，从而替代原先低效的逐帧拒绝采样。

## 1. 背景：原有固定 FBC 仿真方式

在现有 C 模型中，固定 FBC 的典型做法是：

- 按标准 AWGN 信道生成一帧软数据；
- 经过硬判决得到该帧的错误比特数 $E$（原始 FBC）；
- 只保留满足 $E = k = \lfloor N\cdot\mathrm{RBER}\rfloor$ 的帧，其它帧全部丢弃（`continue`）。

这实际上是在估计条件帧误码率
$P(\text{fail}\mid E=k)$，而不是整体 FER。问题在于：

- 在 $E\sim\mathrm{Binomial}(N,p)$ 下，事件 $E=k$ 的概率往往不高；
- 大量模拟的帧被直接丢弃，仿真效率很低。

目标：在不改变统计意义（仍然估计 $P(\text{fail}\mid E=k)$）的前提下，构造一个新的信道分支，使得生成的每一帧天然满足 $E=k$，且符号幅度分布仍然符合 AWGN 模型的条件分布。

## 2. BPSK+AWGN 模型与条件分布

假设标准 BPSK+AWGN：

- 比特 $b_i \in \{0,1\}$，映射为 $x_i = - (2b_i - 1)$；于是 $b=0\to x=+1$，$b=1\to x=-1$；
- 噪声 $n_i \sim \mathcal{N}(0,\sigma^2)$，相互独立；
- 接收符号 $y_i = x_i + n_i$；
- 硬判决 $\hat b_i = \mathbf{1}\{y_i<0\}$（阈值 0）；
- 错误指示 $E_i = 1\{\hat b_i\ne b_i\}$，总错误数 $E = \sum_{i=1}^N E_i$。

我们希望直接从条件分布
$\{y_i\}_{i=1}^N \mid E = k$ 中采样，其中 $k$ 是指定的错误比特数（例如近似为 $N\cdot\mathrm{RBER}$）。

### 2.1 正确侧与错误侧：截断高斯视角

以 $b=0\Rightarrow x=+1$ 为例：

- 正确事件：$y>0$。此时 $y$ 服从 $\mathcal{N}(+1,\sigma^2)$ 在 $(0,+\infty)$ 上的截断高斯分布；
- 错误事件：$y<0$。此时 $y$ 服从 $\mathcal{N}(+1,\sigma^2)$ 在 $(-\infty,0)$ 上的截断高斯分布。

对于 $b=1\Rightarrow x=-1$：

- 正确：$y<0$，即 $\mathcal{N}(-1,\sigma^2)$ 在 $(-\infty,0)$ 上截断；
- 错误：$y>0$，即 $\mathcal{N}(-1,\sigma^2)$ 在 $(0,+\infty)$ 上截断。

如果给定某个错误位置集合 $S\subset\{1,\dots,N\}$，且 $|S| = k$（即这些位置上 $E_i=1$），在 AWGN + 独立假设下：

- 当 $i\in S$ 时，$y_i$ 是“错误侧”的截断高斯；
- 当 $i\notin S$ 时，$y_i$ 是“正确侧”的截断高斯；
- 集合 $S$ 在所有大小为 $k$ 的子集上均匀分布。

因此，只要：

1. 在 $\{1,\dots,N\}$ 中均匀随机选择大小为 $k$ 的集合 $S$；
2. 对 $i\in S$，从错误侧截断高斯中独立采样 $y_i$；
3. 对 $i\notin S$，从正确侧截断高斯中独立采样 $y_i$；

便可得到与“在 AWGN 信道下条件于 $E=k$”完全等价的 $\{y_i\}$ 样本。

## 3. 在 `ch_transmit` 中的 C 级实现思路

在实现上，我们不必显式写解析截断高斯采样器，而可以用“单符号级拒绝采样”来实现：

- 仍按 $y = x + n$ 生成候选样本；
- 若当前期望“正确”，就只接受硬判决正确的 $y$；若期望“错误”，就只接受硬判决错误的 $y$；
- 由于这是在“符号级”进行的拒绝，开销远小于逐帧拒绝。

### 3.1 单符号截断 AWGN 采样函数（伪代码）

```c
static float sample_trunc_awgn(int tx_bit, float sigma, int want_error)
{
    // tx_bit: 0 或 1
    // want_error: 1 -> 希望硬判决错误；0 -> 希望硬判决正确
    float noiseless = -1.0f * (tx_bit * 2.0f - 1.0f); // 0->+1, 1->-1
    while (1) {
        float n = sigma * rand_gaussian();
        float y = noiseless + n;
        int hard_dec = (y >= 0.0f) ? 0 : 1;
        int is_error = (hard_dec != tx_bit);
        if ((want_error && is_error) || (!want_error && !is_error)) {
            return y;
        }
    }
}
```

### 3.2 新增信道分支（示意，不替换原 AWGN）

在 `enum ch_model` 中增加一个新模式（示例）：

- `TRUNC_AWGN`：截断 AWGN，固定硬错误数的信道。

在 `ch_config` 中：

- 与 AWGN 模式相同，根据 SNR 设置 `awgn_sigma`；
- 根据理论公式计算 RBER，用于决定目标错误数 $k$；
- 打印提示该模式用于固定 FBC 仿真。

在 `ch_transmit` 中新增分支（伪代码）：

```c
else if (ch_sel == TRUNC_AWGN)
{
    int k = (int)(blk_len * rber);  // 目标硬错误数

    int *err_mask = (int*)calloc(blk_len, sizeof(*err_mask));
    int *err_pos  = (int*)calloc(k, sizeof(*err_pos));

    // 均匀选择 k 个错误位置
    randomBinError(err_mask, err_pos, blk_len, k);

    for (int i = 0; i < blk_len; i++)
    {
        int want_error = err_mask[i];  // 1: 强制错误, 0: 强制正确
        rx_blk[i] = sample_trunc_awgn(tx_blk[i], awgn_sigma, want_error);
    }

    raw_err_num  = k;
    raw_err_awgn = k;

    free(err_mask);
    free(err_pos);
}
```

注意：

- 原有 `AWGN` 分支保持不变，继续用于标准 FER 仿真；
- `TRUNC_AWGN` 是额外的通道模式，只在需要固定 FBC 的条件 FER 仿真时启用；
- 上述代码仅为思路示例，实际整合需要考虑 HRE 注入统计等细节。

## 4. 与旧有 `if (raw_err_num == N*rber)` 过滤的关系

旧有写法类似：

```c
// 仅保留硬错误数等于期望值的帧
if (raw_err_num != (int)(blk_len * rber)) {
    sim_cnt--;
    continue;
}
```

- 它通过拒绝采样实现 $P(\text{fail}\mid E=k)$ 的估计；
- 但当 $P(E=k)$ 很小时，接受率很低，导致仿真极慢。

截断 AWGN 分支的方案：

- 直接在生成时保证 $E=k$，不丢弃任何帧；
- 每个符号的幅度分布仍然是 AWGN 在“正确/错误侧”的截断；
- 从概率论角度看，与“先跑 AWGN，再筛掉 $E\ne k$ 的帧”完全等价，但计算效率高得多。

因此，当目标是研究“固定硬错误数下的条件 FER”时，截断 AWGN 分支是对旧有 `if (raw_err_num == N*rber)` 过滤逻辑的数学等价、但计算上更高效的替代方案。

## 5. 适用范围与注意事项

- 截断 AWGN 信道适用于研究 $P(\text{fail}\mid E=k)$：即在硬错误数固定为 $k$ 时 LDPC 的鲁棒性；
- 它不适合作为标准 AWGN FER 的替代：真实 FER 必须对 $E$ 的全分布求和（包括均值附近和尾部）；
- 与 HRE 注入的组合要谨慎：
  - 若只关注 AWGN 固定 FBC，建议在 `TRUNC_AWGN` 模式下关闭 HRE（`hre_bit=0`）；
  - 若要同时考虑 HRE+AWGN，则需明确 $k$ 所代表的错误来源（纯 AWGN 还是 AWGN+HRE），以及如何在统计中分别记账。

