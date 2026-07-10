# auto_fer_eval 优化分析与修复方案

## 一、拟合点位过深 -- 根因分析

### 现象

`fit_target_fer` 设为 `1e-6`，实际扫描点位深至 `1e-12` 甚至更低。

### 执行链路

```
fit_target_fer 来源
  ├─ TOML 未配 → 默认取 fer_lo (L484)
  ├─ loop 覆盖 fer_lo 但未覆盖 fit_target_fer → 跟 loop.fer_lo (L539-541)
  └─ TOML 或 loop 显式配置 → 直接使用

拟合外推 (predict_stop_axis_from_manifest, L3930-4050)
  ├─ 筛选 FER_eff ∈ [fit_fer_lo, fit_fer_hi]  (默认 [1e-4, 1e-2])
  ├─ _linfit: log10(FER_eff) = a * log10(RAW_BER) + b  (L3993)
  ├─ 外推: y_target = log10(fit_target_fer), x_target = (y-b)/a  (L4000-4001)
  ├─ 反推 SNR/K → est_stop
  └─ deep_stop = min(x_stop, est_stop)  (L1009)

deep 点位生成 (adaptive_deep_scan, L1903-1934)
  ├─ 从 deep_start 到 deep_stop 按 step_low 一次性全部生成
  └─ 并发提交 (max_in_flight 限制)

硬停止 (L2137-2164)
  ├─ 每完成一个点检查 fer_eff <= fer_lo
  ├─ 满足 → cancel 剩余 in-flight → return
  └─ fer_effective(fer=0, N) = 3/N，FER=0 不会立即触发

completion (L3031-3304)
  └─ 扫描所有已有点位，不区分深度，不完整就重提
```

### 根因

| # | 原因 | 位置 | 严重度 |
|---|------|------|--------|
| **R1** | **拟合杠杆效应** -- 拟合窗口仅 2 decade (`1e-4`~`1e-2`)，外推到 `1e-6` 需要跨 2+ decade。waterfall 深区斜率通常比浅区更陡（曲线凸性），浅区拟合斜率低估深区实际斜率 → est_stop 远超实际需要 | L3982, L3993 | **主因** |
| **R2** | **`_linfit` 无约束** -- 无斜率方向/范围约束、无 R² 检验、无外推距离限制 | L3830-3847 | **主因** |
| **R3** | **adaptive 不用 `fit_max_extend`** -- README 明确说是旧接口。sequential `_run_case` 用了 (L748-764)，adaptive 路径完全没有 | L1000-1011 | **主因** |
| **R4** | **并发预提交过冲** -- deep_stop 本身已经过深，再加上 max_in_flight 个点的并发窗口，达到 fer_lo 时已有多个更深点在跑 | L1992, L2137 | 放大器 |
| **R5** | **completion 不过滤深度** -- `adaptive_finalize_incomplete_logs` 重提所有不完整点，包括 deep 过冲的超深点残留 | L3063-3094 | 放大器 |
| **R6** | **历史残留** -- 同一输出目录的历史 jobs.json/manifest 中的超深点被 completion 捡回来 | L3061 | 放大器 |

### 问题链条

```
拟合窗口浅 (1e-4~1e-2, 2 decade)
  → 斜率 a 低估深区实际斜率 (waterfall 凸性)
  → est_stop 远超实际需要 (可能偏差 1~2 dB)
  → deep_stop 过深 (无 fit_max_extend 限制)
  → 一次性生成大量超深点并发提交
  → fer_lo 硬停有 max_in_flight 个点的并发过冲
  → completion 把过冲点的残留日志补全
  → 最终点位深度远超 fit_target_fer
```

---

## 二、拟合过深修复方案

> **结论：只需 Fix 1 + Fix 2 即可解决核心问题。Fix 3/4/5 在 Fix 1+2 到位后边际收益可忽略。**

### Fix 1: adaptive 路径加入 `fit_max_extend` 限制 (**必做**)

**位置:** L1006-1011

**当前:**
```python
if est_stop is not None:
    if main_case.direction > 0:
        deep_stop = min(float(main_case.x_stop), float(est_stop))
    else:
        deep_stop = max(float(main_case.x_stop), float(est_stop))
```

**改为:**
```python
if est_stop is not None:
    if main_case.direction > 0:
        deep_stop = min(float(main_case.x_stop), float(est_stop))
        if main_case.fit_max_extend > 0 and deep_start is not None:
            deep_stop = min(deep_stop, float(deep_start) + main_case.fit_max_extend)
    else:
        deep_stop = max(float(main_case.x_stop), float(est_stop))
        if main_case.fit_max_extend > 0 and deep_start is not None:
            deep_stop = max(deep_stop, float(deep_start) - main_case.fit_max_extend)
```

**同时修改默认值** (`_parse_case` L487):
```python
fit_max_extend = float(d.get("fit_max_extend", 0.5))  # 从 1.0 改为 0.5
```

**效果:** 从 trigger 点到 deep_stop 最多 0.5 dB。即使拟合外推到 3 dB 远，实际只扫 trigger + 0.5 dB。
以 step_low=0.025 计，0.5 dB 内最多 20 个点，不会浪费资源。

### Fix 2: `_linfit` 增加拟合质量守卫 (**必做**)

**位置:** L3993-3998

**在拟合后加入:**
```python
a, b = fit1

# Guard 1: 斜率方向 -- waterfall 中 FER 应随 RAW_BER 同向变化 (a > 0)
if a <= 0:
    return None

# Guard 2: R² 检验 -- 外推置信度
ss_res = sum((y - (a * x + b)) ** 2 for x, y in zip(xs, ys))
ss_tot = sum((y - sum(ys) / len(ys)) ** 2 for y in ys)
r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
FIT_R_SQ_MIN = 0.90
if r_sq < FIT_R_SQ_MIN:
    return None

# Guard 3: 外推距离限制 -- 不允许外推超过拟合窗口跨度的 N 倍
FIT_EXTRAP_RATIO_MAX = 3.0
fit_span = max(ys) - min(ys)  # 拟合窗口内的 y 跨度 (decade)
extrap_dist = abs(y_target - min(ys))
if fit_span > 0 and extrap_dist / fit_span > FIT_EXTRAP_RATIO_MAX:
    return None
```

**效果:**
- 斜率方向错误 → 直接放弃拟合
- 拟合质量差 (R² < 0.9) → 放弃外推
- 外推距离超过拟合窗口 3 倍 → 放弃（2 decade 窗口最多外推到 6 decade 深度）

### Fix 3/4/5: 边际收益分析 (可不做)

| Fix | 解决什么 | 为何 Fix 1+2 后不需要 |
|-----|---------|----------------------|
| Fix 3 渐进式提交 | fer_lo 触发时减少并发过冲 | Fix 1 限制总跨度 0.5 dB，即使全部预提交也只有 20 个 0.025 dB 间隔的点 |
| Fix 4 completion 过滤 | 历史残留超深点被重提 | 新 case 不再产生超深点，无残留；旧数据一次性手动清理即可 |
| Fix 5 fer_effective 边界 | fer=0,N=None 返回 0.0 | 拟合函数 L3980 已有 `fer_eff<=0` 守卫，不会进入拟合 |

---

## 三、其他优化建议

### 3.1 可靠性 (HIGH)

| # | 问题 | 位置 | 建议 |
|---|------|------|------|
| **H1** | manifest 仅在阶段退出时持久化，崩溃丢失全部运行时进度 | `save_manifest` 仅在 L702/784/905/1076 | poll 循环内每完成一点或周期性刷盘 |
| **H2** | JobDB 无文件锁，并发写入可损坏 jobs.json | `job_db.py` L95-100 | atomic write (先写 tmp 再 rename) + flock |
| **H3** | 固定轮询间隔，无自适应退避 | L1838/2172/2763/3288/3716 | exponential backoff (5s→10s→20s→...→max) |
| **H4** | `_in_flight_has_axis` O(n) 全量扫描 | L2846-2848 | 统一 key 后用 dict O(1) 查找 |

### 3.2 性能 (MEDIUM)

| # | 问题 | 位置 | 建议 |
|---|------|------|------|
| **M1** | `_job_identity` 闭包重复定义 5 次 | L1402/1896/2603/3121/3539 | 提取为模块级函数 |
| **M2** | `pending.pop(0)` 在 list 上 O(n) | L1993/2644/3162/3642 | 改用 `collections.deque.popleft()` |
| **M3** | `find_run` 线性扫描 manifest runs | `manifest.py` L52-56 | 维护 `{(axis_type, axis_key): index}` 索引 |
| **M4** | log 文件在同一决策路径中被读取 2-3 次 | `_log_has_complete_statistics` + `_parse_done_log` | 让检查函数接受已读取的 text 参数 |
| **M5** | LSF bjobs 逐 job 调用 | `executors.py` L371-434 | 添加 `poll_batch()` 一次查多个 job |
| **M6** | 函数体内冗余 `import re` | L623, L2491 | 删除（顶部 L7 已导入） |

### 3.3 代码结构 (LOW，不影响功能)

| # | 问题 | 建议 |
|---|------|------|
| **L1** | `auto_fer_eval.py` 4458 行 | 按职责拆分: scan/deep/backfill/finalize/fitting/export/config/utils |
| **L2** | `adaptive_main_scan` 约 500 行 | 拆分为类或若干专职函数 |
| **L3** | deep/backfill/finalize/completion 四处高度重复的 submit-poll-parse 循环 | 抽象为通用 `parallel_submit_and_poll()` 调度器 |
| **L4** | 配置解析无合理性校验 | `_parse_adaptive` 中添加值域检查 |
| **L5** | `_log_has_complete_statistics` 与 `LsfExecutor._log_has_final_statistics` 功能重复 | 统一实现 |
| **L6** | dead code: `predict_stop_snr_from_manifest` (L4053), `_record_from_manifest` (L4241) | 确认无外部使用后删除 |
| **L7** | 测试覆盖率极低 (仅 2 个测试文件，4 个用例) | 补充 fitting/trigger/timeout 等核心逻辑的单元测试 |

### 3.4 auto_throughput_eval 共享问题

| # | 问题 | 建议 |
|---|------|------|
| **T1** | `fail_fast=false` 时无最大重试次数，失败点无限循环重提 | 加 `max_retries` 参数 |
| **T2** | `run_case_local` 与 `run_case_lsf` 大量重复样板代码 | 与 auto_fer_eval 共享通用 submit-poll 框架 |

---

## 四、建议实施顺序

```
Phase 1 (解决拟合过深 — 必做):
  Fix 1 — adaptive 加 fit_max_extend 限制 (默认 0.5 dB)
  Fix 2 — _linfit 加质量守卫 (斜率方向 + R² + 外推距离)

Phase 2 (可靠性 — 按需):
  H1 — manifest 周期性刷盘
  H2 — atomic write + flock

Phase 3 (性能 — 按需):
  M2 — deque 替代 list.pop(0)
  M5 — bjobs 批量查询
  H4 — in_flight O(1) 查找

Phase 4 (代码结构 — 可选):
  L1~L3 — 拆分大文件、抽象通用调度
```
