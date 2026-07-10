# QC192 Batch Simulation 脚本使用说明

本目录提供脚本：

- `batch_sim.py`
- `batch_sim_simple.py`

用于按尺寸批量仿真 LDPC 矩阵（目录结构类似 `IBEX/ibex_matrix`），流程如下：
 
1. 基于模板 `.cnfg` 按 `M` 自动生成配置（仅修改第 5 行 `parity size = qc/8 * M`）。
2. 每个尺寸先用一个代表矩阵（默认 `matrix_id=1`）调用 `auto_fer_eval` 确定点位。
3. 对同尺寸所有矩阵按同一组点位 `bsub` 提交仿真（不使用 `-cwd`）。
4. 按最深点位做 FER 排名。
5. 每个尺寸输出一张 `rber-fer` 图。

`batch_sim_simple.py` 是简化流程，适合你最新这套“快速选点+排名”需求：

1. 每个尺寸先调用 `auto_fer_eval adaptive`，从控制台解析 `pilot start`（即 `main_start`）后立刻停止。
2. 固定仿真点位：`[main_start, main_start+0.4]`（步长 `0.1`）并追加 `main_start+0.35`、`main_start+0.45`。
3. 排名点位固定为最深三个：`main_start+0.35` / `+0.4` / `+0.45`。
4. 批量阶段实时排名：按 `--bpeek-interval-sec` 周期（默认 `3600` 秒）用 `bpeek` 更新 `summary/realtime/*.csv`。

---

## 1. 前置条件

- Python 3.11+（脚本用到标准库，无额外强依赖）
- 可执行仿真程序（如 `IBEX/ssd_fc_src2`）
- 可用的 `bsub / bjobs / bkill` 命令（真实 LSF 或 shim）
- `scripts/auto_fer_eval/auto_fer_eval.py` 存在
- 矩阵目录结构类似：
  - `<matrix_root>/<M>x<N>/matrix/...`
  - `<matrix_root>/<M>x<N>/fade_matrix/...`
  - `<matrix_root>/<M>x<N>/occupied_matrix/...`

---

## 2. 快速开始（QC192）

```bash
python3 scripts/qc192_batch_sim/batch_sim.py \
  --template-cnfg <template.cnfg> \
  --matrix-root <qc192_matrix_root> \
  --exe <sim_exe> \
  --out-dir <output_dir> \
  --m-range 19-45 \
  --n-matrices 5 \
  --probe-id 1 \
  --qc 192 \
  --queue regr_q
```

---

## 3. 点位规则（已固化在脚本中）

从 probe 阶段 `auto_fer_eval` 的 `main/manifest.json` 读取 `(snr, ldpc_fer)`：

- 若存在 `FER ∈ [1e-5, 1e-4]` 的点：
  - 取该区间最深 SNR 为 `deepest`
  - 最终点位 = `main` 全部点位 + `deepest + 0.05`
- 若不存在：
  - 取 `main` 最深 SNR 为 `deepest`
  - 最终点位 = `main` 全部点位 + `deepest + 0.1`

---

## 4. 已完成判定与重提

正式批量阶段对每个 `(matrix_id, snr)` 的日志做判定：

- 日志存在且包含最终 `[STATISTICS] LDPC FER`：视为已完成
- 否则：视为未完成，重新 `bsub`

通过 `--max-resubmit-rounds` 控制最多重提轮数（默认 3）。

---

## 5. 作业命名与日志命名

脚本提交命名固定为：

- `-J`: `qc{qc}_M{M}_id{id}_snr{snrTag}`
- `-o`: `<out>/<M>x<N>/logs/qc{qc}_M{M}_id{id}_snr{snrTag}.log`

其中 `snrTag` 例如：

- `3.95 -> 3p95`
- `5.0 -> 5`

`bsub` 命令中不包含 `-cwd`。

---

## 6. 主要参数

### 必填参数

- `--template-cnfg`: 模板 `.cnfg`
- `--matrix-root`: 矩阵根目录
- `--exe`: 仿真可执行文件

### 任务范围

- `--m-range`: `M` 范围，支持 `19-45` 或 `19,20,21`（默认 `19-45`）
- `--n-matrices`: 每个尺寸矩阵数量 `n`（默认 `5`）
- `--probe-id`: probe 使用的矩阵 id（默认 `1`）
- `--qc`: QC 尺寸（默认 `192`）

### LSF 与调度

- `--queue`: LSF 队列（默认 `regr_q`）
- `--batch-max-in-flight`: 批量阶段并发上限（默认 `20`）
- `--batch-poll-sec`: 批量轮询间隔秒数（默认 `30`）
- `--max-resubmit-rounds`: 最大重提轮数（默认 `3`）

### probe（auto_fer_eval）相关

- `--auto-fer-eval`: `auto_fer_eval.py` 路径
- `--python-bin`: 调用 auto_fer_eval 的 python（默认当前解释器）
- `--snr-min / --snr-max`: probe 搜索区间
- `--main-step`: probe main 扫描步长（默认 `0.1`）
- `--pilot-max-sim`: probe pilot 最大仿真包数（默认 `100`）
- `--main-max-sim`: probe main 最大仿真包数（默认 `1000`）
- `--probe-max-in-flight`: probe 并发（默认 `20`）
- `--probe-poll-sec`: probe 轮询间隔（默认 `30`）
- `--probe-timeout-sec`: probe main watchdog 超时（默认 `0`，关闭）
- `--probe-log-grace-sec`: probe 日志缓冲等待（默认 `10`）
- `--force-probe`: 即使已有 probe manifest 也强制重跑

### 其他

- `--out-dir`: 输出目录（默认 `output/qc192_batch_out`）
- `--dry-run`: 只打印提交，不实际跑

---

## 7. 输出目录结构

每个尺寸目录示例（`5x72`）：

```text
<out-dir>/5x72/
  configs/
    qc192_M5.cnfg
  toml/
    probe_M5.toml
  probe_autofer/
    M5_probe_id1/adaptive/...
  logs/
    qc192_M5_id1_snr3p5.log
    ...
  snr_points.json
  summary/
    fer_ranking.csv
    rber_fer_points.csv
  plots/
    rber_fer_M5.png
```

---

## 8. 排名与图的定义

- 排名：按“最终最深点位”的 `LDPC FER` 升序排序，输出 `fer_ranking.csv`
- 图：每个尺寸一张图，横轴 `RAW BER`，纵轴 `LDPC FER`（双对数）

---

## 9. 烟测建议（QC512）

如果当前环境没有 qc192 仿真资源，可先做 qc512 烟测：

```bash
python3 scripts/qc192_batch_sim/batch_sim.py \
  --template-cnfg IBEX/config/Ibex_hd_row5.cnfg \
  --matrix-root IBEX/ibex_matrix \
  --exe IBEX/ssd_fc_src2 \
  --out-dir /tmp/qc_batch_smoke \
  --m-range 5-5 \
  --n-matrices 1 \
  --qc 512 \
  --queue regr_q
```

---

## 10. 常见问题

- probe 报 `never reached FER<=0.9 within range`：
  - 扩大 `--snr-min/--snr-max` 区间
- 批量任务长时间不结束：
  - 这是仿真参数导致（如 `max_err_num` 较大 + 高 SNR）
  - 可先用更小仿真预算做烟测，再跑正式参数
- 未生成图片：
  - 环境缺少 `matplotlib`，脚本会给出 warning 并跳过绘图

---

## 11. 简化版脚本快速开始

```bash
python3 scripts/qc192_batch_sim/batch_sim_simple.py \
  --template-cnfg <template.cnfg> \
  --matrix-root <qc192_matrix_root> \
  --exe <sim_exe> \
  --out-dir <output_dir_simple> \
  --m-range 19-45 \
  --n-matrices 5 \
  --probe-id 1 \
  --qc 192 \
  --queue regr_q \
  --parallel-m 4 \
  --batch-max-in-flight 20 \
  --bpeek-interval-sec 3600
```

简化版输出重点：

- 每个尺寸会生成 `snr_points_simple.json`（记录 `main_start`、仿真点位、排名点位）。
- 实时排名文件：
  - `summary/realtime/ranking_snr*.csv`
- 最终排名文件：
  - `summary/ranking_snr*.csv`（三份，对应 `+0.35/+0.4/+0.45`）
- 点位导出与作图（非 dry-run）：
  - `summary/rber_fer_points.csv`
  - `plots/rber_fer_M<M>.png`

---

## 12. 一键整理已有 summary CSV

`organize_summary_csv.py` 用于整理已经跑完的 batch 结果。它不会重新提交仿真，而是把各个 size 目录下已有的 `summary/*.csv` 汇总、重排，并可选地把 best matrix 文件一并挑出来。

适用输入目录：

- `batch_sim.py` 或 `batch_sim_simple.py` 的输出根目录，例如：
  - `output/qc192_batch_out`
  - `output/qc192_batch_simple_out`
- 目录下需要存在按 size 划分的子目录，例如：
  - `<input-root>/19x197/summary/*.csv`
  - `<input-root>/19x197/summary/realtime/*.csv`
  - `<input-root>/20x198/summary/*.csv`

最常用命令：

```bash
python3 scripts/qc192_batch_sim/organize_summary_csv.py \
  --input-root <batch_output_root> \
  --out-dir <organized_out_dir>
```

如果还希望额外导出每个位点对应的 best matrix 文件，需要再提供 `--matrix-root`：

```bash
python3 scripts/qc192_batch_sim/organize_summary_csv.py \
  --input-root <batch_output_root> \
  --out-dir <organized_out_dir> \
  --matrix-root <matrix_root>
```

完整示例：

```bash
python3 scripts/qc192_batch_sim/organize_summary_csv.py \
  --input-root output/qc192_batch_simple_out \
  --out-dir output/qc192_batch_simple_out/summary_organized \
  --matrix-root IBEX/ibex_matrix \
  --collect-deepest-per-size
```

参数说明：

- `--input-root`
  - 必填含义上的主参数
  - 指向 batch 输出根目录，脚本会递归搜索其中所有 `<size>/summary/*.csv`
  - 默认值：`output/qc192_batch_simple_out`
- `--out-dir`
  - 整理后的输出目录
  - 缺省时默认写到 `<input-root>/summary_organized`
- `--matrix-root`
  - 可选
  - 只有在需要导出 `firstPointMatrix/secondPointMatrix/thirdPointMatrix` 或 deepest-point matrix 时才需要提供
  - 目录结构应为 `<matrix-root>/<size>/{matrix,fade_matrix,occupied_matrix}/...`
- `--best-matrix-dir`
  - 可选
  - 自定义 best matrix 输出目录
  - 缺省时默认写到 `<out-dir>/best_point_matrices`
- `--collect-deepest-per-size`
  - 可选
  - 除了三个位点外，再为每个 size 额外导出“当前最深有效点”的 best matrix
- `--no-realtime`
  - 不纳入 `summary/realtime/*.csv`
- `--no-copy`
  - 不复制原始 CSV 到 `by_size/`，只生成 `merged/*.csv`

运行后输出内容：

- `by_size/<size>/...`：按尺寸复制后的原始 CSV（保留相对目录结构）
- `merged/summary_file_index.csv`：所有 summary 文件清单
- `merged/all_summary_rows.csv`：所有 summary 行汇总（统一加了 size/m/n/source 等元数据）
- `merged/ranking_final_all.csv`：汇总 `summary/ranking_snr*.csv`
- `merged/ranking_realtime_all.csv`：汇总 `summary/realtime/ranking_snr*.csv`
- `merged/rber_fer_points_all.csv`：汇总 `summary/rber_fer_points.csv`
- `merged/fer_ranking_all.csv`：汇总 `summary/fer_ranking.csv`
- `best_point_matrices/firstPointMatrix|secondPointMatrix|thirdPointMatrix/...`：
  - 按点位聚合所有尺寸（不再按 size 分目录）：
    - `matrix/`
    - `fade_matrix/`
    - `occupied_matrix/`
  - 不同尺寸拷贝到同一目录时，文件名会自动加 `size/id` 前缀，避免覆盖
- `merged/best_matrix_by_point.csv`：记录三个位点最终选中的 size / matrix_id / snr / fer 与拷贝结果
  - `ranking_kind` 字段表示来源：`ranking_final`（优先）或 `ranking_realtime`（final 缺失时回退）
  - `ldpc_fer <= 0` 的行会被过滤，不参与 best 评选
- 可选（`--collect-deepest-per-size`）：
  - `best_point_matrices/deepestPointMatrix/{matrix,fade_matrix,occupied_matrix}`：
    - 每个 size 取“当前已有结果的最深 SNR 点”上的 best matrix，并聚合到一起
  - `merged/best_matrix_deepest_by_size.csv`：记录每个 size 的 deepest 点选取结果

可选参数：

- `--no-realtime`：不纳入 `summary/realtime/*.csv`
- `--no-copy`：不复制原始 CSV，只生成 `merged/*.csv`
- `--best-matrix-dir`：自定义 `firstPointMatrix/secondPointMatrix/thirdPointMatrix` 输出目录
- `--collect-deepest-per-size`：额外导出“每个 size 在当前最深点的 best matrix”

注意：

- 如果不传 `--matrix-root`，脚本仍然会完成 CSV 汇总，但不会导出任何 matrix 文件。
- `best matrix` 的选择基于 `ranking_snr*.csv`；若同一 size 同一 snr 同时存在 final 和 realtime，优先使用 final。
- `deepestPointMatrix` 不是固定第三个点，而是“当前已有有效结果中 SNR 最大的点”。
