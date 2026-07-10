# auto_throughput_eval：固定网格的 throughput（aver\_iter）评估与调度（SNR/K）

本目录提供一个 **不依赖第三方库** 的 Python 脚本 `auto_throughput_eval/auto_throughput_eval.py`，用于在给定的
网格区间 $[x\\_{low}, x\\_{high}]$ 上按步长 $x\\_{step}$ 逐点运行仿真（`AWGN/TAWGN` 时 $x$ 为 SNR，`ERR_INJ` 时 $x$ 为 $K$），并从日志 `[STATISTICS]` 中提取：

- `RAW BER`（在此脚本中记为 `RBER`）
- `Decoder average iteration(s)`（在此脚本中记为 `aver_iter`，兼容 `Retry/Fast/IBEX` 等变体）

与 `auto_fer_eval` 的区别：不进行区间定位或自适应下探，仅做固定网格扫描。

另提供一版实时进度脚本 `auto_throughput_eval/auto_throughput_eval_realtime.py`：
- 保留原有固定网格与补跑逻辑
- 运行期间实时刷新 `throughput_progress.csv`
- 字段包含：`SNR/K, RBER, LDPC_FER, Fail_CW, Packets, AvgIter, JobState, JobID`
- `AvgIter` 仅在任务最终完成且日志可完整解析后更新
- 当运行中统计达到 `Packets>=1000` 且 `LDPC_FER>=0.99` 时，会自动 kill 该点，并按 cnfg 中的最大 iteration 写入 `AvgIter`
- `lsf` 模式下会优先读取当前日志，并在需要时用 `bpeek` 补抓运行中 stdout

已实现：
- `local`：本地并发运行 + 日志解析 + `manifest.json` 记录
- `lsf`：`bsub/bjobs/bkill` 提交与轮询
- `lsf` 断点重续：`jobs.json` 持久化提交记录；脚本重启后会优先接管仍处于 `PEND/RUN` 的旧 job，避免重复提交同一点
- 完整性补跑：若任务结束但日志不完整（缺 `[STATISTICS]` 或关键字段缺失），会自动重提该点，直到该点日志完整
- 结果导出：从 `manifest.json`/日志汇总并写入 `csv`

---

## 1. 运行方式（local）

命令：

`python3 auto_throughput_eval/auto_throughput_eval.py run --config <your_case>.toml`

若需要实时进度表：

`python3 auto_throughput_eval/auto_throughput_eval_realtime.py run --config <your_case>.toml`

只打印计划、不执行：

`python3 auto_throughput_eval/auto_throughput_eval.py run --config <your_case>.toml --dry-run`

## 1.1 运行方式（lsf）

在 TOML 中设置 `runner.executor = "lsf"` 后，使用同样命令运行即可。

建议先用 `--dry-run` 核对将要提交的 `bsub` 命令（含 `-J` job name 与 `-o` 日志路径）。

---

## 1.2 仅导出 csv（不提交/不运行）

当你已经手动提交并确认所有日志已生成，可直接导出：

`python3 auto_throughput_eval/auto_throughput_eval.py export --config <your_case>.toml`

---

## 2. 配置文件（TOML）

要点：
- 轴选择：
  - `AWGN/TAWGN`：默认 `axis_type="snr"`，使用 `snr_low/snr_high/snr_step`
  - `ERR_INJ`：默认 `axis_type="k"`，使用 `k_low/k_high/k_step`（也兼容复用 `snr_*` 作为别名）
- `max_sim_num`：每个点固定跑的 CW 数（会写入 `.cnfg` 的 `maximum simulation number`）。
- `max_err_num`：建议固定为 0（会写入 `.cnfg` 的 `maximum error number`）。
- `config_tag`、`out_dir`：由你在 TOML 中定义（后续用于 LSF job name 与结果组织）。
- `runner.fail_fast`：默认 `false`。若设为 `true`，则任一任务 `EXIT` 或 `DONE` 但日志解析失败时，会立即取消其它 in-flight 任务并退出；吞吐测试通常仿真代价很小，建议保持 `false` 以避免因共享盘刷盘延迟导致的误判。
- `lsf.use_cwd`：是否在 `bsub` 命令里附带 `-cwd <workdir>`（默认 `true`）。若你的集群策略不允许或不需要 `-cwd`，可设为 `false`。

---

## 3. 输出

对每个 case：
- 日志：`<out_dir>/<case>/ <log_prefix>_snrX.log` 或 `<log_prefix>_kX.log`
- 记录：`<out_dir>/<case>/manifest.json`
- 作业数据库：`<out_dir>/<case>/jobs.json`
- 表格：`<out_dir>/<case>/throughput.csv`
- 实时进度表（realtime 版本）：`<out_dir>/<case>/throughput_progress.csv`

说明：
- `manifest.json` 只接受带最终 `[STATISTICS]` 的完整日志；被 kill 或未完整落盘的 partial log 不会再被当成完成点复用。
- `run` 模式下会持续重提不完整点，直到所有网格点都有完整日志后才导出 `throughput.csv` 并结束（`fail_fast=true` 时仍按快速失败策略中止）。
- `jobs.json` 仅对 `lsf` 模式生效；脚本重启后会先查询旧记录并尝试接管活 job，再决定是否重提。
