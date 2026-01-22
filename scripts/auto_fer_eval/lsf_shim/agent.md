# LSF Shim（`bsub`/`bjobs`/`bkill` 仿真器）设计说明（agent.md）

## 0. 背景与目标

我们希望在**家用电脑**上复现 `auto_fer_eval` 在 LSF 环境下的最小调度闭环，用于验证：

- 提交作业（`bsub`）→ 轮询状态（`bjobs`）→ 取消作业（`bkill`）
- 生成日志（`-o`）→ `auto_fer_eval` 解析日志（至少包含 `[STATISTICS] LDPC FER`）→ 进入后续调度阶段

该 shim 的核心是：**输出格式严格匹配 `auto_fer_eval/executors.py` 的解析假设**，而不是完整复刻 IBM LSF 的全部行为。

## 1. 与 `auto_fer_eval` 的接口契约（必须满足）

### 1.1 `bsub`（提交）

`auto_fer_eval` 典型调用形态（顺序可能略有差异）：

`bsub -o <log_path> -cwd <cwd> [-q <queue>] [-J <job_name>] [bsub_extra...] <cmd...>`

shim 必须做到：

- 将 `<cmd...>` 在 `<cwd>` 下运行（等价于 LSF 的 `-cwd`）。
- 将 stdout/stderr 合并写入 `<log_path>`（等价于 `-o`）。
- `bsub` 的 stdout 必须包含可解析 job id 的片段：`Job <123>`（脚本用正则提取）。
- **允许存在未知的 `bsub` 额外参数**：shim 不应因未知参数直接退出失败。

### 1.2 `bjobs`（查询）

`auto_fer_eval` 典型调用形态：

`bjobs [bjobs_extra...] <job_id>`

shim 必须做到：

- 输出至少两行：header 行 + job 行。
- job 行需满足：第 1 列为 `<job_id>`；**第 3 列为 `STAT`**。
- `STAT` 至少支持：`PEND` / `RUN` / `DONE` / `EXIT`。

### 1.3 `bkill`（取消）

`auto_fer_eval` 典型调用形态：

`bkill [bkill_extra...] <job_id>`

shim 必须做到：

- 尽力终止该 job 对应的进程（推荐终止整个进程组）。
- 状态标记为 `EXIT`（你已确认不需要区分 KILLED/EXIT）。

## 2. 状态模型（最小集合）

| 状态 | 含义 | 触发条件 |
|---|---|---|
| `PEND` | 排队等待 | 默认模拟 $10s$（见 3.2） |
| `RUN` | 正在运行 | 进入执行阶段 |
| `DONE` | 成功结束 | 退出码 `rc==0` |
| `EXIT` | 失败/被 kill | 退出码 `rc!=0` 或 `bkill` |

说明：该状态集合已经覆盖 `auto_fer_eval` 的 `LsfExecutor.poll()` 对 `STAT` 的最小需求。

## 3. Registry（作业注册表）设计

### 3.1 存储位置

使用 JSON 文件持久化 registry，路径规则：

- 优先读取环境变量：`LSF_SHIM_DB=/path/to/jobs.json`
- 未设置则默认：`~/.lsf_shim/jobs.json`

### 3.2 PEND 默认 10 秒

你已确认需要模拟排队：

- 默认 `PEND` 持续 $10s$，随后进入 `RUN`。
- 允许未来扩展：通过环境变量 `LSF_SHIM_PEND_SEC` 覆盖默认值（可选实现点）。

### 3.3 作业记录字段（建议）

每个 job 一条记录（键：`job_id`），至少包含：

- `job_id`：递增整数（全局唯一）
- `state`：`PEND/RUN/DONE/EXIT`
- `rc`：退出码（DONE/EXIT 时可用）
- `submit_time` / `start_time` / `end_time`：时间戳
- `cwd`：工作目录
- `cmd`：命令数组（便于审计）
- `log_path`：日志路径（`-o`）
- `queue` / `job_name`：仅用于记录（不实现真实调度策略）
- `pid` / `pgid`：用于 `bkill` 定位并终止进程（建议使用进程组）

### 3.4 并发一致性

registry 需要支持并发读写（例如多个 `bsub` 同时提交）：

- 写入必须加锁（文件锁），避免 JSON 被并发写坏。
- `bjobs` 可在读 registry 后做一次“pid 是否仍存在”的校验，以修正状态。

## 4. Runner 策略（为何需要额外进程）

为了满足：

- `bsub` **立即返回 job_id**
- 同时又要在后台持续跟踪 `<cmd...>` 的退出并更新 registry

建议实现一个独立的 runner 进程：

1. `bsub`：解析参数 → 写入 registry（state=`PEND`）→ 启动 runner（后台）→ 输出 `Job <id>` → 退出
2. runner：等待 `PEND` 时间 → 启动 `<cmd...>`（stdout/stderr 重定向至 log）→ 等待退出 → 写回 `DONE/EXIT` + `rc`

## 5. `bkill` 语义（实现建议）

- 首选对进程组发送 `SIGTERM`，等待短时间（例如 2s），未退出再 `SIGKILL`。
- 更新 registry：`state=EXIT`，写入 `end_time`，可记录 `cancel_reason="bkill"`（可选）。

## 6. 非目标（明确不做）

shim 不追求：

- 资源表达式/队列策略/公平性（`-R/-W/-n` 等）真实调度
- LSF 全量输出字段与行为一致（仅保证 `auto_fer_eval` 需要的最小子集）

## 7. 集成方式（两种选一）

### 7.1 PATH 覆盖（推荐）

把 shim 的可执行文件目录加入 `PATH`，让 `auto_fer_eval` 默认调用到 shim：

- `export PATH="/abs/path/to/auto_fer_eval/lsf_shim/bin:$PATH"`

### 7.2 TOML 显式指定（可选）

在配置中指定 `bsub_cmd/bjobs_cmd/bkill_cmd`（如果脚本支持显式路径；否则使用 PATH 覆盖）。

## 8. 验收标准（最小可用）

在 `executor=lsf` 配置下，跑通至少一个 dummy case：

- `bsub` 能返回形如 `Job <123>` 的输出
- `bjobs 123` 的第三列状态能从 `PEND/RUN` 变为 `DONE`
- log 文件按 `-o` 路径生成，并包含可解析的 `[STATISTICS] LDPC FER`（或至少包含 `[SIM] LDPC FER` 以便 fallback）
- `bkill` 能终止任务，`bjobs` 显示 `EXIT`

## 9. 已对齐的确认项（来自你的答复）

1. 未知参数：允许（shim 尽量忽略，不因未知参数崩溃）
2. 排队：需要，默认 `PEND=10s`
3. kill 语义：标记为 `EXIT` 可以
4. registry 默认路径：`~/.lsf_shim/jobs.json` 可接受（也支持 `LSF_SHIM_DB` 覆盖）

