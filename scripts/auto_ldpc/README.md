# auto\_ldpc 自动化工具

## 功能概述
- 基于 `submit_ldpc_bsub.py` 的扩展，自动执行矩阵生成、SNR 仿真、ERR\_INJ 检查、深挖高 SNR 以及结果归档。
- 通过 SQLite 记录矩阵/任务状态，实现断点恢复与多轮迭代。
- 解析日志生成 XLSX 报表，并绘制最佳矩阵与基准的 SNR-FER 对比曲线。
- 支持配置化参数、批量生成、停止条件和可选的基准数据。

## 目录结构
- `automation.py`：主控脚本，驱动全流程。
- `config.py`：加载 YAML/JSON 配置并应用默认值。
- `state.py`：SQLite 状态管理。
- `log_parser.py`：读取仿真日志提取统计值。
- `jobs.py`：封装 LSF 提交命令。
- `commands.py`：统一命令执行入口。
- `plotting.py`：绘制最佳矩阵的 SNR-FER 曲线。
- `xlsx_export.py`：导出仿真记录为 XLSX。
- `sample_config.yaml`：示例配置，可复制后按需修改。

## 环境要求
- Python 3.8+（建议在 WSL 中运行）。
- 依赖包：
  - `pip install pyyaml openpyxl matplotlib`
- 可访问 LSF 集群命令（`bsub`、`bjobs`、`bwait`）。
- 确保已有 `submit_ldpc_bsub.py`、`GenLDPC`、`LDPC_Sim_Gen4_dq` 等项目可用。

## 配置说明
1. 复制示例：`cp scripts/auto_ldpc/sample_config.yaml scripts/auto_ldpc/my_config.yaml`
2. 关键字段：
   - `paths`：指定执行文件、矩阵目录、结果目录等。支持 `output_root` + `${OUTPUT_ROOT}` 占位符，便于集中切换输出前缀。
   - `generation`：矩阵生成参数与批次大小（`batch_size` 为持续提交数量）。
   - `simulation`：SNR 列表、锚定点、深挖点、队列等。
   - `selection`：`top_n`（锚定 SNR 最优矩阵数量）、`err_inj_threshold`。
   - `err_inj`： ERR\_INJ 模型配置与命令模板。
   - `limits`：最大矩阵数、运行时间、停止标志文件等。
   - `log_parser`：不同模型的 STATISTICS 关键字与 FER 字段名。
   - `xlsx`：导出的报表文件名。
   - `plot`：最佳矩阵曲线输出路径；若需要基准线，设置 `baseline.path` 指向包含 `SNR`/`FER` 列的 XLSX。

## 运行方式
- 干跑检查：  
  `python3 scripts/auto_ldpc/automation.py --config scripts/auto_ldpc/my_config.yaml --dry-run --once`
- 实际运行（无限循环，直到命中停止条件）：  
  `python3 scripts/auto_ldpc/automation.py --config scripts/auto_ldpc/my_config.yaml`
- 单轮执行（用于手动调度或测试）：  
  `python3 scripts/auto_ldpc/automation.py --config scripts/auto_ldpc/my_config.yaml --once`

## 输出内容
- 状态数据库：`database.path`（可用 `${OUTPUT_ROOT}/state/...` 设置）
- 仿真日志：`paths.sim_out_root` 下自动生成 `awgn/`、`err_inj/`
- 导出的 XLSX：`paths.perf_dir` + `xlsx.filename`
- 绘图文件：`plot.output`（默认 `${OUTPUT_ROOT}/perf/best_plot.png`）
- 连续运行时，可通过创建 `limits.stop_flag_file` 指向的文件来触发停机。

## 调试建议
- `--dry-run` 时会打印所有将要执行的命令，便于核对。
- 若需清空状态，可删除 `state.db`，同时确保历史日志与矩阵目录一致。
- 如果日志键名有调整，务必同步更新 `log_parser` 的 `stat_keys` 与 `fer_key`。
