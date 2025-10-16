# submit_ldpc_bsub.sh 使用指南

## 脚本概述

- `scripts/submit_ldpc_bsub.sh` 提供“一键式”流水：可选矩阵生成（LSF 或本地）、标准重命名、批量提交 `ssd_fc` / `ssd_fc_dq` 仿真任务。
- 默认执行体：`gen4_ldpc_sim/ssd_fc_dq`；默认配置文件：`gen4_ldpc_sim/config/sdec_dq.cnfg`。脚本会把所有路径转换为工作区根目录 `/workspaces/LDPC-codes` 下的绝对路径。
- 所有传入路径（生成输出、矩阵目录、仿真输出等）都会经过 `to_abs_path` 处理，可使用仓库内相对路径或绝对路径。

## 快速上手

```bash
# 预览命令，不执行
bash scripts/submit_ldpc_bsub.sh \
  --gen-n 149 --gen-k 129 --gen-count 1 \
  --snr "4.9" \
  --dry-run

# 确认无误后去掉 --dry-run 真正提交
```

**执行顺序（非 dry-run）：**

1. 若指定 `--gen-n / --gen-k / --gen-count`，先提交 `ldpc_gen` 数组任务 `[1-GEN_COUNT]` 并等待全部完成。
2. 调用 `scripts/rename_ldpc_matrices.sh M N --apply` 重命名矩阵文件（与 `test_ldpc_sim.sh` 保持一致），记录新增矩阵编号。
3. 对新增编号（若无新增则回退扫描目录全部矩阵）与 SNR 列表生成 `bsub` 命令，并通过第 6 参数 + `LDPC_MATRIX_DIR` 传递矩阵目录给仿真器。

## 核心选项

| 选项 | 说明 |
| --- | --- |
| `--exec PATH` | 仿真器可执行文件，默认 `gen4_ldpc_sim/ssd_fc_dq`。切换 `ssd_fc` 时写 `--exec gen4_ldpc_sim/ssd_fc`。 |
| `--config PATH` | 仿真配置；`ssd_fc_dq` 默认 `gen4_ldpc_sim/config/sdec_dq.cnfg`，`ssd_fc` 常用 `config/sdec.cnfg`。 |
| `--gen-bin PATH` | 矩阵生成器，默认 `GenLDPC/ldpc_gen`。 |
| `--gen-n / --gen-k / --gen-count` | 传递给 `ldpc_gen` 的 N、K、phase（job array 大小）。 |
| `--gen-mode MODE` | `lsf`（默认）/`local`，决定生成阶段的执行方式。 |
| `--gen-out-base DIR` | 生成输出根目录，默认 `GenLDPC/output`，同时导出 `GENLDPC_OUT_DIR`。 |
| `--matrix-dir DIR` | 仿真用矩阵目录；未指定时若执行过生成，默认 `GenLDPC/output/<m>x<n>/matrix`。 |
| `--sim-matrix-dir DIR` | 额外指定传给仿真器的目录，并导出 `LDPC_MATRIX_DIR`。 |
| `--out-dir DIR` | 仿真输出目录（默认 `runs/ldpc_batch`），每个矩阵单独落在 `matrix<ID>/snr<SNR>.{out,err}`。 |
| `--rename-enable / --no-rename` | 是否执行标准重命名（默认开启，与 `test_ldpc_sim.sh` 一致）。 |
| `--build` | 提交前自动执行 `make -C GenLDPC` 与 `make -C gen4_ldpc_sim`。 |
| `--snr "…" / --snr-seq A:B:C` | 指定 SNR 列表或区间（二选一），区间由 `awk` 生成浮点序列。 |
| `--dry-run` | 仅打印命令，不执行。 |

## 两类仿真程序接口

| 项目 | `ssd_fc` | `ssd_fc_dq` |
| --- | --- | --- |
| 默认 exec | `gen4_ldpc_sim/ssd_fc` | `gen4_ldpc_sim/ssd_fc_dq` |
| 默认 config | `gen4_ldpc_sim/config/sdec.cnfg` | `gen4_ldpc_sim/config/sdec_dq.cnfg` |
| 调用参数 | `./ssd_fc LDPC <config> <channel> <snr> <matrix_id> [matrix_dir]` | 同左；额外读取 `LDPC_MATRIX_DIR` 环境变量 |
| 典型模式 | 完整 FC 仿真（含 MCRC/LBA） | DQ 版本仿真，需要 `make ssd_fc_dq` |
| 输出 | 常规统计 | 增加 DQ 指标与迭代统计 |

脚本生成的 `inner_cmd` 会：

1. `cd` 到可执行文件所在目录；
2. 如设置 `--sim-matrix-dir` 或 `--matrix-dir`，导出 `LDPC_MATRIX_DIR` 并传递为第 6 参数；
3. 执行 `./<exec> LDPC <config> AWGN <SNR> <matrix_id> [matrix_dir]`。

因此只需替换 `--exec / --config` 即可在 `ssd_fc` 与 `ssd_fc_dq` 之间切换。

## 路径对齐示例

| 变量 | 默认值（相对仓库） | 转换后示例 |
| --- | --- | --- |
| `EXEC` | `gen4_ldpc_sim/ssd_fc_dq` | `/workspaces/LDPC-codes/gen4_ldpc_sim/ssd_fc_dq` |
| `CONFIG` | `gen4_ldpc_sim/config/sdec_dq.cnfg` | `/workspaces/LDPC-codes/gen4_ldpc_sim/config/sdec_dq.cnfg` |
| `GEN_BIN` | `GenLDPC/ldpc_gen` | `/workspaces/LDPC-codes/GenLDPC/ldpc_gen` |
| `GEN_OUT_BASE` | `GenLDPC/output` | `/workspaces/LDPC-codes/GenLDPC/output` |
| `MATRIX_DIR`（默认） | `GenLDLC/output/<m>x<n>/matrix` | `/workspaces/LDPC-codes/GenLDPC/output/20x149/matrix` |
| `SIM_MATRIX_DIR` | 自定义 | `/workspaces/LDPC-codes/...` |
| `OUT_DIR` | `runs/ldpc_batch` | `/workspaces/LDPC-codes/runs/ldpc_batch` |
| `rename_ldpc_matrices.sh` | `scripts/rename_ldpc_matrices.sh` | `/workspaces/LDPC-codes/scripts/rename_ldpc_matrices.sh` |

## 日常使用建议

- **先预览再提交**：添加 `--dry-run` 检查 job array、重命名及 `bsub` 细节，确认无误后再真正执行。  
- **多节点并行**：建议为每个任务指定唯一 `--gen-out-base`（如 `GenLDPC/output/job_${LSB_JOBID}`），同时使用 `--no-rename`，在汇总节点统一重命名与提交，避免冲突。  
- **复用既有矩阵**：无需再生成，直接提供 `--matrix-dir` 或 `--sim-matrix-dir` 即可。  
- **多 SNR 批量**：`--snr` 列出或 `--snr-seq` 生成浮点序列，避免手动输入大量值。  
- **切换仿真器**：`ssd_fc`/`ssd_fc_dq` 仅需调整 `--exec` 与 `--config`；其余流程保持一致。  

只要按照上述规则设置路径和参数，即可在 LSF 环境中快速完成矩阵生成、重命名与仿真提交流水。如需扩展到其它 channel（如 BSC/ERR_INJ），只需调整脚本内 `inner_cmd` 的 `<channel>` 和 `<snr>` 即可。
