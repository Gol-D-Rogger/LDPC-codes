# auto\_ldpc\_local

`auto_ldpc_local` 复用了 `scripts/auto_ldpc` 的配置、状态与日志解析模块，但把矩阵生成、AWGN/ERR_INJ 仿真、深挖等步骤改为 **本地直接执行**，方便在没有 LSF/bsub 的环境快速调试。

## 快速开始
1. 复制示例配置：
   ```bash
   cp scripts/auto_ldpc_local/sample_config.yaml scripts/auto_ldpc_local/my_config.yaml
   ```
2. 根据实际环境修改 `my_config.yaml`：
   - `paths.output_root` 设为 WSL 可写路径；
   - 如在 Windows + WSL 场景运行，可把 `system.command_prefix` 改成 `["wsl"]`；
   - 其他字段与 `auto_ldpc` 相同（SNR 列表、top_n、err_inj 阈值等）。
3. 确保 `scripts/rename_ldpc_matrices.sh` 为 LF 换行：
   ```bash
   dos2unix scripts/rename_ldpc_matrices.sh
   ```
4. 运行：
   ```bash
   python scripts/auto_ldpc_local/automation.py --config scripts/auto_ldpc_local/my_config.yaml
   ```
   如需仅执行一轮，可追加 `--once`；若想查看命令但不实际执行，加 `--dry-run`。

## 注意事项
- 流程完全沿用 `auto_ldpc`：矩阵生成 → AWGN → ERR_INJ → 深挖 → 日志解析 → 绘图 → XLSX 导出；
- 生成的矩阵默认仍会通过 `rename_ldpc_matrices.sh` 重命名，确保 `_QC_H_0.txt` 等文件存在；
- 所有结果（矩阵、日志、state.db、plot 等）都会写到配置中的 `output_root` 下，便于清理或迁移；
- 如需要与远程 LSF 流程对齐，只需将配置与 `scripts/auto_ldpc/sample_config.yaml` 一致即可。*** End Patch
