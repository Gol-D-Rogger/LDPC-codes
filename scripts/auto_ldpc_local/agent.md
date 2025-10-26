# auto_ldpc_local - 本地自动化LDPC矩阵探索系统

## 📋 概述

`auto_ldpc_local` 是 `auto_ldpc` 的本地测试版本，完整复用了 auto_ldpc 的逻辑和链路，但**不使用 bsub/LSF 提交任务**，而是在本地直接运行仿真程序。适用于：

- ✅ 本地开发和调试
- ✅ 功能验证和测试
- ✅ 小规模快速实验
- ✅ Windows/WSL 环境测试

## 🏗️ 架构设计

### 核心差异

| 特性 | auto_ldpc (生产环境) | auto_ldpc_local (本地测试) |
|------|---------------------|--------------------------|
| 任务提交 | bsub (LSF) | subprocess (本地进程) |
| 并行执行 | 集群调度 | 本地串行/简单并行 |
| 日志管理 | LSF日志系统 | 直接文件输出 |
| 资源管理 | 集群队列 | 本地CPU/内存 |
| 代码复用 | - | 继承 auto_ldpc 所有模块 |

### 模块继承关系

```
auto_ldpc/
  ├── config.py         ← 完全复用
  ├── state.py          ← 完全复用
  ├── log_parser.py     ← 完全复用
  ├── xlsx_export.py    ← 完全复用
  ├── plotting.py       ← 完全复用
  ├── jobs.py           ← 完全复用（数据结构）
  └── commands.py       ← 替换为 local_runner.py

auto_ldpc_local/
  ├── local_runner.py   ← 新增：本地任务执行器
  ├── automation.py     ← 修改：使用 local_runner
  ├── sample_config.yaml
  └── agent.md
```

## 📦 文件结构

```
scripts/auto_ldpc_local/
├── agent.md              # 本文档
├── automation.py         # 主自动化脚本（修改版）
├── local_runner.py       # 本地任务执行器（新增）
├── sample_config.yaml    # 示例配置文件
└── README.md             # 快速开始指南
```

## 🔧 实现细节

### 1. local_runner.py - 本地任务执行器

**功能：** 替代 bsub 提交，在本地直接执行仿真程序

```python
class LocalRunner:
    """本地任务执行器，替代 LSF bsub"""
    
    def run_job(self, job: SimulationJob) -> bool:
        """
        直接在本地运行仿真任务
        
        参数:
            job: SimulationJob - 任务信息
            
        返回:
            bool - 是否成功启动
        """
        # 1. 解析命令
        # 2. 设置环境变量
        # 3. 启动子进程
        # 4. 重定向日志
        # 5. 后台运行或等待完成
```

**关键特性：**
- ✅ 支持后台运行（异步）
- ✅ 支持同步等待（串行测试）
- ✅ 自动创建日志文件
- ✅ 进程管理和清理
- ✅ 错误处理和重试

### 2. automation.py - 修改版主脚本

**修改点：**

```python
# 原版 (auto_ldpc)
from scripts.auto_ldpc.commands import CommandRunner
runner = CommandRunner(bsub_prefix, dry_run=dry_run)

# 本地版 (auto_ldpc_local)
from scripts.auto_ldpc_local.local_runner import LocalRunner
runner = LocalRunner(max_parallel=2, dry_run=dry_run)
```

**保持不变的部分：**
- ✅ 所有数据结构（SimulationJob, ConfigBundle, StateStore）
- ✅ 日志解析逻辑
- ✅ Top-N 选择算法
- ✅ 深挖触发逻辑
- ✅ XLSX 导出和绘图
- ✅ 状态管理和数据库

### 3. sample_config.yaml - 本地配置模板

**关键配置差异：**

```yaml
# 本地版本的特殊配置
system:
  command_prefix: null           # 不需要 bsub 前缀
  max_parallel: 2                # 限制并行任务数（本地资源有限）
  
generation:
  mode: "local"                  # 本地生成模式
  disable_rename: false          # 启用重命名（重要！）
  
simulation:
  queue: null                    # 不使用队列
  job_prefix: "ldpc_local"       # 标记为本地任务
  
limits:
  max_matrices: 4                # 本地测试建议小数量
  sleep_seconds: 5               # 更快的检查间隔
```

## 🚀 使用方法

### 快速开始

```bash
# 1. 进入脚本目录
cd scripts/auto_ldpc_local

# 2. 编辑配置文件
vim sample_config.yaml

# 3. 运行自动化脚本
python automation.py --config sample_config.yaml

# 4. 单次测试（推荐首次运行）
python automation.py --config sample_config.yaml --once

# 5. 查看结果
ls -l ../../test_local/runs/
```

### 本地测试工作流

```
1. 启动自动化脚本
   ↓
2. 生成 batch_size 个矩阵
   ↓
3. 本地运行 AWGN 仿真（串行或简单并行）
   ↓
4. 解析日志，提取 FER
   ↓
5. 选择 Top-N 矩阵
   ↓
6. 本地运行 ERR_INJ 测试
   ↓
7. 通过阈值的矩阵进行深挖
   ↓
8. 导出 XLSX 和绘图
   ↓
9. 等待 sleep_seconds 秒后重复（或达到 max_matrices 停止）
```

## 🔍 关键实现代码

### local_runner.py 核心逻辑

```python
import subprocess
import threading
from pathlib import Path
from typing import Optional, List
import time

class LocalRunner:
    def __init__(self, max_parallel: int = 1, dry_run: bool = False):
        self.max_parallel = max_parallel
        self.dry_run = dry_run
        self.running_processes: List[subprocess.Popen] = []
        self.lock = threading.Lock()
    
    def run_job_async(self, job: SimulationJob) -> bool:
        """异步运行任务（后台）"""
        if self.dry_run:
            print(f"[DRY-RUN] Would run: {job.command}")
            return True
        
        # 等待空闲槽位
        self._wait_for_slot()
        
        # 启动子进程
        process = subprocess.Popen(
            job.command,
            shell=True,
            stdout=open(job.log_path, 'w'),
            stderr=open(job.err_path, 'w'),
            cwd=job.work_dir
        )
        
        with self.lock:
            self.running_processes.append(process)
        
        # 后台监控线程
        threading.Thread(
            target=self._monitor_process,
            args=(process, job),
            daemon=True
        ).start()
        
        return True
    
    def run_job_sync(self, job: SimulationJob) -> int:
        """同步运行任务（等待完成）"""
        if self.dry_run:
            print(f"[DRY-RUN] Would run: {job.command}")
            return 0
        
        with open(job.log_path, 'w') as log_f, \
             open(job.err_path, 'w') as err_f:
            result = subprocess.run(
                job.command,
                shell=True,
                stdout=log_f,
                stderr=err_f,
                cwd=job.work_dir
            )
        
        return result.returncode
```

### automation.py 修改点

```python
# 导入本地运行器
from scripts.auto_ldpc_local.local_runner import LocalRunner

class AutomationRunner:
    def __init__(self, bundle: ConfigBundle, dry_run: bool = False):
        # ... 其他初始化代码保持不变 ...
        
        # 使用本地运行器替代 CommandRunner
        system_cfg = self.cfg.get("system", {})
        max_parallel = system_cfg.get("max_parallel", 1)
        self.runner = LocalRunner(max_parallel=max_parallel, dry_run=dry_run)
    
    # submit_awgn_runs, submit_errinj_runs 等方法保持不变
    # 因为它们都通过 submit_jobs(self.runner, jobs) 统一接口
```

## ⚙️ 配置说明

### system 配置块

```yaml
system:
  command_prefix: null           # 本地不需要 bsub 前缀
  max_parallel: 2                # 最大并行任务数
  run_mode: "async"              # "async" 或 "sync"
```

### 本地路径配置

```yaml
paths:
  workspace_root: "."
  output_root: "../../test_local"            # 本地测试输出目录
  exec_path: "../../gen4_ldpc_sim/ssd_fc_dq" # 仿真程序路径（相对或绝对）
  config_path: "../../gen4_ldpc_sim/config/ldpc_test.cnfg"
  gen_bin: "../../GenLDPC/ldpc_gen"
  matrix_dir: "../../GenLDPC/output"
```

### 本地测试建议配置

```yaml
limits:
  max_matrices: 4        # 小数量快速测试
  sleep_seconds: 5       # 快速检查
  
generation:
  batch_size: 2          # 每批2个矩阵
  disable_rename: false  # 必须启用重命名！
  
simulation:
  snr_list: [3.7]        # 只测试 anchor SNR
  deep_snr_list: [3.2, 4.2]  # 深挖2个点
  
selection:
  top_n: 2               # 只选前2个
```

## 🐛 调试和问题排查

### 常见问题

#### 1. 矩阵文件找不到

**症状：**
```
[LDPC] Error opening file .../LDPC_..._1.txt
```

**原因：** 文件名是 `*_1_1.txt` 但期望 `*_1.txt`

**解决：** 确保配置中 `disable_rename: false`

#### 2. 进程卡住

**症状：** 脚本运行但没有输出

**排查：**
```python
# 在 local_runner.py 中添加调试输出
print(f"[DEBUG] 启动进程: {job.command}")
print(f"[DEBUG] 当前运行进程数: {len(self.running_processes)}")
```

#### 3. 日志路径错误

**症状：** 找不到日志文件

**排查：**
```python
# 检查日志路径
print(f"[DEBUG] 日志路径: {job.log_path}")
print(f"[DEBUG] 父目录存在: {job.log_path.parent.exists()}")
```

### 调试模式

```bash
# 使用 dry-run 模式
python automation.py --config sample_config.yaml --dry-run

# 单次迭代测试
python automation.py --config sample_config.yaml --once

# Python unbuffered 模式（实时输出）
python -u automation.py --config sample_config.yaml
```

## 📊 性能对比

### 本地 vs 集群

| 指标 | auto_ldpc (集群) | auto_ldpc_local (本地) |
|------|-----------------|----------------------|
| 并行能力 | 数百任务 | 1-4任务 |
| 启动延迟 | LSF调度延迟 | 立即启动 |
| 适用场景 | 大规模探索 | 快速验证 |
| 资源限制 | 集群配额 | 本地CPU/内存 |
| 调试友好度 | ⭐⭐ | ⭐⭐⭐⭐⭐ |

### 推荐使用场景

**使用 auto_ldpc_local：**
- ✅ 开发新功能时
- ✅ 修复 bug 需要快速验证
- ✅ 测试配置文件
- ✅ 小规模算法验证

**使用 auto_ldpc：**
- ✅ 大规模矩阵探索（数百个）
- ✅ 生产环境长期运行
- ✅ 需要高并行度

## 🔄 迁移到生产环境

当本地测试完成后，迁移到集群很简单：

```bash
# 1. 复制配置文件
cp scripts/auto_ldpc_local/sample_config.yaml scripts/auto_ldpc/my_config.yaml

# 2. 修改配置
vim scripts/auto_ldpc/my_config.yaml
# 修改:
#   - system.command_prefix: "bsub -q regr_q"
#   - limits.max_matrices: 100
#   - generation.batch_size: 10
#   - simulation.snr_list: [2.7, 3.2, 3.7, 4.2]

# 3. 使用生产版本运行
python scripts/auto_ldpc/automation.py --config scripts/auto_ldpc/my_config.yaml
```

配置文件完全兼容，只需调整资源相关参数！

## 📝 开发指南

### 添加新功能

如果要在本地版本添加新功能：

1. **修改 local_runner.py** - 如果涉及任务执行
2. **修改 automation.py** - 如果涉及流程逻辑
3. **复用 auto_ldpc 模块** - 尽量使用现有模块
4. **同步到生产版本** - 验证后合并回 auto_ldpc

### 测试流程

```bash
# 1. 本地小规模测试
python scripts/auto_ldpc_local/automation.py --config test_config.yaml --once

# 2. 本地中等规模测试
# max_matrices: 10, batch_size: 5

# 3. 集群小规模验证
# 使用 auto_ldpc 运行相同配置

# 4. 集群生产运行
# 放大参数到生产规模
```

## 🎯 最佳实践

### 本地测试配置模板

```yaml
# 快速功能验证（5-10分钟）
limits:
  max_matrices: 4
  sleep_seconds: 5

generation:
  batch_size: 2

simulation:
  snr_list: [3.7]
  
# 完整流程测试（30-60分钟）
limits:
  max_matrices: 10
  sleep_seconds: 10

generation:
  batch_size: 5

simulation:
  snr_list: [3.2, 3.7, 4.2]
  deep_snr_list: [2.7, 3.0, 3.2, 3.7, 4.0, 4.2]
```

### 资源管理

```python
# local_runner.py 中的资源控制
class LocalRunner:
    def __init__(self, max_parallel: int = 1, max_memory_mb: int = 4096):
        self.max_parallel = max_parallel
        self.max_memory_mb = max_memory_mb
        # 监控系统资源
        self.check_system_resources()
```

## 🔗 相关链接

- [auto_ldpc 主文档](../auto_ldpc/README.md)
- [配置文件说明](../auto_ldpc/CONFIG.md)
- [LSF vs 本地对比](./COMPARISON.md)
- [故障排查指南](./TROUBLESHOOTING.md)

## 📞 支持

遇到问题？
1. 查看本文档的"调试和问题排查"章节
2. 检查日志文件：`test_local/runs/awgn/*/snr*.log`
3. 使用 `--dry-run` 模式检查命令
4. 使用 `--once` 模式单步调试

---

**版本:** 1.0.0  
**更新日期:** 2025-10-24  
**维护者:** LDPC Team
# auto\_ldpc\_local 说明

1. **准备环境**
   - 进入仓库根目录（WSL 下建议使用 `/mnt/<盘符>/...` 路径）。
   - 将 `scripts/rename_ldpc_matrices.sh` 转为 LF：
     ```bash
     dos2unix scripts/rename_ldpc_matrices.sh
     ```

2. **配置文件**
   - 复制示例：`cp scripts/auto_ldpc_local/sample_config.yaml scripts/auto_ldpc_local/my_config.yaml`
   - 根据环境修改：
     - `paths.output_root` 指向可写目录；
     - 如在 Windows 终端执行，可把 `system.command_prefix` 改成 `["wsl"]`；
     - 其他参数（SNR、top_n、err_inj 阈值等）与 `auto_ldpc` 相同。

3. **运行脚本**
   ```bash
   python scripts/auto_ldpc_local/automation.py \
       --config scripts/auto_ldpc_local/my_config.yaml
   ```
   - 仅执行一轮：追加 `--once`
   - 仅打印命令：追加 `--dry-run`

4. **流程与输出**
   - 自动执行矩阵生成 → AWGN → ERR\_INJ → 深挖；
   - 矩阵重命名仍使用 `scripts/rename_ldpc_matrices.sh`；
   - 结果写入 `output_root`：包含矩阵、日志、`state.db`、XLSX 报表和 SNR-FER 曲线图；
   - 配置中的 `limits` 与 `batch_size` 控制持续运行策略。

这样即可在本地完整复现 `auto_ldpc` 的自动化流程，便于离线调试。*** End Patch
