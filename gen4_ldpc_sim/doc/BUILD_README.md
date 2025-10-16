# 构建说明

## 概述

本项目现在支持构建**两套仿真程序**：

1. **原始版本** (`ssd_fc`) - 使用 `SSD_FC.cpp` + `ldpc_codec.cpp`
2. **DQ版本** (`ssd_fc_dq`) - 使用 `SSD_FC_dq.cpp` + `ldpc_codec dq.cpp`，带 `-DDQ_SIM` 宏定义

## 关键差异

### 编译宏定义
- **原始版本**：不定义特殊宏
- **DQ版本**：定义 `-DDQ_SIM` 宏，用于条件编译（如 `transceiver.cpp` 中的不同行为）

### 源文件
- **共享文件**：所有 `.c` 文件、`fc_dsp.cpp`、`transceiver.cpp`
- **原始版本特定**：`SSD_FC.cpp`、`ldpc_codec.cpp`
- **DQ版本特定**：`SSD_FC_dq.cpp`、`ldpc_codec dq.cpp`

## 编译命令

### 编译所有版本
```bash
make all
# 或者直接
make
```
这将编译生成两个可执行文件：`ssd_fc` 和 `ssd_fc_dq`

### 只编译原始版本
```bash
make ssd_fc
```

### 只编译DQ版本
```bash
make ssd_fc_dq
```

## 清理命令

### 清理所有编译文件
```bash
make clean
```

### 只清理原始版本
```bash
make clean-orig
```

### 只清理DQ版本
```bash
make clean-dq
```

## 运行命令

### 运行原始版本
```bash
make run
# 或者直接运行
./ssd_fc LDPC config/sdec.cnfg AWGN 5.30
```

### 运行DQ版本
```bash
make run-dq
# 或者直接运行
./ssd_fc_dq LDPC config/sdec.cnfg AWGN 5.30
```

## 目录结构

```
gen4_ldpc_sim/
├── build/          # 原始版本编译中间文件
├── build_dq/       # DQ版本编译中间文件
├── ssd_fc          # 原始版本可执行文件 (331KB)
├── ssd_fc_dq       # DQ版本可执行文件 (315KB)
├── src/
│   ├── SSD_FC.cpp          # 原始版本主文件
│   ├── SSD_FC_dq.cpp       # DQ版本主文件
│   ├── ldpc_codec.cpp      # 原始版本codec
│   ├── ldpc_codec dq.cpp   # DQ版本codec
│   └── ...                 # 其他共享源文件
├── config/
└── output/
```

## 技术细节

### 宏定义说明
在 Makefile 中：
- **原始版本** 使用：`CPPFLAGS = -I./src -MMD -MP -g`
- **DQ版本** 使用：`CPPFLAGS_DQ = -I./src -MMD -MP -g -DDQ_SIM`

`-DDQ_SIM` 等价于在代码中定义 `#define DQ_SIM`，使得以下代码生效：
```cpp
#ifdef DQ_SIM
    real_len = blk_len;  // DQ版本的行为
#else
    real_len = x;        // 原始版本的行为
#endif
```

### 文件名空格处理
`ldpc_codec dq.cpp` 文件名中包含空格，在 Makefile 中需要特殊处理：
- 使用反斜杠转义：`ldpc_codec\ dq.cpp`
- 或使用引号：`"ldpc_codec dq.cpp"`
- 目标文件名映射为：`ldpc_codec_dq.cpp.o`（下划线代替空格）

## 注意事项

1. 两个版本使用独立的编译目录（`build/` 和 `build_dq/`），互不干扰
2. 共享的源文件（如 `.c` 文件和其他 `.cpp` 文件）会分别编译到各自的目录，使用不同的编译选项
3. DQ版本的所有文件都带有 `-DDQ_SIM` 宏定义
4. 日志文件默认输出到 `perf/test.log` (原始版本) 和 `perf/test_dq.log` (DQ版本)
5. 编译成功后，两个可执行文件大小略有不同，这是正常的
