# ibex_matrix_family

当前目录下与 matrix family 相关的主流程只保留两段：

1. [`generate_family_artifacts.py`](/Users/roggerzwl/Project/MaxioProject/LDPC-codes/scripts/ibex_matrix_family/generate_family_artifacts.py)
   - 输入：`IBEX/ibex_matrix` 一类的源矩阵目录
   - 输出：
     - `matrix/`
     - `fade_matrix/`
     - `occupied_matrix/`
     - `ens/`
     - `bm_schematic/`
     - `lut/ldpc_matrix_assignments_maxk67.svh`
   - 不再导出：
     - `rdec_sched/`
     - 旧版 LUT (`ldpc_matrix_lut*.svh`)

2. [`generate_family_artifacts_from_occ_fade_rtl.py`](/Users/roggerzwl/Project/MaxioProject/LDPC-codes/scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py)
   - 输入：`ldpc_matrix_assignments_maxk67.svh`
   - 额外输入：固定 `wrap_base`
   - 语义：
     - 先按 `(M,K)` 裁剪 `occupied/fade`
     - 再基于 `wrap_base + delta` 从左到右重建 shift 链
   - 输出：
     - `ldpc_matrix_wrap_base_tmp_adj_ord`
     - `ldpc_matrix_wrap_base_deltasX_adj_ord`
     - `ldpc_matrix_first_mask_shiftMM`
     - `reconstructed_matrix/*.csv`

`RDEC_sched` 另由当前目录下现有的 `generate_cropped_rdec_sched*.py` 脚本生成，不属于这两段主流程的一部分。

## 第一段：源矩阵整理成统一 family 产物

脚本：

- [`generate_family_artifacts.py`](/Users/roggerzwl/Project/MaxioProject/LDPC-codes/scripts/ibex_matrix_family/generate_family_artifacts.py)

输入：

- `IBEX/ibex_matrix` 这类按码率分目录存放的源矩阵目录

输出：

- `matrix/`
- `fade_matrix/`
- `occupied_matrix/`
- `ens/`
- `bm_schematic/`
- `lut/ldpc_matrix_assignments_maxk67.svh`
- `manifest.json`

注意：

- `assignments` 是统一输出到一个文件：
  `lut/ldpc_matrix_assignments_maxk67.svh`
- 不是按不同 `M`、`K` 分成多个 assignments 文件
- 第一段不再导出：
  - `rdec_sched/`
  - 旧版 LUT

命令：

```bash
python3 scripts/ibex_matrix_family/generate_family_artifacts.py \
  --input-root IBEX/ibex_matrix \
  --output-root /tmp/ibex_family_stage1
```

## 第二段：从 assignments 生成参数

脚本：

- [`generate_family_artifacts_from_occ_fade_rtl.py`](/Users/roggerzwl/Project/MaxioProject/LDPC-codes/scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py)

输入：

- 第一段生成的统一 assignments 文件：
  `lut/ldpc_matrix_assignments_maxk67.svh`
- 固定 `wrap_base`

输出：

- `lut/ldpc_matrix_lut_fixed_wrap_base.svh`
- `lut/matrix_lut_summary_fixed_wrap_base.json`
- `lut/wrap_base_delta_table_fixed_wrap_base.csv`
- `lut/first_mask_shift_table_fixed_wrap_base.csv`
- `lut/reconstructed_matrix/*.csv`
- `manifest.json`

### `--wrap-base` 的输入格式

当前 `wrap_base` 通过命令行参数 `--wrap-base` 传入，格式是：

- 逗号分隔的整数串
- 可以带空格，也可以不带空格

等价示例：

```bash
--wrap-base "0,270,415,337,182,301,445,70,490,65,85,75,216"
```

```bash
--wrap-base "0, 270, 415, 337, 182, 301, 445, 70, 490, 65, 85, 75, 216"
```

规则：

- `--max-m 13` 时，至少需要 13 个元素
- `--max-m 17` 时，至少需要 17 个元素
- 如果不传 `--wrap-base`，默认使用 13 行版本：
  `0,270,415,337,182,301,445,70,490,65,85,75,216`

### 13 行视图示例

如果输入 assignments 来自 17 行 family，也可以只导出 13 行视图：

```bash
python3 scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py \
  --assignments-svh /tmp/ibex_family_stage1/lut/ldpc_matrix_assignments_maxk67.svh \
  --output-root /tmp/ibex_family_stage2 \
  --max-m 13
```

### 17 行视图示例

17 行视图必须显式传入 17 个 `wrap_base`：

```bash
python3 scripts/ibex_matrix_family/generate_family_artifacts_from_occ_fade_rtl.py \
  --assignments-svh /tmp/ibex_family_stage1/lut/ldpc_matrix_assignments_maxk67.svh \
  --output-root /tmp/ibex_family_stage2_17 \
  --max-m 17 \
  --wrap-base "b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16"
```

## 当前参数语义

- `wrap_base_deltasX_adj_ord` 的逻辑顺序是：
  - `M` 从小到大
  - 同一 `M` 内，`extra_user_data_col` 从小到大，也就是 `K=64,65,66,67`
- `wrap_base_deltasX_adj_ord` 每一行会补齐到相同长度：
  - `max-m=13` 时统一为 `6*36`
  - `max-m=17` 时统一为 `6*52`
  - 行尾不足部分补 `6'd0`
- `wrap_base_deltasX_adj_ord` 当前采用 Definition A：
  `min k >= 0, s.t. (last_element + k * delta) mod Z == 0`
- `first_mask_shiftMM` 会随 `max-m` 范围生成：
  - `max-m=13` 时生成到 `ldpc_matrix_first_mask_shift13`
  - `max-m=17` 时生成到 `ldpc_matrix_first_mask_shift17`
