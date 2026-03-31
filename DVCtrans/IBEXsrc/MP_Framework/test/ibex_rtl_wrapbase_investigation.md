# IBEX RTL-CN `wrap_base` 调查记录

## 1. 目的

这轮实验要回答 3 个问题：

1. 当前 `BF_IBEX_RTL_CN` 在 `M=10, K=67, ERR_INJ=200` 下到底能不能过。
2. “基于 `wrap_base` 从左到右重建矩阵”之后，新旧矩阵是否真的只是“每行一个常量偏移”。
3. 如果答案是“是”，那为什么还会出现“原矩阵编码出的码字，在新矩阵下不再合法”的现象。

这次实验不改仓库里的 `ldpc_codec.cpp`，只用临时文件 `/tmp/repro_mp_fer_ibex_rtl.cpp` 做验证。

---

## 2. 测试条件

- 译码器：`BF_IBEX_RTL_CN`
- 矩阵配置：`m=10, k=67, n=77, sc=512`
- 数据：全 1 payload
- 信道：`ERR_INJ`
- 输入形式：`MANUAL, sd_num=1, nand_strobes=1`
- FER 统计：`1000 packets`
- 调试统计：`1 packet`

补充说明：

- 这是当前 `ERR_INJ + 1bit/hard-like` 入口，不是多级 soft-input。
- 所以本文结论严格对应这条测试口径。

---

## 3. 相关语义

### 3.1 `WrapBase` 左到右重建

`IBEX/doc/IBEX_Matrix_LUT.md` 里给出的 RTL 语义是：

```text
S_r(n) = (WrapBase[r] + (n + 1) * delta[r]) mod Z
```

也就是：

- 每一行从 `WrapBase + delta` 开始
- 沿列方向每遇到一个 active CPM，就继续 `+delta`

这对应文档的：

- Section 4.2
- Section 4.3

### 3.2 `FirstMaskShift / mask`

同一份文档还说明了，当 `extra_bits_of_parity > 0` 时，最后一行不是完整的 `Z=512` parity 列，`mask/fade` 会把最后一行拆成：

- `occupied + mask==1` 的有效区
- `fade + mask==0` 的补区

文档对应：

- Section 6
- Section 11.1

关键点是：

```text
mask[col][k] = 1 iff (k - element[last_row][col]) mod Z < L
```

这里的 `L = extra_bits_of_parity`。

所以在 `M=10` 这类有 fractional parity 的配置下，最后一行的“行平移不变性”不是完全自由的。

---

## 4. 临时实验方法

临时 harness 在不改仓库代码的前提下做了 4 件事：

1. 从当前 `ldpc_config()` 读出原始 `h_matrix`
2. 生成两份重建矩阵：
   - `current_wrap_base`：用当前 `h_matrix.wrap_base` 左到右重建
   - `lut_wrap_base`：用 RTL 固定 LUT 的 `wrap_base` 左到右重建
3. 对两份重建矩阵都同步重建最后一行 `mask`
4. 用原矩阵编码出的码字，直接在重建矩阵下计算 syndrome

这一步是刻意这样设计的：

- 先不让“编码器自己的实现细节”干扰主结论
- 先直接回答“原码字在新矩阵下还是不是合法码字”

---

## 5. 编译命令

```bash
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I DVCtrans/IBEXsrc/MP_Framework/test/include \
  -I DVCtrans/IBEXsrc/MP_Framework \
  -I IBEX/src \
  -Dpafy_column=parity_column \
  /tmp/repro_mp_fer_ibex_rtl.cpp \
  DVCtrans/IBEXsrc/MP_Framework/test/build/ldpc_codec.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/ibex_rtl_engine.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/transceiver.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/alloc.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/finite_lib.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/intio.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/mod2convert.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/mod2dense.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/mod2sparse.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/rand.o \
  DVCtrans/IBEXsrc/MP_Framework/test/build/vec_op.o \
  -fsanitize=address,undefined \
  -o /tmp/repro_mp_fer_ibex_rtl
```

---

## 6. 实验 1：当前 `BF_IBEX_RTL_CN`，`ERR_INJ=200`

运行：

```bash
REPRO_PACKETS=1000 REPRO_ERR_BITS=200 ASAN_OPTIONS=detect_leaks=0 /tmp/repro_mp_fer_ibex_rtl
```

结果：

```text
[MP-IBEX-RTL-FER] decoder=6 rebuild=0 rebuild_lut=0 m=10 k=67 err_bits=200 packets=1000 fail_count=1000 fer=1 mismatch_frames=1000 avg_dec_mis=200 avg_raw_err=200
```

结论：

- 当前 `ibex_rtl` 在这条 case 下 `FER=1`
- `200 bit` 也完全过不了

---

## 7. 实验 2：用“当前 `h_matrix.wrap_base`”左到右重建

运行：

```bash
REPRO_PACKETS=1000 REPRO_ERR_BITS=200 REPRO_REBUILD_FROM_WRAP_BASE=1 \
ASAN_OPTIONS=detect_leaks=0 /tmp/repro_mp_fer_ibex_rtl
```

关键输出：

```text
[MP-IBEX-RTL-MATRIX] rebuilt_from_wrap_base diff_active_cpm=0 use_lut_wrap_base=0
[MP-IBEX-RTL-SYND] sw_original_under_original=0 sw_original_under_current=0 ...
[MP-IBEX-RTL-FER] decoder=6 rebuild=1 rebuild_lut=0 m=10 k=67 err_bits=200 packets=1000 fail_count=1000 fer=1 mismatch_frames=1000 avg_dec_mis=200 avg_raw_err=200
```

结论：

- 当前 `ldpc_config()` 产出的 `element`，本来就和当前 `h_matrix.wrap_base` 自洽
- 用当前 `wrap_base` 再左到右重建，不会改掉任何 active CPM
- 所以“当前矩阵不是按自身 `wrap_base` 左到右重建出来的”这个怀疑可以排掉
- 但即便这样，`ibex_rtl` 仍然是 `FER=1`

这说明：

- 当前 RTL-CN 失败，不是因为“当前 `h_matrix.element` 和当前 `wrap_base` 不一致”

---

## 8. 实验 3：用 LUT `wrap_base` 左到右重建，并同步重建 `mask`

### 8.1 先看结构关系

运行：

```bash
REPRO_MATRIX_ANALYSIS=1 REPRO_PACKETS=1 REPRO_ERR_BITS=0 \
ASAN_OPTIONS=detect_leaks=0 /tmp/repro_mp_fer_ibex_rtl
```

关键输出：

```text
[MP-IBEX-RTL-ROWSHIFT] tag=lut_wrap_base row=1 ... const_offset=1 offset=148 varying=0
[MP-IBEX-RTL-ROWSHIFT] tag=lut_wrap_base row=2 ... const_offset=1 offset=473 varying=0
[MP-IBEX-RTL-ROWSHIFT] tag=lut_wrap_base row=3 ... const_offset=1 offset=183 varying=0
...
[MP-IBEX-RTL-ROWSHIFT] tag=lut_wrap_base const_rows=10 total_rows=10
```

这说明：

- 用 LUT `wrap_base` 重建后，10 行全部都满足：
  - active CPM 支持集不变
  - 每一行的 shift 都只是加了一个常量偏移

也就是说，你说的这句判断是对的：

> “这种重建方法和直接读矩阵，最多在每行有一个偏移量”

从 `shift` 链本身看，结论确实如此。

### 8.2 但原码字在新矩阵下已经不是合法码字

同一轮输出里还有一条更关键：

```text
[MP-IBEX-RTL-SYND] sw_original_under_original=0 sw_original_under_current=0 sw_original_under_lut=448 mask_diff_current=0 mask_diff_lut=480 current_diff_active=0 lut_diff_active=268 original_fms=259 current_fms=259 lut_fms=259
```

这条的含义是：

- 原矩阵编码出的码字，在原矩阵下 syndrome = 0
- 原矩阵编码出的码字，在“当前 wrap_base 重建矩阵”下 syndrome = 0
- 原矩阵编码出的码字，在“LUT wrap_base 重建矩阵”下 syndrome = 448

同时：

- `mask_diff_lut = 480`
- 说明同步重建 `mask` 以后，最后一行的有效窗口确实发生了变化

这说明：

- 虽然 `shift` 层面是“每行常量偏移”
- 但在 `M=10` 这类 `extra_bits_of_parity > 0` 的配置里，整套 H 语义并不只由“行内相对 shift 差分”决定
- 最后一行的 `occupied/fade + mask` 还绑定了绝对 CN 索引
- 因此这里的“行变换”不是自由的等价变换

换句话说：

- “每行常量偏移”是真的
- “所以码空间一定不变”在这条配置上不成立

---

## 9. 实验 4：LUT `wrap_base` 重建矩阵下，`0 error` 是否通过

运行：

```bash
REPRO_PACKETS=1 REPRO_ERR_BITS=0 REPRO_DUMP_FIRST=1 REPRO_REBUILD_FROM_LUT_WRAP_BASE=1 \
ASAN_OPTIONS=detect_leaks=0 /tmp/repro_mp_fer_ibex_rtl
```

结果：

```text
[MP-IBEX-RTL-DBG] raw_err=0 dec_mis=0 cw_fail=1 decoder=6 init_sw=448 fina_sw=448 cnvg_itr=49 cnvg_lyr=0 iter_out=50 early_term=0 sw_before=448 sw_after=448 col_cnt=-1
[MP-IBEX-RTL-FER] decoder=6 rebuild=1 rebuild_lut=1 m=10 k=67 err_bits=0 packets=1 fail_count=1 fer=1 mismatch_frames=0 avg_dec_mis=0 avg_raw_err=0
```

结论：

- 即使 `0 error`
- 只要把 decoder 侧矩阵切成“LUT wrap_base + 左到右重建 + 同步 mask”
- 当前原矩阵编码出的码字就已经不再满足这个新矩阵

这和上面的 syndrome=448 是完全一致的。

---

## 10. 实验 5：LUT `wrap_base` 重建矩阵下，`ERR_INJ=200`

运行：

```bash
REPRO_PACKETS=1000 REPRO_ERR_BITS=200 REPRO_REBUILD_FROM_LUT_WRAP_BASE=1 \
ASAN_OPTIONS=detect_leaks=0 /tmp/repro_mp_fer_ibex_rtl
```

结果：

```text
[MP-IBEX-RTL-FER] decoder=6 rebuild=1 rebuild_lut=1 m=10 k=67 err_bits=200 packets=1000 fail_count=1000 fer=1 mismatch_frames=1000 avg_dec_mis=200 avg_raw_err=200
```

结论：

- 这条“基于 LUT `wrap_base` 左到右重建矩阵 + ibex_rtl”的临时方案也不能过
- 而且失败并不是“200bit 太难”
- 更直接的原因是：原码字在这个重建矩阵下本来就不是合法码字

---

## 11. 最终结论

### 11.1 当前 `ibex_rtl` 在 `M=10, ERR_INJ=200` 下确实过不了

结果明确是：

```text
FER = 1
```

### 11.2 “当前矩阵没有按自身 `wrap_base` 左到右重建”不是根因

因为：

- 用当前 `h_matrix.wrap_base` 重建后 `diff_active_cpm=0`
- syndrome 仍然是 `0`

说明当前 `ldpc_config()` 产出的矩阵，本来就和它自己的 `wrap_base` 一致。

### 11.3 LUT `wrap_base` 重建后，shift 确实只是“每行常量偏移”

这个判断也是对的：

- 10 行全部满足 `const_offset=1`
- 支持集不变
- 每行只是整体平移

### 11.4 但在 `M=10` 这类 fractional parity 配置里，这种“行平移”不再保持码空间

原因是：

- `extra_bits_of_parity > 0`
- 最后一行 `occupied/fade + mask` 把绝对 CN 索引绑死了
- 文档 `IBEX_Matrix_LUT.md` Section 11.1 已经说明，这时“行平移不变性”不完全成立

实验上对应的直接证据就是：

```text
sw_original_under_original = 0
sw_original_under_lut      = 448
```

也就是说：

- 用 LUT `wrap_base` 左到右重建出的矩阵，不只是“与原矩阵等价的行变换”
- 至少在当前 `M=10, K=67` 这条配置上，它已经改变了码空间语义

### 11.5 因此，这条线的结论不是“RTL_CN 天生不能过”

更准确的说法是：

- 当前 `BF_IBEX_RTL_CN` 这条实现本身在基线下就 `FER=1`
- 另外，如果再把矩阵强行切到 LUT `wrap_base` 语义，问题不会变好，反而会先变成“矩阵与原码字不匹配”

---

## 12. 本轮结论的边界

这轮实验只回答了：

- 当前基线 `ibex_rtl` 为什么 `200bit` 也过不了
- “LUT `wrap_base` 左到右重建”是不是一个可行修复方向

结论是：

- 不是一个可直接成立的修复方向

但这轮还没有继续证明：

- 当前基线 `BF_IBEX_RTL_CN` 的根因究竟更偏向
  - `wrap_base_delta`
  - `CN row rotation`
  - `mask` 使用时机
  - 还是 `sd_num=1` 这条 hard-like 入口本身

如果继续往下查，下一步建议直接盯：

1. `ibex_rtl_engine.cpp` 里的 `wrap_base_delta` 和末尾 wrap 旋转
2. `build_rtl_view()` 对 LUT `wrap_base` / `h_matrix.element` 的混合使用
3. `BF_IBEX_RTL_CN` 在 `ERR_INJ + sd_num=1` 下是否本来就没有有效 correction gain

---

## 13. 后续补充实验：强制重建 `rtl_view.shift/mask`

在 `ibex_rtl_engine.cpp` 里额外加了两个临时调试开关：

- `IBEX_RTL_FORCE_REBUILT_SHIFT=1`
- `IBEX_RTL_FORCE_REBUILT_MASK=1`

语义是：

- `shift` 不再优先取 `h_matrix.element`
- 而是强制使用 `wrap_base + delta` 左到右重建值
- `mask` 也可以选择不再覆盖回 `h_matrix.mask`

### 13.1 `ERR_INJ=200, 1000 packets`

三组结果如下：

```text
baseline:
[MP-IBEX-RTL-FER] ... fail_count=1000 fer=1 avg_dec_mis=200

force rebuilt shift:
[MP-IBEX-RTL-FER] ... fail_count=1000 fer=1 avg_dec_mis=200

force rebuilt shift + mask:
[MP-IBEX-RTL-FER] ... fail_count=1000 fer=1 avg_dec_mis=200
```

结论：

- 只修 `rtl_view.shift`
- 或同时修 `rtl_view.shift + rtl_view.mask`

都不足以让 `BF_IBEX_RTL_CN` 恢复纠错。

### 13.2 单包调试结果

`ERR_INJ=200, 1 packet` 时：

```text
baseline:                    init_sw=708 fina_sw=708
force rebuilt shift:         init_sw=692 fina_sw=692
force rebuilt shift + mask:  init_sw=681 fina_sw=681
```

这说明：

- 强制重建 `shift/mask` 的确会影响初始 syndrome 权重
- 也就是说，它不是“完全没起作用”
- 但 50 轮之后 `fina_sw` 仍然和 `init_sw` 一样，说明迭代里仍然没有形成有效 correction

### 13.3 `0 error` 边界测试

`ERR_INJ=0, 1 packet` 时，上面三组都能过：

```text
init_sw=0 fina_sw=0 cw_fail=0
```

这点很重要，它说明：

- 强制重建 `rtl_view.shift/mask` 并没有把无错码字本身搞坏
- 当前问题更像是“有错时的迭代更新链不工作”
- 而不是“静态视图一上来就不自洽”

### 13.4 普通 `BF_IBEX` 参考线

同样口径下我又跑了普通 `BF_IBEX`：

```text
[MP-IBEX-RTL-FER] decoder=5 ... err_bits=200 packets=100 fail_count=100 fer=1 avg_dec_mis=200
```

这说明：

- 当前问题已经不只是在 `RTL_CN` 的视图重建上
- 至少在你现在这条 `ERR_INJ + nand_strobes=1` 入口下
- 普通 `BF_IBEX` 和 `BF_IBEX_RTL_CN` 都没有产生有效纠错

所以这轮新增实验把结论又往前收了一步：

- “`view->shift` 应该重建”这个判断是对的
- 但它不是当前 `FER=1` 的唯一主因
- 真正还要继续追的是 BF 共享的 likelihood / flip / syndrome 更新链
