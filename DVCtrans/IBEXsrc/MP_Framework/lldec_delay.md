# LLDEC RTL 延迟链路详细分析
## 1. 概述
本文档详细列出 RTL 中 lldec 模块从 syndrome 更新到 syndrome_weight 输出的完整延迟链路,
重点分析 `synd_vn` 寄存器 - `lldec_core_cn_syndw` 输出之间的流水线延迟,
以及 WRAP 周期对迭代边界处 pipeline 时序的影响.
## 2. 延迟链路涉及的模块层次
lldec_cnode.v
- ← syndrome 寄存器(synd_vn)+ bit_sum 实例化
- - lldec_bit_sum_512_x_13.v-← 13行聚合+ 第4级寄存器
lldec_bit_sum_512.v-← 单行512bit popcount,3级寄存器
lldec_sch.v
- ← 使用 lldec_core_cn_syndw 进行 post_process 判断
lldec_likelihood_update.v ← 使用 lldec_core_flip_post_process 进行 likelihood 更新
## 3. 第一级:lldec_cnode.v - synd_vn 寄存器更新
### 3.1 synd_vn 寄存器定义
**文件**:`design/rtl/ldec/dec_core/lldec/lldec_cnode.v`第197行
```verilog
reg [LDPC_P_SIZE*LDPC_M_MAX-1:0] synd_vn;
```
### 3.2 synd_vn 更新逻辑(第210-227行)
```verilog
generate
for (i=0;i<LDPC_M_MAX;i=i+1) begin: gen_syndrome
always @(posedge clk or negedge rst_n)
begin
if (!rst_n) begin
synd_vn[i*LDPC_P_SIZE +: LDPC_P_SIZE] <= {LDPC_P_SIZE{1'd0}};
end
else if(lldec_core_cn_synd_clr_req & lldec_core_cn_synd_clr_ack) begin
synd_vn[i*LDPC_P_SIZE +: LDPC_P_SIZE] <= {LDPC_P_SIZE{1'd0}};
end
// DECODE 周期更新: flip_data XOR synd_vn_delta_shift
else if(lldec_core_cn_synd_upd_en & lldec_core_cn_synd_upd_row_enable[i]
& lldec_core_flip2cn_data_vld & lldec_core_flip2cn_data_rdy) begin
synd_vn[i*LDPC_P_SIZE +:LDPC_P_SIZE]<=
(lldec_core_flip2cn_data & lldec_core_cn_synd_upd_mask_for_row[i*LDPC_P_SIZE +: LDPC_P_SIZE])
^ synd_vn_delta_shift[i*LDPC_P_SIZE +: LDPC_P_SIZE];
end
// WRAP 周期更新:wrap_shift
else if(lldec_core_cn_synd_wrap_req && lldec_core_cn_synd_wrap_ack) begin
synd_vn[i*LDPC_P_SIZE +: LDPC_P_SIZE] <= synd_vn_wrap_shift[i*LDPC_P_SIZE +: LDPC_P_SIZE];
end
end
end
endgenerate
```
### 3.3 synd_vn_delta_shift -组合逻辑(delta_shifter 输出)
**文件**:`lldec_cnode.v`第230-253行
```verilog
generate
for (i=0;i<LDPC_M_MAX;i=i+1) begin:gen_shift
lldec_delta_shifter u_lldec_delta_shifter(
.dout(synd_vn_delta_shift[i*LDPC_P_SIZE +:LDPC_P_SIZE]),// 输出:组合逻辑
.din(synd_vn[i*LDPC_P_SIZE +: LDPC_P_SIZE]),
//输入:synd_vn 寄存器
.delta(lldec_core_cn_synd_upd_matrix_delta[i*W_SHIFT_DELTA +: W_SHIFT_DELTA])

end
endgenerate
```
**关键**:`synd_vn_delta_shift`是 `synd_vn` 经过 delta_shifter 的组合逻辑输出,
没有额外寄存器延迟.
## 3.4 bit_sum 输入连接
**文件**:`lldec_cnode.v`第366-398行(P_SIZE=512 路径)
```verilog
// P_SIZE=512的实例化
lldec_bit_sum_512_x_13 u_lldec_bit_sum_512_x_13(
.o_hamming_weight(lldec_core_cn_syndw),
.o_hamming_weight_valid(),
.o_hamming_weight_valid_col_end(lldec_core_cn_syndw_col_end_vld),
.clk(clk),
.rst_n(rst_n),
.i_en(1'd1),
//始终使能
.i_cn(synd_vn),
// 直连 synd_vn 寄存器输出!
.i_col_end(lldec_core_cn_synd_upd_col_end)
```
**重要**:`.i_cn(synd_vn)` -pipeline 输入直接连接 `synd_vn` 寄存器输出,不是 `synd_vn_delta_shift`!
## 3.5 lldec_core_cn_syndw 输出声明
**文件**:`lldec_cnode.v`第155行
```verilog
output [DNEINF_SYNDW_WID-1:0] lldec_core_cn_syndw; // 4 cycle delay

```
## 4. 第二级:lldec_bit_sum_512_x_13.v-13行聚合+第4级寄存器
### 4.1 模块声明
**文件**:`design/rtl/ldec/dec_core/lldec/lldec_bit_sum_512_x_13.v`
```verilog
module lldec_bit_sum_512_x_13(
output [13:0] o_hamming_weight,
output
o_hamming_weight_valid,
output
o_hamming_weight_valid_col_end,
input
clk, rst_n, i_en, i_col_end,
input [512*LDPC_M_MAX-1:0] i_cn
);

// Description: Computes the sum of 13 512-bit syndromes. Latency: 4
```
### 4.2 内部13个 lldec_bit_sum_512实例
**文件**:`lldec_bit_sum_512_x_13.v`第38-52行
```verilog

generate
for(g = 0; g < LDPC_M_MAX; g = g + 1) begin: weight
lldec_bit_sum_512 u_lldec_bit_sum_512(
.clk(clk),
.rst_n(rst_n),
.en(1'd1),
//始终使能
.i_data(i_cn[g*512 +:512]),//每行512bit输入
.o_hamming_weight(hamming_weight[g]),
.o_hamming_weight_vld(),
.col_end_en(1'd1)

end

endgenerate
```
### 4.3 13行权重求和一组合逻辑
**文件**:`lldec_bit_sum_512_x_13.v` 第54-70行
```verilog
assign hamming_weight_total = hamming_weight[0] + hamming_weight[1] + ... + hamming_weight[16];
```
**注意**:`hamming_weight_total`是组合逻辑输出,直接由17个`hamming_weight[g]`相加.
### 4.4 第4级寄存器一最终输出锁存
**文件**:`lldec_bit_sum_512_x_13.v`第72-89行
```verilog
always @(posedge clk or negedge rst_n) begin
if (!rst_n)
hamming_weight_total_r <= 0;
else if (hamming_weight_valid_r[2]) hamming_weight_total_r <= hamming_weight_total;
end
always @(posedge clk or negedge rst_n) begin
if (!rst_n) hamming_weight_valid_r <= 0;
else
hamming_weight_valid_r <= {hamming_weight_valid_r[2:0], i_en};
end
assign o_hamming_weight = hamming_weight_total_r;
assign o_hamming_weight_valid = hamming_weight_valid_r[3];
assign o_hamming_weight_valid_col_end = hamming_weight_valid_r col_end[4];
```
**关键**:
- `hamming_weight_total_r`在 `hamming_weight_valid_r[2]`为高时锁存
- `hamming_weight_valid_r`是 4bit 移位寄存器:`{r[2:0],i_en}`
-有效信号延迟4 拍:`i_en→r[0]→r[1]→r[2]→r[3]`
-数据在第3拍(`r[2]`)锁存到`total_r`,第4拍通过`r[3]` 输出有效
## 5. 第三级:lldec_bit_sum_512.v -单行512bit popcount,3级寄存器
### 5.1 模块声明
**文件**:`design/rtl/ldec/dec_core/lldec/lldec_bit_sum_512.v`
```verilog
module lldec_bit_sum_512(
input
clk, rst_n, en, col_end_en,
input [511:0] i_data,
output [9:0] o_hamming_weight,
output
o_hamming_weight_vld,
output
o_hamming_weight_vld_col_end

/// Description: The sum of 512 bit data in 3 clock cycles
```
### 5.2 第1级寄存器(posedge clk)- s04→s05
**文件**:`lldec_bit_sum_512.v`第190-224行
```verilog
always @(posedge clk or negedge rst_n) begin
if (rst_n == 0) begin
s0500<=0;s0501<=0;... s05_05<=0;
hamming_weight_vld[0] <= 1'd0;
hamming_weight_vld_col_end[0] <= 1'd0;
end
else begin
//组合逻辑计算 s04 + popcount →s05
for (i=0;i<2;i=i+1)
{s05_01[i],s05_00[i]} <= s0400[3*i] + s04_00[3*i+1]+ s04_00[3*i+2];
//...更多并行加法...
hamming_weight_vld[0] <= en;
hamming_weight_vld_col_end[0] <= col_end_en;
end
end
```
### 5.3 第2级寄存器(posedge clk)- s09 →s10
**文件**:`lldec_bit_sum_512.v`第304-332行
```verilog
always @(posedge clk or negedge rst_n) begin
if (rst_n == 0) begin
s1000<= 0;s1001<=0;... s10_08 <=0;
hamming_weight_vld[1] <= 1'd0;
hamming_weight_vld_col_end[1] <= 1'd0;
end
else begin
s10_00 <= s09_00;
s10_01 <= s09_01;
//...更多传递和加法...
hamming_weight_vld[1] <= hamming_weight_vld[0];
hamming_weight_vld_col_end[1] <= hamming_weight_vld_col_end[0];
end
end
```
### 5.4 第3级寄存器(posedge clk)- s15 →hamming_weight
**文件**:`lldec_bit_sum_512.v`第403-423行
```verilog
always @(posedge clk or negedge rst_n) begin
if (rst_n == 0) begin
hamming_weight <= 0;
hamming_weight_vld[2] <= 1'd0;
hamming_weight_vld_col_end[2] <= 1'd0;
end
else begin
hamming_weight[0] <= s15_00;
hamming_weight[1] <= s15_01;
//..更多传递..
hamming_weight[9] <= s15_09[0]+ s15_09[1];
hamming_weight_vld[2] <= hamming_weight_vld[1];
hamming_weight_vld_col_end[2] <= hamming_weight_vld_col_end[1];
end
end

```
### 5.5 输出
```verilog
assign o_hamming_weight = hamming_weight;

assign o_hamming_weight_vld = hamming_weight_vld[2];

```
## 6. 完整4级流水线延迟总结
| 级别 | 模块 | 寄存器 | 时钟周期 | 说明 |
|---|---|---|---|---|
| Stage 0 | lldec_cnode | synd_vn | T+0 | synd_vn 在 posedge clk 更新 (DECODE 周期: flip XOR delta_shift; WRAP 周期: wrap_shift) |
| Stage 1 | lldec_bit_sum_512 | s05_xxx, vld[0] | T+1 | pipeline 第1级锁存:捕获 synd_vn 的组合逻辑 popcount 结果 |
| Stage 2 | lldec_bit_sum_512 | s10_xxx, vld[1] | T+2 | pipeline 第2级锁存|
| Stage 3 | lldec_bit_sum_512 | hamming_weight, vld[2] | T+3 | pipeline 第3级锁存:输出 per-row hamming_weight |
| Stage 4 | lldec_bit_sum_512_x_13 | hamming_weight_total_r, vld[3] | T+4 | 13行求和锁存:输出 lldec_core_cn_syndw |
**总延迟**:从 `i_cn = synd_vn`输入到 `o_hamming_weight = lldec_core_cn_syndw` 输出=**4 个时钟周期**

## 7. 时钟周期级精确时序分析
## 7.1 speed_mode_full(mode 0)下的时序
在 speed_mode_full 下,每个时钟周期处理一列(col), `col_itr_en`每周期拉高.
## 7.2 关键时序关系
时钟周期:
T
T+1
T+2
T+3
T+4
DECODE状态:
col_itr_en: col=J col=J+1 col=J+2 col=J+3 col=J+4
synd_vn:
[J更新][J+1更新][J+2更新][J+3更新][J+4更新]
(posedge)


bit_sum Stage1:-[捕获J][捕获J+1][捕获J+2][捕获J+3]
bit_sum Stage2:---[传递J][传递J+1][传递J+2]
bit_sum Stage3: --
[输出J][输出J+1]
bit_sum_x_13: ---

[锁存J→syndw=SW(J)]
结论:lldec_core_cn_syndw 在 col=J+4 时输出 SW(J)
即:syndw(K)=SW(K-4)(K>=4,同一 iteration内)
## 7.3 Pipeline 输入:捕获的是 synd_vn 的哪个值?
**关键细节**:`i_cn = synd_vn` 是 wire 直连.
在 posedge clk T 时刻:
1.`synd_vn`被**更新**为 col J的 flip XOR结果(非阻塞赋值,T时刻后生效)
2.同时,bit_sum Stage1 捕获的 `i_cn`是 **T-1 时刻** 的 `synd_vn`值
**因此**:
- 在 col J 的 posedge clk, Stage1 捕获的是 col J-1 更新后的 `synd_vn`
- 在 col J+4 的 posedge clk, Stage4 输出的是 col J-1 更新后的 `synd_vn`对应的 SW
等价地说:
- **syndw 在 col K 输出 = SW(col K-4 对应的 synd_vn 状态)**
- 即 **syndw(K) = SW(K-4)**,其中 SW(K-4) 是 col K-4 处理完毕后 synd_vn 对应的 syndrome_weight

## 8. WRAP 周期对Pipeline 的影响
### 8.1 FSM状态转换
**文件**:`design/rtl/ldec/dec_core/lldec/lldec_sch.v`第894-906行
```verilog
ST_SCH_WRAP:
begin
if(!lldec_max_itr_awon && lldec_core_cn_syndw is zero)
sch_fsm_p = ST_SCH_WAIT_SYNDW;
else if(cmd_type === SOFT_CMD && (itr_cnt==0) && lldec_sd_ll_mode & lldec_core_cn_syndw_col_end_vld)
sch_fsm_p = ST_SCH_INIT_LL;
else if(cmd_type == SOFT_CMD && (itr_cnt==0) & lldec_sd_ll_mode && !lldec_core_cn_syndw_col_end_vld)
sch_fsm_p = ST_SCH_WRAP;
else if (itr_cnt == 0)
sch_fsm_p = ST_SCH_INIT_LL;
else
sch_fsm_p = ST_SCH_DECODE;
//←正常路径:WRAP →DECODE
end
```
## 8.2 迭代间WRAP周期的时序
当 `one_itr_end & speed_mode_en`时,FSM 进入 WRAP 状态(1个时钟周期),
然后转回 DECODE 状态.
迭代N 最后一列(col=N_COL-1):
时钟 T:DECODE,col=N_COL-1,synd_vn 更新
时钟 T+1:WRAP, synd_vn 被 wrap_shift 更新←1个额外时钟周期
迭代N+1的前几列:
时钟 T+2:DECODE, col=0, synd_vn 更新(col 0 的 flip XOR)
时钟 T+3:DECODE,col=1,synd_vn 更新
时钟 T+4:DECODE,col=2,synd_vn 更新
时钟 T+5:DECODE,col=3,synd_vn 更新
## 8.3 WRAP 期间pipeline的行为
**WRAP 周期(T+1)发生什么?**
- `lldec_core_cn_synd_upd_en = 0`(因为 sch_fsm != ST_SCH_DECODE,参见第1474行)
-但是`i_en=1'd1`(在 bit_sum_512_x_13实例化中硬连接为1'd1)
-所以 pipeline **不会停止**,它继续处理!
WRAP 周期的 pipeline 行为:
时钟T+1(WRAP):
- synd_vn 被 wrap_shift 更新(WRAP 操作)
- bit_sum Stage1 捕获:时钟 T 处理 col=N_COL-1 后的 synd_vn(DECODE 最后一次更新)
- bit_sum Stage2:传递上一级的数据
-bit_sum Stage3:传递上上级的数据
- bit_sum_x_13 Stage4: 锁存并输出 →syndw = SW(N_COL-4 时的 synd_vn)
时钟T+2(新迭代 col=0):
- synd_vn 被 col=0 的 flip XOR 更新
- bit_sum Stage1 捕获: WRAP后的 synd_vn (wrap_shift 的结果)
时钟T+3(新迭代 col=1):
- synd_vn 被 col=1 的 flip XOR 更新
- bit_sum Stage1 捕获: col=0 更新后的 synd_vn
- bit_sum Stage2:传递 WRAP 后 synd_vn 对应的 hamming_weight
时钟T+4(新迭代 col=2):
- synd_vn 被 col=2 的 flip XOR 更新
- bit_sum Stage1 捕获: col=1 更新后的 synd_vn
- bit_sum Stage2:传递 col=0 更新后的结果
时钟T+5(新迭代col=3):
- synd_vn 被 col=3 的 flip XOR 更新
- bit_sum Stage1 捕获: col=2 更新后的 synd_vn
- bit_sum Stage2:传递 col=1 更新后的结果
- bit_sum Stage3:传递 col=0 更新后的结果
...
时钟T+6(新迭代col=4):
- synd_vn 被 col=4 的 flip XOR 更新
- bit_sum Stage1 捕获: col=3 更新后的 synd_vn
- bit_sum Stage2:传递 col=2 更新后的结果
- bit_sum Stage3:传递 col=1 更新后的结果
- bit_sum_x_13 Stage4: 锁存并输出 →syndw = SW(WRAP后的 synd_vn) ←不是 SW(col -4)!
### 8.4 关键结论:WRAP周期使得迭代边界处延迟变为5拍
由于 WRAP 周期占据1个时钟周期,且此期间 pipeline 不停 (i_en=1),
pipeline 中多"消化"了一个 WRAP 后的 synd_vn 状态:
| 迭代 N+1 的 col | syndw 输出对应的 synd_vn 状态 | 延迟拍数 |
|---|---|---|
| col 0 | SW(迭代N col=N_COL-4) | - |
| col 1 | SW(迭代N col=N_COL-3) | - |
| col 2 | SW(迭代N col=N_COL-2) | - |
| col 3 | SW(迭代N col=N_COL-1) | - |
| **col 4** | **SW(WRAP后,即迭代N+1 col=0之前的 wrap 状态)** | **5** |
| col 5 | SW(迭代N+1 col=0)| 5 |
| col 6 | SW(迭代N+1 col=1)| 5 |
| ... | ... | ... |
| col N_COL-1 | SW(迭代N+1 col=N_COL-5) | 4 |
**等等,需要重新分析**.让我们更精确地追踪:
## 9. 精确时钟周期追踪表
### 9.1 假设
- 迭代N 最后一列 col=72 (data len=73, col 0~72)
- speed_mode full,每时钟一列
- Pipeline 始终运行(i_en=1'd1 硬连接)
### 9.2 迭代N-N+1边界处的逐周期追踪
| 周期编号 | FSM状态 | col | synd_vn 内容(posedge 后) | syndw 输出(posedge 后) |
|---|---|---|---|---|
| T-3 | DECODE | 69 | synd_vn[69] | SW(65) |
| T-2 | DECODE | 70 | synd_vn[70] | SW(66) |
| T-1 | DECODE | 71 | synd_vn[71] | SW(67) |
| T | DECODE | 72 | synd_vn[72](最后一列) | SW(68) |
| T+1 | WRAP | -- | synd_vn_wrap(rotate) | SW(69) |
| T+2 | DECODE | 0 | synd_vn[0](新iter) | SW(70) |
| T+3 | DECODE | 1 | synd_vn[1] | SW(71) |
| T+4 | DECODE | 2 | synd_vn[2] | SW(72) |
| T+5 | DECODE | 3 | synd_vn[3] | SW(WRAP后)←关键 |
| T+6 | DECODE | 4 | synd_vn[4] | SW(0)←关键 |
| T+7 | DECODE | 5 | synd_vn[5] | SW(1) |
## 9.3 关键分析
- **周期 T+1(WRAP)**: `synd_vn` 被 wrap_shift 更新.pipeline Stage1 捕获的是 T 时刻(col=72 处理后) 的 synd_vn.syndw 输出 SW(69).
- **周期 T+2(新iter col=0)**: pipeline Stage1 捕获 WRAP 后的 synd_vn.syndw 输出 SW(70).
- **周期 T+3(新iter col=1)**: Stage1 捕获 col=0 更新后的 synd_vn.syndw 输出 SW(71).
**周期 T+4(新iter col=2)**: Stage1 捕获 col=1 更新后的 synd_vn.syndw 输出 SW(72).
- **周期 T+5(新iter col=3)**: Stage1 捕获 col=2 更新后的 synd_vn.syndw 输出 **SW(WRAP后)**.
- **周期 T+6(新iter col=4)**: Stage1 捕获 col=3 更新后的 synd_vn.syndw 输出 **SW(col=0)**.
## 9.4 syndw 与当前 col 的关系总结

| 当前处理 col | syndw 输出对应 | 延迟拍数(相对当前 col) | 说明 |
|---|---|---|---|
| 新iter col=0 | SW(旧iter col=70) | 4+2=6? | 不对.. |
**等等,让我重新精确分析**.
syndw 在 posedge clk 输出,与 synd_vn 在 posedge clk 更新是同时的.

所以在 col=0的 posedge clk 时,syndw 输出的值是 *上一拍* Stage4 锁存的结果.
让我重新按"每个 posedge clk 瞬间各信号状态"来追踪:

```text
posedge clk 时刻 | FSM    | col_itr | synd_vn 变为   | Stage1 捕获             | Stage4 输出(syndw)
=================|========|=========|================|=========================|====================
T                | DECODE | col=69  | synd_vn[69]    | synd_vn[68](T-1 的)     | SW(65)
T+1              | DECODE | col=70  | synd_vn[70]    | synd_vn[69](T 的)       | SW(66)
T+2              | DECODE | col=71  | synd_vn[71]    | synd_vn[70](T+1 的)     | SW(67)
T+3              | DECODE | col=72  | synd_vn[72]    | synd_vn[71](T+2 的)     | SW(68)
T+4              | WRAP   | --      | synd_vn_wrap   | synd_vn[72](T+3 的)     | SW(69)
T+5              | DECODE | col=0   | synd_vn[0]     | synd_vn_wrap(T+4 的)    | SW(70)
T+6              | DECODE | col=1   | synd_vn[1]     | synd_vn[0](T+5 的)      | SW(71)
T+7              | DECODE | col=2   | synd_vn[2]     | synd_vn[1](T+6 的)      | SW(72)
T+8              | DECODE | col=3   | synd_vn[3]     | synd_vn[2](T+7 的)      | SW(WRAP)
T+9              | DECODE | col=4   | synd_vn[4]     | synd_vn[3](T+8 的)      | SW(col=0)
T+10             | DECODE | col=5   | synd_vn[5]     | synd_vn[4](T+9 的)      | SW(col=1)
```

## 9.5 正确的延迟对应关系
从上表可以清楚地看到:

```text
| 当前 col       | syndw 对应         | 延迟(时钟周期) | 相当于             |
|:---------------|:-------------------|:---------------|:-------------------|
| col=0(新iter) | SW(旧iter col=70)  | 4+WRAP偏移     | 旧iter col=70      |
| col=1          | SW(旧iter col=71)  | 4+WRAP偏移     | 旧iter col=71      |
| col=2          | SW(旧iter col=72)  | 4+WRAP偏移     | 旧iter col=72      |
| col=3          | SW(WRAP后)         | 4+WRAP偏移     | WRAP后、col=0前    |
| col=4          | SW(col=0)          | 5              | 新iter col=0       |
| col=5          | SW(col=1)          | 5              | 新iter col=1       |
| ...            | ...                | ...            | ...                |
| col=72         | SW(col=67)         | 5              | 新iter col=67      |
```

**但这里有个关键问题**:在col=0~2 时,syndw 对应的是旧迭代末尾的 SW,
这些 SW 值在 C 模型中来自上一迭代的 syndrome_weight_r[]移位寄存器.
## 10. lldec_sch.v -lldec_core_cn_syndw 的使用
## 10.1 post_process信号生成
**文件**:`design/rtl/ldec/dec_core/lldec/lldec_sch.v`第1340-1345行
```verilog
always @* begin
for(ii = 0; ii < LDPC_P_SIZE; ii = ii + 1) begin
lldec_core_flip_post_process[ii] = (st_decode &&(lldec_core_cn_syndw < lldec_post_weight_thr0)
&& (itr_cnt[3:0] < lldec_post_ratio))
& (itr_cnt >= lldec_post_trigger) && lldec_post_process_en
? prng[ii] : 1'd0;
lldec_core_flip_post_process[LDPC_P_SIZE + ii] = (st_decode && (lldec_core_cn_syndw < lldec_post_weight_thr1)
&& (itr_cnt[3:0] >= lldec_post_ratio))
& (itr_cnt >= lldec_post_trigger) && lldec_post_process_en
? prng[ii] : 1'd0;
end
end
```
**关键**:
- `st_decode = (sch_fsm == ST_SCH_DECODE)` - 只在 DECODE 状态下 post_process 生效
- `lldec_core_cn_syndw`是组合逻辑直接使用(`always @*`),没有额外寄存器
### 10.2 post_process条件总结
post_process[0] = st_decode
&&(syndw < synd thr_qc) // synd_thr_qc= 16
&& (itr_cnt[3:0] < post_ratio) // post_ratio = 12
&& (itr_cnt >= post_trigger) // post_trigger = post iter = 512
&& post_process_en
post_process[1] = st_decode
&&(syndw < synd_thr_post) // synd_thr_post = 48
&(itr_cnt[3:0] >= post_ratio)
&& (itr_cnt >= post_trigger)
& post_process_en
### 10.3 WRAP 周期时 post_process = 0
WRAP 周期中 `sch_fsm = ST_SCH_WRAP`,所以 `st_decode = 0`,
`lldec_core_flip_post_process`全部为 0.**post_process 在 WRAP 周期不生效**.
### 10.4 PRNG
更新
**文件**:`lldec_sch.v`第1287-1289行
```verilog
assign prng_init = (sch_fsm != ST_SCH_WRAP) && (sch_fsm_p == ST_SCH_WRAP)
&(itr_cnt == lldec_post_trigger - 1'd1);
assign prng_update = (sch_fsm == ST_SCH_DECODE) && col_itr_en && speed_mode_en;
assign prng_seed = 512'h1fe0_1fe0_..._1fe0;
```
**关键**:
- `prng_init`:在 WRAP-DECODE 转变且 `itr_cnt == post_trigger - 1`时初始化
- `prng_update`:只在 DECODE 状态且 `col_itr_en && speed_mode_en`时更新
- **WRAP 周期 PRNG 不更新**(sch_fsm != ST_SCH_DECODE)
## 11. lldec_likelihood_update.v -post_process 的使用
### 11.1 输入信号
**文件**:`design/rtl/ldec/dec_core/lldec/lldec_likelihood_update.v`第135行
```verilog
input [2-1:0] lldec_core_flip_post_process; // 2bit: [0] = post_process1, [1] = post_process2
```
### 11.2 post_process 在likelihood 更新中的使用
**文件**:`lldec_likelihood_update.v`第253-271行
```verilog
always @(*) begin
/// post_process[0]: flipped bit, likelihood == flip thr →likelihood = flip thr+1
if (lldec_core_flip_post_process[0]
&& (likelihood_new pre itr org == likelihood_flip thr)
&& !likelihood_new pre itr org sign) begin
likelihood_new post itr org = likelihood_flip thr + 1'd1;
likelihood_new post itr org sign = 1'd0;
end
/// post_process[0]: unflipped bit→likelihood = flip thr-1
else if (lldec_core_flip_post_process[0]
&& ((likelihood_new pre_itr org < likelihood_flip_thr)
likelihood_new pre itr org sign)) begin
likelihood_new post itr org = likelihood_flip_thr - 1'd1;
likelihood_new post itr org sign = 1'd0;
end
// post_process[1]: special case bump
else if (lldec_core_flip_post_process[1]
&(!likelihood_new_pre_itr_org_sign
&(likelihood_new_pre_itr_org == likelihood_flip_thr - 1'd1))
&& !flipped_old
&&(num_unsatisfy chk eq[W_LDPC_COL_WGT_MAX:0] == 1)) begin
likelihood_new post itr org = likelihood_new pre itr_org + 1'd1;
likelihood_new post itr org sign = 1'd0;
end
else begin
likelihood_new post itr org = likelihood_new pre itr org;
likelihood_new post itr org sign = likelihood_new pre itr org sign;
end
end

```
## 12. C模型当前延迟实现
## 12.1 旧函数 ldpc_dec_bf_ibex (line 3553)
```cpp
int syndrome_weight_r[5];// 5级移位寄存器
```
## 12.2 新函数 ldpc_dec_bf_ibex_rtl(line 4212)
```cpp
int syndrome_weight_r[5] = {0,0,0,0,0};
```
## 12.3 当前延迟链代码(line 4444-4449)
```cpp
syndrome_weight_r[4] = syndrome_weight_r[3];
syndrome_weight_r[3] = syndrome_weight_r[2];
syndrome_weight_r[2] = syndrome_weight_r[1];
syndrome_weight_r[1] = syndrome_weight_r[0];
syndrome_weight_r[0] = syndrome_weight;
syndrome_weight_delayed = (j < 4) ? syndrome_weight_r[4] : syndrome_weight_r[3];
```
### 12.4 C 模型列循环处理顺序(ldpc_dec_bf_ibex_rtl)
for each column j:
1. PRNG update (if post iteration)
←行4422-4431
2. shift syndrome_weight_r[], set delayed
←行4444-4449
3. syndrome_weight = f_check_node_weight(cn)←行4452
4. post_trigger based on delayed
←行4454-4470
5. per-bit: likelihood update, flip, cn update ←行4505-4669
6. rotate cn rows (delta_shift)
←行4690-4695
7. syndrome_weight = f_check_node_weight(cn) ←行4697
### 12.5 C模型 vs RTL延迟对应关系分析
C 模型的延迟链是在**每个 col 的开头**移位,然后使用 `syndrome_weight_delayed`来判断 post_process.
但在RTL中:
1. **syndw 是在当前 col 处理时同时输出的**(pipeline 组合逻辑直连)
2. **post_process 判断也是组合逻辑**(`always @*`块),直接使用 syndw 输出
3. **post_process 信号在同一时钟周期内直接影响 likelihood 更新**(组合逻辑路径)
这意味着RTL中:
- col K 处理时,使用的 syndw = pipeline 输出 = 对应 col K-4(或更早)的 SW
-这个 syndw 在**同一时钟周期**内就决定了 post_process,并影响 likelihood 更新
C模型中:
- col J 处理时,使用的 `syndrome_weight_delayed`来自移位寄存器
-移位寄存器在 col J开头移位, `r[0]`获取的是**上一个 col 结束时**的 `syndrome_weight`
**时序差异的关键**:
- RTL pipeline 捕获的是 synd_vn 寄存器输出,是 **col K 的 flip 更新前** 的 synd_vn
-C 模型的 `syndrome_weight`是在**每个 col 结束时**(flip + rotate 之后)重新计算的
-RTL 的 syndw 比当前 col 的 SW 延迟 4~5 拍
## 13. 完整信号流路径图

```text
posedge clk (col J)
        |
        |-----------------------------------------------|
        |                                               |
        v                                               |
synd_vn 更新                         Stage1 捕获       |
(flip_data XOR                       i_cn = synd_vn_old |
 delta_shift)                        (J-1 更新后的)     |
        |                                      |         |
        |                                      v         |
        |                         (posedge clk J+1)      |
        |                              Stage2 捕获       |
        |                                      |         |
        |                                      v         |
        |                         (posedge clk J+2)      |
        |                              Stage3 捕获       |
        |                                      |         |
        |                                      v         |
        |                         (posedge clk J+3)      |
        |                              Stage4 锁存       |
        |                              hamming_weight_total_r
        |                                      |
        |                                      v
        |                         lldec_core_cn_syndw = SW(J-1)
        |                                      |
        |                                      v
        |                         lldec_sch.v (always @*)
        |                         post_process = f(syndw, itr_cnt, ...)
        |                                      |
        |                                      v
        |                         lldec_likelihood_update.v (always @*)
        |                         likelihood_new = f(post_process, ...)
        |                                      |
        |                                      v
        |                         lldec_core_flip_new_ll  (组合逻辑输出)
        |                         lldec_core_flip_flipped (组合逻辑输出)
        |                         lldec_core_flip_toggle  (组合逻辑输出)
        |                                      |
        |--------------------------------------|
                                               |
                                               v
                                  (posedge clk J+4)
                                  flip_data 更新, synd_vn 更新 ...
```

## 14. 迭代边界完整周期追踪(以73列为例)

```text
                                      RTL syndw 输出
                         ┌───────────────────────────────────┐
时钟   FSM      col      synd_vn             syndw           | post_process?
---------------------------------------------------------------------------
T-73   DECODE   0        sv[0]               SW(旧iter)      | 否(iter 不够)
...
T-1    DECODE   71       sv[71]              SW(67)          | 否
T      DECODE   72       sv[72]              SW(68)          | 否
T+1    WRAP     --       sv_wrap             SW(69)          | 否(st_decode=0)
T+2    DECODE   0        sv[0]               SW(70)          | 否(旧iter 数据)
T+3    DECODE   1        sv[1]               SW(71)          | 否
T+4    DECODE   2        sv[2]               SW(72)          | 否
T+5    DECODE   3        sv[3]               SW(wrap后)      | 可能(取决于 SW 值)
T+6    DECODE   4        sv[4]               SW(col=0)       | 可能
T+7    DECODE   5        sv[5]               SW(col=1)       | 可能
...
T+N    DECODE   N-2      sv[N-2]             SW(N-6)         | 是(如果 iter>=512)
```

**核心问题**:C 模型中 col0~3的`syndrome_weight_delayed`应该对应 RTL的哪个值?
## 15. C模型延迟链的正确映射
## 15.1 C 模型列循环中syndrome_weight的更新时序
C模型每个colJ:
(1) shift r[]:r[4]=r[3],r[3]=r[2],r[2]=r[1],r[1]=r[0],r[0]=syndrome_weight
(2) syndrome_weight_delayed = (j<4) ? r[4] : r[3]
(3) syndrome_weight =f_check_node_weight(cn)←cn 是更新前的
(4)使用 syndrome_weight_delayed 做 post_trigger 判断
(5) per-bit: likelihood update, flip, cn XOR
(6) rotate cn rows (delta_shift)
(7) syndrome_weight = f_check_node_weight(cn) ←cn 是更新+rotate后的
**关键**:步骤(1)中的`syndrome_weight`是**上一列的步骤(7)**计算出来的,
即上一列处理完毕后(flip + rotate 后)的 SW.
### 15.2 RTL vs C 模型 syndw 对应关系
在RTL中:
- `i_cn = synd_vn` 直连, pipeline 捕获的是 **当前 col 的 flip 更新前** 的 synd_vn
- syndw 在 col K 输出 = SW(捕获时的 synd_vn) = SW(col K-4 更新后的 synd_vn)
在C模型中:
- `syndrome_weight` 是 **上一列结束** (flip + rotate 后)的 SW
- `r[0] = syndrome_weight`记录的是上一列结束后的 SW
- `r[N]`对应 N 列之前的 SW
所以 RTL syndw 在 col K的输出,等价于 C 模型中 col K-4 结束后的 SW.
C 模型 `r[3]`=3 列之前的 SW, `r[4]`=4 列之前的 SW.
### 15.3 迭代边界处的映射
在迭代边界:
- RTL 在 col 0~3 时, syndw 仍然输出旧迭代末尾+WRAP后的 SW
- C 模型在 col 0~3 时, `r[]`中也存储着旧迭代末尾的 SW
**但差异在于**:RTL WRAP 周期占据1个额外时钟周期,pipeline 在此期间不停.
这意味着 pipeline 中"消化"了 WRAP 周期的 synd_vn 状态.
具体地:
- 新迭代 col=0 时,RTL syndw = SW(旧iter col=70)(延迟4+WRAP偏移)
- 新迭代 col=1 时,RTL syndw = SW(旧iter col=71)
- 新迭代 col=2 时,RTL syndw = SW(旧iter col=72)
- 新迭代 col=3 时,RTL syndw = SW(WRAP后)
新迭代 col=4 时,RTL syndw = SW(col=0)
- 新迭代 col=5 时,RTL syndw = SW(col=1)

**C模型映射**:
- 新迭代 col=0 时, `syndrome_weight_delayed` = r[3] = SW(旧iter col=69)
或 r[4]= SW(旧iter col=68)
- 新迭代 col=4 时, `syndrome_weight_delayed` = r[3] = SW(旧iter col=72)
- 新迭代 col=5 时, `syndrome_weight_delayed` = r[3] = SW(新iter col=0)
**差异**:在col 4 及以后,RTL syndw 比C 模型 r[3] **多延迟1拍**(因为 WRAP 占1拍),
等价于C模型应该用 r[4]而不是 r[3].
但在 col 0~3,RTL syndw 对应旧迭代末尾的值,而C 模型 r[3]/r[4]也对应旧迭代末尾的值,
只是偏移不同.
## 16. 迭代首列的延迟链初始状态
### 16.1 首次迭代(iter=0)的 pipeline 初始状态
RTL 中 pipeline 在解码开始前(ST_SCH IDLE)被 rst_n 复位为 0.
ST_SCH SYND 阶段会向 synd_vn 载入初始 syndrome 数据,
然后进入 WRAP →INIT_LL →DECODE.
首次进入DECODE时
- synd_vn 已被初始化(SYND 阶段 + WRAP 阶段)
- pipeline 中仍有残留的0值(需要4拍才能填满)
所以首次迭代 col 0~3 的 syndw 输出可能为0 或无效值.
C 模型也应同样处理:`syndrome_weight_r[]`初始化为 0.
## 17. 附录:原始注释
### 17.1 lldec_bit_sum_512_x_13.v注释
// Description: Computes the sum of 13 512-bit syndromes. Latency: 4
### 17.2 lldec_bit_sum_512.v 注释
// Description: The sum of 512 bit data in 3 clock cycles
### 17.3 lldec_cnode.v 注释
```verilog
output [DNEINF_SYNDW_WID-1:0] lldec_core_cn_syndw; // 4 cycle delay
```
### 17.4 C 模型原始注释(ldpc_dec_bf_ibex 函数,~line 3766)
```cpp
// To match verilog resource sharing and pipelines:
// syndrome_weight module has a latency of 4,
// but because of the "wrap" cycle at the end of the matrix,
// column 0 weight has a latency of 1 fewer.
```
**注**:这个原始注释说"column 0 weight has a latency of 1 fewer",
意味着原始C模型作者认为 col0的延迟是3而不是4.
但实际上 pipeline 始终运行(i_en=1),WRAP 周期不会减少延迟,
反而会增加1拍偏移.这个注释可能是错误的.
