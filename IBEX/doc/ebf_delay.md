![image](https://cdn-mineru.openxlab.org.cn/result/2026-07-08/a9fa2029-a4f0-474b-a35c-0a38557e0f7e/625e153d9c9bbf6f1b04185ca1fda7d8689aef4fe2e37222cbe04d35555926aa.jpg)


- 为了Q1T1 latency 考虑，raw data memory 本体不在core 内部，图中用虚框表示

- Flip memory 位宽要做成P_SIZEx2，并且写优先级要比读高，这样写端不会有反压，设计会简单一些；flip 每个cycle 都需要读，只有当前列有toggle 时，才进行写

- Memory wrapper 内部都支持data gating

- 打拍HW考虑打1拍和打2拍的情况，可以使用参数进行配置，默认打一拍（如图中实线DFF），如果打2拍，图中实线和虚线DFF都存在

- Syndrome weight 计算:

Syndrome weight new= Syndrome weight old +/- delat (delta 为当前迭代列非0 cirluant对应的 syndrome_weight_old - syndrome_weight_new 之和)

- DV checkpoint: 

- Task out 

Data out 

- flipped_new 

Energy 

- sk(不满足校验方程的个数)

Syndrome 

- Flip memory设计从功耗面积评估3种情况的搭配：

- 1片 2P，数据位宽为P_SIZE

- 2片 1P，数据位宽为P_SIZE/2

- 1片 1P，数据位宽为P_SIZE x 2

假双端口？

- 关于Byte align mask:

对于第0行，mask={{byte_len_rem*8{1'd1}},{(P_SIZE-byte_len_rem*8){1'd0}};

- 对于fade处，mask=~mask_0th

因为第一行是单位矩阵，mask的是第一行，因此在非P_SIZE对齐时，也只需要列重个circulant



![image-20260708235538817](/Users/roggerzwl/Library/Application Support/typora-user-images/image-20260708235538817.png)



## Pipeline:

![image-20260708235550916](/Users/roggerzwl/Library/Application Support/typora-user-images/image-20260708235550916.png)