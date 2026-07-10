# RTX 3090 row11 hard-decode runtime

Server: `root@connect.nmb2.seetacloud.com:27270`
GPU: NVIDIA GeForce RTX 3090 24GB
Code path: `/root/autodl-tmp/LDPC-codes/GPUcode/0410_run_on_server`
Config matched to `QC512/Ibex_hd_row11.cnfg`: `bytes_of_userdata=4112`, `bytes_of_parity=672`, generated `11x76` matrix.

The time below is the program's `Decoder Elapsed time` from `execute_gpu_kernel`, i.e. accumulated decode-loop time including kernel launch/sync/result copy per batch. It excludes compilation, program setup, and RNG initialization.

| SNR | RBER | requested CW | actual CW | failures | FER | avg_itr | decoder elapsed |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4.70 | 0.012148429736067569 | 10000 | 10250 | 1965 | 0.1917073171 | 16.8715 | 31.4661 s |
| 4.80 | 0.011348881149113716 | 10000 | 10250 | 139 | 0.0135609756 | 9.94878 | 24.2327 s |
| 4.85 | 0.010963277989231404 | 30000 | 30340 | 57 | 0.0018787080 | 8.44252 | 60.3859 s |

Remote result logs:

- `result_row11_3090/snr_4.70_cfg4112_total10000.txt`
- `result_row11_3090/snr_4.80_cfg4112_total10000.txt`
- `result_row11_3090/snr_4.85_cfg4112_total30000.txt`
