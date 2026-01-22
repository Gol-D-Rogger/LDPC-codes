#pragma once
enum fc_sim_mode {FC_SIM=0, LDPC_SIM=1};

// Simulation
enum fc_sim_mode sim_mode = LDPC_SIM;
int max_sim_num = 1000;   // 默认值：1000 个数据包
int max_err_num = 10;     // 默认值：10 个错误
int sim_step = 100;        // 默认值：每 100 个包更新一次

char *config_file;
char pchk_file[500];
char matrix_dir[500];
int h_m;
int h_n;
int h_sc;
int h_st;
int h_dense;
int h_wt;

// LDPC decoder
char dec_sel[10];
enum fc_dec_mode dec_mode;
int fdec_max_itr;
int ldec_max_itr;
int fdec_col_skip_itr;
float alpha;
int sd_num;
int rd_num;
float *vref;
float hd0_llr, hd1_llr;
int synd_wt_thrshd;

int hre_bit;
int hre_mode;
int hre_dec;
int target_fbc;  // Target total FBC for fixed-FBC mode (0=disabled)

// LDPC quantization
int finite_mode;
int finite_q_num;   
int finite_r_num;
int finite_f_num;
int finite_llr_num;
int finite_llr_f_num;

// MPCRC and MCRC
long crc32_poly_val;

// Randomizer
long rand32_poly_val;

// Simulation CH
enum ch_model ch_mode;
float ch_para;
char sd_sel[10];
enum sd_mode sd_type;

// Data format
int dsp_lba_size;
int dsp_lba_num;
int dsp_meta_size;
int dsp_lba_len;
int dsp_src_len;
int dsp_info_len;
int dsp_pad_len;
int dsp_blk_len;

// statistics

void read_arg(int, char **);
void print_usage();
void read_config_file();
