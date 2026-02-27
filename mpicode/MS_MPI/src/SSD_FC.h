// Simulation
int max_sim_num;   
int max_err_num;   
int sim_step;      

char *config_file;
int h_m;
int h_n;
int h_sc;

// LDPC decoder
int ldec_max_itr;
float alpha;
int sd_num;
int rd_num;
float *vref;
float hd0_llr, hd1_llr;

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

// Data format
int dsp_lba_size;
int dsp_lba_num;
int dsp_meta_size;
int dsp_lba_len;
int dsp_src_len;
int dsp_info_len;
int dsp_pad_len;
int dsp_blk_len;

int bytes_of_userdata;
int bytes_of_parity;

// statistics

void read_arg(int, char **);
void print_usage();
void read_config_file();
