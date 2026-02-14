#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "ldpc_decoder.h"

using namespace std;

ldpc_decoder::ldpc_decoder()
{
    int i;
}

ldpc_decoder::~ldpc_decoder()
{

}

extern "C"
void* ldpc_decoder_new()
{
    return new ldpc_decoder();
}

extern "C"
void ldpc_decoder_inv (void* sv_inst,
                       svOpenArrayHandle sv_hardbit_data,
                       svOpenArrayHandle sv_softbit_data,
                       svOpenArrayHandle sv_softbit_data1,
                       const svOpenArrayHandle sv_decodedBits,
                       const svOpenArrayHandle sv_decoding_failure,
                       const svOpenArrayHandle sv_early_termination,
                       const svOpenArrayHandle sv_num_iter,
                       const svOpenArrayHandle sv_col_cnt,
                       const svOpenArrayHandle sv_num_errors_in_codeword,
                       const svOpenArrayHandle sv_num_errors_in_userdata,
                       const svOpenArrayHandle sv_syndrome_weight_before,
                       const svOpenArrayHandle sv_syndrome_weight_after,
                       int user_data_bytes = 0,
                       int parity_bytes = 0,
                       int post_iter = 0,
                       int max_iter = 0,
                       int nand_strobes = 1,
                       int codeword_4k_8k = 0,
                       int syndrome_weight_thr_qc = 0,
                       int syndrome_weight_thr_post = 0,
                       int early_termination_dis = 0,
                       int post_process_en = 1,
                       int ldpc_decoder_control_likelihood_0 = 0x00000000,
                       int ldpc_decoder_control_likelihood_1 = 0x00000000,
                       int ldpc_decoder_control_likelihood_2 = 0x00000000,
                       int ldpc_decoder_control_likelihood_3 = 0x00000000,
                       int ldpc_decoder_control_post = 0x00000000,
                       unsigned int ldpc_early_term_0 = 0x00000000,
                       unsigned int ldpc_early_term_1 = 0x00000000,
                       unsigned int ldpc_early_term_2 = 0x00000000,
                       unsigned int ldpc_early_term_3 = 0x00000000,
                       unsigned int ldpc_early_term_4 = 0x00000000,
                       unsigned int ldpc_early_term_5 = 0x00000000,
                       unsigned int ldpc_early_term_6 = 0x00000000
                    )
{
    int soft_bits = 0;
    int size = svHigh(sv_hardbit_data, 1);
    unsigned char* hardbit_data = new unsigned char[size + 1];
    unsigned char* softbit_data = new unsigned char[size + 1];
    unsigned char* softbit_data1 = new unsigned char[size + 1];

    int soft_bit_control_0 = 0, soft_bit_control_1 = 0;

    soft_bits = (nand_strobes == 0) ? 0 : (nand_strobes == 1) ? 1 : 2;
    printf("### LDPC DECODER C++: USER DATA BYTES = %4d, PARITY BYTES = %3d ###\n", user_data_bytes, parity_bytes);
    for (int i = 0; i <= size; i++)
    {
        hardbit_data[i] = *((unsigned char*)svGetArrElemPtr(sv_hardbit_data, i));
        if (soft_bits > 0)
            softbit_data[i] = *((unsigned char*)svGetArrElemPtr(sv_softbit_data, i));
        if (soft_bits > 1)
            softbit_data1[i] = *((unsigned char*)svGetArrElemPtr(sv_softbit_data1, i));
    }

    size = svHigh(sv_decodedBits, 1);
    unsigned char* decodedBits = new unsigned char[size + 1];
    ldpc_decoder* inst =  (ldpc_decoder*) sv_inst;
    inst->h_matrix = inst->f_h_matrix(user_data_bytes, parity_bytes);
    inst->ldpc_decoder_input.post_iteration = post_iter;
    inst->ldpc_decoder_input.iteration_limit = max_iter;
    inst->ldpc_decoder_input.nand_strobes = nand_strobes;
    inst->ldpc_decoder_input.soft_bits = soft_bits;
    inst->ldpc_decoder_input.syndrome_cal_only = syndrome_cal_only;
    inst->ldpc_decoder_parameters.post_process_en = post_process_en;
    inst->ldpc_decoder_parameters.syndrome_weight_thr_qc = syndrome_weight_thr_qc;
    inst->ldpc_decoder_parameters.syndrome_weight_thr_post = syndrome_weight_thr_post;
    inst->ldpc_decoder_parameters.likelihood_thr = (unsigned int)(ldpc_decoder_control_post >> 16)%256; // 64
    inst->ldpc_decoder_parameters.post_ratio = (unsigned int)(ldpc_decoder_control_post >> 24)%16;      // 4

    inst->ldpc_decoder_parameters.likelihood_init_coef_all[0][0] = (ldpc_decoder_control_likelihood_0 >> 0 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[0][1] = (ldpc_decoder_control_likelihood_0 >> 4 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[0][2] = (ldpc_decoder_control_likelihood_0 >> 8 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[0][3] = (ldpc_decoder_control_likelihood_0 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[1][0] = (ldpc_decoder_control_likelihood_0 >> 16) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[1][1] = (ldpc_decoder_control_likelihood_0 >> 20) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[1][2] = (ldpc_decoder_control_likelihood_0 >> 24) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[1][3] = (ldpc_decoder_control_likelihood_0 >> 28) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[2][0] = (ldpc_decoder_control_likelihood_1 >> 0 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[2][1] = (ldpc_decoder_control_likelihood_1 >> 4 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[2][2] = (ldpc_decoder_control_likelihood_1 >> 8 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[2][3] = (ldpc_decoder_control_likelihood_1 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[3][0] = (ldpc_decoder_control_likelihood_1 >> 16) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[3][1] = (ldpc_decoder_control_likelihood_1 >> 20) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[3][2] = (ldpc_decoder_control_likelihood_1 >> 24) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[3][3] = (ldpc_decoder_control_likelihood_1 >> 28) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[4][0] = (ldpc_decoder_control_likelihood_2 >> 0 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[4][1] = (ldpc_decoder_control_likelihood_2 >> 4 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[4][2] = (ldpc_decoder_control_likelihood_2 >> 8 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[4][3] = (ldpc_decoder_control_likelihood_2 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[5][0] = (ldpc_decoder_control_likelihood_2 >> 16) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[5][1] = (ldpc_decoder_control_likelihood_2 >> 20) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[5][2] = (ldpc_decoder_control_likelihood_2 >> 24) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[5][3] = (ldpc_decoder_control_likelihood_2 >> 28) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[6][0] = (ldpc_decoder_control_likelihood_3 >> 0 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[6][1] = (ldpc_decoder_control_likelihood_3 >> 4 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[6][2] = (ldpc_decoder_control_likelihood_3 >> 8 ) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[6][3] = (ldpc_decoder_control_likelihood_3 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[7][0] = (ldpc_decoder_control_likelihood_3 >> 16) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[7][1] = (ldpc_decoder_control_likelihood_3 >> 20) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[7][2] = (ldpc_decoder_control_likelihood_3 >> 24) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_coef_all[7][3] = (ldpc_decoder_control_likelihood_3 >> 28) % 16;

    inst->ldpc_decoder_parameters.likelihood_init_fraction[1] = (ldpc_early_term_0 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_init_fraction[2] = (ldpc_early_term_0 >> 28) % 16;

    inst->ldpc_decoder_parameters.likelihood_map[0] = (ldpc_early_term_1 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_map[1] = (ldpc_early_term_1 >> 28) % 16;
    inst->ldpc_decoder_parameters.likelihood_map[2] = (ldpc_early_term_2 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_map[3] = (ldpc_early_term_2 >> 28) % 16;
    inst->ldpc_decoder_parameters.likelihood_map[4] = (ldpc_early_term_3 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_map[5] = (ldpc_early_term_3 >> 28) % 16;
    inst->ldpc_decoder_parameters.likelihood_map[6] = (ldpc_early_term_4 >> 12) % 16;
    inst->ldpc_decoder_parameters.likelihood_map[7] = (ldpc_early_term_4 >> 28) % 16;

    inst->ldpc_decoder_parameters.early_terminate_dis = early_termination_dis;
    inst->ldpc_decoder_parameters.early_terminate_thr[0][0] = 600;
    inst->ldpc_decoder_parameters.early_terminate_thr[0][1] = 750;
    inst->ldpc_decoder_parameters.early_terminate_thr[0][2] = 900;
    inst->ldpc_decoder_parameters.early_terminate_thr[0][3] = 1050;
    inst->ldpc_decoder_parameters.early_terminate_thr[0][4] = 1200;
    inst->ldpc_decoder_parameters.early_terminate_thr[0][5] = 1350;
    inst->ldpc_decoder_parameters.early_terminate_thr[0][6] = 1500;
    inst->ldpc_decoder_parameters.early_terminate_thr[1][0] = 1200;
    inst->ldpc_decoder_parameters.early_terminate_thr[1][1] = 1500;
    inst->ldpc_decoder_parameters.early_terminate_thr[1][2] = 1800;
    inst->ldpc_decoder_parameters.early_terminate_thr[1][3] = 2100;
    inst->ldpc_decoder_parameters.early_terminate_thr[1][4] = 2400;
    inst->ldpc_decoder_parameters.early_terminate_thr[1][5] = 2700;
    inst->ldpc_decoder_parameters.early_terminate_thr[1][6] = 3000;

    // soft_bit_control_0
    unsigned int *c_soft_data_table = new unsigned int[8];
    for (int i = 0; i < 8; i++)
    {
        int shift_bit = i * 4;
        c_soft_data_table[i] = (soft_bit_control_0 >> shift_bit) & 0xF;
    }

    // soft_bit_control_1
    unsigned int *c_soft_data_scale = new unsigned int[8];
    for (int i = 0; i < 8; i++)
    {
        int shift_bit = i * 4;
        c_soft_data_scale[i] = (soft_bit_control_1 >> shift_bit) & 0xF;
    }

    inst->ldpc_decoder_parameters.soft_bit_table[0] = c_soft_data_table[0];
    inst->ldpc_decoder_parameters.soft_bit_table[1] = c_soft_data_table[1];
    inst->ldpc_decoder_parameters.soft_bit_table[2] = c_soft_data_table[2];
    inst->ldpc_decoder_parameters.soft_bit_table[3] = c_soft_data_table[3];
    inst->ldpc_decoder_parameters.soft_bit_table[4] = c_soft_data_table[4];
    inst->ldpc_decoder_parameters.soft_bit_table[5] = c_soft_data_table[5];
    inst->ldpc_decoder_parameters.soft_bit_table[6] = c_soft_data_table[6];
    inst->ldpc_decoder_parameters.soft_bit_table[7] = c_soft_data_table[7];

    int fractional_bytes_of_userdata = user_data_bytes - ((user_data_bytes >> 6) << 6);
    int fractional_bytes_of_parity = parity_bytes - ((parity_bytes >> 6) << 6);
    int fractional_bytes_of_total = fractional_bytes_of_userdata + fractional_bytes_of_parity;
    int unused_bytes_of_userdata = (fractional_bytes_of_userdata == 0) ? 0 : (64 - fractional_bytes_of_userdata);
    int unused_bytes_of_parity = (fractional_bytes_of_parity == 0) ? 0 : (64 - fractional_bytes_of_parity);
    int unused_bits_of_userdata = unused_bytes_of_userdata << 3;
    int unused_bits_of_parity = unused_bytes_of_parity << 3;
    int bit_index = 0;
    int bit_index_padded = 0;
    int bytes_of_padding = 0;
    if (fractional_bytes_of_userdata > 0)
        bytes_of_padding += (64 - fractional_bytes_of_userdata);
    if (fractional_bytes_of_parity > 0)
        bytes_of_padding += (64 - fractional_bytes_of_parity);
    int bits_of_padding = bytes_of_padding << 3;

    printf("### LDPC DECODER C++: USER DATA BYTES = %4d, PARITY BYTES = %3d FRACTIONAL BYTES: USER DATA: %3d PARITY: %3d TOTAL: %3d PADDING BIT: %4d UNUSED BITS USERDATA: %3d PARITY: %3d\n",
           user_data_bytes, parity_bytes,
           fractional_bytes_of_userdata,
           fractional_bytes_of_parity,
           fractional_bytes_of_total,
           bits_of_padding,
           unused_bits_of_userdata,
           unused_bits_of_parity);

    for (int j = 0; j < inst->h_matrix.cols; j++)
    {
        for (int k = 0; k < inst->h_matrix.bits; k++)
        {
            if ((j == (inst->h_matrix.cols - inst->h_matrix.rows - 1)) && (k >= (inst->h_matrix.bits - unused_bits_of_userdata)))
            {
                inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = 0;
                if (soft_bits >= 1)
                    inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = inst->ldpc_decoder_parameters.likelihood_map[3];
                if (soft_bits >= 2)
                    inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = inst->ldpc_decoder_parameters.likelihood_map[3];
            }
            else if ((j == (inst->h_matrix.cols - inst->h_matrix.rows)) && (k >= (inst->h_matrix.bits - unused_bits_of_parity)))
            {
                inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = 0;
                if (soft_bits >= 1)
                    inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = inst->ldpc_decoder_parameters.likelihood_map[3];
                if (soft_bits >= 2)
                    inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = inst->ldpc_decoder_parameters.likelihood_map[3];
            }
            else
            {
                inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = hardbit_data[bit_index];
                if (soft_bits >= 1)
                    inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = softbit_data[bit_index];
                if (soft_bits >= 2)
                    inst->ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = softbit_data1[bit_index];
                bit_index++;
            }
            bit_index_padded++;
        }
    }

    inst->ldpc_decoder_output = inst->f_ldpc_decoder(inst->ldpc_decoder_input, inst->ldpc_decoder_parameters, inst->h_matrix);
    printf("### LDPC DECODER RESULTS C++: FAILURE: %1d ERRORS: %4d %4d ITERATIONS: %4d SYNDROME WEIGHT: %5d %5d ###\n",
           inst->ldpc_decoder_output.decoding_failure,
           inst->ldpc_decoder_output.num_errors_in_codeword,
           inst->ldpc_decoder_output.num_errors_in_userdata,
           inst->ldpc_decoder_output.iterations,
           inst->ldpc_decoder_output.syndrome_weight_before,
           inst->ldpc_decoder_output.syndrome_weight_after);
    
    const int decoded_len = svHigh(sv_decodedBits, 1) + 1;
    const int stream_len  = (user_data_bytes + parity_bytes) * 8;

    if (decoded_len != stream_len) {
    fprintf(stderr, "[WARN] decodedBits length mismatch: sv=%d stream=%d\n", decoded_len, stream_len);
    }

    const int last_payload_col = inst->h_matrix.cols - inst->h_matrix.rows - 1;
    const int first_parity_col = inst->h_matrix.cols - inst->h_matrix.rows;
    const int unused_bits_userdata = inst->h_matrix.unused_bytes_of_userdata * 8;
    const int unused_bits_parity   = inst->h_matrix.unused_bytes_of_parity   * 8;

    int out_bit = 0;
    for (int j = 0; j < inst->h_matrix.cols; j++) {
    for (int k = 0; k < inst->h_matrix.bits; k++) {
        const bool is_padding =
            ((unused_bits_userdata > 0) && (j == last_payload_col) && (k >= inst->h_matrix.bits - unused_bits_userdata)) ||
            ((unused_bits_parity   > 0) && (j == first_parity_col) && (k >= inst->h_matrix.bits - unused_bits_parity));

        if (is_padding) continue;

        if (out_bit >= decoded_len) {
        fprintf(stderr,
                "[ERR] stream write exceeds decodedBits: out_bit=%d decoded_len=%d (j=%d k=%d)\n",
                out_bit, decoded_len, j, k);
        fflush(stderr);
        return;
        }

        unsigned char* p = (unsigned char*)svGetArrElemPtr(sv_decodedBits, out_bit);
        *p = inst->ldpc_decoder_output.corrected_codeword.c[j].b[k];
        out_bit++;
    }
    }

    if (out_bit != stream_len) {
    fprintf(stderr, "[ERR] stream_len mismatch after pack: out_bit=%d stream_len=%d\n", out_bit, stream_len);
    fflush(stderr);
    }

    unsigned int *c_decoding_failure = (unsigned int*)svGetArrElemPtr(sv_decoding_failure, 0);
    *c_decoding_failure = inst->ldpc_decoder_output.decoding_failure;
    unsigned int *c_early_termination = (unsigned int*)svGetArrElemPtr(sv_early_termination, 0);
    *c_early_termination = inst->ldpc_decoder_output.early_termination;
    unsigned int *c_num_errors_in_codeword = (unsigned int*)svGetArrElemPtr(sv_num_errors_in_codeword, 0);
    *c_num_errors_in_codeword = inst->ldpc_decoder_output.num_errors_in_codeword;
    unsigned int *c_num_errors_in_userdata = (unsigned int*)svGetArrElemPtr(sv_num_errors_in_userdata, 0);
    *c_num_errors_in_userdata = inst->ldpc_decoder_output.num_errors_in_userdata;
    unsigned int *c_num_iter = (unsigned int*)svGetArrElemPtr(sv_num_iter, 0);
    *c_num_iter = inst->ldpc_decoder_output.iterations;
    int *c_col_cnt = (int*)svGetArrElemPtr(sv_col_cnt, 0);
    *c_col_cnt = inst->ldpc_decoder_output.col_cnt;
    unsigned int *c_syndrome_weight_before = (unsigned int*)svGetArrElemPtr(sv_syndrome_weight_before, 0);
    *c_syndrome_weight_before = inst->ldpc_decoder_output.syndrome_weight_before;
    unsigned int *c_syndrome_weight_after = (unsigned int*)svGetArrElemPtr(sv_syndrome_weight_after, 0);
    *c_syndrome_weight_after = inst->ldpc_decoder_output.syndrome_weight_after;
}                    
