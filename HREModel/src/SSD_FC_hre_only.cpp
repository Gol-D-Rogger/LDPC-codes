#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rand.h"
#include "transceiver_hre_only.h"
#include "SSD_FC_hre_only.h"

int h_matrix_id = 0;  // 保留全局编号但不生成矩阵

int main(int argc, char **argv)
{
    int sim_cnt;
    long raw_err_tot = 0;
    long raw_err_awgn_tot = 0;
    long hre_like_tot = 0;

    // 读取参数与配置
    read_arg(argc, argv);
    read_config_file();

    int h_k = h_n - h_m;
    int H_M = h_m * h_sc;
    int H_N = h_n * h_sc;
    int H_K = h_k * h_sc;
    dsp_info_len = dsp_src_len + 32;
    dsp_pad_len = H_K - dsp_info_len;
    dsp_blk_len = dsp_info_len + H_M;

    printf("[HRE ONLY] CW len %d, info %d, parity %d, rate %f\n",
           dsp_blk_len, dsp_info_len, H_M, dsp_info_len*1.0/dsp_blk_len);

    // 配置信道
    ch_packet sim_ch;
    rand_seed(0);
    sim_ch.ch_config(dsp_info_len, dsp_blk_len, ch_mode, ch_para, hre_bit, hre_mode, hre_dec);
    sim_ch.ch_pckt_alloc();
    sim_ch.ch_llr_alloc(sd_type, sd_num, vref);
    sim_ch.ch_llr_gen(hd0_llr, hd1_llr, finite_llr_num, finite_llr_f_num);

    printf("[HRE ONLY] Start simulation for %d packets\n", max_sim_num);
    for (sim_cnt = 0; (sim_cnt < max_sim_num) || (max_sim_num == 0); sim_cnt++)
    {
        // 随机生成 TX
        for (int i = 0; i < dsp_blk_len; i++)
        {
            sim_ch.tx_blk[i] = rand_int(2);
        }

        sim_ch.ch_transmit();
        sim_ch.ch_detector();

        raw_err_tot += sim_ch.raw_err_num;
        raw_err_awgn_tot += sim_ch.raw_err_awgn;
        hre_like_tot += sim_ch.hre_like_cnt;

        if ((sim_cnt != 0) && (sim_cnt % sim_step == 0))
        {
            printf("[HRE ONLY] %d pkts | RAW BER %e | FBC(bit) %f | FBC(no HRE) %f | HRE-like %f\n",
                   sim_cnt,
                   raw_err_tot*1.0/sim_cnt/dsp_blk_len,
                   raw_err_tot*1.0/sim_cnt,
                   raw_err_awgn_tot*1.0/sim_cnt,
                   hre_like_tot*1.0/sim_cnt);
        }
    }

    printf("--------------------------------------------------------\n");
    printf("[HRE ONLY] Total packets simulated: %d\n", sim_cnt);
    printf("[HRE ONLY] RAW BER       : %e\n", raw_err_tot*1.0/sim_cnt/dsp_blk_len);
    printf("[HRE ONLY] FBC(bit)      : %f\n", raw_err_tot*1.0/sim_cnt);
    printf("[HRE ONLY] FBC(no HRE)   : %f\n", raw_err_awgn_tot*1.0/sim_cnt);
    printf("[HRE ONLY] HRE-like(cnt) : %f\n", hre_like_tot*1.0/sim_cnt);
    printf("--------------------------------------------------------\n");

    sim_ch.ch_llr_clean();
    sim_ch.ch_pckt_clean();
    free(vref);

    return 0;
}

// 以下保留原有配置读取逻辑，便于与主仿真一致
void print_usage()
{
    printf("Usage: ./ssd_fc_hre_only.out sim_mode config_file ch_model ch_para(OPTIONAL)\n");
    printf("  SIM_MODE: FC or LDPC\n");
    printf("  CH_MODE: CLEAN, AWGN, BSC, ERR_INJ, MAX_ERR\n");
    printf("  CH_PARA: not needed for CLEAN channel\n");
    printf("           SNR for AWGN\n");
    printf("           RBER for BSC channel\n");
    printf("           Number of injected bit errors for ERR_INJ channel\n");
    printf("           Max number of injected bit errors for MAX_ERR channel\n");
    exit(1);
}

void read_arg(int argc, char **argv)
{
    if (argc < 3)
        print_usage();

    if (strcmp(argv[1], "FC") == 0)
    {
        sim_mode = FC_SIM;
    }
    else if (strcmp(argv[1], "LDPC") == 0)
    {
        sim_mode = LDPC_SIM;
    }
    else
    {
        print_usage();
    }

    config_file = argv[2];

    if (strcmp(argv[3], "CLEAN") == 0)
    {
        ch_mode = CLEAN;
    }
    else if(strcmp(argv[3], "AWGN")==0)
    {
        ch_mode = AWGN;
        if (argc > 4)
        {
            ch_para = atof(argv[4]);
        }
        else
            ch_para = 0;
    }
    else if(strcmp(argv[3], "TRUNC_AWGN")==0 || strcmp(argv[3], "TAWGN")==0)
    {
        ch_mode = TRUNC_AWGN;
        if (argc > 4)
        {
            ch_para = atof(argv[4]);
        }
        else
            ch_para = 0;
    }
    else if(strcmp(argv[3], "BSC")==0)
    {
        ch_mode = BSC;
        if (argc > 4)
        {
            ch_para = atof(argv[4]);
        }
        else
            ch_para = 3e-3;
    }
    else if(strcmp(argv[3], "ERR_INJ")==0)
    {
        ch_mode = ERR_INJ;
        if (argc > 4)
        {
            ch_para = atof(argv[4]);
        }
        else
            ch_para = 10;
    }
    else if(strcmp(argv[3], "MAX_ERR")==0)
    {
        ch_mode = MAX_ERR;
        if (argc > 4)
        {
            ch_para = atof(argv[4]);
        }
        else
            ch_para = 10;
    }
    else
        print_usage();

    if (argc >= 6)
    {
        // 支持可配置矩阵目录
        const char *md_env = getenv("LDPC_MATRIX_DIR");
        if (argv[5][0] != '-')
            strncpy(matrix_dir, argv[5], sizeof(matrix_dir)-1), matrix_dir[sizeof(matrix_dir)-1]='\0';
        else if (md_env)
            strncpy(matrix_dir, md_env, sizeof(matrix_dir)-1), matrix_dir[sizeof(matrix_dir)-1]='\0';
        else
            strcpy(matrix_dir, "./matrix");
    }
}

void read_config_file()
{
    FILE *fp;
    char str_tmp[800];

    fp = fopen(config_file, "r");

    // 0 data format
    fscanf(fp, "%d", &dsp_meta_size);
    fgets(str_tmp, 800, fp);

    fscanf(fp, "%d", &dsp_lba_size);
    fgets(str_tmp, 800, fp);

    fscanf(fp, "%d", &dsp_lba_num);
    fgets(str_tmp, 800, fp);

    dsp_meta_size = dsp_meta_size*8; // byte --> bit
    dsp_lba_size = dsp_lba_size*8; // byte --> bit
    dsp_lba_len = dsp_lba_size*dsp_lba_num;
    dsp_src_len = dsp_meta_size + dsp_lba_len;

    // 1. h_m (row number of base matrix)
    fscanf(fp, "%d", &h_m);
    fgets(str_tmp, 800, fp);

    // 2. h_n (column number of base matrix)
    fscanf(fp, "%d", &h_n);
    fgets(str_tmp, 800, fp);

    // 3. h_sc (expansion factor of base matrix)
    fscanf(fp, "%d", &h_sc);
    fgets(str_tmp, 800, fp);

    // 4. h_dense (rows of E matrix in base matrix)
    fscanf(fp, "%d", &h_dense);
    fgets(str_tmp, 800, fp);
    h_st = h_m - h_dense;

    // 5. h_wt (column weight)
    fscanf(fp, "%d", &h_wt);
    fgets(str_tmp, 800, fp);

    // generate pchk file（支持可配置矩阵目录）
    sprintf(pchk_file, "%s/LDPC_%dx%dex%d_w%d_dense%d_QC_H.txt", matrix_dir, h_m, h_n, h_sc, h_wt, h_dense);

    // 6. max_sim_num
    fscanf(fp, "%d", &max_sim_num);
    fgets(str_tmp, 800, fp);

    sim_step = (max_sim_num>=1000) ? max_sim_num/10 : 10000;

    // 7. max_err_num
    fscanf(fp, "%d", &max_err_num);
    fgets(str_tmp, 800, fp);

    // 8. LDPC decoder selection
    fscanf(fp, "%s", dec_sel);
    fgets(str_tmp, 800, fp);

    if (strcmp(dec_sel, "FC_SKIP")==0)
        dec_mode = FC_SKIP;
    else if (strcmp(dec_sel, "FC_FDEC")==0)
        dec_mode = FC_FDEC;
    else if (strcmp(dec_sel, "FC_FDEC2")==0)
        dec_mode = FC_FDEC_G2;
    else if (strcmp(dec_sel, "FC_RDEC")==0)
        dec_mode = FC_RDEC;
    else if (strcmp(dec_sel, "FC_MIX")==0)
        dec_mode = FC_MIX;
    else if (strcmp(dec_sel, "FC_MIX2")==0)
        dec_mode = FC_MIX_G2;

    // HRE Bit Count
    fscanf(fp, "%d", &hre_bit);
    fgets(str_tmp, 800, fp);

    // HRE Bit Set
    fscanf(fp, "%d", &hre_mode);
    fgets(str_tmp, 800, fp);

    // HRE Bit LLR Set
    fscanf(fp, "%d", &hre_dec);
    fgets(str_tmp, 800, fp);

    // 9. Maximum LDPC BF decoding iteration
    fscanf(fp, "%d", &fdec_max_itr);
    fgets(str_tmp, 800, fp);

    // 9.1 FDEC iteration number to turn on column skip feature
    fscanf(fp, "%d", &fdec_col_skip_itr);
    fgets(str_tmp, 800, fp);
    
    // 10. Maximum LDPC Layer decoding iteration
    fscanf(fp, "%d", &ldec_max_itr);
    fgets(str_tmp, 800, fp);

    // 11. alpha
    fscanf(fp, "%f", &alpha);
    fgets(str_tmp, 800, fp);

    // 12. finite_mode: 0-float point 1-quantization
    fscanf(fp, "%d", &finite_mode);
    fgets(str_tmp, 800, fp);

    // 13. Q/APP message bit width
    fscanf(fp, "%d", &finite_q_num);
    fgets(str_tmp, 800, fp);

    // 14. R message bit width
    fscanf(fp, "%d", &finite_r_num);
    fgets(str_tmp, 800, fp);

    // 15. Fraction bit width
    fscanf(fp, "%d", &finite_f_num);
    fgets(str_tmp, 800, fp);

    // 16. LLR bit width
    fscanf(fp, "%d", &finite_llr_num);
    fgets(str_tmp, 800, fp);

    // 17. LLR fraction bit width
    fscanf(fp, "%d", &finite_llr_f_num);
    fgets(str_tmp, 800, fp);

    // 18. NAND LLR mode
    fscanf(fp, "%s", sd_sel);
    fgets(str_tmp, 800, fp);
    if (strcmp(sd_sel, "MANUAL")==0)
        sd_type = MANUAL;
    else if (strcmp(sd_sel, "VENDOR0")==0)
        sd_type = VENDOR0;
    else if (strcmp(sd_sel, "VENDOR1")==0)
        sd_type = VENDOR1;
    else if (strcmp(sd_sel, "DIRECT")==0)
        sd_type = DIRECT;

    // 19. NAND read number
    fscanf(fp, "%d", &sd_num);
    fgets(str_tmp, 800, fp);
    if (sd_type == MANUAL)
        rd_num = sd_num;
    else if (sd_type == VENDOR0)
        rd_num = 2*sd_num - 1;
    else if (sd_type == VENDOR1)
        rd_num = 2*sd_num - 1;
    else if (sd_type == DIRECT)
        rd_num = sd_num;

    // 20. read reference voltage
    vref = (float*)calloc(rd_num, sizeof(*vref));
    for (int i=0; i<rd_num; i++)
    {
        if (sd_type != DIRECT)
            fscanf(fp, "%f", vref+i);
        else
            vref[i] = 0;
    }
    fgets(str_tmp, 800, fp);

    if (dec_mode!=FC_RDEC)
    {
        sd_num = 1;
        rd_num = 1;
        vref[0] = 0;
    }

    // 21. HD LLRs for retry decoder
    fscanf(fp, "%f", &hd0_llr);
    fscanf(fp, "%f", &hd1_llr);
    fgets(str_tmp, 800, fp);

    // 22. Initial syndrome weight threshold to skip FDEC
    fscanf(fp, "%d", &synd_wt_thrshd);
    fgets(str_tmp, 800, fp);

    // 23. CRC-32 polynomial value for MCRC and MPCRC
    fscanf(fp, "%lx", &crc32_poly_val);
    fgets(str_tmp, 800, fp);

    // 24. 32-bit LFSR polynomial value for randomizer
    fscanf(fp, "%lx", &rand32_poly_val);
    fgets(str_tmp, 800, fp);

    fclose(fp);

    // print out simulation configuration (简化输出)
    printf("---------------HRE ONLY Configuration----------------\n");
    printf("LDPC matrix %dx%dex%dwt%d-T(%d)\n", h_m, h_n, h_sc, h_wt, h_st);
    if (ch_mode==CLEAN)
        printf("Clean channel\n");
    else if (ch_mode==AWGN)
        printf("AWGN channel with SNR=%f\n", ch_para);
    else if (ch_mode==TRUNC_AWGN)
        printf("Truncated AWGN channel with SNR=%f (fixed FBC)\n", ch_para);
    else if (ch_mode==BSC)
        printf("BSC channel with RBER=%f\n", ch_para);
    else if (ch_mode==ERR_INJ)
        printf("Error injection channel with %d random bit errors\n", (int)ch_para);
    else if (ch_mode==MAX_ERR)
        printf("Max Error injection channel with up to %d random bit errors\n", (int)ch_para);
    printf("HRE bit num %d, mode %d, LLR proc %d\n", hre_bit, hre_mode, hre_dec);
    printf("-----------------------------------------------------\n");
}
