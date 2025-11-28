#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fc_dsp.h"
#include "rand.h"
#include "transceiver.h"
#include "utility.h"
#include "vec_op.h"
#include "SSD_FC.h"

int h_matrix_id = 0;  // 矩阵编号


int main (int argc, char **argv)
{
    FILE *fp;
    int h_k;
    int H_M;
    int H_N;
    int H_K;

    int sim_cnt;
    int sim_err = 0;
    long raw_err_tot = 0;
    long dec_err_tot = 0;
    long cw_fail_tot = 0;
    long cw_misc_tot = 0;
    long cw_mcrc_tot = 0;
    long lba_err_tot = 0;
    long meta_err_tot = 0;
    long fdec_itr_tot = 0;
    long ldec_itr_tot = 0;
    long fdec_cyc_tot = 0;
    long fdec_cyc_org = 0;
    int ldec_tot = 0;
    int fdec_false_skip_cnt = 0;
    int fdec_false_used_cnt = 0;
    int fdec_corr_skip_cnt = 0;
    long *fdec_cnvg_itr;
    long *ldec_cnvg_itr;
    long *err_num_dist;
    struct dsp_packet dsp_pckt;
    struct dsp_packet *sim_pckt = &dsp_pckt;

#ifdef _SIM_DUMP
    fp=fopen("./output/err_num_dist.txt", "w");
#endif

// initialization
    read_arg(argc, argv);
    read_config_file();
    h_k = h_n - h_m;
    H_M = h_m*h_sc;
    H_N = h_n*h_sc;
    H_K = h_k*h_sc;
    dsp_info_len = dsp_src_len + 32;
    dsp_pad_len = H_K - dsp_info_len;
    printf("[DEBUG GEN MATRIX] h_n%d, h_m%d, h_k%d, H_N%d, H_M%d, H_K%d\n", h_n, h_m, h_k, H_N, H_M, H_K);
    dsp_blk_len = dsp_info_len + H_M;
    fdec_cnvg_itr = (long*)calloc((fdec_max_itr>0) ? fdec_max_itr : -1*fdec_max_itr, sizeof(*fdec_cnvg_itr));
    ldec_cnvg_itr = (long*)calloc((ldec_max_itr>0) ? ldec_max_itr : -1*ldec_max_itr, sizeof(*ldec_cnvg_itr));
    err_num_dist = (long*)calloc(H_N/64+1, sizeof(*err_num_dist));

    printf("------------------Data Format--------------------\n");
    printf("META data       : %dB\n", dsp_meta_size/8);
    printf("LBA data        : %dB\n", dsp_lba_size/8);
    printf("LBA # per CW    : %d\n", dsp_lba_num);
    printf("DSP in data     : %dB\n", dsp_src_len/8);
    printf("MCRC size       : 4B\n");
    printf("ECC user data   : %dB\n", dsp_info_len/8);
    printf("ECC padding     : %dB\n", dsp_pad_len/8);
    printf("ECC parity      : %dB\n", H_M/8);
    printf("ECC CW size     : %dB\n", dsp_blk_len/8);
    printf("ECC CW rate     : %f\n", dsp_info_len*1.0/dsp_blk_len);
    printf("--------------------------------------------------\n");

    rand_seed(0);
#ifdef _SIM_DEBUG    
    printf("[SIM DEBUG] Configuring CH ...\n");
#endif
    sim_pckt->ch_config(dsp_info_len, dsp_blk_len, ch_mode, ch_para);
    sim_pckt->ch_llr_alloc(sd_type, sd_num, vref);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Generating LLR tables ...\n");
#endif
    sim_pckt->ch_llr_gen(hd0_llr, hd1_llr, finite_llr_num, finite_llr_f_num);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Configuring data format ...\n");
#endif
    sim_pckt->dfmt_config(dsp_meta_size, dsp_lba_size, dsp_lba_num);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Configuring LDPC ...\n");
#endif    
    sim_pckt->ldpc_config(h_m, h_n, h_sc, h_st, h_wt, pchk_file);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Configuring LDPC decoder ...\n");
#endif
    sim_pckt->ldpc_dec_config(fdec_max_itr, fdec_col_skip_itr, ldec_max_itr, alpha, finite_mode, finite_q_num, finite_r_num, finite_f_num);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Configuring randomizer ...\n");
#endif
    sim_pckt->rand_config(rand32_poly_val);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Configuring MCRC ...\n");
#endif
    sim_pckt->mcrc_config(crc32_poly_val);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Allocating simulation packet ...\n");
#endif
    sim_pckt->dsp_pckt_alloc();

    // transceiver
    printf("######################################");
    printf("[SIM] Start simulation @ ");
    print_time();
    
    for (sim_cnt=0; (((sim_cnt<max_sim_num)||(max_sim_num==0)) || sim_err<max_err_num); sim_cnt++)
    {
#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] Packet %d\n", sim_cnt);
#endif

        if ((sim_cnt!=0) && (sim_cnt%sim_step==0))
        {
            printf("[SIM] Finish %d packets @ ", sim_cnt);
            print_time();

            printf("[SIM] Statistical result of %d packets simulated: \n", sim_cnt);
            printf("[SIM] SNR       : %f\n", ch_para);
            printf("[SIM] FAIL CW   : %ld\n", cw_fail_tot);
            printf("[SIM] RAW BER   : %e\n", raw_err_tot*1.0/sim_cnt/dsp_blk_len);
            printf("[SIM] LDPC BER  : %e\n", dec_err_tot*1.0/sim_cnt/dsp_blk_len);
            printf("[SIM] LDPC FER  : %e\n", cw_fail_tot*1.0/sim_cnt);
            printf("[SIM] LDPC MIS  : %e\n", cw_misc_tot*1.0/sim_cnt);
            printf("[SIM] MCRC FER  : %e\n", cw_mcrc_tot*1.0/sim_cnt);
            printf("[SIM] DATA FER  : %e\n", (lba_err_tot+meta_err_tot)*1.0/sim_cnt);
            if ((dec_mode == FC_MIX) || (dec_mode==FC_MIX_G2))
                printf("[SIM] ECC Retry Rate: %e\n", ldec_tot*1.0/sim_cnt);
        }

        // generate random data
        if ((sim_mode==FC_SIM) && (!skip_ecc_encoder))
        {
            for (int i=0; i<dsp_lba_len; i++)
            {
                sim_pckt->wr_lba_blk[i] = rand_int(2);
            }

            for (int i=0; i<dsp_meta_size; i++)
            {
                sim_pckt->wr_meta_blk[i] = rand_int(2);
            }
        }
        else if ((sim_mode==LDPC_SIM) && (!skip_ecc_encoder))
        {
            for (int i=0; i<dsp_src_len; i++)
            {
                sim_pckt->wr_rand_blk[i] = rand_int(2);
            }
        }
        
        if ((sim_mode==FC_SIM) && (!skip_ecc_encoder))
        {
#ifdef _SIM_DEBUG
            printf("[SIM DEBUG] Insert META for packet %d\n", sim_cnt);
#endif
            sim_pckt->meta_insert();
        }

        if ((sim_mode==FC_SIM) && (!skip_ecc_encoder))
        {
#ifdef _SIM_DEBUG
            printf("[SIM DEBUG] Scramble packet %d\n", sim_cnt);
#endif
            sim_pckt->rand_encoder(1, 0, 0);
        }

#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] encoding packet %d\n", sim_cnt);
#endif
        if (!skip_ecc_encoder)
        {
            sim_pckt->ecc_encoder();
        }
        else
        {
            vec_clr(sim_pckt->tx_blk, dsp_blk_len);
        }

#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] trasmitting packet %d\n", sim_cnt);
#endif
        sim_pckt->ch_transmit();        

#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] detecting packet %d\n", sim_cnt);
#endif
        sim_pckt->ch_detector();        

#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] decoding packet %d\n", sim_cnt);
#endif
        sim_pckt->ecc_decoder(dec_mode);

        fdec_cyc_tot += sim_pckt->fdec_cyc_num;
        fdec_cyc_org += sim_pckt->fdec_cyc_org;

        if ((sim_pckt->rdec_used==0)&&(sim_pckt->init_synd_wt>=synd_wt_thrshd))
            fdec_false_skip_cnt++;
        if ((sim_pckt->rdec_used==1)&&(sim_pckt->init_synd_wt<synd_wt_thrshd))
            fdec_false_used_cnt++;
        if ((sim_pckt->rdec_used==1)&&(sim_pckt->init_synd_wt>=synd_wt_thrshd))
            fdec_corr_skip_cnt++;

        if (sim_pckt->rdec_used==1)
        {
            ldec_tot++;
            fdec_cnvg_itr[sim_pckt->fdec_max_itr-1]++;
            ldec_cnvg_itr[sim_pckt->cnvg_itr]++;
#ifdef _SIM_DEBUG
            printf("[SIM DEBUG] Packet %d TRIGGERED RETRY DECODER!\n", sim_cnt);
#endif            
        }
        else
        {
            fdec_cnvg_itr[sim_pckt->cnvg_itr]++;
        }

        if ((sim_mode==FC_SIM) && (!skip_ecc_encoder))
        {
#ifdef _SIM_DEBUG
            printf("[SIM DEBUG] Descramble packet %d\n", sim_cnt);
#endif
            sim_pckt->rand_decoder(1,0,0);
        }

        if ((sim_mode==FC_SIM) && (!skip_ecc_encoder))
        {
#ifdef _SIM_DEBUG
            printf("[SIM DEBUG] META extract for packet %d\n", sim_cnt);
#endif
            sim_pckt->meta_extract();
        }            

        if ((sim_mode==FC_SIM) && (!skip_ecc_encoder))
        {
#ifdef _SIM_DEBUG
            printf("[SIM DEBUG] LBA data compare packet %d\n", sim_cnt);
#endif
            sim_pckt->lba_err = vec_cmp(sim_pckt->wr_lba_blk, sim_pckt->rd_lba_blk, 0, 0, dsp_lba_len);

#ifdef _SIM_DEBUG
            printf("[SIM DEBUG] META data compare packet %d\n", sim_cnt);
#endif
            sim_pckt->meta_err = vec_cmp(sim_pckt->wr_meta_blk, sim_pckt->rd_meta_blk, 0, 0, dsp_meta_size);
        }
        
        if ((sim_mode == FC_SIM) && (!skip_ecc_encoder))
        {
#ifdef _SIM_DUMP
            sim_pckt->packet_dump(sim_cnt);
#endif                        
        }
        else if (sim_mode == LDPC_SIM)
        {
#ifdef _UNC_DUMP
            if ((sim_pckt->mcrc_err == 1) || (sim_pckt->cw_fail == 1))
            {
                sim_pckt->packet_dump(sim_err);
            }
#endif
        }

        raw_err_tot += sim_pckt->raw_err_num;
        dec_err_tot += sim_pckt->dec_err_num;

#ifdef _SIM_DUMP
        if (sim_pckt->raw_err_num < (H_N/64))
        {
            err_num_dist[sim_pckt->raw_err_num]++;
        }
        else
        {
            err_num_dist[H_N/64]++;
        }
#endif

        if ((sim_pckt->cw_fail == 1) || (sim_pckt->cw_miscorr == 1) ||
            (sim_pckt->mcrc_err == 1) || (sim_pckt->lba_err == 1) || (sim_pckt->meta_err == 1))
        {
            sim_err++;
            printf("[SIM] Packet %d error (%d): %d-bit errors! \n", sim_cnt, sim_err, sim_pckt->raw_err_num);
        }

        if (sim_pckt->cw_fail==1)
        {
            cw_fail_tot++;
            printf("[SIM] - DECODING FAILED(%ld)\n", cw_fail_tot);
        }
        if (sim_pckt->cw_miscorr==1)
        {
            cw_misc_tot++;
            printf("[SIM] - DECODING MISCORRECTED(%ld)\n", cw_misc_tot);
        }
        if (sim_pckt->mcrc_err==1)
        {
            cw_mcrc_tot++;
            printf("[SIM] - MCRC ERROR(%ld)\n", cw_mcrc_tot);
        }
        if (sim_pckt->lba_err==1)
        {
            lba_err_tot++;
            printf("[SIM] - LBA ERROR(%ld)\n", lba_err_tot);
        }
        if (sim_pckt->meta_err==1)
        {
            meta_err_tot++;
            printf("[SIM] - META ERROR(%ld)\n", meta_err_tot);
        }

        if ((sim_pckt->cw_fail == 1) || (sim_pckt->cw_miscorr==1) ||
            (sim_pckt->mcrc_err==1) || (sim_pckt->lba_err==1) || (sim_pckt->meta_err==1))
        {
            printf("[SIM] Statistical result of %d packets simulated: \n", (sim_cnt+1));
            printf("[SIM] SNR       : %f\n", ch_para);
            printf("[SIM] FAIL CW   : %ld\n", cw_fail_tot);
            printf("[SIM] RAW BER   : %e\n", raw_err_tot*1.0/(sim_cnt+1)/dsp_blk_len);
            printf("[SIM] LDPC BER  : %e\n", dec_err_tot*1.0/(sim_cnt+1)/dsp_blk_len);
            printf("[SIM] LDPC FER  : %e\n", cw_fail_tot*1.0/(sim_cnt+1));
            printf("[SIM] LDPC MIS  : %e\n", cw_misc_tot*1.0/(sim_cnt+1));
            printf("[SIM] MCRC FER  : %e\n", cw_mcrc_tot*1.0/(sim_cnt+1));
            printf("[SIM] DATA FER  : %e\n", (lba_err_tot+meta_err_tot)*1.0/(sim_cnt+1));
            if ((dec_mode == FC_MIX) || (dec_mode==FC_MIX_G2))
                printf("[SIM] ECC Retry Rate: %e\n", ldec_tot*1.0/(sim_cnt+1));
        } 
    }

    printf("[SIM] Finish simulation @ ");
    print_time();

    printf("--------------------------------------------------------\n");
    printf("[STATISTICS] Total packets simulated: %d\n", sim_cnt);
    printf("[STATISTICS] RAW BER   : %e\n", raw_err_tot*1.0/sim_cnt/dsp_blk_len);
    printf("[STATISTICS] LDPC BER  : %e\n", dec_err_tot*1.0/sim_cnt/dsp_blk_len);
    printf("[STATISTICS] LDPC FER  : %e\n", cw_fail_tot*1.0/sim_cnt);
    printf("[STATISTICS] LDPC MIS  : %e\n", cw_misc_tot*1.0/sim_cnt);
    printf("[STATISTICS] MCRC FER  : %e\n", cw_mcrc_tot*1.0/sim_cnt);
    printf("[STATISTICS] DATA FER  : %e\n", (lba_err_tot+meta_err_tot)*1.0/sim_cnt);

    if ((dec_mode == FC_MIX) || (dec_mode==FC_MIX_G2))
    {
        printf("[STATISTICS] ECC Retry Rate: %e\n", ldec_tot*1.0/sim_cnt);
        printf("[STATISTICS] Fast decoder falsely skipped   #: %d\n", fdec_false_skip_cnt);
        printf("[STATISTICS] Fast decoder correctly skipped #: %d (%d)\n", fdec_corr_skip_cnt, fdec_false_used_cnt);
        printf("[STATISTICS] LDPC initial syndrome weight    : [%d, %d]\n", sim_pckt->init_synd_wt_min, sim_pckt->init_synd_wt_max);
    }
    printf("---------------------------------------------------------\n");
    
    for (int i=0; i<sim_pckt->fdec_max_itr; i++)
    {
        fdec_itr_tot += (i+1)*fdec_cnvg_itr[i];
    }
    for (int i=0; i<sim_pckt->ldec_max_itr; i++)
    {
        ldec_itr_tot += (i+1)*ldec_cnvg_itr[i];
    }

    if ((dec_mode == FC_FDEC) || (dec_mode==FC_FDEC_G2) || (dec_mode==FC_MIX) || (dec_mode==FC_MIX_G2))
    {
        printf("[STATISTICS] Fast decoder average iteration: %f\n", fdec_itr_tot*1.0/sim_cnt);
        printf("[STATISTICS] Fast decoder saves %f cycles from column skip\n", (100-100.0*fdec_cyc_tot/fdec_cyc_org));

        printf("[STATISTICS] Iteration distribution as below: \n");
        for (int i=0; i<sim_pckt->fdec_max_itr; i=i+8)
        {
            printf("Itr: ");
            for (int j=0; j<8; j++)
            {
                if ((i+j) < sim_pckt->fdec_max_itr)
                    printf("%8d |", i+j+1);
            }
            printf("\n");
            printf("Num: ");
            for (int j=0; j<8; j++)
            {
                if ((i+j) < sim_pckt->fdec_max_itr)
                    printf("%8ld |", fdec_cnvg_itr[i+j]);
            }
            printf("\n");
            printf("------------------------------------------------------------\n");
        }
    }

    if ((dec_mode==FC_RDEC) || ((dec_mode==FC_MIX) && ldec_tot>0) || ((dec_mode==FC_MIX_G2) && ldec_tot>0))
    {
        printf("[STATISTICS] Retry decoder average iteration: %f\n", ldec_itr_tot*1.0/ldec_tot);

        printf("[STATISTICS] Retry decoder iteration distribution as below: \n");
        for (int i=0; i<sim_pckt->ldec_max_itr; i=i+8)
        {
            printf("Itr: ");
            for (int j=0; j<8; j++)
            {
                if ((i+j) < sim_pckt->ldec_max_itr)
                    printf("%8d |", i+j+1);
            }
            printf("\n");
            printf("Num: ");
            for (int j=0; j<8; j++)
            {
                if ((i+j) < sim_pckt->ldec_max_itr)
                    printf("%8ld |", ldec_cnvg_itr[i+j]);
            }
            printf("\n");
            printf("------------------------------------------------------------\n");
        }
    }

#ifdef _SIM_DUMP
    for (int i=0; i<=H_N/64; i++)
    {
        if (err_num_dist[i]!=0)
            fprintf(fp, "ERR # = %3d: %6e\n", i, err_num_dist[i]*1.0/sim_cnt);
    }
#endif

    // clean up
    sim_pckt->ldpc_clean();
    sim_pckt->mcrc_clean();
    sim_pckt->rand_clean();
    sim_pckt->dsp_pckt_clean();
    sim_pckt->ch_llr_clean();
    free(fdec_cnvg_itr);
    free(ldec_cnvg_itr);
    free(err_num_dist);
    free(vref);

#ifdef _SIM_DUMP
    fclose(fp);
#endif    

    return 0;
}



void read_arg(int argc, char **argv)
{
    char *sim_sel;
    char *ch_sel;
    char junk;
    int matrix_id = 0;  // 添加矩阵编号变量

    if (!(sim_sel=argv[1]) || !(config_file=argv[2]) || !(ch_sel=argv[3]))
        print_usage();

    if (strcmp(sim_sel, "FC")==0)
        sim_mode = FC_SIM;
    else if (strcmp(sim_sel, "LDPC")==0)
        sim_mode = LDPC_SIM;
    else
        print_usage();

    if (strcmp(ch_sel, "CLEAN") == 0)
        ch_mode = CLEAN;
    else if (strcmp(ch_sel, "AWGN") == 0)
        ch_mode = AWGN;
    else if (strcmp(ch_sel, "BSC") == 0)
        ch_mode = BSC;
    else if (strcmp(ch_sel, "ERR_INJ") == 0)
        ch_mode = ERR_INJ;
    else if (strcmp(ch_sel, "MAX_ERR") == 0)
        ch_mode = MAX_ERR;
    else if (strcmp(ch_sel, "ALL_ZERO") == 0)
        ch_mode = ALL_ZERO;
    else
        print_usage();

    if (ch_mode != CLEAN)
    {
        if ((!argv[4] || sscanf(argv[4], "%f%c", &ch_para, &junk) != 1))
            print_usage();
        
        // 添加矩阵编号参数检查
        if (argv[5] && sscanf(argv[5], "%d%c", &matrix_id, &junk) == 1)
        {
            // 将矩阵编号设置为全局变量或传递给配置函数
            h_matrix_id = matrix_id;
        }
        else
        {
            h_matrix_id = 0;  // 默认编号为0
        }

        // 矩阵目录（可选第6参），否则环境变量 LDPC_MATRIX_DIR，否则默认 ./matrix
        const char *md_env = getenv("LDPC_MATRIX_DIR");
        if (argv[6])
            strncpy(matrix_dir, argv[6], sizeof(matrix_dir)-1), matrix_dir[sizeof(matrix_dir)-1]='\0';
        else if (md_env)
            strncpy(matrix_dir, md_env, sizeof(matrix_dir)-1), matrix_dir[sizeof(matrix_dir)-1]='\0';
        else
            strcpy(matrix_dir, "./matrix");
    }
    else
    {
        // CLEAN模式也可以接收矩阵编号
        if (argv[4] && sscanf(argv[4], "%d%c", &matrix_id, &junk) == 1)
        {
            h_matrix_id = matrix_id;
        }
        else
        {
            h_matrix_id = 0;
        }

        // 矩阵目录：第5参（CLEAN 少一个 ch_para），或环境变量，或默认
        const char *md_env = getenv("LDPC_MATRIX_DIR");
        if (argv[5])
            strncpy(matrix_dir, argv[5], sizeof(matrix_dir)-1), matrix_dir[sizeof(matrix_dir)-1]='\0';
        else if (md_env)
            strncpy(matrix_dir, md_env, sizeof(matrix_dir)-1), matrix_dir[sizeof(matrix_dir)-1]='\0';
        else
            strcpy(matrix_dir, "./matrix");
    }
}

void print_usage()
{
    printf("Usage: ./ssd_fc.out sim_mode config_file ch_model ch_para(OPTIONAL)\n");
    printf("  SIM_MODE: FC or LDPC\n");
    printf("  CH_MODE: CLEAN, AWGN, BSC, ERR_INJ, MAX_ERR, ALL_ZERO\n");
    printf("  CH_PARA: not needed for CLEAN channel\n");
    printf("           SNR for AWGN/ALL_ZERO\n");
    printf("           RBER for BSC channel\n");
    printf("           Number of injected bit errors for ERR_INJ channel\n");
    printf("           Max number of injected bit errors for MAX_ERR channel\n");
    printf("           ALL_ZERO forces zero codeword with given SNR (AWGN)\n");
    exit(1);
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

    // 0. Maximum LDPC BF decoding iteration
    fscanf(fp, "%d", &fdec_max_itr);
    fgets(str_tmp, 800, fp);

    // 9.1 FDEC iteration number to turn on column skip featur
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

    // 25. Skip ECC encoder flag (optional, default 0)
    if (fscanf(fp, "%d", &skip_ecc_encoder) == 1)
    {
        fgets(str_tmp, 800, fp);
    }
    else
    {
        skip_ecc_encoder = 0;
    }

    fclose(fp);

    // print out simulation configuration
    printf("---------------Simulation Configuration----------------\n");
    printf("LDPC matrix %dx%dex%dwt%d-T(%d)\n", h_m, h_n, h_sc, h_wt, h_st);
    printf("H matrix file: %s\n", pchk_file);
    if (ch_mode==CLEAN)
        printf("Clean channel\n");
    else if (ch_mode==AWGN)
        printf("AWGN channel with SNR=%f\n", ch_para);
    else if (ch_mode==BSC)
        printf("BSC channel with RBER=%f\n", ch_para);
    else if (ch_mode==ERR_INJ)
        printf("Error injection channel with %d random bit errors\n", (int)ch_para);
    else if (ch_mode==MAX_ERR)
        printf("Max Error injection channel with up to %d random bit errors\n", (int)ch_para);
    else if (ch_mode==ALL_ZERO)
        printf("ALL_ZERO channel forcing zero codeword with SNR=%f\n", ch_para);

    if (skip_ecc_encoder)
        printf("ECC encoder: SKIPPED (tx_blk uses custom/all-zero pattern)\n");
    else
        printf("ECC encoder: enabled\n");
    
    if (dec_mode==FC_SKIP)
        printf("Skip ECC decoder\n");
    else if (dec_mode==FC_FDEC)
        printf("Fast decoder only (HD only)\n");
    else if (dec_mode==FC_RDEC)
        printf("Retry decoder only (%d-bit soft data)\n", sd_num);
    else if (dec_mode==FC_MIX)
        printf("LDPC decoder: Fast + Retry decoder (HD only)\n");
    else if (dec_mode==FC_FDEC_G2)
        printf("Fast decoder (GEN2) only (HD only)\n");
    else if (dec_mode==FC_MIX_G2)
        printf("Fast (GEN2) + retry decoder (HD only)\n");

    printf("Fast decoder: max_iter = %d\n", fdec_max_itr);
    if (fdec_col_skip_itr == 0)
        printf("Column skip feature OFF\n");
    else
        printf("Column skip @ iterations %d\n", fdec_col_skip_itr);
    printf("Retry decoder: max_iter = %d\n", ldec_max_itr);
    printf("alpha = %f\n", alpha);

    if (finite_mode==0)
        printf("Floating simulation \n");
    else
    {
        printf("Q/APP (%d-bit tot/%d-bit frac)\n", finite_q_num, finite_f_num);
        printf("R   (%d-bit tot/%d-bit frac)\n", finite_r_num, finite_f_num);
    }

    if (max_sim_num!=0)
        printf("Simulation stops @ %d packets or %d errors\n", max_sim_num, max_err_num);
    else
        printf("Simulation stops @ %d errors\n", max_err_num);
    printf("Update simulation status per %d packets\n", sim_step);
    printf("--------------------------------------------------------\n");
}
