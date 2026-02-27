#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <iostream>
#include "fc_dsp.h"
#include "rand.h"
#include "transceiver.h"
#include "utility.h"
#include "vec_op.h"
#include "SSD_FC.h"
#include <mpi.h>


int main (int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0)
        printf("There are %d MPI processes\n", size);

    FILE *fp;
    int h_k;
    int H_M;
    int H_N;
    int H_K;
    bool use_zero_cw = true;
    long sim_cnt;
    long sim_err = 0;
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
    long ldec_tot = 0;
    int fdec_false_skip_cnt = 0;
    int fdec_false_used_cnt = 0;
    int fdec_corr_skip_cnt = 0;
    long *fdec_cnvg_itr;
    long *ldec_cnvg_itr;
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

    dsp_info_len =bytes_of_userdata*8;
    dsp_pad_len = H_K - dsp_info_len;
    if (rank == 0)
        printf("[DEBUG GEN MATRIX] h_n%d, h_m%d, h_k%d, H_N%d, H_M%d, H_K%d\n", h_n, h_m, h_k, H_N, H_M, H_K);
    dsp_blk_len = dsp_info_len + bytes_of_parity*8;
    
    ldec_cnvg_itr = (long*)calloc((ldec_max_itr>0) ? ldec_max_itr : -1*ldec_max_itr, sizeof(*ldec_cnvg_itr));

    if (rank == 0)
    {
        printf("------------------Data Format--------------------\n");
        printf("META data       : %dB\n", dsp_meta_size/8);
        printf("LBA data        : %dB\n", dsp_lba_size/8);
        printf("LBA # per CW    : %d\n", dsp_lba_num);
        printf("DSP in data     : %dB\n", dsp_src_len/8);
        printf("MCRC size       : 4B\n");
        printf("ECC user data   : %dB\n", dsp_info_len/8);
        printf("ECC padding     : %dB\n", dsp_pad_len/8);
        printf("ECC parity      : %dB\n", bytes_of_parity);
        printf("ECC CW size     : %dB\n", dsp_blk_len/8);
        printf("ECC CW rate     : %f\n", dsp_info_len*1.0/dsp_blk_len);
        printf("--------------------------------------------------\n");
    }

    int seed = static_cast<int>(time(nullptr)) + rank*19;
    rand_seed(seed);
    printf("rank %d seed %d\n", rank, seed);

#ifdef _SIM_DEBUG    
    printf("[SIM DEBUG] Configuring CH ...\n");
#endif
    sim_pckt->ch_config(dsp_info_len, dsp_blk_len, ch_mode, ch_para);
    sim_pckt->ch_llr_alloc(sd_num, vref);

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
    sim_pckt->ldpc_config(h_m, h_n, h_sc);

#ifdef _SIM_DEBUG
    printf("[SIM DEBUG] Configuring LDPC decoder ...\n");
#endif
    sim_pckt->ldpc_dec_config(ldec_max_itr, alpha, finite_mode, finite_q_num, finite_r_num, finite_f_num);

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
    
    for (sim_cnt = 0; sim_cnt < max_sim_num; sim_cnt++)
    {
#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] Packet %ld\n", sim_cnt);
#endif

        if ((sim_cnt!=0) && (sim_cnt%sim_step==0))
        {
            printf("[SIM] Finish %ld packets @ ", sim_cnt);
            print_time();

            printf("[SIM] Statistical result of %ld packets simulated: \n", sim_cnt);
            printf("[SIM] SNR       : %f\n", ch_para);
            printf("[SIM] FAIL CW   : %ld\n", cw_fail_tot);
            printf("[SIM] RAW BER   : %e\n", raw_err_tot*1.0/sim_cnt/dsp_blk_len);
            printf("[SIM] LDPC BER  : %e\n", dec_err_tot*1.0/sim_cnt/dsp_blk_len);
            printf("[SIM] LDPC FER  : %e\n", cw_fail_tot*1.0/sim_cnt);
            printf("[SIM] LDPC MIS  : %e\n", cw_misc_tot*1.0/sim_cnt);
            printf("[SIM] MCRC FER  : %e\n", cw_mcrc_tot*1.0/sim_cnt);
            printf("[SIM] DATA FER  : %e\n", (lba_err_tot+meta_err_tot)*1.0/sim_cnt);
        }

        if (!use_zero_cw)
        {
            for (int i=0; i<dsp_meta_size; i++)
            {
                sim_pckt->wr_rand_blk[i] = rand_int(2);
            }
            for (int i=0; i<dsp_lba_size*dsp_lba_num; i++)
            {
                sim_pckt->wr_rand_blk[i+dsp_meta_size] = rand_int(2);
            }
            sim_pckt->ecc_encoder();
        }

#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] transmitting packet %d\n", sim_cnt);
#endif
        sim_pckt->ch_transmit();        

#ifdef _SIM_DEBUG
        printf("[SIM DEBUG] detecting packet %d\n", sim_cnt);
#endif
        sim_pckt->ch_detector();        

        sim_pckt->ecc_decoder();

        if (sim_pckt->rdec_used==1)
        {
            ldec_tot++;
            ldec_cnvg_itr[sim_pckt->cnvg_itr]++;
        }

        raw_err_tot += sim_pckt->raw_err_num;
        dec_err_tot += sim_pckt->dec_err_num;

        if ((sim_pckt->cw_fail == 1) || (sim_pckt->cw_miscorr == 1) ||
            (sim_pckt->mcrc_err == 1) || (sim_pckt->lba_err == 1) || (sim_pckt->meta_err == 1))
        {
            sim_err++;
            printf("[SIM] Packet %ld error (%ld): %d-bit errors! at rank (%d)\n", sim_cnt, sim_err, sim_pckt->raw_err_num, rank);
        }
        
        if (sim_pckt->cw_fail==1)
        {
            cw_fail_tot++;
            printf("[SIM] - DECODING FAILED(%ld) at rank (%d)\n", cw_fail_tot, rank);
        }
        if (sim_pckt->cw_miscorr==1)
        {
            cw_misc_tot++;
            printf("[SIM] - DECODING MISCORRECTION(%ld)\n", cw_misc_tot);
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

        // if ((sim_pckt->cw_fail == 1) || (sim_pckt->cw_miscorr==1) ||
        //     (sim_pckt->mcrc_err==1) || (sim_pckt->lba_err==1) || (sim_pckt->meta_err==1))
        // {
        //     printf("[SIM] Statistical result of %d packets simulated: \n", sim_cnt+1);
        //     printf("[SIM] SNR       : %f\n", ch_para);
        //     printf("[SIM] FAIL CW   : %ld\n", cw_fail_tot);
        //     printf("[SIM] RAW BER   : %e\n", raw_err_tot*1.0/(sim_cnt+1)/dsp_blk_len);
        //     printf("[SIM] LDPC BER  : %e\n", dec_err_tot*1.0/(sim_cnt+1)/dsp_blk_len);
        //     printf("[SIM] LDPC FER  : %e\n", cw_fail_tot*1.0/(sim_cnt+1));
        //     printf("[SIM] LDPC MIS  : %e\n", cw_misc_tot*1.0/(sim_cnt+1));
        //     printf("[SIM] MCRC FER  : %e\n", cw_mcrc_tot*1.0/(sim_cnt+1));
        //     printf("[SIM] DATA FER  : %e\n", (lba_err_tot+meta_err_tot)*1.0/(sim_cnt+1));
        //     if ((dec_mode == FC_MIX) || (dec_mode==FC_MIX_G2))
        //         printf("[SIM] ECC Retry Rate: %e\n", ldec_tot*1.0/(sim_cnt+1));
        // } 

        long int total_fail_cnt = 0;
        if (sim_cnt%500==29)
        {
            MPI_Allreduce(&sim_err, &total_fail_cnt, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (total_fail_cnt >= max_err_num)
                break;
            if (rank == 0)
            {
                std::cout << "Completed " << (sim_cnt+1)*size
                  << " codewords for all MPI processes, and total " << total_fail_cnt
                  << " errors @";
                print_time();
            }
        }
    }

    printf("[SIM] rank %d Finish simulation @ ", rank);
    print_time();
    MPI_Barrier(MPI_COMM_WORLD);

    printf("--------------------------------------------------------\n");
    printf("[STATISTICS] rank : %d, Total simulated number: %ld\n", rank, sim_cnt+1);
    printf("[STATISTICS] rank : %d, RAW  BER: %e\n", rank, raw_err_tot*1.0/(sim_cnt+1)/dsp_blk_len);
    printf("[STATISTICS] rank : %d, LDPC BER: %e\n", rank, dec_err_tot*1.0/(sim_cnt+1)/dsp_blk_len);
    printf("[STATISTICS] rank : %d, LDPC FER: %e\n", rank, cw_fail_tot*1.0/(sim_cnt+1));
    printf("[STATISTICS] rank : %d, LDPC MIS: %e\n", rank, cw_misc_tot*1.0/(sim_cnt+1));
    printf("[STATISTICS] rank : %d, MCRC FER: %e\n", rank, cw_mcrc_tot*1.0/(sim_cnt+1));
    printf("[STATISTICS] rank : %d, DATA FER: %e\n", rank, (lba_err_tot+meta_err_tot)*1.0/(sim_cnt+1));
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i=0; i<sim_pckt->ldec_max_itr; i++)
        ldec_itr_tot += (i+1)*ldec_cnvg_itr[i];

    if (rank == 0)
        MPI_Reduce(MPI_IN_PLACE, ldec_cnvg_itr, sim_pckt->ldec_max_itr, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    else
        MPI_Reduce(ldec_cnvg_itr, NULL, sim_pckt->ldec_max_itr, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    printf("before allreduce,   ldec_itr_tot = %ld at rank = %d \n", ldec_itr_tot, rank);
    printf("before allreduce,   sim_cnt = %ld at rank = %d \n", sim_cnt+1, rank);
    printf("before allreduce,   ldec_tot = %ld at rank = %d \n", ldec_tot+1, rank);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &ldec_itr_tot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &sim_cnt, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &ldec_tot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &cw_fail_tot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &raw_err_tot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("---------------------------------------------------------\n");
        printf("[STATISTICS] Total simulated number for all ranks: %ld\n", sim_cnt+size);
        printf("[STATISTICS] ALL RANK LDPC RBER: %e\n", raw_err_tot*1.0/(sim_cnt+size)/dsp_blk_len);
        printf("[STATISTICS] ALL RANK LDPC FER: %e\n", cw_fail_tot*1.0/(sim_cnt+size));
        printf("[STATISTICS] ALL RANK Retry Decoder average iterations: %f\n", ldec_itr_tot*1.0/ldec_tot);
        printf("[STATISTICS] ALL RANK Iteration distribution as below: \n");

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
            printf("-------------------------------------------------------------------------------------\n");
        }
        printf("[SIM] ALL RANK Finish simulation @ ");
        print_time();
    }

    // clean up
    sim_pckt->mcrc_clean();
    sim_pckt->rand_clean();
    sim_pckt->dsp_pckt_clean();
    sim_pckt->ch_llr_clean();
    free(ldec_cnvg_itr);
    free(vref);

    MPI_Finalize();
    return 0;
}



void read_arg(int argc, char **argv)
{
    char *ch_sel;
    char junk;

    if ( !(config_file=argv[1]) || !(ch_sel=argv[2]))
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
        if ((!argv[3] || sscanf(argv[3], "%f%c", &ch_para, &junk) != 1))
            print_usage();
    }
}

void print_usage()
{
    printf("Usage: ./ssd_fc.out sim_mode config_file ch_model ch_para(OPTIONAL)\n");
    printf("  CH_MODE: CLEAN, AWGN, BSC, ERR_INJ, MAX_ERR, ALL_ZERO\n");
    printf("  CH_PARA: not needed for CLEAN channel\n");
    printf("           SNR for AWGN/ALL_ZERO\n");
    printf("           RBER for BSC channel\n");
    printf("           Number of injected bit errors for ERR_INJ channel\n");
    printf("           Max number of injected bit errors for MAX_ERR channel\n");
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

    fscanf(fp, "%d", &bytes_of_userdata);
    fgets(str_tmp, 800, fp);

    fscanf(fp, "%d", &bytes_of_parity);
    fgets(str_tmp, 800, fp);

    dsp_meta_size = dsp_meta_size*8; // byte --> bit
    dsp_lba_size = dsp_lba_size*8; // byte --> bit
    dsp_lba_len = dsp_lba_size*dsp_lba_num;
    dsp_src_len = dsp_meta_size + dsp_lba_len;
    // MCRC generator always appends 4 bytes (32 bits) after the DSP source data.
    // For IBEX sc=512 configs, `bytes_of_userdata` must therefore be >= (DSP in data bytes + 4).
    // Some cnfg files mistakenly set `bytes_of_userdata == DSP in data bytes` and will corrupt heap in `mcrc_gen()`.
    // const int dsp_src_bytes = dsp_src_len / 8;
    // const int mcrc_bytes = 4;
    // if (bytes_of_userdata < dsp_src_bytes + mcrc_bytes) {
    //     printf("[CFG Warning] bytes_of_userdata=%dB is smaller than DSP in data=%dB + %dB MCRC; auto-fix to %dB\n",
    //            bytes_of_userdata, dsp_src_bytes, mcrc_bytes, dsp_src_bytes + mcrc_bytes);
    //     bytes_of_userdata = dsp_src_bytes + mcrc_bytes;
    // }

    // 1. h_m (row number of base matrix)
    fscanf(fp, "%d", &h_m);
    fgets(str_tmp, 800, fp);

    // 2. h_n (column number of base matrix)
    fscanf(fp, "%d", &h_n);
    fgets(str_tmp, 800, fp);

    // 3. h_sc (expansion factor of base matrix)
    fscanf(fp, "%d", &h_sc);
    fgets(str_tmp, 800, fp);

    h_m = (bytes_of_parity+63)/64;
    h_n = (bytes_of_userdata+63)/64 + h_m;
    printf("LDPC DEBUG %d %d\n", h_m, h_n);

    // 6. max_sim_num
    fscanf(fp, "%d", &max_sim_num);
    fgets(str_tmp, 800, fp);

    sim_step = (max_sim_num>=1000) ? max_sim_num/10 : 10000;

    // 7. max_err_num
    fscanf(fp, "%d", &max_err_num);
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

    // 19. NAND read number
    fscanf(fp, "%d", &sd_num);
    fgets(str_tmp, 800, fp);

    rd_num = sd_num;

    // 20. read reference voltage
    vref = (float*)calloc(rd_num, sizeof(*vref));
    for (int i=0; i<rd_num; i++)
    {
        fscanf(fp, "%f", vref+i);
    }
    fgets(str_tmp, 800, fp);

    // 21. HD LLRs for retry decoder
    fscanf(fp, "%f", &hd0_llr);
    fscanf(fp, "%f", &hd1_llr);
    fgets(str_tmp, 800, fp);

    // 23. CRC-32 polynomial value for MCRC and MPCRC
    fscanf(fp, "%lx", &crc32_poly_val);
    fgets(str_tmp, 800, fp);

    // 24. 32-bit LFSR polynomial value for randomizer
    fscanf(fp, "%lx", &rand32_poly_val);
    fgets(str_tmp, 800, fp);

    fclose(fp);

    // print out simulation configuration
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        printf("---------------Simulation Configuration----------------\n");
        printf("Userdata %dB, Parity %dB, CW %dB (Info %dB + Pad %dB)\n", bytes_of_userdata, bytes_of_parity, dsp_blk_len/8, dsp_info_len/8, dsp_pad_len/8);
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

        printf(" Retry decoder: max_iter = %d\n", ldec_max_itr);
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

}
