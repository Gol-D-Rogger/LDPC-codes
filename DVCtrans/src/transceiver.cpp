#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rand.h"
#include "finite_lib.h"
#include "transceiver.h"

//############################
// w--> Length of information block
// x--> Length of data block
// y --> CH selection
// z --> CH parameter
//############################
void ch_packet::ch_config(int w, int x, ch_model y, float z)
{
    info_len = w;
    blk_len = x;
    ch_sel = y;
    
    if (ch_sel == CLEAN)
    {
        printf("[CH_TRX] Clean channel selected.\n");
    }
    else if (ch_sel == AWGN)
    {
        snr = z;
        snr_code = snr - 10 * log10(blk_len * 1.0 / info_len);
        awgn_sigma = pow(10, -snr_code * 1.0 / 20) / sqrt(2);
        printf("[CH_TRX] AWGN channel selected with SNR = %f dB.\n", snr);
    }
    else if (ch_sel == BSC)
    {
        ber = z;
        printf("[CH_TRX] BSC channel selected with BER = %f.\n", ber);
    }
    else if (ch_sel == ERR_INJ)
    {
        err_num = (int)z;
        printf("[CH_TRX] Error Injection channel selected with %d errors.\n", err_num);
    }
    else if (ch_sel == MAX_ERR)
    {
        err_num = (int)z;
        printf("[CH_TRX] Maximum Error Injection channel selected with %d errors.\n", err_num);
    }
    else if (ch_sel == ALL_ZERO)
    {
        snr = z;
        snr_code = snr - 10 * log10(blk_len * 1.0 / info_len);
        awgn_sigma = pow(10, -snr_code * 1.0 / 20) / sqrt(2);
        printf("[CH_TRX] ALL_ZERO channel selected (forcing zero codeword) with SNR = %f dB.\n", snr);
    }
}

// Allocate memory for the channel packet
void ch_packet::ch_pckt_alloc()
{
    tx_blk = (char*)calloc(blk_len, sizeof(*tx_blk));
    rx_blk = (float*)calloc(blk_len, sizeof(*rx_blk));
    det_blk = (char*)calloc(blk_len, sizeof(*det_blk));
    rd_blk = (char*)calloc(blk_len*128, sizeof(*rd_blk));
}

void ch_packet::ch_pckt_clean()
{
    free(rx_blk);
    rx_blk = NULL;
    free(det_blk);
    det_blk = NULL;
    free(tx_blk);
    tx_blk = NULL;
    free(rd_blk);
    rd_blk = NULL;

#ifdef _ASIC_DUMP
    free(sd_blk);
    sd_blk = NULL;
#endif
}

// Allocate AWGN related memory
void ch_packet::ch_llr_alloc()
{
    llr_tbl = (float *)calloc(128,sizeof(*llr_tbl));
    split_bin = (char *)calloc(128,sizeof(*split_bin));
    bin_id = (char *)calloc(128,sizeof(*bin_id));
    llr_bin = (float *)calloc(128,sizeof(*llr_bin));
    vref_bin = (float *)calloc(128,sizeof(*vref_bin));
    bin_distr = (int *)calloc(128,sizeof(*bin_distr));
}

void ch_packet::ch_llr_clean()
{
    free(llr_tbl);
    free(split_bin);
    free(bin_id);
    free(llr_bin);
    free(vref_bin);
    free(bin_distr);
}

// transmit data with selected channel model
void ch_packet::ch_transmit()
{
    int *err_vec;
    int *err_pos;
    if (ch_sel == ALL_ZERO)
    {
        for (int i = 0; i < blk_len; i++)
        {
            tx_blk[i] = 0;
        }
    }

    if ((ch_sel == ERR_INJ) || (ch_sel == MAX_ERR))
    {
        if (ch_sel == ERR_INJ)
            raw_err_num = err_num;
        else if (ch_sel == MAX_ERR)
            raw_err_num = rand_int(err_num + 1);

        err_vec = (int *)calloc(blk_len, sizeof(*err_vec));
        err_pos = (int *)calloc(raw_err_num, sizeof(*err_pos));

        randomBinError(err_vec, err_pos, blk_len, raw_err_num);
    }

    for (int i = 0; i < blk_len; i++)
    {
        if (ch_sel == CLEAN)
        {
            rx_blk[i] = -1 * (tx_blk[i] * 2.0 - 1);
        }
        else if ((ch_sel == AWGN) || (ch_sel == ALL_ZERO))
        {
            rx_blk[i] = -1 * (tx_blk[i] * 2.0 - 1) + awgn_sigma * rand_gaussian();
        }
        else if (ch_sel == BSC)
        {
            rx_blk[i] = (rand_uniform() >= ber) ? -1 * (tx_blk[i] * 2.0 - 1) : (tx_blk[i] * 2.0 - 1);
        }
        else if (ch_sel == ERR_INJ || ch_sel == MAX_ERR)
        {
            rx_blk[i] = (err_vec[i] == 0) ? -1 * (tx_blk[i] * 2.0 - 1) : (tx_blk[i] * 2.0 - 1);
        }
    }

    if ((ch_sel == ERR_INJ) || (ch_sel == MAX_ERR))
    {
        free(err_vec);
        free(err_pos);
    }
}

void ch_packet::ch_llr_gen(int rd_num, float *vref, float llr0, float llr1, int tot_num, int frac_num)
{
    float awgn_sigma_sqaure = awgn_sigma * awgn_sigma;
    float llr_max = (pow(2, tot_num - 1)/pow(2, frac_num)); 
    float llr_min = -(pow(2, tot_num - 1)/pow(2, frac_num));

    // Step 1. Decide the bin to be split and also determine the bin index in ascending order
    // > Vref, keep the original order
    // <= Vref, new bin index, i.e., rd_index + 1
    for (int rd_indx = 0; rd_indx < rd_num; rd_indx++)
    {
        // 1st bin split
        if (rd_indx == 0)
        {
            split_bin[0] = 0;
            vref_bin[0] = vref[0];
            bin_id[0] = 1;
            bin_id[1] = 0;
        }
        // split the rightest bin with largest Vref
        else if (vref[rd_indx] > vref_bin[rd_indx - 1])
        {
            split_bin[rd_indx] = bin_id[rd_indx];
            bin_id[rd_indx + 1] = bin_id[rd_indx];
            bin_id[rd_indx] = rd_indx + 1;
            vref_bin[rd_indx] = vref[rd_indx];
        }
        // other cases
        else
        {
            for (int i = 0; i < rd_indx; i++)
            {
                // search from the leftest bin with smallest Vref
                if (vref[rd_indx] < vref_bin[i])
                {
                    split_bin[rd_indx] = bin_id[i];

                    for (int j = rd_indx; j > i; j--)
                        bin_id[j + 1] = bin_id[j];
                    bin_id[i] = rd_indx + 1;

                    for (int j=(rd_indx -1); j>=i; j--)
                        vref_bin[j +1] = vref_bin[j];
                    vref_bin[i] = vref[rd_indx];

                    break;
                }
            }
        }
    }

    if (rd_num == 1)
    {
        llr_bin[0] = llr1;
        llr_bin[1] = llr0;
    }
    // LDPC soft decoding
    // Use LLR of middile point for Vref bin
    else if (rd_num > 1)
    {
        llr_bin[0] = (float)Sat_Quan((double)(2*(-0.1 + vref_bin[0]) / awgn_sigma_sqaure), 
                                    llr_max, llr_min, tot_num, frac_num);
        llr_bin[rd_num] = (float)Sat_Quan((double)(2*(0.1 + vref_bin[rd_num - 1]) / awgn_sigma_sqaure), 
                                    llr_max, llr_min, tot_num, frac_num);

        for (int i = 1; i < rd_num; i++)
        {
            llr_bin[i] = (float)Sat_Quan((double)((vref_bin[i] + vref_bin[i - 1]) / awgn_sigma_sqaure), 
                                    llr_max, llr_min, tot_num, frac_num);
        }
    }

    printf("[LLRGEN] VREF <%5.2f: BIN ID %2d -- %7.4f(%d)\n", 
            vref_bin[0], bin_id[0], llr_bin[0], (int)(llr_bin[0] * pow(2, frac_num)));

    for (int i = 1; i < rd_num; i++)
    {
        printf("[LLRGEN] VREF <%5.2f: BIN ID %2d -- %7.4f(%d)\n", 
            vref_bin[i], bin_id[i], llr_bin[i], (int)(llr_bin[i] * pow(2, frac_num)));        
    }
    printf("[LLRGEN] VREF <%5.2f: BIN ID %2d -- %7.4f(%d)\n", 
            vref_bin[rd_num-1], bin_id[rd_num], llr_bin[rd_num], (int)(llr_bin[rd_num] * pow(2, frac_num)));


    // gen LLR table
    for (int i = 0; i < rd_num; i++)
    {
        for (int j = 0; j <= rd_num; j++)
        {
            if( bin_id[j] == i)
            {
                llr_tbl[i] = llr_bin[j];
            }
        }
    }
}

void ch_packet::ch_detector(int rd_num, float *vref)
{
    int err_cnt = 0;
    int split_bin_id = 0;

    // only bin ID 0 without any read
    for (int i = 0; i < blk_len; i++)
        det_blk[i] = 0;

    for (int rd_indx = 0; rd_indx < rd_num; rd_indx++)
    {
        // read data and update bin
        for (int i = 0; i < blk_len; i++)
        {
            rd_blk[rd_indx * blk_len + i] = (rx_blk[i] >= vref[rd_indx]) ? 0 : 1;
            if ((det_blk[i]==split_bin[rd_indx]) && (rd_blk[rd_indx * blk_len + i]==1))
            {
                det_blk[i] = rd_indx+1;
            }
        }
    }

    // Bin distribution
    for (int i = 0; i < 128; i++)
        bin_distr[i] = 0;

    for (int i = 0; i < blk_len; i++)
        bin_distr[det_blk[i]]++;

#ifdef _CH_DEBUG
    printf("[SIM] Bin Distribution:\n");
    printf("[SIM] ");
    for (int i = 0; i <= rd_num; i++)
    {
        printf("%d(bin %d) ", bin_distr[bin_id[i]], bin_id[i]);
    }
    printf("\n");
#endif

    // check error count
    for (int i = 0; i < blk_len; i++)
    {
        if (tx_blk[i] != ((rx_blk[i] >= 0) ? 0 : 1))
        {
            err_cnt++;
#ifdef _CH_DEBUG
        printf("[SIM]    CH error @ bit %d\n", i);
#endif
        }
    }

#ifdef _CH_DEBUG
    printf("[SIM DEBUG] Total %d CH errors detected.\n", err_cnt);
    if ((ch_sel == ERR_INJ) & (err_cnt != raw_err_num))
    {
        printf("[SIM]!!!! Error numbers mismatch (%d/%d) !!!!\n", err_cnt, raw_err_num);
    }
    if ((ch_sel == MAX_ERR) & (err_cnt != raw_err_num))
    {
        printf("[SIM]!!!! Error numbers mismatch (%d/%d) !!!!\n", err_cnt, raw_err_num);
    }
#endif
    raw_err_num = err_cnt;
}


