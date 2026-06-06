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
void ch_packet::ch_config(int w, int x, ch_model y, float z, int width)
{
    info_len = w;
    blk_len = x;
    ch_sel = y;
    VN_BITS = width;
    
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
#ifdef _ASIC_DUMP
    sd_blk = (char*)calloc(blk_len, sizeof(*sd_blk));
#endif
}

void ch_packet::ch_pckt_clean()
{
    free(rx_blk);
    rx_blk = NULL;
    free(det_blk);
    det_blk = NULL;
    free(tx_blk);
    tx_blk = NULL;

#ifdef _ASIC_DUMP
    free(sd_blk);
    sd_blk = NULL;
#endif
}

// Allocate AWGN related memory
void ch_packet::ch_llr_alloc(enum sd_mode sd_sel, int sd_bit, float *vref_in)
{
    sd_type = sd_sel;
    sd_num = sd_bit;

    if (sd_type == MANUAL)
    {    
        rd_num = sd_num;
        bin_num = rd_num + 1;
    }
    else if (sd_type == DIRECT)
    {
        rd_num = sd_num;
        bin_num = pow(2, sd_num) - 1;
    }
    else if (sd_type == VENDOR0 || sd_type == VENDOR1)
    {
        rd_num = 2 * sd_num - 1;
        bin_num = pow(2, sd_num);
    }

    vref = (float *)calloc(rd_num, sizeof(*vref));
    for (int i = 0; i < rd_num; i++)
    {
        vref[i] = vref_in[i];
    }
    
    llr_tbl = (float *)calloc(bin_num, sizeof(*llr_tbl));
    bin_split = (char *)calloc(rd_num, sizeof(*bin_split));
    bin_asc_ord = (char *)calloc(rd_num+1, sizeof(*bin_asc_ord));
    llr_asc_ord = (float *)calloc(rd_num+1, sizeof(*llr_asc_ord));
    vref_asc_ord = (float *)calloc(rd_num, sizeof(*vref_asc_ord));
    sd_asc_ord = (char *)calloc(rd_num+1, sizeof(*sd_asc_ord));
    bin_distr = (int *)calloc(bin_num, sizeof(*bin_distr));
}

void ch_packet::ch_llr_clean()
{
    free(vref);
    free(llr_tbl);
    free(bin_split);
    free(bin_asc_ord);
    free(llr_asc_ord);
    free(vref_asc_ord);
    free(sd_asc_ord);
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

void ch_packet::ch_llr_gen(float hd0_llr, float hd1_llr, int llr_tot_bit, int llr_frac_bit)
{
    float awgn_sigma_sq = awgn_sigma * awgn_sigma;
    llr_tot_num = llr_tot_bit;
    llr_frac_num = llr_frac_bit;

    llr_max = (pow(2, llr_tot_num - 1) - 1) / pow(2, llr_frac_num);
    llr_min = -llr_max;

    if (sd_type == DIRECT)
    {
        for (int i = (1 - pow(2, sd_num -1)); i < pow(2, sd_num - 1); i++)
        {
            llr_tbl[(int)(i + pow(2, sd_num - 1) - 1)] = (float)(1.0 * i / pow(2, llr_frac_bit));
        }

        max_llr_bin = bin_num - 1;

        printf("[LLRGEN] Direct mode LLR table:\n");
        for (int i = 0; i < bin_num; i++)
        {
            printf("[LLRGEN] BIN ID %2d -- %7.4f\n", i, llr_tbl[i]);
        }
    }
    else // LLR table for real case
    {
        // Step 1. Decide the bin to be split and also determine the bin index in ascending order
        // > Vref, keep the original order
        // <= Vref, new bin index, i.e., rd_index + 1
        for (int rd_indx = 0; rd_indx < rd_num; rd_indx++)
        {
            // 1st bin split
            if (rd_indx == 0)
            {
                bin_split[0] = 0;
                vref_asc_ord[0] = vref[0];
                bin_asc_ord[0] = 1;
                bin_asc_ord[1] = 0;
            }
            // split the rightest bin with largest Vref
            else if (vref[rd_indx] > vref_asc_ord[rd_indx - 1])
            {
                bin_split[rd_indx] = bin_asc_ord[rd_indx];
                bin_asc_ord[rd_indx + 1] = bin_asc_ord[rd_indx];
                bin_asc_ord[rd_indx] = rd_indx + 1;
                vref_asc_ord[rd_indx] = vref[rd_indx];
            }
            // other cases
            else
            {
                for (int i = 0; i < rd_indx; i++)
                {
                    // search from the leftest bin with smallest Vref
                    if (vref[rd_indx] < vref_asc_ord[i])
                    {
                        bin_split[rd_indx] = bin_asc_ord[i];

                        // reorder bin and verf
                        for(int j = rd_indx; j >= i; j--)
                        {
                            bin_asc_ord[j + 1] = bin_asc_ord[j];
                        }
                        bin_asc_ord[i] = rd_indx + 1;

                        for (int j = (rd_indx - 1); j >= i; j--)
                        {
                            vref_asc_ord[j + 1] = vref_asc_ord[j];
                        }
                        vref_asc_ord[i] = vref[rd_indx];
                        break;
                    }
                }
            }
        }

        for (int i = 0; i < rd_num; i++)
        {
            printf("[LLRGEN] Vref[%d] = %f\n", i, vref_asc_ord[i]);
        }

        // Step 2. Determine the LLR in ascending order of Vref
        // LDPC hard decoding
        if (rd_num == 1)
        {
            llr_asc_ord[0] = hd1_llr;
            llr_asc_ord[1] = hd0_llr;
        }
        // LDPC soft decoding
        // Use LLR of middile point for Vref bin
        else 
        {
            llr_asc_ord[0] = (float)Sat_Quan((double)((-1+vref_asc_ord[0]) / awgn_sigma_sq), 
                                            llr_max, llr_min, llr_tot_num, llr_frac_num);
            llr_asc_ord[rd_num] = (float)Sat_Quan((double)((1+vref_asc_ord[rd_num - 1]) / awgn_sigma_sq),
                                            llr_max, llr_min, llr_tot_num, llr_frac_num);

            for (int i = 1; i < rd_num; i++)
            {
                llr_asc_ord[i] = (float)Sat_Quan((double)((vref_asc_ord[i] + vref_asc_ord[i - 1]) / (awgn_sigma_sq)),
                                            llr_max, llr_min, llr_tot_num, llr_frac_num);
            }
        }

        printf("\n");
        for (int i = 0; i < rd_num; i++)
        {
            printf("[LLRGEN] VREF <%5.2f: BIN ID %2d -- %7.4f(%d)\n", 
                vref_asc_ord[i], bin_asc_ord[i], llr_asc_ord[i], (int)(llr_asc_ord[i] * pow(2, llr_frac_num)));
        }

        // Step 3. Determine soft data from NAND
        int **nand_read;

        nand_read = (int **)calloc(rd_num+1, sizeof(*nand_read));
        for (int i = 0; i <= rd_num; i++)
        {
            nand_read[i] = (int *)calloc(rd_num, sizeof(*nand_read[i]));
        }

        if (sd_type == MANUAL)
        {
            for (int i = 0; i < rd_num; i++)
            {
                for (int j = 0; j < rd_num; j++)
                {
                    nand_read[i][j] = (vref_asc_ord[i] <= vref[j]) ? 1 : 0;
                }
            }

            for (int j = 0; j < rd_num; j++)
            {
                nand_read[rd_num][j] = 0;
            }

            for (int i = 0; i <= rd_num; i++)
            {
                sd_asc_ord[i] = nand_read[i][0];
                for (int j = 1; j < rd_num; j++)
                {
                    sd_asc_ord[i] = sd_asc_ord[i] * 2 + nand_read[i][j];
                }
            }
        }
        else
        {
            for (int i = 0; i < rd_num; i++)
            {
                for (int j = 0; j < rd_num; j++)
                {
                    nand_read[i][j] = (vref_asc_ord[i] <= vref_asc_ord[j]) ? 1 : 0;
                }
            }

            for (int j = 0; j < rd_num; j++)
            {
                nand_read[rd_num][j] = 0;
            }
            
            for (int i = 0; i <= rd_num; i++)
            {
                int hard_bit = nand_read[i][sd_num - 1];
                int soft_bits = 0;

                for (int j = 0; j < sd_num - 1; j++)
                {
                    soft_bits = soft_bits * 2;

                    if (nand_read[i][j] == nand_read[i][rd_num -1 -j])
                    {
                        if (sd_type == VENDOR1)
                        {
                            soft_bits += 1;
                        }
                    }
                    else {
                        if (sd_type == VENDOR0)
                        {
                            soft_bits += 1;
                        }
                    }
                }

                sd_asc_ord[i] = soft_bits * 2 + hard_bit;
                bin_asc_ord[i] = sd_asc_ord[i];
            }
        }

        for (int i = 0; i < rd_num; i++)
        {
            free(nand_read[i]);
        }
        free(nand_read);


        // Step 4. Determine the LLR Table
        if (sd_num == 1)
        {
            printf("[LLRGEN] LLR table for hard decoding:\n");
        }
        else if (sd_type == MANUAL)
        {
            printf("[LLRGEN] LLR table for manual mode:\n");
        }
        else if (sd_type == VENDOR0)
        {
            printf("[LLRGEN] LLR table for vendor0 mode:\n");
        }
        else if (sd_type == VENDOR1)
        {
            printf("[LLRGEN] LLR table for vendor1 mode:\n");
        }

        printf("[LLRGEN] VREF <%5.2f: BIN ID %2d -- %7.4f(%d)\n", 
            vref_asc_ord[0], bin_asc_ord[0], llr_asc_ord[0], (int)(llr_asc_ord[0] * pow(2, llr_frac_num)));

        for (int i = 1; i < rd_num; i++)
        {
            printf("[LLRGEN] VREF <%5.2f: BIN ID %2d -- %7.4f(%d)\n", 
                vref_asc_ord[i], bin_asc_ord[i], llr_asc_ord[i], (int)(llr_asc_ord[i] * pow(2, llr_frac_num)));
        }

        printf("[LLRGEN] VREF <%5.2f: BIN ID %2d -- %7.4f(%d)\n", 
            vref_asc_ord[rd_num - 1], bin_asc_ord[rd_num], llr_asc_ord[rd_num], (int)(llr_asc_ord[rd_num] * pow(2, llr_frac_num)));

        // gen LLR table
        for (int i = 0; i < bin_num; i++)
        {
            for (int j = 0; j <= rd_num; j++)
            {
                if( bin_asc_ord[j] == i)
                {
                    llr_tbl[i] = llr_asc_ord[j];
                }
            }
        }

        max_llr_bin = bin_asc_ord[rd_num];
    }

#ifdef _CH_DEBUG
    printf("[LLRGEN] Bin %d has the max LLR=%f\n", max_llr_bin, llr_tbl[max_llr_bin]);
#endif
}

void ch_packet::ch_detector()
{
    int err_cnt = 0;

    for (int i = 0; i < blk_len; i++)
    {
        if (sd_type == DIRECT)
        {
            det_blk[i] = (int)((Sat_Quan((double)(2 * rx_blk[i] / (awgn_sigma * awgn_sigma)), 
                            llr_max, llr_min, llr_tot_num, llr_frac_num) - llr_min) * pow(2, llr_frac_num));
        }
        else
        {
#ifdef _ASIC_DUMP
            sd_blk[i] = sd_asc_ord[rd_num];
#endif            
            det_blk[i] = bin_asc_ord[rd_num];

            for (int rd_indx = 0; rd_indx < rd_num; rd_indx++)
            {
                if (rx_blk[i] < vref_asc_ord[rd_indx])
                {
#ifdef _ASIC_DUMP
                    sd_blk[i] = sd_asc_ord[rd_indx];
#endif
                    det_blk[i] = bin_asc_ord[rd_indx];
                    break;
                }           
            }
        }

        // check error count
        if (tx_blk[i] != ((rx_blk[i] >= 0) ? 0 : 1))
        {
            err_cnt++;
#ifdef _CH_DEBUG
            printf("[SIM] CH error @ bit %d\n", i);
#endif
        }  
    }

    raw_err_num = err_cnt;

    if (sd_type != DIRECT)
    {
#ifdef _CH_DEBUG
        // bin distribution
        for (int i = 0; i < bin_num; i++)
        {
            bin_distr[i] = 0;
        }

        for (int i = 0; i < blk_len; i++)
        {
            bin_distr[det_blk[i]]++;
        }
#endif        

#ifdef _CH_DEBUG
        printf("[SIM] Bin distribution:\n");
        printf("[SIM] ");
        for (int i = 0; i < bin_num; i++)
        {
            printf("%d(bin %d) ", bin_distr[bin_asc_ord[i]], bin_asc_ord[i]);
        }
        printf("\n");
#endif
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
}
