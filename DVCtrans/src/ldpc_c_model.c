#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rand.h"
#include "alloc.h"
#include "vec_op.h"
#include "mod2sparse.h"
#include "transceiver.h"
#include "ldpc_codec.h"

struct ldpc_packet *sim_pckt;

extern "C"
void ldpc_config(int h_m,
                 int h_n,
                 int h_sc,
                 int h_st,
                 int h_wt,
                 int info_num,
                 int ch_mode,
                 int ch_para,
                 int fdec_max_itr,
                 int ldec_max_itr,
                 int llr0,
                 int llr1,
                 int alpha,
                 int finite_q_num,
                 int finite_r_num,
                 int sd_num,
                 svOpenArrayHandle v_ref_sv,
                 int debug,
                 int sdlite_llr_config,
                 int sdlite_llr0,
                 int sdlite_llr1,
                 int sdlite_llr2,
                 int sdlite_llr3
                )
{
    sim_pckt = &ldpc_pckt;
    char pchk_file[200];
    float m_ch_para;
    float m_alpha;
    float m_llr0;
    float m_llr1;
    int i;
    int *v_ref;

    int h_k = h_n - h_m;
    int pad_num = h_k*h_sc-info_num;
    int blk_num = h_n*h_sc-pad_num;

    sprintf(pchk_file, "ldpc_h_%d_%d_%d_%d_%d.txt", h_m, h_n, h_sc, h_wt, h_wt);
    m_ch_para = (float)ch_para/100.0;
    m_alpha = (float)alpha/100.0;
    m_llr0 = (float)llr0/16;
    m_llr1 = (float)llr1/16;

    //for soft
    if(sd_num >= 2)
    {
        v_ref = (int *)svGetArrayPtr(v_ref_sv);
        for(i=0; i<sd_num; i++)
            vref[i] = (float)v_ref[i]/1000;
    }
    else {
        vref[0] = 0;
    }

    // CH configuration
    if (ch_mode==0)
        sim_pckt->ch_config(info_num, blk_num, CLEAN,   m_ch_para);
    else if (ch_mode==1)
        sim_pckt->ch_config(info_num, blk_num, AWGN,    m_ch_para);
    else if (ch_mode==2)
        sim_pckt->ch_config(info_num, blk_num, BSC,     m_ch_para);
    else if (ch_mode==3)
        sim_pckt->ch_config(info_num, blk_num, ERR_INJ, m_ch_para);
    else if (ch_mode==4)
        sim_pckt->ch_config(info_num, blk_num, MAX_ERR, m_ch_para);

    // LDPC configuration
    sim_pckt->ldpc_config(h_m, h_n, h_sc, h_st, h_wt, pchk_file);

    sim_pckt->out_len = blk_num;
    // packet allocation
    sim_pckt->ldpc_pckt_alloc();
    sim_pckt->ch_llr_alloc();

    //1st reference: sd_num, 1: for hard, others for soft, at most 128 reads
    // vref: -1 to 1; usually use -0.1 to 0.1
    sim_pckt->ch_llr_gen(sd_num, vref, m_llr0, m_llr1, finite_q_num-1, 4);
    sim_pckt->ldpc_dec_config(fdec_max_itr, ldec_max_itr, m_alpha, 1, finite_q_num, finite_r_num, 4,
                              sdlite_llr_config, sdlite_llr0, sdlite_llr1, sdlite_llr2, sdlite_llr3);
}                

extern "C"
void ldpc_cleanup()
{
    sim_pckt->ldpc_clean();
    sim_pckt->ldpc_pckt_clean();
    sim_pckt->ch_llr_clean();
}

extern "C"
void ldpc_enc(svOpenArrayHandle usr_data_sv,
              svOpenArrayHandle enc_data_sv,
              int debug)
{
    int i, j, k;
    int *usr_data;
    int *enc_data;
    unsigned int tmp;

    // get data pointer
    usr_data = (int *)svGetArrayPtr(usr_data_sv);
    enc_data = (int *)svGetArrayPtr(enc_data_sv);

    j = 0;
    for (i = 0; i < sim_pckt->info_len; i++)
    {
        if (j==0)
            tmp = (unsigned int)(usr_data[i/32]);

        sim_pckt->usr_blk[i] = (tmp >> (31-j))%2;
        j++;
        j = j%32;
    }

    // encoder
    sim_pckt->ldpc_encoder();

    // data out
    j = 0;
    k = 0;
    tmp = 0;

    for (i = 0; i < sim_pckt->out_len; i++)
    {
        tmp = tmp*2 + sim_pckt->tx_blk[i];
        j++;
        j=j%32;

        if ((j==0) || (i==sim_pckt->out_len-1))
        {
            enc_data[k] = tmp;
            k++;
            tmp = 0;
        }
    }
}              

extern "C"
void ldpc_dec(svOpenArrayHandle det_data_sv,
              svOpenArrayHandle dec_data_sv,
              int *dec_unc_sv,
              int *init_synd_wt_sv,
              int *dec_cnvg_itr_sv,
              int *dec_cnvg_col_sv,
              int *fina_synd_wt_sv,
              int dec_mode,
              int sd_num,
              int debug,
              int h_n,
              int h_sc
)
{
    int i, j, k;
    int *det_data;
    int *dec_data;
    unsigned int tmp;

    det_data = (int *)svGetArrayPtr(det_data_sv);
    dec_data = (int *)svGetArrayPtr(dec_data_sv);

    // data in
    if(sd_num < 2)
    {
        j = 0;
        for (i = 0; i < sim_pckt->blk_len; i++)
        {
            if (j==0)
                tmp = (unsigned int)(det_data_sv[i/32]);

            sim_pckt->det_blk[i] = (tmp >> (31-j))%2;
            j++;
            j = j%32;
        }
    }
    else
    {
        for (i = 0; i < sim_pckt->blk_len; i++)
            sim_pckt->det_blk[i] = (unsigned int)(det_data[i]);
    }

    // decoder
    if (dec_mode==0)
        sim_pckt->ldpc_decoder(BF_P3);
    if (dec_mode==1)
        sim_pckt->ldpc_decoder(LAYER);
    if (dec_mode==2)
        sim_pckt->ldpc_decoder(TBFDEC);

    *dec_unc_sv = sim_pckt->cw_fail;
    *init_synd_wt_sv = sim_pckt->init_synd_wt;
    *dec_cnvg_itr_sv = sim_pckt->cnvg_itr;
    *dec_cnvg_col_sv = sim_pckt->cnvg_lyr;
    *fina_synd_wt_sv = sim_pckt->fina_synd_wt;
    vpi_printf("c_debug %d, %d, %d, %d, %d\n", *dec_unc_sv, *init_synd_wt_sv, *dec_cnvg_itr_sv, *dec_cnvg_col_sv, *fina_synd_wt_sv);

    // data out
    j = 0;
    k = 0;
    tmp = 0;
    for (i = 0; i < sim_pckt->out_len; i++)
    {
        tmp = tmp*2 + sim_pckt->dec_blk[i];
        j++;
        j=j%32;

        if ((j==0) || (i==sim_pckt->out_len-1))
        {
            dec_data[k] = tmp;
            tmp = 0;
            k++;
        }
    }
    while (i%32 !=0)
        i++;

    for (;i<h_n*h_sc;i+=32)
    {
        dec_data[k] = 0;
        k++;
    }
}

// CH update
extern "C"
void ch_update(int ch_mode,
               int ch_para,
               int debug)
{
    float m_ch_para;

    m_ch_para = (float)ch_para/100.0;

    // CH configuration
    if (ch_mode==0)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->out_len, CLEAN,   m_ch_para);
    else if (ch_mode==1)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->out_len, AWGN,    m_ch_para);
    else if (ch_mode==2)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->out_len, BSC,     m_ch_para);
    else if (ch_mode==3)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->out_len, ERR_INJ, m_ch_para);
    else if (ch_mode==4)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->out_len, MAX_ERR, m_ch_para);
}               

// error injection
extern "C"
void ch_err_inj(svOpenArrayHandle tx_data_sv,
                svOpenArrayHandle rx_data_sv,
                int debug)
{
    float vref=0;
    int i, j, k;
    unsigned int tmp;
    int *tx_data;
    int *det_data;

    // get data pointer
    tx_data = (int *)svGetArrayPtr(tx_data_sv);
    det_data = (int *)svGetArrayPtr(rx_data_sv);

    // data in
    j = 0;
    for (i = 0; i < sim_pckt->blk_len; i++)
    {
        if (j==0)
            tmp = (unsigned int)(tx_data[i/32]);

        sim_pckt->tx_blk[i] = (tmp >> (31-j))%2;
        j++;
        j = j%32;
    }

    // error injection
    sim_pckt->ch_transmit();
    sim_pckt->ch_detector(1, &vref);

    // data out
    j = 0;
    k = 0;
    tmp = 0;
    for (i = 0; i < sim_pckt->out_len; i++)
    {
        tmp = tmp*2 + sim_pckt->det_blk[i];
        j++;
        j=j%32;

        if ((j==0) || (i==sim_pckt->out_len-1))
        {
            det_data[k] = tmp;
            tmp = 0;
            k++;
        }
    }
}                

// this function should combine with ldpc_config, which sd_num>=2
// soft error injection
// input data: original data
// output data: soft read data(user mode now)
//              bin(6 bits for each CW's bit)
//              llr table(can get from ldpc_config, but get here is also okay)
//              split bin: same as vref, at most 127
//              bin_id: after vref, bin_id value from left to right
extern "C"
void sd_err_inj(svOpenArrayHandle tx_data_sv,
                   svOpenArrayHandle rx_data_sv,    // bin, CW_size*6
                   svOpenArrayHandle rd_data_sv,    // rd_data, CW_size*rd_num
                   svOpenArrayHandle llr_tbl_sv,    // [128]
                   svOpenArrayHandle split_bin_sv,  // splited bin, for user mode
                   svOpenArrayHandle bin_id_sv,     // BIN ID, from left to right
                   int sd_num,
                   svOpenArrayHandle v_ref_sv,      // at most 127
                   int debug)
{
    float vref[128];
    int i, j, k;
    unsigned int tmp;
    int *tx_data;
    int *bin_data;
    int *rd_data;
    int *llr_tbl;
    int *v_ref;
    int *split_bin;
    int *bin_id;

    tx_data = (int *)svGetArrayPtr(tx_data_sv);
    bin_data = (int *)svGetArrayPtr(rx_data_sv);
    rd_data = (int *)svGetArrayPtr(rd_data_sv);
    llr_tbl = (int *)svGetArrayPtr(llr_tbl_sv);
    split_bin = (int *)svGetArrayPtr(split_bin_sv);
    bin_id = (int *)svGetArrayPtr(bin_id_sv);
    v_ref = (int *)svGetArrayPtr(v_ref_sv);

    if (sd_num >= 2)
    {
        for(i=0; i<sd_num; i++)
            vref[i] = (float)v_ref[i]/1000;
    }
    else {
        return;
    }

    // data in
    j = 0;
    for (i = 0; i < sim_pckt->out_len; i++)
    {
        if (j==0)
            tmp = (unsigned int)(tx_data[i/32]);

        sim_pckt->tx_blk[i] = (tmp >> (31-j))%2;
        j++;
        j = j%32;
    }

    // error injection
    sim_pckt->ch_transmit();
    sim_pckt->ch_detector(sd_num, vref);

    // bin out
    j = 0;
    k = 0;
    tmp = 0;
    for (i = 0; i < sim_pckt->out_len; i++)
    {
        tmp = sim_pckt->det_blk[i];
        bin_data[k] = tmp;
        tmp = 0;
        k++;
    }

    // rd data out
    j = 0;
    k = 0;
    tmp = 0;
    for (i = 0; i < sim_pckt->out_len*sd_num; i++)
    {
        tmp = tmp*2 + sim_pckt->rd_blk[i];
        j++;
        j=j%32;

        if ((j==0) || (i==sim_pckt->out_len*sd_num-1))
        {
            rd_data[k] = tmp;
            tmp = 0;
            k++;
        }
    }

    // get llr_tbl
    for(i=0; i<128; i++)
        llr_tbl[i] = sim_pckt->llr_tbl[i] * 16;

    // get split_bin
    for(i=0; i<128; i++)
        split_bin[i] = sim_pckt->split_bin[i];

    // get bin_id
    for(i=0; i<128; i++)
        bin_id[i] = sim_pckt->bin_id[i];
}           

void sd_err_inj_for_tbdec(svOpenArrayHandle tx_data_sv,
                          svOpenArrayHandle rx_data_sv_hd,    // bin, CW_size*6
                          svOpenArrayHandle rd_data_sv_sf,    // rd_data, CW_size*rd_num
                          svOpenArrayHandle llr_tbl_sv,    // [128]
                          svOpenArrayHandle rx_data_c,      // bin, CW_size*6
                          svOpenArrayHandle v_ref_sv,      // at most 127
                          svOpenArrayHandle bin_map_sv,
                          int data2fdec_core,
                          int debug)
{
    float vref[3];
    int i, j, k;
    unsigned int tmp, tmp0, tmp1;
    int *tx_data, *det_data_hd, *det_data_sf, *det_data_c, *llr_tbl, *v_ref, *bin_map;

    tx_data = (int *)svGetArrayPtr(tx_data_sv);
    det_data_hd = (int *)svGetArrayPtr(det_data_sv_hd);
    det_data_sf = (int *)svGetArrayPtr(rd_data_sv_sf);
    det_data_c = (int *)svGetArrayPtr(rx_data_c);
    llr_tbl = (int *)svGetArrayPtr(llr_tbl_sv);
    bin_map = (int *)svGetArrayPtr(bin_map_sv);
    v_ref = (int *)svGetArrayPtr(v_ref_sv);

    vref[0] = (float) v_ref[0] / (float) 1000;
    vref[1] = (float) v_ref[1] / (float) 1000;
    vref[2] = (float) v_ref[2] / (float) 1000;

    for (i=0; i<128; i++)
        llr_tbl[i] = sim_pckt->llr_tbl[i] * 16;

    int llr_tbl_tmp[4];
    llr_tbl_tmp[0] = llr_tbl[0];
    llr_tbl_tmp[1] = llr_tbl[1];
    llr_tbl_tmp[2] = llr_tbl[2];
    llr_tbl_tmp[3] = llr_tbl[3];

    // min at right([3])
    for (i=0; i<3; i++) {
        for (j=0; j<4-1-i; j++) {
            if (llr_tbl_tmp[j] < llr_tbl_tmp[j+1]) {
                k = llr_tbl_tmp[j];
                llr_tbl_tmp[j] = llr_tbl_tmp[j+1];
                llr_tbl_tmp[j+1] = k;
            }
        }
    }

    for (i=0; i<4; i++) {
        if (llr_tbl[i] == llr_tbl_tmp[0])
            bin_map[i] = 1; // 01
        else if (llr_tbl[i] == llr_tbl_tmp[1])
            bin_map[i] = 0; // 00
        else if (llr_tbl[i] == llr_tbl_tmp[2])
            bin_map[i] = 2; // 10
        else if (llr_tbl[i] == llr_tbl_tmp[3])
            bin_map[i] = 3; // 11
    }

    vpi_printf("llr_tbl_tmp[0~3]: %d   %d   %d   %d\n", llr_tbl_tmp[0], llr_tbl_tmp[1], llr_tbl_tmp[2], llr_tbl_tmp[3]);
    vpi_printf("bin_map[0~3]: %d   %d   %d   %d\n", bin_map[0], bin_map[1], bin_map[2], bin_map[3]);

    // data in
    j = 0;
    for (i = 0; i < sim_pckt->out_len; i++)
    {
        if (j==0)
            tmp = (unsigned int)(tx_data[i/32]);

        sim_pckt->tx_blk[i] = (tmp >> (31-j))%2;
        j++;
        j = j%32;
    }

    sim_pckt->ch_transmit();
    sim_pckt->ch_detector(3, vref, 0);

    // data out
    j = 0;
    k = 0;
    tmp = 0;
    for (i = 0; i < sim_pckt->out_len; i++)
    {
        tmp = sim_pckt->det_blk[i];
        if (i>=sim_pckt->blk_len)
            det_data_c[k] = 0;
        else
            det_data_c[k] = tmp;
        tmp = 0;
        k++;
    }

    j = 0;
    k = 0;
    tmp0 = 0;
    tmp1 = 0;
    int det_blk_bin_map = 0;
    vpi_printf("parameter data2fdec_core is: %d\n", data2fdec_core);
    for (i = 0; i < sim_pckt->out_len; i++)
    {
        if (data2fdec_core == 1)
            det_blk_bin_map = bin_map[sim_pckt->det_blk[i]];
        else
            det_blk_bin_map = sim_pckt->det_blk[i];

        if (i >= sim_pckt->blk_len)
        {
            tmp0=tmp0*2;
            tmp1=tmp1*2;
        }
        else
        {
            if (det_blk_bin_map == 0)
            {
                tmp0 = tmp0 * 2;
                tmp1 = tmp1 * 2;
            }
            else if (det_blk_bin_map == 1)
            {
                tmp0 = tmp0 * 2;
                tmp1 = tmp1 * 2+1;
            }
            else if (det_blk_bin_map == 2)
            {
                tmp0 = tmp0 * 2+1;
                tmp1 = tmp1 * 2;
            }
            else if (det_blk_bin_map == 3)
            {
                tmp0 = tmp0 * 2+1;
                tmp1 = tmp1 * 2+1;
            }
            else {
                vpi_printf("det_blk error\n");
            }
        }

        j++;
        j=j%32;

        if ((j==0) || (i==(sim_pckt->out_len-1)))
        {
            det_data_hd[k] = tmp0;
            det_data_sf[k] = tmp1;
            tmp0 = 0;
            tmp1 = 0;
            k++;
        }
    }
}           
