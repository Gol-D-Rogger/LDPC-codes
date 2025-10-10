#pragma once
#ifndef _CH_PACKET_H
#define _CH_PACKET_H

enum sd_mode {DIRECT=0, MANUAL=1, VENDOR0=2, VENDOR1=3};
enum ch_model {CLEAN=0, AWGN=1, BSC=2, ERR_INJ=3, MAX_ERR=4};

struct ch_packet{
    enum ch_model ch_sel;
    int    info_len;
    int    blk_len;
    int    real_len;
    float  snr;
    float  ber;
    int    raw_err_num;
    char  *tx_blk;
    float *rx_blk;
    char  *det_blk;
    char  *sd_blk;

    // CH
    float  snr_code;
    float  awgn_sigma;
    int    err_num;

    // detector
    enum sd_mode sd_type;   // soft decision type, select LLR GEN Mode
    int    sd_num;          // soft decision bit num
    int    rd_num;          // read num or Vref num
    int    bin_num;         // 量化区间的num
    float *vref;
    float *llr_tbl;
    char  *bin_split;       // 标记哪个旧的Bin被新的Vref分割了
    char  *bin_asc_ord;
    float *vref_asc_ord;
    float *llr_asc_ord;
    char  *sd_asc_ord;
    int   *bin_distr;

    int   max_llr_bin;

    int   llr_tot_num;
    int   llr_frac_num;
    float llr_max;
    float llr_min;

    void ch_config(int, int, int, enum ch_model, float);
    void ch_pckt_alloc();
    void ch_pckt_clean();
    void ch_llr_alloc(enum sd_mode, int, float*);
    void ch_llr_clean();
    void ch_transmit();
    void ch_detector();
    void ch_llr_gen(float, float, int, int);
};

#endif // _CH_PACKET_H