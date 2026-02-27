#pragma once
#include <cstddef>
#ifndef _CH_PACKET_H
#define _CH_PACKET_H

enum ch_model {CLEAN=0, AWGN=1, BSC=2, ERR_INJ=3, MAX_ERR=4, ALL_ZERO=5};

struct ch_packet{
    ch_packet()
        : ch_sel(CLEAN),
          info_len(0),
          blk_len(0),
          snr(0.0f),
          ber(0.0f),
          raw_err_num(0),
          tx_blk(NULL),
          rx_blk(NULL),
          det_blk(NULL),
          sd_blk(NULL),
          snr_code(0.0f),
          awgn_sigma(0.0f),
          err_num(0),
          sd_num(0),
          rd_num(0),
          bin_num(0),
          vref(NULL),
          llr_tbl(NULL),
          bin_split(NULL),
          bin_asc_ord(NULL),
          vref_asc_ord(NULL),
          llr_asc_ord(NULL),
          sd_asc_ord(NULL),
          bin_distr(NULL),
          max_llr_bin(0),
          llr_tot_num(0),
          llr_frac_num(0),
          llr_max(0.0f),
          llr_min(0.0f)
    {}
    enum ch_model ch_sel;
    int    info_len;
    int    blk_len;
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

    void ch_config(int, int, enum ch_model, float);
    void ch_pckt_alloc();
    void ch_pckt_clean();
    void ch_llr_alloc(int, float*);
    void ch_llr_clean();
    void ch_transmit();
    void ch_detector();
    void ch_llr_gen(float, float, int, int);
};

#endif // _CH_PACKET_H
