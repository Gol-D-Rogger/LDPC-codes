#pragma once
#include <cstddef>
#ifndef _CH_PACKET_H
#define _CH_PACKET_H

enum sd_mode {DIRECT=0, MANUAL=1, VENDOR0=2, VENDOR1=3};
enum ch_model {CLEAN=0, AWGN=1, BSC=2, ERR_INJ=3, MAX_ERR=4, ALL_ZERO=5};

struct ch_packet
{
    enum ch_model ch_sel;
    int    info_len;
    int    blk_len;
    float  snr;
    float  ber;
    int    raw_err_num;
    char  *tx_blk;
    float *rx_blk;
    char  *det_blk;
    char  *rd_blk;

    // CH
    float  snr_code;
    float  awgn_sigma;
    int    err_num;

    // detector
    float *llr_tbl;
    char  *bin_id;
    char  *split_bin;
    float *vref_bin;
    float *llr_bin;
    int   *bin_distr;

    void ch_config(int, int, enum ch_model, float);
    void ch_pckt_alloc();
    void ch_pckt_clean();
    void ch_llr_alloc();
    void ch_llr_clean();
    void ch_transmit();
    void ch_detector(int, float*);
    void ch_llr_gen(int, float*, float, float, int, int);
};

#endif // _CH_PACKET_H
