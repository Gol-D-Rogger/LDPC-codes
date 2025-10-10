// Description:
// 1. DSP packet allocation/clean
// 2. Data Randomization
// 3. Media CRC

#include "ldpc_codec.h"

#ifndef _FC_DSP_H
#define _FC_DSP_H

// FC_SKIP: skip LDPC decoding
// FC_FDEC: fast LDPC decoding
// FC_MIX:  hybrid LDPC decoding: fast + retry
// FC_RDEC: retry LDPC decoding only
enum fc_dec_mode
{
    FC_SKIP = 0,
    FC_FDEC,
    FC_RDEC,
    FC_MIX,
    FC_FDEC_G2,
    FC_MIX_G2,
    FC_TBFDEC,
};

struct dsp_packet : ldpc_packet
{
    // randminzer & MCRC configuration
    int *rand_poly;
    int *mcrc_poly;
    int rseed_num;
    unsigned int *rseed_tbl; // randomizer seed table
    int agit_num;
    unsigned int *agit_tbl; // agitator table
    unsigned int agit_page_cycle;
    int remap_ofst_init;
    int remap_ofst_diff_init;
    int cycle_init;
    int cycle_diff_init;
    int d1, d2;
    int page_cw_num;

    // data block
    int meta_size;
    int lba_size;   // Logical Block Addressing size
    int lba_num;
    int lba_len;    // CW total LBA length
    int src_len;   // CW Source block length

    char *wr_meta_blk; // meta data block
    char *wr_lba_blk;  // LBA data block
    char *src_blk;    // source data block
    char *wr_rand_blk; // write path random data block
    char *rd_rand_blk; // read path random data block
    char *rcv_blk; // recovered data block
    char *rd_meta_blk; // read path meta data block
    char *rd_lba_blk; 

    // CW status
    int mcrc_err;
    int lba_err;
    int meta_err;
    int rdec_used;

    // alloc/cleanup packet
    void dsp_pckt_alloc();
    void dsp_pckt_clean();

    // QC LDPC config & clean up
    void dfmt_config(int msize, int lsize, int lnum);
    void rand_config(long);
    void rand_clean();
    void mcrc_config(long);
    void mcrc_clean();

    // META data handling
    void meta_insert();
    void meta_extract();
    // randomizer 
    unsigned int rand_seed_sel(int, int, int);
    void rand_encoder(int, int, int);
    void rand_decoder(int, int, int);
    void rand_eng(unsigned int, char*, char*, int);
    
    // MCRC
    void mcrc_gen();
    void mcrc_chk();
    void crc32_eng(char*, char*, int);

    // ECC
    void ecc_encoder();
    void ecc_decoder(enum fc_dec_mode fc_dec_mode);

    // dump all data
    void packet_dump(int pckt_id);
};

#endif // _FC_DSP_H