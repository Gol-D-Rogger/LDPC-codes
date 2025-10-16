#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "vec_op.h"
#include "finite_lib.h"
#include "ldpc_codec.h"
#include "fc_dsp.h"

// randomizer cofiguration
void dsp_packet::rand_config(long rand_p_val)
{
    if (rand_p_val >= pow(2, 33))
    {
        printf("[DSP ERROR] Randomizer only supports degree 32 polynomial! (%lx)\n", rand_p_val);
        exit(1);
    }

    rand_poly = (int *)calloc(32, sizeof(*rand_poly));

    //generate randomizer polynomial
    int i = 0;
    long tmp = rand_p_val;
    while (tmp!=1)
    {
        rand_poly[31-i] = tmp%2;
        tmp = tmp / 2;
        i++;
    }

    if (i!=32)
    {
        printf("[DSP ERROR] Randomizer only supports degree 32 polynomial! (%lx)\n", rand_p_val);
        exit(1);
    }

    // randomizer seed selection
    remap_ofst_init = 0;
    remap_ofst_diff_init = 1;
    cycle_init = 4;
    cycle_diff_init = 0;
    d1 = 0;
    d2 = 0;
    page_cw_num = 8;

    // randomizer seed table
    rseed_num = 67;
    rseed_tbl = (unsigned int *)calloc(rseed_num, sizeof(*rseed_tbl));
    rseed_tbl[0] = 0xab11cc9b;
    rseed_tbl[1] = 0x2e1da61b;
    rseed_tbl[2] = 0xf8372be7;
    rseed_tbl[3] = 0xf2cdf49f;
    rseed_tbl[4] = 0xd3d22f7f;
    rseed_tbl[5] = 0x5be09d1c;
    rseed_tbl[6] = 0xf98e137e;
    rseed_tbl[7] = 0x6f670345;
    rseed_tbl[8] = 0xd91b5d72;
    rseed_tbl[9] = 0xedec5375;
    rseed_tbl[10]= 0x8a50e1d2;
    rseed_tbl[11]= 0x167637db;
    rseed_tbl[12]= 0x3988baed;
    rseed_tbl[13]= 0x5075b36e;
    rseed_tbl[14]= 0x27ea9ab7;
    rseed_tbl[15]= 0x49d5a266;
    rseed_tbl[16]= 0x948e296e;
    rseed_tbl[17]= 0xbfd518f3;
    rseed_tbl[18]= 0xdf01ec4a;
    rseed_tbl[19]= 0x686e1bb1;
    rseed_tbl[20]= 0xfa7ee2c9;
    rseed_tbl[21]= 0x7397c076;
    rseed_tbl[22]= 0xd275da40;
    rseed_tbl[23]= 0xcf428cef;
    rseed_tbl[24]= 0xe1ab245b;
    rseed_tbl[25]= 0x57e26a5e;
    rseed_tbl[26]= 0x2ac04837;
    rseed_tbl[27]= 0x99a115a6;
    rseed_tbl[28]= 0xd3c92d47;
    rseed_tbl[29]= 0x196d8a26;
    rseed_tbl[30]= 0x9a1a1aa9;
    rseed_tbl[31]= 0xb8259e70;
    rseed_tbl[32]= 0x547badc1;
    rseed_tbl[33]= 0xd3ed8ca8;
    rseed_tbl[34]= 0xa34a4e94;
    rseed_tbl[35]= 0xe8bc32f0;
    rseed_tbl[36]= 0x3063210f;
    rseed_tbl[37]= 0x4a5808d8;
    rseed_tbl[38]= 0x94772728;
    rseed_tbl[39]= 0x7d832e2b;
    rseed_tbl[40]= 0xbe1626a2;
    rseed_tbl[41]= 0xf046bea7;
    rseed_tbl[42]= 0x8c451487;
    rseed_tbl[43]= 0xb8bcd15f;
    rseed_tbl[44]= 0x555fdde8;
    rseed_tbl[45]= 0x86006513;
    rseed_tbl[46]= 0x4df771ed;
    rseed_tbl[47]= 0xe391fbee;
    rseed_tbl[48]= 0xdbc1fd98;
    rseed_tbl[49]= 0x164cfd33;
    rseed_tbl[50]= 0x2eae0831;
    rseed_tbl[51]= 0x57077929;
    rseed_tbl[52]= 0x72cfee1f;
    rseed_tbl[53]= 0x7b2a7fcc;
    rseed_tbl[54]= 0x39d47304;
    rseed_tbl[55]= 0xafa49758;
    rseed_tbl[56]= 0xc9c5e5a3;
    rseed_tbl[57]= 0x45a5fa49;
    rseed_tbl[58]= 0xe544e572;
    rseed_tbl[59]= 0x7e8f78a6;
    rseed_tbl[60]= 0x42074860;
    rseed_tbl[61]= 0x6a183bbc;
    rseed_tbl[62]= 0x4d797600;
    rseed_tbl[63]= 0x03be9b49;
    rseed_tbl[64]= 0xb80dfafb;
    rseed_tbl[65]= 0xf2dad720;
    rseed_tbl[66]= 0x9aaaa2c7;

    // Agitation value
    agit_num = 11;
    agit_page_cycle = 72;
    agit_tbl = (unsigned int *)calloc(agit_num, sizeof(*agit_tbl));
    agit_tbl[0] = 0x00000000;
    agit_tbl[1] = 0x84c099ab;
    agit_tbl[2] = 0x77498833;
    agit_tbl[3] = 0xc3072608;
    agit_tbl[4] = 0xa615813c;
    agit_tbl[5] = 0x7bdfed05;
    agit_tbl[6] = 0xd1dedc6f;
    agit_tbl[7] = 0xfed76b76;
    agit_tbl[8] = 0xf3593c12;
    agit_tbl[9] = 0x2ace4776;
    agit_tbl[10]= 0x44d46da8;
} // rand_config

// randomizer clean up
void dsp_packet::rand_clean()
{
    free(rand_poly);
    free(rseed_tbl);
    free(agit_tbl);
}

// MCRC clean up
void dsp_packet::mcrc_clean()
{
    free(mcrc_poly);
}

// MCRC config
void dsp_packet::mcrc_config(long mcrc_p_val)
{
    if (mcrc_p_val >= pow(2, 33))
    {
        printf("[DSP ERROR] MCRC only supports degree 32 polynomial! (%lx)\n", mcrc_p_val);
        exit(1);
    }

    mcrc_poly = (int *)calloc(32, sizeof(*mcrc_poly));

    // MCRC polynomial
    int i = 0;
    long tmp = mcrc_p_val;
    while (tmp!=1)
    {
        mcrc_poly[31-i] = tmp%2;
        tmp = tmp / 2;
        i++;
    }

    if (i!=32)
    {
        printf("[DSP ERROR] MCRC only supports degree 32 polynomial! (%lx)\n", mcrc_p_val);
        exit(1);
    }
}

// DFT config
// msize --> Meta data size
// lsize --> LBA size
// lnum --> Number of LBAs in a CW
void dsp_packet::dfmt_config(int msize, int lsize, int lnum)
{
    meta_size = msize;
    lba_size = lsize;
    lba_num = lnum;

    lba_len = lba_size * lba_num;
    src_len = meta_size + lba_len;
}

void dsp_packet::dsp_pckt_alloc()
{
    ldpc_pckt_alloc();

    wr_meta_blk = (char *)calloc(meta_size, sizeof(*wr_meta_blk));
    wr_lba_blk = (char *)calloc(lba_len, sizeof(*wr_lba_blk));
    rd_meta_blk = (char *)calloc(meta_size, sizeof(*rd_meta_blk));
    rd_lba_blk = (char *)calloc(lba_len, sizeof(*rd_lba_blk));
    src_blk = (char *)calloc(src_len, sizeof(*src_blk));
    wr_rand_blk = (char *)calloc(src_len, sizeof(*wr_rand_blk));
    rd_rand_blk = (char *)calloc(src_len, sizeof(*rd_rand_blk));
    rcv_blk = (char *)calloc(src_len, sizeof(*rcv_blk));

    lba_err = 0;
    meta_err = 0;
    mcrc_err = 0;
}

void dsp_packet::dsp_pckt_clean()
{
    ldpc_pckt_clean();
    free(wr_meta_blk);
    wr_meta_blk = NULL;
    free(wr_lba_blk);
    wr_lba_blk = NULL;
    free(rd_meta_blk);
    rd_meta_blk = NULL;
    free(rd_lba_blk);
    rd_lba_blk = NULL;
    free(src_blk);
    src_blk = NULL;
    free(wr_rand_blk);
    wr_rand_blk = NULL;
    free(rd_rand_blk);
    rd_rand_blk = NULL;
    free(rcv_blk);
    rcv_blk = NULL;
}

void dsp_packet::meta_insert()
{
    vec_copy(wr_lba_blk, src_blk, 0, 0, lba_len);
    vec_copy(wr_meta_blk, src_blk, 0, lba_len, meta_size);
}

void dsp_packet::meta_extract()
{
    vec_copy(rcv_blk, rd_lba_blk, 0, 0, lba_len);
    vec_copy(rcv_blk, rd_meta_blk, lba_len, 0, meta_size);
}

// randomizer seed selector
// paddr --> page address
// cw_indx --> codeword index within a page
// pe_cnt --> PE cycle
unsigned int dsp_packet::rand_seed_sel(int paddr, int cw_indx, int pe_cnt)
{
    // remap calculate offset
    int cycle = cycle_init;
    int cycle_diff = cycle_diff_init;
    int remap_ofst = remap_ofst_init;
    int remap_ofst_diff = remap_ofst_diff_init;
    int j = 0;

    for (int i=0; i<= paddr; i++)
    {
        if (j== cycle)
        {
            cycle_diff += d2;
            cycle += cycle_diff;
            remap_ofst += remap_ofst_diff;
            remap_ofst_diff += d1;
            j = 0;
        }
        
        j++;
    }

    // rand seed table select
    int rseed_tbl_sel;
    int rseed_tbl_val;

    rseed_tbl_sel = (paddr * page_cw_num + cw_indx + remap_ofst + pe_cnt) % rseed_num;
    rseed_tbl_val = rseed_tbl[rseed_tbl_sel];

    // agitation table select
    int agit_tbl_val;
    agit_tbl_val = agit_tbl[paddr%agit_page_cycle];

    // randomizer seed
    return (agit_tbl_val ^ rseed_tbl_val);
}

// randomizer engine (LFSR)
// rseed --> randomizer seed
// din_blk --> input data block
// dout_blk --> output data block
// blk_len --> block length
void dsp_packet::rand_eng(unsigned int rseed, char *din_blk, char *dout_blk, int blk_len)
{
    char *lfsr = (char *)calloc(32, sizeof(*lfsr));
    
    for (int i = 0; i < 32; i++)
    {
        lfsr[31-i] = rseed % 2;
        rseed = rseed / 2;
        i++;
    }

    for (int i = 0; i < blk_len; i = i + 8)
    {
        dout_blk[i] = (din_blk[i] + lfsr[15]) % 2;      // X17
        dout_blk[i + 1] = (din_blk[i + 1] + lfsr[19]) % 2; // X13
        dout_blk[i + 2] = (din_blk[i + 2] + lfsr[20]) % 2; // X12
        dout_blk[i + 3] = (din_blk[i + 3] + lfsr[23]) % 2; // X9
        dout_blk[i + 4] = (din_blk[i + 4] + lfsr[26]) % 2; // X6
        dout_blk[i + 5] = (din_blk[i + 5] + lfsr[27]) % 2; // X5
        dout_blk[i + 6] = (din_blk[i + 6] + lfsr[29]) % 2; // X3
        dout_blk[i + 7] = (din_blk[i + 7] + lfsr[30]) % 2; // X2

        // move LFSR
        for (int i = 0; i < 32 - 1; i++)
        {
            if (rand_poly[i] == 1)
            {
                lfsr[i] = (lfsr[0] + lfsr[i + 1]) % 2;
            }
            else
            {
                lfsr[i] = lfsr[i + 1];
            }
        }
        lfsr[31] = lfsr[0];
    }
}

// randomizer
void dsp_packet::rand_encoder(int paddr, int cw_indx, int pe_cnt)
{
    unsigned int rseed = rand_seed_sel(paddr, cw_indx, pe_cnt);
    rand_eng(rseed, src_blk, wr_rand_blk, src_len);
}

// randomizer decoder
void dsp_packet::rand_decoder(int paddr, int cw_indx, int pe_cnt)
{
    unsigned int rseed = rand_seed_sel(paddr, cw_indx, pe_cnt);
    rand_eng(rseed, rd_rand_blk, rcv_blk, src_len);
}

// CRC32 generator
void dsp_packet::crc32_eng(char *din_blk, char *crc_blk, int blk_len)
{
    int feedback;

    for (int i = 0; i < blk_len; i++)
    {
        feedback = (crc_blk[0] + din_blk[i]) % 2;

        for (int j = 0; j < 31; j++)
        {
            crc_blk[j] = (crc_blk[j+1] + feedback * mcrc_poly[j]) % 2;
        }

        crc_blk[31] = feedback;
    }
}

// 4Byte  MCRC generator
// wr_rand_blk --> usr_blk
void dsp_packet::mcrc_gen()
{
    char *crc_poly = (char *)calloc(32, sizeof(*crc_poly));

    crc32_eng(wr_rand_blk, crc_poly, src_len);

    vec_copy(crc_poly, usr_blk, 0, src_len, 32);
    vec_copy(wr_rand_blk, usr_blk, 0, 0, src_len);

    free(crc_poly);
}

// 4Byte MCRC checker
// dec_blk --> rd_rand_blk
void dsp_packet::mcrc_chk()
{
    char *crc_poly = (char *)calloc(32, sizeof(*crc_poly));

    crc32_eng(dec_blk, crc_poly, src_len);

    vec_copy(dec_blk, rd_rand_blk, 0, 0, src_len);

    // checker
    mcrc_err = 0;
    for (int i = 0; i < 32; i++)
    {
        if (dec_blk[i+src_len] != crc_poly[i])
        {
            mcrc_err = 1;
        }
    }

    free(crc_poly);
}

void dsp_packet::ecc_encoder()
{
    mcrc_gen();
    ldpc_encoder();
}

void dsp_packet::ecc_decoder(enum fc_dec_mode fc_mode)
{
    enum dec_model dec_mode;
    
    rdec_used = 0;

    // decoding skip
    if (fc_mode == FC_SKIP)
    {
        dec_mode = SKIP;
        ldpc_decoder(dec_mode);
        mcrc_chk();
    }

    // fast decoding
    if ((fc_mode == FC_FDEC) || (fc_mode == FC_MIX))
    {
        dec_mode = BF_P3;
        ldpc_decoder(dec_mode);
        mcrc_chk();
    }
    if (fc_mode == FC_TBFDEC)
    {
        dec_mode = TBFDEC;
        ldpc_decoder(dec_mode);
        mcrc_chk();
    }
    if ((fc_mode == FC_FDEC_G2) || (fc_mode == FC_MIX_G2))
    {
        dec_mode = BF_G2;
        ldpc_decoder(dec_mode);
        mcrc_chk();
    }

    // retry decoding only
    if (fc_mode == FC_RDEC)
    {
        rdec_used = 1;
        dec_mode = LAYER;
        ldpc_decoder(dec_mode);
        mcrc_chk();
    }

    // mix decoding after fast decoding failed
    if ((fc_mode == FC_MIX) || (fc_mode == FC_MIX_G2))
    {
        if ((cw_fail == 1) || (mcrc_err == 1))
        {
#ifdef _DSP_DEBUG
            if (cw_fail == 1)
            {
                printf("[DSP DEBUG] BF decoding failed, start layer decoding...\n");
            }
            else if (mcrc_err == 1)
            {
                printf("[DSP DEBUG] BF mis-correction and MCRC check failed, start layer decoding...\n");
            }
#endif
            rdec_used = 1;

            dec_mode = LAYER;
            ldpc_decoder(dec_mode);
            mcrc_chk();
        }
    }
}

// Dump the whole packet
void dsp_packet::packet_dump(int pckt_id)
{
    FILE *fp;

    // meta data
    char WR_META_DATA[50];
    sprintf(WR_META_DATA, "./output/wr_meta_data_%d.txt", pckt_id);
    fp = fopen(WR_META_DATA, "w");
    vec_print(wr_meta_blk, meta_size, fp);
    fclose(fp);

    char RD_META_DATA[50];
    sprintf(RD_META_DATA, "./output/rd_meta_data_%d.txt", pckt_id);
    fp = fopen(RD_META_DATA, "w");
    vec_print(rd_meta_blk, meta_size, fp);
    fclose(fp);

    // LBA
    char WR_LBA_DATA[50];
    sprintf(WR_LBA_DATA, "./output/wr_lba_data_%d.txt", pckt_id);
    fp = fopen(WR_LBA_DATA, "w");
    vec_print(wr_lba_blk, lba_len, fp);
    fclose(fp);

    char RD_LBA_DATA[50];
    sprintf(RD_LBA_DATA, "./output/rd_lba_data_%d.txt", pckt_id);
    fp = fopen(RD_LBA_DATA, "w");
    vec_print(rd_lba_blk, lba_len, fp);
    fclose(fp);

    // source data
    char SRC_DATA[50];
    sprintf(SRC_DATA, "./output/src_data_%d.txt", pckt_id);
    fp = fopen(SRC_DATA, "w");
    vec_print(src_blk, src_len, fp);
    fclose(fp);

    // write path randomized data
    char WR_RAND_DATA[50];
    sprintf(WR_RAND_DATA, "./output/wr_rand_data_%d.txt", pckt_id);
    fp = fopen(WR_RAND_DATA, "w");
    vec_print(wr_rand_blk, src_len, fp);
    fclose(fp);

    // ECC user data
    char USER_DATA[50];
    sprintf(USER_DATA, "./output/user_data_%d.txt", pckt_id);
    fp = fopen(USER_DATA, "w");
    vec_print(usr_blk, info_len, fp);
    fclose(fp);

    // ECC encoded data
    char ENC_DATA[50];
    sprintf(ENC_DATA, "./output/enc_data_%d.txt", pckt_id);
    fp = fopen(ENC_DATA, "w");
    vec_print(tx_blk, real_len, fp);
    fclose(fp);

    // ECC decoded data
    char DEC_DATA[50];
    sprintf(DEC_DATA, "./output/dec_data_%d.txt", pckt_id);
    fp = fopen(DEC_DATA, "w");
    vec_print(dec_blk, real_len, fp);
    fclose(fp);

    // read path randomized data
    char RD_RAND_DATA[50];
    sprintf(RD_RAND_DATA, "./output/rd_rand_data_%d.txt", pckt_id);
    fp = fopen(RD_RAND_DATA, "w");
    vec_print(rd_rand_blk, src_len, fp);
    fclose(fp);

    // recovered data
    char RCV_DATA[50];
    sprintf(RCV_DATA, "./output/rcv_data_%d.txt", pckt_id);
    fp = fopen(RCV_DATA, "w");
    vec_print(rcv_blk, src_len, fp);
    fclose(fp);
} // packet_dump
