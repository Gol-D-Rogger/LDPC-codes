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

// 定义全局指针，指向 ldpc_packet 实例
struct ldpc_packet *sim_pckt;

// 声明 DQ 特有的配置函数，如果头文件中未声明
// 假设 ldpc_packet 结构体定义在编译时会包含 ldpc_config_dq
// void ldpc_config_dq(int m, int n, int sc, int st, int wt, int pad_bit, char *pchk_file, char *mask_file);

// 长度记录
static int g_info_len = 0;   // 用户信息长度
static int g_blk_len  = 0;   // 传输块长度（信息+校验，不含 padding 0）
static int g_hm_k     = 0;   // 系统信息长度（补零后）
static int g_hm_m     = 0;   // parity 长度
static int g_pad_num  = 0;   // 补零长度

extern "C"
void ldpc_config(int h_m,
                 int h_n,
                 int h_sc,
                 int h_st,
                 int h_wt,
                 int info_num,
                 int pad_bit,       // 新增参数: pad_bit
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

    // Gen4 DQ 维度（简化版）：扩 1 行 1 列，保留 pad_bit，屏蔽 mask_len
    int bm_m = h_m + 1;
    int bm_n = h_n + 1;
    int mask_len = h_sc - pad_bit;
    g_info_len = info_num;
    g_hm_m = bm_m * h_sc - mask_len;          // parity 长 = h_m*h_sc + pad_bit
    g_hm_k = (bm_n - bm_m) * h_sc;            // 系统信息长度
    g_pad_num = g_hm_k - g_info_len;          // 补零长度
    g_blk_len = g_info_len + g_hm_m;          // 传输块长度（信息 + 校验）

    // 文件名构造 (使用物理维度 bm_m, bm_n)
    sprintf(pchk_file, "ldpc_h_%d_%d_%d_%d_%d.txt", bm_m, bm_n, h_sc, h_wt, h_wt);

    m_ch_para = (float)ch_para/100.0;
    m_alpha = (float)alpha/100.0;
    m_llr0 = (float)llr0/16;
    m_llr1 = (float)llr1/16;

    // 软判决参考电压处理
    if(sd_num >= 2)
    {
        v_ref = (int *)svGetArrayPtr(v_ref_sv);
        for(i=0; i<sd_num; i++)
            vref[i] = (float)v_ref[i]/1000;
    }
    else {
        vref[0] = 0;
    }

    // 信道配置: 使用 blk_len（信息+校验）
    if (ch_mode==0)
        sim_pckt->ch_config(g_info_len, g_blk_len, CLEAN,   m_ch_para);
    else if (ch_mode==1)
        sim_pckt->ch_config(g_info_len, g_blk_len, AWGN,    m_ch_para);
    else if (ch_mode==2)
        sim_pckt->ch_config(g_info_len, g_blk_len, BSC,     m_ch_para);
    else if (ch_mode==3)
        sim_pckt->ch_config(g_info_len, g_blk_len, ERR_INJ, m_ch_para);
    else if (ch_mode==4)
        sim_pckt->ch_config(g_info_len, g_blk_len, MAX_ERR, m_ch_para);

    // LDPC 配置: 当前 DVCtrans 仅有标准接口，使用扩展后的尺寸近似配置
    sim_pckt->ldpc_config(bm_m, bm_n, h_sc, h_st, h_wt, pchk_file);

    // 缓冲分配
    sim_pckt->ldpc_pckt_alloc();
    sim_pckt->ch_llr_alloc();

    // LLR 生成
    sim_pckt->ch_llr_gen(sd_num, vref, m_llr0, m_llr1, finite_q_num-1, 4);
    
    // 解码器配置 (DQ 版参数可能略有不同，这里沿用通用配置)
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

    // 获取数据指针
    usr_data = (int *)svGetArrayPtr(usr_data_sv);
    enc_data = (int *)svGetArrayPtr(enc_data_sv);

    // 1. 解包用户数据 (长度 info_len)
    j = 0;
    for (i = 0; i < sim_pckt->info_len; i++)
    {
        if (j==0)
            tmp = (unsigned int)(usr_data[i/32]);

        sim_pckt->usr_blk[i] = (tmp >> (31-j))%2;
        j++;
        j = j%32;
    }

    // 显式清零 padding 区，防止残留
    for (i = 0; i < sim_pckt->pad_len; i++)
    {
        sim_pckt->usr_blk[sim_pckt->info_len + i] = 0;
    }

    // 注意: padding 填充通常在 ldpc_encoder 内部根据 pad_len 自动处理
    // Gen4 encoder 会处理 enc_di_blk 的 padding

    // 2. 执行编码
    sim_pckt->ldpc_encoder();

    // 3. 打包编码数据
    // 长度使用 blk_len (在 DQ 模式下即为 hm_n, 有效物理长度)
    j = 0;
    k = 0;
    tmp = 0;

    for (i = 0; i < sim_pckt->blk_len; i++)
    {
        tmp = tmp*2 + sim_pckt->tx_blk[i];
        j++;
        j=j%32;

        if ((j==0) || (i==sim_pckt->blk_len-1))
        {
            enc_data[k] = tmp;
            k++;
            tmp = 0;
        }
    }

    if ((sim_pckt->blk_len % 32) == 16) 
    {
        unsigned int last = enc_data[k-1];
        last = last << 16;
        enc_data[k-1] = last;
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

    // 1. 数据输入解包（长度 = blk_len = info + parity）
    if(sd_num < 2)
    {
        j = 0;
        for (i = 0; i < sim_pckt->blk_len; i++)
        {
            if (j==0)
                tmp = (unsigned int)(det_data[i/32]);

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

    // 2. 重构 [Info | Pad | Parity] 到 det_blk 前 hm_k+hm_m 位置
    int info_len = sim_pckt->info_len;
    int hm_k     = sim_pckt->hm_k;     // = info_len + pad_len
    int hm_m     = sim_pckt->hm_m;
    int pad_len  = sim_pckt->pad_len;

    // 使用临时缓冲防止覆盖，并避免超出 det_blk 原有容量
    char *tmp_blk = (char *)calloc(hm_k + hm_m, sizeof(*tmp_blk));
    for (i = 0; i < info_len; i++)
        tmp_blk[i] = sim_pckt->det_blk[i];
    for (i = 0; i < pad_len; i++)
        tmp_blk[info_len + i] = 0;
    for (i = 0; i < hm_m; i++)
        tmp_blk[hm_k + i] = sim_pckt->det_blk[info_len + i];

    // 暂时替换 det_blk 指针和 blk_len 供译码使用
    char *old_det = sim_pckt->det_blk;
    sim_pckt->det_blk = tmp_blk;

    // 3. 执行解码
    if (dec_mode==0)
        sim_pckt->ldpc_decoder(BF_P3);
    if (dec_mode==1)
        sim_pckt->ldpc_decoder(LAYER);
    if (dec_mode==2)
        sim_pckt->ldpc_decoder(TBFDEC);

    // 恢复指针
    sim_pckt->det_blk = old_det;
    free(tmp_blk);

    // 4. 状态输出
    *dec_unc_sv = sim_pckt->cw_fail;
    *init_synd_wt_sv = sim_pckt->init_synd_wt;
    *dec_cnvg_itr_sv = sim_pckt->cnvg_itr;
    *dec_cnvg_col_sv = sim_pckt->cnvg_lyr;
    *fina_synd_wt_sv = sim_pckt->fina_synd_wt;
    vpi_printf("c_debug %d, %d, %d, %d, %d\n", *dec_unc_sv, *init_synd_wt_sv, *dec_cnvg_itr_sv, *dec_cnvg_col_sv, *fina_synd_wt_sv);

    // 4. 数据输出打包
    // Gen4 通常只输出 Info 部分 (info_len)
    // 检查 Gen3 行为: Gen3 输出 blk_len (Info+Parity)
    // 为保持一致性，这里暂且输出 blk_len。如果需要仅输出 Info，需改为循环 info_len。
    // 鉴于 DVC 通常用于 verify 用户数据，且 DQ Parity 很长，建议确认需求。
    // 但根据现有代码结构，使用 sim_pckt->dec_blk (已去除了 padding)，其有效部分长度视 dec_blk 分配而定。
    // 原始代码: dec_data[k]... 循环 blk_len。
    // DQ ldpc_decoder: vec_copy(dec_do_blk, dec_blk, 0, 0, info_len); vec_copy(..., hm_k, info_len, hm_m);
    // 所以 dec_blk 里存的是 [Info, Parity] (无 padding)，总长 info_len + hm_m = blk_len (hm_n)。
    // 所以循环 blk_len 是正确的。
    
    j = 0;
    k = 0;
    tmp = 0;
    for (i = 0; i < sim_pckt->blk_len; i++)
    {
        tmp = tmp*2 + sim_pckt->dec_blk[i];
        j++;
        j=j%32;

        if ((j==0) || (i==sim_pckt->blk_len-1))
        {
            dec_data[k] = tmp;
            tmp = 0;
            k++;
        }
    }

    if ((sim_pckt->blk_len % 32) == 16) 
    {
        unsigned int last = dec_data[k-1];
        last = last << 16;
        dec_data[k-1] = last;
    }

    // 补齐对齐
    while (i%32 !=0)
        i++;

    // 原始代码还有一段补 0 到 h_n*h_sc 的逻辑，DQ 模式下 h_n 传入值可能已变，
    // 或者如果 SV 端 buffer 很大，可能需要清零。保留原逻辑结构但注意 bounds。
    // 这里的 h_n, h_sc 是 ldpc_dec 的入参，可能还是 Config 时的 logical dimension。
    // 为安全起见，防止 SV 越界读取，可以保留补零。
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
    // 使用 sim_pckt->blk_len (已在 config 中设为 hm_n)
    if (ch_mode==0)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->blk_len, CLEAN,   m_ch_para);
    else if (ch_mode==1)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->blk_len, AWGN,    m_ch_para);
    else if (ch_mode==2)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->blk_len, BSC,     m_ch_para);
    else if (ch_mode==3)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->blk_len, ERR_INJ, m_ch_para);
    else if (ch_mode==4)
        sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->blk_len, MAX_ERR, m_ch_para);
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
    // j = 0;
    // k = 0;
    // tmp = 0;
    // for (i = 0; i < sim_pckt->blk_len; i++)
    // {
    //     tmp = tmp*2 + sim_pckt->det_blk[i];
    //     j++;
    //     j=j%32;

    //     if ((j==0) || (i==sim_pckt->blk_len-1))
    //     {
    //         det_data[k] = tmp;
    //         tmp = 0;
    //         k++;
    //     }
    // }

    int n_words = (sim_pckt->blk_len + 31) / 32;
    for (int w = 0; w < n_words; w++) 
    {
      unsigned int word = (unsigned int)tx_data[w];

        for (int b = 0; b < 32; b++) 
        {
            int bit_index = w * 32 + (31 - b);  // 与解包时相反的 MSB 顺序
            if (bit_index >= sim_pckt->blk_len)
                break;  // 超出码字长度的位保持 TX 原值

            unsigned int bit  = (unsigned int)(sim_pckt->det_blk[bit_index] & 1u);
            unsigned int mask = 1u << b;

            if (bit)
                word |= mask;
            else
                word &= ~mask;
        }
        det_data[w] = word;
    }
}                

// soft error injection
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
    sim_pckt->ch_detector(sd_num, vref);

    // bin out
    j = 0;
    k = 0;
    tmp = 0;
    for (i = 0; i < sim_pckt->blk_len; i++)
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
    for (i = 0; i < sim_pckt->blk_len*sd_num; i++)
    {
        tmp = tmp*2 + sim_pckt->rd_blk[i];
        j++;
        j=j%32;

        if ((j==0) || (i==sim_pckt->blk_len*sd_num-1))
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
