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

static ldpc_packet ldpc_pckt;
struct ldpc_packet *sim_pckt;

static int g_configured = 0;

static int g_vn_bits = 8;

// 长度记录（Gen4 口径：DV 侧物理长度域）
static int g_info_len = 0; // 用户信息长度
static int g_blk_len = 0;  // 传输块长度（信息 + parity，不含 padding 0）
static int g_hm_k = 0;     // 系统信息长度（补零后）
static int g_hm_m = 0;     // parity 长度（DV 侧物理）
static int g_hm_n = 0;     // 码长（DV 侧物理）
static int g_pad_len = 0;  // padding 0 长度（shortening）
static int g_bm_m = 0;
static int g_bm_n = 0;

// IBEX BF decoder DV-side parameters (ldpc_decoder_inv style)
static int g_post_iter = 0;
static int g_max_iter = 0;
static int g_nand_strobes = 1;

extern "C" void ldpc_cleanup();

static unsigned int dvc_mask_valid_bits_msb(int total_bits)
{
    const int rem = total_bits % 32;
    if (rem == 0)
        return 0xFFFFFFFFu;
    return 0xFFFFFFFFu << (32 - rem);
}

static void dvc_pack_bits_msb32(const char *bits, int n_bits, unsigned int *words, int n_words)
{
    unsigned int tmp = 0;
    int j = 0;
    int k = 0;
    for (int i = 0; i < n_bits; i++) {
        tmp = tmp * 2u + (unsigned int)(bits[i] & 1);
        j++;
        j = j % 32;
        if ((j == 0) || (i == (n_bits - 1))) {
            if (k < n_words)
                words[k] = tmp;
            k++;
            tmp = 0;
        }
    }
    if (k < n_words) {
        for (; k < n_words; k++)
            words[k] = 0;
    }
}

static void dvc_check_words_against_bits_msb(const char *tag,
                                             const int *words,
                                             const char *bits,
                                             int n_bits,
                                             int debug)
{
    if (debug < 2)
        return;
    const int n_words = (n_bits + 31) / 32;
    if (n_words <= 0)
        return;
    if (n_words > 128) {
        printf("[DVC IBEX WARN] %s pack-check skipped: n_words=%d > 128\n", tag, n_words);
        return;
    }
    unsigned int repack[128] = {0};
    dvc_pack_bits_msb32(bits, n_bits, repack, n_words);
    const unsigned int last_mask = dvc_mask_valid_bits_msb(n_bits);
    for (int w = 0; w < n_words; w++) {
        const unsigned int mask = (w == (n_words - 1)) ? last_mask : 0xFFFFFFFFu;
        const unsigned int a = ((unsigned int)words[w]) & mask;
        const unsigned int b = repack[w] & mask;
        if (a != b) {
            printf("[DVC IBEX WARN] %s unpack/pack mismatch at word=%d mask=%08x in=%08x repack=%08x\n", tag, w, mask,
                   a, b);
            return;
        }
    }
}

static int dvc_report_bit_mismatches_msb(const char *tag,
                                        const char *a_bits,
                                        const char *b_bits,
                                        int n_bits,
                                        int max_list)
{
    int err = 0;
    int listed = 0;
    for (int i = 0; i < n_bits; i++) {
        const int a = a_bits[i] & 1;
        const int b = b_bits[i] & 1;
        if (a != b) {
            err++;
            if (max_list > 0 && listed < max_list) {
                const int w = i / 32;
                const int j = i % 32;
                const int bitpos = 31 - j; // MSB-first
                printf("%s mismatch i=%d (word=%d bit=%d): a=%d b=%d\n", tag, i, w, bitpos, a, b);
                listed++;
            }
        }
    }
    printf("%s mismatch_count=%d / %d bits\n", tag, err, n_bits);
    return err;
}

static int dvc_qc_syndrome_weight_bm(mod2sparse *qc_bm, const char *hard_bits, int bm_m, int cir_sz)
{
    if (!qc_bm || !hard_bits || bm_m <= 0 || cir_sz <= 0)
        return 0;

    int synd_wt = 0;
    char *layer_synd = (char *)calloc(cir_sz, sizeof(*layer_synd));
    char *vn_dec_hd = (char *)calloc(cir_sz, sizeof(*vn_dec_hd));
    char *cn_dec_hd = (char *)calloc(cir_sz, sizeof(*cn_dec_hd));
    if (!layer_synd || !vn_dec_hd || !cn_dec_hd) {
        free(layer_synd);
        free(vn_dec_hd);
        free(cn_dec_hd);
        return 0;
    }

    for (int layer = 0; layer < bm_m; layer++) {
        vec_clr(layer_synd, cir_sz);
        for (mod2entry *e = mod2sparse_first_in_row(qc_bm, layer); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e)) {
            vec_copy((char *)hard_bits, vn_dec_hd, e->col * cir_sz, 0, cir_sz);
            vec_shift(vn_dec_hd, cn_dec_hd, cir_sz, -1 * e->shift);
            vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);
        }
        synd_wt += vec_sum(layer_synd, cir_sz);
    }

    free(layer_synd);
    free(vn_dec_hd);
    free(cn_dec_hd);
    return synd_wt;
}

static unsigned int dvc_sig_fnv1a_bits01(const char *bits, int start, int len)
{
    unsigned int h = 2166136261u;
    for (int i = 0; i < len; i++) {
        h ^= (unsigned int)(bits[start + i] & 1);
        h *= 16777619u;
    }
    return h;
}

static unsigned int dvc_sig_fnv1a_rx_bins01(const ldpc_packet *pckt, int start, int len)
{
    unsigned int h = 2166136261u;
    if (!pckt || !pckt->det_blk) {
        return h;
    }

    const int bin_num = pckt->bin_num;
    const int max_bin = (bin_num > 0) ? (bin_num - 1) : 0;
    for (int i = 0; i < len; i++) {
        int bit = 0;
        if (pckt->llr_tbl && (bin_num > 0)) {
            int bin = (int)(unsigned char)pckt->det_blk[start + i];
            if (bin > max_bin)
                bin = max_bin;
            bit = (pckt->llr_tbl[bin] < 0.0f) ? 1 : 0;
        } else {
            bit = (pckt->det_blk[start + i] & 1);
        }
        h ^= (unsigned int)(bit & 1);
        h *= 16777619u;
    }
    return h;
}

static void dvc_rotate_segment_left(char *buf, int start, int len, int sh)
{
    if (!buf || len <= 0)
        return;
    sh %= len;
    if (sh < 0)
        sh += len;
    if (sh == 0)
        return;

    char *tmp = (char *)calloc(len, sizeof(*tmp));
    if (!tmp)
        return;
    for (int i = 0; i < len; i++)
        tmp[i] = buf[start + ((i + sh) % len)];
    for (int i = 0; i < len; i++)
        buf[start + i] = tmp[i];
    free(tmp);
}

extern "C"
void ldpc_config(int h_m,
                 int h_n,
                 int h_sc,
                 int h_st,
                 int h_wt,
                 int info_num,
                 int pad_bit,
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
                 int sdlite_llr3,
                 int dv_user_data_bytes,
                 int dv_parity_bytes,
                 int post_iter,
                 int max_iter,
                 int nand_strobes,
                 int codeword_4k_8k,
                 int syndrome_weight_thr_qc,
                 int syndrome_weight_thr_post,
                 int early_termination_dis,
                 int post_process_en,
                 int ldpc_decoder_control_likelihood_0,
                 int ldpc_decoder_control_likelihood_1,
                 int ldpc_decoder_control_likelihood_2,
                 int ldpc_decoder_control_likelihood_3,
                 int ldpc_decoder_control_post,
                 unsigned int ldpc_early_term_0,
                 unsigned int ldpc_early_term_1,
                 unsigned int ldpc_early_term_2,
                 unsigned int ldpc_early_term_3,
                 unsigned int ldpc_early_term_4,
                 unsigned int ldpc_early_term_5,
                 unsigned int ldpc_early_term_6)
{
    if (g_configured)
        ldpc_cleanup();

    sim_pckt = &ldpc_pckt;

    const float m_ch_para = (float)ch_para / 1000.0f;
    const float m_alpha = (float)alpha / 1000.0f;
    const float m_llr0 = (float)llr0 / 16.0f;
    const float m_llr1 = (float)llr1 / 16.0f;

    g_bm_m = h_m + 1;
    g_bm_n = h_n + 1;
    const int mask_len = h_sc - pad_bit;
    g_info_len = info_num;
    g_hm_m = g_bm_m * h_sc - mask_len;
    g_hm_n = g_bm_n * h_sc - mask_len;
    g_hm_k = (g_bm_n - g_bm_m) * h_sc;
    g_pad_len = g_hm_k - g_info_len;
    g_blk_len = g_info_len + g_hm_m;

    const int user_data_bytes = g_info_len / 8;
    const int parity_bytes = (g_blk_len - g_info_len) / 8;
    const int blk_bytes = g_blk_len / 8;

    g_post_iter = post_iter;
    g_max_iter = max_iter;
    g_nand_strobes = nand_strobes;

    // MP_Framework: only support IBEX internal configuration path (cir_sz==512).
    // External matrix-file based configuration (cir_sz==256) is removed by design.
    if (h_sc != 512) {
        printf("[DVC IBEX ERROR] Unsupported cir_sz=%d in MP_Framework. Only sc=512 is supported.\n", h_sc);
        return;
    }

    printf("[DVC IBEX] ldpc_config: h_m=%d h_n=%d sc=%d pad_bit=%d => bm_m=%d bm_n=%d\n", h_m, h_n, h_sc, pad_bit,
           g_bm_m, g_bm_n);
    printf("[DVC IBEX] len(bits): info=%d hm_k=%d hm_m=%d hm_n=%d pad_len=%d blk_len=%d\n", g_info_len, g_hm_k, g_hm_m,
           g_hm_n, g_pad_len, g_blk_len);
    printf("[DVC IBEX] len(bytes): user=%d parity=%d blk=%d (VN_BITS=%d)\n", user_data_bytes, parity_bytes, blk_bytes,
           g_vn_bits);
    printf("[DVC IBEX] dv_cfg(bytes): user=%d parity=%d | post_iter=%d max_iter=%d nand_strobes=%d codeword_4k_8k=%d\n",
           dv_user_data_bytes, dv_parity_bytes, post_iter, max_iter, nand_strobes, codeword_4k_8k);

    if ((g_info_len % 8) || ((g_blk_len - g_info_len) % 8) || (g_blk_len % 8)) {
        printf("[DVC IBEX WARN] byte alignment: info%%8=%d parity%%8=%d blk%%8=%d\n", g_info_len % 8,
               (g_blk_len - g_info_len) % 8, g_blk_len % 8);
    }

    if ((dv_user_data_bytes > 0) && (dv_user_data_bytes != user_data_bytes)) {
        printf("[DVC IBEX WARN] dv_user_data_bytes mismatch: dv=%d calc=%d (g_info_len=%d)\n", dv_user_data_bytes,
               user_data_bytes, g_info_len);
    }
    if ((dv_parity_bytes > 0) && (dv_parity_bytes != parity_bytes)) {
        printf("[DVC IBEX WARN] dv_parity_bytes mismatch: dv=%d calc=%d (hm_m(bits)=%d)\n", dv_parity_bytes, parity_bytes,
               g_hm_m);
    }

    int rd_num = sd_num;
    if (rd_num < 1)
        rd_num = 1;
    if (rd_num > 127)
        rd_num = 127;
    float vref[128] = {0.0f};
    if (sd_num >= 2) {
        int *v_ref = (int *)svGetArrayPtr(v_ref_sv);
        for (int i = 0; i < rd_num; i++)
            vref[i] = (float)v_ref[i] / 1000.0f;
    } else {
        vref[0] = 0.0f;
    }

    // CH configuration (Gen4/DVC style mapping)
    enum ch_model ch_sel = CLEAN;
    if (ch_mode == 0)
        ch_sel = CLEAN;
    else if (ch_mode == 1)
        ch_sel = AWGN;
    else if (ch_mode == 2)
        ch_sel = BSC;
    else if (ch_mode == 3)
        ch_sel = ERR_INJ;
    else if (ch_mode == 4)
        ch_sel = MAX_ERR;
    else
        printf("[DVC IBEX WARN] unknown ch_mode=%d, fallback CLEAN\n", ch_mode);

    sim_pckt->ch_config(g_info_len, g_blk_len, ch_sel, m_ch_para, g_vn_bits);
    if (debug >= 1) {
        if (ch_sel == BSC)
            printf("[DVC IBEX] ch_config: mode=%d(BSC) ch_para=%d => ber=%f\n", ch_mode, ch_para, sim_pckt->ber);
        else if ((ch_sel == ERR_INJ) || (ch_sel == MAX_ERR))
            printf("[DVC IBEX] ch_config: mode=%d(ERR_INJ/MAX_ERR) ch_para=%d => err_num=%d\n", ch_mode, ch_para,
                   sim_pckt->err_num);
        else if ((ch_sel == AWGN) || (ch_sel == ALL_ZERO))
            printf("[DVC IBEX] ch_config: mode=%d(AWGN/ALL_ZERO) ch_para=%d => snr=%f\n", ch_mode, ch_para, sim_pckt->snr);
        else
            printf("[DVC IBEX] ch_config: mode=%d(CLEAN) ch_para=%d\n", ch_mode, ch_para);
    }

    // LDPC configuration
    sim_pckt->ldpc_config(g_bm_m, g_bm_n, h_sc, h_st, h_wt);

    // packet allocation
    sim_pckt->ldpc_pckt_alloc();

    // LLR tables (reuse IBEX transceiver implementation)
    sim_pckt->ch_llr_alloc(MANUAL, rd_num, vref);
    sim_pckt->ch_llr_gen(m_llr0, m_llr1, finite_q_num - 1, 4);

    // decoder config
    sim_pckt->ldpc_dec_config(fdec_max_itr, 0, ldec_max_itr, m_alpha, 1, finite_q_num, finite_r_num, 4, sdlite_llr_config,
                              sdlite_llr0, sdlite_llr1, sdlite_llr2, sdlite_llr3);

    // IBEX BF decoder parameter passthrough (DV-side ldpc_decoder_inv style).
    // NOTE: input-related fields (e.g. corrupted codeword) will be set per-decode via `ldpc_ibex_input`.
    sim_pckt->ldpc_ibex_parameters(post_process_en, syndrome_weight_thr_qc, syndrome_weight_thr_post, early_termination_dis,
                                   ldpc_decoder_control_likelihood_0, ldpc_decoder_control_likelihood_1,
                                   ldpc_decoder_control_likelihood_2, ldpc_decoder_control_likelihood_3,
                                   ldpc_decoder_control_post, ldpc_early_term_0, ldpc_early_term_1, ldpc_early_term_2,
                                   ldpc_early_term_3, ldpc_early_term_4, ldpc_early_term_5, ldpc_early_term_6);

    g_configured = 1;
}

extern "C"
void ldpc_cleanup()
{
    if (!sim_pckt || !g_configured)
        return;

    sim_pckt->ldpc_clean();
    sim_pckt->ldpc_pckt_clean();
    sim_pckt->ch_llr_clean();

    g_configured = 0;
}

extern "C"
void ldpc_enc(svOpenArrayHandle usr_data_sv, svOpenArrayHandle enc_data_sv, int debug)
{
    if (!sim_pckt || !g_configured) {
        printf("[DVC IBEX ERROR] ldpc_enc called before ldpc_config.\n");
        return;
    }

    int *usr_data = (int *)svGetArrayPtr(usr_data_sv);
    int *enc_data = (int *)svGetArrayPtr(enc_data_sv);

    // 1) Unpack user bits (MSB-first) into usr_blk[0..info_len-1]
    unsigned int tmp = 0;
    int j = 0;
    for (int i = 0; i < sim_pckt->info_len; i++) {
        if (j == 0)
            tmp = (unsigned int)usr_data[i / 32];
        sim_pckt->usr_blk[i] = (char)((tmp >> (31 - j)) & 1u);
        j++;
        j = j % 32;
    }

    // 2) Encode (IBEX path depends on cir_sz)
    if (sim_pckt->cir_sz == 256) {
        if (sim_pckt->blk_len != (sim_pckt->info_len + sim_pckt->hm_m)) {
            printf("[DVC IBEX ERROR] sc=256 encoder requires blk_len==info_len+hm_m. blk_len=%d info_len=%d hm_m=%d\n",
                   sim_pckt->blk_len, sim_pckt->info_len, sim_pckt->hm_m);
            return;
        }
        sim_pckt->ldpc_encoder();
    } else if (sim_pckt->cir_sz == 512) {
        sim_pckt->ldpc_ibex_encoder();
    } else {
        printf("[DVC IBEX ERROR] U@nsupported cir_sz=%d in ldpc_enc.\n", sim_pckt->cir_sz);
        return;
    }

    // 3) Pack tx bits (MSB-first) from tx_blk[0..blk_len-1] into 32-bit words
    tmp = 0;
    j = 0;
    int k = 0;
    for (int i = 0; i < sim_pckt->blk_len; i++) {
        tmp = tmp * 2u + (unsigned int)(sim_pckt->tx_blk[i] & 1);
        j++;
        j = j % 32;
        if ((j == 0) || (i == (sim_pckt->blk_len - 1))) {
            enc_data[k] = (int)tmp;
            k++;
            tmp = 0;
        }
    }

    // MSB-left align for partial last word. With 1B-step padding, blk_len%32 may be 8/16/24.
    const int enc_rem = sim_pckt->blk_len % 32;
    if ((enc_rem != 0) && (k > 0)) {
        unsigned int last = (unsigned int)enc_data[k - 1];
        last = last << (32 - enc_rem);
        enc_data[k - 1] = (int)last;
    }
}

extern "C" void ldpc_dec(svOpenArrayHandle det_data_sv,
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
                         int h_sc)
{
    if (!sim_pckt || !g_configured) {
        printf("[DVC IBEX ERROR] ldpc_dec called before ldpc_config.\n");
        return;
    }

    int *dec_data = (int *)svGetArrayPtr(dec_data_sv);
    int *det_data = (int *)svGetArrayPtr(det_data_sv);
    const int n_words = (sim_pckt->blk_len + 31) / 32;
    int det_is_llr_bin = (sd_num >= 2);
    int det_hard_ones = 0;

    // 1) Unpack/copy detector input into det_blk[0..blk_len-1]
    if (sd_num < 2) {
        unsigned int tmp = 0;
        int j = 0;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            if (j == 0)
                tmp = (unsigned int)det_data[i / 32];
            sim_pckt->det_blk[i] = (char)((tmp >> (31 - j)) & 1u);
            det_hard_ones += (sim_pckt->det_blk[i] & 1);
            j++;
            j = j % 32;
        }
        dvc_check_words_against_bits_msb("ldpc_dec DET(hard)", det_data, sim_pckt->det_blk, sim_pckt->blk_len, debug);
    } else {
        const int max_bin = (sim_pckt->bin_num > 0) ? (sim_pckt->bin_num - 1) : 0;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            int v = det_data[i];
            if (v < 0)
                v = 0;
            if (v > max_bin)
                v = max_bin;
            sim_pckt->det_blk[i] = (char)v;
        }
        if (debug >= 2) {
            printf("[DVC IBEX DBG] ldpc_dec DET(soft) first bins:");
            const int n = (sim_pckt->blk_len > 16) ? 16 : sim_pckt->blk_len;
            for (int i = 0; i < n; i++)
                printf(" %d", (int)(unsigned char)sim_pckt->det_blk[i]);
            if (n < sim_pckt->blk_len)
                printf(" ...");
            printf("\n");
        }
    }

    // 2) Decode mode mapping (keep gen4 interface behavior, but allow IBEX-native values too)
    enum dec_model model = LAYER;
    if (dec_mode == 0)
        model = BF_P3;
    else if (dec_mode == 1)
        model = LAYER;
    else if (dec_mode == 2)
        model = (sd_num >= 2) ? BF_IBEX : BF_G2;
    else if (dec_mode == (int)BF_P3)
        model = BF_P3;
    else if (dec_mode == (int)BF_G2)
        model = BF_G2;
    else if (dec_mode == (int)BF_IBEX)
        model = BF_IBEX;
    else if (dec_mode == (int)LAYER)
        model = LAYER;
    else if (dec_mode == (int)LAYER_G2)
        model = LAYER_G2;
    else if (dec_mode == (int)SKIP)
        model = SKIP;
    else
        printf("[DVC IBEX WARN] unknown dec_mode=%d, fallback LAYER\n", dec_mode);

    if (debug >= 1) {
        printf("[DVC IBEX] ldpc_dec: dec_mode=%d -> model=%d sd_num=%d blk_len=%d n_words=%d bin_num=%d\n", dec_mode,
               (int)model, sd_num, sim_pckt->blk_len, n_words, sim_pckt->bin_num);
        printf("[DVC IBEX] ldpc_dec sizes: info_len=%d pad_len=%d hm_k=%d hm_m=%d hm_n=%d parity_det=%d\n",
               sim_pckt->info_len, sim_pckt->pad_len, sim_pckt->hm_k, sim_pckt->hm_m, sim_pckt->hm_n,
               sim_pckt->blk_len - sim_pckt->info_len);
    }

    // For LAYER/LAYER_G2 decoding, dec_di_blk is used as an index into llr_tbl (see ldpc_dec_layer/ldpc_dec_layer2).
    // If hard bits (0/1) are passed with sd_num<2, they must be mapped to "strong-0" / "strong-1" bins first.
    if (model == LAYER || model == LAYER_G2) {
        if (!sim_pckt->llr_tbl || sim_pckt->bin_num <= 0) {
            printf("[DVC IBEX ERROR] LAYER requires llr_tbl/bin_num configured.\n");
            return;
        }
        int bin_min = 0;
        int bin_max = 0;
        float llr_min = sim_pckt->llr_tbl[0];
        float llr_max = sim_pckt->llr_tbl[0];
        for (int b = 1; b < sim_pckt->bin_num; b++) {
            const float v = sim_pckt->llr_tbl[b];
            if (v < llr_min) {
                llr_min = v;
                bin_min = b;
            }
            if (v > llr_max) {
                llr_max = v;
                bin_max = b;
            }
        }

        if ((sd_num < 2) && !det_is_llr_bin) {
            // Map hard bits: 0 -> max-LLR bin (strong 0), 1 -> min-LLR bin (strong 1).
            int ones = 0;
            for (int i = 0; i < sim_pckt->blk_len; i++) {
                const char bit = sim_pckt->det_blk[i] & 1;
                ones += (bit != 0);
                sim_pckt->det_blk[i] = (bit == 0) ? (char)bin_max : (char)bin_min;
            }
            det_is_llr_bin = 1;
            if (debug >= 1) {
                printf("[DVC IBEX] ldpc_dec LAYER hard->bin: ones_in=%d/%d bit0->bin%d(llr=%f) bit1->bin%d(llr=%f)\n",
                       ones, sim_pckt->blk_len, bin_max, llr_max, bin_min, llr_min);
                if (sim_pckt->max_llr_bin != bin_max) {
                    printf("[DVC IBEX WARN] max_llr_bin mismatch: cfg=%d argmax=%d (LLR cfg=%f argmax=%f)\n",
                           sim_pckt->max_llr_bin, bin_max, sim_pckt->llr_tbl[sim_pckt->max_llr_bin], llr_max);
                }
                printf("[DVC IBEX] ldpc_dec LAYER cfg: ldec_max_itr=%d early_term=%d alpha=%f finite(mode=%d q=%d r=%d f=%d)\n",
                       sim_pckt->ldec_max_itr, sim_pckt->ldec_early_term_en, sim_pckt->alpha, sim_pckt->finite_mode,
                       sim_pckt->finite_q_num, sim_pckt->finite_r_num, sim_pckt->finite_f_num);
            }
        }

        if (debug >= 1) {
            int neg = 0;
            int bin0_cnt = 0;
            int bin1_cnt = 0;
            int bin_other_cnt = 0;
            for (int i = 0; i < sim_pckt->blk_len; i++) {
                const int bin = (int)(unsigned char)sim_pckt->det_blk[i];
                if (bin == 0)
                    bin0_cnt++;
                else if (bin == 1)
                    bin1_cnt++;
                else
                    bin_other_cnt++;
                if ((bin >= 0) && (bin < sim_pckt->bin_num) && (sim_pckt->llr_tbl[bin] < 0.0f))
                    neg++;
            }
            if (sim_pckt->bin_num == 2) {
                printf("[DVC IBEX] ldpc_dec LAYER input LLR neg=%d/%d (bin_num=2 bins0=%d bins1=%d other=%d llr0=%f llr1=%f)\n",
                       neg, sim_pckt->blk_len, bin0_cnt, bin1_cnt, bin_other_cnt,
                       sim_pckt->llr_tbl[0], sim_pckt->llr_tbl[1]);
            } else {
                printf("[DVC IBEX] ldpc_dec LAYER input LLR neg=%d/%d (bin_num=%d max_llr_bin=%d llr[max]=%f llr[min]=%f)\n",
                       neg, sim_pckt->blk_len, sim_pckt->bin_num, sim_pckt->max_llr_bin,
                       sim_pckt->llr_tbl[sim_pckt->max_llr_bin], llr_min);
            }

            if (sim_pckt->tx_blk && sim_pckt->bin_num == 2) {
                int tx_ones = 0;
                int rx_ones = 0;
                int rx_vs_tx = 0;
                for (int i = 0; i < sim_pckt->blk_len; i++) {
                    const int txb = sim_pckt->tx_blk[i] & 1;
                    const int bin = (int)(unsigned char)sim_pckt->det_blk[i];
                    const int rxb = (sim_pckt->llr_tbl[bin] < 0.0f) ? 1 : 0;
                    tx_ones += txb;
                    rx_ones += rxb;
                    rx_vs_tx += (txb != rxb);
                }
                printf("[DVC IBEX] ldpc_dec LLR check: tx_ones=%d rx_ones=%d rx_vs_tx_mismatch=%d\n", tx_ones, rx_ones,
                       rx_vs_tx);
            }
        }
    }

    // If a hard-decision decoder is selected while sd_num>=2, derive hard bits from LLR table sign.
    if ((sd_num >= 2) && (model != LAYER) && (model != BF_IBEX) && (model != SKIP)) {
        printf("[DVC IBEX WARN] sd_num=%d with hard decoder (mode=%d): converting bins to hard bits via llr_tbl sign.\n",
               sd_num, (int)model);
        const int max_bin = (sim_pckt->bin_num > 0) ? (sim_pckt->bin_num - 1) : 0;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            int bin = (int)(unsigned char)sim_pckt->det_blk[i];
            if (bin > max_bin)
                bin = max_bin;
            sim_pckt->det_blk[i] = (sim_pckt->llr_tbl && sim_pckt->llr_tbl[bin] >= 0.0f) ? 0 : 1;
        }
    }

    // BF_IBEX requires its dedicated input structure (DVsrc ldpc_decoder_inv style).
    if (model == BF_IBEX) {
        if (sd_num < 2) {
            printf("[DVC IBEX ERROR] BF_IBEX requires sd_num>=2 (soft decision). sd_num=%d\n", sd_num);
            return;
        }
        // NOTE: `ldpc_ibex_input()` derives hard bits from `rx_blk` sign, not from `det_blk`.
        // When DV calls `ldpc_dec()` directly (bypassing `sd_err_inj/ch_err_inj`), `rx_blk` may be stale.
        // Reconstruct a consistent `rx_blk` from detector bins and the configured LLR table so that:
        //   bit_hard = (llr < 0) ? 1 : 0
        if (sim_pckt->rx_blk && sim_pckt->llr_tbl && sim_pckt->bin_num > 0) {
            const int max_bin = sim_pckt->bin_num - 1;
            for (int i = 0; i < sim_pckt->blk_len; i++) {
                int bin = (int)(unsigned char)sim_pckt->det_blk[i];
                if (bin > max_bin)
                    bin = max_bin;
                sim_pckt->rx_blk[i] = sim_pckt->llr_tbl[bin];
            }
        } else if (debug >= 1) {
            printf("[DVC IBEX WARN] BF_IBEX rx_blk/llr_tbl unavailable: rx_blk=%p llr_tbl=%p bin_num=%d\n",
                   sim_pckt->rx_blk, sim_pckt->llr_tbl, sim_pckt->bin_num);
        }
        if ((g_nand_strobes > 0) && (g_nand_strobes != sd_num)) {
            printf("[DVC IBEX WARN] nand_strobes mismatch: cfg=%d call_sd_num=%d (soft_bits uses cfg)\n", g_nand_strobes,
                   sd_num);
        }
        sim_pckt->ldpc_ibex_input(/*syndrome_cal_only=*/0, g_max_iter, g_post_iter, g_nand_strobes);
    }

    // 3) Run decode (IBEX codec handles [Info|Pad|Parity] reconstruction internally)
    sim_pckt->ldpc_decoder(model);

    if ((debug >= 1) && (model == LAYER || model == LAYER_G2)) {
        int c0 = 0;
        int c1 = 0;
        int cother = 0;
        int first_i = -1;
        int first_v = 0;
        for (int bi = 0; bi < sim_pckt->blk_len; bi++) {
            const int v = (int)(unsigned char)sim_pckt->dec_blk[bi];
            if (v == 0)
                c0++;
            else if (v == 1)
                c1++;
            else {
                if (first_i < 0) {
                    first_i = bi;
                    first_v = v;
                }
                cother++;
            }
        }
        printf("[DVC IBEX] ldpc_dec LAYER bits: in_ones=%d out_ones=%d (blk_len=%d)\n",
               (sd_num < 2) ? det_hard_ones : -1, c1, sim_pckt->blk_len);
        printf("[DVC IBEX] ldpc_dec LAYER dec_blk distribution: 0=%d 1=%d other=%d\n", c0, c1, cother);
        if (cother > 0)
            printf("[DVC IBEX WARN] ldpc_dec LAYER dec_blk has non-binary value at i=%d v=%d\n", first_i, first_v);
    }

    // 4) Status outputs
    *dec_unc_sv = sim_pckt->cw_fail;
    *init_synd_wt_sv = sim_pckt->init_synd_wt;
    *dec_cnvg_itr_sv = sim_pckt->cnvg_itr;
    *dec_cnvg_col_sv = sim_pckt->cnvg_lyr;
    *fina_synd_wt_sv = sim_pckt->fina_synd_wt;

    if (debug >= 1) {
        printf("[DVC IBEX] ldpc_dec status: cw_fail=%d init_synd=%d fina_synd=%d cnvg_itr=%d cnvg_col=%d\n",
               sim_pckt->cw_fail, sim_pckt->init_synd_wt, sim_pckt->fina_synd_wt, sim_pckt->cnvg_itr, sim_pckt->cnvg_lyr);
    }

    if ((debug >= 1) && sim_pckt->tx_blk) {
        int mis_total = 0;
        int mis_info = 0;
        int mis_parity = 0;
        for (int bi = 0; bi < sim_pckt->blk_len; bi++) {
            const int txb = sim_pckt->tx_blk[bi] & 1;
            const int dcb = sim_pckt->dec_blk[bi] & 1;
            if (txb != dcb) {
                mis_total++;
                if (bi < sim_pckt->info_len)
                    mis_info++;
                else
                    mis_parity++;
            }
        }
        printf("[DVC IBEX] ldpc_dec tx vs dec mismatch: total=%d info=%d parity=%d\n", mis_total, mis_info, mis_parity);

        // Minimal parity "shift probe": check if decoded parity matches TX parity better after a small fixed shift.
        if (debug >= 2) {
            const int parity_start = sim_pckt->info_len;
            const int parity_len = sim_pckt->blk_len - parity_start;
            if (parity_len > 0) {
                const int shifts[4] = {-512, -4, 4, 512};
                int best_shift = 0;
                int best_m = mis_parity;
                int best_n = parity_len;
                for (int si = 0; si < 4; si++) {
                    const int sh = shifts[si];
                    int m = 0;
                    int n = 0;
                    for (int i = 0; i < parity_len; i++) {
                        const int j = i + sh;
                        if ((j < 0) || (j >= parity_len))
                            continue;
                        const int a = sim_pckt->tx_blk[parity_start + j] & 1;
                        const int b = sim_pckt->dec_blk[parity_start + i] & 1;
                        n++;
                        m += (a != b);
                    }
                    if ((n > 0) && ((long long)m * (long long)best_n < (long long)best_m * (long long)n)) {
                        best_shift = sh;
                        best_m = m;
                        best_n = n;
                    }
                }
                printf("[DVC IBEX] parity_shift_probe: base=%d/%d best_shift=%d best=%d/%d\n", mis_parity, parity_len,
                       best_shift, best_m, best_n);
            }
        }

        // (1) Parity column signatures (TX/RX/DEC) for quick shift detection.
        if ((debug >= 2) && (model == LAYER_G2)) {
            const int z = sim_pckt->cir_sz;
            const int parity_start = sim_pckt->info_len;
            const int parity_len = sim_pckt->blk_len - parity_start;
            const int cols = (z > 0) ? (parity_len / z) : 0;
            if ((cols > 0) && sim_pckt->dec_blk) {
                printf("[DVC IBEX] parity_col_sig z=%d cols=%d tx=", z, cols);
                for (int c = 0; c < cols; c++) {
                    const int off = parity_start + c * z;
                    const unsigned int s = dvc_sig_fnv1a_bits01(sim_pckt->tx_blk, off, z);
                    printf("%s%08x", (c == 0) ? "" : ",", s);
                }
                printf(" rx=");
                for (int c = 0; c < cols; c++) {
                    const int off = parity_start + c * z;
                    const unsigned int s = dvc_sig_fnv1a_rx_bins01(sim_pckt, off, z);
                    printf("%s%08x", (c == 0) ? "" : ",", s);
                }
                printf(" dec=");
                for (int c = 0; c < cols; c++) {
                    const int off = parity_start + c * z;
                    const unsigned int s = dvc_sig_fnv1a_bits01(sim_pckt->dec_blk, off, z);
                    printf("%s%08x", (c == 0) ? "" : ",", s);
                }
                printf("\n");
            }
        }

        // (2) Syndrome-weight after rotating parity by ±Z inside full QC view (dec_do_blk).
        if ((debug >= 2) && (model == LAYER_G2) && sim_pckt->qc_bm && sim_pckt->dec_do_blk && (sim_pckt->cir_sz > 0) &&
            (sim_pckt->hm_m > 0) && (sim_pckt->hm_n > 0) && (sim_pckt->hm_k >= 0) &&
            ((sim_pckt->hm_k + sim_pckt->hm_m) <= sim_pckt->hm_n)) {
            const int z = sim_pckt->cir_sz;
            const int synd0 = dvc_qc_syndrome_weight_bm(sim_pckt->qc_bm, sim_pckt->dec_do_blk, sim_pckt->bm_m, z);
            int synd_m1 = -1;
            int synd_p1 = -1;
            char *tmp = (char *)calloc(sim_pckt->hm_n, sizeof(*tmp));
            if (tmp) {
                memcpy(tmp, sim_pckt->dec_do_blk, sim_pckt->hm_n);
                dvc_rotate_segment_left(tmp, sim_pckt->hm_k, sim_pckt->hm_m, z);
                synd_m1 = dvc_qc_syndrome_weight_bm(sim_pckt->qc_bm, tmp, sim_pckt->bm_m, z);

                memcpy(tmp, sim_pckt->dec_do_blk, sim_pckt->hm_n);
                dvc_rotate_segment_left(tmp, sim_pckt->hm_k, sim_pckt->hm_m, sim_pckt->hm_m - z);
                synd_p1 = dvc_qc_syndrome_weight_bm(sim_pckt->qc_bm, tmp, sim_pckt->bm_m, z);
                free(tmp);
            }
            printf("[DVC IBEX] dec_synd_rot z=%d synd=%d rot(-1)=%d rot(+1)=%d\n", z, synd0, synd_m1, synd_p1);
        }
    }

    if ((sd_num < 2) && !det_is_llr_bin && (debug >= 1)) {
        dvc_report_bit_mismatches_msb("[DVC IBEX] ldpc_dec det vs dec", sim_pckt->det_blk, sim_pckt->dec_blk,
                                      sim_pckt->blk_len, (debug >= 2) ? 4 : 0);
    }

    // 5) Pack decoded bits (MSB-first) from dec_blk[0..blk_len-1]
    unsigned int tmp = 0;
    int j = 0;
    int k = 0;
    int i;
    for (i = 0; i < sim_pckt->blk_len; i++) {
        tmp = tmp * 2u + (unsigned int)(sim_pckt->dec_blk[i] & 1);
        j++;
        j = j % 32;
        if ((j == 0) || (i == (sim_pckt->blk_len - 1))) {
            dec_data[k] = (int)tmp;
            tmp = 0;
            k++;
        }
    }

    // MSB-left align for partial last word. With 1B-step padding, blk_len%32 may be 8/16/24.
    const int dec_rem = sim_pckt->blk_len % 32;
    if ((dec_rem != 0) && (k > 0)) {
        unsigned int last = (unsigned int)dec_data[k - 1];
        last = last << (32 - dec_rem);
        dec_data[k - 1] = (int)last;
    }

    if ((sd_num < 2) && (debug >= 1)) {
        const unsigned int last_mask = dvc_mask_valid_bits_msb(sim_pckt->blk_len);
        int same_words = 1;
        for (int w = 0; w < n_words; w++) {
            const unsigned int mask = (w == (n_words - 1)) ? last_mask : 0xFFFFFFFFu;
            const unsigned int a = ((unsigned int)det_data[w]) & mask;
            const unsigned int b = ((unsigned int)dec_data[w]) & mask;
            if (a != b) {
                same_words = 0;
                break;
            }
        }
        printf("[DVC IBEX] ldpc_dec hard in/out %s\n", same_words ? "SAME" : "DIFF");
    }

    // Pad remaining words up to h_n*h_sc bits with zeros (gen4 behavior)
    while (i % 32 != 0)
        i++;
    for (; i < h_n * h_sc; i += 32) {
        dec_data[k] = 0;
        k++;
    }
}

extern "C"
void ch_update(int ch_mode, int ch_para, int debug)
{
    if (!sim_pckt || !g_configured) {
        printf("[DVC IBEX ERROR] ch_update called before ldpc_config.\n");
        return;
    }

    const float m_ch_para = (float)ch_para / 100.0f;

    enum ch_model ch_sel = CLEAN;
    if (ch_mode == 0)
        ch_sel = CLEAN;
    else if (ch_mode == 1)
        ch_sel = AWGN;
    else if (ch_mode == 2)
        ch_sel = BSC;
    else if (ch_mode == 3)
        ch_sel = ERR_INJ;
    else if (ch_mode == 4)
        ch_sel = MAX_ERR;
    else
        printf("[DVC IBEX WARN] unknown ch_mode=%d, fallback CLEAN\n", ch_mode);

    sim_pckt->ch_config(sim_pckt->info_len, sim_pckt->blk_len, ch_sel, m_ch_para, g_vn_bits);
    if (debug >= 1) {
        if (ch_sel == BSC)
            printf("[DVC IBEX] ch_update: mode=%d(BSC) ch_para=%d => ber=%f\n", ch_mode, ch_para, sim_pckt->ber);
        else if ((ch_sel == ERR_INJ) || (ch_sel == MAX_ERR))
            printf("[DVC IBEX] ch_update: mode=%d(ERR_INJ/MAX_ERR) ch_para=%d => err_num=%d\n", ch_mode, ch_para,
                   sim_pckt->err_num);
        else if ((ch_sel == AWGN) || (ch_sel == ALL_ZERO))
            printf("[DVC IBEX] ch_update: mode=%d(AWGN/ALL_ZERO) ch_para=%d => snr=%f\n", ch_mode, ch_para, sim_pckt->snr);
        else
            printf("[DVC IBEX] ch_update: mode=%d(CLEAN) ch_para=%d\n", ch_mode, ch_para);
    }
}

extern "C"
void ch_err_inj(svOpenArrayHandle tx_data_sv, svOpenArrayHandle rx_data_sv, int debug)
{
    if (!sim_pckt || !g_configured) {
        printf("[DVC IBEX ERROR] ch_err_inj called before ldpc_config.\n");
        return;
    }

    int *tx_data = (int *)svGetArrayPtr(tx_data_sv);
    int *det_data = (int *)svGetArrayPtr(rx_data_sv);

    // 1) Unpack TX (MSB-first) into tx_blk[0..blk_len-1]
    unsigned int tmp = 0;
    int j = 0;
    for (int i = 0; i < sim_pckt->blk_len; i++) {
        if (j == 0)
            tmp = (unsigned int)tx_data[i / 32];
        sim_pckt->tx_blk[i] = (char)((tmp >> (31 - j)) & 1u);
        j++;
        j = j % 32;
    }

    dvc_check_words_against_bits_msb("ch_err_inj TX", tx_data, sim_pckt->tx_blk, sim_pckt->blk_len, debug);

    // Sanity: the transmitted codeword should satisfy H (syndrome weight == 0).
    // Use IBEX's own H representation to avoid relying on qc_hm (may be uninitialized on sc=512 path).
    if ((debug >= 1) && (sim_pckt->hm_n > 0) && (sim_pckt->hm_m > 0) && (sim_pckt->hm_k > 0) &&
        (sim_pckt->info_len >= 0) && (sim_pckt->pad_len >= 0) && (sim_pckt->blk_len > 0) && (sim_pckt->cir_sz == 512)) {
        char *cw_full = (char *)calloc(sim_pckt->hm_n, sizeof(*cw_full));
        if (!cw_full) {
            printf("[DVC IBEX WARN] ch_err_inj tx_synd_wt_ibex skipped: alloc failed (hm_n=%d)\n", sim_pckt->hm_n);
        } else {
            // blk layout: [info_len][parity(bits_in_det)]
            // hm  layout: [info_len][pad_len zeros][parity(hm_m full)]
            for (int i = 0; i < sim_pckt->info_len && i < sim_pckt->blk_len; i++)
                cw_full[i] = (char)(sim_pckt->tx_blk[i] & 1);

            // Insert parity. For non-shortened case parity_bits_in_det == hm_m, this is a simple copy.
            // For shortened parity tail in the first parity column, fill missing tail bits with 0.
            const int parity_bits_in_det = sim_pckt->blk_len - sim_pckt->info_len;
            if ((sim_pckt->h_matrix.extra_bits_of_parity > 0) && (parity_bits_in_det < sim_pckt->hm_m)) {
                int first_parity_bits = sim_pckt->h_matrix.extra_bits_of_parity;
                if (first_parity_bits > parity_bits_in_det)
                    first_parity_bits = parity_bits_in_det;
                for (int i = 0; i < first_parity_bits; i++)
                    cw_full[sim_pckt->hm_k + i] = (char)(sim_pckt->tx_blk[sim_pckt->info_len + i] & 1);
                for (int i = first_parity_bits; i < sim_pckt->cir_sz; i++)
                    cw_full[sim_pckt->hm_k + i] = 0;
                const int remaining = parity_bits_in_det - first_parity_bits;
                for (int i = 0; i < remaining; i++)
                    cw_full[sim_pckt->hm_k + sim_pckt->cir_sz + i] =
                        (char)(sim_pckt->tx_blk[sim_pckt->info_len + first_parity_bits + i] & 1);
            } else {
                for (int i = 0; i < parity_bits_in_det; i++) {
                    const int src = sim_pckt->info_len + i;
                    const int dst = sim_pckt->hm_k + i;
                    if ((src >= 0) && (src < sim_pckt->blk_len) && (dst >= 0) && (dst < sim_pckt->hm_n))
                        cw_full[dst] = (char)(sim_pckt->tx_blk[src] & 1);
                }
            }

            // Pad zeros (shortening) already zero-initialized in cw_full[info_len..hm_k-1].
            s_hard_codeword vn;
            memset(&vn, 0, sizeof(vn));
            for (int j = 0; j < sim_pckt->h_matrix.cols; j++) {
                for (int k = 0; k < sim_pckt->h_matrix.bits; k++) {
                    vn.c[j].b[k] = (cw_full[j * sim_pckt->h_matrix.bits + k] & 1) ? true : false;
                }
            }
            s_check_nodes cn = sim_pckt->f_check_nodes(sim_pckt->h_matrix, vn);
            const int wt_ibex = sim_pckt->f_check_node_weight(sim_pckt->h_matrix, cn);
            const int wt_qc =
                (sim_pckt->qc_bm != NULL) ? dvc_qc_syndrome_weight_bm(sim_pckt->qc_bm, cw_full, sim_pckt->bm_m, sim_pckt->cir_sz) : -1;
            printf("[DVC IBEX] ch_err_inj tx_synd_wt: ibex=%d qc_bm=%d (expect 0)\n", wt_ibex, wt_qc);
        }
        if (cw_full)
            free(cw_full);
    }

    // 2) Channel transmit
    sim_pckt->ch_transmit();
    if (debug >= 1) {
        if (sim_pckt->ch_sel == BSC)
            printf("[DVC IBEX] ch_err_inj: ch_sel=BSC ber=%f\n", sim_pckt->ber);
        else if ((sim_pckt->ch_sel == ERR_INJ) || (sim_pckt->ch_sel == MAX_ERR))
            printf("[DVC IBEX] ch_err_inj: ch_sel=%s err_num=%d raw_err_num=%d\n",
                   (sim_pckt->ch_sel == ERR_INJ) ? "ERR_INJ" : "MAX_ERR", sim_pckt->err_num, sim_pckt->raw_err_num);
        else if ((sim_pckt->ch_sel == AWGN) || (sim_pckt->ch_sel == ALL_ZERO))
            printf("[DVC IBEX] ch_err_inj: ch_sel=%s snr=%f\n", (sim_pckt->ch_sel == AWGN) ? "AWGN" : "ALL_ZERO",
                   sim_pckt->snr);
        else
            printf("[DVC IBEX] ch_err_inj: ch_sel=CLEAN\n");
    }

    // 3) Hard decision (threshold 0), consistent with DVCtrans hard detector behavior
    for (int i = 0; i < sim_pckt->blk_len; i++) {
        sim_pckt->det_blk[i] = (sim_pckt->rx_blk[i] >= 0.0f) ? 0 : 1;
    }

    if (debug >= 1) {
        dvc_report_bit_mismatches_msb("[DVC IBEX] ch_err_inj tx vs det", sim_pckt->tx_blk, sim_pckt->det_blk,
                                      sim_pckt->blk_len, (debug >= 2) ? 4 : 0);
    }

    // 4) Output: keep original TX words for bits beyond blk_len (avoid tail-word mismatch when blk_len%32!=0)
    const int n_words = (sim_pckt->blk_len + 31) / 32;
    for (int w = 0; w < n_words; w++) {
        unsigned int word = (unsigned int)tx_data[w];
        for (int b = 31; b >= 0; b--) {
            const int bit_index = w * 32 + (31 - b);
            if (bit_index >= sim_pckt->blk_len)
                break;

            const unsigned int bit = (unsigned int)(sim_pckt->det_blk[bit_index] & 1u);
            const unsigned int mask = 1u << b;
            if (bit)
                word |= mask;
            else
                word &= ~mask;
        }
        det_data[w] = (int)word;
    }
}

extern "C" void sd_err_inj(svOpenArrayHandle tx_data_sv,
                           svOpenArrayHandle rx_data_sv,
                           svOpenArrayHandle rd_data_sv,
                           svOpenArrayHandle llr_tbl_sv,
                           svOpenArrayHandle split_bin_sv,
                           svOpenArrayHandle bin_id_sv,
                           int sd_num,
                           svOpenArrayHandle v_ref_sv,
                           int debug)
{
    if (!sim_pckt || !g_configured) {
        printf("[DVC IBEX ERROR] sd_err_inj called before ldpc_config.\n");
        return;
    }
    if (sd_num < 2) {
        return;
    }

    int rd_num = sd_num;
    if (rd_num > 127)
        rd_num = 127;

    int *tx_data = (int *)svGetArrayPtr(tx_data_sv);
    int *bin_data = (int *)svGetArrayPtr(rx_data_sv);
    int *rd_data = (int *)svGetArrayPtr(rd_data_sv);
    int *llr_tbl = (int *)svGetArrayPtr(llr_tbl_sv);
    int *split_bin = (int *)svGetArrayPtr(split_bin_sv);
    int *bin_id = (int *)svGetArrayPtr(bin_id_sv);
    int *v_ref = (int *)svGetArrayPtr(v_ref_sv);

    float vref[128] = {0.0f};
    for (int i = 0; i < rd_num; i++)
        vref[i] = (float)v_ref[i] / 1000.0f;

    if (sim_pckt->rd_num != rd_num) {
        printf("[DVC IBEX WARN] sd_err_inj sd_num mismatch: cfg rd_num=%d, call rd_num=%d\n", sim_pckt->rd_num, rd_num);
    }
    if (sim_pckt->vref && (sim_pckt->rd_num >= rd_num)) {
        for (int i = 0; i < rd_num; i++) {
            const float diff = (float)fabs((double)(sim_pckt->vref[i] - vref[i]));
            if (diff > 1e-6f) {
                printf("[DVC IBEX WARN] sd_err_inj vref mismatch at i=%d: cfg=%f call=%f\n", i, sim_pckt->vref[i], vref[i]);
                break;
            }
        }
    }

    // 1) Unpack TX (MSB-first) into tx_blk[0..blk_len-1]
    unsigned int tmp = 0;
    int j = 0;
    for (int i = 0; i < sim_pckt->blk_len; i++) {
        if (j == 0)
            tmp = (unsigned int)tx_data[i / 32];
        sim_pckt->tx_blk[i] = (char)((tmp >> (31 - j)) & 1u);
        j++;
        j = j % 32;
    }

    // 2) Channel transmit + bin detector (bin index output)
    sim_pckt->ch_transmit();
    if (debug >= 1) {
        if (sim_pckt->ch_sel == BSC)
            printf("[DVC IBEX] sd_err_inj: ch_sel=BSC ber=%f\n", sim_pckt->ber);
        else if ((sim_pckt->ch_sel == ERR_INJ) || (sim_pckt->ch_sel == MAX_ERR))
            printf("[DVC IBEX] sd_err_inj: ch_sel=%s err_num=%d raw_err_num=%d\n",
                   (sim_pckt->ch_sel == ERR_INJ) ? "ERR_INJ" : "MAX_ERR", sim_pckt->err_num, sim_pckt->raw_err_num);
        else if ((sim_pckt->ch_sel == AWGN) || (sim_pckt->ch_sel == ALL_ZERO))
            printf("[DVC IBEX] sd_err_inj: ch_sel=%s snr=%f\n", (sim_pckt->ch_sel == AWGN) ? "AWGN" : "ALL_ZERO",
                   sim_pckt->snr);
        else
            printf("[DVC IBEX] sd_err_inj: ch_sel=CLEAN\n");
    }

    // 3) Gen4/DVC style soft detector:
    //    - det_blk init to 0
    //    - for each read rd_indx: rd_bit = (rx>=vref)?0:1
    //      if det_blk == split_bin[rd_indx] and rd_bit==1 => det_blk = rd_indx+1
    for (int i = 0; i < sim_pckt->blk_len; i++)
        sim_pckt->det_blk[i] = 0;

    // 4) Read data out (packed bits), rd_data order: rd_indx major, then bit index
    tmp = 0;
    j = 0;
    int k = 0;
    for (int rd_indx = 0; rd_indx < rd_num; rd_indx++) {
        const int split = (sim_pckt->bin_split && (rd_indx < sim_pckt->rd_num)) ? (int)(unsigned char)sim_pckt->bin_split[rd_indx]
                                                                                : 0;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            const unsigned int bit = (sim_pckt->rx_blk[i] >= vref[rd_indx]) ? 0u : 1u;
            if (((int)(unsigned char)sim_pckt->det_blk[i] == split) && (bit == 1u))
                sim_pckt->det_blk[i] = (char)(rd_indx + 1);

            tmp = tmp * 2u + bit;
            j++;
            j = j % 32;
            if ((j == 0) || ((rd_indx == (rd_num - 1)) && (i == (sim_pckt->blk_len - 1)))) {
                rd_data[k] = (int)tmp;
                tmp = 0;
                k++;
            }
        }
    }

    // 5) Bin out (per-bit integer)
    for (int i = 0; i < sim_pckt->blk_len; i++) {
        bin_data[i] = (int)(unsigned char)sim_pckt->det_blk[i];
    }

    // 6) LLR table (scaled by 16 as fixed-point with frac=4), fill up to 128 elements
    for (int i = 0; i < 128; i++) {
        if (sim_pckt->llr_tbl && (i < sim_pckt->bin_num))
            llr_tbl[i] = (int)(sim_pckt->llr_tbl[i] * 16.0f);
        else
            llr_tbl[i] = 0;
    }

    // 7) split_bin / bin_id outputs (fill rest with 0)
    for (int i = 0; i < 128; i++) {
        if (sim_pckt->bin_split && (i < sim_pckt->rd_num))
            split_bin[i] = (int)(unsigned char)sim_pckt->bin_split[i];
        else
            split_bin[i] = 0;
    }
    for (int i = 0; i < 128; i++) {
        if (sim_pckt->bin_asc_ord && (i <= sim_pckt->rd_num))
            bin_id[i] = (int)(unsigned char)sim_pckt->bin_asc_ord[i];
        else
            bin_id[i] = 0;
    }
}
