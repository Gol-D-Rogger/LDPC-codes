#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Forward-declare vpi_printf so we don't have to depend on vpi_user.h being on
// the IDE/clangd include path; the simulator provides the real symbol at link
// time (same approach as DVCtrans/src/ldpc_c_model.c).
extern "C" int vpi_printf(const char *format, ...);

#include "rand.h"
#include "alloc.h"
#include "vec_op.h"
#include "mod2sparse.h"
#include "transceiver.h"
#include "ldpc_codec.h"

static ldpc_packet ldpc_pckt;
struct ldpc_packet *sim_pckt;

static int g_configured = 0;

static int g_vn_bits = 3;

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
static int g_sd_err_inj_dump_seq = 0;

static const int kFiniteFracBits = 3;
static const float kFiniteFracScale = 8.0f;

extern "C" void ldpc_cleanup();

static unsigned int dvc_mask_valid_bits_msb(int total_bits)
{
    const int rem = total_bits % 32;
    if (rem == 0)
        return 0xFFFFFFFFu;
    return 0xFFFFFFFFu << (32 - rem);
}

static unsigned int dvc_bitreverse32(unsigned int x)
{
    x = ((x & 0x55555555u) << 1) | ((x >> 1) & 0x55555555u);
    x = ((x & 0x33333333u) << 2) | ((x >> 2) & 0x33333333u);
    x = ((x & 0x0F0F0F0Fu) << 4) | ((x >> 4) & 0x0F0F0F0Fu);
    x = ((x & 0x00FF00FFu) << 8) | ((x >> 8) & 0x00FF00FFu);
    x = (x << 16) | (x >> 16);
    return x;
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

static int dvc_parity_stream_to_hm_index(const ldpc_packet *pckt, int parity_idx)
{
    if (!pckt || (parity_idx < 0))
        return -1;

    if ((pckt->h_matrix.unused_bytes_of_parity > 0) && (pckt->cir_sz > 0)) {
        const int first_valid_bits = pckt->h_matrix.extra_bits_of_parity;
        if (parity_idx < first_valid_bits)
            return pckt->hm_k + parity_idx;
        return pckt->hm_k + pckt->cir_sz + (parity_idx - first_valid_bits);
    }

    return pckt->hm_k + parity_idx;
}

static int dvc_fill_sd_bin_data_max_llr_view(const ldpc_packet *pckt, int *bin_data, int bin_elems)
{
    if (!pckt || !bin_data || (bin_elems <= 0))
        return 0;

    for (int i = 0; i < bin_elems; i++)
        bin_data[i] = 0;

    const int use_hm_layout =
        (pckt->hm_n > pckt->blk_len) &&
        (bin_elems >= pckt->hm_n) &&
        (pckt->hm_k >= pckt->info_len) &&
        ((pckt->hm_k + pckt->hm_m) <= pckt->hm_n);
    if (!use_hm_layout) {
        for (int i = 0; i < pckt->blk_len && i < bin_elems; i++)
            bin_data[i] = (int)(unsigned char)pckt->det_blk[i];
        return 0;
    }

    // DV soft-buffer export view:
    //   [User Data][User Padding=max][Valid Parity][PMSK]
    //   0         info_len          hm_k            hm_k+valid_parity hm_n
    int bit_wd = 0;
    if (rd_num >= 2 && rd_num < 4)
        bit_wd = 2;
    else if (rd_num >= 4 && rd_num < 8)
        bit_wd = 3;
    else if (rd_num >= 8 && rd_num < 16)
        bit_wd = 4;
    else if (rd_num >= 16 && rd_num < 32)
        bit_wd = 5;
    else if (rd_num >= 32 && rd_num < 64)
        bit_wd = 6;
    else if (rd_num >= 64 && rd_num < 128)
        bit_wd = 7;

    int parity_mask_bin = pckt->max_llr_bin;
    for (int i = 0; i < pckt->hm_n; i++)
    {
        if (sim_pckt->sd_type == MANUAL)
            parity_mask_bin = sim_pckt->max_llr_bin;
        else if (sim_pckt->sd_type == VENDOR1)
            parity_mask_bin = int(pow(2, bit_wd) - 2);
        else if (sim_pckt->sd_type == VENDOR0)
            parity_mask_bin = 0;
    }

    for (int i = 0; i < pckt->hm_n; i++)
        bin_data[i] = pckt->max_llr_bin;

    for (int i = 0; i < pckt->info_len && i < pckt->blk_len; i++)
        bin_data[i] = (int)(unsigned char)pckt->det_blk[i];

    const int parity_bits_in_det = pckt->blk_len - pckt->info_len;
    for (int i = 0; i < parity_bits_in_det; i++) {
        const int src = pckt->info_len + i;
        const int dst = pckt->hm_k + i;
        if ((src >= 0) && (src < pckt->blk_len) && (dst >= 0) && (dst < pckt->hm_n))
            bin_data[dst] = (int)(unsigned char)pckt->det_blk[src];
    }

    const int pmsk_start = pckt->hm_k + parity_bits_in_det;
    for (int i = pmsk_start; i < pckt->hm_n; i++) {
        if ((i >= 0) && (i < bin_elems))
            bin_data[i] = parity_mask_bin;
    }

    return 1;
}

static void dvc_dump_codec_overview(const ldpc_packet *pckt,
                                    int dec_mode_in,
                                    enum dec_model model,
                                    int sd_num,
                                    int debug)
{
    if (!pckt || debug < 1)
        return;

    const char *ch_name = "UNKNOWN";
    switch (pckt->ch_sel) {
    case CLEAN:
        ch_name = "CLEAN";
        break;
    case AWGN:
        ch_name = "AWGN";
        break;
    case BSC:
        ch_name = "BSC";
        break;
    case ERR_INJ:
        ch_name = "ERR_INJ";
        break;
    default:
        break;
    }

    const char *model_name = "UNKNOWN";
    switch (model) {
    case SKIP:
        model_name = "SKIP";
        break;
    case BF_P0:
        model_name = "BF_P0";
        break;
    case BF_P3:
        model_name = "BF_P3";
        break;
    case BF_G2:
        model_name = "BF_G2";
        break;
    case BF_IBEX:
        model_name = "BF_IBEX";
        break;
    case BF_IBEX_RTL_CN:
        model_name = "BF_IBEX_RTL_CN";
        break;
    case LAYER:
        model_name = "LAYER";
        break;
    case LAYER_G2:
        model_name = "LAYER_G2";
        break;
    default:
        break;
    }

    int det_bin_hist[128] = {0};
    int det_ones = 0;
    int tx_ones = 0;
    int dec_ones = 0;
    int tx_vs_det_total = 0;
    int tx_vs_det_info = 0;
    int tx_vs_det_parity = 0;
    int tx_vs_dec_total = 0;
    int tx_vs_dec_info = 0;
    int tx_vs_dec_parity = 0;
    int det_vs_dec_total = 0;
    int det_vs_dec_info = 0;
    int det_vs_dec_parity = 0;
    const int max_bin = (pckt->bin_num > 0) ? (pckt->bin_num - 1) : 0;
    const int rows = pckt->h_matrix.rows;
    const int cols = pckt->h_matrix.cols;
    const int bits = pckt->h_matrix.bits;

    int row_occ_cnt[LDPC_MAX_ROWS] = {0};
    int row_fade_cnt[LDPC_MAX_ROWS] = {0};
    int col_occ_cnt[LDPC_MAX_COLS] = {0};
    int col_fade_cnt[LDPC_MAX_COLS] = {0};
    int col_mask_bits[LDPC_MAX_COLS] = {0};

    int ibex_in_hard_ones = 0;
    int ibex_in_q1_ones = 0;
    int ibex_in_q2_ones = 0;
    int ibex_in_corrected = 0;
    int ibex_in_error = 0;
    int ibex_in_sw_nonzero = 0;
    int ibex_in_sw_sum = 0;
    int ibex_in_sw_max = 0;
    int ibex_in_level_hist[8] = {0};
    int ibex_out_hard_ones = 0;
    unsigned int ibex_in_hard_sig = 2166136261u;
    unsigned int ibex_in_level_sig = 2166136261u;
    unsigned int ibex_in_sw_sig = 2166136261u;
    unsigned int ibex_out_hard_sig = 2166136261u;

    for (int i = 0; i < pckt->blk_len; i++) {
        int det_hard = 0;
        if (sd_num >= 2) {
            int bin = (int)(unsigned char)pckt->det_blk[i];
            if (bin < 0)
                bin = 0;
            if (bin > max_bin)
                bin = max_bin;
            if (bin < 128)
                det_bin_hist[bin]++;
            det_hard = (pckt->llr_tbl && (pckt->llr_tbl[bin] < 0.0f)) ? 1 : 0;
        } else {
            det_hard = pckt->det_blk ? (pckt->det_blk[i] & 1) : 0;
        }

        const int txb = pckt->tx_blk ? (pckt->tx_blk[i] & 1) : 0;
        const int dcb = pckt->dec_blk ? (pckt->dec_blk[i] & 1) : 0;
        const int is_info = (i < pckt->info_len);

        det_ones += det_hard;
        tx_ones += txb;
        dec_ones += dcb;

        if (pckt->tx_blk) {
            if (txb != det_hard) {
                tx_vs_det_total++;
                if (is_info)
                    tx_vs_det_info++;
                else
                    tx_vs_det_parity++;
            }
            if (txb != dcb) {
                tx_vs_dec_total++;
                if (is_info)
                    tx_vs_dec_info++;
                else
                    tx_vs_dec_parity++;
            }
        }

        if (pckt->dec_blk) {
            if (det_hard != dcb) {
                det_vs_dec_total++;
                if (is_info)
                    det_vs_dec_info++;
                else
                    det_vs_dec_parity++;
            }
        }
    }

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (pckt->h_matrix.occupied[r][c]) {
                row_occ_cnt[r]++;
                col_occ_cnt[c]++;
            }
            if (pckt->h_matrix.fade[r][c]) {
                row_fade_cnt[r]++;
                col_fade_cnt[c]++;
            }
        }
    }
    for (int c = 0; c < cols; c++) {
        for (int k = 0; k < bits; k++) {
            if (pckt->h_matrix.mask[c][k])
                col_mask_bits[c]++;
        }
    }

    for (int c = 0; c < cols; c++) {
        for (int k = 0; k < bits; k++) {
            const s_codeword_bit *inb = &pckt->ldpc_decoder_input.corrupted_codeword.c[c].b[k];
            const int hb = inb->bit_hard ? 1 : 0;
            const int q1 = inb->bit_questionable ? 1 : 0;
            const int q2 = inb->bit_questionable2 ? 1 : 0;
            const int lv = (int)inb->level;
            const int sw = (int)inb->syndrome_weight;
            const int outb = pckt->ldpc_decoder_output.corrected_codeword.c[c].b[k] ? 1 : 0;

            ibex_in_hard_ones += hb;
            ibex_in_q1_ones += q1;
            ibex_in_q2_ones += q2;
            ibex_in_corrected += inb->bit_corrected ? 1 : 0;
            ibex_in_error += inb->bit_is_error ? 1 : 0;
            ibex_in_sw_sum += sw;
            if (sw != 0)
                ibex_in_sw_nonzero++;
            if (sw > ibex_in_sw_max)
                ibex_in_sw_max = sw;
            if ((lv >= 0) && (lv < 8))
                ibex_in_level_hist[lv]++;

            ibex_in_hard_sig ^= (unsigned int)hb;
            ibex_in_hard_sig *= 16777619u;
            ibex_in_level_sig ^= (unsigned int)(lv & 0xFF);
            ibex_in_level_sig *= 16777619u;
            ibex_in_sw_sig ^= (unsigned int)(sw & 0xFFFF);
            ibex_in_sw_sig *= 16777619u;

            ibex_out_hard_ones += outb;
            ibex_out_hard_sig ^= (unsigned int)outb;
            ibex_out_hard_sig *= 16777619u;
        }
    }

    const unsigned int tx_sig =
        (pckt->tx_blk && (pckt->blk_len > 0)) ? dvc_sig_fnv1a_bits01(pckt->tx_blk, 0, pckt->blk_len) : 0u;
    const unsigned int det_sig =
        (pckt->det_blk && (pckt->blk_len > 0)) ? dvc_sig_fnv1a_rx_bins01(pckt, 0, pckt->blk_len) : 0u;
    const unsigned int dec_sig =
        (pckt->dec_blk && (pckt->blk_len > 0)) ? dvc_sig_fnv1a_bits01(pckt->dec_blk, 0, pckt->blk_len) : 0u;

    printf("[DVC OVW] ==================== LDPC OVERVIEW BEGIN ====================\n");
    printf("[DVC OVW] meta  dec_mode_in=%d model=%d(%s) sd_num=%d debug=%d cfg=%d vn_bits=%d\n",
           dec_mode_in, (int)model, model_name, sd_num, debug, g_configured, g_vn_bits);
    printf("[DVC OVW] size  info_len=%d pad_len=%d blk_len=%d hm_k=%d hm_m=%d hm_n=%d bm_m=%d bm_n=%d cir_sz=%d\n",
           pckt->info_len, pckt->pad_len, pckt->blk_len, pckt->hm_k, pckt->hm_m, pckt->hm_n, pckt->bm_m, pckt->bm_n,
           pckt->cir_sz);
    printf("[DVC OVW] hmat  rows=%d cols=%d bits=%d extra_ud_bits=%d extra_pa_bits=%d unused_ud_bytes=%d unused_pa_bytes=%d\n",
           pckt->h_matrix.rows, pckt->h_matrix.cols, pckt->h_matrix.bits, pckt->h_matrix.extra_bits_of_userdata,
           pckt->h_matrix.extra_bits_of_parity, pckt->h_matrix.unused_bytes_of_userdata,
           pckt->h_matrix.unused_bytes_of_parity);
    printf("[DVC OVW] ch    sel=%d(%s) snr=%f snr_code=%f sigma=%f ber=%f err_num=%d raw_err_num=%d\n",
           (int)pckt->ch_sel, ch_name, pckt->snr, pckt->snr_code, pckt->awgn_sigma, pckt->ber, pckt->err_num,
           pckt->raw_err_num);
    printf("[DVC OVW] soft  sd_type=%d rd_num=%d bin_num=%d max_llr_bin=%d llr_tot=%d llr_frac=%d llr_max=%f llr_min=%f\n",
           (int)pckt->sd_type, pckt->rd_num, pckt->bin_num, pckt->max_llr_bin, pckt->llr_tot_num, pckt->llr_frac_num,
           pckt->llr_max, pckt->llr_min);
    printf("[DVC OVW] dec   ldec_max_itr=%d early_term=%d alpha=%f finite(mode=%d q=%d r=%d f=%d) max_llr_bin=%d\n",
           pckt->ldec_max_itr, pckt->ldec_early_term_en, pckt->alpha, pckt->finite_mode, pckt->finite_q_num,
           pckt->finite_r_num, pckt->finite_f_num, pckt->max_llr_bin);
    printf("[DVC OVW] pms   point1=%f point2=%f alpha_pms=[%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f] beta_pms=[%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f]\n",
           pckt->point1, pckt->point2, pckt->alpha_pms[0], pckt->alpha_pms[1], pckt->alpha_pms[2], pckt->alpha_pms[3],
           pckt->alpha_pms[4], pckt->alpha_pms[5], pckt->beta_pms[0], pckt->beta_pms[1], pckt->beta_pms[2],
           pckt->beta_pms[3], pckt->beta_pms[4], pckt->beta_pms[5]);
    printf("[DVC OVW] bf.in nand=%d soft_bits=%d post_iter=%d itr_limit=%d synd_only=%d err_ud=%d err_cw=%d rber=%f verbose=%d\n",
           pckt->ldpc_decoder_input.nand_strobes, pckt->ldpc_decoder_input.soft_bits,
           pckt->ldpc_decoder_input.post_iteration, pckt->ldpc_decoder_input.iteration_limit,
           pckt->ldpc_decoder_input.syndrome_cal_only ? 1 : 0,
           pckt->ldpc_decoder_input.errors_in_userdata, pckt->ldpc_decoder_input.errors,
           pckt->ldpc_decoder_input.rber, pckt->ldpc_decoder_input.verbose ? 1 : 0);
    printf("[DVC OVW] bf.out failure=%d err_ud=%d err_cw=%d iterations=%d cycles=%d synd_before=%d synd_after=%d early_term=%d col_cnt=%d\n",
           pckt->ldpc_decoder_output.failure, pckt->ldpc_decoder_output.errors_in_userdata,
           pckt->ldpc_decoder_output.errors_in_codeword, pckt->ldpc_decoder_output.iterations,
           pckt->ldpc_decoder_output.clock_cycles, pckt->ldpc_decoder_output.syndrome_weight_before,
           pckt->ldpc_decoder_output.syndrome_weight_after, pckt->ldpc_decoder_output.early_termination,
           pckt->ldpc_decoder_output.col_cnt);
    printf("[DVC OVW] bf.par post_en=%d synd_thr_qc=%d synd_thr_post=%d likelihood_thr=%d post_ratio=%d questionable_sense=%d early_term_dis=%d\n",
           pckt->ldpc_decoder_parameters.post_process_en, pckt->ldpc_decoder_parameters.syndrome_weight_thr_qc,
           pckt->ldpc_decoder_parameters.syndrome_weight_thr_post, pckt->ldpc_decoder_parameters.likelihood_thr,
           pckt->ldpc_decoder_parameters.post_ratio, pckt->ldpc_decoder_parameters.questionable_sense ? 1 : 0,
           pckt->ldpc_decoder_parameters.early_terminate_dis ? 1 : 0);
    printf("[DVC OVW] bf.tbl soft_bit_table=");
    for (int i = 0; i < 8; i++)
        printf("%s%d", (i == 0) ? "" : ",", pckt->ldpc_decoder_parameters.soft_bit_table[i]);
    printf(" likelihood_map=");
    for (int i = 0; i < 8; i++)
        printf("%s%d", (i == 0) ? "" : ",", pckt->ldpc_decoder_parameters.likelihood_map[i]);
    printf(" init_frac=%d,%d,%d init_coef=%d,%d,%d,%d\n",
           pckt->ldpc_decoder_parameters.likelihood_init_fraction[0],
           pckt->ldpc_decoder_parameters.likelihood_init_fraction[1],
           pckt->ldpc_decoder_parameters.likelihood_init_fraction[2],
           pckt->ldpc_decoder_parameters.likelihood_init_coef[0],
           pckt->ldpc_decoder_parameters.likelihood_init_coef[1],
           pckt->ldpc_decoder_parameters.likelihood_init_coef[2],
           pckt->ldpc_decoder_parameters.likelihood_init_coef[3]);
    for (int i = 0; i < 8; i++) {
        printf("[DVC OVW] bf.coef row=%d all=%d,%d,%d,%d ethr0=%d ethr1=%d\n",
               i,
               pckt->ldpc_decoder_parameters.likelihood_init_coef_all[i][0],
               pckt->ldpc_decoder_parameters.likelihood_init_coef_all[i][1],
               pckt->ldpc_decoder_parameters.likelihood_init_coef_all[i][2],
               pckt->ldpc_decoder_parameters.likelihood_init_coef_all[i][3],
               pckt->ldpc_decoder_parameters.early_terminate_thr[0][i],
               pckt->ldpc_decoder_parameters.early_terminate_thr[1][i]);
    }

    if (pckt->llr_tbl && (pckt->bin_num > 0)) {
        printf("[DVC OVW] llr   bin_num=%d max_bin=%d llr_max=%f llr_min=%f tbl:",
               pckt->bin_num, pckt->max_llr_bin, pckt->llr_max, pckt->llr_min);
        for (int i = 0; i < pckt->bin_num; i++)
            printf(" [%d]=%0.3f", i, pckt->llr_tbl[i]);
        printf("\n");
    }

    if (pckt->vref && (pckt->rd_num > 0)) {
        printf("[DVC OVW] vref  rd_num=%d vals:", pckt->rd_num);
        for (int i = 0; i < pckt->rd_num; i++)
            printf(" [%d]=%0.3f", i, pckt->vref[i]);
        printf("\n");
    }
    if (pckt->vref_asc_ord && (pckt->rd_num > 0)) {
        printf("[DVC OVW] vasc  vals:");
        for (int i = 0; i < pckt->rd_num; i++)
            printf(" [%d]=%0.3f", i, pckt->vref_asc_ord[i]);
        printf("\n");
    }
    if (pckt->llr_asc_ord && (pckt->bin_num > 0)) {
        printf("[DVC OVW] lasc  vals:");
        for (int i = 0; i < pckt->bin_num; i++)
            printf(" [%d]=%0.3f", i, pckt->llr_asc_ord[i]);
        printf("\n");
    }

    if (pckt->bin_split && (pckt->rd_num > 0)) {
        printf("[DVC OVW] split split_bin:");
        for (int i = 0; i < pckt->rd_num; i++)
            printf(" [%d]=%d", i, (int)(unsigned char)pckt->bin_split[i]);
        printf("\n");
    }

    if (pckt->bin_asc_ord && (pckt->rd_num >= 0)) {
        printf("[DVC OVW] binid bin_asc_ord:");
        for (int i = 0; i <= pckt->rd_num; i++)
            printf(" [%d]=%d", i, (int)(unsigned char)pckt->bin_asc_ord[i]);
        printf("\n");
    }
    if (pckt->sd_asc_ord && (pckt->bin_num > 0)) {
        printf("[DVC OVW] sasc  sd_asc_ord:");
        for (int i = 0; i < pckt->bin_num; i++)
            printf(" [%d]=%d", i, (int)(unsigned char)pckt->sd_asc_ord[i]);
        printf("\n");
    }

    printf("[DVC OVW] stat  cw_fail=%d init_synd=%d fina_synd=%d cnvg_itr=%d cnvg_col=%d\n",
           pckt->cw_fail, pckt->init_synd_wt, pckt->fina_synd_wt, pckt->cnvg_itr, pckt->cnvg_lyr);
    printf("[DVC OVW] bits  tx_ones=%d det_ones=%d dec_ones=%d blk_len=%d\n",
           tx_ones, det_ones, dec_ones, pckt->blk_len);
    printf("[DVC OVW] mis   tx_vs_det=%d(info=%d parity=%d) tx_vs_dec=%d(info=%d parity=%d) det_vs_dec=%d(info=%d parity=%d)\n",
           tx_vs_det_total, tx_vs_det_info, tx_vs_det_parity, tx_vs_dec_total, tx_vs_dec_info, tx_vs_dec_parity,
           det_vs_dec_total, det_vs_dec_info, det_vs_dec_parity);
    printf("[DVC OVW] sig   tx=%08x det=%08x dec=%08x\n", tx_sig, det_sig, dec_sig);
    printf("[DVC OVW] bfsum in_hard_ones=%d q1_ones=%d q2_ones=%d corrected=%d error=%d sw_nonzero=%d sw_sum=%d sw_max=%d out_hard_ones=%d\n",
           ibex_in_hard_ones, ibex_in_q1_ones, ibex_in_q2_ones, ibex_in_corrected, ibex_in_error,
           ibex_in_sw_nonzero, ibex_in_sw_sum, ibex_in_sw_max, ibex_out_hard_ones);
    printf("[DVC OVW] bfsig in_hard=%08x in_level=%08x in_sw=%08x out_hard=%08x level_hist=%d,%d,%d,%d,%d,%d,%d,%d\n",
           ibex_in_hard_sig, ibex_in_level_sig, ibex_in_sw_sig, ibex_out_hard_sig,
           ibex_in_level_hist[0], ibex_in_level_hist[1], ibex_in_level_hist[2], ibex_in_level_hist[3],
           ibex_in_level_hist[4], ibex_in_level_hist[5], ibex_in_level_hist[6], ibex_in_level_hist[7]);
    for (int lv = 0; lv < 4; lv++) {
        printf("[DVC OVW] bflw  lv=%d err=%u,%u,%u,%u,%u cor=%u,%u,%u,%u,%u prob=%0.6f,%0.6f,%0.6f,%0.6f,%0.6f\n",
               lv,
               pckt->ldpc_decoder_output.corrected_codeword.errors_at_level_and_weight[lv][0],
               pckt->ldpc_decoder_output.corrected_codeword.errors_at_level_and_weight[lv][1],
               pckt->ldpc_decoder_output.corrected_codeword.errors_at_level_and_weight[lv][2],
               pckt->ldpc_decoder_output.corrected_codeword.errors_at_level_and_weight[lv][3],
               pckt->ldpc_decoder_output.corrected_codeword.errors_at_level_and_weight[lv][4],
               pckt->ldpc_decoder_output.corrected_codeword.correct_at_level_and_weight[lv][0],
               pckt->ldpc_decoder_output.corrected_codeword.correct_at_level_and_weight[lv][1],
               pckt->ldpc_decoder_output.corrected_codeword.correct_at_level_and_weight[lv][2],
               pckt->ldpc_decoder_output.corrected_codeword.correct_at_level_and_weight[lv][3],
               pckt->ldpc_decoder_output.corrected_codeword.correct_at_level_and_weight[lv][4],
               pckt->ldpc_decoder_output.corrected_codeword.probability_of_error_at_level_and_weight[lv][0],
               pckt->ldpc_decoder_output.corrected_codeword.probability_of_error_at_level_and_weight[lv][1],
               pckt->ldpc_decoder_output.corrected_codeword.probability_of_error_at_level_and_weight[lv][2],
               pckt->ldpc_decoder_output.corrected_codeword.probability_of_error_at_level_and_weight[lv][3],
               pckt->ldpc_decoder_output.corrected_codeword.probability_of_error_at_level_and_weight[lv][4]);
    }

    if (sd_num >= 2 && (pckt->bin_num > 0)) {
        printf("[DVC OVW] dist  det_bin_hist:");
        for (int i = 0; i < pckt->bin_num; i++)
            printf(" b%d=%d", i, det_bin_hist[i]);
        printf("\n");
    }

    if ((rows > 0) && (cols > 0)) {
        printf("[DVC OVW] rowwt");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, pckt->h_matrix.row_weight[r]);
        printf("\n");
        printf("[DVC OVW] rowfe");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, pckt->h_matrix.first_element[r]);
        printf("\n");
        printf("[DVC OVW] rowle");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, pckt->h_matrix.last_element[r]);
        printf("\n");
        printf("[DVC OVW] rowwr");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, pckt->h_matrix.wraparound[r]);
        printf("\n");
        printf("[DVC OVW] rowwb");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, pckt->h_matrix.wrap_base[r]);
        printf("\n");
        printf("[DVC OVW] rownd");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, pckt->h_matrix.wrap_num_deltas[r]);
        printf("\n");
        printf("[DVC OVW] rowoc");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, row_occ_cnt[r]);
        printf("\n");
        printf("[DVC OVW] rowfd");
        for (int r = 0; r < rows; r++)
            printf(" [%d]=%d", r, row_fade_cnt[r]);
        printf("\n");

        printf("[DVC OVW] colwt");
        for (int c = 0; c < cols; c++)
            printf(" [%d]=%d", c, pckt->h_matrix.col_weight[c]);
        printf("\n");
        printf("[DVC OVW] colpc");
        for (int c = 0; c < cols; c++)
            printf(" [%d]=%d", c, pckt->h_matrix.parity_column[c] ? 1 : 0);
        printf("\n");
        printf("[DVC OVW] coloc");
        for (int c = 0; c < cols; c++)
            printf(" [%d]=%d", c, col_occ_cnt[c]);
        printf("\n");
        printf("[DVC OVW] colfd");
        for (int c = 0; c < cols; c++)
            printf(" [%d]=%d", c, col_fade_cnt[c]);
        printf("\n");
        printf("[DVC OVW] colmk");
        for (int c = 0; c < cols; c++)
            printf(" [%d]=%d", c, col_mask_bits[c]);
        printf("\n");
    }

    printf("[DVC OVW] ===================== LDPC OVERVIEW END =====================\n");
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

static unsigned int dvc_get_bit_from_words_msb32(const int *words, int bit_index)
{
    if (!words || (bit_index < 0))
        return 0u;
    const int w = bit_index / 32;
    const int j = bit_index % 32;
    return ((unsigned int)words[w] >> (31 - j)) & 1u;
}

static void dvc_dump_words_hex_rows(const char *tag, const int *words, int n_words)
{
    if (!tag || !words || (n_words <= 0))
        return;
    printf("%s\n", tag);
    for (int w = 0; w < n_words; w++)
        printf("%04d %08x\n", w, (unsigned int)words[w]);
}

static void dvc_dump_words_hex_rows_fp(FILE *fp, const char *tag, const int *words, int n_words)
{
    if (!fp || !tag || !words || (n_words <= 0))
        return;
    fprintf(fp, "%s\n", tag);
    for (int w = 0; w < n_words; w++)
        fprintf(fp, "%04d %08x\n", w, (unsigned int)words[w]);
}

static void dvc_dump_sd_rd_planes_hex_rows(const int *rd_data, int rd_num, int bits_per_plane)
{
    if (!rd_data || (rd_num <= 0) || (bits_per_plane <= 0))
        return;

    const int words_per_plane = (bits_per_plane + 31) / 32;
    char tag[96];
    for (int rd_indx = 0; rd_indx < rd_num; rd_indx++) {
        snprintf(tag, sizeof(tag), "[DVC IBEX] sd_err_inj rd_plane[%d] bits=%d words=%d", rd_indx, bits_per_plane,
                 words_per_plane);
        dvc_dump_words_hex_rows(tag, rd_data + rd_indx * words_per_plane, words_per_plane);
    }
}

static void dvc_dump_bins_hex_rows_fp(FILE *fp, const char *tag, const int *bins, int n_bins)
{
    if (!fp || !tag || !bins || (n_bins <= 0))
        return;
    fprintf(fp, "%s\n", tag);
    for (int base = 0; base < n_bins; base += 32) {
        fprintf(fp, "%05d ", base);
        for (int i = 0; (i < 32) && ((base + i) < n_bins); i++)
            fprintf(fp, "%x", ((unsigned int)bins[base + i]) & 0xFu);
        fprintf(fp, "\n");
    }
}

static void dvc_dump_bins_hex_range_fp(FILE *fp, const char *tag, const int *bins, int start, int len)
{
    if (!fp || !tag || !bins)
        return;
    if (start < 0)
        start = 0;
    if (len < 0)
        len = 0;

    fprintf(fp, "%s range=[%d,%d) len=%d\n", tag, start, start + len, len);
    for (int base = 0; base < len; base += 32) {
        fprintf(fp, "%05d ", start + base);
        for (int i = 0; (i < 32) && ((base + i) < len); i++)
            fprintf(fp, "%x", ((unsigned int)bins[start + base + i]) & 0xFu);
        fprintf(fp, "\n");
    }
}

static void dvc_dump_int_array_fp(FILE *fp, const char *tag, const int *vals, int n_vals)
{
    if (!fp || !tag || !vals || (n_vals <= 0))
        return;
    fprintf(fp, "%s", tag);
    for (int i = 0; i < n_vals; i++)
        fprintf(fp, "%s%d", (i == 0) ? "" : " ", vals[i]);
    fprintf(fp, "\n");
}

static void dvc_dump_float_array_fp(FILE *fp, const char *tag, const float *vals, int n_vals)
{
    if (!fp || !tag || !vals || (n_vals <= 0))
        return;
    fprintf(fp, "%s", tag);
    for (int i = 0; i < n_vals; i++)
        fprintf(fp, "%s%.6f", (i == 0) ? "" : " ", vals[i]);
    fprintf(fp, "\n");
}

static void dvc_dump_sd_rd_planes_hex_rows_fp(FILE *fp, const int *rd_data, int rd_num, int bits_per_plane)
{
    if (!fp || !rd_data || (rd_num <= 0) || (bits_per_plane <= 0))
        return;

    const int words_per_plane = (bits_per_plane + 31) / 32;
    char tag[96];
    for (int rd_indx = 0; rd_indx < rd_num; rd_indx++) {
        snprintf(tag, sizeof(tag), "[RD_PLANE %d] bits=%d words=%d", rd_indx, bits_per_plane, words_per_plane);
        dvc_dump_words_hex_rows_fp(fp, tag, rd_data + rd_indx * words_per_plane, words_per_plane);
    }
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
                 int ch_llr_mode,
                 int fdec_max_itr,
                 int ldec_max_itr,
                 int llr0,
                 int llr1,
                 int alpha,
                 svOpenArrayHandle alpha_sv,
                 svOpenArrayHandle beta_sv,
                 int point1,
                 int point2,
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
    const float m_point1 = (float)point1 / 1000.0f;
    const float m_point2 = (float)point2 / 1000.0f;
    float alpha_pms[LDPC_PMS_LUT_SIZE];
    float beta_pms[LDPC_PMS_LUT_SIZE];
    for (int i = 0; i < LDPC_PMS_LUT_SIZE; i++)
        alpha_pms[i] = m_alpha;
    for (int i = 0; i < LDPC_PMS_LUT_SIZE; i++)
        beta_pms[i] = 0.0f;
    const int alpha_len = svHigh(alpha_sv, 1) + 1;
    int *alpha_in = (int *)svGetArrayPtr(alpha_sv);
    if ((alpha_in != NULL) && (alpha_len > 0)) {
        int copy_cnt = alpha_len;
        if (copy_cnt > LDPC_PMS_LUT_SIZE)
            copy_cnt = LDPC_PMS_LUT_SIZE;
        for (int i = 0; i < copy_cnt; i++)
            alpha_pms[i] = (float)alpha_in[i] / 1000.0f;
    }
    const int beta_len = svHigh(beta_sv, 1) + 1;
    int *beta_in = (int *)svGetArrayPtr(beta_sv);
    if ((beta_in != NULL) && (beta_len > 0)) {
        int copy_cnt = beta_len;
        if (copy_cnt > LDPC_PMS_LUT_SIZE)
            copy_cnt = LDPC_PMS_LUT_SIZE;
        for (int i = 0; i < copy_cnt; i++)
            beta_pms[i] = (float)beta_in[i] / 1000.0f;
    }
    const float m_llr0 = (float)llr0 / kFiniteFracScale;
    const float m_llr1 = (float)llr1 / kFiniteFracScale;

    g_bm_m = h_m;
    g_bm_n = h_n;
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
    else
        printf("[DVC IBEX WARN] unknown ch_mode=%d, fallback CLEAN\n", ch_mode);

    sim_pckt->ch_config(g_info_len, g_blk_len, ch_sel, m_ch_para, g_vn_bits);

    // LDPC configuration
    sim_pckt->ldpc_config(g_bm_m, g_bm_n, h_sc, h_st, h_wt);

    // packet allocation
    sim_pckt->ldpc_pckt_alloc();

    // LLR tables (reuse IBEX transceiver implementation)
    if (ch_llr_mode == 0)
        sim_pckt->ch_llr_alloc(MANUAL, rd_num, vref);
    else if (ch_llr_mode == 1)
        sim_pckt->ch_llr_alloc(VENDOR0, log2(rd_num + 1), vref);
    else if (ch_llr_mode == 2)
        sim_pckt->ch_llr_alloc(VENDOR1, log2(rd_num + 1), vref);

    printf("[LLR DBG] before ch_llr_gen: ch_sel=%d snr=%f snr_code=%f awgn_sigma=%f sigma2=%f llr_tot_arg=%d llr_frac_arg=%d\n",
           (int)sim_pckt->ch_sel, sim_pckt->snr, sim_pckt->snr_code, sim_pckt->awgn_sigma,
           sim_pckt->awgn_sigma * sim_pckt->awgn_sigma, finite_r_num, kFiniteFracBits);
    sim_pckt->ch_llr_gen(m_llr0, m_llr1, finite_r_num, kFiniteFracBits);

    // decoder config
    sim_pckt->ldpc_dec_config(fdec_max_itr, 0, ldec_max_itr, m_alpha, alpha_pms, beta_pms, m_point1, m_point2, LDPC_PMS_LUT_SIZE, 1, finite_q_num, finite_r_num, kFiniteFracBits, sdlite_llr_config,
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

    sim_pckt->ldpc_ibex_encoder();

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
    const int det_elems = svHigh(det_data_sv, 1) + 1;
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
    }

    // 2) Decode mode mapping (keep gen4 interface behavior, but allow IBEX-native values too)
    enum dec_model model = LAYER;
    if (dec_mode == (int)BF_IBEX)
        model = BF_IBEX;
    else if (dec_mode == (int)LAYER_G2)
        model = LAYER_G2;

    if (debug >= 1) {
        printf("[DVC IBEX] ldpc_dec: dec_mode=%d -> model=%d sd_num=%d blk_len=%d n_words=%d bin_num=%d\n", dec_mode,
               (int)model, sd_num, sim_pckt->blk_len, n_words, sim_pckt->bin_num);
        printf("[DVC IBEX] ldpc_dec sizes: info_len=%d pad_len=%d hm_k=%d hm_m=%d hm_n=%d parity_det=%d\n",
               sim_pckt->info_len, sim_pckt->pad_len, sim_pckt->hm_k, sim_pckt->hm_m, sim_pckt->hm_n,
               sim_pckt->blk_len - sim_pckt->info_len);
    }

    // For LAYER_G2 decoding, dec_di_blk is used as an index into llr_tbl (see ldpc_dec_layer/ldpc_dec_layer2).
    // If hard bits (0/1) are passed with sd_num<2, they must be mapped to "strong-0" / "strong-1" bins first.
    if (model == LAYER_G2) {
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

            if (sim_pckt->tx_blk) {
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
    if ((sd_num >= 2) && (model != LAYER_G2)) {
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
        int max_bin = sim_pckt->bin_num - 1;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            int bin = (int)(unsigned char)sim_pckt->det_blk[i];
            if (bin > max_bin)
                bin = max_bin;
            sim_pckt->rx_blk[i] = sim_pckt->llr_tbl[bin];
        }
        sim_pckt->ldpc_ibex_input(/*syndrome_cal_only=*/0, g_max_iter, g_post_iter, g_nand_strobes);
    }

    // 3) Run decode (IBEX codec handles [Info|Pad|Parity] reconstruction internally)
    sim_pckt->ldpc_decoder(model);

    if ((debug >= 1) && (model == LAYER_G2)) {
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
        dvc_dump_codec_overview(sim_pckt, dec_mode, model, sd_num, debug);
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
        if (debug >= 1) {
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
        if ((debug >= 1) && (model == LAYER_G2) && sim_pckt->qc_bm && sim_pckt->dec_do_blk && (sim_pckt->cir_sz > 0) &&
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

    const float m_ch_para = (float)ch_para / 1000.0f;

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

    // 1) Unpack TX (MSB-first) into tx_blk[0..blk_len-1].
    //    Support two caller layouts:
    //    - BLK layout: [info][parity_det] (len=blk_len)
    //    - HM  layout: [info][pad][parity_full] (len=hm_n), then we skip pad/mask.
    const int tx_words = svHigh(tx_data_sv, 1) + 1; 
    const int tx_bits_capacity = (tx_words > 0) ? (tx_words * 32) : 0;
    const int parity_bits_in_det = sim_pckt->blk_len - sim_pckt->info_len;
    const int use_hm_layout =
        (sim_pckt->pad_len > 0) && (sim_pckt->hm_n > sim_pckt->blk_len) &&
        (tx_bits_capacity >= sim_pckt->hm_n) && (sim_pckt->hm_k >= sim_pckt->info_len) &&
        (sim_pckt->hm_k + parity_bits_in_det <= sim_pckt->hm_n);

    if (use_hm_layout) {
        for (int i = 0; i < sim_pckt->info_len; i++)
            sim_pckt->tx_blk[i] = (char)dvc_get_bit_from_words_msb32(tx_data, i);
        for (int i = 0; i < parity_bits_in_det; i++) {
            sim_pckt->tx_blk[sim_pckt->info_len + i] = 
                (char)dvc_get_bit_from_words_msb32(tx_data, sim_pckt->hm_k+i);
        }
        if (debug >= 1) {
            printf("[DVC IBEX] ch_err_inj TX layout=HM_N (skip user pad/parity mask): tx_bits=%d hm_n=%d blk_len=%d info=%d hm_k=%d\n",
                   tx_bits_capacity, sim_pckt->hm_n, sim_pckt->blk_len, sim_pckt->info_len, sim_pckt->hm_k);
        }
    } else {
        unsigned int tmp = 0;
        int j = 0;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            if (j == 0)
                tmp = (unsigned int)tx_data[i / 32];
            sim_pckt->tx_blk[i] = (char)((tmp >> (31 - j)) & 1u);
            j++;
            j = j % 32;
        }
    }

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
            if ((sim_pckt->h_matrix.unused_bytes_of_parity > 0) && (parity_bits_in_det < sim_pckt->hm_m)) {
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
                        (char)(sim_pckt->tx_blk[sim_pckt->info_len + first_parity_bits + i] & 1); //(char)(sim_pckt->tx_blk[sim_pckt->info_len + i] & 1);
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
        else if (sim_pckt->ch_sel == AWGN)
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
    const int tx_words = svHigh(tx_data_sv, 1) + 1;
    const int rd_elems = svHigh(rd_data_sv, 1) + 1;

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

    // 1) Unpack TX (MSB-first) into tx_blk[0..blk_len-1].
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

    for (int i = 0; i < sim_pckt->blk_len; i++)
        sim_pckt->det_blk[i] = 0;

    // 3) Read data out (packed bits), rd_data order: rd_indx major, then bit index.
    // Keep old BLK-length behavior unless DV provides a larger buffer that can
    // hold full/HM_N per-plane data, in which case zero-fill the tail bits in
    // each read plane. RDEC soft-bin input is driven by bin_data/vdr_dat.
    const int rd_words_blk = (sim_pckt->blk_len + 31) / 32;
    const int rd_words_hmn = (sim_pckt->hm_n + 31) / 32;
    const int use_rd_hm_layout =
        (sim_pckt->hm_n > sim_pckt->blk_len) &&
        (rd_elems >= rd_num * rd_words_hmn);
    const int rd_bits_per_plane = use_rd_hm_layout ? sim_pckt->hm_n : sim_pckt->blk_len;
    const int rd_words_per_plane = use_rd_hm_layout ? rd_words_hmn : rd_words_blk;
    for (int i = 0; i < rd_elems; i++)
        rd_data[i] = 0;

    for (int rd_indx = 0; rd_indx < rd_num; rd_indx++) {
        tmp = 0;
        j = 0;
        int k = rd_indx * rd_words_per_plane;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            const unsigned int bit = (sim_pckt->rx_blk[i] >= vref[rd_indx]) ? 0u : 1u;
            tmp = tmp * 2u + bit;
            j++;
            j = j % 32;
            if ((j == 0) || (i == (sim_pckt->blk_len - 1))) {
                if (j == 0)
                    rd_data[k] = (int)tmp;
                else
                    rd_data[k] |= (int)tmp << (32 - j);
                tmp = 0;
                k++;
            }
        }
    }

    // 4) Detector bin out.
    // MANUAL uses the split-bin sequence from vref input order.
    // VENDOR0/VENDOR1 use the sorted-threshold interval mapping generated by
    // ch_llr_gen(), matching ch_packet::ch_detector().
    if (sim_pckt->sd_type == MANUAL) {
        for (int i = 0; i < sim_pckt->blk_len; i++)
            sim_pckt->det_blk[i] = 0;

        for (int rd_indx = 0; rd_indx < rd_num; rd_indx++) {
            const int split = (sim_pckt->bin_split && (rd_indx < sim_pckt->rd_num))
                                  ? (int)(unsigned char)sim_pckt->bin_split[rd_indx]
                                  : 0;
            for (int i = 0; i < sim_pckt->blk_len; i++) {
                const unsigned int bit = (sim_pckt->rx_blk[i] >= vref[rd_indx]) ? 0u : 1u;
                if (((int)(unsigned char)sim_pckt->det_blk[i] == split) && (bit == 1u))
                    sim_pckt->det_blk[i] = (char)(rd_indx + 1);
            }
        }
    } else if ((sim_pckt->sd_type == VENDOR0) || (sim_pckt->sd_type == VENDOR1)) {
        const int cfg_rd_num = sim_pckt->rd_num;
        for (int i = 0; i < sim_pckt->blk_len; i++) {
            int bin = (sim_pckt->bin_asc_ord && (cfg_rd_num >= 0))
                          ? (int)(unsigned char)sim_pckt->bin_asc_ord[cfg_rd_num]
                          : 0;
            for (int rd_indx = 0; rd_indx < cfg_rd_num; rd_indx++) {
                if (sim_pckt->rx_blk[i] < sim_pckt->vref_asc_ord[rd_indx]) {
                    bin = (int)(unsigned char)sim_pckt->bin_asc_ord[rd_indx];
                    break;
                }
            }
            sim_pckt->det_blk[i] = (char)bin;
        }
    }

    // 5) Bin out (per-bit integer).
    // Keep transmitted-block order regardless of rx_data_sv capacity:
    //   [User Data][User Padding=max][Valid Parity][PMSK_tail]
    const int bin_elems = svHigh(rx_data_sv, 1) + 1;
    for (int i = 0; i < bin_elems; i++)
        bin_data[i] = 0;

    const int use_bin_hm_layout = dvc_fill_sd_bin_data_max_llr_view(sim_pckt, bin_data, bin_elems);

    const int bin_export_bits = use_bin_hm_layout ? sim_pckt->hm_n : sim_pckt->blk_len;
    const int bin_tail_zero_bits = (bin_elems > bin_export_bits) ? (bin_elems - bin_export_bits) : 0;

    if (debug >= 1) {
        printf("[DVC IBEX] sd_err_inj BIN layout=%s: bin_elems=%d hm_n=%d blk_len=%d pad_bits=%d parity_mask_bits=%d tail_zero_bits=%d max_llr_bin=%d\n",
               use_bin_hm_layout ? "HM_N+tail_pmsk" : "BLK+tail0",
               bin_elems, sim_pckt->hm_n, sim_pckt->blk_len, sim_pckt->pad_len,
               sim_pckt->h_matrix.unused_bytes_of_parity * 8, bin_tail_zero_bits, sim_pckt->max_llr_bin);
        printf("[DVC IBEX] sd_err_inj RD layout=%s: rd_elems=%d words_per_plane=%d bits_per_plane=%d tail_zero_bits=%d\n",
               use_rd_hm_layout ? "HM_N+tail0" : "BLK",
               rd_elems, rd_words_per_plane, rd_bits_per_plane,
               (rd_bits_per_plane > sim_pckt->blk_len) ? (rd_bits_per_plane - sim_pckt->blk_len) : 0);
        printf("[DVC IBEX] sd_err_inj final bin_data is a single per-bit bin-index array; per-read bit planes are in rd_data[0..sd_num-1].\n");
        dvc_dump_sd_rd_planes_hex_rows(rd_data, rd_num, rd_bits_per_plane);
    }

    // 6) LLR table (scaled by kFiniteFracScale for fixed-point export), fill up to 128 elements
    for (int i = 0; i < 128; i++) {
        if (sim_pckt->llr_tbl && (i < sim_pckt->bin_num))
            llr_tbl[i] = (int)(sim_pckt->llr_tbl[i] * kFiniteFracScale);
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
