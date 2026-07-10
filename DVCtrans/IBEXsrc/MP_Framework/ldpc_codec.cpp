#include <algorithm>
#include <cctype>
#include <stdint.h>
#include <cstring>
#include <errno.h>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "finite_lib.h"
#include "ldpc_codec.h"
#include "ibex_rtl_engine.h"
#include "mod2convert.h"
#include "mod2dense.h"
#include "mod2sparse.h"
#include "vec_op.h"

static void dvc_bins_to_hard_bits(const char *bins, char *hard_bits, int n_bits,
                                  const float *llr_tbl, int bin_num) {
  if (!bins || !hard_bits || n_bits <= 0)
    return;
  if (!llr_tbl || bin_num <= 0) {
    for (int i = 0; i < n_bits; i++)
      hard_bits[i] = (char)(bins[i] & 1);
    return;
  }

  const int max_bin = (bin_num > 0) ? (bin_num - 1) : 0;
  for (int i = 0; i < n_bits; i++) {
    int bin = (int)(unsigned char)bins[i];
    if (bin > max_bin)
      bin = max_bin;
    hard_bits[i] = (llr_tbl[bin] < 0.0f) ? 1 : 0;
  }
}

static void dvc_mkdir_if_missing(const char *path) {
  if (!path || !*path)
    return;
  struct stat st;
  if (stat(path, &st) == 0) {
    if (!S_ISDIR(st.st_mode))
      printf("[LDPC WARN] path exists but is not a directory: %s\n", path);
    return;
  }
  if (mkdir(path, 0777) != 0 && errno != EEXIST)
    printf("[LDPC WARN] mkdir failed: %s (errno=%d)\n", path, errno);
}

static void dvc_ensure_output_c_code_dir() {
  static int done = 0;
  if (done)
    return;
  done = 1;
  dvc_mkdir_if_missing("./output");
  dvc_mkdir_if_missing("./output/c_code");
}

static int dvc_qc_syndrome_weight(mod2sparse *qc_bm, const char *hard_bits,
                                  int bm_m, int cir_sz) {
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
    for (mod2entry *e = mod2sparse_first_in_row(qc_bm, layer);
         !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e)) {
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

static int dvc_ibex_syndrome_weight(ldpc_packet *packet,
                                    const char *hard_bits) {
  if (!packet || !hard_bits)
    return 0;

  s_hard_codeword hard_codeword;
  memset(&hard_codeword, 0, sizeof(hard_codeword));

  const int cols = packet->h_matrix.cols;
  const int bits = packet->h_matrix.bits;
  for (int j = 0; j < cols; j++) {
    for (int k = 0; k < bits; k++) {
      hard_codeword.c[j].b[k] = ((hard_bits[j * bits + k] & 1) != 0);
    }
  }

  s_check_nodes cn = packet->f_check_nodes(packet->h_matrix, hard_codeword);
  return packet->f_check_node_weight(packet->h_matrix, cn);
}

static FILE *dvc_open_bf_ibex_sw_delta_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  dvc_mkdir_if_missing("./output");
  snprintf(path, sizeof(path), "./output/ldpc_dbg_bf_ibex_sw_delta_cw%04u.txt",
           cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open BF IBEX sw-delta dump: %s\n", path);
  return fp;
}

static FILE *dvc_open_bf_ibex_row_sw_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  dvc_mkdir_if_missing("./output");
  snprintf(path, sizeof(path), "./output/ldpc_dbg_bf_ibex_row_sw_cw%04u.txt",
           cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open BF IBEX row-sw dump: %s\n", path);
  return fp;
}

static void dvc_log_bf_ibex_sw(FILE *fp, int iteration, const char *phase,
                               int col, int syndrome_weight) {
  if (!fp || !phase)
    return;
  fprintf(fp, "%8d %12s %6d %16d\n", iteration, phase, col, syndrome_weight);
}

static void dvc_log_bf_ibex_row_sw_header(FILE *fp, int rows) {
  if (!fp)
    return;
  fprintf(fp, "%8s %12s %6s %16s", "iter", "phase", "col", "syndrome_weight");
  for (int row = 0; row < rows; row++)
    fprintf(fp, " row%02d_sw", row);
  fputc('\n', fp);
}

static void dvc_log_bf_ibex_row_sw(FILE *fp, int iteration, const char *phase,
                                   int col, const s_h_matrix &h_matrix,
                                   const s_check_nodes &cn) {
  if (!fp || !phase)
    return;

  int row_sw[LDPC_MAX_ROWS];
  int total_sw = 0;
  for (int row = 0; row < h_matrix.rows; row++) {
    row_sw[row] = 0;
    for (int bit = 0; bit < h_matrix.bits; bit++)
      row_sw[row] += cn.r[row].b[bit] ? 1 : 0;
    total_sw += row_sw[row];
  }

  fprintf(fp, "%8d %12s %6d %16d", iteration, phase, col, total_sw);

  for (int row = 0; row < h_matrix.rows; row++) {
    fprintf(fp, " %8d", row_sw[row]);
  }
  fputc('\n', fp);
}

#if defined(_LDPC_DBG_DUMP) || defined(_LDPC_DEBUG_DUMP)
static char dvc_hex_digit(unsigned int val) {
  static const char kHex[] = "0123456789ABCDEF";
  return kHex[val & 0xF];
}

static void dvc_dump_hard_bits_hex128(const char *path, const char *bits,
                                      int nbits) {
  if (!path || !bits || nbits <= 0)
    return;

  dvc_ensure_output_c_code_dir();
  FILE *fp = fopen(path, "w");
  if (!fp) {
    printf("[LDPC WARN] failed to open hard-decision dump: %s\n", path);
    return;
  }

  const int bits_per_line = 128;
  for (int base = 0; base < nbits; base += bits_per_line) {
    const int line_bits = std::min(bits_per_line, nbits - base);
    unsigned int nibble = 0;
    int nibble_bits = 0;

    for (int i = 0; i < line_bits; i++) {
      nibble = (nibble << 1) | (unsigned int)(bits[base + i] & 1);
      nibble_bits++;
      if (nibble_bits == 4) {
        fputc(dvc_hex_digit(nibble), fp);
        nibble = 0;
        nibble_bits = 0;
      }
    }

    if (nibble_bits != 0) {
      nibble <<= (4 - nibble_bits);
      fputc(dvc_hex_digit(nibble), fp);
    }
    fputc('\n', fp);
  }

  fclose(fp);
}
#endif

#ifdef _LDPC_DBG_DUMP
static FILE *dvc_open_bf_ibex_like_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  snprintf(path, sizeof(path), "./output/ldpc_dbg_bf_ibex_like_cw%04u.txt",
           cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open BF IBEX likelihood dump: %s\n", path);
  return fp;
}

static void dvc_dump_bf_ibex_col_likelihood(FILE *fp,
                                            const s_variable_nodes &vn,
                                            int itr_cnt, int col_cnt,
                                            int cir_bits) {
  if (!fp || (col_cnt < 0) || (cir_bits <= 0))
    return;

  fprintf(fp, "itr_cnt:%03d col_cnt:%03d:", itr_cnt, col_cnt);
  for (int k = cir_bits - 1; k >= 0; k--) {
    fprintf(fp, "%X", (unsigned int)vn.c[col_cnt].b[k].likelihood);
  }
  fputc('\n', fp);
}

static void dvc_fprint_bool_bits_hex(FILE *fp, const bool *bits, int nbits,
                                     int hex_wrap, const char *prefix_first,
                                     const char *prefix_cont) {
  if (!fp || !bits || (nbits <= 0) || (hex_wrap <= 0))
    return;

  const int total_hex = (nbits + 3) / 4;
  int hex_count = 0;
  if (prefix_first)
    fputs(prefix_first, fp);

  unsigned int nibble = 0;
  int nibble_bits = 0;
  for (int bit = nbits - 1; bit >= 0; bit--) {
    nibble = (nibble << 1) | (bits[bit] ? 1u : 0u);
    nibble_bits++;
    if (nibble_bits == 4) {
      fputc(dvc_hex_digit(nibble), fp);
      hex_count++;
      if ((hex_count < total_hex) && ((hex_count % hex_wrap) == 0)) {
        fputc('\n', fp);
        if (prefix_cont)
          fputs(prefix_cont, fp);
      }
      nibble = 0;
      nibble_bits = 0;
    }
  }

  if (nibble_bits != 0) {
    nibble <<= (4 - nibble_bits);
    fputc(dvc_hex_digit(nibble), fp);
    hex_count++;
    if ((hex_count < total_hex) && ((hex_count % hex_wrap) == 0)) {
      fputc('\n', fp);
      if (prefix_cont)
        fputs(prefix_cont, fp);
    }
  }
}

static void dvc_fprint_shifted_bool_bits_hex(FILE *fp, const bool *bits,
                                             int total_bits, int bit_lo,
                                             int nbits, int shift,
                                             bool right_shift) {
  if (!fp || !bits || (total_bits <= 0) || (nbits <= 0) || (bit_lo < 0) ||
      ((bit_lo + nbits) > total_bits))
    return;

  shift %= total_bits;
  if (shift < 0)
    shift += total_bits;

  const int total_hex = (nbits + 3) / 4;
  int hex_count = 0;
  unsigned int nibble = 0;
  int nibble_bits = 0;
  for (int rel = nbits - 1; rel >= 0; rel--) {
    const int bit_idx = bit_lo + rel;
    const int src_idx = right_shift
                            ? (bit_idx + total_bits - shift) % total_bits
                            : (bit_idx + shift) % total_bits;
    nibble = (nibble << 1) | (bits[src_idx] ? 1u : 0u);
    nibble_bits++;
    if (nibble_bits == 4) {
      fputc(dvc_hex_digit(nibble), fp);
      hex_count++;
      if ((hex_count < total_hex) && ((hex_count % 1024) == 0))
        fputc('\n', fp);
      nibble = 0;
      nibble_bits = 0;
    }
  }

  if (nibble_bits != 0) {
    nibble <<= (4 - nibble_bits);
    fputc(dvc_hex_digit(nibble), fp);
    hex_count++;
    if ((hex_count < total_hex) && ((hex_count % 1024) == 0))
      fputc('\n', fp);
  }
}

static void dvc_dump_cn_synd_file(const char *path, const s_check_nodes &cn,
                                  const s_h_matrix &h_matrix) {
  FILE *fp = fopen(path, "w");
  if (!fp) {
    printf("[LDPC WARN] failed to open BF IBEX syndrome dump: %s\n", path);
    return;
  }

  const int syndrome_dump_bits = (h_matrix.bits >= 512) ? 512 : h_matrix.bits;
  const int synd_chunk_bits = 128;
  for (int rr = 0; rr < h_matrix.rows; rr++) {
    for (int local_hi = syndrome_dump_bits; local_hi > 0;
         local_hi -= synd_chunk_bits) {
      const int local_lo = std::max(0, local_hi - synd_chunk_bits);
      const int chunk_len = local_hi - local_lo;
      char line_prefix[32];
      snprintf(line_prefix, sizeof(line_prefix), "row%02d ", rr);
      fputs(line_prefix, fp);
      dvc_fprint_bool_bits_hex(fp, cn.r[rr].b + local_lo, chunk_len, 1024,
                               NULL, NULL);
      fputc('\n', fp);
    }
  }

  fclose(fp);
}

static void dvc_dump_cn_synd_shifted_file(const char *path, const s_check_nodes &cn,
                                          const s_h_matrix &h_matrix,
                                          bool right_shift) {
  FILE *fp = fopen(path, "w");
  if (!fp) {
    printf("[LDPC WARN] failed to open BF IBEX shifted syndrome dump: %s\n",
           path);
    return;
  }

  const int syndrome_dump_bits = (h_matrix.bits >= 512) ? 512 : h_matrix.bits;
  const int synd_chunk_bits = 128;
  for (int col = 0; col < h_matrix.cols; col++) {
    fprintf(fp, "col%02d\n", col);
    for (int rr = 0; rr < h_matrix.rows; rr++) {
      const int shift = h_matrix.element[rr][col];

      for (int local_hi = syndrome_dump_bits; local_hi > 0;
           local_hi -= synd_chunk_bits) {
        const int local_lo = std::max(0, local_hi - synd_chunk_bits);
        const int chunk_len = local_hi - local_lo;
        char line_prefix[32];
        snprintf(line_prefix, sizeof(line_prefix), "row%02d s=%03d ", rr,
                 shift);
        fputs(line_prefix, fp);
        if (shift >= 0)
          dvc_fprint_shifted_bool_bits_hex(fp, cn.r[rr].b, h_matrix.bits,
                                           local_lo, chunk_len, shift,
                                           right_shift);
        else
          dvc_fprint_bool_bits_hex(fp, cn.r[rr].b + local_lo, chunk_len, 1024,
                                   NULL, NULL);
        fputc('\n', fp);
      }
    }
    fputc('\n', fp);
  }

  fclose(fp);
}

static void dvc_dump_cn_synd_shifted_col_file(const char *path,
                                              const s_check_nodes &cn,
                                              const s_h_matrix &h_matrix,
                                              int col, bool right_shift) {
  FILE *fp = fopen(path, "w");
  if (!fp) {
    printf("[LDPC WARN] failed to open BF IBEX shifted syndrome dump: %s\n",
           path);
    return;
  }

  const int syndrome_dump_bits = (h_matrix.bits >= 512) ? 512 : h_matrix.bits;
  const int synd_chunk_bits = 128;
  fprintf(fp, "col%02d\n", col);
  for (int rr = 0; rr < h_matrix.rows; rr++) {
    const int shift = h_matrix.element[rr][col];

    for (int local_hi = syndrome_dump_bits; local_hi > 0;
         local_hi -= synd_chunk_bits) {
      const int local_lo = std::max(0, local_hi - synd_chunk_bits);
      const int chunk_len = local_hi - local_lo;
      char line_prefix[32];
      snprintf(line_prefix, sizeof(line_prefix), "row%02d s=%03d ", rr, shift);
      fputs(line_prefix, fp);
      if (shift >= 0)
        dvc_fprint_shifted_bool_bits_hex(fp, cn.r[rr].b, h_matrix.bits,
                                         local_lo, chunk_len, shift,
                                         right_shift);
      else
        dvc_fprint_bool_bits_hex(fp, cn.r[rr].b + local_lo, chunk_len, 1024,
                                 NULL, NULL);
      fputc('\n', fp);
    }
  }
  fputc('\n', fp);
  fclose(fp);
}
#endif

#ifdef _LDPC_DEBUG_DUMP
static void dvc_fprint_bits_hex(FILE *fp, const char *bits, int nbits) {
  if (!fp || !bits || (nbits <= 0))
    return;

  unsigned int nibble = 0;
  int nibble_bits = 0;
  for (int bit = nbits - 1; bit >= 0; bit--) {
    nibble = (nibble << 1) | (unsigned int)(bits[bit] & 0x1);
    nibble_bits++;
    if (nibble_bits == 4) {
      fputc(dvc_hex_digit(nibble), fp);
      nibble = 0;
      nibble_bits = 0;
    }
  }

  if (nibble_bits != 0) {
    nibble <<= (4 - nibble_bits);
    fputc(dvc_hex_digit(nibble), fp);
  }
}

static void dvc_dump_col_major_codeword(FILE *fp, const char *bits, int cols,
                                        int cir_bits, const char *tag) {
  if (!fp || !bits || (cols <= 0) || (cir_bits <= 0) || !tag)
    return;

  fprintf(fp, "[DVC] %s\n", tag);
  for (int col = 0; col < cols; col++) {
    fprintf(fp, "COL%02d:", col);
    dvc_fprint_bits_hex(fp, bits + col * cir_bits, cir_bits);
    fputc('\n', fp);
  }
}
#endif

void ldpc_packet::ldpc_ibex_phck(s_h_matrix h_matrix) {
  mod2entry *e;
  total_cir = 0;

  qc_bm = mod2sparse_allocate(bm_m, bm_n);
  std::vector<int> row_wt(static_cast<size_t>(bm_m), 0);
  std::vector<int> col_wt(static_cast<size_t>(bm_n), 0);

  for (int i = 0; i < bm_m; i++) {
    for (int j = 0; j < bm_n; j++) {
      const bool should_insert = ((h_matrix.extra_bits_of_parity > 0) &&
                                  (h_matrix.element[i][j] >= 0)) ||
                                 ((h_matrix.extra_bits_of_parity == 0) &&
                                  (h_matrix.occupied[i][j] == 1));
      if (should_insert) {
        e = mod2sparse_insert(qc_bm, i, j);
        e->shift = h_matrix.element[i][j];
        total_cir++;
        row_wt[static_cast<size_t>(i)]++;
        col_wt[static_cast<size_t>(j)]++;
      }
    }
  }
  int max_row_wt = 0;
  for (int i = 0; i < bm_m; i++)
    max_row_wt = std::max(max_row_wt, row_wt[static_cast<size_t>(i)]);
  int max_col_wt = 0;
  for (int j = 0; j < bm_n; j++)
    max_col_wt = std::max(max_col_wt, col_wt[static_cast<size_t>(j)]);

  printf("[LDPC] H matrix porting ready! Total circulants: %d\n", total_cir);
  printf("[LDPC] Base-matrix max row weight: %d, max col weight: %d\n",
         max_row_wt, max_col_wt);
}

void ldpc_packet::print_hm() {
  FILE *fp, *fp1;
  FILE *fp_bm;
  mod2entry *e, *e_pre;
  int tmp;
  int tmp_val;
  uint32_t lr_sch;
  char ENS1[50];
  char ENS2[50];
  char BMC[50];
  char BMR[50];
  char BFS[50];
  char LRS[50];
  char BMS[60];

  snprintf(BMC, sizeof(BMC), "./output/HM_COL_order_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  snprintf(BMR, sizeof(BMR), "./output/HM_ROW_order_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  snprintf(ENS1, sizeof(ENS1), "./output/enc1_sched_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  snprintf(ENS2, sizeof(ENS2), "./output/enc2_sched_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  snprintf(BFS, sizeof(BFS), "./output/fdec_sched_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  snprintf(LRS, sizeof(LRS), "./output/rdec_sched_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  snprintf(BMS, sizeof(BMS), "./output/bm_schematic_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  
  // NOTE:
  // - For `cir_sz == 256`, encoder/FDEC schedules are available (qc_fi is
  // constructed in `ldpc_gen_gm()`).
  // - For `cir_sz == 512` (IBEX), these matrices/schedules are not constructed
  // in this code path; dumping them here can
  //   overflow buffers or dereference null pointers. We therefore only dump the
  //   layer-decoder (RDEC) scheduler.
  if (cir_sz == 256) {
    printf("[LDPC] Dump encoder scheduler 1 to file %s\n", ENS1);
    int *enc_sch = (int *)calloc(2 * col_wt, sizeof(*enc_sch));
    unsigned int **enc_fi;
    enc_fi = (unsigned int **)calloc((hm_m - tm_sz), sizeof(*enc_fi));
    for (int i = 0; i < (bm_m - tm_sz); i++)
      enc_fi[i] =
          (unsigned int *)calloc(cir_sz * (bm_m - tm_sz), sizeof(*enc_fi[i]));

    fp = fopen(ENS1, "w");
    fp1 = fopen(ENS2, "w");

    for (int i = 0; i < bm_n; i++) {
      for (int j = 0; j < col_wt; j++)
        enc_sch[j] = 0xFF;
      for (int j = col_wt; j < 2 * col_wt; j++)
        enc_sch[j] = 0x1F;

      tmp = 0;

      for (e = mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e)) {
        if (tmp >= col_wt) {
          printf("[LDPC Warning] ENC-scheduler: column %d has > col_wt(%d) "
                 "entries; extra entries are dropped.\n",
                 i, col_wt);
          break;
        }
        if (i < (bm_n - bm_m)) {
          enc_sch[tmp] = e->shift;
          enc_sch[tmp + col_wt] = e->row;
          tmp++;
        } else if (i < (bm_n - tm_sz)) {
          if (e->row < tm_sz) {
            enc_sch[tmp] = e->shift;
            enc_sch[tmp + col_wt] = e->row;
            tmp++;
          }
        } else {
          if (e->row >= tm_sz) {
            enc_sch[tmp] = e->shift;
            enc_sch[tmp + col_wt] = e->row;
            tmp++;
          }
        }
      }

      tmp = 0;
      for (int j = 0; j < col_wt; j++)
        tmp = (tmp << 5) + enc_sch[j + col_wt];

      if (tmp < pow(2, 4))
        fprintf(fp, "0000000%1X", tmp);
      else if (tmp < pow(2, 8))
        fprintf(fp, "000000%2X", tmp);
      else if (tmp < pow(2, 12))
        fprintf(fp, "00000%3X", tmp);
      else if (tmp < pow(2, 16))
        fprintf(fp, "0000%4X", tmp);
      else if (tmp < pow(2, 20))
        fprintf(fp, "000%5X", tmp);
      else if (tmp < pow(2, 24))
        fprintf(fp, "00%6X", tmp);
      else if (tmp < pow(2, 28))
        fprintf(fp, "0%7X", tmp);
      else
        fprintf(fp, "%8X", tmp);

      for (int j = 0; j < col_wt; j++) {
        if (enc_sch[j] < pow(2, 4))
          fprintf(fp, "0%1X", enc_sch[j]);
        else
          fprintf(fp, "%2X", enc_sch[j]);
      }

      fprintf(fp, "\n");
    }

    printf("[LDPC] Dump encoder scheduler 2 to file %s\n", ENS2);

    for (int i = 0; i < (bm_m - tm_sz); i++) {
      for (e = mod2sparse_first_in_col(qc_fi, i * cir_sz);
           !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e)) {
        enc_fi[i][e->row] = 1;
      }
    }

    unsigned int fi_tmp;
    for (int i = 0; i < (bm_m - tm_sz); i++) {
      for (int j = 0; j < (bm_m - tm_sz); j++) {
        for (int k = 0; k < (cir_sz / 32); k++) {
          fi_tmp = 0;
          for (int l = 0; l < 32; l++)
            fi_tmp += enc_fi[i][j * cir_sz + k * 32 + l] << l;

          if (fi_tmp < pow(2, 4))
            fprintf(fp1, "0000000%1X", fi_tmp);
          else if (fi_tmp < pow(2, 8))
            fprintf(fp1, "000000%2X", fi_tmp);
          else if (fi_tmp < pow(2, 12))
            fprintf(fp1, "00000%3X", fi_tmp);
          else if (fi_tmp < pow(2, 16))
            fprintf(fp1, "0000%4X", fi_tmp);
          else if (fi_tmp < pow(2, 20))
            fprintf(fp1, "000%5X", fi_tmp);
          else if (fi_tmp < pow(2, 24))
            fprintf(fp1, "00%6X", fi_tmp);
          else if (fi_tmp < pow(2, 28))
            fprintf(fp1, "0%7X", fi_tmp);
          else
            fprintf(fp1, "%8X", fi_tmp);

          fprintf(fp1, "\n");
        }
      }
    }

    fclose(fp);
    fclose(fp1);
    free(enc_sch);
    for (int i = 0; i < (hm_m - tm_sz); i++)
      free(enc_fi[i]);
    free(enc_fi);

    printf("[LDPC] Dump fast decoder scheduler to file %s\n", BFS);
    int *bf_sch = (int *)calloc(6, sizeof(*bf_sch));

    fp = fopen(BMC, "w");
    fp1 = fopen(BFS, "w");

    for (int i = 0; i < bm_n; i++) {
      fprintf(fp, "COL %3d: ", i);

      tmp = 0;
      for (int j = 0; j < 6; j++)
        bf_sch[j] = 0;

      for (e = mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e)) {
        fprintf(fp, "%2d(%3d) ", e->row, e->shift);
        tmp++;

        for (int j = 0; j < 6; j++)
          bf_sch[j] = (bf_sch[j] << 13);
        bf_sch[0] += (e->row + (e->shift << 5));
        for (int j = 0; j < 5; j++) {
          if (bf_sch[j] >= (1 >> 16)) {
            tmp_val = bf_sch[j] >> (16);
            bf_sch[j] -= tmp_val << (16);
            bf_sch[j + 1] += tmp_val;
          }
        }
      }

      for (int j = 0; j < col_wt; j++) {
        for (int j = 0; j < 6; j++)
          bf_sch[j] = (bf_sch[j] << 13);
        bf_sch[0] += 0x1FFF;
        for (int j = 0; j < 5; j++) {
          if (bf_sch[j] >= (1 >> 16)) {
            tmp_val = bf_sch[j] >> (16);
            bf_sch[j] -= tmp_val << (16);
            bf_sch[j + 1] += tmp_val;
          }
        }
      }

      fprintf(fp, "(col_wt = %d)\n", tmp);
      for (int j = 4; j >= 0; j--) {
        if (bf_sch[j] >= pow(2, 12))
          fprintf(fp1, "%4X", bf_sch[j]);
        else if (bf_sch[j] >= pow(2, 8))
          fprintf(fp1, "0%3X", bf_sch[j]);
        else if (bf_sch[j] >= pow(2, 4))
          fprintf(fp1, "00%2X", bf_sch[j]);
        else
          fprintf(fp1, "000%1X", bf_sch[j]);
      }
      fprintf(fp1, "\n");
    }

    fclose(fp);
    fclose(fp1);
    free(bf_sch);
  } else {
    printf("[LDPC] Skip ENC/FDEC scheduler dump for cir_sz=%d (IBEX mode)\n",
           cir_sz);
  }

  // Layer decoder scheduler
  printf("[LDPC] Dump layer decoder scheduler to file %s\n", LRS);
  printf("[LDPC] Dump base-matrix schematic to file %s\n", BMS);

  fp_bm = fopen(BMS, "w");
  if (fp_bm != NULL) {
    const int payload_cols_total = bm_n - bm_m;
    const int base_userdata_cols = (cir_sz == 512) ? 64 : payload_cols_total;
    int base_payload_cols = payload_cols_total;
    if (base_payload_cols > base_userdata_cols)
      base_payload_cols = base_userdata_cols;

    fprintf(fp_bm, "# Base-matrix schematic (bm_m=%d, bm_n=%d, Z=%d)\n", bm_m,
            bm_n, cir_sz);
    fprintf(fp_bm, "# Legend: 1=non-zero CPM, 0=zero CPM, X=zero CPM in "
                   "extra-userdata col (placeholder/skip)\n");
    fprintf(fp_bm,
            "# Groups: base_userdata_cols|extra_userdata_cols|parity_cols\n");
    fprintf(fp_bm,
            "# base_userdata_cols=%d, extra_userdata_cols=%d, parity_cols=%d\n",
            base_payload_cols, payload_cols_total - base_payload_cols, bm_m);
    fprintf(fp_bm,
            "# NOTE: extra_bits_of_userdata=%d (lane-level skip in last "
            "payload column, not shown per-CPM)\n\n",
            h_matrix.extra_bits_of_userdata);

    for (int r = 0; r < bm_m; r++) {
      fprintf(fp_bm, "ROW %2d: ", r);

      // base user-data columns
      for (int c = 0; c < base_payload_cols; c++) {
        const bool nz = (h_matrix.extra_bits_of_parity > 0)
                            ? (h_matrix.element[r][c] >= 0)
                            : (h_matrix.occupied[r][c] == 1);
        fputc(nz ? '1' : '0', fp_bm);
      }
      fputc('|', fp_bm);

      // extra user-data columns
      for (int c = base_payload_cols; c < payload_cols_total; c++) {
        const bool nz = (h_matrix.extra_bits_of_parity > 0)
                            ? (h_matrix.element[r][c] >= 0)
                            : (h_matrix.occupied[r][c] == 1);
        fputc(nz ? '1' : 'X', fp_bm);
      }
      fputc('|', fp_bm);

      // parity columns
      for (int c = payload_cols_total; c < bm_n; c++) {
        const bool nz = (h_matrix.extra_bits_of_parity > 0)
                            ? (h_matrix.element[r][c] >= 0)
                            : (h_matrix.occupied[r][c] == 1);
        fputc(nz ? '1' : '0', fp_bm);
      }
      fputc('\n', fp_bm);
    }
    fclose(fp_bm);
  } else {
    printf("[LDPC Warning] Failed to open %s for write\n", BMS);
  }

  int *sch_col;
  int last_in_row;
  int row_wt, row_dis;
  fp = fopen(BMR, "w");
  fp1 = fopen(LRS, "w");

  int cir_cnt = 0;
  int sch_out_cnt = 0;
  // Extra user-data columns (beyond the first 64 payload columns) are streamed
  // out even if the base-matrix entry is zero. For a zero circulant, we output
  // a 44-bit all-ones word as a placeholder marker. NOTE: This is only
  // meaningful for Z=512 (64 bytes per payload column).
  const int base_userdata_cols = 64;
  const int payload_cols_total = bm_n - bm_m;
  const int extra_userdata_col_start = base_userdata_cols;
  const int extra_userdata_col_end = payload_cols_total; // exclusive
  const int extra_userdata_col_cnt =
      (cir_sz == 512)
          ? std::max(0, extra_userdata_col_end - extra_userdata_col_start)
          : 0;
  int qc_bm_nnz = 0;
  for (int r = 0; r < bm_m; r++)
    for (e = mod2sparse_first_in_row(qc_bm, r); !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e))
      qc_bm_nnz++;
  if (qc_bm_nnz <= 0) {
    qc_bm_nnz = bm_n * col_wt;
    printf(
        "[LDPC Warning] qc_bm_nnz not detected, fallback to bm_n*col_wt=%d\n",
        qc_bm_nnz);
  }
  sch_col = (int *)calloc(qc_bm_nnz, sizeof(*sch_col));

#if 0
  // ---------------------------------------------------------------------------
  // Original RDEC scheduler export (kept for reference).
  //
  // Row-internal order:
  // - Group by pre_row distance via `for (j=1; j<bm_m; j++)`.
  // - Within each group, scan entries in `col` order.
  //
  // extra-userdata:
  // - Missing extra-userdata columns are emitted as trailing placeholders
  //   (`44'hFFFFFFFFFFF`) after all real circulants in the row.
  // ---------------------------------------------------------------------------
  for (int i=0; i<bm_m; i++)
  {
    fprintf(fp, "ROW %3d: ", i);

    tmp=0;
    row_wt = 0;
    last_in_row = 0;

    std::vector<uint8_t> extra_userdata_present;
    if (extra_userdata_col_cnt > 0)
      extra_userdata_present.assign(static_cast<size_t>(extra_userdata_col_cnt), 0);
    for (e = mod2sparse_first_in_row(qc_bm, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e)) {
      row_wt++;
      if ((extra_userdata_col_cnt > 0) && (e->col >= extra_userdata_col_start) && (e->col < extra_userdata_col_end)) {
        extra_userdata_present[static_cast<size_t>(e->col - extra_userdata_col_start)] = 1;
      }
    }

    for (int j=1; j<bm_m; j++)
    {
      for (e = mod2sparse_first_in_row(qc_bm, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
      {
        // search the previous circulant in the column
        e_pre = mod2sparse_prev_in_col(e);
        if (mod2sparse_at_end(e_pre))
          e_pre = mod2sparse_last_in_col(qc_bm, e->col);

        if (e_pre->row == (i+j)%bm_m) // non overlapped
        {
          fprintf(fp, "%3d(%3d/%2d) ", e->col, e->shift, e_pre->row);
          tmp++;

          // column log
          if (cir_cnt < qc_bm_nnz) {
            sch_col[cir_cnt] = e->col;
            cir_cnt++;
          } else {
            printf("[LDPC Error] RDEC-scheduler overflow: cir_cnt=%d >= qc_bm_nnz=%d\n", cir_cnt, qc_bm_nnz);
            break;
          }

          // delta shift
          tmp_val = (e->shift - e_pre->shift + cir_sz) % cir_sz;

          // last in the row (1 cycle in advance)
          if (tmp == (row_wt-1))
            last_in_row = 1;
          else
            last_in_row = 0;

          // RDEC scheduler packing (32-bit)
          // [6:0]   col (7)
          // [15:7]  shift (9)
          // [19:16] pre_cir_row (4)
          // [28:20] shift_delta (9)
          // [29]    last_in_row (1)
          // [30]    flag_64_extra_userdata (1): set when the circulant is in extra user-data column(s)
          // [31]    mask_flag (1): this circulant requires lane mask (MASK or INVMASK)
          // [40:32] mask_shift (9): signed shift of the last-row circulant in the same column (two's complement, [-256,255])

          // Defensive overflow checks (printing only; scheduler still truncates by design)
          if ((e->col >= (1 << 7)) || (e->shift >= (1 << 9)) || (e_pre->row >= (1 << 4)) || (tmp_val >= (1 << 9))) {
            printf("[LDPC Warning] RDEC-scheduler field overflow: row=%d col=%d shift=%d pre_row=%d shift_delta=%d\n", i, e->col,
                   e->shift, e_pre->row, tmp_val);
          }

          // Determine mask_flag (MASK/INVMASK) and per-column last-row shift for hardware.
          uint64_t mask_flag = 0;
          int signed_mask_shift = 0;
          // flag_64_extra_userdata is repurposed: mark entries in the extra user-data column(s).
          // Extra user-data columns are payload columns beyond the first 64 columns: [64, payload_cols_total).
          const uint64_t flag_64_extra_userdata =
              ((extra_userdata_col_cnt > 0) && (e->col >= extra_userdata_col_start) && (e->col < extra_userdata_col_end))
                  ? 1u : 0u;
          if ((cir_sz == 512) && (h_matrix.bits == 512) && (h_matrix.extra_bits_of_parity > 0) && (e->col >= 0) &&
              (e->col < h_matrix.cols) && (i >= 0) && (i < h_matrix.rows)) {
            if (h_matrix.fade[i][e->col])
              mask_flag = 1;
            else if (h_matrix.occupied[i][e->col] && (i == (h_matrix.rows - 1)))
              mask_flag = 1;

            const int last_row = h_matrix.rows - 1;
            const int last_row_shift = h_matrix.element[last_row][e->col];
            if (last_row_shift >= 0) {
              // Convert 0..511 to signed [-256,255] for RTL convenience.
              signed_mask_shift = (last_row_shift >= (cir_sz / 2)) ? (last_row_shift - cir_sz) : last_row_shift;
            }
          }

          // Pack a scheduler word with total 41 bits:
          // - [31:0] legacy fields + flags
          // - [40:32] signed mask_shift
          // We output it as 11 hex digits (44 bits) with the top 3 bits always 0.
          uint64_t sch64 = 0;
          sch64 |= (uint64_t(e->col) & 0x7Fu);
          sch64 |= (uint64_t(e->shift) & 0x1FFu) << 7;
          sch64 |= (uint64_t(e_pre->row) & 0x0Fu) << 16;
          sch64 |= (uint64_t(tmp_val) & 0x1FFu) << 20;
          sch64 |= (uint64_t(last_in_row) & 0x1u) << 29;
          sch64 |= (flag_64_extra_userdata & 0x1u) << 30;
          sch64 |= (mask_flag & 0x1u) << 31;
          sch64 |= (uint64_t(uint32_t(signed_mask_shift) & 0x1FFu)) << 32;

	          const unsigned long long mmem_word = (unsigned long long)(sch64 & ((1ull << 41) - 1));
	          const int sch_idx = sch_out_cnt++;
	          if (sch_idx >= (1 << 10)) {
	            printf("[LDPC Warning] RDEC-scheduler address overflow: sch_idx=%d (needs >10 bits)\n", sch_idx);
	          }
	          fprintf(fp1, "10'd%-6d:mmem_rdt=44'h%011llX;\n", sch_idx, mmem_word);
	        }
	      }

        if (j==bm_m-2)
        {
          row_dis = tmp;

          if (row_dis <= rdec_cmem_cont_thrshd)
            printf("[LDPC Error] RDEC-scheduler violation of C-MEM constraint!\n");
        }
	    }

	    // Export placeholders for zero-circulant entries in extra user-data columns.
	    // For each row, ensure we stream out one entry per extra user-data payload column, even if it is a zero circulant.
	    if (extra_userdata_col_cnt > 0) {
	      for (int c = extra_userdata_col_start; c < extra_userdata_col_end; c++) {
	        if (extra_userdata_present[static_cast<size_t>(c - extra_userdata_col_start)] == 0) {
	          const int sch_idx = sch_out_cnt++;
	          const unsigned long long mmem_word = 0xFFFFFFFFFFFULL; // 44-bit all ones
	          if (sch_idx >= (1 << 10)) {
	            printf("[LDPC Warning] RDEC-scheduler address overflow: sch_idx=%d (needs >10 bits)\n", sch_idx);
	          }
	          fprintf(fp1, "10'd%-6d:mmem_rdt=44'h%011llX;\n", sch_idx, mmem_word);
	        }
	      }
		    }

		    fprintf(fp, "(row_wt = %d, row distance = %d)\n", row_wt, row_dis);
		  }
#endif

  // ---------------------------------------------------------------------------
  // Modified RDEC scheduler export (P1 + contiguous extra-userdata block).
  //
  // Goal:
  // - Keep the original output order of all NON extra-userdata circulants.
  // - For extra-userdata payload columns `col ∈ [extra_userdata_col_start,
  // extra_userdata_col_end)`:
  //   emit a contiguous block in `col` ascending order, filling missing columns
  //   with the placeholder marker `44'hFFFFFFFFFFF`.
  // - Insert that block at P1: the location where the first extra-userdata
  // circulant
  //   would have appeared in the original per-row schedule.
  //
  // Notes:
  // - Placeholders are NOT counted as real circulants for HD-MEM/C-MEM checks
  // (conservative).
  // - `last_in_row` semantics remain tied to the number of real circulants
  // (`row_wt`).
  // ---------------------------------------------------------------------------
  for (int i = 0; i < bm_m; i++) {
    fprintf(fp, "ROW %3d: ", i);

    tmp = 0;
    row_wt = 0;
    last_in_row = 0;
    row_dis = -1;

    const int cmem_hazard_pre_row = (i + bm_m - 1) % bm_m; // i-1 (mod bm_m)

    // Map each extra-userdata payload column to its (optional) non-zero entry
    // in this row.
    std::vector<mod2entry *> extra_entry_by_col;
    if (extra_userdata_col_cnt > 0)
      extra_entry_by_col.assign(static_cast<size_t>(extra_userdata_col_cnt),
                                (mod2entry *)0);

    for (e = mod2sparse_first_in_row(qc_bm, i); !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e)) {
      row_wt++;
      if ((extra_userdata_col_cnt > 0) &&
          (e->col >= extra_userdata_col_start) &&
          (e->col < extra_userdata_col_end)) {
        extra_entry_by_col[static_cast<size_t>(e->col -
                                               extra_userdata_col_start)] = e;
      }
    }

    bool extra_block_emitted = false;

// Emit one REAL (non-placeholder) circulant entry to the row dump + scheduler
// ROM. `tmp` counts real emitted circulants in this row.
#define RDEC_EMIT_REAL_ENTRY(EE)                                               \
  do {                                                                         \
    e_pre = mod2sparse_prev_in_col((EE));                                      \
    if (mod2sparse_at_end(e_pre))                                              \
      e_pre = mod2sparse_last_in_col(qc_bm, (EE)->col);                        \
                                                                               \
    if ((row_dis < 0) && (e_pre->row == cmem_hazard_pre_row))                  \
      row_dis = tmp;                                                           \
                                                                               \
    fprintf(fp, "%3d(%3d/%2d) ", (EE)->col, (EE)->shift, e_pre->row);          \
    tmp++;                                                                     \
                                                                               \
    /* column log (used for HD-MEM constraint check only) */                   \
    if (cir_cnt < qc_bm_nnz) {                                                 \
      sch_col[cir_cnt] = (EE)->col;                                            \
      cir_cnt++;                                                               \
    } else {                                                                   \
      printf("[LDPC Error] RDEC-scheduler overflow: cir_cnt=%d >= "            \
             "qc_bm_nnz=%d\n",                                                 \
             cir_cnt, qc_bm_nnz);                                              \
    }                                                                          \
                                                                               \
    /* delta shift */                                                          \
    tmp_val = ((EE)->shift - e_pre->shift + cir_sz) % cir_sz;                  \
                                                                               \
    /* last in the row (1 cycle in advance) */                                 \
    if (tmp == (row_wt - 1))                                                   \
      last_in_row = 1;                                                         \
    else                                                                       \
      last_in_row = 0;                                                         \
                                                                               \
    /* Defensive overflow checks (printing only; scheduler still truncates by  \
     * design) */                                                              \
    if (((EE)->col >= (1 << 7)) || ((EE)->shift >= (1 << 9)) ||                \
        (e_pre->row >= (1 << 4)) || (tmp_val >= (1 << 9))) {                   \
      printf("[LDPC Warning] RDEC-scheduler field overflow: row=%d col=%d "    \
             "shift=%d pre_row=%d shift_delta=%d\n",                           \
             i, (EE)->col, (EE)->shift, e_pre->row, tmp_val);                  \
    }                                                                          \
                                                                               \
    /* Determine mask_flag (MASK/INVMASK) and per-column last-row shift for    \
     * hardware. */                                                            \
    uint64_t mask_flag = 0;                                                    \
    int signed_mask_shift = 0;                                                 \
    const uint64_t flag_64_extra_userdata =                                    \
        ((extra_userdata_col_cnt > 0) &&                                       \
         ((EE)->col >= extra_userdata_col_start) &&                            \
         ((EE)->col < extra_userdata_col_end))                                 \
            ? 1u                                                               \
            : 0u;                                                              \
    if ((cir_sz == 512) && (h_matrix.bits == 512) &&                           \
        (h_matrix.extra_bits_of_parity > 0) && ((EE)->col >= 0) &&             \
        ((EE)->col < h_matrix.cols) && (i >= 0) && (i < h_matrix.rows)) {      \
      if (h_matrix.fade[i][(EE)->col])                                         \
        mask_flag = 1;                                                         \
      else if (h_matrix.occupied[i][(EE)->col] && (i == (h_matrix.rows - 1)))  \
        mask_flag = 1;                                                         \
                                                                               \
      const int last_row = h_matrix.rows - 1;                                  \
      const int last_row_shift = h_matrix.element[last_row][(EE)->col];        \
      if (last_row_shift >= 0) {                                               \
        /* Convert 0..511 to signed [-256,255] for RTL convenience. */         \
        signed_mask_shift = (last_row_shift >= (cir_sz / 2))                   \
                                ? (last_row_shift - cir_sz)                    \
                                : last_row_shift;                              \
      }                                                                        \
    }                                                                          \
                                                                               \
    /* Pack a scheduler word with total 41 bits (written as 44'hXXXXXXXXXXX).  \
     */                                                                        \
    uint64_t sch64 = 0;                                                        \
    sch64 |= (uint64_t((EE)->col) & 0x7Fu);                                    \
    sch64 |= (uint64_t((EE)->shift) & 0x1FFu) << 7;                            \
    sch64 |= (uint64_t(e_pre->row) & 0x0Fu) << 16;                             \
    sch64 |= (uint64_t(tmp_val) & 0x1FFu) << 20;                               \
    sch64 |= (uint64_t(last_in_row) & 0x1u) << 29;                             \
    sch64 |= (flag_64_extra_userdata & 0x1u) << 30;                            \
    sch64 |= (mask_flag & 0x1u) << 31;                                         \
    sch64 |= (uint64_t(uint32_t(signed_mask_shift) & 0x1FFu)) << 32;           \
                                                                               \
    const unsigned long long mmem_word =                                       \
        (unsigned long long)(sch64 & ((1ull << 41) - 1));                      \
    const int sch_idx = sch_out_cnt++;                                         \
    if (sch_idx >= (1 << 10)) {                                                \
      printf("[LDPC Warning] RDEC-scheduler address overflow: sch_idx=%d "     \
             "(needs >10 bits)\n",                                             \
             sch_idx);                                                         \
    }                                                                          \
    fprintf(fp1, "10'd%-6d:mmem_rdt=44'h%011llX;\n", sch_idx, mmem_word);      \
  } while (0)

    // Main schedule generation in original order, with a one-time insertion of
    // the extra block at P1.
    for (int j = 1; j < bm_m; j++) {
      for (e = mod2sparse_first_in_row(qc_bm, i); !mod2sparse_at_end(e);
           e = mod2sparse_next_in_row(e)) {
        // search the previous circulant in the column
        e_pre = mod2sparse_prev_in_col(e);
        if (mod2sparse_at_end(e_pre))
          e_pre = mod2sparse_last_in_col(qc_bm, e->col);

        if (e_pre->row == (i + j) % bm_m) { // non overlapped
          const bool is_extra = (extra_userdata_col_cnt > 0) &&
                                (e->col >= extra_userdata_col_start) &&
                                (e->col < extra_userdata_col_end);

          if (is_extra) {
            if (!extra_block_emitted) {
#ifdef _LDPC_RDEC_EXTRA_CONTIG_DBG
              fprintf(stderr,
                      "[LDPC DBG] RDEC extra-block emit: row=%d j=%d "
                      "first_extra_col=%d\n",
                      i, j, e->col);
#endif
              // Emit extra-userdata payload columns as a fixed-length
              // contiguous block (col ascending).
              for (int c = extra_userdata_col_start; c < extra_userdata_col_end;
                   c++) {
                mod2entry *e_extra = extra_entry_by_col[static_cast<size_t>(
                    c - extra_userdata_col_start)];
                if (e_extra) {
                  RDEC_EMIT_REAL_ENTRY(e_extra);
                } else {
                  const int sch_idx = sch_out_cnt++;
                  const unsigned long long mmem_word =
                      0xFFFFFFFFFFFULL; // 44-bit all ones
                  if (sch_idx >= (1 << 10)) {
                    printf("[LDPC Warning] RDEC-scheduler address overflow: "
                           "sch_idx=%d (needs >10 bits)\n",
                           sch_idx);
                  }
                  fprintf(fp1, "10'd%-6d:mmem_rdt=44'h%011llX;\n", sch_idx,
                          mmem_word);
                }
              }

              extra_block_emitted = true;
            }

            // Skip this extra entry in the original walk; it has already been
            // emitted in the extra block.
            continue;
          }

          // Non-extra entry: keep original order.
          RDEC_EMIT_REAL_ENTRY(e);
        }
      }
    }

    // If this row contains no extra circulants at all, fall back to emitting
    // the extra block at the row end (this matches the original behavior in
    // this corner-case).
    if (!extra_block_emitted && (extra_userdata_col_cnt > 0)) {
      for (int c = extra_userdata_col_start; c < extra_userdata_col_end; c++) {
        mod2entry *e_extra = extra_entry_by_col[static_cast<size_t>(
            c - extra_userdata_col_start)];
        if (e_extra) {
          RDEC_EMIT_REAL_ENTRY(e_extra);
        } else {
          const int sch_idx = sch_out_cnt++;
          const unsigned long long mmem_word =
              0xFFFFFFFFFFFULL; // 44-bit all ones
          if (sch_idx >= (1 << 10)) {
            printf("[LDPC Warning] RDEC-scheduler address overflow: sch_idx=%d "
                   "(needs >10 bits)\n",
                   sch_idx);
          }
          fprintf(fp1, "10'd%-6d:mmem_rdt=44'h%011llX;\n", sch_idx, mmem_word);
        }
      }
    }

    if (row_dis < 0)
      row_dis = row_wt;
    if (row_dis <= rdec_cmem_cont_thrshd)
      printf("[LDPC Error] RDEC-scheduler violation of C-MEM constraint!\n");

    fprintf(fp, "(row_wt = %d, row distance = %d)\n", row_wt, row_dis);

#undef RDEC_EMIT_REAL_ENTRY
  }

  // check HD-MEM constraint
  if (cir_cnt > 0) {
    for (int i = 0; i < cir_cnt; i++)
      for (int j = 1; j <= rdec_hdmem_cont_thrshd; j++)
        if (sch_col[i] == sch_col[(i + cir_cnt - j) % cir_cnt])
          printf(
              "[LDPC Error] RDEC-scheduler violation of HD-MEM constraint!\n");
  }

  free(sch_col);
  fclose(fp);
  fclose(fp1);
}

/*
void ldpc_packet::print_hm()
{
  FILE *fp, *fp1;
  mod2entry *e, *e_pre;
  int tmp;
  int tmp_val;
  unsigned int lr_sch;
  char ENS1[50];
  char ENS2[50];
  char BMC[50];
  char BMR[50];
  char BFS[50];
  char LRS[50];

  sprintf(BMC, "./output/HM_COL_order_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz,
col_wt); sprintf(BMR, "./output/HM_ROW_order_%dx%dex%d_w%d.txt", bm_m, bm_n,
cir_sz, col_wt); sprintf(ENS1, "./output/enc1_sched_%dx%dex%d_w%d.txt", bm_m,
bm_n, cir_sz, col_wt); sprintf(ENS2, "./output/enc2_sched_%dx%dex%d_w%d.txt",
bm_m, bm_n, cir_sz, col_wt); sprintf(BFS,
"./output/fdec_sched_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz, col_wt);
  sprintf(LRS, "./output/rdec_sched_%dx%dex%d_w%d.txt", bm_m, bm_n, cir_sz,
col_wt);

  printf("[LDPC] Dump encoder scheduler 1 to file %s\n", ENS1);
  int *enc_sch = (int *)calloc(2*col_wt, sizeof(*enc_sch));
  unsigned int **enc_fi;
  enc_fi = (unsigned int **)calloc((hm_m-tm_sz), sizeof(*enc_fi));
  for (int i = 0; i < (bm_m - tm_sz); i++)
    enc_fi[i] = (unsigned int *)calloc(cir_sz*(bm_m-tm_sz), sizeof(*enc_fi[i]));

  fp = fopen(ENS1, "w");
  fp1 = fopen(ENS2, "w");

  for (int i = 0; i < bm_n; i++)
  {
    for (int j=0; j<col_wt; j++)
      enc_sch[j] = 0xFF;
    for (int j=col_wt; j<2*col_wt; j++)
      enc_sch[j] = 0x1F;

    tmp = 0;

    for (e= mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e); e =
mod2sparse_next_in_col(e))
    {
      if (i<(bm_n-bm_m))
      {
        enc_sch[tmp] = e->shift;
        enc_sch[tmp+col_wt] = e->row;
        tmp++;
      } else if (i<(bm_n-tm_sz))
      {
        if (e->row<tm_sz)
        {
          enc_sch[tmp] = e->shift;
          enc_sch[tmp+col_wt] = e->row;
          tmp++;
        }
      } else {
        if (e->row>=tm_sz)
        {
          enc_sch[tmp] = e->shift;
          enc_sch[tmp+col_wt] = e->row;
          tmp++;
        }
      }
    }

    tmp=0;
    for (int j=0; j<col_wt; j++)
      tmp = (tmp<<5)+enc_sch[j+col_wt];

    if (tmp<pow(2, 4))
      fprintf(fp, "0000000%1X", tmp);
    else if (tmp<pow(2, 8))
      fprintf(fp, "000000%2X", tmp);
    else if (tmp<pow(2, 12))
      fprintf(fp, "00000%3X", tmp);
    else if (tmp<pow(2, 16))
      fprintf(fp, "0000%4X", tmp);
    else if (tmp<pow(2, 20))
      fprintf(fp, "000%5X", tmp);
    else if (tmp<pow(2, 24))
      fprintf(fp, "00%6X", tmp);
    else if (tmp<pow(2, 28))
      fprintf(fp, "0%7X", tmp);
    else
      fprintf(fp, "%8X", tmp);

    for (int j=0; j<col_wt; j++)
    {
      if (enc_sch[j]<pow(2, 4))
        fprintf(fp, "0%1X", enc_sch[j]);
      else
        fprintf(fp, "%2X", enc_sch[j]);
    }

    fprintf(fp, "\n");
  }

  printf("[LDPC] Dump encoder scheduler 2 to file %s\n", ENS2);

  for (int i = 0; i < (bm_m - tm_sz); i++)
  {
    for (e = mod2sparse_first_in_col(qc_fi, i*cir_sz); !mod2sparse_at_end(e); e
= mod2sparse_next_in_col(e))
    {
      enc_fi[i][e->row] = 1;
    }
  }

  unsigned int fi_tmp;
  for (int i=0; i<(bm_m-tm_sz); i++)
  {
    for (int j=0; j<(bm_m-tm_sz); j++)
    {
      for (int k=0; k<(cir_sz/32); k++)
      {
        fi_tmp = 0;
        for (int l=0; l<32; l++)
          fi_tmp += enc_fi[i][j*cir_sz + k*32 + l]<<l;

        if (fi_tmp<pow(2, 4))
          fprintf(fp1, "0000000%1X", fi_tmp);
        else if (fi_tmp<pow(2, 8))
          fprintf(fp1, "000000%2X", fi_tmp);
        else if (fi_tmp<pow(2, 12))
          fprintf(fp1, "00000%3X", fi_tmp);
        else if (fi_tmp<pow(2, 16))
          fprintf(fp1, "0000%4X", fi_tmp);
        else if (fi_tmp<pow(2, 20))
          fprintf(fp1, "000%5X", fi_tmp);
        else if (fi_tmp<pow(2, 24))
          fprintf(fp1, "00%6X", fi_tmp);
        else if (fi_tmp<pow(2, 28))
          fprintf(fp1, "0%7X", fi_tmp);
        else
          fprintf(fp1, "%8X", fi_tmp);

        fprintf(fp1, "\n");
      }
    }
  }

  fclose(fp);
  fclose(fp1);
  free(enc_sch);
  for (int i = 0; i < (hm_m - tm_sz); i++)
    free(enc_fi[i]);
  free(enc_fi);

  printf("[LDPC] Dump fast decoder scheduler to file %s\n", BFS);
  int *bf_sch = (int *)calloc(col_wt, sizeof(*bf_sch));

  fp = fopen(BMC, "w");
  fp1 = fopen(BFS, "w");

  for (int i = 0; i < bm_n; i++)
  {
    fprintf(fp, "COL %3d: ", i);

    tmp=0;
    for (int j=0; j<col_wt; j++)
      bf_sch[j] = 0;

    for (e= mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e); e =
mod2sparse_next_in_col(e))
    {
      fprintf(fp, "%2d(%3d) ", e->row, e->shift);
      tmp++;

      for (int j=0; j<6; j++)
        bf_sch[j] = (bf_sch[j]<<13);
      bf_sch[0] += (e->row+(e->shift<<5));
      for (int j=0; j<5; j++){
        if (bf_sch[j]>=(1>>16))
        {
          tmp_val = bf_sch[j]>>(16);
          bf_sch[j] -= tmp_val<<(16);
          bf_sch[j+1] += tmp_val;
        }
      }
    }

    for (int j=0; j<col_wt; j++)
    {
      for (int j=0; j<6; j++)
        bf_sch[j] = (bf_sch[j]<<13);
      bf_sch[0] += 0x1FFF;
      for (int j=0; j<5; j++){
        if (bf_sch[j]>=(1>>16))
        {
          tmp_val = bf_sch[j]>>(16);
          bf_sch[j] -= tmp_val<<(16);
          bf_sch[j+1] += tmp_val;
        }
      }
    }

    fprintf(fp, "(col_wt = %d)\n", tmp);
    for (int j=4; j>=0; j--)
    {
      if (bf_sch[j]>=pow(2, 12))
        fprintf(fp1, "%4X", bf_sch[j]);
      else if (bf_sch[j]>=pow(2, 8))
        fprintf(fp1, "0%3X", bf_sch[j]);
      else if (bf_sch[j]>=pow(2, 4))
        fprintf(fp1, "00%2X", bf_sch[j]);
      else
        fprintf(fp1, "000%1X", bf_sch[j]);
    }
    fprintf(fp1, "\n");
  }

  fclose(fp);
  fclose(fp1);
  free(bf_sch);

  // Layer decoder scheduler
  printf("[LDPC] Dump layer decoder scheduler to file %s\n", LRS);

  int *sch_col;
  int last_in_row;
  int row_wt, row_dis;
  fp = fopen(BMR, "w");
  fp1 = fopen(LRS, "w");

  int cir_cnt = 0;
  sch_col = (int *)calloc(bm_n*col_wt-1, sizeof(*sch_col));

  for (int i=0; i<bm_m; i++)
  {
    fprintf(fp, "ROW %3d: ", i);

    tmp=0;
    row_wt = 0;
    last_in_row = 0;

    for (e = mod2sparse_first_in_row(qc_bm, i); !mod2sparse_at_end(e); e =
mod2sparse_next_in_row(e)) row_wt++;

    for (int j=1; j<bm_m; j++)
    {
      for (e = mod2sparse_first_in_row(qc_bm, i); !mod2sparse_at_end(e); e =
mod2sparse_next_in_row(e))
      {
        // search the previous circulant in the column
        e_pre = mod2sparse_prev_in_col(e);
        if (mod2sparse_at_end(e_pre))
          e_pre = mod2sparse_last_in_col(qc_bm, e->col);

        if (e_pre->row == (i+j)%bm_m) // non overlapped
        {
          fprintf(fp, "%3d(%3d/%2d) ", e->col, e->shift, e_pre->row);
          tmp++;

          // column log
          sch_col[cir_cnt] = e->col;
          cir_cnt++;

          // delta shift
          tmp_val = (e->shift - e_pre->shift + cir_sz) % cir_sz;

          // last in the row (1 cycle in advance)
          if (tmp == (row_wt-1))
            last_in_row = 1;
          else
            last_in_row = 0;

          // {last_in_row, shift_delta[7:0], pre_cir_row[4:0], shift[7:0],
col_index[7:0]} lr_sch = last_in_row*pow(2, 29) + tmp_val*pow(2, 21) +
e_pre->row*pow(2, 16) + e->shift*pow(2, 8) + e->col;

          if (lr_sch < pow(2, 4))
            fprintf(fp1, "0000000%1X\n", lr_sch);
          else if (lr_sch < pow(2, 8))
            fprintf(fp1, "000000%2X\n", lr_sch);
          else if (lr_sch < pow(2, 12))
            fprintf(fp1, "00000%3X\n", lr_sch);
          else if (lr_sch < pow(2, 16))
            fprintf(fp1, "0000%4X\n", lr_sch);
          else if (lr_sch < pow(2, 20))
            fprintf(fp1, "000%5X\n", lr_sch);
          else if (lr_sch < pow(2, 24))
            fprintf(fp1, "00%6X\n", lr_sch);
          else if (lr_sch < pow(2, 28))
            fprintf(fp1, "0%7X\n", lr_sch);
          else
            fprintf(fp1, "%8X\n", lr_sch);
        }
      }

      if (j==bm_m-2)
      {
        row_dis = tmp;

        if (row_dis <= rdec_cmem_cont_thrshd)
          printf("[LDPC Error] RDEC-scheduler violation of C-MEM
constraint!\n");
      }
    }

    fprintf(fp, "(row_wt = %d, row distance = %d)\n", row_wt, row_dis);
  }

  // check HD-MEM constraint
  for (int i=0; i<(bm_n*col_wt-1); i++)
    for (int j=1; j<=rdec_hdmem_cont_thrshd; j++)
      if (sch_col[i] == sch_col[(i+cir_cnt-j)%cir_cnt])
        printf("[LDPC Error] RDEC-scheduler violation of HD-MEM constraint!\n");

  free(sch_col);
  fclose(fp);
  fclose(fp1);
}
*/

// Generate G matrices from H
void ldpc_packet::ldpc_gen_gm() {
  FILE *fp;
  mod2sparse *qc_ac, *qc_bd, *qc_te;
  mod2sparse *qc_exb, *qc_f;
  mod2dense *qc_f_d, *qc_fi_d;
  int *qc_a_cols, *qc_b_cols, *qc_t_cols;
  int *qc_a_rows, *qc_c_rows;

  qc_ac = mod2sparse_allocate(hm_m, hm_k);
  qc_bd = mod2sparse_allocate(hm_m, (bm_m - tm_sz) * cir_sz);
  qc_te = mod2sparse_allocate(hm_m, tm_sz * cir_sz);

  qc_a_cols = (int *)calloc((bm_n - bm_m) * cir_sz, sizeof(*qc_a_cols));
  qc_b_cols = (int *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*qc_b_cols));
  qc_t_cols = (int *)calloc(tm_sz * cir_sz, sizeof(*qc_t_cols));
  qc_a_rows = (int *)calloc(tm_sz * cir_sz, sizeof(*qc_a_rows));
  qc_c_rows = (int *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*qc_c_rows));

  printf("[LDPC] Generating A/B/C/D/E matrices from H ...\n");

  for (int i = 0; i < (bm_n - bm_m) * cir_sz; i++)
    qc_a_cols[i] = i;
  for (int i = 0; i < (bm_m - tm_sz) * cir_sz; i++)
    qc_b_cols[i] = i + (bm_n - bm_m) * cir_sz;
  for (int i = 0; i < tm_sz * cir_sz; i++)
    qc_t_cols[i] = i + (bm_n - tm_sz) * cir_sz;

  for (int i = 0; i < tm_sz * cir_sz; i++)
    qc_a_rows[i] = i;
  for (int i = 0; i < (bm_m - tm_sz) * cir_sz; i++)
    qc_c_rows[i] = i + tm_sz * cir_sz;

  // split column first
  mod2sparse_copycols(qc_hm, qc_ac, qc_a_cols);
  mod2sparse_copycols(qc_hm, qc_bd, qc_b_cols);
  mod2sparse_copycols(qc_hm, qc_te, qc_t_cols);

  // then split row
  mod2sparse_copyrows(qc_ac, qc_a, qc_a_rows);
  mod2sparse_copyrows(qc_ac, qc_c, qc_c_rows);
  mod2sparse_copyrows(qc_bd, qc_b, qc_a_rows);
  mod2sparse_copyrows(qc_bd, qc_d, qc_c_rows);
  mod2sparse_copyrows(qc_te, qc_e, qc_c_rows);
  printf("[LDPC] A/B/C/D/E matrices ready!\n");

  /*
  #ifdef _LDPC_DUMP
    char MH[50] = "H_matrix.txt";
    char MA[50] = "A_matrix.txt";
    char MB[50] = "B_matrix.txt";
    char MC[50] = "C_matrix.txt";
    char MD[50] = "D_matrix.txt";
    char ME[50] = "E_matrix.txt";

    printf("[LDPC] Dump H matrix to file %s\n", MH);
    fp = fopen(MH, "w");
    mod2sparse_print(fp, qc_hm);
    fclose(fp);
    printf("[LDPC] Dump A matrix to file %s\n", MA);
    fp = fopen(MA, "w");
    mod2sparse_print(fp, qc_a);
    fclose(fp);
    printf("[LDPC] Dump B matrix to file %s\n", MB);
    fp = fopen(MB, "w");
    mod2sparse_print(fp, qc_b);
    fclose(fp);
    printf("[LDPC] Dump C matrix to file %s\n", MC);
    fp = fopen(MC, "w");
    mod2sparse_print(fp, qc_c);
    fclose(fp);
    printf("[LDPC] Dump D matrix to file %s\n", MD);
    fp = fopen(MD, "w");
    mod2sparse_print(fp, qc_d);
    fclose(fp);
    printf("[LDPC] Dump E matrix to file %s\n", ME);
    fp = fopen(ME, "w");
    mod2sparse_print(fp, qc_e);
    fclose(fp);
  #endif
  */

  // generate inverse F matrix (F=E*B+D)
  printf("[LDPC] Generating inverse F matrix from H matrix ...\n");

  qc_exb =
      mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);
  qc_f = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);

  mod2sparse_multiply(qc_e, qc_b, qc_exb);
  mod2sparse_add(qc_exb, qc_d, qc_f);

  qc_f_d = mod2dense_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);
  qc_fi_d =
      mod2dense_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);

  mod2sparse_to_dense(qc_f, qc_f_d);
  mod2dense_invert(qc_f_d, qc_fi_d);
  mod2dense_to_sparse(qc_fi_d, qc_fi);

  printf("[LDPC] Inverse F matrix generated!\n");

#ifdef _LDPC_FI_DUMP
  char MFi[50] = "./output/Fi_matrix.txt";
  printf("[LDPC] Dump inverse F matrix to file %s\n", MFi);
  fp = fopen(MFi, "w");
  mod2sparse_print(fp, qc_fi);
  fclose(fp);
#endif

  // free up all internals
  mod2sparse_free(qc_ac);
  mod2sparse_free(qc_bd);
  mod2sparse_free(qc_te);
  mod2sparse_free(qc_exb);
  mod2sparse_free(qc_f);
  mod2dense_free(qc_f_d);
  mod2dense_free(qc_fi_d);
  free(qc_a_cols);
  free(qc_b_cols);
  free(qc_t_cols);
  free(qc_a_rows);
  free(qc_c_rows);
} // ldpc_gen_gm

void ldpc_packet::f_print_h_matrix(s_h_matrix h_matrix) {
  int i;
  int j;
  for (i = 0; i < h_matrix.rows; i++) {
    printf("%02d  ", h_matrix.row_weight[i]);
    for (j = 0; j < h_matrix.cols; j++) {
      if (h_matrix.bits > 256) {
        if (h_matrix.occupied[i][j])
          printf("%03d ", h_matrix.element[i][j]);
        else
          printf("... ");
      } else {
        if (h_matrix.occupied[i][j])
          printf("%02X ", h_matrix.element[i][j]);
        else
          printf(".. ");
      }
    }
    printf("\n");
  }
  printf("\n");

  for (i = 0; i < h_matrix.rows; i++) {
    printf("%02d  ", h_matrix.row_weight[i]);
    for (j = 0; j < h_matrix.cols; j++) {
      if (h_matrix.occupied[i][j])
        printf("x,");
      else if (h_matrix.fade[i][j])
        printf("y,");
      else
        printf(" ,");
    }
    printf("\n");
  }
  printf("\n");
  printf("OCCUPIED\n");

  for (i = 0; i < h_matrix.rows; i++) {
    printf("%02d  ", h_matrix.row_weight[i]);
    for (j = 0; j < h_matrix.cols; j++) {
      if (h_matrix.occupied[i][j])
        printf("1 ");
      else
        printf(". ");
    }
    printf("\n");
  }
  printf("\n");
  printf("FADE\n");

  for (i = 0; i < h_matrix.rows; i++) {
    printf("%02d  ", h_matrix.row_weight[i]);
    for (j = 0; j < h_matrix.cols; j++) {
      if (h_matrix.fade[i][j])
        printf("1 ");
      else
        printf(". ");
    }
    printf("\n");
  }
  printf("\n");
}

// Function: IBEX Hard Codeword print
void ldpc_packet::f_print_hard_codeword(s_hard_codeword data, int cols,
                                        int bits) {
  int j;
  int k;
  int b;
  for (j = 0; j < cols; j++) {
    printf("%3d  0x", j);
    for (k = bits - 1; k >= 0; k--) {
      if ((k % 4) == 3)
        b = data.c[j].b[k] ? 8 : 0;
      else if ((k % 4) == 2)
        b += data.c[j].b[k] ? 4 : 0;
      else if ((k % 4) == 1)
        b += data.c[j].b[k] ? 2 : 0;
      else
        b += data.c[j].b[k] ? 1 : 0;

      if ((k % 4) == 0) {
        if (b == 0)
          printf(".");
        else
          printf("%1x", b);
      }

      if ((k % 16) == 0)
        printf("_");
    }
    printf("\n");
  }
}

// Function: IBEX CN print
void ldpc_packet::f_print_check_nodes(s_check_nodes data, int rows, int bits) {
  int i;
  int k;
  int b;

  for (i = 0; i < rows; i++) {
    printf("%3d  0x", i);
    for (k = bits - 1; k >= 0; k--) {
      if ((k % 4) == 3)
        b = data.r[i].b[k] ? 8 : 0;
      else if ((k % 4) == 2)
        b += data.r[i].b[k] ? 4 : 0;
      else if ((k % 4) == 1)
        b += data.r[i].b[k] ? 2 : 0;
      else
        b += data.r[i].b[k] ? 1 : 0;

      if (((k % 4) == 0) && (b == 0))
        printf(".");
      else if ((k % 4) == 0)
        printf("%1x", b);
    }
    printf("\n");
  }
}

void ldpc_packet::f_print_check_nodes_shifted(s_check_nodes data,
                                              s_h_matrix h_matrix, int column) {
  int i;
  int k;
  int b;
  s_check_nodes data_shifted;

  for (i = 0; i < h_matrix.rows; i++) {
    if (h_matrix.occupied[i][column] || h_matrix.fade[i][column])
      for (k = 0; k < h_matrix.bits; k++)
        data_shifted.r[i].b[k] =
            data.r[i].b[(h_matrix.bits + k - h_matrix.element[i][column]) %
                        h_matrix.bits];
    else
      for (k = 0; k < h_matrix.bits; k++)
        data_shifted.r[i].b[k] = 0;
  }

  for (i = 0; i < h_matrix.rows; i++) {
    if (h_matrix.occupied[i][column] || h_matrix.fade[i][column]) {
      if (h_matrix.occupied[i][column])
        printf("C++ CN SHIFTED COL %2d ROW %2d OCCUPIED  0x", column, i);
      else if (h_matrix.fade[i][column])
        printf("C++ CN SHIFTED COL %2d ROW %2d FADED     0x", column, i);

      for (k = h_matrix.bits - 1; k >= 0; k--) {
        if ((k % 4) == 3)
          b = data_shifted.r[i].b[k] ? 8 : 0;
        else if ((k % 4) == 2)
          b += data_shifted.r[i].b[k] ? 4 : 0;
        else if ((k % 4) == 1)
          b += data_shifted.r[i].b[k] ? 2 : 0;
        else
          b += data_shifted.r[i].b[k] ? 1 : 0;

        if (((k % 4) == 0) && (b == 0))
          printf(".");
        else if ((k % 4) == 0)
          printf("%1x", b);
      }
      printf("\n");
    }
  }
}

void ldpc_packet::f_print_s_256_bits(s_256_bits s) {
  int k;
  int b;
  printf("### PRNG C  ");

  for (k = 255; k >= 0; k--) {
    if ((k % 4) == 3)
      b = s.b[k] ? 8 : 0;
    else if ((k % 4) == 2)
      b += s.b[k] ? 4 : 0;
    else if ((k % 4) == 1)
      b += s.b[k] ? 2 : 0;
    else
      b += s.b[k] ? 1 : 0;

    if ((k % 4) == 0)
      printf("%1x", b);
    if ((k % 16) == 0)
      printf("_");
  }
  printf("\n");
}

void ldpc_packet::f_print_s_512_bits(s_512_bits s) {
  int k;
  int b;
  printf("### PRNG C  ");

  for (k = 512; k >= 0; k--) {
    if ((k % 4) == 3)
      b = s.b[k] ? 8 : 0;
    else if ((k % 4) == 2)
      b += s.b[k] ? 4 : 0;
    else if ((k % 4) == 1)
      b += s.b[k] ? 2 : 0;
    else
      b += s.b[k] ? 1 : 0;

    if ((k % 4) == 0)
      printf("%1x", b);
    if ((k % 16) == 0)
      printf("_");
  }
  printf("\n");
}

s_256_bits ldpc_packet::f_256_bit_lfsr(s_256_bits data_in) {
  s_256_bits data_out;
  int i;
  bool bit0;
  bit0 = data_in.b[0];

  for (i = 0; i < 256; i++) {
    if (i == 255)
      data_out.b[i] = bit0;
    else if ((i == 253) || (i == 250) || (i == 245))
      data_out.b[i] = data_in.b[i + 1] ^ bit0;
    else
      data_out.b[i] = data_in.b[i + 1];
  }
  return data_out;
}

s_512_bits ldpc_packet::f_512_bit_lfsr(s_512_bits data_in) {
  s_512_bits data_out;
  int i;
  bool bit0;
  bit0 = data_in.b[0];

  for (i = 0; i < 512; i++) {
    if (i == 511)
      data_out.b[i] = bit0;
    else if ((i == 509) || (i == 506) || (i == 503))
      data_out.b[i] = data_in.b[i + 1] ^ bit0;
    else
      data_out.b[i] = data_in.b[i + 1];
  }
  return data_out;
}

// Function: IBEX Syndrome calculate
s_check_nodes ldpc_packet::f_check_nodes(s_h_matrix h_matrix,
                                         s_hard_codeword vn) {
  s_check_nodes cn;
  int i, j, k, m;

  // initialize check nodes
  for (i = 0; i < h_matrix.rows; i++)
    for (k = 0; k < h_matrix.bits; k++)
      cn.r[i].b[k] = 0;

  // calculate check nodes
  for (j = 0; j < h_matrix.cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      if (vn.c[j].b[k]) {
        for (i = 0; i < h_matrix.rows; i++) {
          m = (k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits;
          if (h_matrix.extra_bytes_of_parity == 0) {
            if (h_matrix.occupied[i][j])
              cn.r[i].b[m] = !cn.r[i].b[m];
          } else {
            if (h_matrix.occupied[i][j] && (i < h_matrix.rows - 1))
              cn.r[i].b[m] = !cn.r[i].b[m];
            if (h_matrix.occupied[i][j] && (i == h_matrix.rows - 1) &&
                h_matrix.mask[j][k])
              cn.r[i].b[m] = !cn.r[i].b[m];
            if (h_matrix.fade[i][j] && !h_matrix.mask[j][k])
              cn.r[i].b[m] = !cn.r[i].b[m];
          }
        }
      }
    }
  }
  return cn;
}

// Function: IBEX Syndrome weight calculate
int ldpc_packet::f_check_node_weight(s_h_matrix h_matrix, s_check_nodes cn) {
  int check_node_weight;
  int row, bit;
  check_node_weight = 0;
  for (row = 0; row < h_matrix.rows; row++) {
    for (bit = 0; bit < h_matrix.bits; bit++) {
      if (cn.r[row].b[bit])
        check_node_weight++;
    }
  }
  return check_node_weight;
}

// Function: IBEX Likelihood level Initial
s_likelihood_levels ldpc_packet::f_likelihood_levels(
    int strobes, s_ldpc_decoder_parameters ldpc_decoder_parameters,
    int syndrome_weight, int rows) {
  s_likelihood_levels likelihood_levels;
  int address;
  int delta[4];
  int delta_sum;
  int delta_total;
  int coef[4];
  int weak_minus_strong;
  int coef_index = 0;
  likelihood_levels.max = (1 << VN_BITS) - 1;
  likelihood_levels.min = 0;

  // VN 位宽自适应的翻转阈值与强弱档初始化
  if (VN_BITS <= 2) {
    likelihood_levels.flip_thr = 2; // 3
    likelihood_levels.weak = 0;     // 2
    likelihood_levels.strong = 0;   // 1
  } else if (VN_BITS == 3) {
    likelihood_levels.flip_thr = likelihood_levels.max - 3;
    likelihood_levels.weak = likelihood_levels.flip_thr - 4;
    likelihood_levels.strong = likelihood_levels.weak;
  } else {
    likelihood_levels.flip_thr =
        std::max(likelihood_levels.min + 1, likelihood_levels.max - 7);
    likelihood_levels.weak =
        std::min(likelihood_levels.min, likelihood_levels.flip_thr - 4);
    likelihood_levels.strong = likelihood_levels.weak;
  }

  likelihood_levels.level[0] = likelihood_levels.weak;
  likelihood_levels.level[1] = likelihood_levels.weak;
  likelihood_levels.level[2] = likelihood_levels.weak;
  likelihood_levels.level[3] = likelihood_levels.weak;
  if (strobes > 0) {
    address = syndrome_weight >> 5;
    if (address >= 64)
      address = 63;
    int coef_index = 0;
    if (rows == 7)
      coef_index = 1;
    else if (rows == 8)
      coef_index = 2;
    else if (rows == 9)
      coef_index = 3;
    else if (rows == 10)
      coef_index = 4;
    else if (rows == 11)
      coef_index = 5;
    else if (rows == 12)
      coef_index = 6;
    else if (rows >= 13)
      coef_index = 7;
    coef[0] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][0];
    coef[1] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][1];
    coef[2] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][2];
    coef[3] = ldpc_decoder_parameters.likelihood_init_coef_all[coef_index][3];
    delta[0] = (address >= 0) ? (coef[0] * (address - 0)) : 0;
    delta[1] = (address >= 8) ? (coef[1] * (address - 8)) : 0;
    delta[2] = (address >= 16) ? (coef[2] * (address - 16)) : 0;
    delta[3] = (address >= 24) ? (coef[3] * (address - 24)) : 0;
    delta_sum = delta[0] + delta[1] + delta[2] + delta[3];
    if (VN_BITS == 2) {
      delta_total = delta_sum >> 8;
      likelihood_levels.strong = likelihood_levels.min + delta_total;
      if (likelihood_levels.strong > likelihood_levels.weak)
        likelihood_levels.strong = likelihood_levels.weak;
    } else {
      delta_total = (VN_BITS == 8) ? ((delta_sum >> 1) + (delta_sum >> 2))
                                   : ((delta_sum >> 2) + (delta_sum >> 3));
      likelihood_levels.strong = likelihood_levels.weak - delta_total;
      if (likelihood_levels.strong < likelihood_levels.min)
        likelihood_levels.strong = likelihood_levels.min;
    }

    likelihood_levels.level[0] = likelihood_levels.strong;
    likelihood_levels.level[1] = likelihood_levels.level[0];
    likelihood_levels.level[2] = likelihood_levels.level[3];
  }
  if (strobes > 1) {
    likelihood_levels.level[1] = likelihood_levels.level[0];
    likelihood_levels.level[2] = likelihood_levels.level[3];
    weak_minus_strong = likelihood_levels.level[3] - likelihood_levels.level[0];
    int split_threshold = (VN_BITS <= 2) ? 1 : 8;
    if (weak_minus_strong >= split_threshold) {
      likelihood_levels.level[1] =
          likelihood_levels.level[3] -
          ((weak_minus_strong *
            ldpc_decoder_parameters.likelihood_init_fraction[1]) >>
           4);
      likelihood_levels.level[2] =
          likelihood_levels.level[3] -
          ((weak_minus_strong *
            ldpc_decoder_parameters.likelihood_init_fraction[2]) >>
           4);
    }
  }
  likelihood_levels.min = likelihood_levels.strong;
  return likelihood_levels;
}

int ldpc_packet::f_update_vn_post(int likelihood, int weight,
                                  int min_likelihood, int max_likelihood,
                                  bool post_process, bool post_process2,
                                  bool be_aggressive, int flip_threshold,
                                  bool pushing) {
  int likelihood_new;
  bool clamp_to_min;
  bool do_post_flipped;
  bool do_post_unflipped;
  bool do_aggr;
  bool flipped = (likelihood >= flip_threshold);

  int delta;
  if (VN_BITS <= 2)
    delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
  else
    delta = weight;

  likelihood_new =
      flipped ? (likelihood - delta)
              : (likelihood + delta - 1); // when clamp!=1, same as strong
  // likelihood_new = flipped ? (likelihood - weight) : (likelihood + weight -
  // 1); // when clamp!=1, same as strong

  bool flipped_new = (likelihood_new >= flip_threshold);
  do_post_flipped = post_process && (likelihood_new == flip_threshold);
  do_post_unflipped = post_process && (likelihood_new < flip_threshold);
  if (do_post_unflipped)
    likelihood_new = flip_threshold - 1;
  else if (do_post_flipped)
    likelihood_new = flip_threshold + 1;

  if (post_process2 && !flipped && (weight == 1) &&
      (likelihood_new == (flip_threshold - 1)))
    likelihood_new++;
  if (likelihood_new <= min_likelihood)
    likelihood_new = min_likelihood;
  if (likelihood_new >= max_likelihood)
    likelihood_new = max_likelihood;

  return likelihood_new;
}

// int ldpc_packet::f_update_vn_post(int likelihood, int weight, int
// min_likelihood, int max_likelihood, bool post_process,
//                                   bool post_process2, bool be_aggressive, int
//                                   flip_threshold, bool look, bool scale2x) {
//   const int scale = (scale2x && VN_BITS <= 2) ? 1 : 0;
//   const int unit = 1;

//   int like = (likelihood << scale) + ((scale) ? 1 : 0);
//   int thr;
//   if (scale) thr = 4; // 3-bit 域的经典阈值
//   else thr = flip_threshold;
//   int mn = min_likelihood << scale;
//   int mx = (max_likelihood << scale) + ((scale) ? 1 : 0);
//   int w_ext = weight;

//   bool flipped = (like >= thr);
//   int like_new = flipped ? like - w_ext : like + w_ext - unit;
//   bool do_post_flipped = post_process && (like_new == thr);
//   bool do_post_unflipped = post_process && (like_new < thr);
//   if (do_post_unflipped)
//     like_new = thr - unit;
//   else if (do_post_flipped)
//     like_new = thr + unit;

//   if (post_process2 && !flipped && (weight == 1) && (like_new == thr - unit))
//     like_new += unit;
//   if (like_new <= mn)
//     like_new = mn;
//   if (like_new >= mx)
//     like_new = mx;

//   return like_new >> scale;
// }

// LDPC code configuration
// void ldpc_packet::ldpc_config(int m, int n, int sc, int st, int wt) {
//   bm_m = m;
//   bm_n = n;
//   bm_k = n - m;
//   tm_sz = st;
//   cir_sz = sc;

//   col_wt = wt;
//   hm_m = m * sc;
//   hm_n = n * sc;
//   hm_k = hm_n - hm_m;
//   pad_len = hm_k - info_len;

//   printf(
//       "[LDPC] Configuring LDPC code with m=%d, n=%d, sc=%d, st=%d, col_wt=%d\n",
//       bm_m, bm_n, cir_sz, tm_sz, col_wt);

//   if (sc != 512) {
//     printf("[LDPC ERROR] Unsupported cir_sz=%d in MP_Framework. Please use sc=512.\n",
//            sc);
//     return;
//   }

//   {
//     int VERBOSITY = 0;
//     int i;
//     int j;
//     int k;
//     int bit;
//     int occupied_flag;
//     int fade_flag;
//     int n_full;
//     int h_matrix_index[LDPC_MAX_ROWS];
//     int init_base;
//     int shift_base;
//     int wrap_found;
//     int payload_cols;
//     int source_ok = 0;
//     char matrix_path[512] = "";
//     char occupied_path[512] = "";
//     char fade_path[512] = "";
//     char source_desc[512] = "";
//     FILE *fp_h = NULL;
//     FILE *fp_occupied = NULL;
//     FILE *fp_fade = NULL;
//     const char *env_root = getenv("DVC_IBEX_MATRIX_ROOT");
//     const char *default_roots[] = {
//         "IBEX/ibex_matrix_flat_13rate",
//         "../IBEX/ibex_matrix_flat_13rate",
//         "../../IBEX/ibex_matrix_flat_13rate",
//         "../../../IBEX/ibex_matrix_flat_13rate",
//         "../../../../IBEX/ibex_matrix_flat_13rate",
//     };

//     h_matrix.rows = m;
//     h_matrix.cols = n;
//     h_matrix.bits = sc;
//     h_matrix.bytes_of_userdata = info_len / 8;
//     h_matrix.bytes_of_parity = (blk_len - info_len) / 8;
//     h_matrix.min_rows = 5;
//     h_matrix.max_rows = LDPC_MAX_ROWS;

//     if ((h_matrix.rows <= 0) || (h_matrix.rows > LDPC_MAX_ROWS) ||
//         (h_matrix.cols <= 0) || (h_matrix.cols > LDPC_MAX_COLS)) {
//       printf("[LDPC ERROR] Unsupported matrix size bm_m=%d bm_n=%d (max rows=%d, max cols=%d)\n",
//              h_matrix.rows, h_matrix.cols, LDPC_MAX_ROWS, LDPC_MAX_COLS);
//       return;
//     }
//     if ((h_matrix.cols - h_matrix.rows) < 0 ||
//         (h_matrix.cols - h_matrix.rows) > LDPC_MAX_PAYLOAD_COLS) {
//       printf("[LDPC ERROR] Unsupported payload cols=%d (max supported=%d)\n",
//              h_matrix.cols - h_matrix.rows, LDPC_MAX_PAYLOAD_COLS);
//       return;
//     }

//     h_matrix.unused_bytes_of_parity =
//         (h_matrix.rows * (h_matrix.bits >> 3)) - h_matrix.bytes_of_parity;
//     if (h_matrix.unused_bytes_of_parity != 0)
//       h_matrix.extra_bytes_of_parity =
//           (h_matrix.bits >> 3) - h_matrix.unused_bytes_of_parity;
//     else
//       h_matrix.extra_bytes_of_parity = 0;

//     h_matrix.unused_bytes_of_userdata =
//         ((h_matrix.cols - h_matrix.rows) * (h_matrix.bits >> 3)) -
//         h_matrix.bytes_of_userdata;
//     if (h_matrix.unused_bytes_of_userdata != 0)
//       h_matrix.extra_bytes_of_userdata =
//           (h_matrix.bits >> 3) - h_matrix.unused_bytes_of_userdata;
//     else
//       h_matrix.extra_bytes_of_userdata = 0;

//     h_matrix.extra_bits_of_parity = h_matrix.extra_bytes_of_parity << 3;
//     h_matrix.extra_bits_of_userdata = h_matrix.extra_bytes_of_userdata << 3;

//     if (VERBOSITY > 0)
//       printf("### MATRIX: BITS: %4d ROWS: %2d COLS: %3d BYTES_OF_USERDATA: %5d "
//              "BYTES_OF_PARITY: %4d EXTRA_BYTES_OF_USERDATA: %5d EXTRA_BYTES_OF_PARITY: %4d\n",
//              h_matrix.bits, h_matrix.rows, h_matrix.cols,
//              h_matrix.bytes_of_userdata, h_matrix.bytes_of_parity,
//              h_matrix.extra_bytes_of_userdata, h_matrix.extra_bytes_of_parity);

//     n_full = LDPC_MAX_PAYLOAD_COLS + m;
//     if (n > n_full) {
//       printf("[LDPC ERROR] target n=%d > n_full=%d, cannot trim\n", n, n_full);
//       return;
//     }

//     for (i = 0; i < LDPC_MAX_ROWS; i++) {
//       h_matrix.delta[i] = 0;
//       h_matrix.first_element[i] = 0;
//       h_matrix.last_element[i] = 0;
//       h_matrix.wraparound[i] = 0;
//       h_matrix.wrap_base[i] = 0;
//       h_matrix.wrap_num_deltas[i] = 0;
//       h_matrix.row_weight[i] = 0;
//       h_matrix_index[i] = 0;
//       for (j = 0; j < LDPC_MAX_COLS; j++) {
//         h_matrix.element[i][j] = -1;
//         h_matrix.occupied[i][j] = 0;
//         h_matrix.fade[i][j] = 0;
//       }
//     }

//     for (j = 0; j < LDPC_MAX_COLS; j++) {
//       h_matrix.col_weight[j] = 0;
//       h_matrix.parity_column[j] = 0;
//       for (k = 0; k < h_matrix.bits; k++)
//         h_matrix.mask[j][k] = 0;
//     }

//     for (int root_index = -1;
//          root_index < (int)(sizeof(default_roots) / sizeof(default_roots[0]));
//          root_index++) {
//       const char *root =
//           (root_index < 0) ? env_root : default_roots[root_index];

//       if (!root || !root[0])
//         continue;

//       snprintf(matrix_path, sizeof(matrix_path),
//                "%s/matrix/LDPC_%dx%dex%d_w4_dense5_QC_H_1.txt", root,
//                h_matrix.rows, n_full, h_matrix.bits);
//       snprintf(occupied_path, sizeof(occupied_path),
//                "%s/occupied_matrix/LDPC_%dx%dex%d_w4_dense5_occupied_1.txt",
//                root, h_matrix.rows, n_full, h_matrix.bits);
//       snprintf(fade_path, sizeof(fade_path),
//                "%s/fade_matrix/LDPC_%dx%dex%d_w4_dense5_fade_1.txt", root,
//                h_matrix.rows, n_full, h_matrix.bits);

//       fp_h = fopen(matrix_path, "r");
//       fp_occupied = fopen(occupied_path, "r");
//       fp_fade = fopen(fade_path, "r");

//       if (!fp_h || !fp_occupied || !fp_fade) {
//         if (fp_h)
//           fclose(fp_h);
//         if (fp_occupied)
//           fclose(fp_occupied);
//         if (fp_fade)
//           fclose(fp_fade);
//         fp_h = NULL;
//         fp_occupied = NULL;
//         fp_fade = NULL;
//         continue;
//       }

//       snprintf(source_desc, sizeof(source_desc), "%s", matrix_path);
//       source_ok = 1;
//       break;
//     }

//     if (!source_ok) {
//       printf("[LDPC ERROR] Failed to load matrix/fade/occupied files for bm_m=%d bm_n=%d. "
//              "Set DVC_IBEX_MATRIX_ROOT or run from a tree containing IBEX/ibex_matrix_flat_13rate.\n",
//              h_matrix.rows, h_matrix.cols);
//       return;
//     }

//     h_matrix.delta[0] = 0;
//     h_matrix.delta[1] = 13;
//     h_matrix.delta[2] = 19;
//     h_matrix.delta[3] = 29;
//     h_matrix.delta[4] = 41;
//     h_matrix.delta[5] = 67;
//     h_matrix.delta[6] = 73;
//     h_matrix.delta[7] = 79;
//     h_matrix.delta[8] = 91;
//     h_matrix.delta[9] = 97;
//     h_matrix.delta[10] = 103;
//     h_matrix.delta[11] = 111;
//     h_matrix.delta[12] = 119;
//     h_matrix.delta[13] = 127;
//     h_matrix.delta[14] = 131;
//     h_matrix.delta[15] = 137;
//     h_matrix.delta[16] = 149;

//     h_matrix_index[0] = 0 + (3 * h_matrix.delta[0]);
//     h_matrix_index[1] = 0 + (4 * h_matrix.delta[1]);
//     h_matrix_index[2] = 0 + (4 * h_matrix.delta[2]);
//     h_matrix_index[3] = 0 + (4 * h_matrix.delta[3]);
//     h_matrix_index[4] = 0 + (3 * h_matrix.delta[4]);
//     for (i = 5; i < h_matrix.rows; i++)
//       h_matrix_index[i] = 0;

//     for (i = 0; i < h_matrix.rows; i++)
//       h_matrix.last_element[i] = h_matrix_index[i];

//     payload_cols = h_matrix.cols - h_matrix.rows;
//     for (i = 0; i < h_matrix.rows; i++) {
//       for (j = 0; j < n_full; j++) {
//         if (fscanf(fp_fade, "%d", &fade_flag) != 1)
//           source_ok = 0;
//         if (fscanf(fp_occupied, "%d", &occupied_flag) != 1)
//           source_ok = 0;

//         if ((j >= payload_cols) && (j < LDPC_MAX_PAYLOAD_COLS))
//           continue;

//         if ((j < payload_cols) && (occupied_flag == 1))
//           h_matrix.occupied[i][j] = 1;
//         if ((j < payload_cols) && (fade_flag == 1))
//           h_matrix.fade[i][j] = 1;
//         if ((j >= LDPC_MAX_PAYLOAD_COLS) && (occupied_flag == 1))
//           h_matrix.occupied[i][payload_cols + j - LDPC_MAX_PAYLOAD_COLS] = 1;
//         if ((j >= LDPC_MAX_PAYLOAD_COLS) && (fade_flag == 1))
//           h_matrix.fade[i][payload_cols + j - LDPC_MAX_PAYLOAD_COLS] = 1;
//       }
//     }

//     fclose(fp_h);
//     fclose(fp_occupied);
//     fclose(fp_fade);
//     fp_h = NULL;
//     fp_occupied = NULL;
//     fp_fade = NULL;

//     if (!source_ok) {
//       printf("[LDPC ERROR] Failed while parsing occupied/fade files for bm_m=%d bm_n=%d from %s\n",
//              h_matrix.rows, h_matrix.cols, source_desc);
//       return;
//     }

//     for (j = 0; j < h_matrix.cols; j++) {
//       const int column_index = h_matrix.cols - 1 - j;

//       k = column_index;
//       h_matrix.parity_column[j] = (j < h_matrix.rows);
//       h_matrix.bits_in_last_column = 8 * 20;
//       for (i = 0; i < h_matrix.rows; i++) {
//         if (h_matrix.occupied[i][column_index] || h_matrix.fade[i][column_index]) {
//           h_matrix.element[i][k] = h_matrix_index[i];
//           h_matrix.first_element[i] = h_matrix.element[i][k];
//           h_matrix_index[i] =
//               (h_matrix_index[i] + h_matrix.bits - h_matrix.delta[i]) %
//               h_matrix.bits;
//           if (h_matrix.occupied[i][column_index]) {
//             h_matrix.row_weight[i]++;
//             h_matrix.col_weight[k]++;
//           }
//         } else
//           h_matrix.element[i][k] = -1;
//       }
//     }

//     for (i = 0; i < h_matrix.rows; i++) {
//       h_matrix.wraparound[i] =
//           (h_matrix.bits + h_matrix.first_element[i] - h_matrix.last_element[i]) %
//           h_matrix.bits;
//       init_base =
//           (h_matrix.first_element[i] - h_matrix.delta[i] + h_matrix.bits) %
//           h_matrix.bits;
//       wrap_found = 0;
//       for (j = 0; j < h_matrix.bits; j++) {
//         shift_base = init_base + j;
//         for (k = 0; k < 50; k++) {
//           if ((h_matrix.last_element[i] + shift_base +
//                k * h_matrix.delta[i] + j) %
//                   h_matrix.bits ==
//               shift_base) {
//             h_matrix.wrap_num_deltas[i] = k;
//             h_matrix.wrap_base[i] = shift_base;
//             wrap_found = 1;
//             break;
//           }
//         }
//         if (wrap_found)
//           break;
//       }
//     }

//     h_matrix.bits_in_last_column = 8 * 20;

//     printf("[LDPC] Loaded external matrix from %s (source_cols=%d, crop_payload=%d)\n",
//            source_desc, n_full, h_matrix.cols - h_matrix.rows);

//     for (j = 0; j < h_matrix.cols; j++) {
//       if (h_matrix.occupied[h_matrix.rows - 1][j]) {
//         for (k = 0; k < h_matrix.bits; k++) {
//           bit = (k + h_matrix.bits -
//                  h_matrix.element[h_matrix.rows - 1][j]) %
//                 h_matrix.bits;
//           if ((bit < h_matrix.extra_bits_of_parity) ||
//               (h_matrix.extra_bits_of_parity == 0))
//             h_matrix.mask[j][k] = 1;
//         }
//       }
//     }

// #ifdef _LDPC_DEBUG_DUMP
//     printf("[DVC][LDEC_MATRIX] source=%s source_cols=%d rows=%d cols=%d hm_k=%d info_len=%d blk_len=%d extra_ud_bits=%d extra_pa_bits=%d\n",
//            source_desc, n_full, h_matrix.rows, h_matrix.cols, hm_k,
//            info_len, blk_len, h_matrix.extra_bits_of_userdata,
//            h_matrix.extra_bits_of_parity);
//     for (i = 0; i < h_matrix.rows; i++) {
//       printf("[DVC][LDEC_MATRIX] row=%02d FE=%03d LE=%03d WA=%03d RW=%02d\n", i,
//              h_matrix.first_element[i], h_matrix.last_element[i],
//              h_matrix.wraparound[i], h_matrix.row_weight[i]);
//     }
//     f_print_h_matrix(h_matrix);
// #endif

//     if (VERBOSITY > 0)
//       f_print_h_matrix(h_matrix);

//     ldpc_ibex_phck(h_matrix);
//   }

// #ifdef _LDPC_DUMP
//   rdec_cmem_cont_thrshd = 10;
//   rdec_hdmem_cont_thrshd = 6;

//   print_hm();
// #endif
// } // ldpc_config

void ldpc_packet::ldpc_config(int m, int n, int sc, int st, int wt) {
  bm_m = m;
  bm_n = n;
  bm_k = n - m;
  tm_sz = st;
  cir_sz = sc;

  col_wt = wt;
  hm_m = m * sc;
  hm_n = n * sc;
  hm_k = hm_n - hm_m;
  pad_len = hm_k - info_len;

  printf(
      "[LDPC] Configuring LDPC code with m=%d, n=%d, sc=%d, st=%d, col_wt=%d\n",
      bm_m, bm_n, cir_sz, tm_sz, col_wt);

  if (sc != 512) {
    printf("[LDPC ERROR] Unsupported cir_sz=%d in MP_Framework. Please use sc=512.\n",
           sc);
    return;
  }

  {
    int VERBOSITY = 0;
    int i;
    int j;
    int k;
    int bit;
    int col_shift;
    int occupied_flag;
    int fade_flag;
    int src_col;
    int n_full;
    int h_matrix_index[LDPC_MAX_ROWS];
    int row_nz_cnt[LDPC_MAX_ROWS];
    int full_shift[LDPC_MAX_ROWS][LDPC_MAX_COLS];
    int full_occupied[LDPC_MAX_ROWS][LDPC_MAX_COLS];
    int full_fade[LDPC_MAX_ROWS][LDPC_MAX_COLS];
    int source_cols = 0;
    int source_ok = 0;
    char matrix_path[512] = "";
    char occupied_path[512] = "";
    char fade_path[512] = "";
    char source_desc[512] = "";
    char source_occupied_desc[512] = "";
    char source_fade_desc[512] = "";
    char source_root_desc[512] = "";
    FILE *fp_h = NULL;
    FILE *fp_occupied = NULL;
    FILE *fp_fade = NULL;
    const char *env_root = getenv("DVC_IBEX_MATRIX_ROOT");
    const char *default_roots[] = {
        "IBEX/ibex_matrix_flat_13rate",
        "../IBEX/ibex_matrix_flat_13rate",
        "../../IBEX/ibex_matrix_flat_13rate",
        "../../../IBEX/ibex_matrix_flat_13rate",
        "../../../../IBEX/ibex_matrix_flat_13rate",
    };

    h_matrix.rows = m;
    h_matrix.cols = n;
    h_matrix.bits = sc;
    h_matrix.bytes_of_userdata = info_len / 8;
    h_matrix.bytes_of_parity = (blk_len - info_len) / 8;
    h_matrix.min_rows = 5;
    h_matrix.max_rows = LDPC_MAX_ROWS;

    if ((h_matrix.rows <= 0) || (h_matrix.rows > LDPC_MAX_ROWS) ||
        (h_matrix.cols <= 0) || (h_matrix.cols > LDPC_MAX_COLS)) {
      printf("[LDPC ERROR] Unsupported matrix size bm_m=%d bm_n=%d (max rows=%d, max cols=%d)\n",
             h_matrix.rows, h_matrix.cols, LDPC_MAX_ROWS, LDPC_MAX_COLS);
      return;
    }
    if ((h_matrix.cols - h_matrix.rows) < 0 ||
        (h_matrix.cols - h_matrix.rows) > LDPC_MAX_PAYLOAD_COLS) {
      printf("[LDPC ERROR] Unsupported payload cols=%d (max supported=%d)\n",
             h_matrix.cols - h_matrix.rows, LDPC_MAX_PAYLOAD_COLS);
      return;
    }

    h_matrix.unused_bytes_of_parity =
        (h_matrix.rows * (h_matrix.bits >> 3)) - h_matrix.bytes_of_parity;
    if (h_matrix.unused_bytes_of_parity != 0)
      h_matrix.extra_bytes_of_parity =
          (h_matrix.bits >> 3) - h_matrix.unused_bytes_of_parity;
    else
      h_matrix.extra_bytes_of_parity = 0;

    h_matrix.unused_bytes_of_userdata =
        ((h_matrix.cols - h_matrix.rows) * (h_matrix.bits >> 3)) -
        h_matrix.bytes_of_userdata;
    if (h_matrix.unused_bytes_of_userdata != 0)
      h_matrix.extra_bytes_of_userdata =
          (h_matrix.bits >> 3) - h_matrix.unused_bytes_of_userdata;
    else
      h_matrix.extra_bytes_of_userdata = 0;

    h_matrix.extra_bits_of_parity = h_matrix.extra_bytes_of_parity << 3;
    h_matrix.extra_bits_of_userdata = h_matrix.extra_bytes_of_userdata << 3;

    if (VERBOSITY > 0)
      printf("### MATRIX: BITS: %4d ROWS: %2d COLS: %3d BYTES_OF_USERDATA: %5d "
             "BYTES_OF_PARITY: %4d EXTRA_BYTES_OF_USERDATA: %5d EXTRA_BYTES_OF_PARITY: %4d\n",
             h_matrix.bits, h_matrix.rows, h_matrix.cols,
             h_matrix.bytes_of_userdata, h_matrix.bytes_of_parity,
             h_matrix.extra_bytes_of_userdata, h_matrix.extra_bytes_of_parity);

    n_full = LDPC_MAX_PAYLOAD_COLS + m;
    if (n > n_full) {
      printf("[LDPC ERROR] target n=%d > n_full=%d, cannot trim\n", n, n_full);
      return;
    }

    for (i = 0; i < LDPC_MAX_ROWS; i++) {
      h_matrix.delta[i] = 0;
      h_matrix.first_element[i] = 0;
      h_matrix.last_element[i] = 0;
      h_matrix.wraparound[i] = 0;
      h_matrix.wrap_base[i] = 0;
      h_matrix.wrap_num_deltas[i] = 0;
      h_matrix.row_weight[i] = 0;
      h_matrix_index[i] = 0;
      row_nz_cnt[i] = 0;
      for (j = 0; j < LDPC_MAX_COLS; j++) {
        h_matrix.element[i][j] = -1;
        h_matrix.occupied[i][j] = 0;
        h_matrix.fade[i][j] = 0;
        full_shift[i][j] = -1;
        full_occupied[i][j] = 0;
        full_fade[i][j] = 0;
      }
    }

    for (j = 0; j < LDPC_MAX_COLS; j++) {
      h_matrix.col_weight[j] = 0;
      h_matrix.parity_column[j] = 0;
      for (k = 0; k < h_matrix.bits; k++)
        h_matrix.mask[j][k] = 0;
    }

    for (int root_index = -1;
         root_index < (int)(sizeof(default_roots) / sizeof(default_roots[0]));
         root_index++) {
      const char *root =
          (root_index < 0) ? env_root : default_roots[root_index];
      const int source_candidates[2] = {h_matrix.cols, n_full};

      if (!root || !root[0])
        continue;

      for (int candidate_index = 0; candidate_index < 2; candidate_index++) {
        int read_ok = 1;
        int read_cols = source_candidates[candidate_index];

        if ((read_cols < h_matrix.cols) || (read_cols > LDPC_MAX_COLS))
          continue;
        if ((candidate_index == 1) && (read_cols == source_candidates[0]))
          continue;

        snprintf(matrix_path, sizeof(matrix_path),
                 "%s/matrix/LDPC_%dx%dex%d_w4_dense5_QC_H_1.txt", root,
                 h_matrix.rows, read_cols, h_matrix.bits);
        snprintf(occupied_path, sizeof(occupied_path),
                 "%s/occupied_matrix/LDPC_%dx%dex%d_w4_dense5_occupied_1.txt",
                 root, h_matrix.rows, read_cols, h_matrix.bits);
        snprintf(fade_path, sizeof(fade_path),
                 "%s/fade_matrix/LDPC_%dx%dex%d_w4_dense5_fade_1.txt", root,
                 h_matrix.rows, read_cols, h_matrix.bits);

        fp_h = fopen(matrix_path, "r");
        fp_occupied = fopen(occupied_path, "r");
        fp_fade = fopen(fade_path, "r");

        if (!fp_h || !fp_occupied || !fp_fade) {
          if (fp_h)
            fclose(fp_h);
          if (fp_occupied)
            fclose(fp_occupied);
          if (fp_fade)
            fclose(fp_fade);
          fp_h = NULL;
          fp_occupied = NULL;
          fp_fade = NULL;
          continue;
        }

        for (i = 0; i < h_matrix.rows; i++) {
          for (j = 0; j < read_cols; j++) {
            if (fscanf(fp_h, "%d", &full_shift[i][j]) != 1)
              read_ok = 0;
            if (fscanf(fp_occupied, "%d", &full_occupied[i][j]) != 1)
              read_ok = 0;
            if (fscanf(fp_fade, "%d", &full_fade[i][j]) != 1)
              read_ok = 0;
          }
        }

        fclose(fp_h);
        fclose(fp_occupied);
        fclose(fp_fade);
        fp_h = NULL;
        fp_occupied = NULL;
        fp_fade = NULL;

        if (!read_ok)
          continue;

        source_cols = read_cols;
        snprintf(source_desc, sizeof(source_desc), "%s", matrix_path);
        snprintf(source_occupied_desc, sizeof(source_occupied_desc), "%s",
                 occupied_path);
        snprintf(source_fade_desc, sizeof(source_fade_desc), "%s", fade_path);
        snprintf(source_root_desc, sizeof(source_root_desc), "%s", root);
        source_ok = 1;
        break;
      }

      if (source_ok)  
        break;
    }

    if (!source_ok) {
      printf("[LDPC ERROR] Failed to load matrix/fade/occupied files for bm_m=%d bm_n=%d. "
             "Set DVC_IBEX_MATRIX_ROOT or run from a tree containing IBEX/ibex_matrix_flat_13rate.\n",
             h_matrix.rows, h_matrix.cols);
      return;
    }

    if (source_cols != h_matrix.cols)
      printf("[LDPC] Matrix trimming: read %dx%d, trim to %dx%d (K %d->%d)\n",
             h_matrix.rows, source_cols, h_matrix.rows, h_matrix.cols,
             source_cols - h_matrix.rows, h_matrix.cols - h_matrix.rows);
    printf("[LDPC] Matrix source root: %s\n", source_root_desc);
    printf("[LDPC] Matrix path: %s\n", source_desc);
    printf("[LDPC] Occupied path: %s\n", source_occupied_desc);
    printf("[LDPC] Fade path: %s\n", source_fade_desc);

    for (i = 0; i < h_matrix.rows; i++) {
      for (j = 0; j < h_matrix.cols; j++) {
        src_col = (j < (h_matrix.cols - h_matrix.rows))
                      ? j
                      : (source_cols - h_matrix.cols + j);
        col_shift = full_shift[i][src_col];
        occupied_flag = full_occupied[i][src_col];
        fade_flag = full_fade[i][src_col];
        h_matrix.element[i][j] = -1;

        if (fade_flag == 0) {
          h_matrix.fade[i][j] = 0;
          if ((occupied_flag != 0) && (col_shift >= 0)) {
            h_matrix.occupied[i][j] = 1;
            h_matrix.last_element[i] = col_shift;
          } else {
            h_matrix.occupied[i][j] = 0;
          }
        } else {
          h_matrix.occupied[i][j] = 0;
          h_matrix.fade[i][j] = 1;
        }

        if ((col_shift >= 0) && (h_matrix.occupied[i][j] || h_matrix.fade[i][j]))
          h_matrix.last_element[i] = col_shift;
      }
    }

    h_matrix.delta[0] = 0;
    h_matrix.delta[1] = 13;
    h_matrix.delta[2] = 19;
    h_matrix.delta[3] = 29;
    h_matrix.delta[4] = 41;
    h_matrix.delta[5] = 67;
    h_matrix.delta[6] = 73;
    h_matrix.delta[7] = 79;
    h_matrix.delta[8] = 91;
    h_matrix.delta[9] = 97;
    h_matrix.delta[10] = 103;
    h_matrix.delta[11] = 111;
    h_matrix.delta[12] = 119;
    h_matrix.delta[13] = 127;
    h_matrix.delta[14] = 131;
    h_matrix.delta[15] = 137;
    h_matrix.delta[16] = 149;

    for (i = 0; i < h_matrix.rows; i++)
      h_matrix_index[i] = h_matrix.last_element[i];

    for (j = 0; j < h_matrix.cols; j++) {
      src_col = h_matrix.cols - 1 - j;
      k = src_col;
      h_matrix.parity_column[j] = (j < h_matrix.rows);
      for (i = 0; i < h_matrix.rows; i++) {
        if (h_matrix.occupied[i][src_col] || h_matrix.fade[i][src_col]) {
          h_matrix.element[i][k] = h_matrix_index[i];
          h_matrix.first_element[i] = h_matrix.element[i][k];
          h_matrix_index[i] =
              (h_matrix_index[i] + h_matrix.bits - h_matrix.delta[i]) %
              h_matrix.bits;
          row_nz_cnt[i]++;
          if (h_matrix.occupied[i][src_col]) {
            h_matrix.row_weight[i]++;
            h_matrix.col_weight[k]++;
          }
        } else
          h_matrix.element[i][k] = -1;
      }
    }

    for (i = 0; i < h_matrix.rows; i++) {
      h_matrix.wraparound[i] =
          (h_matrix.bits + h_matrix.first_element[i] - h_matrix.last_element[i]) %
          h_matrix.bits;
      h_matrix.wrap_base[i] =
          (h_matrix.delta[i] == 0)
              ? h_matrix.first_element[i]
              : (h_matrix.first_element[i] - h_matrix.delta[i] +
                 h_matrix.bits) %
                    h_matrix.bits;
      h_matrix.wrap_num_deltas[i] =
          (row_nz_cnt[i] > 0) ? (row_nz_cnt[i] - 1) : 0;
    }

    printf("[LDPC] Loaded external matrix from %s (source_cols=%d, crop_payload=%d)\n",
           source_desc, source_cols, h_matrix.cols - h_matrix.rows);

    if (h_matrix.extra_bits_of_parity > 0) {
      for (j = 0; j < h_matrix.cols; j++) {
        if (h_matrix.occupied[h_matrix.rows - 1][j]) {
          for (k = 0; k < h_matrix.bits; k++) {
            bit = (k + h_matrix.bits -
                   h_matrix.element[h_matrix.rows - 1][j]) %
                  h_matrix.bits;
            if (bit < h_matrix.extra_bits_of_parity)
              h_matrix.mask[j][k] = 1;
          }
        }
      }
    }

#ifdef _LDPC_DEBUG_DUMP
    printf("[DVC][LDEC_MATRIX] source=%s source_cols=%d rows=%d cols=%d hm_k=%d info_len=%d blk_len=%d extra_ud_bits=%d extra_pa_bits=%d\n",
           source_desc, source_cols, h_matrix.rows, h_matrix.cols, hm_k,
           info_len, blk_len, h_matrix.extra_bits_of_userdata,
           h_matrix.extra_bits_of_parity);
    for (i = 0; i < h_matrix.rows; i++) {
      printf("[DVC][LDEC_MATRIX] row=%02d FE=%03d LE=%03d WA=%03d RW=%02d\n", i,
             h_matrix.first_element[i], h_matrix.last_element[i],
             h_matrix.wraparound[i], h_matrix.row_weight[i]);
    }
    f_print_h_matrix(h_matrix);
#endif

    if (VERBOSITY > 0)
      f_print_h_matrix(h_matrix);

    ldpc_ibex_phck(h_matrix);
  }

#ifdef _LDPC_DUMP
  rdec_cmem_cont_thrshd = 10;
  rdec_hdmem_cont_thrshd = 6;

  print_hm();
#endif
} // ldpc_config




void ldpc_packet::ldpc_ibex_input(int syndrome_cal_only, int max_iter,
                                  int post_iter, int nand_strobes) {
  ldpc_decoder_input.post_iteration = post_iter;
  ldpc_decoder_input.iteration_limit = max_iter;
  ldpc_decoder_input.nand_strobes = nand_strobes;
  int soft_bits = (nand_strobes <= 0) ? 0 : (nand_strobes == 1) ? 1 : 2;
  ldpc_decoder_input.soft_bits = soft_bits;
  ldpc_decoder_input.syndrome_cal_only = syndrome_cal_only;
  int bit_index = 0;
  int bin_ibex;
  int err_cnt = 0;

  for (int j = 0; j < h_matrix.cols; j++) {
    for (int k = 0; k < h_matrix.bits; k++) {
      if ((j == (h_matrix.cols - h_matrix.rows - 1)) &&
          (k >= (h_matrix.bits - 8 * h_matrix.unused_bytes_of_userdata))) {
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = 0;
        if (soft_bits >= 1)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable =
              ldpc_decoder_parameters.likelihood_map[3];
        if (soft_bits >= 2)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
              ldpc_decoder_parameters.likelihood_map[3];
        else if (soft_bits >= 1)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
              ldpc_decoder_input.corrupted_codeword.c[j]
                  .b[k]
                  .bit_questionable; // replicate
      } else if ((j == (h_matrix.cols - h_matrix.rows)) &&
                 (k >= (h_matrix.bits - 8 * h_matrix.unused_bytes_of_parity))) {
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = 0;
        if (soft_bits >= 1)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable =
              ldpc_decoder_parameters.likelihood_map[3];
        if (soft_bits >= 2)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
              ldpc_decoder_parameters.likelihood_map[3];
        else if (soft_bits >= 1)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
              ldpc_decoder_input.corrupted_codeword.c[j]
                  .b[k]
                  .bit_questionable; // replicate
      } else {
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard =
            (rx_blk[bit_index] >= 0) ? 0 : 1;
        // debug
        err_cnt += tx_blk[bit_index] != ((rx_blk[bit_index] >= 0) ? 0 : 1);
        if (rx_blk[bit_index] >= 0) {
          bin_ibex = (det_blk[bit_index] + rd_num - 1) % (rd_num + 1);
          if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable =
                (bin_ibex / 2) % 2;
          if (soft_bits >= 2)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
                (bin_ibex / 4) % 2;
          else if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
                ldpc_decoder_input.corrupted_codeword.c[j]
                    .b[k]
                    .bit_questionable; // replicate
        } else {
          bin_ibex = det_blk[bit_index];
          if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable =
                (bin_ibex / 2) % 2;
          if (soft_bits >= 2)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
                (bin_ibex / 4) % 2;
          else if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
                ldpc_decoder_input.corrupted_codeword.c[j]
                    .b[k]
                    .bit_questionable; // replicate
        }
        bit_index++;
      }
    }
  }

  if (err_cnt != raw_err_num)
    printf("### LDPC IBEX INPUT ERROR: ERR_CNT %d != RAW_ERR_NUM %d\n", err_cnt,
           raw_err_num);
}

void ldpc_packet::ldpc_ibex_parameters(
    int post_process_en, int syndrome_weight_thr_qc,
    int syndrome_weight_thr_post, int early_termination_dis,
    int ldpc_decoder_control_likelihood_0,
    int ldpc_decoder_control_likelihood_1,
    int ldpc_decoder_control_likelihood_2,
    int ldpc_decoder_control_likelihood_3, int ldpc_decoder_control_post,
    unsigned int ldpc_early_term_0, unsigned int ldpc_early_term_1,
    unsigned int ldpc_early_term_2, unsigned int ldpc_early_term_3,
    unsigned int ldpc_early_term_4, unsigned int ldpc_early_term_5,
    unsigned int ldpc_early_term_6) {
  ldpc_decoder_parameters.post_process_en = post_process_en;
  ldpc_decoder_parameters.syndrome_weight_thr_qc = syndrome_weight_thr_qc;
  ldpc_decoder_parameters.syndrome_weight_thr_post = syndrome_weight_thr_post;

  const unsigned int post_cfg =
      static_cast<unsigned int>(ldpc_decoder_control_post);
  ldpc_decoder_parameters.likelihood_thr = (post_cfg >> 16) & 0xff;
  ldpc_decoder_parameters.post_ratio = (post_cfg >> 24) & 0xf;

  const unsigned int like0 =
      static_cast<unsigned int>(ldpc_decoder_control_likelihood_0);
  const unsigned int like1 =
      static_cast<unsigned int>(ldpc_decoder_control_likelihood_1);
  const unsigned int like2 =
      static_cast<unsigned int>(ldpc_decoder_control_likelihood_2);
  const unsigned int like3 =
      static_cast<unsigned int>(ldpc_decoder_control_likelihood_3);

  ldpc_decoder_parameters.likelihood_init_coef_all[0][0] = (like0 >> 0) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[0][1] = (like0 >> 4) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[0][2] = (like0 >> 8) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[0][3] = (like0 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][0] = (like0 >> 16) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][1] = (like0 >> 20) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][2] = (like0 >> 24) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][3] = (like0 >> 28) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][0] = (like1 >> 0) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][1] = (like1 >> 4) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][2] = (like1 >> 8) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][3] = (like1 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][0] = (like1 >> 16) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][1] = (like1 >> 20) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][2] = (like1 >> 24) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][3] = (like1 >> 28) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][0] = (like2 >> 0) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][1] = (like2 >> 4) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][2] = (like2 >> 8) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][3] = (like2 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][0] = (like2 >> 16) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][1] = (like2 >> 20) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][2] = (like2 >> 24) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][3] = (like2 >> 28) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][0] = (like3 >> 0) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][1] = (like3 >> 4) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][2] = (like3 >> 8) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][3] = (like3 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][0] = (like3 >> 16) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][1] = (like3 >> 20) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][2] = (like3 >> 24) & 0xf;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][3] = (like3 >> 28) & 0xf;

  ldpc_decoder_parameters.likelihood_init_fraction[1] =
      (ldpc_early_term_0 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_init_fraction[2] =
      (ldpc_early_term_0 >> 28) & 0xf;

  ldpc_decoder_parameters.likelihood_map[0] = (ldpc_early_term_1 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_map[1] = (ldpc_early_term_1 >> 28) & 0xf;
  ldpc_decoder_parameters.likelihood_map[2] = (ldpc_early_term_2 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_map[3] = (ldpc_early_term_2 >> 28) & 0xf;
  ldpc_decoder_parameters.likelihood_map[4] = (ldpc_early_term_3 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_map[5] = (ldpc_early_term_3 >> 28) & 0xf;
  ldpc_decoder_parameters.likelihood_map[6] = (ldpc_early_term_4 >> 12) & 0xf;
  ldpc_decoder_parameters.likelihood_map[7] = (ldpc_early_term_4 >> 28) & 0xf;

  ldpc_decoder_parameters.early_terminate_dis = early_termination_dis;
  ldpc_decoder_parameters.early_terminate_thr[0][0] = 600;
  ldpc_decoder_parameters.early_terminate_thr[0][1] = 750;
  ldpc_decoder_parameters.early_terminate_thr[0][2] = 900;
  ldpc_decoder_parameters.early_terminate_thr[0][3] = 1050;
  ldpc_decoder_parameters.early_terminate_thr[0][4] = 1200;
  ldpc_decoder_parameters.early_terminate_thr[0][5] = 1350;
  ldpc_decoder_parameters.early_terminate_thr[0][6] = 1500;
  ldpc_decoder_parameters.early_terminate_thr[1][0] = 1200;
  ldpc_decoder_parameters.early_terminate_thr[1][1] = 1500;
  ldpc_decoder_parameters.early_terminate_thr[1][2] = 1800;
  ldpc_decoder_parameters.early_terminate_thr[1][3] = 2100;
  ldpc_decoder_parameters.early_terminate_thr[1][4] = 2400;
  ldpc_decoder_parameters.early_terminate_thr[1][5] = 2700;
  ldpc_decoder_parameters.early_terminate_thr[1][6] = 3000;

  int soft_bit_control_0 = 0;
  unsigned int c_soft_data_table[8]; // 32 bit size / 4 bits each field = 8
  for (int i = 0; i < 8; i++) {
    const int shift_bit = i * 4; // 4 bit at a time
    c_soft_data_table[i] =
        (static_cast<unsigned int>(soft_bit_control_0) >> shift_bit) & 0xf;
  }
  for (int i = 0; i < 8; i++)
    ldpc_decoder_parameters.soft_bit_table[i] = c_soft_data_table[i];

  (void)ldpc_early_term_5;
  (void)ldpc_early_term_6;
}

void ldpc_packet::ldpc_apply_sdlite_llr_override() {
  if (!reg_sdlite_llr_config_en)
    return;

  int cfg_len = 0;
  if (rd_num == 3)
    cfg_len = 4;
  else if (rd_num == 5)
    cfg_len = 6;
  else if (rd_num == 7)
    cfg_len = 8;
  else
    return;

  if ((llr_tbl == NULL) || (bin_num < cfg_len))
    return;

  const double scale = pow(2, -1 * finite_f_num);
  for (int i = 0; i < cfg_len; i++)
    llr_tbl[i] = (float)(reg_sdlite_llr_table[i] * scale);
}

// LDPC decoder config
void ldpc_packet::ldpc_dec_config(int max_fdec_itr, int fdec_col_skip,
                                  int max_ldec_itr, float dec_alpha,
                                  const float *dec_alpha_pms,
                                  const float *dec_beta_pms,
                                  float dec_point1, float dec_point2,
                                  int dec_pms_lut_len, int fin_mode,
                                  int fin_q_num, int fin_r_num,
                                  int fin_f_num) {
  // syndrome weight
  init_synd_wt_min = hm_m;
  init_synd_wt_max = 0;

  reg_sdlite_llr_config_en = 0;
  for (int i = 0; i < 8; i++)
    reg_sdlite_llr_table[i] = 0;

  // layer config
  if (max_ldec_itr > 0) {
    ldec_max_itr = max_ldec_itr;
    ldec_early_term_en = 1;
  } else {
    ldec_max_itr = -1 * max_ldec_itr;
    ldec_early_term_en = 0;
  }

  alpha = dec_alpha;
  for (int i = 0; i < LDPC_PMS_LUT_SIZE; i++) {
    alpha_pms[i] = dec_alpha;
    beta_pms[i] = 0.0f;
  }
  if ((dec_alpha_pms != nullptr) && (dec_pms_lut_len > 0)) {
    int valid_cnt = dec_pms_lut_len;
    if (valid_cnt > LDPC_PMS_LUT_SIZE)
      valid_cnt = LDPC_PMS_LUT_SIZE;
    for (int i = 0; i < valid_cnt; i++)
      alpha_pms[i] = dec_alpha_pms[i];
  }
  if ((dec_beta_pms != nullptr) && (dec_pms_lut_len > 0)) {
    int valid_cnt = dec_pms_lut_len;
    if (valid_cnt > LDPC_PMS_LUT_SIZE)
      valid_cnt = LDPC_PMS_LUT_SIZE;
    for (int i = 0; i < valid_cnt; i++)
      beta_pms[i] = dec_beta_pms[i];
  }
  point1 = dec_point1;
  point2 = dec_point2;

  // initialize quantization
  finite_mode = fin_mode;
  finite_q_num = fin_q_num;
  finite_r_num = fin_r_num;
  finite_c_num = finite_r_num;
  finite_f_num = fin_f_num;

  int tmp = pow(2, finite_f_num);
  finite_q_max = (pow(2, finite_q_num - 1) - 1) / tmp;
  finite_q_min = -(pow(2, finite_q_num - 1) - 1) / tmp;
  finite_r_max = (pow(2, finite_r_num - 1) - 1) / tmp;
  finite_r_min = -(pow(2, finite_r_num - 1) - 1) / tmp;
  finite_c_max = (pow(2, finite_c_num - 1) - 1) / tmp;
  finite_c_min = 0;

#ifdef _LDPC_DEBUG
  printf("[LDPC DEBUG] Q MSG: %d/%d, max %f, min %f \n", finite_q_num,
         fin_f_num, finite_q_max, finite_q_min);
  printf("[LDPC DEBUG] R MSG: %d/%d, max %f, min %f \n", finite_r_num,
         fin_f_num, finite_r_max, finite_r_min);
  printf("[LDPC DEBUG] C MSG: %d/%d, max %f, min %f \n", finite_c_num - 1,
         fin_f_num, finite_c_max, finite_c_min);
#endif

  // BF config
  if (max_fdec_itr > 0) {
    fdec_max_itr = max_fdec_itr;
    fdec_early_term_en = 1;
  } else {
    fdec_max_itr = -1 * max_fdec_itr;
    fdec_early_term_en = 0;
  }

  col_skip_itr = fdec_col_skip;

  flp_thrshd0 = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd0));
  flp_thrshd1 = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd1));

  flp_thrshd0[0] = 6;
  flp_thrshd0[1] = 5;
  flp_thrshd0[2] = 5;
  flp_thrshd0[3] = 5;
  flp_thrshd0[4] = 5;
  flp_thrshd0[5] = 5;
  flp_thrshd0[6] = 5;
  flp_thrshd0[7] = 5;
  for (int i = 8; i < fdec_max_itr; i++)
    flp_thrshd0[i] = 4;

  flp_thrshd1[0] = 3;
  flp_thrshd1[1] = 1;
  flp_thrshd1[2] = 1;
  flp_thrshd1[3] = 1;
  flp_thrshd1[4] = 1;
  flp_thrshd1[5] = 1;
  flp_thrshd1[6] = 2;
  flp_thrshd1[7] = 4;
  for (int i = 8; i < fdec_max_itr; i++)
    flp_thrshd1[i] = 2;

  // special iteration
  // break some deadlocks that would degrade performance
  // especially for matrix 143x14wt6

  flp_thrshd0_s = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd0_s));
  flp_thrshd0_w = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd0_w));
  flp_thrshd1_s = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd1_s));
  flp_thrshd1_w = (int *)calloc(fdec_max_itr, sizeof(*flp_thrshd1_w));

  flp_thrshd0_s[0] = 6;
  flp_thrshd0_s[1] = 6;
  flp_thrshd0_s[2] = 6;
  flp_thrshd0_s[3] = 6;
  flp_thrshd0_s[4] = 6;
  flp_thrshd0_s[5] = 6;
  flp_thrshd0_s[6] = 6;
  flp_thrshd0_s[7] = 6;
  for (int i = 8; i < fdec_max_itr; i++)
    flp_thrshd0_s[i] = 5;

  flp_thrshd0_w[0] = 6;
  flp_thrshd0_w[1] = 5;
  flp_thrshd0_w[2] = 5;
  flp_thrshd0_w[3] = 5;
  flp_thrshd0_w[4] = 5;
  flp_thrshd0_w[5] = 5;
  flp_thrshd0_w[6] = 5;
  flp_thrshd0_w[7] = 5;
  for (int i = 8; i < fdec_max_itr; i++)
    flp_thrshd0_w[i] = 4;

  flp_thrshd1_s[0] = 3;
  flp_thrshd1_s[1] = 2;
  flp_thrshd1_s[2] = 2;
  flp_thrshd1_s[3] = 2;
  flp_thrshd1_s[4] = 2;
  flp_thrshd1_s[5] = 2;
  flp_thrshd1_s[6] = 2;
  flp_thrshd1_s[7] = 4;
  for (int i = 8; i < fdec_max_itr; i++)
    flp_thrshd1_s[i] = 3;

  flp_thrshd1_w[0] = 3;
  flp_thrshd1_w[1] = 1;
  flp_thrshd1_w[2] = 1;
  flp_thrshd1_w[3] = 1;
  flp_thrshd1_w[4] = 1;
  flp_thrshd1_w[5] = 1;
  flp_thrshd1_w[6] = 2;
  flp_thrshd1_w[7] = 4;
  for (int i = 8; i < fdec_max_itr; i++)
    flp_thrshd1_w[i] = 2;

  sb_thrshd0_s0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_s0));
  sb_thrshd0_s1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_s1));
  sb_thrshd0_w0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_w0));
  sb_thrshd0_w1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd0_w1));
  sb_thrshd1_s0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_s0));
  sb_thrshd1_s1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_s1));
  sb_thrshd1_w0 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_w0));
  sb_thrshd1_w1 = (int *)calloc(fdec_max_itr, sizeof(*sb_thrshd1_w1));

  sb_thrshd0_s0[0] = 6;
  sb_thrshd0_s0[1] = 4;
  sb_thrshd0_s0[2] = 4;
  sb_thrshd0_s0[3] = 4;
  sb_thrshd0_s0[4] = 4;
  sb_thrshd0_s0[5] = 4;
  sb_thrshd0_s0[6] = 4;
  sb_thrshd0_s0[7] = 4;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd0_s0[i] = 3;

  sb_thrshd0_s1[0] = 7;
  sb_thrshd0_s1[1] = 7;
  sb_thrshd0_s1[2] = 7;
  sb_thrshd0_s1[3] = 7;
  sb_thrshd0_s1[4] = 7;
  sb_thrshd0_s1[5] = 7;
  sb_thrshd0_s1[6] = 7;
  sb_thrshd0_s1[7] = 7;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd0_s1[i] = 7;

  sb_thrshd0_w0[0] = 3;
  sb_thrshd0_w0[1] = 3;
  sb_thrshd0_w0[2] = 3;
  sb_thrshd0_w0[3] = 3;
  sb_thrshd0_w0[4] = 3;
  sb_thrshd0_w0[5] = 3;
  sb_thrshd0_w0[6] = 3;
  sb_thrshd0_w0[7] = 3;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd0_w0[i] = 2;

  sb_thrshd0_w1[0] = 6;
  sb_thrshd0_w1[1] = 6;
  sb_thrshd0_w1[2] = 6;
  sb_thrshd0_w1[3] = 6;
  sb_thrshd0_w1[4] = 6;
  sb_thrshd0_w1[5] = 6;
  sb_thrshd0_w1[6] = 6;
  sb_thrshd0_w1[7] = 6;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd0_w1[i] = 6;

  sb_thrshd1_s0[0] = 0;
  sb_thrshd1_s0[1] = 0;
  sb_thrshd1_s0[2] = 0;
  sb_thrshd1_s0[3] = 0;
  sb_thrshd1_s0[4] = 0;
  sb_thrshd1_s0[5] = 0;
  sb_thrshd1_s0[6] = 0;
  sb_thrshd1_s0[7] = 0;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd1_s0[i] = 0;

  sb_thrshd1_s1[0] = 6;
  sb_thrshd1_s1[1] = 6;
  sb_thrshd1_s1[2] = 6;
  sb_thrshd1_s1[3] = 6;
  sb_thrshd1_s1[4] = 6;
  sb_thrshd1_s1[5] = 6;
  sb_thrshd1_s1[6] = 6;
  sb_thrshd1_s1[7] = 6;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd1_s1[i] = 6;

  sb_thrshd1_w0[0] = 0;
  sb_thrshd1_w0[1] = 0;
  sb_thrshd1_w0[2] = 0;
  sb_thrshd1_w0[3] = 0;
  sb_thrshd1_w0[4] = 0;
  sb_thrshd1_w0[5] = 0;
  sb_thrshd1_w0[6] = 0;
  sb_thrshd1_w0[7] = 0;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd1_w0[i] = 0;

  sb_thrshd1_w1[0] = 5;
  sb_thrshd1_w1[1] = 5;
  sb_thrshd1_w1[2] = 5;
  sb_thrshd1_w1[3] = 5;
  sb_thrshd1_w1[4] = 5;
  sb_thrshd1_w1[5] = 5;
  sb_thrshd1_w1[6] = 5;
  sb_thrshd1_w1[7] = 5;
  for (int i = 8; i < fdec_max_itr; i++)
    sb_thrshd1_w1[i] = 5;
} // ldpc_dec_config

void ldpc_packet::ldpc_clean() {
  if (qc_bm) {
    mod2sparse_free(qc_bm);
    qc_bm = NULL;
  }
  if (qc_hm) {
    mod2sparse_free(qc_hm);
    qc_hm = NULL;
  }
  if (qc_a) {
    mod2sparse_free(qc_a);
    qc_a = NULL;
  }
  if (qc_b) {
    mod2sparse_free(qc_b);
    qc_b = NULL;
  }
  if (qc_c) {
    mod2sparse_free(qc_c);
    qc_c = NULL;
  }
  if (qc_d) {
    mod2sparse_free(qc_d);
    qc_d = NULL;
  }
  if (qc_e) {
    mod2sparse_free(qc_e);
    qc_e = NULL;
  }
  if (qc_fi) {
    mod2sparse_free(qc_fi);
    qc_fi = NULL;
  }

  free(flp_thrshd0);
  flp_thrshd0 = NULL;
  free(flp_thrshd1);
  flp_thrshd1 = NULL;
  free(flp_thrshd0_s);
  flp_thrshd0_s = NULL;
  free(flp_thrshd1_s);
  flp_thrshd1_s = NULL;
  free(flp_thrshd0_w);
  flp_thrshd0_w = NULL;
  free(flp_thrshd1_w);
  flp_thrshd1_w = NULL;
  free(sb_thrshd0_s0);
  sb_thrshd0_s0 = NULL;
  free(sb_thrshd0_s1);
  sb_thrshd0_s1 = NULL;
  free(sb_thrshd1_s0);
  sb_thrshd1_s0 = NULL;
  free(sb_thrshd1_s1);
  sb_thrshd1_s1 = NULL;
  free(sb_thrshd0_w0);
  sb_thrshd0_w0 = NULL;
  free(sb_thrshd0_w1);
  sb_thrshd0_w1 = NULL;
  free(sb_thrshd1_w0);
  sb_thrshd1_w0 = NULL;
  free(sb_thrshd1_w1);
  sb_thrshd1_w1 = NULL;
}

void ldpc_packet::ldpc_pckt_alloc() {
  ch_pckt_alloc();

  usr_blk = (char *)calloc(info_len, sizeof(*usr_blk));
  enc_di_blk = (char *)calloc(hm_k, sizeof(*enc_di_blk));
  enc_do_blk = (char *)calloc(hm_n, sizeof(*enc_do_blk));
  dec_di_blk = (char *)calloc(hm_n, sizeof(*dec_di_blk));
  dec_do_blk = (char *)calloc(hm_n, sizeof(*dec_do_blk));
  dec_blk = (char *)calloc(blk_len, sizeof(*dec_blk));

  cw_fail = 0;
  cw_miscorr = 0;
  cor_err_num = 0;
  dec_err_num = 0;
}

void ldpc_packet::ldpc_pckt_clean() {
  ch_pckt_clean();

  free(usr_blk);
  usr_blk = NULL;
  free(dec_blk);
  dec_blk = NULL;
  free(enc_di_blk);
  enc_di_blk = NULL;
  free(enc_do_blk);
  enc_do_blk = NULL;
  free(dec_di_blk);
  dec_di_blk = NULL;
  free(dec_do_blk);
  dec_do_blk = NULL;

} // ldpc_pckt_free

void ldpc_packet::ldpc_encoder() {
  char *az1, *eaz1, *cz1, *sumz1, *z2, *bz2, *z3;

  az1 = (char *)calloc(tm_sz * cir_sz, sizeof(*az1));
  eaz1 = (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*eaz1));
  cz1 = (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*cz1));
  sumz1 = (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*sumz1));
  z2 = (char *)calloc((bm_m - tm_sz) * cir_sz, sizeof(*z2));
  bz2 = (char *)calloc(tm_sz * cir_sz, sizeof(*bz2));
  z3 = (char *)calloc(tm_sz * cir_sz, sizeof(*z3));

  // padding 0s
  vec_copy(usr_blk, enc_di_blk, 0, 0, info_len);
  for (int i = 0; i < pad_len; i++)
    enc_di_blk[hm_k - pad_len + i] = 0;

  // A*Z1
  mod2sparse_mulvec(qc_a, enc_di_blk, az1);
  // E*(A*Z1)
  mod2sparse_mulvec(qc_e, az1, eaz1);
  // C*Z1
  mod2sparse_mulvec(qc_c, enc_di_blk, cz1);
  // E*(A*Z1)+C*Z1
  vec_mod2_add(eaz1, cz1, sumz1, (bm_m - tm_sz) * cir_sz);
  // Z2 = F_inv*[E*(A*Z1)+C*Z1]
  mod2sparse_mulvec(qc_fi, sumz1, z2);
  // B*Z2
  mod2sparse_mulvec(qc_b, z2, bz2);
  // Z3 = BZ2 +AZ1
  vec_mod2_add(az1, bz2, z3, tm_sz * cir_sz);

  // encoded data
  vec_copy(enc_di_blk, enc_do_blk, 0, 0, hm_k);
  vec_copy(z2, enc_do_blk, 0, hm_k, (bm_m - tm_sz) * cir_sz);
  vec_copy(z3, enc_do_blk, 0, (bm_n - tm_sz) * cir_sz, tm_sz * cir_sz);

  // removing 0 padding
  vec_copy(enc_do_blk, tx_blk, 0, 0, info_len);
  vec_copy(enc_do_blk, tx_blk, hm_k, info_len, hm_m);

  // free space
  free(az1);
  free(eaz1);
  free(cz1);
  free(sumz1);
  free(z2);
  free(bz2);
  free(z3);
}

void ldpc_packet::ldpc_ibex_encoder() {
  // #include "ldpc_matrix.h"
  // #include "ldpc_matrix_inverse.h"

  int VERBOSITY = 0;
  int DEBUG_MODE = 0;

  int i;
  int j;
  int k;

  s_hard_codeword ldpc_encoder_input;
  int user_data_bits = h_matrix.bytes_of_userdata * 8;
  int parity_bits = h_matrix.bytes_of_parity * 8;
  int skip_parity_bits_in_first_parity_column =
      (h_matrix.rows * h_matrix.bits) - parity_bits;
  int last_parity_bits_in_first_parity_column =
      h_matrix.bits - 1 - skip_parity_bits_in_first_parity_column;
  int first_parity = ((h_matrix.cols * h_matrix.bits) - parity_bits - 1);
  int bit_location = 0;
  int unpacked_bit_location = 0;

  for (j = 0; j < h_matrix.cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      if (bit_location < user_data_bits)
        ldpc_encoder_input.c[j].b[k] = usr_blk[bit_location++];
      else
        ldpc_encoder_input.c[j].b[k] = 0;

      if (DEBUG_MODE)
        ldpc_encoder_input.c[j].b[k] = 0;
    }
  }

  s_hard_codeword ldpc_encoder_output =
      f_ldpc_encode(ldpc_encoder_input, h_matrix);

  bit_location = 0;

  char c;

  for (j = 0; j < h_matrix.cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      if (bit_location < user_data_bits) {
        tx_blk[bit_location++] = ldpc_encoder_output.c[j].b[k];
        c = ldpc_encoder_output.c[j].b[k] ? '1' : '0';
        unpacked_bit_location++;
      } else if (j == (h_matrix.cols - h_matrix.rows)) {
        if (k < (512 - skip_parity_bits_in_first_parity_column)) {
          tx_blk[bit_location++] = ldpc_encoder_output.c[j].b[k];
          c = ldpc_encoder_output.c[j].b[k] ? '1' : '0';
        }
        unpacked_bit_location++;
      } else if (j > (h_matrix.cols - h_matrix.rows)) {
        tx_blk[bit_location++] = ldpc_encoder_output.c[j].b[k];
        c = ldpc_encoder_output.c[j].b[k] ? '1' : '0';
        unpacked_bit_location++;
      } else {
        unpacked_bit_location++;
      }
    }
  }
}

s_hard_codeword ldpc_packet::f_ldpc_encode(s_hard_codeword ldpc_encoder_input,
                                           s_h_matrix h_matrix) {
  int VERBOSITY = 0;
  int i;
  int j;
  int k;
  int m;
  int cv;
  int ch;
  s_hard_codeword vn;
  s_hard_codeword vnpf;
  s_check_nodes cn;
  int matrix_sel;

  s_hard_codeword codeword_payload;

  bool payload[LDPC_MAX_PAYLOAD_COLS][LDPC_MAX_CIRC_BITS];
  bool parity[LDPC_MAX_ROWS][LDPC_MAX_CIRC_BITS];
  bool check_node[LDPC_MAX_ROWS][LDPC_MAX_CIRC_BITS];
  bool ldpc_matrix_occupied[LDPC_MAX_ROWS][LDPC_MAX_COLS];
  bool ldpc_matrix_fade[LDPC_MAX_ROWS][LDPC_MAX_COLS];
  unsigned int ldpc_matrix[LDPC_MAX_ROWS][LDPC_MAX_COLS];
  int ldpc_encoder_failure;
  int bit_rotated;
  int syndrome_weight;

#include "ldpc_matrix_inverse.h"

  int num_bytes = h_matrix.bits >> 3;
  int num_payload_cols =
      int((h_matrix.bytes_of_userdata + num_bytes - 1) / num_bytes);
  int num_payload_bits = 8 * h_matrix.bytes_of_userdata;
  int num_parity_cols = h_matrix.rows;
  int extra_payload_cols = num_payload_cols - 64;
  int unused_parity_bytes =
      (h_matrix.rows * num_bytes) - h_matrix.bytes_of_parity;
  int unused_parity_bits = unused_parity_bytes * 8;
  int matrix_element[LDPC_MAX_ROWS];

  int bit_location;
  int byte_data;
  int index;

  // Set the initial VN to the input data (payload part will be correct, parity
  // part will be 0)
  vn = ldpc_encoder_input;
  // compute the CN for the payload
  cn = f_check_nodes(h_matrix, vn);

  // Apply the first part of the parity matrix, between the user data and the
  // last 5 columns

  for (j = h_matrix.cols - h_matrix.rows; j < h_matrix.cols - 5; j++) {
    m = h_matrix.cols - j - 1;
    for (k = 0; k < h_matrix.bits; k++) {
      vn.c[j].b[k] = cn.r[m].b[k];
    }
    // int shift = h_matrix.element[m][j];
    // for (k = 0; k < h_matrix.bits; k++)
    //   vn.c[j].b[k] = cn.r[m].b[(k + shift) % h_matrix.bits];

    for (k = 0; k < h_matrix.bits; k++) {
      if (vn.c[j].b[k]) {
        for (i = 0; i < h_matrix.rows; i++) {
          m = (k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits;
          if (h_matrix.extra_bytes_of_parity == 0) {
            if (h_matrix.occupied[i][j])
              cn.r[i].b[m] = !cn.r[i].b[m];
          } else {
            if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)))
              cn.r[i].b[m] = !cn.r[i].b[m];
            if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
                h_matrix.mask[j][k])
              cn.r[i].b[m] = !cn.r[i].b[m];
            if (h_matrix.fade[i][j] && !h_matrix.mask[j][k])
              cn.r[i].b[m] = !cn.r[i].b[m];
          }
        }
      }
    }
  }

  // Apply the 5x5 parity matrix to get the last 5 columns of parity
  for (j = 0; j < 5; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      parity[h_matrix.rows - 5 + j][k] = 0;
      for (i = 0; i < 5; i++) {
        for (m = 0; m < h_matrix.bits; m++) {
          parity[h_matrix.rows - 5 + j][k] ^=
              (cn.r[i].b[m] &
               encoder_matrix_first_col[j][i][(h_matrix.bits + m - k) %
                                              h_matrix.bits]);
        }
      }
    }
  }

  for (j = num_parity_cols - 5; j < num_parity_cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      if ((j > 0) || (unused_parity_bytes == 0) ||
          ((j == 0) && (k >= unused_parity_bits)))
        vn.c[j + h_matrix.cols - h_matrix.rows].b[k] = parity[j][k];
    }
  }

  cn = f_check_nodes(h_matrix, vn);
  syndrome_weight = f_check_node_weight(h_matrix, cn);

  dvc_ensure_output_c_code_dir();
  if (FILE *swfp = fopen("./output/c_code/enc_synd_wt.txt", "a")) {
    fprintf(swfp, "%d\n", syndrome_weight);
    fclose(swfp);
  } else {
    printf("[LDPC WARN] fopen failed: %s (errno=%d)\n",
           "./output/c_code/enc_synd_wt.txt", errno);
  }

  ldpc_encoder_failure = (syndrome_weight > 0);
  if (ldpc_encoder_failure != 0) {
    printf("ERROR: Encoding! Syndrome weight = %d > 0\n", syndrome_weight);
    exit(1);
  }

  return vn;
}

// [DVC IBEX PROOF] Dump a bit buffer (`bits[0..len-1]`, each char !=0 treated as 1)
// matching DV SV `$displayh(logic[N-1:0] data)` big-integer convention:
//   - bits[i] goes to logical bit `i` of the integer (LSB-first: word bit(i&31))
//   - The full row is printed MSB → LSB, so word0 ends up RIGHTMOST, word15 LEFTMOST.
//     Example: bits[0..3]=1,1,1,0 and the rest = 0 → row = "00000000...00000007".
// One row = 16 words = 512 bits = one cir_sz=512 circulant column.
static void dvc_proof_dump_bits_hex_rows(const char *tag, const char *bits,
                                         int len) {
  if (!bits || len <= 0) {
    printf("[DVC IBEX PROOF] %s: <empty> len=%d\n", tag ? tag : "?", len);
    return;
  }
  unsigned int words[16] = {0};
  int word_idx = 0;
  int row_idx = 0;
  int start_bit = 0;
  int i;
  for (i = 0; i < len; i++) {
    if (bits[i] & 1)
      words[word_idx] |= (1u << (i & 31));
    if (((i + 1) & 31) == 0) {
      word_idx++;
      if (word_idx == 16) {
        char row[16 * 8 + 1];
        int pos = 0;
        // word0 is rightmost in the big-integer view → emit in reverse.
        for (int w = 15; w >= 0; w--)
          pos += snprintf(row + pos, sizeof(row) - pos, "%08x", words[w]);
        row[pos] = '\0';
        printf("[DVC IBEX PROOF] %s row%03d [bit %05d..%05d]: %s\n", tag,
               row_idx, start_bit, i, row);
        row_idx++;
        for (int w = 0; w < 16; w++)
          words[w] = 0;
        word_idx = 0;
        start_bit += 512;
      }
    }
  }
  // Tail: include any partial word, then flush remaining words in reverse order.
  if ((i & 31) != 0)
    word_idx++;
  if (word_idx > 0) {
    char row[16 * 8 + 1];
    int pos = 0;
    for (int w = word_idx - 1; w >= 0; w--)
      pos += snprintf(row + pos, sizeof(row) - pos, "%08x", words[w]);
    row[pos] = '\0';
    printf("[DVC IBEX PROOF] %s row%03d [bit %05d..%05d]: %s\n", tag, row_idx,
           start_bit, len - 1, row);
  }
}

void ldpc_packet::ldpc_decoder(enum dec_model dec_mode) {
  // add 0 padding
  if ((dec_mode != BF_IBEX) && (dec_mode != BF_IBEX_RTL_CN) &&
      (dec_mode != LAYER_G2)) {
    vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
    for (int i = 0; i < pad_len; i++)
      dec_di_blk[info_len + i] = max_llr_bin;
  }
  if (dec_mode == LAYER_G2) {
    vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
    for (int i = 0; i < pad_len; i++)
      dec_di_blk[info_len + i] = max_llr_bin;
    // `LAYER_G2` supports IBEX shortening mode where the first parity column is
    // fractional
    // (`extra_bits_of_parity > 0`). For a full-parity codeword
    // (`extra_bits_of_parity == 0`), parity is contiguous `hm_m` bits and must
    // start from the first parity column.
    if (h_matrix.extra_bits_of_parity > 0) {
      vec_copy(det_blk, dec_di_blk, info_len, hm_k,
               h_matrix.extra_bits_of_parity);
      for (int i = 0; i < h_matrix.unused_bytes_of_parity * 8; i++)
        dec_di_blk[hm_k + h_matrix.extra_bits_of_parity + i] = max_llr_bin;
      vec_copy(det_blk, dec_di_blk, info_len + h_matrix.extra_bits_of_parity,
               hm_k + cir_sz, (bm_m - 1) * cir_sz);
    } else {
      vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);
    }
  }

  // NOTE: for IBEX (cir_sz==512), the transmitted block `det_blk` does not
  // contain the padded tail bits in the first parity column when
  // `unused_bytes_of_parity > 0`. The decoder input `dec_di_blk` is always
  // sized to `hm_n` (full QC length), so we must explicitly insert those
  // missing bits instead of blindly copying `hm_m` bits (which would read past
  // `det_blk`).
  if ((cir_sz == 512) && (h_matrix.unused_bytes_of_parity > 0)) {
    const int skip_parity_bits_in_first_parity_column =
        h_matrix.unused_bytes_of_parity * 8;
    const int valid_parity_bits_in_first_parity_column =
        h_matrix.bits - skip_parity_bits_in_first_parity_column;
    const int parity_bits_in_det = blk_len - info_len;

    if (valid_parity_bits_in_first_parity_column > 0) {
      vec_copy(det_blk, dec_di_blk, info_len, hm_k,
               valid_parity_bits_in_first_parity_column);
    }

    for (int i = 0; i < skip_parity_bits_in_first_parity_column; i++)
      dec_di_blk[hm_k + valid_parity_bits_in_first_parity_column + i] =
          max_llr_bin;

    const int remaining_parity_bits =
        parity_bits_in_det - valid_parity_bits_in_first_parity_column;
    if (remaining_parity_bits > 0) {
      vec_copy(det_blk, dec_di_blk,
               info_len + valid_parity_bits_in_first_parity_column,
               hm_k + h_matrix.bits, remaining_parity_bits);
    }
  } else {
    vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);
  }

  // [DVC IBEX PROOF] Right before dispatch, dump what the decoder is actually
  // about to see. Compare three logs end-to-end to pinpoint where data drifts:
  //   (1) DV SV-side TX/err_data print
  //   (2) ldpc_c_model.c ldpc_dec entry "err_data row..." dump
  //   (3) THIS dump (det_blk + dec_di_blk right before the decoder runs)
  // All three should be bit-identical on the payload bits; any divergence
  // localizes the bug to the segment between consecutive prints.
  {
    printf("[DVC IBEX PROOF] ldpc_decoder dispatch: dec_mode=%d info_len=%d "
           "pad_len=%d blk_len=%d hm_k=%d hm_m=%d hm_n=%d max_llr_bin=%d\n",
           (int)dec_mode, info_len, pad_len, blk_len, hm_k, hm_m, hm_n,
           max_llr_bin);
    dvc_proof_dump_bits_hex_rows("det_blk", det_blk, blk_len);
    dvc_proof_dump_bits_hex_rows("dec_di_blk", dec_di_blk, hm_n);
    if ((dec_mode == BF_IBEX) || (dec_mode == BF_IBEX_RTL_CN)) {
      printf("[DVC IBEX PROOF] NOTE: BF_IBEX path reads corrupted_codeword "
             "(populated by ldpc_ibex_input from rx_blk/det_blk), not "
             "dec_di_blk; the dec_di_blk view above is the stale "
             "pre-dispatch content and is NOT what BF_IBEX decodes.\n");
    }
  }

  if (dec_mode == SKIP)
    ldpc_dec_skip();
  else if (dec_mode == BF_P0)
    ldpc_dec_bf(0, col_skip_itr);
  else if (dec_mode == BF_P3)
    ldpc_dec_bf(3, col_skip_itr);
  else if (dec_mode == BF_G2)
    ldpc_dec_bf2(3, col_skip_itr);
  else if (dec_mode == LAYER)
    ldpc_dec_layer();
  else if (dec_mode == LAYER_G2)
    ldpc_dec_layer2();
  else if (dec_mode == BF_IBEX)
    ldpc_dec_bf_ibex(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);
  else if (dec_mode == BF_IBEX_RTL_CN)
    ldpc_dec_bf_ibex_rtl_cn(ldpc_decoder_input, ldpc_decoder_parameters,
                            h_matrix);

  // remove padding
  // NOTE:
  // - `dec_do_blk` is always the full-length codeword view (hm_n bits, includes
  // pad bits and shortened parity tail bits).
  // - `dec_blk` is the DV-visible transmitted block (blk_len bits, excludes pad
  // bits and shortened parity tail bits). For IBEX (cir_sz==512) with
  // shortening, `blk_len - info_len` is smaller than `hm_m`, so we must not
  // copy `hm_m` parity bits into `dec_blk` (would overflow and corrupt memory,
  // causing seemingly random decoded output).
  vec_copy(dec_do_blk, dec_blk, 0, 0, info_len);
  const int parity_bits_in_det = blk_len - info_len;
  if ((cir_sz == 512) && (h_matrix.extra_bits_of_parity > 0)) {
    int first_parity_bits = h_matrix.extra_bits_of_parity;
    if (first_parity_bits > parity_bits_in_det)
      first_parity_bits = parity_bits_in_det;
    vec_copy(dec_do_blk, dec_blk, hm_k, info_len, first_parity_bits);
    const int remaining_parity_bits = parity_bits_in_det - first_parity_bits;
    if (remaining_parity_bits > 0) {
      vec_copy(dec_do_blk, dec_blk, hm_k + cir_sz, info_len + first_parity_bits,
               remaining_parity_bits);
    }
  } else {
    vec_copy(dec_do_blk, dec_blk, hm_k, info_len, parity_bits_in_det);
  }

  // check error bit number
  dec_err_num = 0;
  cor_err_num = 0;
  for (int i = 0; i < blk_len; i++) {
    if (dec_blk[i] != tx_blk[i])
      dec_err_num++;
    if (dec_blk[i] != det_blk[i])
      cor_err_num++;
  }

  // mis-correction case
  if ((cw_fail == 0) && (dec_err_num != 0))
    cw_miscorr = 1;
  else
    cw_miscorr = 0;

#ifdef _LDPC_DEBUG
  if (cw_fail == 0) {
    if (dec_err_num == 0)
      printf("[LDPC DEBUG] Decoding success, no error bits found.\n");
    else
      printf("[LDPC DEBUG] MIS-CORRECTION!!! (%d)\n", dec_err_num);
  } else {
    printf("[LDPC DEBUG] Decoding failed, (%d)!\n", dec_err_num);
  }
#endif
}

void ldpc_packet::ldpc_dec_bf(int p_num, int col_skip_itr) {
  mod2entry *e;
  int synd_wt;
  int col_updt;
  int itr_updt;
  bool col_skip;
  char *cn_synd_mem; // syndrome memory in CN order
  char *cn_synd_sel; // selected syndrome in CN order
  char *vn_synd_sel; // selected syndrome in VN order
  char *vn_synd_cnt; // syndrome weight of select columns in VN order
  char
      *vn_hd_sel; // current HD of selected column (from dec_do_blk) in VN order
  char *vn_raw_sel; // raw data of selected column (from dec_di_blk) in VN order
  char **vn_flp_sel; // flip flag of selected column in VN order
  int *vn_flp_col;   // column index of the vn_flp_sel
  int *vn_flp_itr;   // iteration of the vn_flp_sel
  char *cn_flp_sel;  // flip flag of selected column in CN order
  char *cn_synd_new; // new syndrome in CN order
  char *cn_synd_old;

  // allocate memory
  cn_synd_mem = (char *)calloc(hm_m, sizeof(*cn_synd_mem));
  cn_synd_sel = (char *)calloc(cir_sz, sizeof(*cn_synd_sel));
  vn_synd_sel = (char *)calloc(cir_sz, sizeof(*vn_synd_sel));
  vn_synd_cnt = (char *)calloc(cir_sz, sizeof(*vn_synd_cnt));
  vn_hd_sel = (char *)calloc(cir_sz, sizeof(*vn_hd_sel));
  vn_raw_sel = (char *)calloc(cir_sz, sizeof(*vn_raw_sel));
  vn_flp_sel = (char **)calloc(p_num + 1, sizeof(*vn_flp_sel));
  for (int i = 0; i <= p_num; i++)
    vn_flp_sel[i] = (char *)calloc(cir_sz, sizeof(*vn_flp_sel[i]));
  vn_flp_col = (int *)calloc(p_num + 1, sizeof(*vn_flp_col));
  vn_flp_itr = (int *)calloc(p_num + 1, sizeof(*vn_flp_itr));
  cn_flp_sel = (char *)calloc(cir_sz, sizeof(*cn_flp_sel));
  cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));
  cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));

  cw_fail = 1;
  cw_miscorr = 0;
  fina_synd_wt = 0;

  vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);
  vec_clr(cn_synd_mem, hm_m);
  for (int i = 0; i < p_num; i++)
    vec_clr(vn_flp_sel[i], cir_sz);
  fdec_cyc_num = 0;
  fdec_cyc_org = 0;

  for (int itr = 0;
       (itr <= fdec_max_itr) && ((fdec_early_term_en == 0) || (cw_fail == 1));
       itr++) {
    for (int i = 0; i < (itr == fdec_max_itr ? p_num : bm_n); i++) {
      if (itr == 0) {
        vec_copy(dec_di_blk, vn_flp_sel[p_num], i * cir_sz, 0, cir_sz);
        col_updt = i;
        itr_updt = itr;
        col_skip = false;
      } else {
        vec_clr(vn_synd_cnt, cir_sz);
        for (e = mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e);
             e = mod2sparse_next_in_col(e)) {
          // read syndrome
          vec_copy(cn_synd_mem, cn_synd_sel, e->row * cir_sz, 0, cir_sz);
          // barrel shift
          vec_shift(cn_synd_sel, vn_synd_sel, cir_sz, e->shift);
          // increment vn_synd_cnt
          vec_incr(vn_synd_cnt, vn_synd_sel, cir_sz);
        }

        // previous column is skipped
        if (col_skip)
          col_skip = false; // Column skip feature OFF
        else if (col_skip_itr == 0)
          col_skip = false; // non-skip iterations
        else if ((col_skip_itr > 0) && (itr < col_skip_itr))
          col_skip = false;
        else {
          col_skip = true;
          // make a skip decision
          for (int j = 0; j < cir_sz; j++) {
            if (vn_synd_cnt[j] >= flp_thrshd1[itr - 1])
              col_skip = false;
          }
        }

        // flip logics
        if (col_skip == false) {
          // read raw and current HD
          vec_copy(dec_di_blk, vn_raw_sel, i * cir_sz, 0, cir_sz);
          vec_copy(dec_do_blk, vn_hd_sel, i * cir_sz, 0, cir_sz);

          // pipelines
          for (int j = p_num; j > 0; j--)
            vec_copy(vn_flp_sel[j - 1], vn_flp_sel[j], 0, 0, cir_sz);
          for (int j = p_num; j > 0; j--) {
            vn_flp_col[j] = vn_flp_col[j - 1];
            vn_flp_itr[j] = vn_flp_itr[j - 1];
          }
          vn_flp_col[0] = i;
          vn_flp_itr[0] = itr;

          for (int j = 0; j < cir_sz; j++) {
            if (((vn_raw_sel[j] == vn_hd_sel[j]) &&
                 (vn_synd_cnt[j] >= flp_thrshd0[itr - 1])) ||
                ((vn_raw_sel[j] != vn_hd_sel[j]) &&
                 (vn_synd_cnt[j] >= flp_thrshd1[itr - 1]))) {
              vn_flp_sel[0][j] = 1;
            } else {
              vn_flp_sel[0][j] = 0;
            }
          }

          col_updt = vn_flp_col[p_num];
          itr_updt = vn_flp_itr[p_num];

          for (int j = 0; j < cir_sz; j++) {
            if (vn_flp_sel[p_num][j] == 1) {
              dec_do_blk[col_updt * cir_sz + j] =
                  (dec_do_blk[col_updt * cir_sz + j] + 1) % 2;
            }
          }
        } // non-skipped columns(flip logic)
      }   // non-1st iteration columns
      fdec_cyc_org++;

      if (col_skip == false) {
        fdec_cyc_num++;

        for (e = mod2sparse_first_in_col(qc_bm, col_updt);
             !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e)) {
          // barrel shift
          vec_shift(vn_flp_sel[p_num], cn_flp_sel, cir_sz, -1 * e->shift);
          // read old syndrome
          vec_copy(cn_synd_mem, cn_synd_old, e->row * cir_sz, 0, cir_sz);
          // update new syndrome
          vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
          // update syndrome memory
          vec_copy(cn_synd_new, cn_synd_mem, 0, e->row * cir_sz, cir_sz);
        }

        if ((itr > 0) || (i == (bm_n - 1))) {
          synd_wt = vec_sum(cn_synd_mem, hm_m);
          if (synd_wt == 0) {
            cw_fail = 0;
            cnvg_itr = itr_updt;
            cnvg_lyr = col_updt;
            if (fdec_early_term_en == 1)
              break;
          }
        }
      }
    }
  }

  if ((fdec_early_term_en == 0) || (cw_fail == 1)) {
    cnvg_itr = fdec_max_itr - 1;
    cnvg_lyr = bm_n - 1;
    fina_synd_wt = vec_sum(cn_synd_mem, hm_m);
  }

  free(cn_synd_mem);
  free(cn_synd_sel);
  free(vn_synd_sel);
  free(vn_synd_cnt);
  free(vn_hd_sel);
  free(vn_raw_sel);
  for (int i = 0; i <= p_num; i++)
    free(vn_flp_sel[i]);
  free(vn_flp_sel);
  free(vn_flp_col);
  free(vn_flp_itr);
  free(cn_flp_sel);
  free(cn_synd_new);
  free(cn_synd_old);
}

void ldpc_packet::ldpc_dec_bf2(int p_num, int col_skip_itr) {
  mod2entry *e;
  int synd_wt;
  int col_updt;
  int itr_updt;

  char *dec_sb_blk;  // CW soft bi
  char *cn_synd_mem; // syndrome memory in CN order
  char *cn_synd_sel; // selected syndrome in CN order
  char *vn_synd_sel; // selected syndrome in VN order
  char *vn_synd_cnt; // syndrome weight of select columns in VN order
  char *vn_sb_sel; // soft bit of selected column (from dec_sb_blk) in VN order
  char
      *vn_hd_sel; // current HD of selected column (from dec_do_blk) in VN order
  char *vn_raw_sel; // raw data of selected column (from dec_di_blk) in VN order
  char *vn_flp_sel; // flip flag of selected column in VN order
  char *cn_flp_sel; // flip flag of selected column in CN order
  char *cn_synd_new; // new syndrome in CN order
  char *cn_synd_old;

  // allocate memory
  dec_sb_blk = (char *)calloc(hm_n, sizeof(*dec_sb_blk));
  cn_synd_mem = (char *)calloc(hm_m, sizeof(*cn_synd_mem));
  cn_synd_sel = (char *)calloc(cir_sz, sizeof(*cn_synd_sel));
  vn_synd_sel = (char *)calloc(cir_sz, sizeof(*vn_synd_sel));
  vn_synd_cnt = (char *)calloc(cir_sz, sizeof(*vn_synd_cnt));
  vn_sb_sel = (char *)calloc(cir_sz, sizeof(*vn_sb_sel));
  vn_hd_sel = (char *)calloc(cir_sz, sizeof(*vn_hd_sel));
  vn_raw_sel = (char *)calloc(cir_sz, sizeof(*vn_raw_sel));
  vn_flp_sel = (char *)calloc(cir_sz, sizeof(*vn_flp_sel));
  cn_flp_sel = (char *)calloc(cir_sz, sizeof(*cn_flp_sel));
  cn_synd_old = (char *)calloc(cir_sz, sizeof(*cn_synd_old));
  cn_synd_new = (char *)calloc(cir_sz, sizeof(*cn_synd_new));

  // Initialize decoder
  cw_fail = 1;
  cw_miscorr = 0;
  fina_synd_wt = 0;
  vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);
  vec_clr(cn_synd_mem, hm_m);
  vec_clr(vn_flp_sel, cir_sz);
  vec_set(dec_sb_blk, hm_n);

  for (int itr = 0;
       (itr <= fdec_max_itr) && ((fdec_early_term_en == 0) || (cw_fail == 1));
       itr++) {
    for (int i = 0; i < bm_n; i++) {
      if (itr == 0) {
        vec_copy(dec_di_blk, vn_flp_sel, i * cir_sz, 0, cir_sz);
        col_updt = i;
        itr_updt = itr;
      } else {
        vec_clr(vn_synd_cnt, cir_sz);
        for (e = mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e);
             e = mod2sparse_next_in_col(e)) {
          vec_copy(cn_synd_mem, cn_synd_sel, e->row * cir_sz, 0, cir_sz);
          vec_shift(cn_synd_sel, vn_synd_sel, cir_sz, e->shift);
          vec_incr(vn_synd_cnt, vn_synd_sel, cir_sz);
        }

        // read raw and current HD and sb
        vec_copy(dec_di_blk, vn_raw_sel, i * cir_sz, 0, cir_sz);
        vec_copy(dec_do_blk, vn_hd_sel, i * cir_sz, 0, cir_sz);
        vec_copy(dec_sb_blk, vn_sb_sel, i * cir_sz, 0, cir_sz);

        for (int j = 0; j < cir_sz; j++) {
          if (vn_raw_sel[j] == vn_hd_sel[j]) // same
          {
            if (vn_sb_sel[j] == 1) // strong
            {
              if (vn_synd_cnt[j] >= flp_thrshd0_s[itr - 1])
                vn_flp_sel[j] = 1;
              else
                vn_flp_sel[j] = 0;

              if ((vn_synd_cnt[j] >= sb_thrshd0_s0[itr - 1]) &&
                  (vn_synd_cnt[j] < sb_thrshd0_s1[itr - 1]))
                vn_sb_sel[j] = 0;
            } else // weak
            {
              if (vn_synd_cnt[j] >= flp_thrshd0_w[itr - 1])
                vn_flp_sel[j] = 1;
              else
                vn_flp_sel[j] = 0;

              if ((vn_synd_cnt[j] < sb_thrshd0_w0[itr - 1]) ||
                  (vn_synd_cnt[j] >= sb_thrshd0_w1[itr - 1]))
                vn_sb_sel[j] = 1;
            }
          } else {
            if (vn_sb_sel[j] == 1) // strong
            {
              if (vn_synd_cnt[j] >= flp_thrshd1_s[itr - 1])
                vn_flp_sel[j] = 1;
              else
                vn_flp_sel[j] = 0;

              if ((vn_synd_cnt[j] >= sb_thrshd1_s0[itr - 1]) &&
                  (vn_synd_cnt[j] < sb_thrshd1_s1[itr - 1]))
                vn_sb_sel[j] = 0;
            } else // weak
            {
              if (vn_synd_cnt[j] >= flp_thrshd1_w[itr - 1])
                vn_flp_sel[j] = 1;
              else
                vn_flp_sel[j] = 0;

              if ((vn_synd_cnt[j] < sb_thrshd1_w0[itr - 1]) ||
                  (vn_synd_cnt[j] >= sb_thrshd1_w1[itr - 1]))
                vn_sb_sel[j] = 1;
            }
          }
        }

        itr_updt = itr;
        col_updt = i;

        for (int j = 0; j < cir_sz; j++) {
          if (vn_flp_sel[j] == 1) {
            dec_do_blk[col_updt * cir_sz + j] =
                (dec_do_blk[col_updt * cir_sz + j] + 1) % 2;
          }
        }

        vec_copy(vn_sb_sel, dec_sb_blk, 0, i * cir_sz, cir_sz);
      }

      // update syndrome memory
      // multiple times (col_wt)
      for (e = mod2sparse_first_in_col(qc_bm, col_updt); !mod2sparse_at_end(e);
           e = mod2sparse_next_in_col(e)) {
        // barrel shift
        vec_shift(vn_flp_sel, cn_flp_sel, cir_sz, -1 * e->shift);
        // read old syndrome
        vec_copy(cn_synd_mem, cn_synd_old, e->row * cir_sz, 0, cir_sz);
        // update new syndrome
        vec_mod2_add(cn_flp_sel, cn_synd_old, cn_synd_new, cir_sz);
        // update syndrome memory
        vec_copy(cn_synd_new, cn_synd_mem, 0, e->row * cir_sz, cir_sz);
      }

      if ((itr > 0) || (i == (bm_n - 1))) {
        synd_wt = vec_sum(cn_synd_mem, hm_m);
        if (synd_wt == 0) {
          cw_fail = 0;
          cnvg_itr = itr_updt;
          cnvg_lyr = col_updt;
          if (fdec_early_term_en == 1)
            break;
        }
      }
    }
  }

  if ((fdec_early_term_en == 0) || (cw_fail == 1)) {
    cnvg_itr = fdec_max_itr - 1;
    cnvg_lyr = bm_n - 1;
    fina_synd_wt = vec_sum(cn_synd_mem, hm_m);
  }

  free(dec_sb_blk);
  free(cn_synd_mem);
  free(cn_synd_sel);
  free(vn_synd_sel);
  free(vn_synd_cnt);
  free(vn_sb_sel);
  free(vn_hd_sel);
  free(vn_raw_sel);
  free(vn_flp_sel);
  free(cn_flp_sel);
  free(cn_synd_new);
  free(cn_synd_old);
} // ldpc_dec_bf2

// gen3 layer decoder
void ldpc_packet::ldpc_dec_layer() {
#ifdef _LDPC_DEBUG_DUMP
  FILE *cfp, *sfp, *hdfp, *lfp;
  int stmp, vtmp;
  dvc_ensure_output_c_code_dir();
  char cmem_dump[50] = "./output/c_code/rdec_cmem_dump.txt";
  char stot_dump[50] = "./output/c_code/rdec_stot_dump.txt";
  char hdmem_dump[50] = "./output/c_code/rdec_hdmem_dump.txt";
  char log_dump[50] = "./output/c_code/rdec_log_dump.txt";
  cfp = fopen(cmem_dump, "w");
  sfp = fopen(stot_dump, "w");
  hdfp = fopen(hdmem_dump, "w");
  lfp = fopen(log_dump, "w");
#endif

  mod2entry *e, *e_pre;
  char *dec_init;
  int shift_val1;
  int shift_val2;
  int cir_cnt;
  int sign_tmp;
  float val_tmp;
  int hd_init;

  struct cn_msg **cn_c_mem;
  struct cn_msg *cn_c_updt_cur; // current layer check node msg to be updt
  struct cn_msg *cn_c_sel_cur;  // current layer check node msg
  struct cn_msg *cn_c_sel_pre;  // previous layer check node msg
  float **cn_q_mem;             // Q mem in CN order of previous layer
  float *cn_q_sel_pre; // Q msg of the select circulant from previous layer
  float *cn_r_new_pre; // New R msg in CN order of previous layer
  float *cn_app_pre;   // APP = Q + R_new in CN order of previous layer
  float *cn_app_cur;   // APP = Q + R_new in CN order of current layer
  float *cn_q_sel_cur; // Q msg of the select circulant from current layer
  float *cn_r_old_cur; // old R msg in CN order of previous layer
  float
      *cn_q_updt_cur; // Updated Q msg of the select circulant in current layer
  int **cn_q_sign;    // Q sign

  char *layer_synd;
  char *cn_dec_hd;
  char *vn_dec_hd;
  int hd_updated;
  int layer_synd_wt;
  int synd_pass_cnt = 0;
  int hd_stable_cnt = 0;

  // allocation
  dec_init = (char *)calloc(bm_n, sizeof(*dec_init));
  vec_set(dec_init, bm_n);

  cn_c_mem = (struct cn_msg **)calloc(bm_m, sizeof(*cn_c_mem));
  for (int i = 0; i < bm_m; i++)
    cn_c_mem[i] = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_mem[i]));
  cn_c_updt_cur = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_updt_cur));

  cn_q_mem = (float **)calloc(bm_n, sizeof(*cn_q_mem));
  for (int i = 0; i < bm_n; i++)
    cn_q_mem[i] = (float *)calloc(cir_sz, sizeof(*cn_q_mem[i]));

  cn_r_new_pre = (float *)calloc(cir_sz, sizeof(*cn_r_new_pre));
  cn_app_pre = (float *)calloc(cir_sz, sizeof(*cn_app_pre));
  cn_app_cur = (float *)calloc(cir_sz, sizeof(*cn_app_cur));
  cn_q_sel_cur = (float *)calloc(cir_sz, sizeof(*cn_q_sel_cur));
  cn_r_old_cur = (float *)calloc(cir_sz, sizeof(*cn_r_old_cur));
  cn_q_updt_cur = (float *)calloc(cir_sz, sizeof(*cn_q_updt_cur));

  cn_q_sign = (int **)calloc(bm_n * col_wt, sizeof(*cn_q_sign));
  for (int i = 0; i < bm_n * col_wt; i++)
    cn_q_sign[i] = (int *)calloc(cir_sz, sizeof(*cn_q_sign[i]));

  layer_synd = (char *)calloc(cir_sz, sizeof(*layer_synd));
  vn_dec_hd = (char *)calloc(cir_sz, sizeof(*vn_dec_hd));
  cn_dec_hd = (char *)calloc(cir_sz, sizeof(*cn_dec_hd));

  // initialize decoder
  cw_fail = 1;
  cw_miscorr = 0;
  vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);

  ldpc_apply_sdlite_llr_override();

  {
    char *hard_init = (char *)calloc(hm_n, sizeof(*hard_init));
    if (hard_init) {
      dvc_bins_to_hard_bits(dec_di_blk, hard_init, hm_n, llr_tbl, bin_num);
      const int init_synd_qc =
          dvc_qc_syndrome_weight(qc_bm, hard_init, bm_m, cir_sz);
      const int init_synd_ibex = dvc_ibex_syndrome_weight(this, hard_init);
      init_synd_wt = std::max(init_synd_qc, init_synd_ibex);
      init_synd_wt_min = std::min(init_synd_wt_min, init_synd_wt);
      init_synd_wt_max = std::max(init_synd_wt_max, init_synd_wt);
      free(hard_init);
    }
  }

  for (int i = 0; i < bm_n; i++)
    for (int j = 0; j < cir_sz; j++)
      cn_q_mem[i][j] = (float)llr_tbl[dec_di_blk[i * cir_sz + j]];

  // iterative decoding
  for (int itr = 0;
       (itr <= ldec_max_itr) && ((ldec_early_term_en == 0) || (cw_fail == 1));
       itr++) {
    // Q sign mem index
    cir_cnt = 0;

    // layer decoding
    for (int layer = 0;
         layer < bm_m && ((ldec_early_term_en == 0) || (cw_fail == 1));
         layer++) {
      // initilize HD mem
      hd_init = (vec_sum(dec_init, bm_n) != 0);

#ifdef _LDPC_DEBUG_DUMP
      printf("[LDPC DEBUG] Layer decoding @ iteration %d, layer %d ...\n", itr,
             layer);
#endif

      // init current layer C-MSG
      // C-MSG of previous iteration

      cn_c_sel_cur = cn_c_mem[layer];
      // C-MSG to be updt
      for (int i = 0; i < cir_sz; i++) {
        cn_c_updt_cur[i].min1_val = 100000;
        cn_c_updt_cur[i].min2_val = 100000;
        cn_c_updt_cur[i].min1_pos = 0;
        cn_c_updt_cur[i].sign_tot = 1;
      }

      // Earlier termination init
      hd_updated = 0;
      vec_clr(layer_synd, cir_sz);

      // Per circulant of the layer
      for (e = mod2sparse_first_in_row(qc_bm, layer);
           !mod2sparse_at_end(e) &&
           ((ldec_early_term_en == 0) || (cw_fail == 1));
           e = mod2sparse_next_in_row(e)) {
        // read Q from the previous layer of the selected column
        cn_q_sel_pre = cn_q_mem[e->col];

        // find out the previous layer of select column
        e_pre = mod2sparse_prev_in_col(e);
        if (mod2sparse_at_end(e_pre))
          e_pre = mod2sparse_last_in_col(qc_bm, e->col);

        cn_c_sel_pre = cn_c_mem[e_pre->row]; // read previous layer C msg

        // cal Rnew and APP
        for (int i = 0; i < cir_sz; i++) {
          // Qmsg sign
          sign_tmp = (cn_q_sel_pre[i] >= 0) ? 1 : -1;

          // Rnew
          if (cn_c_sel_pre[i].min1_pos == e->col)
            cn_r_new_pre[i] =
                cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
          else
            cn_r_new_pre[i] =
                cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;

          // APP in CN order of previous layer
          cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];

          // Quantization
          if (finite_mode == 1) {
            cn_app_pre[i] =
                (float)Sat_Quan((double)cn_app_pre[i], finite_q_max,
                                finite_q_min, finite_q_num, finite_f_num);
          }
        }

        // APP shift
        // when decoder initilized, Q msg are in VN order
        if (dec_init[e->col] == 1) {
          shift_val1 = e->shift;
          shift_val2 = 0;
          dec_init[e->col] = 0;
        } else {
          shift_val1 = -1 * e_pre->shift + e->shift;
          shift_val2 = -1 * e_pre->shift;
        }

        for (int i = 0; i < cir_sz; i++) {
          cn_app_cur[i] = cn_app_pre[(i + shift_val1 + cir_sz) % cir_sz];
          vn_dec_hd[i] =
              cn_app_pre[(i + shift_val2 + cir_sz) % cir_sz] >= 0 ? 0 : 1;
        }

        // CW converge check logic per circulant
        // 1. check if HD updated
        if (hd_updated == 0)
          if (vec_cmp(dec_do_blk, vn_dec_hd, e->col * cir_sz, 0, cir_sz) == 1)
            hd_updated = 1;

        vec_copy(vn_dec_hd, dec_do_blk, 0, e->col * cir_sz, cir_sz);

        // 2. accumulate syndrome
        vec_shift(vn_dec_hd, cn_dec_hd, cir_sz, -1 * e->shift);
        vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);

        // calculate R_old and current Q, update current layer C and Q
        for (int i = 0; i < cir_sz; i++) {
          if (cn_c_sel_cur[i].min1_pos == e->col)
            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val *
                              cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
          else
            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val *
                              cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

          // Q -= Rold
          cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          // Quantization
          if (finite_mode == 1) {
            cn_q_updt_cur[i] =
                (float)Sat_Quan((double)cn_q_updt_cur[i], finite_q_max,
                                finite_q_min, finite_q_num, finite_f_num);
          }

          // update C
          sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
          val_tmp = cn_q_updt_cur[i] * sign_tmp;
          cn_c_updt_cur[i].sign_tot *= sign_tmp;

          if (val_tmp < cn_c_updt_cur[i].min1_val) {
            cn_c_updt_cur[i].min2_val = cn_c_updt_cur[i].min1_val;
            cn_c_updt_cur[i].min1_val = val_tmp;
            cn_c_updt_cur[i].min1_pos = e->col;
          } else if (val_tmp < cn_c_updt_cur[i].min2_val) {
            cn_c_updt_cur[i].min2_val = val_tmp;
          }

          cn_q_sign[cir_cnt][i] = sign_tmp;
        }

        // update Q memory
        for (int i = 0; i < cir_sz; i++) {
          cn_q_mem[e->col][i] = cn_q_updt_cur[i];
        }

        cir_cnt++;
      } // per circulant

      // update C_MSG per layer
      for (int i = 0; i < cir_sz; i++) {
        const float min1_pre_q = cn_c_updt_cur[i].min1_val * alpha;
        const float min2_pre_q = cn_c_updt_cur[i].min2_val * alpha;
        cn_c_mem[layer][i].min1_val = min1_pre_q;
        cn_c_mem[layer][i].min2_val = min2_pre_q;

        if (finite_mode == 1) {
          cn_c_mem[layer][i].min1_val =
              (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, finite_c_max,
                              finite_c_min, finite_c_num, finite_f_num);
          cn_c_mem[layer][i].min2_val =
              (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, finite_c_max,
                              finite_c_min, finite_c_num, finite_f_num);
        }

        cn_c_mem[layer][i].min1_pos = cn_c_updt_cur[i].min1_pos;
        cn_c_mem[layer][i].sign_tot = cn_c_updt_cur[i].sign_tot;
      }

      // check converage checking
      layer_synd_wt = vec_sum(layer_synd, cir_sz);
      if (hd_init == 1) {
        hd_stable_cnt = 0;
        synd_pass_cnt = 0;
      } else if ((hd_updated == 0) && (layer_synd_wt == 0)) {
        hd_stable_cnt++;
        synd_pass_cnt++;
      } else {
        hd_stable_cnt = 0;
        synd_pass_cnt = 0;
      }

      if ((synd_pass_cnt >= bm_m) && (hd_stable_cnt >= bm_m - 1)) {
        cw_fail = 0;
        cnvg_itr = itr;
        cnvg_lyr = layer;
      }
    }
  }

  if ((cw_fail == 1) || (ldec_early_term_en == 0)) {
    cnvg_itr = ldec_max_itr - 1;
    cnvg_lyr = bm_m - 1;
  }

  {
    const int fina_synd_qc =
        dvc_qc_syndrome_weight(qc_bm, dec_do_blk, bm_m, cir_sz);
    const int fina_synd_ibex = dvc_ibex_syndrome_weight(this, dec_do_blk);
    fina_synd_wt = std::max(fina_synd_qc, fina_synd_ibex);
  }

  // free all
  free(dec_init);
  for (int i = 0; i < bm_m; i++)
    free(cn_c_mem[i]);
  free(cn_c_mem);
  free(cn_c_updt_cur);
  for (int i = 0; i < bm_n; i++)
    free(cn_q_mem[i]);
  free(cn_q_mem);
  free(cn_r_new_pre);
  free(cn_app_pre);
  free(cn_app_cur);
  free(cn_q_sel_cur);
  free(cn_r_old_cur);
  free(cn_q_updt_cur);
  for (int i = 0; i < bm_n * col_wt; i++)
    free(cn_q_sign[i]);
  free(cn_q_sign);
  free(layer_synd);
  free(cn_dec_hd);
  free(vn_dec_hd);
} // ldpc_dec_layer

void ldpc_packet::ldpc_dec_layer2() {
#ifdef _ldpc_debug_dump
  file *cfp, *sfp, *hdfp, *lfp;
  int stmp, vtmp;
  dvc_ensure_output_c_code_dir();
  char cmem_dump[50] = "./output/c_code/rdec_cmem_dump.txt";
  char stot_dump[50] = "./output/c_code/rdec_stot_dump.txt";
  char hdmem_dump[50] = "./output/c_code/rdec_hdmem_dump.txt";
  char log_dump[50] = "./output/c_code/rdec_log_dump.txt";
  cfp = fopen(cmem_dump, "w");
  sfp = fopen(stot_dump, "w");
  hdfp = fopen(hdmem_dump, "w");
  lfp = fopen(log_dump, "w");

  // one-time probe to help align dump vs rtl when a specific column seems
  // missing. keep it short: list whether target_col exists in each layer's
  // schedule.
  if (lfp) {
    const int target_col = 73;
    fprintf(lfp, "[dvc] qc_bm col probe: target_col=%d bm_m=%d bm_n=%d\n",
            target_col, bm_m, bm_n);
    for (int rr = 0; rr < bm_m; rr++) {
      int found = 0;
      for (mod2entry *ee = mod2sparse_first_in_row(qc_bm, rr);
           !mod2sparse_at_end(ee); ee = mod2sparse_next_in_row(ee)) {
        if (ee->col == target_col) {
          fprintf(lfp, "[dvc] layer=%d has col=%d shift=%d\n", rr, ee->col,
                  ee->shift);
          found = 1;
          break;
        }
      }
      if (!found)
        fprintf(lfp, "[dvc] layer=%d missing col=%d\n", rr, target_col);
    }
    fflush(lfp);
  }
#endif

  mod2entry *e, *e_pre;
  char *dec_init;
  int shift_val1;
  int shift_val2;
  int cir_cnt;
  int sign_tmp;
  float val_tmp;
  int hd_init;

  struct cn_msg **cn_c_mem;
  struct cn_msg *cn_c_updt_cur; // current layer check node msg to be updt
  struct cn_msg *cn_c_sel_cur;  // current layer check node msg
  struct cn_msg *cn_c_sel_pre;  // previous layer check node msg
  float **cn_q_mem;             // q mem in cn order of previous layer
  float *cn_q_sel_pre; // q msg of the select circulant from previous layer
  float *cn_r_new_pre; // new r msg in cn order of previous layer
  float *cn_app_pre;   // app = q + r_new in cn order of previous layer
  float *cn_app_cur;   // app = q + r_new in cn order of current layer
  float *cn_q_sel_cur; // q msg of the select circulant from current layer
  float *cn_r_old_cur; // old r msg in cn order of previous layer
  float
      *cn_q_updt_cur; // updated q msg of the select circulant in current layer
  int **cn_q_sign;    // q sign

  char *layer_synd;
  char *cn_dec_hd;
  char *vn_dec_hd;
  int hd_updated;
  int layer_synd_wt;
  int synd_pass_cnt = 0;
  int hd_stable_cnt = 0;
#ifdef _ldpc_debug_dump
  int dbg_converged_printed = 0;
#endif

  // allocation
  dec_init = (char *)calloc(bm_n, sizeof(*dec_init));
  vec_set(dec_init, bm_n);

  cn_c_mem = (struct cn_msg **)calloc(bm_m, sizeof(*cn_c_mem));
  for (int i = 0; i < bm_m; i++)
    cn_c_mem[i] = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_mem[i]));
  cn_c_updt_cur = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_updt_cur));

  cn_q_mem = (float **)calloc(bm_n, sizeof(*cn_q_mem));
  for (int i = 0; i < bm_n; i++)
    cn_q_mem[i] = (float *)calloc(cir_sz, sizeof(*cn_q_mem[i]));

  cn_r_new_pre = (float *)calloc(cir_sz, sizeof(*cn_r_new_pre));
  cn_app_pre = (float *)calloc(cir_sz, sizeof(*cn_app_pre));
  cn_app_cur = (float *)calloc(cir_sz, sizeof(*cn_app_cur));
  cn_q_sel_cur = (float *)calloc(cir_sz, sizeof(*cn_q_sel_cur));
  cn_r_old_cur = (float *)calloc(cir_sz, sizeof(*cn_r_old_cur));
  cn_q_updt_cur = (float *)calloc(cir_sz, sizeof(*cn_q_updt_cur));

  cn_q_sign = (int **)calloc(total_cir, sizeof(*cn_q_sign));
  for (int i = 0; i < total_cir; i++)
    cn_q_sign[i] = (int *)calloc(cir_sz, sizeof(*cn_q_sign[i]));

  layer_synd = (char *)calloc(cir_sz, sizeof(*layer_synd));
  vn_dec_hd = (char *)calloc(cir_sz, sizeof(*vn_dec_hd));
  cn_dec_hd = (char *)calloc(cir_sz, sizeof(*cn_dec_hd));

  // initialize decoder
  cw_fail = 1;
  cw_miscorr = 0;
  vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);

  ldpc_apply_sdlite_llr_override();

  // {
  //   char *hard_init = (char *)calloc(hm_n, sizeof(*hard_init));
  //   if (hard_init) {
  //     dvc_bins_to_hard_bits(dec_di_blk, hard_init, hm_n, llr_tbl, bin_num);
  //     init_synd_wt = dvc_ibex_syndrome_weight(this, hard_init);
  //     init_synd_wt_min = std::min(init_synd_wt_min, init_synd_wt);
  //     init_synd_wt_max = std::max(init_synd_wt_max, init_synd_wt);
  //     free(hard_init);
  //   }
  // }


  for (int i = 0; i < bm_n; i++) {
    for (int j = 0; j < cir_sz; j++) {
      cn_q_mem[i][j] = (float)llr_tbl[dec_di_blk[i * cir_sz + j]];
      if (finite_mode == 1)
        cn_q_mem[i][j] =
            (float)Sat_Quan((double)cn_q_mem[i][j], finite_q_max, finite_q_min,
                            finite_q_num, finite_f_num);
    }
  }

  // iterative decoding starts from iter1 (iter0 is init syndrome only).
  for (int itr = 0;
       (itr <= ldec_max_itr) && ((ldec_early_term_en == 0) || (cw_fail == 1));
       itr++) {
    // q sign mem index
    cir_cnt = 0;

    // layer decoding
    for (int layer = 0;
         layer < bm_m && ((ldec_early_term_en == 0) || (cw_fail == 1));
         layer++) {
      int last_nz_col = -1;
      int last_nz_shift = -1;
      // initilize hd mem
      hd_init = (vec_sum(dec_init, bm_n) != 0);

#ifdef _ldpc_debug_dump
      printf("[ldpc debug] layer decoding @ iteration %d, layer %d ...\n", itr,
             layer);
#endif

      // init current layer c-msg
      // c-msg of previous iteration

      cn_c_sel_cur = cn_c_mem[layer];
      // c-msg to be updt
      for (int i = 0; i < cir_sz; i++) {
        cn_c_updt_cur[i].min1_val = 100000;
        cn_c_updt_cur[i].min2_val = 100000;
        cn_c_updt_cur[i].min1_pos = 0;
        cn_c_updt_cur[i].sign_tot = 1;
      }

      // earlier termination init
      hd_updated = 0;
      vec_clr(layer_synd, cir_sz);

      // per circulant of the layer
      for (e = mod2sparse_first_in_row(qc_bm, layer);
           !mod2sparse_at_end(e) &&
           ((ldec_early_term_en == 0) || (cw_fail == 1));
           e = mod2sparse_next_in_row(e)) {
        last_nz_col = e->col;
        last_nz_shift = e->shift;
        // read q from the previous layer of the selected column
        cn_q_sel_pre = cn_q_mem[e->col];

        // find out the previous layer of select column
        e_pre = mod2sparse_prev_in_col(e);
        if (mod2sparse_at_end(e_pre))
          e_pre = mod2sparse_last_in_col(qc_bm, e->col);

        cn_c_sel_pre = cn_c_mem[e_pre->row]; // read previous layer c msg

        // cal rnew and app
        for (int i = 0; i < cir_sz; i++) {
          // qmsg sign
          sign_tmp = (cn_q_sel_pre[i] >= 0) ? 1 : -1;

          // rnew
          if (cn_c_sel_pre[i].min1_pos == e->col)
            cn_r_new_pre[i] =
                cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
          else
            cn_r_new_pre[i] =
                cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;

          // app in cn order of previous layer
          if (h_matrix.extra_bytes_of_parity == 0) {
            cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
          } else {
            if (h_matrix.occupied[e_pre->row][e_pre->col] &&
                (e_pre->row < (h_matrix.rows - 1)))
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else if (h_matrix.occupied[e_pre->row][e_pre->col] &&
                     (e_pre->row == (h_matrix.rows - 1)) &&
                     h_matrix.mask[e_pre->col][(i + e_pre->shift) % cir_sz])
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else if (h_matrix.fade[e_pre->row][e_pre->col] &&
                     !h_matrix.mask[e_pre->col][(i + e_pre->shift) % cir_sz])
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else
              cn_app_pre[i] = cn_q_sel_pre[i];
          }

          // quantization
          if (finite_mode == 1) {
            cn_app_pre[i] =
                (float)Sat_Quan((double)cn_app_pre[i], finite_q_max,
                                finite_q_min, finite_q_num, finite_f_num);
          }
        }

        // app shift
        // when decoder initilized, q msg are in vn order
        if (dec_init[e->col] == 1) {
          shift_val1 = e->shift;
          shift_val2 = 0;
          dec_init[e->col] = 0;
        } else {
          shift_val1 = -1 * e_pre->shift + e->shift;
          shift_val2 = -1 * e_pre->shift;
        }

        for (int i = 0; i < cir_sz; i++) {
          cn_app_cur[i] = cn_app_pre[(i + shift_val1 + cir_sz) % cir_sz];
          vn_dec_hd[i] =
              cn_app_pre[(i + shift_val2 + cir_sz) % cir_sz] >= 0 ? 0 : 1;
        }

#ifdef _ldpc_debug_dump
        if (lfp) {
          fprintf(lfp, "itr%2d/layer%2d/col%2d: \n", itr, layer, e->col);
          const int dump_groups_8 = (cir_sz + 7) / 8;
          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "q pre msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_q_sel_pre[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "r new msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_r_new_pre[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "app-c msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_app_pre[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "app-s msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_app_cur[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }
        }
#endif

        // cw converge check logic per circulant
        // 1. check if hd updated
        if (hd_updated == 0)
          if (vec_cmp(dec_do_blk, vn_dec_hd, e->col * cir_sz, 0, cir_sz) == 1)
            hd_updated = 1;

#ifdef _ldpc_debug_dump
        if (lfp) {
          for (int i = 0; i < cir_sz; i++) {
            if (dec_do_blk[e->col * cir_sz + i] != vn_dec_hd[i]) {
              fprintf(lfp, "itr%2d/layer%2d/col%2d: flip bit %d (%d-->%d)\n",
                      itr, layer, e->col, i, dec_do_blk[e->col * cir_sz + i],
                      vn_dec_hd[i]);
            }
          }
        }
#endif
        vec_copy(vn_dec_hd, dec_do_blk, 0, e->col * cir_sz, cir_sz);

#ifdef _ldpc_debug_dump
        if (hdfp) {
          fprintf(hdfp, "itr%2d/layer%2d/col%2d: ", itr, layer, e->col);
          for (int i = cir_sz / 4 - 1; i >= 0; i--) {
            stmp = 0;
            for (int j = 3; j >= 0; j--)
              stmp = stmp * 2 + vn_dec_hd[i * 4 + j];
            fprintf(hdfp, "%lx", stmp);
          }
          fprintf(hdfp, "\n");
        }
#endif

        // 2. accumulate syndrome
        vec_shift(vn_dec_hd, cn_dec_hd, cir_sz, -1 * e->shift);
        if (h_matrix.extra_bytes_of_parity == 0) {
          vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);
        } else {
          if (h_matrix.occupied[e->row][e->col] &&
              (e->row == (h_matrix.rows - 1))) {
            for (int i = 0; i < cir_sz; i++) {
              if (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])
                cn_dec_hd[i] = 0;
            }
          }
          if (h_matrix.fade[e->row][e->col]) {
            for (int i = 0; i < cir_sz; i++) {
              if (h_matrix.mask[e->col][(i + e->shift) % cir_sz])
                cn_dec_hd[i] = 0;
            }
          }
          vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);
        }

        // calculate r_old and current q, update current layer c and q
        for (int i = 0; i < cir_sz; i++) {
          if (cn_c_sel_cur[i].min1_pos == e->col)
            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val *
                              cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
          else
            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val *
                              cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

          // q -= rold
          if (h_matrix.extra_bytes_of_parity == 0) {
            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          } else {
            if (h_matrix.occupied[e->row][e->col] &&
                (e->row == (h_matrix.rows - 1)) &&
                (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              cn_r_old_cur[i] = 0;
            }
            if (h_matrix.fade[e->row][e->col] &&
                (h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              cn_r_old_cur[i] = 0;
            }
            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          }

          // quantization
          if (finite_mode == 1) {
            cn_q_updt_cur[i] =
                (float)Sat_Quan((double)cn_q_updt_cur[i], finite_q_max,
                                finite_q_min, finite_q_num, finite_f_num);
          }

          // update c
          if (h_matrix.extra_bytes_of_parity == 0) {
            sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
            val_tmp = cn_q_updt_cur[i] * sign_tmp;
          } else {
            if (h_matrix.occupied[e->row][e->col] &&
                (e->row == (h_matrix.rows - 1)) &&
                (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              sign_tmp = 1;
              val_tmp = 100000;
            }
            if (h_matrix.fade[e->row][e->col] &&
                (h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              sign_tmp = 1;
              val_tmp = 100000;
            } else {
              sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
              val_tmp = cn_q_updt_cur[i] * sign_tmp;
            }
          }

          cn_c_updt_cur[i].sign_tot *= sign_tmp;

          if (val_tmp < cn_c_updt_cur[i].min1_val) {
            cn_c_updt_cur[i].min2_val = cn_c_updt_cur[i].min1_val;
            cn_c_updt_cur[i].min1_val = val_tmp;
            cn_c_updt_cur[i].min1_pos = e->col;
          } else if (val_tmp < cn_c_updt_cur[i].min2_val) {
            cn_c_updt_cur[i].min2_val = val_tmp;
          }

          cn_q_sign[cir_cnt][i] = sign_tmp;
        }

        // update q memory
        for (int i = 0; i < cir_sz; i++) {
          cn_q_mem[e->col][i] = cn_q_updt_cur[i];
        }

#ifdef _ldpc_debug_dump
        if (lfp) {
          const int dump_groups_8 = (cir_sz + 7) / 8;
          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "r old msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_r_old_cur[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "q new msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_q_updt_cur[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }
        }
#endif

        cir_cnt++;
      } // per circulant

      // update c_msg per layer
      for (int i = 0; i < cir_sz; i++) {
        const float min1_pre_q = cn_c_updt_cur[i].min1_val * alpha;
        const float min2_pre_q = cn_c_updt_cur[i].min2_val * alpha;
        cn_c_mem[layer][i].min1_val = min1_pre_q;
        cn_c_mem[layer][i].min2_val = min2_pre_q;

        if (finite_mode == 1) {
          cn_c_mem[layer][i].min1_val =
              (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, finite_c_max,
                              finite_c_min, finite_c_num, finite_f_num);
          cn_c_mem[layer][i].min2_val =
              (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, finite_c_max,
                              finite_c_min, finite_c_num, finite_f_num);
        }

        cn_c_mem[layer][i].min1_pos = cn_c_updt_cur[i].min1_pos;
        cn_c_mem[layer][i].sign_tot = cn_c_updt_cur[i].sign_tot;
#ifdef _ldpc_debug_dump
        const int cmem_dump_scale = 1 << finite_f_num;
        fprintf(cfp, "itr%d/l%d/c%d: min1 %x, min2 %x, min1 pos %d, sign_tot %d\n",
          itr, layer, i,
          int(cmem_dump_scale * cn_c_mem[layer][i].min1_val),
          int(cmem_dump_scale * cn_c_mem[layer][i].min2_val),
          cn_c_mem[layer][i].min1_pos,
          (1-cn_c_mem[layer][i].sign_tot)/2 );
#endif
      }

#ifdef _ldpc_debug_dump
      if (sfp) {
        for (int i = (cir_sz / 4 - 1); i >= 0; i--) {
          stmp = 0;
          for (int j = 3; j >= 0; j--)
            stmp = stmp * 2 + (1 - cn_c_mem[layer][i * 4 + j].sign_tot) / 2;
          fprintf(sfp, "%lx", stmp);
        }
        fprintf(sfp, "\n");
      }
#endif

      // check converage checking
      layer_synd_wt = vec_sum(layer_synd, cir_sz);
      int synd_pass_cnt_old = synd_pass_cnt;
      int hd_stable_cnt_old = hd_stable_cnt;
      int conv_now;

      if (itr == 0)
        init_synd_wt += layer_synd_wt;
      if (hd_init == 1) {
        hd_stable_cnt = 0;
        synd_pass_cnt = 0;
      } else if ((hd_updated == 0) && (layer_synd_wt == 0)) {
        hd_stable_cnt++;
        synd_pass_cnt++;
      } else {
        hd_stable_cnt = 0;
        synd_pass_cnt = 0;
      }
      conv_now = ((synd_pass_cnt >= bm_m) && (hd_stable_cnt >= bm_m - 1));

#ifdef _ldpc_debug_dump
      if (lfp) {
        fprintf(lfp,
                "[dvc] layer summary: itr=%d layer=%d tail_col=%d tail_shift=%d layer_synd_wt=%d hd_updated=%d synd_pass_cnt=%d->%d hd_stable_cnt=%d->%d conv_now=%d cw_fail=%d\n",
                itr, layer, last_nz_col, last_nz_shift, layer_synd_wt,
                hd_updated, synd_pass_cnt_old, synd_pass_cnt,
                hd_stable_cnt_old, hd_stable_cnt, conv_now, cw_fail);
      }
#endif

      if (conv_now) {
        cw_fail = 0;
        cnvg_itr = itr;
        cnvg_lyr = layer;
#ifdef _ldpc_debug_dump
        if (lfp && !dbg_converged_printed) {
          fprintf(lfp,
                  "[dvc] converged latch: itr=%d layer=%d synd_pass_cnt=%d hd_stable_cnt=%d (early_term=%d)\n",
                  itr, layer, synd_pass_cnt, hd_stable_cnt, ldec_early_term_en);
          fflush(lfp);
          dbg_converged_printed = 1;
        }
#endif
      }
    }
  }

  if ((cw_fail == 1) || (ldec_early_term_en == 0)) {
    cnvg_itr = ldec_max_itr;
    cnvg_lyr = bm_m - 1;
  }

  fina_synd_wt = dvc_ibex_syndrome_weight(this, dec_do_blk);

#ifdef _ldpc_debug_dump
  if (lfp) {
    fprintf(lfp,
            "[dvc] cw summary: init_synd_wt=%d final_synd_wt=%d init_ok=%d final_ok=%d cw_fail=%d cnvg_itr=%d cnvg_lyr=%d\n",
            init_synd_wt, fina_synd_wt, (init_synd_wt == 0),
            (fina_synd_wt == 0), cw_fail, cnvg_itr, cnvg_lyr);
    dvc_dump_col_major_codeword(lfp, dec_do_blk, bm_n, cir_sz, "dec_do");
    fflush(lfp);
  }
#endif

  // free all
  free(dec_init);
  for (int i = 0; i < bm_m; i++)
    free(cn_c_mem[i]);
  free(cn_c_mem);
  free(cn_c_updt_cur);
  for (int i = 0; i < bm_n; i++)
    free(cn_q_mem[i]);
  free(cn_q_mem);
  free(cn_r_new_pre);
  free(cn_app_pre);
  free(cn_app_cur);
  free(cn_q_sel_cur);
  free(cn_r_old_cur);
  free(cn_q_updt_cur);
  for (int i = 0; i < total_cir; i++)
    free(cn_q_sign[i]);
  free(cn_q_sign);
  free(layer_synd);
  free(cn_dec_hd);
  free(vn_dec_hd);

#ifdef _ldpc_debug_dump
  if (cfp)
    fclose(cfp);
  if (sfp)
    fclose(sfp);
  if (hdfp)
    fclose(hdfp);
  if (lfp)
    fclose(lfp);
#endif
} // ldpc_dec_layer2

void ldpc_packet::ldpc_dec_pms() {
#ifdef _ldpc_debug_dump
  file *cfp, *sfp, *hdfp, *lfp;
  int stmp, vtmp;
  dvc_ensure_output_c_code_dir();
  char cmem_dump[50] = "./output/c_code/rdec_cmem_dump.txt";
  char stot_dump[50] = "./output/c_code/rdec_stot_dump.txt";
  char hdmem_dump[50] = "./output/c_code/rdec_hdmem_dump.txt";
  char log_dump[50] = "./output/c_code/rdec_log_dump.txt";
  cfp = fopen(cmem_dump, "w");
  sfp = fopen(stot_dump, "w");
  hdfp = fopen(hdmem_dump, "w");
  lfp = fopen(log_dump, "w");

  // one-time probe to help align dump vs rtl when a specific column seems
  // missing. keep it short: list whether target_col exists in each layer's
  // schedule.
  if (lfp) {
    const int target_col = 73;
    fprintf(lfp, "[dvc] qc_bm col probe: target_col=%d bm_m=%d bm_n=%d\n",
            target_col, bm_m, bm_n);
    for (int rr = 0; rr < bm_m; rr++) {
      int found = 0;
      for (mod2entry *ee = mod2sparse_first_in_row(qc_bm, rr);
           !mod2sparse_at_end(ee); ee = mod2sparse_next_in_row(ee)) {
        if (ee->col == target_col) {
          fprintf(lfp, "[dvc] layer=%d has col=%d shift=%d\n", rr, ee->col,
                  ee->shift);
          found = 1;
          break;
        }
      }
      if (!found)
        fprintf(lfp, "[dvc] layer=%d missing col=%d\n", rr, target_col);
    }
    fflush(lfp);
  }
#endif

  mod2entry *e, *e_pre;
  char *dec_init;
  int shift_val1;
  int shift_val2;
  int cir_cnt;
  int sign_tmp;
  float val_tmp;
  int hd_init;

  struct cn_msg **cn_c_mem;
  struct cn_msg *cn_c_updt_cur; // current layer check node msg to be updt
  struct cn_msg *cn_c_sel_cur;  // current layer check node msg
  struct cn_msg *cn_c_sel_pre;  // previous layer check node msg
  float **cn_q_mem;             // q mem in cn order of previous layer
  float *cn_q_sel_pre; // q msg of the select circulant from previous layer
  float *cn_r_new_pre; // new r msg in cn order of previous layer
  float *cn_app_pre;   // app = q + r_new in cn order of previous layer
  float *cn_app_cur;   // app = q + r_new in cn order of current layer
  float *cn_q_sel_cur; // q msg of the select circulant from current layer
  float *cn_r_old_cur; // old r msg in cn order of previous layer
  float *cn_q_updt_cur; // updated q msg of the select circulant in current layer
  int **cn_q_sign;    // q sign
  float *min3_val;

  char *layer_synd;
  char *cn_dec_hd;
  char *vn_dec_hd;
  int hd_updated;
  int layer_synd_wt;
  int synd_pass_cnt = 0;
  int hd_stable_cnt = 0;
  int pl1, pl2;
#ifdef _ldpc_debug_dump
  int dbg_converged_printed = 0;
#endif

  // allocation
  dec_init = (char *)calloc(bm_n, sizeof(*dec_init));
  vec_set(dec_init, bm_n);

  cn_c_mem = (struct cn_msg **)calloc(bm_m, sizeof(*cn_c_mem));
  for (int i = 0; i < bm_m; i++)
    cn_c_mem[i] = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_mem[i]));
  cn_c_updt_cur = (struct cn_msg *)calloc(cir_sz, sizeof(*cn_c_updt_cur));

  cn_q_mem = (float **)calloc(bm_n, sizeof(*cn_q_mem));
  for (int i = 0; i < bm_n; i++)
    cn_q_mem[i] = (float *)calloc(cir_sz, sizeof(*cn_q_mem[i]));

  cn_r_new_pre = (float *)calloc(cir_sz, sizeof(*cn_r_new_pre));
  cn_app_pre = (float *)calloc(cir_sz, sizeof(*cn_app_pre));
  cn_app_cur = (float *)calloc(cir_sz, sizeof(*cn_app_cur));
  cn_q_sel_cur = (float *)calloc(cir_sz, sizeof(*cn_q_sel_cur));
  cn_r_old_cur = (float *)calloc(cir_sz, sizeof(*cn_r_old_cur));
  cn_q_updt_cur = (float *)calloc(cir_sz, sizeof(*cn_q_updt_cur));

  cn_q_sign = (int **)calloc(total_cir, sizeof(*cn_q_sign));
  for (int i = 0; i < total_cir; i++)
    cn_q_sign[i] = (int *)calloc(cir_sz, sizeof(*cn_q_sign[i]));

  min3_val = (float*)calloc(cir_sz, sizeof(*min3_val));

  layer_synd = (char *)calloc(cir_sz, sizeof(*layer_synd));
  vn_dec_hd = (char *)calloc(cir_sz, sizeof(*vn_dec_hd));
  cn_dec_hd = (char *)calloc(cir_sz, sizeof(*cn_dec_hd));

  // initialize decoder
  cw_fail = 1;
  cw_miscorr = 0;
  vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);

  ldpc_apply_sdlite_llr_override();

  {
    char *hard_init = (char *)calloc(hm_n, sizeof(*hard_init));
    if (hard_init) {
      dvc_bins_to_hard_bits(dec_di_blk, hard_init, hm_n, llr_tbl, bin_num);
#ifdef _LDPC_DEBUG_DUMP
      {
        static unsigned int pms_hard_cw_seq = 0;
        char hard_dump_path[128];
        snprintf(hard_dump_path, sizeof(hard_dump_path),
                 "./output/c_code/pms_hard_init_cw%04u.txt",
                 pms_hard_cw_seq++);
        dvc_dump_hard_bits_hex128(hard_dump_path, hard_init, hm_n);
      }
#endif
      init_synd_wt = dvc_ibex_syndrome_weight(this, hard_init);
      init_synd_wt_min = std::min(init_synd_wt_min, init_synd_wt);
      init_synd_wt_max = std::max(init_synd_wt_max, init_synd_wt);
      free(hard_init);
    }
  }

  for (int i = 0; i < bm_n; i++) {
    for (int j = 0; j < cir_sz; j++) {
      cn_q_mem[i][j] = (float)llr_tbl[dec_di_blk[i * cir_sz + j]];
      if (finite_mode == 1)
        cn_q_mem[i][j] =
            (float)Sat_Quan((double)cn_q_mem[i][j], finite_q_max, finite_q_min,
                            finite_q_num, finite_f_num);
    }
  }

  // iterative decoding
  for (int itr = 0;
       (itr <= ldec_max_itr) && ((ldec_early_term_en == 0) || (cw_fail == 1));
       itr++) {
    // q sign mem index
    cir_cnt = 0;

    // layer decoding
    for (int layer = 0;
         layer < bm_m && ((ldec_early_term_en == 0) || (cw_fail == 1));
         layer++) {
      // initilize hd mem
      hd_init = (vec_sum(dec_init, bm_n) != 0);

#ifdef _ldpc_debug_dump
      printf("[ldpc debug] layer decoding @ iteration %d, layer %d ...\n", itr,
             layer);
#endif

      // init current layer c-msg
      // c-msg of previous iteration

      cn_c_sel_cur = cn_c_mem[layer];
      // c-msg to be updt
      for (int i = 0; i < cir_sz; i++) {
        cn_c_updt_cur[i].min1_val = 100000;
        cn_c_updt_cur[i].min2_val = 100000;
        cn_c_updt_cur[i].min1_pos = 0;
        cn_c_updt_cur[i].sign_tot = 1;
        min3_val[i] = 100000;
      }

      // earlier termination init
      hd_updated = 0;
      vec_clr(layer_synd, cir_sz);

      // per circulant of the layer
      for (e = mod2sparse_first_in_row(qc_bm, layer);
           !mod2sparse_at_end(e) &&
           ((ldec_early_term_en == 0) || (cw_fail == 1));
           e = mod2sparse_next_in_row(e)) {
        // read q from the previous layer of the selected column
        cn_q_sel_pre = cn_q_mem[e->col];

        // find out the previous layer of select column
        e_pre = mod2sparse_prev_in_col(e);
        if (mod2sparse_at_end(e_pre))
          e_pre = mod2sparse_last_in_col(qc_bm, e->col);

        cn_c_sel_pre = cn_c_mem[e_pre->row]; // read previous layer c msg

        // cal rnew and app
        for (int i = 0; i < cir_sz; i++) {
          // qmsg sign
          sign_tmp = (cn_q_sel_pre[i] >= 0) ? 1 : -1;

          // rnew
          if (cn_c_sel_pre[i].min1_pos == e->col)
            cn_r_new_pre[i] =
                cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
          else
            cn_r_new_pre[i] =
                cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;

          // app in cn order of previous layer
          if (h_matrix.extra_bytes_of_parity == 0) {
            cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
          } else {
            if (h_matrix.occupied[e_pre->row][e_pre->col] &&
                (e_pre->row < (h_matrix.rows - 1)))
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else if (h_matrix.occupied[e_pre->row][e_pre->col] &&
                     (e_pre->row == (h_matrix.rows - 1)) &&
                     h_matrix.mask[e_pre->col][(i + e_pre->shift) % cir_sz])
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else if (h_matrix.fade[e_pre->row][e_pre->col] &&
                     !h_matrix.mask[e_pre->col][(i + e_pre->shift) % cir_sz])
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else
              cn_app_pre[i] = cn_q_sel_pre[i];
          }

          // quantization
          if (finite_mode == 1) {
            cn_app_pre[i] =
                (float)Sat_Quan((double)cn_app_pre[i], finite_q_max,
                                finite_q_min, finite_q_num, finite_f_num);
          }
        }

        // app shift
        // when decoder initilized, q msg are in vn order
        if (dec_init[e->col] == 1) {
          shift_val1 = e->shift;
          shift_val2 = 0;
          dec_init[e->col] = 0;
        } else {
          shift_val1 = -1 * e_pre->shift + e->shift;
          shift_val2 = -1 * e_pre->shift;
        }

        for (int i = 0; i < cir_sz; i++) {
          cn_app_cur[i] = cn_app_pre[(i + shift_val1 + cir_sz) % cir_sz];
          vn_dec_hd[i] =
              cn_app_pre[(i + shift_val2 + cir_sz) % cir_sz] >= 0 ? 0 : 1;
        }

#ifdef _ldpc_debug_dump
        if (lfp) {
          fprintf(lfp, "itr%2d/layer%2d/col%2d: \n", itr, layer, e->col);
          const int dump_groups_8 = (cir_sz + 7) / 8;
          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "q pre msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_q_sel_pre[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "r new msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_r_new_pre[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "app-c msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_app_pre[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "app-s msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_app_cur[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }
        }
#endif

        // cw converge check logic per circulant
        // 1. check if hd updated
        if (hd_updated == 0)
          if (vec_cmp(dec_do_blk, vn_dec_hd, e->col * cir_sz, 0, cir_sz) == 1)
            hd_updated = 1;

#ifdef _ldpc_debug_dump
        if (lfp) {
          for (int i = 0; i < cir_sz; i++) {
            if (dec_do_blk[e->col * cir_sz + i] != vn_dec_hd[i]) {
              fprintf(lfp, "itr%2d/layer%2d/col%2d: flip bit %d (%d-->%d)\n",
                      itr, layer, e->col, i, dec_do_blk[e->col * cir_sz + i],
                      vn_dec_hd[i]);
            }
          }
        }
#endif
        vec_copy(vn_dec_hd, dec_do_blk, 0, e->col * cir_sz, cir_sz);

#ifdef _ldpc_debug_dump
        if (hdfp) {
          fprintf(hdfp, "itr%2d/layer%2d/col%2d: ", itr, layer, e->col);
          for (int i = cir_sz / 4 - 1; i >= 0; i--) {
            stmp = 0;
            for (int j = 3; j >= 0; j--)
              stmp = stmp * 2 + vn_dec_hd[i * 4 + j];
            fprintf(hdfp, "%lx", stmp);
          }
          fprintf(hdfp, "\n");
        }
#endif

        // 2. accumulate syndrome
        vec_shift(vn_dec_hd, cn_dec_hd, cir_sz, -1 * e->shift);
        if (h_matrix.extra_bytes_of_parity == 0) {
          vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);
        } else {
          if (h_matrix.occupied[e->row][e->col] &&
              (e->row == (h_matrix.rows - 1))) {
            for (int i = 0; i < cir_sz; i++) {
              if (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])
                cn_dec_hd[i] = 0;
            }
          }
          if (h_matrix.fade[e->row][e->col]) {
            for (int i = 0; i < cir_sz; i++) {
              if (h_matrix.mask[e->col][(i + e->shift) % cir_sz])
                cn_dec_hd[i] = 0;
            }
          }
          vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, cir_sz);
        }

        // calculate r_old and current q, update current layer c and q
        for (int i = 0; i < cir_sz; i++) {
          if (cn_c_sel_cur[i].min1_pos == e->col)
            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val *
                              cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
          else
            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val *
                              cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

          // q -= rold
          if (h_matrix.extra_bytes_of_parity == 0) {
            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          } else {
            if (h_matrix.occupied[e->row][e->col] &&
                (e->row == (h_matrix.rows - 1)) &&
                (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              cn_r_old_cur[i] = 0;
            }
            if (h_matrix.fade[e->row][e->col] &&
                (h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              cn_r_old_cur[i] = 0;
            }
            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          }

          // quantization
          if (finite_mode == 1) {
            cn_q_updt_cur[i] =
                (float)Sat_Quan((double)cn_q_updt_cur[i], finite_q_max,
                                finite_q_min, finite_q_num, finite_f_num);
          }

          // update c
          if (h_matrix.extra_bytes_of_parity == 0) {
            sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
            val_tmp = cn_q_updt_cur[i] * sign_tmp;
          } else {
            if (h_matrix.occupied[e->row][e->col] &&
                (e->row == (h_matrix.rows - 1)) &&
                (!h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              sign_tmp = 1;
              val_tmp = 100000;
            }
            if (h_matrix.fade[e->row][e->col] &&
                (h_matrix.mask[e->col][(i + e->shift) % cir_sz])) {
              sign_tmp = 1;
              val_tmp = 100000;
            } else {
              sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
              val_tmp = cn_q_updt_cur[i] * sign_tmp;
            }
          }

          cn_c_updt_cur[i].sign_tot *= sign_tmp;

          if (val_tmp < cn_c_updt_cur[i].min1_val) {
            min3_val[i] = cn_c_updt_cur[i].min2_val;
            cn_c_updt_cur[i].min2_val = cn_c_updt_cur[i].min1_val;
            cn_c_updt_cur[i].min1_val = val_tmp;
            cn_c_updt_cur[i].min1_pos = e->col;
          } else if (val_tmp < cn_c_updt_cur[i].min2_val) {
            min3_val[i] = cn_c_updt_cur[i].min2_val;
            cn_c_updt_cur[i].min2_val = val_tmp;
          } else if (val_tmp < min3_val[i]) {
            min3_val[i] = val_tmp;
          }

          cn_q_sign[cir_cnt][i] = sign_tmp;
        }

        // update q memory
        for (int i = 0; i < cir_sz; i++) {
          cn_q_mem[e->col][i] = cn_q_updt_cur[i];
        }

#ifdef _ldpc_debug_dump
        if (lfp) {
          const int dump_groups_8 = (cir_sz + 7) / 8;
          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "r old msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_r_old_cur[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }

          for (int gi = 0; gi < dump_groups_8; gi++) {
            fprintf(lfp, "q new msg:");
            for (int gj = 0; gj < 8; gj++) {
              const int idx = gi * 8 + gj;
              if (idx >= cir_sz)
                break;
              vtmp = int(cn_q_updt_cur[idx] * pow(2, finite_f_num));
              if (vtmp < 0) {
                vtmp = -vtmp;
                fprintf(lfp, " %3d-%2x", idx, vtmp);
              } else {
                fprintf(lfp, " %3d+%2x", idx, vtmp);
              }
            }
            fprintf(lfp, "\n");
          }
        }
#endif

        cir_cnt++;
      } // per circulant

      // update c_msg per layer
      for (int i = 0; i < cir_sz; i++) {
        pl1 = ldpc_pms_ind(cn_c_updt_cur[i].min1_val, cn_c_updt_cur[i].min2_val);
        pl2 = ldpc_pms_ind(cn_c_updt_cur[i].min2_val, cn_c_updt_cur[i].min2_val);

        cn_c_mem[layer][i].min1_val = cn_c_updt_cur[i].min1_val * alpha_pms[pl1] - beta_pms[pl1];
        cn_c_mem[layer][i].min2_val = cn_c_updt_cur[i].min1_val * alpha_pms[pl2] - beta_pms[pl2];

        if (cn_c_mem[layer][i].min1_val < 0)
          cn_c_mem[layer][i].min1_val = 0;
        if (cn_c_mem[layer][i].min2_val < 0)
          cn_c_mem[layer][i].min2_val = 0;

        if (finite_mode == 1) {
          cn_c_mem[layer][i].min1_val =
              (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, finite_c_max,
                              finite_c_min, finite_c_num, finite_f_num);
          cn_c_mem[layer][i].min2_val =
              (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, finite_c_max,
                              finite_c_min, finite_c_num, finite_f_num);
        }

        cn_c_mem[layer][i].min1_pos = cn_c_updt_cur[i].min1_pos;
        cn_c_mem[layer][i].sign_tot = cn_c_updt_cur[i].sign_tot;
#ifdef _ldpc_debug_dump
        const int cmem_dump_scale = 1 << finite_f_num;
        fprintf(cfp, "itr%d/l%d/c%d: min1 %x, min2 %x, min1 pos %d, sign_tot %d\n",
          itr, layer, i,
          int(cmem_dump_scale * cn_c_mem[layer][i].min1_val),
          int(cmem_dump_scale * cn_c_mem[layer][i].min2_val),
          cn_c_mem[layer][i].min1_pos,
          (1-cn_c_mem[layer][i].sign_tot)/2 );
#endif
      }

#ifdef _ldpc_debug_dump
      if (sfp) {
        for (int i = (cir_sz / 4 - 1); i >= 0; i--) {
          stmp = 0;
          for (int j = 3; j >= 0; j--)
            stmp = stmp * 2 + (1 - cn_c_mem[layer][i * 4 + j].sign_tot) / 2;
          fprintf(sfp, "%lx", stmp);
        }
        fprintf(sfp, "\n");
      }
#endif

      // check converage checking
      layer_synd_wt = vec_sum(layer_synd, cir_sz);
      if (hd_init == 1) {
        hd_stable_cnt = 0;
        synd_pass_cnt = 0;
      } else if ((hd_updated == 0) && (layer_synd_wt == 0)) {
        hd_stable_cnt++;
        synd_pass_cnt++;
      } else {
        hd_stable_cnt = 0;
        synd_pass_cnt = 0;
      }

#ifdef _ldpc_debug_dump
      if (lfp) {
        fprintf(lfp,
                "[dvc] layer summary: itr=%d layer=%d layer_synd_wt=%d hd_updated=%d synd_pass_cnt=%d hd_stable_cnt=%d cw_fail=%d\n",
                itr, layer, layer_synd_wt, hd_updated, synd_pass_cnt,
                hd_stable_cnt, cw_fail);
      }
#endif

      if ((synd_pass_cnt >= bm_m) && (hd_stable_cnt >= bm_m - 1)) {
        cw_fail = 0;
        cnvg_itr = itr;
        cnvg_lyr = layer;
#ifdef _ldpc_debug_dump
        if (lfp && !dbg_converged_printed) {
          fprintf(lfp, "[dvc] converged: itr=%d layer=%d (early_term=%d)\n",
                  itr, layer, ldec_early_term_en);
          fflush(lfp);
          dbg_converged_printed = 1;
        }
#endif
      }
    }
  }

  if ((cw_fail == 1) || (ldec_early_term_en == 0)) {
    cnvg_itr = ldec_max_itr - 1;
    cnvg_lyr = bm_m - 1;
  }

  fina_synd_wt = dvc_ibex_syndrome_weight(this, dec_do_blk);

#ifdef _ldpc_debug_dump
  if (lfp) {
    fprintf(lfp,
            "[dvc] cw summary: init_synd_wt=%d final_synd_wt=%d init_ok=%d final_ok=%d cw_fail=%d cnvg_itr=%d cnvg_lyr=%d\n",
            init_synd_wt, fina_synd_wt, (init_synd_wt == 0),
            (fina_synd_wt == 0), cw_fail, cnvg_itr, cnvg_lyr);
    dvc_dump_col_major_codeword(lfp, dec_do_blk, bm_n, cir_sz, "dec_do");
    fflush(lfp);
  }
#endif

  // free all
  free(dec_init);
  for (int i = 0; i < bm_m; i++)
    free(cn_c_mem[i]);
  free(cn_c_mem);
  free(cn_c_updt_cur);
  for (int i = 0; i < bm_n; i++)
    free(cn_q_mem[i]);
  free(cn_q_mem);
  free(cn_r_new_pre);
  free(cn_app_pre);
  free(cn_app_cur);
  free(cn_q_sel_cur);
  free(cn_r_old_cur);
  free(cn_q_updt_cur);
  for (int i = 0; i < total_cir; i++)
    free(cn_q_sign[i]);
  free(cn_q_sign);
  free(layer_synd);
  free(cn_dec_hd);
  free(vn_dec_hd);
  free(min3_val);

#ifdef _ldpc_debug_dump
  if (cfp)
    fclose(cfp);
  if (sfp)
    fclose(sfp);
  if (hdfp)
    fclose(hdfp);
  if (lfp)
    fclose(lfp);
#endif
} // ldpc_dec_layer2

int ldpc_packet::ldpc_pms_ind(float min1, float min2) {
  int ind;
  if (min1 <= point1)
  {
    if (min2 <= point1)
      ind = 0;
    else if (min2 <= point2)
      ind = 1;
    else
      ind = 2;
  }
  else if (min1 <= point2)
  {
    if (min2 <= point2)
      ind = 3;
    else
      ind = 4;
  }
  else
    ind = 5;

  return ind;
}

void ldpc_packet::ldpc_dec_skip() {
  for (int i = 0; i < hm_n; i++)
    dec_do_blk[i] = dec_di_blk[i];

  cw_fail = ldpc_synd(dec_do_blk);
}

int ldpc_packet::ldpc_synd(char *cw) {
  int synd_fail = 0;
  char *synd;

  synd = (char *)calloc(hm_m, sizeof(*synd));

  mod2sparse_mulvec(qc_hm, cw, synd);

  for (int i = 0; i < hm_m; i++) {
    if (synd[i] != 0) {
      synd_fail = 1;
      break;
    }
  }

  free(synd);

  return synd_fail;
}

// 2-bit 专用入口：在不修改原逻辑的情况下，通过调整参数禁用激进加权
// ldpc_packet::ldpc_dec_bf_ibex_2bit(s_ldpc_decoder_input ldpc_decoder_input,
//                                    s_ldpc_decoder_parameters
//                                    ldpc_decoder_parameters, s_h_matrix
//                                    h_matrix) {
//   if (vn_bits <= 2) {
//     // 令 be_aggressive 判定为假，避免 2bit 下权重放大导致抖动
//     ldpc_decoder_parameters.likelihood_thr = 0;
//   }
//   ldpc_dec_bf_ibex(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);
// }

void ldpc_packet::ldpc_dec_bf_ibex(
    s_ldpc_decoder_input ldpc_decoder_input,
    s_ldpc_decoder_parameters ldpc_decoder_parameters, s_h_matrix h_matrix) {
  int VERBOSITY = 0;
  int MAX_ERROR_COUNT = 4095;
  int i;
  int j;
  int k;
  int m;
  int iteration = 0;
  int clock_cycles = 0;
  int syndrome_weight;
  int weight;
  int decode_mode;
  int decode_mode2;
  bool post_process = 0;
  bool post_process_2 = 0;
  bool do_post_flipped = 0;
  bool do_post_unflipped = 0;
  bool give_up = 0;
  bool finished = 0;
  bool post_processing = 0;
  bool qc_parity_en;
  bool disable_update;
  bool hamming_weight_le_circ_thr;
  bool hamming_weight_lt_circ_thr;
  bool hamming_weight_lt_post_thr;
  bool post_trigger;
  bool post_trigger2;
  bool flipped_prev;
  bool flipped;
  bool be_aggressive;
  bool prng_post_process;
  bool prng_post_process2;
  s_hard_codeword hard_codeword;
  s_hard_codeword soft_codeword;
  s_variable_nodes vn;
  s_check_nodes cn;
  // s_check_nodes cn_delay1;
  // s_check_nodes cn_delay2;
  s_check_nodes cn_shifted;
  s_likelihood_levels likelihood_levels;
  s_256_bits prng_256;
  s_512_bits prng_512;
  int prng_init[32];
  int syndrome_weight_delayed;

  int syndrome_weight_r[5];
  bool do_not_use_this_bit;
  bool look;
  FILE *sw_delta_fp = dvc_open_bf_ibex_sw_delta_dump_file();
  FILE *row_sw_fp = dvc_open_bf_ibex_row_sw_dump_file();
  if (sw_delta_fp)
    fprintf(sw_delta_fp, "%8s %12s %6s %16s\n", "iter", "phase", "col",
            "syndrome_weight");
  dvc_log_bf_ibex_row_sw_header(row_sw_fp, h_matrix.rows);

#ifdef _LDPC_DBG_DUMP
  static int dbg_made_dir = 0;
  if (!dbg_made_dir) {
    system("mkdir -p output >/dev/null 2>&1");
    dbg_made_dir = 1;
  }
  FILE *dbg_sw_fp = fopen("./output/ldpc_dbg_sw.txt", "w");
  if (dbg_sw_fp)
    fprintf(dbg_sw_fp, "%6s %6s\n", "iter", "sw");
  FILE *dbg_tr_fp = NULL;
  FILE *dbg_like_fp = dvc_open_bf_ibex_like_dump_file();
#endif

  prng_init[31] = 0x083d;
  prng_init[30] = 0x3214;
  prng_init[29] = 0xa8a1;
  prng_init[28] = 0x5327;
  prng_init[27] = 0x71bc;
  prng_init[26] = 0x3edb;
  prng_init[25] = 0xba50;
  prng_init[24] = 0xc946;
  prng_init[23] = 0x4c9a;
  prng_init[22] = 0x0b73;
  prng_init[21] = 0xef18;
  prng_init[20] = 0x31ee;
  prng_init[19] = 0xf2b2;
  prng_init[18] = 0xd98f;
  prng_init[17] = 0x89f9;
  prng_init[16] = 0x9375;
  prng_init[15] = 0x043d;
  prng_init[14] = 0x6214;
  prng_init[13] = 0xa8a1;
  prng_init[12] = 0x5d27;
  prng_init[11] = 0x71cc;
  prng_init[10] = 0x2edb;
  prng_init[9] = 0xb550;
  prng_init[8] = 0xc646;
  prng_init[7] = 0x8c9a;
  prng_init[6] = 0x1b73;
  prng_init[5] = 0xef08;
  prng_init[4] = 0x30ee;
  prng_init[3] = 0xf7b2;
  prng_init[2] = 0xda8f;
  prng_init[1] = 0x49f9;
  prng_init[0] = 0x9365;

  for (j = 0; j < h_matrix.cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      hard_codeword.c[j].b[k] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard;
      vn.c[j].b[k].bit_hard = hard_codeword.c[j].b[k];
      soft_codeword.c[j].b[k] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable;
      vn.c[j].b[k].flipped = 0;
    }
  }

  cn = f_check_nodes(h_matrix, hard_codeword);
  // cn_delay1 = cn;
  // cn_delay2 = cn;
  if (VERBOSITY > 0) {
    printf("[LDPC DEBUG] Starting BF decoding with max %d iterations.\n",
           fdec_max_itr);
    f_print_hard_codeword(hard_codeword, h_matrix.cols, h_matrix.bits);
    printf("### DECODER C++: H MATRIX:\n");
    f_print_h_matrix(h_matrix);
    printf(
        "### DECODER C++: AFTER SYNDROME CHECK, CODEWORD HAS CHECK NODES:\n");
    f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
  }

  for (i = 0; i < h_matrix.rows; i++)
    for (k = 0; k < h_matrix.bits; k++)
      cn_shifted.r[i].b[k] = 0;
  syndrome_weight = f_check_node_weight(h_matrix, cn);
  clock_cycles = (2 * h_matrix.cols) + 1;
  if (syndrome_weight == 0)
    finished = 1;
  ldpc_decoder_output.syndrome_weight_before = syndrome_weight;
  syndrome_weight_delayed = syndrome_weight;
  for (i = 0; i < 5; i++)
    syndrome_weight_r[i] = syndrome_weight;

  int prev_sw =
      (iteration == 0) ? syndrome_weight_delayed : syndrome_weight_r[3];
  bool pushing = (syndrome_weight_delayed >= prev_sw);

  ldpc_decoder_output.early_termination = 0;

  if (ldpc_decoder_parameters.early_terminate_dis == 0) {
    int early_term_thr;
    if (ldpc_decoder_input.nand_strobes == 0) {
      if (h_matrix.rows <= 6)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][0];
      else if (h_matrix.rows == 7)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][1];
      else if (h_matrix.rows == 8)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][2];
      else if (h_matrix.rows == 9)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][3];
      else if (h_matrix.rows == 10)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][4];
      else if (h_matrix.rows == 11)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][5];
      else
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[0][6];
    } else {
      if (h_matrix.rows <= 6)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][0];
      else if (h_matrix.rows == 7)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][1];
      else if (h_matrix.rows == 8)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][2];
      else if (h_matrix.rows == 9)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][3];
      else if (h_matrix.rows == 10)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][4];
      else if (h_matrix.rows == 11)
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][5];
      else
        early_term_thr = ldpc_decoder_parameters.early_terminate_thr[1][6];
    }

    if (syndrome_weight >= early_term_thr) {
      finished = 1;
      ldpc_decoder_output.early_termination = 1;
      if ((VERBOSITY > 0) && (syndrome_weight >= early_term_thr))
        printf("### LDPC DECODER C++ EARLY TERMINATION: %5d\n",
               syndrome_weight);
    }
  }
  if (ldpc_decoder_input.syndrome_cal_only) {
    finished = 1;
  }

  likelihood_levels = f_likelihood_levels(ldpc_decoder_input.nand_strobes,
                                          ldpc_decoder_parameters,
                                          syndrome_weight, h_matrix.rows);
  if (!finished)
    iteration = 1;
  const int post_start_iteration =
      (ldpc_decoder_input.post_iteration < 1) ? 1
                                              : ldpc_decoder_input.post_iteration;

  // soft_data -> likelihood_level
  bool soft_data[2];
  for (j = 0; j < h_matrix.cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      soft_data[0] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable;
      soft_data[1] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2;
      if ((soft_data[1] == 0) && (soft_data[0] == 0))
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters
                                        .likelihood_map[0]]; // 3 verilog: 3
      if ((soft_data[1] == 0) && (soft_data[0] == 1))
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters
                                        .likelihood_map[1]]; // 2 verilog: 2
      if ((soft_data[1] == 1) && (soft_data[0] == 0))
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters
                                        .likelihood_map[2]]; // 1 verilog: 3
      if ((soft_data[1] == 1) && (soft_data[0] == 1))
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters
                                        .likelihood_map[3]]; // 0 verilog: 2

      if ((VERBOSITY > 0) && (ldpc_decoder_input.soft_bits == 1) &&
          (ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable ==
           ldpc_decoder_parameters.questionable_sense))
        printf("### INITIAL %x %x %x\n", j, k, vn.c[j].b[k].likelihood);
    }
  }

  if (VERBOSITY > 0) {
    printf("### C++ INPUT DATA %2d %3d  %x %x %x  LIKELIHOOD: %d\n", 71, 405,
           ldpc_decoder_input.corrupted_codeword.c[71].b[405].bit_hard,
           ldpc_decoder_input.corrupted_codeword.c[71].b[405].bit_questionable,
           ldpc_decoder_input.corrupted_codeword.c[71].b[405].bit_questionable2,
           vn.c[71].b[405].likelihood);
  }

  // 4x5 statistic table init
  for (j = 0; j <= 3; j++) {
    for (k = 0; k <= 4; k++) {
      hard_codeword.errors_at_level_and_weight[j][k] = 0;
      hard_codeword.correct_at_level_and_weight[j][k] = 0;
    }
  }

  while ((iteration < ldpc_decoder_input.iteration_limit) && (finished == 0) &&
         (give_up == 0)) {
    // Treat syndrome calculation as the standalone iteration before this pass:
    // col0/col1 read this snapshot, and col2 first sees col0's live update.
    // Temporarily disable cn_delay for weight calculation debug.
    // cn_delay1 = cn;
    // cn_delay2 = cn;
    dvc_log_bf_ibex_sw(sw_delta_fp, iteration, "iter_pre", -1, syndrome_weight);
    dvc_log_bf_ibex_row_sw(row_sw_fp, iteration, "iter_pre", -1, h_matrix, cn);
#ifdef _LDPC_DBG_DUMP
    {
      char dbg_tr_path[128];
      snprintf(dbg_tr_path, sizeof(dbg_tr_path),
               "./output/ldpc_dbg_trace_iter%d.txt", iteration);
      dbg_tr_fp = fopen(dbg_tr_path, "w");
      if (dbg_tr_fp)
        fprintf(dbg_tr_fp, "%6s %6s %4s %3s %1s\n", "iter", "idx", "like", "wt",
                "f");
    }
#endif
    for (j = 0; j < h_matrix.cols; j++) {
      clock_cycles++;
      if ((iteration == post_start_iteration) && (j == 0)) {
        for (i = 0; i < 256; i++)
          prng_256.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
        for (i = 0; i < 512; i++)
          prng_512.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
        for (i = 0; i < 512; i++)
          prng_512.b[i] = (0x1fe0 >> (i % 16)) & 1; // to match verilog
      } else if (iteration >= post_start_iteration) {
        prng_256 = f_256_bit_lfsr(prng_256);
        prng_512 = f_512_bit_lfsr(prng_512);
      }

      syndrome_weight_r[4] = syndrome_weight_r[3];
      syndrome_weight_r[3] = syndrome_weight_r[2];
      syndrome_weight_r[2] = syndrome_weight_r[1];
      syndrome_weight_r[1] = syndrome_weight_r[0];
      syndrome_weight_r[0] = syndrome_weight;
      syndrome_weight_delayed = syndrome_weight_r[3];

      if (VERBOSITY > 0)
        printf("### C++ ITERATION %4d, COLUMN %2d, SYNDROME WEIGHT: %4d "
               "DELAYED WEIGHT: %4d ###\n",
               iteration, j, syndrome_weight, syndrome_weight_delayed);

      // distribute different thr for post processing stage
      if (iteration >= ldpc_decoder_input.post_iteration) {
        hamming_weight_lt_post_thr =
            (syndrome_weight_delayed <
             ldpc_decoder_parameters.syndrome_weight_thr_post);
        hamming_weight_lt_circ_thr =
            (syndrome_weight_delayed <
             ldpc_decoder_parameters.syndrome_weight_thr_qc);
        post_trigger =
            hamming_weight_lt_circ_thr &&
            ((iteration % 16) < ldpc_decoder_parameters.post_ratio) &&
            ldpc_decoder_parameters.post_process_en;
        post_trigger2 =
            hamming_weight_lt_post_thr &&
            ((iteration % 16) >= ldpc_decoder_parameters.post_ratio) &&
            ldpc_decoder_parameters.post_process_en;

        if (VERBOSITY > 0) {
          printf("### ITERATION %4d POST PROCESSING   BEFORE UPDATE OF COLUMN: "
                 "%2dWEIGHT: %4d DELAYED WEIGHT: %4d "
                 "POST PROCESS: %d THR_POST: %3d THR_QC: %3d POST_TRIGGER: %d "
                 "POST_TRIGGER2: %d ###\n",
                 iteration, j, syndrome_weight, syndrome_weight_delayed,
                 post_process, ldpc_decoder_parameters.syndrome_weight_thr_post,
                 ldpc_decoder_parameters.syndrome_weight_thr_qc, post_trigger,
                 post_trigger2);
        }
      } else {
        post_process = 0;
        post_trigger = 0;
        post_trigger2 = 0;
      }

      // be_aggressive is different with rtl and not used
      be_aggressive =
          (ldpc_decoder_input.soft_bits > 0) &&
          (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr) &&
          !post_trigger &&
          !post_trigger2; // different with verilog, add two post trigger judge
      // be_aggressive = (ldpc_decoder_input.soft_bits > 0) &&
      // (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr) &&
      // (iteration < (ldpc_decoder_input.post_iteration << 1)); be_aggressive =
      // (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min <
      // ldpc_decoder_parameters.likelihood_thr) && (iteration <
      // (ldpc_decoder_input.post_iteration + 100)); be_aggressive =
      // (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min <
      // ldpc_decoder_parameters.likelihood_thr);
      // const s_check_nodes &cn_for_weight = cn_delay2;
      const s_check_nodes &cn_for_weight = cn;
      for (k = 0; k < h_matrix.bits; k++) {
        look = 0;
        do_not_use_this_bit = 0;
        do_not_use_this_bit |= ((h_matrix.extra_bits_of_parity > 0) &&
                                (j == (h_matrix.cols - h_matrix.rows)) &&
                                (k >= h_matrix.extra_bits_of_parity));
        do_not_use_this_bit |= ((h_matrix.extra_bits_of_userdata > 0) &&
                                (j == (h_matrix.cols - h_matrix.rows - 1)) &&
                                (k >= h_matrix.extra_bits_of_userdata));
        if ((VERBOSITY > 0) && do_not_use_this_bit)
          printf("### DO NOT USE THIS BIT %2d %3d\n", j, k);
        if (do_not_use_this_bit) {
#ifdef _LDPC_DBG_DUMP
          if (dbg_tr_fp) {
            int idx = (j * h_matrix.bits) + k;
            fprintf(dbg_tr_fp, "%6d %6d %4d %3d %1d\n", iteration, idx,
                    (int)vn.c[j].b[k].likelihood, -1,
                    (int)vn.c[j].b[k].flipped);
          }
#endif
        } else {
          weight = 0;
          for (i = 0; i < h_matrix.rows; i++) {
            m = (k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits;
            if (h_matrix.extra_bytes_of_parity == 0) {
              if (h_matrix.occupied[i][j] &&
                  (cn_for_weight.r[i].b[m] == 1))
                weight++;
            } else {
              if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)) &&
                  (cn_for_weight.r[i].b[m] == 1))
                weight++;
              if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
                  (cn_for_weight.r[i].b[m] == 1) && h_matrix.mask[j][k])
                weight++;
              if (h_matrix.fade[i][j] && (cn_for_weight.r[i].b[m] == 1) &&
                  !h_matrix.mask[j][k])
                weight++;
            }
          }
          if (iteration == 1) {
            if (ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_is_error)
              hard_codeword.errors_at_level_and_weight
                  [ldpc_decoder_input.corrupted_codeword.c[j].b[k].level]
                  [weight]++;
            else
              hard_codeword.correct_at_level_and_weight
                  [ldpc_decoder_input.corrupted_codeword.c[j].b[k].level]
                  [weight]++;
          }
          flipped_prev = vn.c[j].b[k].flipped;
          int likelihood_prev = vn.c[j].b[k].likelihood;
          look = (VERBOSITY > 0) && ((j == 33) && (k == 163));

          if (iteration >= ldpc_decoder_input.post_iteration) {
            // normal disurbance
            prng_post_process =
                post_trigger &&
                (syndrome_weight_delayed <
                 ldpc_decoder_parameters.syndrome_weight_thr_qc) &&
                ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
            // radical disturbance
            prng_post_process2 =
                post_trigger2 &&
                (syndrome_weight_delayed <
                 ldpc_decoder_parameters.syndrome_weight_thr_post) &&
                ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
          } else {
            prng_post_process = 0;
            prng_post_process2 = 0;
          }
          if ((VERBOSITY > 0) && look) {
            printf("### C++ LOOK ITERATION %4d, SW: %4d, BIT: %2d %3d  "
                   "ORIGINAL: %x CORRUPTED: %x    WEIGHT: %1d MIN: "
                   "%3d LIKELIHOOD: %3d POST_TRIGGER: %d %d POST: %d %d "
                   "PRNG[%3d]: %d\n",
                   iteration, syndrome_weight_delayed, j, k, 0,
                   ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                   weight, likelihood_levels.min, vn.c[j].b[k].likelihood,
                   post_trigger, post_trigger2, prng_post_process,
                   prng_post_process2, k, prng_512.b[k]);
            f_print_s_512_bits(prng_512);
          }

          // Adjust weight if aggressive mode
          bool aggr =
              (ldpc_decoder_input.soft_bits > 0) &&
              (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr &&
               !flipped_prev);
          int w = weight;
          if (aggr && (weight == 0))
            w = weight + 0;
          if (aggr && (weight == 1))
            w = weight + 0; // =1, verilog delta: 0
          if (aggr && (weight == 2))
            w = weight + 1; // =3, verilog delta: 2
          if (aggr && (weight == 3))
            w = weight + 2; // =5, verilog delta: 4
          if (aggr && (weight == 4))
            w = weight + 3; // =7, verilog delta: 7
          if ((VERBOSITY > 0) && look)
            printf("C++ LOOK   ITERATION: %4d BEFORE LIKELIHOOD UPDATE SW: %4d "
                   "CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x "
                   "CORRUPTED: %x FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: "
                   "%3d\n",
                   iteration, syndrome_weight, j, k, 0,
                   ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                   vn.c[j].b[k].flipped, weight, likelihood_levels.min,
                   vn.c[j].b[k].likelihood);
          vn.c[j].b[k].likelihood = f_update_vn_post(
              vn.c[j].b[k].likelihood, w, likelihood_levels.min,
              likelihood_levels.max, prng_post_process, prng_post_process2,
              aggr, likelihood_levels.flip_thr, pushing);

          if ((VERBOSITY > 0) && look)
            printf("C++ LOOK   ITERATION: %4d AFTER LIKELIHOOD UPDATE SW: %4d "
                   "CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x "
                   "CORRUPTED: %x FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: "
                   "%3d\n",
                   iteration, syndrome_weight, j, k, 0,
                   ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                   vn.c[j].b[k].flipped, weight, likelihood_levels.min,
                   vn.c[j].b[k].likelihood);
          vn.c[j].b[k].flipped =
              (vn.c[j].b[k].likelihood >= likelihood_levels.flip_thr);

#ifdef _LDPC_DBG_DUMP
          if (dbg_tr_fp) {
            int idx = (j * h_matrix.bits) + k;
            fprintf(dbg_tr_fp, "%6d %6d %4d %3d %1d\n", iteration, idx,
                    (int)vn.c[j].b[k].likelihood, weight,
                    (int)vn.c[j].b[k].flipped);
          }
#endif

          if ((VERBOSITY > 0) && look)
            printf(
                "C++ ITERATION: %4d AFTER FLIP CHECK  SW: %4d CURRENT BIT "
                "FLIP: %2d %3d  ORIGINAL: %x CORRUPTED: %x "
                "FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: %3d THR: %3d\n",
                iteration, syndrome_weight, j, k, 0,
                ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                vn.c[j].b[k].flipped, weight, likelihood_levels.min,
                vn.c[j].b[k].likelihood, likelihood_levels.flip_thr);

          // Increment update syndrome if flipped state changed
          if (flipped_prev != vn.c[j].b[k].flipped) {
            for (i = 0; i < h_matrix.rows; i++) {
              m = (k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits;
              if (h_matrix.extra_bytes_of_parity == 0) {
                if (h_matrix.occupied[i][j])
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
              } else {
                if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)))
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
                if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
                    h_matrix.mask[j][k])
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
                if (h_matrix.fade[i][j] && !h_matrix.mask[j][k])
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
              }
            }
            if ((VERBOSITY > 0) && vn.c[j].b[k].flipped)
              printf("### C++ ITERATION %4d   FLIPPED BIT: %2d %3d ###\n",
                     iteration, j, k);
            if ((VERBOSITY > 0) && !vn.c[j].b[k].flipped)
              printf("### C++ ITERATION %4d   UNFLIPPED BIT: %2d %3d ###\n",
                     iteration, j, k);

            if (VERBOSITY > 0) {
              printf(
                  "### DECODER C++: AFTER TOGGLE, CODEWORD HAS CHECK NODES:\n");
              f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
            }
          }
        }
      }

#ifdef _LDPC_DBG_DUMP
      dvc_dump_bf_ibex_col_likelihood(dbg_like_fp, vn, iteration, j,
                                      h_matrix.bits);

      char dbg_cn_col_path[128];
      snprintf(dbg_cn_col_path, sizeof(dbg_cn_col_path),
               "./output/ldpc_dbg_cn_synd_iter%d_col%02d.txt", iteration, j);
      dvc_dump_cn_synd_file(dbg_cn_col_path, cn, h_matrix);

      char dbg_cn_rshift_col_path[128];
      snprintf(dbg_cn_rshift_col_path, sizeof(dbg_cn_rshift_col_path),
               "./output/ldpc_dbg_cn_synd_rshift_iter%d_col%02d.txt",
               iteration, j);
      dvc_dump_cn_synd_shifted_col_file(dbg_cn_rshift_col_path, cn, h_matrix,
                                        j, true);

      char dbg_cn_lshift_col_path[128];
      snprintf(dbg_cn_lshift_col_path, sizeof(dbg_cn_lshift_col_path),
               "./output/ldpc_dbg_cn_synd_lshift_iter%d_col%02d.txt",
               iteration, j);
      dvc_dump_cn_synd_shifted_col_file(dbg_cn_lshift_col_path, cn, h_matrix,
                                        j, false);
#endif

      for (i = 0; i < h_matrix.rows; i++) {
        if (h_matrix.occupied[i][j] || h_matrix.fade[i][j])
          for (k = 0; k < h_matrix.bits; k++)
            cn_shifted.r[i].b[k] =
                cn.r[i].b[(k + h_matrix.bits - h_matrix.element[i][j]) %
                          h_matrix.bits];
      }
      if (VERBOSITY > 0) {
        printf("### DECODER C++: DECODING ITERATION %4d AT COLUMN %2d, "
               "CODEWORD HAS CHECK NODES:\n",
               iteration, j);
        f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
        printf("### DECODER C++: DECODING ITERATION %4d AT COLUMN %2d, "
               "CODEWORD HAS CHECK NODES SHIFTED:\n",
               iteration, j);
        f_print_check_nodes_shifted(cn, h_matrix, j);
      }

      syndrome_weight = f_check_node_weight(h_matrix, cn);
      if (sw_delta_fp) {
        dvc_log_bf_ibex_sw(sw_delta_fp, iteration, "col_post", j,
                           syndrome_weight);
      }
      dvc_log_bf_ibex_row_sw(row_sw_fp, iteration, "col_post", j, h_matrix, cn);
      // cn_delay2 = cn_delay1;
      // cn_delay1 = cn;
      if (syndrome_weight == 0) {
        ldpc_decoder_output.col_cnt = j;
        finished = 1;
        break;
      }
    }

    syndrome_weight = f_check_node_weight(h_matrix, cn);
    finished = (syndrome_weight == 0);
    syndrome_weight_r[4] = syndrome_weight_r[3];
    syndrome_weight_r[3] = syndrome_weight_r[2];
    syndrome_weight_r[2] = syndrome_weight_r[1];
    syndrome_weight_r[1] = syndrome_weight_r[0];
    syndrome_weight_r[0] = syndrome_weight;

#ifdef _LDPC_DBG_DUMP
    if (dbg_sw_fp) {
      fprintf(dbg_sw_fp, "%6d %6d\n", iteration, syndrome_weight);
      fflush(dbg_sw_fp);
    }

    char dbg_cn_path[128];
    snprintf(dbg_cn_path, sizeof(dbg_cn_path),
             "./output/ldpc_dbg_cn_synd_iter%d.txt", iteration);
    dvc_dump_cn_synd_file(dbg_cn_path, cn, h_matrix);

    char dbg_cn_rshift_path[128];
    snprintf(dbg_cn_rshift_path, sizeof(dbg_cn_rshift_path),
             "./output/ldpc_dbg_cn_synd_rshift_iter%d.txt", iteration);
    dvc_dump_cn_synd_shifted_file(dbg_cn_rshift_path, cn, h_matrix, true);

    char dbg_cn_lshift_path[128];
    snprintf(dbg_cn_lshift_path, sizeof(dbg_cn_lshift_path),
             "./output/ldpc_dbg_cn_synd_lshift_iter%d.txt", iteration);
    dvc_dump_cn_synd_shifted_file(dbg_cn_lshift_path, cn, h_matrix, false);

    if (dbg_tr_fp) {
      fclose(dbg_tr_fp);
      dbg_tr_fp = NULL;
    }
#endif

    clock_cycles++;
    iteration++;
  }

  if (VERBOSITY > 0) {
    printf("### AFTER %4d ITERATIONS, CODEWORD HAS CHECK NODES:\n", iteration);
    f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
  }

  if (sw_delta_fp)
    fclose(sw_delta_fp);
  if (row_sw_fp)
    fclose(row_sw_fp);

  ldpc_decoder_output.iterations = iteration;
  ldpc_decoder_output.clock_cycles = clock_cycles;
  ldpc_decoder_output.syndrome_weight_after = f_check_node_weight(h_matrix, cn);
  ldpc_decoder_output.failure =
      (ldpc_decoder_output.syndrome_weight_after != 0);
  ldpc_decoder_output.errors_in_userdata = 0;
  ldpc_decoder_output.errors_in_codeword = 0;
  for (j = 0; j < h_matrix.cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      ldpc_decoder_output.corrected_codeword.c[j].b[k] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard ^
          vn.c[j].b[k].flipped;
      dec_do_blk[j * h_matrix.bits + k] =
          ldpc_decoder_output.corrected_codeword.c[j].b[k];
      if (vn.c[j].b[k].flipped)
        ldpc_decoder_output.errors_in_codeword++;
      if (vn.c[j].b[k].flipped && (j < (h_matrix.cols - h_matrix.rows)))
        ldpc_decoder_output.errors_in_userdata++;
      if ((VERBOSITY > 0) && (vn.c[j].b[k].flipped))
        printf("### C FLIPPED BIT %2d %3d ### %4d ERRORS IN CODEWORD, %4d "
               "ERRORS IN USERDATA\n",
               j, k, ldpc_decoder_output.errors_in_codeword,
               ldpc_decoder_output.errors_in_userdata);
    }
  }

  if (ldpc_decoder_output.early_termination) {
    if (VERBOSITY > 0)
      printf("### C EARLY TERMINATION ###\n");
    ldpc_decoder_output.iterations = 0;
    ldpc_decoder_output.syndrome_weight_after = syndrome_weight;
    ldpc_decoder_output.failure = 1;
    ldpc_decoder_output.errors_in_userdata = 0;
    ldpc_decoder_output.errors_in_codeword = 0;
  } else if (ldpc_decoder_input.syndrome_cal_only) {
    if (VERBOSITY > 0)
      printf("### C SYNDROME CALCULATION ONLY ###\n");
    ldpc_decoder_output.iterations = 0;
    ldpc_decoder_output.syndrome_weight_after = syndrome_weight;
    ldpc_decoder_output.failure = 0;
    ldpc_decoder_output.errors_in_userdata = 0;
    ldpc_decoder_output.errors_in_codeword = 0;
  } else if (ldpc_decoder_output.syndrome_weight_before == 0) {
    if (VERBOSITY > 0)
      printf("### C NO ERRORS ###\n");
    ldpc_decoder_output.iterations = 0;
    ldpc_decoder_output.syndrome_weight_after = 0;
    ldpc_decoder_output.failure = 0;
    ldpc_decoder_output.errors_in_userdata = 0;
    ldpc_decoder_output.errors_in_codeword = 0;
  }
  if (ldpc_decoder_output.errors_in_codeword > MAX_ERROR_COUNT)
    ldpc_decoder_output.errors_in_codeword = MAX_ERROR_COUNT;
  if (ldpc_decoder_output.errors_in_userdata > MAX_ERROR_COUNT)
    ldpc_decoder_output.errors_in_userdata = MAX_ERROR_COUNT;
  if (VERBOSITY > 0) {
    printf("### C++ OUTPUT ###\n");
    f_print_hard_codeword(ldpc_decoder_output.corrected_codeword, h_matrix.cols,
                          h_matrix.bits);
  }

#ifdef _LDPC_DBG_DUMP
  if (dbg_sw_fp)
    fclose(dbg_sw_fp);
  if (dbg_tr_fp)
    fclose(dbg_tr_fp);
  if (dbg_like_fp)
    fclose(dbg_like_fp);
#endif

  init_synd_wt = ldpc_decoder_output.syndrome_weight_before;
  fina_synd_wt = ldpc_decoder_output.syndrome_weight_after;
  // Align DV status semantics with ldpc_dec_layer2():
  // - cnvg_itr is a 0-based "last iteration index" (not a count).
  // - on failure, cnvg_col is forced to 0 per DV expectation.
  /*
  if (ldpc_decoder_output.failure)
    cnvg_itr = ldpc_decoder_input.iteration_limit - 1;
  else
    cnvg_itr = ldpc_decoder_output.iterations - 1;
  */
  if (iteration == 0)
    cnvg_itr = 1;
  else
    cnvg_itr = ldpc_decoder_output.iterations;

  // `cnvg_lyr` is reported to DV as `cnvg_col` for BF-style decoders.
  // ldpc_dec_bf_ibex tracks the convergence column in `ldpc_decoder_output.col_cnt`.
  // - On success, propagate the 0-based convergence column index.
  // - On failure (or if no column index is available), force to 0 to match DV expectation.
  if (ldpc_decoder_output.failure || (ldpc_decoder_output.col_cnt < 0))
    cnvg_lyr = 0;
  else
    cnvg_lyr = ldpc_decoder_output.col_cnt;

  cw_fail = ldpc_decoder_output.failure;
}
