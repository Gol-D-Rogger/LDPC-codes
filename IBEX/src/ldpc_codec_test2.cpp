#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstdint>
#include <vector>

#include "finite_lib.h"
#include "ldpc_codec.h"
#include "mod2convert.h"
#include "mod2dense.h"
#include "mod2sparse.h"
#include "vec_op.h"

void ldpc_packet::ldpc_rd_phck(char *pchk_file) {
  mod2entry *e;
  FILE *fp;
  int col_shift;
  qc_bm = mod2sparse_allocate(bm_m, bm_n);
  qc_hm = mod2sparse_allocate(hm_m, hm_n);

  printf("[LDPC] Porting H matrix from %s ...\n", pchk_file);

  fp = fopen(pchk_file, "r");

  if (fp == NULL) {
    printf("[LDPC] Error opening file %s \n", pchk_file);
    exit(0);
  }

  for (int i = 0; i < bm_m; i++) {
    for (int j = 0; j < bm_n; j++) {
      fscanf(fp, "%d", &col_shift);

     // check submatrix T (identity matrix)
      if ((i < tm_sz) && (j > (bm_n - tm_sz - 1))) {
        if (((i != (j + tm_sz - bm_n)) && (col_shift >= 0)) || ((i == (j + tm_sz - bm_n)) && (col_shift != 0))) {
          printf("[LDPC] Error: Submatrix T is not identity matrix!\n");
          exit(0);
        }
      }

      // insert mod2entry to base matrix
      if (col_shift >= 0) {
        e = mod2sparse_insert(qc_bm, i, j);
        e->shift = col_shift;

        for (int k = 0; k < cir_sz; k++)
          mod2sparse_insert(qc_hm, i * cir_sz + k, j * cir_sz + (k + col_shift) % cir_sz);
      }
    }
  }

  fclose(fp);
  printf("[LDPC] H matrix porting ready!\n");
}

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

  // generate inverse F matrix (F=E*B+D)
  printf("[LDPC] Generating inverse F matrix from H matrix ...\n");

  qc_exb = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);
  qc_f = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);

  mod2sparse_multiply(qc_e, qc_b, qc_exb);
  mod2sparse_add(qc_exb, qc_d, qc_f);

  qc_f_d = mod2dense_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);
  qc_fi_d = mod2dense_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);

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
void ldpc_packet::f_print_hard_codeword(s_hard_codeword data, int cols, int bits) {
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

void ldpc_packet::f_print_check_nodes_shifted(s_check_nodes data, s_h_matrix h_matrix, int column) {
  int i;
  int k;
  int b;
  s_check_nodes data_shifted;

  for (i = 0; i < h_matrix.rows; i++) {
    if (h_matrix.occupied[i][column] || h_matrix.fade[i][column])
      for (k = 0; k < h_matrix.bits; k++)
        data_shifted.r[i].b[k] = data.r[i].b[(h_matrix.bits + k - h_matrix.element[i][column]) % h_matrix.bits];
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
s_check_nodes ldpc_packet::f_check_nodes(s_h_matrix h_matrix, s_hard_codeword vn) {
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
          }
          else {
            if (h_matrix.occupied[i][j] && (i < h_matrix.rows - 1))
              cn.r[i].b[m] = !cn.r[i].b[m];
            if (h_matrix.occupied[i][j] && (i == h_matrix.rows - 1) && h_matrix.mask[j][k])
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
s_likelihood_levels ldpc_packet::f_likelihood_levels(int strobes, s_ldpc_decoder_parameters ldpc_decoder_parameters,
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
    likelihood_levels.flip_thr = 2;      // 3
    likelihood_levels.weak = 0; // 2
    likelihood_levels.strong = 0;        // 1
  } else if (VN_BITS == 3) {
    likelihood_levels.flip_thr = likelihood_levels.max - 3;
    likelihood_levels.weak = likelihood_levels.flip_thr - 4;
    likelihood_levels.strong = likelihood_levels.weak;
  } else {
    likelihood_levels.flip_thr = std::max(likelihood_levels.min + 1, likelihood_levels.max - 7);
    likelihood_levels.weak = std::min(likelihood_levels.min, likelihood_levels.flip_thr - 4);
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
    else if (rows == 13)
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
      delta_total = (VN_BITS == 8) ? ((delta_sum >> 1) + (delta_sum >> 2)) : ((delta_sum >> 2) + (delta_sum >> 3));
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
          likelihood_levels.level[3] - ((weak_minus_strong * ldpc_decoder_parameters.likelihood_init_fraction[1]) >> 4);
      likelihood_levels.level[2] =
          likelihood_levels.level[3] - ((weak_minus_strong * ldpc_decoder_parameters.likelihood_init_fraction[2]) >> 4);
    }
  }
  likelihood_levels.min = likelihood_levels.strong;
  return likelihood_levels;
}
int ldpc_packet::f_update_vn_post(int likelihood, int weight, int min_likelihood, int max_likelihood, bool post_process,
                                  bool post_process2, bool be_aggressive, int flip_threshold, bool pushing) {
  int likelihood_new;
  bool do_post_flipped;
  bool do_post_unflipped;

  // ========== 方案 T：迭代自适应步长 ==========
  // 前期用半步防发散，后期（be_aggressive=true）用全步突破 error floor
  if (VN_BITS <= 2) {
    bool flipped = (likelihood >= flip_threshold);
    int delta;

	    if (!flipped) {
	      if (0 && be_aggressive) {  // 禁用，只靠 aggr 权重放大
	        delta = weight;  // 后期全步进攻
	      } else {
	        delta = (weight + 1) >> 1;  // 前期半步
	      }
	    } else {
	      delta = (weight >> 1) + ((weight & 1) && pushing ? 1 : 0);
	    }

    likelihood_new = flipped ? (likelihood - delta) : (likelihood + delta - 1);

    do_post_flipped = post_process && (likelihood_new == flip_threshold);
		    do_post_unflipped = post_process && (likelihood_new < flip_threshold);
		    if (do_post_unflipped)
		      likelihood_new = flip_threshold - 1;
		    else if (do_post_flipped)
		      likelihood_new = flip_threshold + 1;

    if (likelihood_new < min_likelihood) likelihood_new = min_likelihood;
    if (likelihood_new > max_likelihood) likelihood_new = max_likelihood;

    return likelihood_new;
  }
  // ========== 3-bit 及以上：原逻辑 ==========
  bool flipped = (likelihood >= flip_threshold);
  int delta = weight;

  likelihood_new = flipped ? (likelihood - delta) : (likelihood + delta - 1);

  do_post_flipped = post_process && (likelihood_new == flip_threshold);
  do_post_unflipped = post_process && (likelihood_new < flip_threshold);
  if (do_post_unflipped)
    likelihood_new = flip_threshold - 1;
  else if (do_post_flipped)
    likelihood_new = flip_threshold + 1;

  if (post_process2 && !flipped && (weight == 1) && (likelihood_new == (flip_threshold - 1)))
    likelihood_new++;
  if (likelihood_new <= min_likelihood)
    likelihood_new = min_likelihood;
  if (likelihood_new >= max_likelihood)
    likelihood_new = max_likelihood;

  return likelihood_new;
}
// int ldpc_packet::f_update_vn_post(int likelihood, int weight, int min_likelihood, int max_likelihood, bool post_process,
//                                   bool post_process2, bool be_aggressive, int flip_threshold, bool look,
//                                   bool scale2x) {
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
void ldpc_packet::ldpc_config(int m, int n, int sc, int st, int wt, char *pchk_file) {
  // QC matrix sizes
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

  printf("[LDPC] Configuring LDPC code with m=%d, n=%d, sc=%d, st=%d, col_wt=%d\n", bm_m, bm_n, cir_sz, tm_sz, col_wt);

  // read parity check matrix
  if (sc == 256) {
    ldpc_rd_phck(pchk_file);

    // G matrices
    qc_a = mod2sparse_allocate(tm_sz * cir_sz, (bm_n - bm_m) * cir_sz);
    qc_b = mod2sparse_allocate(tm_sz * cir_sz, (bm_m - tm_sz) * cir_sz);
    qc_c = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_n - bm_m) * cir_sz);
    qc_d = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);
    qc_e = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, tm_sz * cir_sz);
    qc_fi = mod2sparse_allocate((bm_m - tm_sz) * cir_sz, (bm_m - tm_sz) * cir_sz);

    // generate G matrices
    ldpc_gen_gm();
  } else if (sc == 512) {
    int VERBOSITY = 0;
    int i;
    int j;
    int k;
    int bit;
    int h_matrix_index[13];

    h_matrix.rows = m;
    h_matrix.cols = n;
    h_matrix.bits = sc;
    h_matrix.bytes_of_userdata = info_len / 8;
    h_matrix.bytes_of_parity = (blk_len - info_len) / 8;
    h_matrix.min_rows = 5;
    h_matrix.max_rows = 13;

    h_matrix.unused_bytes_of_parity =
        (h_matrix.rows * (h_matrix.bits >> 3)) - h_matrix.bytes_of_parity; // padding bytes
    if (h_matrix.unused_bytes_of_parity != 0)
      h_matrix.extra_bytes_of_parity =
          (h_matrix.bits >> 3) - h_matrix.unused_bytes_of_parity; // fractional_parity_bytes
    else
      h_matrix.extra_bytes_of_parity = 0;

    h_matrix.unused_bytes_of_userdata =
        ((h_matrix.cols - h_matrix.rows) * (h_matrix.bits >> 3)) - h_matrix.bytes_of_userdata;
    if (h_matrix.unused_bytes_of_userdata != 0)
      h_matrix.extra_bytes_of_userdata = (h_matrix.bits >> 3) - h_matrix.unused_bytes_of_userdata;
    else
      h_matrix.extra_bytes_of_userdata = 0;

    h_matrix.extra_bits_of_parity = h_matrix.extra_bytes_of_parity << 3;
    h_matrix.extra_bits_of_userdata = h_matrix.extra_bytes_of_userdata << 3;
    if (VERBOSITY > 0)
      printf("### MATRIX: BITS: %4d ROWS: %2d COLS: %3d BYTES_OF_USERDATA: %5d BYTES_OF_PARITY: %4d "
             "EXTRA_BYTES_OF_USERDATA: %5d EXTRA_BYTES_OF_PARITY: %4d\n",
             h_matrix.bits, h_matrix.rows, h_matrix.cols, h_matrix.bytes_of_userdata, h_matrix.bytes_of_parity,
             h_matrix.extra_bytes_of_userdata, h_matrix.extra_bytes_of_parity);

    int mx[13][80];
    int my[13][80];
    int rw[13][80] = {};
    int rw_max[13 + 1];
    int rw_min[13 + 1];
    float rw_avg[13 + 1];
    int num_reduced;
    bool ldpc_matrix_occupied[13][13][80];
    int ldpc_matrix_fade[13][13][80];
    for (k = 0; k < 13; k++) {
      for (i = 0; i < 13; i++) {
        for (j = 0; j < 80; j++) {
          ldpc_matrix_occupied[k][i][j] = 0;
          ldpc_matrix_fade[k][i][j] = 0;
        }
      }
    }

    for (i = 5 - 1; i < 13; i++) {
      rw_avg[i] = (float)(4 * (67 + i + 1)) / (float)(i + 1);
      rw_min[i] = int(rw_avg[i]);
      rw_max[i] = int(rw_avg[i] + 0.999999);
      if (VERBOSITY > 0)
        printf("### FOR MATRIX %2d, ROW_WEIGHT_AVG: %8.4f = %2d %d/%d   MIN: %2d MAX: %2d\n", i, rw_avg[i],
               int(rw_avg[i]), int(rw_avg[i] * (i + 1)) % (i + 1), (i + 1), rw_min[i], rw_max[i]);
    }

    srand(1);

    bool use_location;
    for (i = 0; i < 5; i++) {
      for (j = 75; j < 80; j++) {
        k = 80 - 1 - j;
        use_location = (k % 5) != i;
        if ((i == (5 - 1)) && (j == (80 - 1)))
          use_location = 0; // Make matrix invertible
        if (use_location) {
          ldpc_matrix_occupied[5 - 1][i][j] = 1;
          rw[5 - 1][i]++;
        }
      }
    } // Generate 5 x 5 sub-matrix, same in all matrix

    for (i = 0; i < 5; i++) {
      for (j = 0; j < 67; j++) {
        k = 80 - 1 - j;
        if (j == 0)
          use_location = (i != 0);
        else if (j == 1)
          use_location = (i != 1);
        else if (j == 2)
          use_location = (i != 2);
        else if (j == 3)
          use_location = (i != 3);
        else if (j == 4)
          use_location = (i != 4);
        else if (j == 5)
          use_location = (i != 1);
        else if (j == 6)
          use_location = (i != 3);
        else if (j == 7)
          use_location = (i != 4);
        else if (j == 8)
          use_location = (i != 2);
        else if (j == 9)
          use_location = (i != 0);
        else if (j == 10)
          use_location = (i != 4);
        else if (j == 11)
          use_location = (i != 3);
        else if (j == 12)
          use_location = (i != 1);
        else if (j == 13)
          use_location = (i != 0);
        else if (j == 14)
          use_location = (i != 2);
        else if (j == 15)
          use_location = (i != 3);
        else if (j == 16)
          use_location = (i != 1);
        else if (j == 17)
          use_location = (i != 2);
        else if (j == 18)
          use_location = (i != 0);
        else if (j == 19)
          use_location = (i != 4);
        else if (j == 20)
          use_location = (i != 0);
        else if (j == 21)
          use_location = (i != 0);
        else if (j == 22)
          use_location = (i != 4);
        else if (j == 23)
          use_location = (i != 4);
        else if (j == 24)
          use_location = (i != 2);
        else if (j == 25)
          use_location = (i != 2);
        else if (j == 26)
          use_location = (i != 3);
        else if (j == 27)
          use_location = (i != 3);
        else if (j == 28)
          use_location = (i != 1);
        else if (j == 29)
          use_location = (i != 1);
        else if (j == 30)
          use_location = (i != 3);
        else if (j == 31)
          use_location = (i != 0);
        else if (j == 32)
          use_location = (i != 3);
        else if (j == 33)
          use_location = (i != 1);
        else if (j == 34)
          use_location = (i != 4);
        else if (j == 35)
          use_location = (i != 0);
        else if (j == 36)
          use_location = (i != 4);
        else if (j == 37)
          use_location = (i != 2);
        else if (j == 38)
          use_location = (i != 2);
        else if (j == 39)
          use_location = (i != 1);
        else if (j == 40)
          use_location = (i != 1);
        else if (j == 41)
          use_location = (i != 4);
        else if (j == 42)
          use_location = (i != 3);
        else if (j == 43)
          use_location = (i != 2);
        else if (j == 44)
          use_location = (i != 0);
        else if (j == 45)
          use_location = (i != 2);
        else if (j == 46)
          use_location = (i != 4);
        else if (j == 47)
          use_location = (i != 0);
        else if (j == 48)
          use_location = (i != 1);
        else if (j == 49)
          use_location = (i != 3);
        else if (j <= 59)
          use_location = ((k + 4) % 5) != i;
        else if (j <= 69)
          use_location = ((k + 2) % 5) != i;
        if (use_location) {
          ldpc_matrix_occupied[5 - 1][i][j] = 1;
          rw[5 - 1][i]++;
        }
      }
    }

    if (VERBOSITY > 0) {
      for (i = 0; i < 13; i++) {
        k = 5 - 1;
        printf("### MATRIX: %2d ROW: %2d WEIGHT: %2d OCCUPIED: ", k, i, rw[k][i]);
        for (j = 0; j < 80; j++)
          printf("%1x ", ldpc_matrix_occupied[k][i][j]);
        printf("\n");
      }
      printf("\n");
    }

    for (i = 0; i < 13; i++)
      if (rw[5 - 1][i] > rw_max[5 - 1])
        rw_max[5 - 1] = rw[5 - 1][i];

    for (k = 5; k < 13; k++) {
      for (i = 0; i < 13; i++) {
        rw[k][i] = 0;
        for (j = 0; j < 80; j++) {
          ldpc_matrix_occupied[k][i][j] = ldpc_matrix_occupied[k - 1][i][j];
          rw[k][i] += ldpc_matrix_occupied[k][i][j];
        }
      }

      for (i = 0; i < 5; i++) {
        ldpc_matrix_occupied[k][k - i][80 - k - 1] = 1;
        rw[k][k - i]++;
      }
      while (rw[k][k] < rw_min[k]) {
        i = (rand() % (67 + 0));
        if (ldpc_matrix_occupied[k][k][i] == 0) {
          ldpc_matrix_occupied[k][k][i] = 1;
          rw[k][k]++;
        }
      }

      for (j = 0; j < 80; j++) {
        if (ldpc_matrix_occupied[k][k][j]) {
          int rw_max_in_col = 0;
          for (i = 0; i < k; i++)
            if (ldpc_matrix_occupied[k][i][j] && (rw[k][i] > rw_max_in_col))
              rw_max_in_col = rw[k][i];

          bool found_it = 0;
          for (i = 0; i < k; i++) {
            if (ldpc_matrix_occupied[k][i][j] && (rw[k][i] == rw_max_in_col) && !found_it) {
              found_it = 1;
              rw[k][i]--;
              ldpc_matrix_occupied[k][i][j] = 0;
              ldpc_matrix_fade[k][i][j] = 1;
            }
          }
        }
      }

      if (VERBOSITY > 0) {
        for (i = 0; i < 13; i++) {
          printf("### MATRIX: %2d ROW: %2d WEIGHT: %2d OCCUPIED: ", k, i, rw[k][i]);
          for (j = 0; j < 80; j++)
            printf("%1x ", ldpc_matrix_occupied[k][i][j]);
          printf("\n");
        }
      }
      printf("\n");
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

    h_matrix_index[0] = 0 + (3 * h_matrix.delta[0]);
    h_matrix_index[1] = 0 + (4 * h_matrix.delta[1]);
    h_matrix_index[2] = 0 + (4 * h_matrix.delta[2]);
    h_matrix_index[3] = 0 + (4 * h_matrix.delta[3]);
    h_matrix_index[4] = 0 + (3 * h_matrix.delta[4]);
    h_matrix_index[5] = 0;
    h_matrix_index[6] = 0;
    h_matrix_index[7] = 0;
    h_matrix_index[8] = 0;
    h_matrix_index[9] = 0;
    h_matrix_index[10] = 0;
    h_matrix_index[11] = 0;
    h_matrix_index[12] = 0;

    for (i = 0; i < h_matrix.rows; i++)
      h_matrix.last_element[i] = h_matrix_index[i];

    for (i = 0; i < h_matrix.rows; i++) {
      for (j = 0; j < h_matrix.cols; j++) {
        if (j < (h_matrix.cols - h_matrix.rows)) {
          h_matrix.occupied[i][j] = ldpc_matrix_occupied[h_matrix.rows - 1][i][j];
          h_matrix.fade[i][j] = ldpc_matrix_fade[h_matrix.rows - 1][i][j];
        } else {
          h_matrix.occupied[i][j] = ldpc_matrix_occupied[h_matrix.rows - 1][i][80 - h_matrix.cols + j];
          h_matrix.fade[i][j] = ldpc_matrix_fade[h_matrix.rows - 1][i][80 - h_matrix.cols + j];
        }
      }
    }

    for (i = 0; i < h_matrix.rows; i++)
      h_matrix.row_weight[i] = 0;

    for (j = 0; j < h_matrix.cols; j++)
      h_matrix.col_weight[j] = 0;

    for (j = 0; j < h_matrix.cols; j++) {
      n = h_matrix.cols - 1 - j;
      k = n;
      h_matrix.parity_column[j] = (j < h_matrix.rows);
      h_matrix.bits_in_last_column = 8 * 20;
      for (i = 0; i < h_matrix.rows; i++) {
        if (h_matrix.occupied[i][n] || h_matrix.fade[i][n]) {
          h_matrix.element[i][k] = h_matrix_index[i];
          h_matrix.first_element[i] = h_matrix.element[i][k];
          h_matrix_index[i] = (h_matrix_index[i] + h_matrix.bits - h_matrix.delta[i]) % h_matrix.bits;
          if (h_matrix.occupied[i][n]) {
            h_matrix.row_weight[i]++;
            h_matrix.col_weight[k]++;
          }
        } else
          h_matrix.element[i][k] = -1;
      }
    }
    for (i = 0; i < h_matrix.rows; i++) {
      // h_matrix.wraparound[i] =
      //     (h_matrix.bits + h_matrix.last_element[i] - h_matrix.first_element[i]) % h_matrix.bits; // should delete?
      h_matrix.wraparound[i] = (h_matrix.bits + h_matrix.first_element[i] - h_matrix.last_element[i]) % h_matrix.bits;
    }

    for (j = 0; j < 80; j++)
      for (k = 0; k < h_matrix.bits; k++)
        h_matrix.mask[j][k] = 0;

    if (h_matrix.extra_bits_of_parity > 0) {
      for (j = 0; j < h_matrix.cols; j++) {
        if (h_matrix.occupied[h_matrix.rows - 1][j]) {
          for (k = 0; k < h_matrix.bits; k++) {
            bit = (k + h_matrix.bits - h_matrix.element[h_matrix.rows - 1][j]) % h_matrix.bits;
            if ((bit < h_matrix.extra_bits_of_parity) || (h_matrix.extra_bits_of_parity == 0))
              h_matrix.mask[j][k] = 1;
          }
        }
      }
    }
    if (VERBOSITY > 0)
      f_print_h_matrix(h_matrix);
  }
} // ldpc_config

void ldpc_packet::ldpc_ibex_input(int syndrome_cal_only = 0, int max_iter = 0, int post_iter = 0) {
  ldpc_decoder_input.post_iteration = post_iter;
  ldpc_decoder_input.iteration_limit = max_iter;
  ldpc_decoder_input.nand_strobes = (rd_num == 1) ? 0 : (rd_num == 3) ? 1 : (rd_num == 5) ? 2 : 3;
  int soft_bits = (rd_num == 1) ? 0 : (rd_num == 3) ? 1 : 2;
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
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = ldpc_decoder_parameters.likelihood_map[3];
        if (soft_bits >= 2)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = ldpc_decoder_parameters.likelihood_map[3];
        else if (soft_bits >= 1)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable; // replicate
      } else if ((j == (h_matrix.cols - h_matrix.rows)) &&
                 (k >= (h_matrix.bits - 8 * h_matrix.unused_bytes_of_parity))) {
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = 0;
        if (soft_bits >= 1)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = ldpc_decoder_parameters.likelihood_map[3];
        if (soft_bits >= 2)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = ldpc_decoder_parameters.likelihood_map[3];
        else if (soft_bits >= 1)
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable; // replicate
      } else {
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = (rx_blk[bit_index] >= 0) ? 0 : 1;
        // debug
        err_cnt += tx_blk[bit_index] != ((rx_blk[bit_index] >= 0) ? 0 : 1);
        if (rx_blk[bit_index] >= 0) {
          bin_ibex = (det_blk[bit_index] + rd_num - 1) % (rd_num + 1);
          if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = (bin_ibex / 2) % 2;
          if (soft_bits >= 2)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = (bin_ibex / 4) % 2;
          else if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
                ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable; // replicate
        } else {
          bin_ibex = det_blk[bit_index];
          if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = (bin_ibex / 2) % 2;
          if (soft_bits >= 2)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = (bin_ibex / 4) % 2;
          else if (soft_bits >= 1)
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 =
                ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable; // replicate
        }
        bit_index++;
      }
    }
  }

  if (err_cnt != raw_err_num)
    printf("### LDPC IBEX INPUT ERROR: ERR_CNT %d != RAW_ERR_NUM %d\n", err_cnt, raw_err_num);
}

void ldpc_packet::ldpc_ibex_parameters(int post_process_en = 1, int syndrome_weight_thr_qc = 0,
                                       int syndrome_weight_thr_post = 0, int early_terminate_dis = 0) {
  ldpc_decoder_parameters.post_process_en = post_process_en;
  ldpc_decoder_parameters.syndrome_weight_thr_qc = syndrome_weight_thr_qc;
  ldpc_decoder_parameters.syndrome_weight_thr_post = syndrome_weight_thr_post;
  ldpc_decoder_parameters.likelihood_thr = 64;
  ldpc_decoder_parameters.post_ratio = 12;

  ldpc_decoder_parameters.likelihood_init_coef_all[0][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[0][1] = 7;
  ldpc_decoder_parameters.likelihood_init_coef_all[0][2] = 6;
  ldpc_decoder_parameters.likelihood_init_coef_all[0][3] = 5;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][1] = 6;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][2] = 5;
  ldpc_decoder_parameters.likelihood_init_coef_all[1][3] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][1] = 5;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][2] = 5;
  ldpc_decoder_parameters.likelihood_init_coef_all[2][3] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][1] = 5;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][2] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[3][3] = 3;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][1] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][2] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[4][3] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][1] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][2] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[5][3] = 3;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][1] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][2] = 3;
  ldpc_decoder_parameters.likelihood_init_coef_all[6][3] = 3;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][0] = 0;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][1] = 4;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][2] = 3;
  ldpc_decoder_parameters.likelihood_init_coef_all[7][3] = 2;

  ldpc_decoder_parameters.likelihood_init_fraction[1] = 8;
  ldpc_decoder_parameters.likelihood_init_fraction[2] = 4;

  ldpc_decoder_parameters.likelihood_map[0] = 3;
  ldpc_decoder_parameters.likelihood_map[1] = 2;
  ldpc_decoder_parameters.likelihood_map[2] = 1;
  ldpc_decoder_parameters.likelihood_map[3] = 0;
  ldpc_decoder_parameters.likelihood_map[4] = 3;
  ldpc_decoder_parameters.likelihood_map[5] = 2;
  ldpc_decoder_parameters.likelihood_map[6] = 1;
  ldpc_decoder_parameters.likelihood_map[7] = 0;

  ldpc_decoder_parameters.early_terminate_dis = early_terminate_dis;
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
  unsigned int *c_soft_data_table = new unsigned int[8]; // 32 bit size / 4 bits each field = 8
  for (int i = 0; i < 8; i++) {
    int shift_bit = i * 4; // 4 bit at a time
    c_soft_data_table[i] = (soft_bit_control_0 >> shift_bit) & 0xf;
  }
  ldpc_decoder_parameters.soft_bit_table[0] = c_soft_data_table[0];
  ldpc_decoder_parameters.soft_bit_table[1] = c_soft_data_table[1];
  ldpc_decoder_parameters.soft_bit_table[2] = c_soft_data_table[2];
  ldpc_decoder_parameters.soft_bit_table[3] = c_soft_data_table[3];
  ldpc_decoder_parameters.soft_bit_table[4] = c_soft_data_table[4];
  ldpc_decoder_parameters.soft_bit_table[5] = c_soft_data_table[5];
  ldpc_decoder_parameters.soft_bit_table[6] = c_soft_data_table[6];
  ldpc_decoder_parameters.soft_bit_table[7] = c_soft_data_table[7];
}

// LDPC decoder config
void ldpc_packet::ldpc_dec_config(int max_fdec_itr, int fdec_col_skip, int max_ldec_itr, float dec_alpha, int fin_mode,
                                  int fin_q_num, int fin_r_num, int fin_f_num) {
  // syndrome weight
  init_synd_wt_min = hm_m;
  init_synd_wt_max = 0;

  // layer config
  if (max_ldec_itr > 0) {
    ldec_max_itr = max_ldec_itr;
    ldec_early_term_en = 1;
  } else {
    ldec_max_itr = -1 * max_ldec_itr;
    ldec_early_term_en = 0;
  }

  alpha = dec_alpha;

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
  printf("[LDPC DEBUG] Q MSG: %d/%d, max %f, min %f \n", finite_q_num, fin_f_num, finite_q_max, finite_q_min);
  printf("[LDPC DEBUG] R MSG: %d/%d, max %f, min %f \n", finite_r_num, fin_f_num, finite_r_max, finite_r_min);
  printf("[LDPC DEBUG] C MSG: %d/%d, max %f, min %f \n", finite_c_num - 1, fin_f_num, finite_c_max, finite_c_min);
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
  if (qc_bm) mod2sparse_free(qc_bm);
  if (qc_hm) mod2sparse_free(qc_hm);
  if (qc_a) mod2sparse_free(qc_a);
  if (qc_b) mod2sparse_free(qc_b);
  if (qc_c) mod2sparse_free(qc_c);
  if (qc_d) mod2sparse_free(qc_d);
  if (qc_e) mod2sparse_free(qc_e);
  if (qc_fi) mod2sparse_free(qc_fi);

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
  int skip_parity_bits_in_first_parity_column = (h_matrix.rows * h_matrix.bits) - parity_bits;
  int last_parity_bits_in_first_parity_column = h_matrix.bits - 1 - skip_parity_bits_in_first_parity_column;
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

  s_hard_codeword ldpc_encoder_output = f_ldpc_encode(ldpc_encoder_input, h_matrix);

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

s_hard_codeword ldpc_packet::f_ldpc_encode(s_hard_codeword ldpc_encoder_input, s_h_matrix h_matrix) {
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

  bool payload[67][512];
  bool parity[13][512];
  bool check_node[13][512];
  bool ldpc_matrix_occupied[13][80];
  bool ldpc_matrix_fade[13][80];
  unsigned int ldpc_matrix[13][80];
  int ldpc_encoder_failure;
  int bit_rotated;
  int syndrome_weight;

#include "ldpc_matrix_inverse.h"

  int num_bytes = h_matrix.bits >> 3;
  int num_payload_cols = int((h_matrix.bytes_of_userdata + num_bytes - 1) / num_bytes);
  int num_payload_bits = 8 * h_matrix.bytes_of_userdata;
  int num_parity_cols = h_matrix.rows;
  int extra_payload_cols = num_payload_cols - 64;
  int unused_parity_bytes = (h_matrix.rows * num_bytes) - h_matrix.bytes_of_parity;
  int unused_parity_bits = unused_parity_bytes * 8;
  int matrix_element[13];

  int bit_location;
  int byte_data;
  int index;

  // Set the initial VN to the input data (payload part will be correct, parity part will be 0)
  vn = ldpc_encoder_input;
  // compute the CN for the payload
  cn = f_check_nodes(h_matrix, vn);

  // Apply the first part of the parity matrix, between the user data and the last 5 columns

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
            if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) && h_matrix.mask[j][k])
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
              (cn.r[i].b[m] & encoder_matrix_first_col[j][i][(h_matrix.bits + m - k) % h_matrix.bits]);
        }
      }
    }
  }

  for (j = num_parity_cols - 5; j < num_parity_cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      if ((j > 0) || (unused_parity_bytes == 0) || ((j == 0) && (k >= unused_parity_bits)))
        vn.c[j + h_matrix.cols - h_matrix.rows].b[k] = parity[j][k];
    }
  }

  cn = f_check_nodes(h_matrix, vn);
  syndrome_weight = f_check_node_weight(h_matrix, cn);
  ldpc_encoder_failure = (syndrome_weight > 0);
  if (ldpc_encoder_failure != 0) {
    printf("ERROR: Encoding! Syndrome weight = %d > 0\n", syndrome_weight);
    exit(1);
  }

  return vn;
}

void ldpc_packet::ldpc_decoder(enum dec_model dec_mode) {
  // add 0 padding
  vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
  for (int i = 0; i < pad_len; i++)
    dec_di_blk[info_len + i] = max_llr_bin;

  vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);

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
  else if (dec_mode == BF_IBEX)
    ldpc_dec_bf_ibex(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);

  // remove padding
  if (dec_mode != BF_IBEX) {
    vec_copy(dec_do_blk, dec_blk, 0, 0, info_len);
    vec_copy(dec_do_blk, dec_blk, hm_k, info_len, hm_m);
  } else {
    vec_copy(dec_do_blk, dec_blk, 0, 0, info_len);
    if (h_matrix.extra_bits_of_parity > 0) {
      vec_copy(dec_do_blk, dec_blk, hm_k, info_len, h_matrix.extra_bits_of_parity);
      vec_copy(dec_do_blk, dec_blk, hm_k + cir_sz, info_len + h_matrix.extra_bits_of_parity, (bm_m - 1) * cir_sz);
    } else {
      vec_copy(dec_do_blk, dec_blk, hm_k, info_len, hm_m);
    }
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
  char *vn_hd_sel;   // current HD of selected column (from dec_do_blk) in VN order
  char *vn_raw_sel;  // raw data of selected column (from dec_di_blk) in VN order
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

  for (int itr = 0; (itr <= fdec_max_itr) && ((fdec_early_term_en == 0) || (cw_fail == 1)); itr++) {
    for (int i = 0; i < (itr == fdec_max_itr ? p_num : bm_n); i++) {
      if (itr == 0) {
        vec_copy(dec_di_blk, vn_flp_sel[p_num], i * cir_sz, 0, cir_sz);
        col_updt = i;
        itr_updt = itr;
        col_skip = false;
      } else {
        vec_clr(vn_synd_cnt, cir_sz);
        for (e = mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e)) {
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
            if (((vn_raw_sel[j] == vn_hd_sel[j]) && (vn_synd_cnt[j] >= flp_thrshd0[itr - 1])) ||
                ((vn_raw_sel[j] != vn_hd_sel[j]) && (vn_synd_cnt[j] >= flp_thrshd1[itr - 1]))) {
              vn_flp_sel[0][j] = 1;
            } else {
              vn_flp_sel[0][j] = 0;
            }
          }

          col_updt = vn_flp_col[p_num];
          itr_updt = vn_flp_itr[p_num];

          for (int j = 0; j < cir_sz; j++) {
            if (vn_flp_sel[p_num][j] == 1) {
              dec_do_blk[col_updt * cir_sz + j] = (dec_do_blk[col_updt * cir_sz + j] + 1) % 2;
            }
          }
        } // non-skipped columns(flip logic)
      }   // non-1st iteration columns
      fdec_cyc_org++;

      if (col_skip == false) {
        fdec_cyc_num++;

        for (e = mod2sparse_first_in_col(qc_bm, col_updt); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e)) {
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
  char *vn_sb_sel;   // soft bit of selected column (from dec_sb_blk) in VN order
  char *vn_hd_sel;   // current HD of selected column (from dec_do_blk) in VN order
  char *vn_raw_sel;  // raw data of selected column (from dec_di_blk) in VN order
  char *vn_flp_sel;  // flip flag of selected column in VN order
  char *cn_flp_sel;  // flip flag of selected column in CN order
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

  for (int itr = 0; (itr <= fdec_max_itr) && ((fdec_early_term_en == 0) || (cw_fail == 1)); itr++) {
    for (int i = 0; i < bm_n; i++) {
      if (itr == 0) {
        vec_copy(dec_di_blk, vn_flp_sel, i * cir_sz, 0, cir_sz);
        col_updt = i;
        itr_updt = itr;
      } else {
        vec_clr(vn_synd_cnt, cir_sz);
        for (e = mod2sparse_first_in_col(qc_bm, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e)) {
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

              if ((vn_synd_cnt[j] >= sb_thrshd0_s0[itr - 1]) && (vn_synd_cnt[j] < sb_thrshd0_s1[itr - 1]))
                vn_sb_sel[j] = 0;
            } else // weak
            {
              if (vn_synd_cnt[j] >= flp_thrshd0_w[itr - 1])
                vn_flp_sel[j] = 1;
              else
                vn_flp_sel[j] = 0;

              if ((vn_synd_cnt[j] < sb_thrshd0_w0[itr - 1]) || (vn_synd_cnt[j] >= sb_thrshd0_w1[itr - 1]))
                vn_sb_sel[j] = 1;
            }
          } else {
            if (vn_sb_sel[j] == 1) // strong
            {
              if (vn_synd_cnt[j] >= flp_thrshd1_s[itr - 1])
                vn_flp_sel[j] = 1;
              else
                vn_flp_sel[j] = 0;

              if ((vn_synd_cnt[j] >= sb_thrshd1_s0[itr - 1]) && (vn_synd_cnt[j] < sb_thrshd1_s1[itr - 1]))
                vn_sb_sel[j] = 0;
            } else // weak
            {
              if (vn_synd_cnt[j] >= flp_thrshd1_w[itr - 1])
                vn_flp_sel[j] = 1;
              else
                vn_flp_sel[j] = 0;

              if ((vn_synd_cnt[j] < sb_thrshd1_w0[itr - 1]) || (vn_synd_cnt[j] >= sb_thrshd1_w1[itr - 1]))
                vn_sb_sel[j] = 1;
            }
          }
        }

        itr_updt = itr;
        col_updt = i;

        for (int j = 0; j < cir_sz; j++) {
          if (vn_flp_sel[j] == 1) {
            dec_do_blk[col_updt * cir_sz + j] = (dec_do_blk[col_updt * cir_sz + j] + 1) % 2;
          }
        }

        vec_copy(vn_sb_sel, dec_sb_blk, 0, i * cir_sz, cir_sz);
      }

      // update syndrome memory
      // multiple times (col_wt)
      for (e = mod2sparse_first_in_col(qc_bm, col_updt); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e)) {
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
  char cmem_dump[50] = "./output/rdec_cmem_dump.txt";
  char stot_dump[50] = "./output/rdec_stot_dump.txt";
  char hdmem_dump[50] = "./output/rdec_hdmem_dump.txt";
  char log_dump[50] = "./output/rdec_log_dump.txt";
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
  float *cn_q_sel_pre;          // Q msg of the select circulant from previous layer
  float *cn_r_new_pre;          // New R msg in CN order of previous layer
  float *cn_app_pre;            // APP = Q + R_new in CN order of previous layer
  float *cn_app_cur;            // APP = Q + R_new in CN order of current layer
  float *cn_q_sel_cur;          // Q msg of the select circulant from current layer
  float *cn_r_old_cur;          // old R msg in CN order of previous layer
  float *cn_q_updt_cur;         // Updated Q msg of the select circulant in current layer
  int **cn_q_sign;              // Q sign

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

  for (int i = 0; i < bm_n; i++)
    for (int j = 0; j < cir_sz; j++)
      cn_q_mem[i][j] = (float)llr_tbl[dec_di_blk[i * cir_sz + j]];

  // iterative decoding
  for (int itr = 0; (itr <= ldec_max_itr) && ((ldec_early_term_en == 0) || (cw_fail == 1)); itr++) {
    // Q sign mem index
    cir_cnt = 0;

    // layer decoding
    for (int layer = 0; layer < bm_m && ((ldec_early_term_en == 0) || (cw_fail == 1)); layer++) {
      // initilize HD mem
      hd_init = (vec_sum(dec_init, bm_n) != 0);

#ifdef _LDPC_DEBUG_DUMP
      printf("[LDPC DEBUG] Layer decoding @ iteration %d, layer %d ...\n", itr, layer);
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
           !mod2sparse_at_end(e) && ((ldec_early_term_en == 0) || (cw_fail == 1)); e = mod2sparse_next_in_row(e)) {
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
            cn_r_new_pre[i] = cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
          else
            cn_r_new_pre[i] = cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;

          // APP in CN order of previous layer
          cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];

          // Quantization
          if (finite_mode == 1) {
            cn_app_pre[i] =
                (float)Sat_Quan((double)cn_app_pre[i], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
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
          vn_dec_hd[i] = cn_app_pre[(i + shift_val2 + cir_sz) % cir_sz] >= 0 ? 0 : 1;
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
            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
          else
            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

          // Q -= Rold
          cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          // Quantization
          if (finite_mode == 1) {
            cn_q_updt_cur[i] =
                (float)Sat_Quan((double)cn_q_updt_cur[i], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
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
        cn_c_mem[layer][i].min1_val = cn_c_updt_cur[i].min1_val * alpha;
        cn_c_mem[layer][i].min2_val = cn_c_updt_cur[i].min2_val * alpha;

        if (finite_mode == 1) {
          cn_c_mem[layer][i].min1_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min1_val, finite_c_max, finite_c_min,
                                                        finite_c_num, finite_f_num);
          cn_c_mem[layer][i].min2_val = (float)Sat_Quan((double)cn_c_mem[layer][i].min2_val, finite_c_max, finite_c_min,
                                                        finite_c_num, finite_f_num);
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
} // ldpc_dec_layer3

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
//                                    s_ldpc_decoder_parameters ldpc_decoder_parameters, s_h_matrix h_matrix) {
//   if (vn_bits <= 2) {
//     // 令 be_aggressive 判定为假，避免 2bit 下权重放大导致抖动
//     ldpc_decoder_parameters.likelihood_thr = 0;
//   }
//   ldpc_dec_bf_ibex(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);
// }

void ldpc_packet::ldpc_dec_bf_ibex(s_ldpc_decoder_input ldpc_decoder_input,
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
  s_check_nodes cn_shifted;
  s_likelihood_levels likelihood_levels;
  s_256_bits prng_256;
  s_512_bits prng_512;
  int prng_init[32];
  int syndrome_weight_delayed;

  int syndrome_weight_r[5];
  bool do_not_use_this_bit;
  bool look;

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
      hard_codeword.c[j].b[k] = ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard;
      vn.c[j].b[k].bit_hard = hard_codeword.c[j].b[k];
      soft_codeword.c[j].b[k] = ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable;
      vn.c[j].b[k].flipped = 0;
    }
  }

  cn = f_check_nodes(h_matrix, hard_codeword);
  if (VERBOSITY > 0) {
    printf("[LDPC DEBUG] Starting BF decoding with max %d iterations.\n", fdec_max_itr);
    f_print_hard_codeword(hard_codeword, h_matrix.cols, h_matrix.bits);
    printf("### DECODER C++: H MATRIX:\n");
    f_print_h_matrix(h_matrix);
    printf("### DECODER C++: AFTER SYNDROME CHECK, CODEWORD HAS CHECK NODES:\n");
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

  // Tunable thresholds for 2-bit aggr gate (optional env override for quick sweeps).
  const char *env_aggr_iter_hi = getenv("IBEX_AGGR_ITER_HI");
  const char *env_aggr_iter_lo = getenv("IBEX_AGGR_ITER_LO");
  const char *env_aggr_synd_th = getenv("IBEX_AGGR_SYND_TH");
  const char *env_aggr_strong_synd_th = getenv("IBEX_AGGR_STRONG_SW_TH");
  const int aggr_iter_hi = env_aggr_iter_hi ? atoi(env_aggr_iter_hi) : 220;
  const int aggr_iter_lo = env_aggr_iter_lo ? atoi(env_aggr_iter_lo) : 110;
  const int aggr_synd_th = env_aggr_synd_th ? atoi(env_aggr_synd_th) : 280;
  const int aggr_strong_synd_th = env_aggr_strong_synd_th ? atoi(env_aggr_strong_synd_th) : (1 << 30);

  // Multi-phase diversity (no per-VN state):
  // Allow retry phases (phase>0) to use different aggr thresholds to explore alternative trajectories.
  // By default, retry phases inherit phase0 thresholds unless explicit env overrides are provided.
  const char *env_phase1_aggr_iter_hi = getenv("IBEX_PHASE1_AGGR_ITER_HI");
  const char *env_phase1_aggr_iter_lo = getenv("IBEX_PHASE1_AGGR_ITER_LO");
  const char *env_phase1_aggr_synd_th = getenv("IBEX_PHASE1_AGGR_SYND_TH");
  const char *env_phase1_aggr_strong_synd_th = getenv("IBEX_PHASE1_AGGR_STRONG_SW_TH");
  int phase1_aggr_iter_hi = env_phase1_aggr_iter_hi ? atoi(env_phase1_aggr_iter_hi) : -1;
  int phase1_aggr_iter_lo = env_phase1_aggr_iter_lo ? atoi(env_phase1_aggr_iter_lo) : -1;
  int phase1_aggr_synd_th = env_phase1_aggr_synd_th ? atoi(env_phase1_aggr_synd_th) : -1;
  int phase1_aggr_strong_synd_th =
      env_phase1_aggr_strong_synd_th ? atoi(env_phase1_aggr_strong_synd_th) : -1;
  if (phase1_aggr_iter_hi < 0) phase1_aggr_iter_hi = -1;
  if (phase1_aggr_iter_lo < 0) phase1_aggr_iter_lo = -1;
  if (phase1_aggr_synd_th < 0) phase1_aggr_synd_th = -1;
  if (phase1_aggr_strong_synd_th < 0) phase1_aggr_strong_synd_th = -1;
  // Optional phase2 overrides (for restart_phases >= 3). If unset, phase2 inherits phase1 behavior.
  const char *env_phase2_aggr_iter_hi = getenv("IBEX_PHASE2_AGGR_ITER_HI");
  const char *env_phase2_aggr_iter_lo = getenv("IBEX_PHASE2_AGGR_ITER_LO");
  const char *env_phase2_aggr_synd_th = getenv("IBEX_PHASE2_AGGR_SYND_TH");
  const char *env_phase2_aggr_strong_synd_th = getenv("IBEX_PHASE2_AGGR_STRONG_SW_TH");
  int phase2_aggr_iter_hi = env_phase2_aggr_iter_hi ? atoi(env_phase2_aggr_iter_hi) : -1;
  int phase2_aggr_iter_lo = env_phase2_aggr_iter_lo ? atoi(env_phase2_aggr_iter_lo) : -1;
  int phase2_aggr_synd_th = env_phase2_aggr_synd_th ? atoi(env_phase2_aggr_synd_th) : -1;
  int phase2_aggr_strong_synd_th =
      env_phase2_aggr_strong_synd_th ? atoi(env_phase2_aggr_strong_synd_th) : -1;
  if (phase2_aggr_iter_hi < 0) phase2_aggr_iter_hi = -1;
  if (phase2_aggr_iter_lo < 0) phase2_aggr_iter_lo = -1;
  if (phase2_aggr_synd_th < 0) phase2_aggr_synd_th = -1;
  if (phase2_aggr_strong_synd_th < 0) phase2_aggr_strong_synd_th = -1;
  const char *env_2bit_mode = getenv("IBEX_2BIT_MODE");
  const int mode_2bit = env_2bit_mode ? atoi(env_2bit_mode) : 0;

  // Pure-2bit framework profiles (global control only; MUST NOT add any per-VN state).
  const char *env_profile = getenv("IBEX_PROFILE");
  int ibex_profile = env_profile ? atoi(env_profile) : 0;
  if (ibex_profile < 0)
    ibex_profile = 0;
  if (ibex_profile > 2)
    ibex_profile = 2;

  // UP-GDBF inspired "active iteration": delay stochastic escapes until later iterations / retry phases.
  const char *env_upgdbf_active_iter = getenv("IBEX_UPGDBF_ACTIVE_ITER");
  int upgdbf_active_iter_base =
      env_upgdbf_active_iter ? atoi(env_upgdbf_active_iter)
                             : ((ibex_profile >= 2) ? (ldpc_decoder_input.post_iteration + 32)
                                                    : ldpc_decoder_input.post_iteration);
  if (upgdbf_active_iter_base < 0)
    upgdbf_active_iter_base = 0;
  if (upgdbf_active_iter_base > ldpc_decoder_input.iteration_limit)
    upgdbf_active_iter_base = ldpc_decoder_input.iteration_limit;
  const char *env_upgdbf_active_phase_bonus = getenv("IBEX_UPGDBF_ACTIVE_PHASE_BONUS");
  int upgdbf_active_phase_bonus = env_upgdbf_active_phase_bonus
                                     ? atoi(env_upgdbf_active_phase_bonus)
                                     : ((ibex_profile >= 2) ? 16 : 0);
  if (upgdbf_active_phase_bonus < 0)
    upgdbf_active_phase_bonus = 0;
  if (upgdbf_active_phase_bonus > ldpc_decoder_input.iteration_limit)
    upgdbf_active_phase_bonus = ldpc_decoder_input.iteration_limit;
  const char *env_upgdbf_stall_early = getenv("IBEX_UPGDBF_STALL_EARLY");
  const bool upgdbf_stall_early =
      env_upgdbf_stall_early ? (atoi(env_upgdbf_stall_early) != 0) : (ibex_profile >= 2);

		  // A-direction knobs: make "pushing" meaningful and add controlled stochasticity for w=2 boost.
		  const char *env_push_dynamic = getenv("IBEX_PUSH_DYNAMIC");
		  const bool push_dynamic =
		      env_push_dynamic ? (atoi(env_push_dynamic) != 0) : (ibex_profile >= 1);
		  const char *env_push_mode = getenv("IBEX_PUSH_MODE");
		  int push_mode = env_push_mode ? atoi(env_push_mode) : 0;
		  if (push_mode < 0) push_mode = 0;
		  if (push_mode > 1) push_mode = 1;
		  const char *env_w2_stoch = getenv("IBEX_W2_STOCH");
		  const bool w2_stoch = env_w2_stoch ? (atoi(env_w2_stoch) != 0) : (ibex_profile >= 1);
		  const char *env_w2_stoch_xor = getenv("IBEX_W2_STOCH_XOR");
		  const bool w2_stoch_xor = env_w2_stoch_xor ? (atoi(env_w2_stoch_xor) != 0) : false;
		  const char *env_w2_boost_only_when_pushing = getenv("IBEX_W2_BOOST_ONLY_WHEN_PUSHING");
		  const bool w2_boost_only_when_pushing =
		      env_w2_boost_only_when_pushing ? (atoi(env_w2_boost_only_when_pushing) != 0) : (ibex_profile >= 1);
		  const char *env_w2_boost_weak_only = getenv("IBEX_W2_BOOST_WEAK_ONLY");
		  const bool w2_boost_weak_only =
		      env_w2_boost_weak_only ? (atoi(env_w2_boost_weak_only) != 0) : false;
		  const char *env_w2_boost_soft_only = getenv("IBEX_W2_BOOST_SOFT_ONLY");
		  const bool w2_boost_soft_only =
		      env_w2_boost_soft_only ? (atoi(env_w2_boost_soft_only) != 0) : false;
		  const char *env_w2_tail_guard = getenv("IBEX_W2_TAIL_GUARD");
		  const bool w2_tail_guard = env_w2_tail_guard ? (atoi(env_w2_tail_guard) != 0) : false;
		  // Syndrome-weight delta gate for w=2 boost (global only; no per-VN state):
		  // In retry phases, if the syndrome weight worsens too much between iterations, temporarily forbid w=2 boosting
		  // to avoid "adding fuel to the fire" on diverging trajectories.
		  const char *env_w2_boost_sw_delta_gate = getenv("IBEX_W2_BOOST_SW_DELTA_GATE");
		  const bool w2_boost_sw_delta_gate =
		      env_w2_boost_sw_delta_gate ? (atoi(env_w2_boost_sw_delta_gate) != 0) : false;
		  const char *env_w2_boost_sw_delta_hi = getenv("IBEX_W2_BOOST_SW_DELTA_HI");
		  int w2_boost_sw_delta_hi = env_w2_boost_sw_delta_hi ? atoi(env_w2_boost_sw_delta_hi) : 8;
		  if (w2_boost_sw_delta_hi < 0) w2_boost_sw_delta_hi = 0;
		  const char *env_w2_boost_sw_delta_min_iter = getenv("IBEX_W2_BOOST_SW_DELTA_MIN_ITER");
		  int w2_boost_sw_delta_min_iter =
		      env_w2_boost_sw_delta_min_iter ? atoi(env_w2_boost_sw_delta_min_iter)
		                                     : ldpc_decoder_input.post_iteration;
		  if (w2_boost_sw_delta_min_iter < 0)
		    w2_boost_sw_delta_min_iter = 0;
		  if (w2_boost_sw_delta_min_iter > ldpc_decoder_input.iteration_limit)
		    w2_boost_sw_delta_min_iter = ldpc_decoder_input.iteration_limit;
		  const char *env_toggle_strong = getenv("IBEX_TOGGLE_STRONG");
		  const bool toggle_strong = env_toggle_strong ? (atoi(env_toggle_strong) != 0) : false;
		  // Temporarily forbid any extra per-VN state beyond the 2-bit likelihood (pure 2bit route).
		  static constexpr bool k_allow_extra_vn_state = false;
		  const char *env_w2_cand_strong = getenv("IBEX_W2_CAND_STRONG");
		  const bool w2_cand_strong =
		      (k_allow_extra_vn_state && env_w2_cand_strong) ? (atoi(env_w2_cand_strong) != 0) : false;
		  const char *env_soft_guard = getenv("IBEX_SOFT_GUARD");
		  const bool soft_guard = env_soft_guard ? (atoi(env_soft_guard) != 0) : false;
		  const char *env_init_soft_bias = getenv("IBEX_INIT_SOFT_BIAS");
		  const bool init_soft_bias = env_init_soft_bias ? (atoi(env_init_soft_bias) != 0) : false;
		  const char *env_w1_stoch = getenv("IBEX_W1_STOCH");
		  const bool w1_stoch = env_w1_stoch ? (atoi(env_w1_stoch) != 0) : false;
		  const char *env_post_only_when_pushing = getenv("IBEX_POST_ONLY_WHEN_PUSHING");
		  const bool post_only_when_pushing =
		      env_post_only_when_pushing ? (atoi(env_post_only_when_pushing) != 0) : false;
		  const char *env_tabu1 = getenv("IBEX_TABU1");
		  const bool tabu1 = (k_allow_extra_vn_state && env_tabu1) ? (atoi(env_tabu1) != 0) : false;
		  const char *env_tabu_rev = getenv("IBEX_TABU_REV");
		  const bool tabu_rev = (k_allow_extra_vn_state && env_tabu_rev) ? (atoi(env_tabu_rev) != 0) : false;
				  const char *env_restart_phases = getenv("IBEX_RESTART_PHASES");
				  int restart_phases = env_restart_phases ? atoi(env_restart_phases) : ((ibex_profile >= 2) ? 3 : 1);
				  if (restart_phases < 1) restart_phases = 1;
				  if (restart_phases > 4) restart_phases = 4;
				  const char *env_phase1_w2_not_pushing = getenv("IBEX_PHASE1_W2_NOT_PUSHING");
				  const bool phase1_w2_not_pushing =
				      env_phase1_w2_not_pushing ? (atoi(env_phase1_w2_not_pushing) != 0) : false;
				  const char *env_phase1_w2_boost_to = getenv("IBEX_PHASE1_W2_BOOST_TO");
				  int phase1_w2_boost_to = env_phase1_w2_boost_to ? atoi(env_phase1_w2_boost_to) : 3;
				  if (phase1_w2_boost_to < 3) phase1_w2_boost_to = 3;
				  if (phase1_w2_boost_to > 7) phase1_w2_boost_to = 7;
				  const char *env_phase2_w2_not_pushing = getenv("IBEX_PHASE2_W2_NOT_PUSHING");
				  int phase2_w2_not_pushing_override = env_phase2_w2_not_pushing ? atoi(env_phase2_w2_not_pushing) : -1;
				  if (phase2_w2_not_pushing_override < 0) phase2_w2_not_pushing_override = -1;
				  const char *env_phase2_w2_boost_to = getenv("IBEX_PHASE2_W2_BOOST_TO");
				  int phase2_w2_boost_to = env_phase2_w2_boost_to ? atoi(env_phase2_w2_boost_to) : -1;
				  if (phase2_w2_boost_to >= 0) {
				    if (phase2_w2_boost_to < 3) phase2_w2_boost_to = 3;
				    if (phase2_w2_boost_to > 7) phase2_w2_boost_to = 7;
				  } else {
				    phase2_w2_boost_to = -1;
				  }
				  // FM-PGDBF style "flip probability gate" (phase2-only; no per-VN state).
				  const char *env_phase2_flip_rand_gate = getenv("IBEX_PHASE2_FLIP_RAND_GATE");
				  const bool phase2_flip_rand_gate =
				      env_phase2_flip_rand_gate ? (atoi(env_phase2_flip_rand_gate) != 0) : false;
				  const char *env_phase2_flip_rand_gate_max_w = getenv("IBEX_PHASE2_FLIP_RAND_GATE_MAX_W");
				  int phase2_flip_rand_gate_max_w = env_phase2_flip_rand_gate_max_w ? atoi(env_phase2_flip_rand_gate_max_w) : 3;
				  if (phase2_flip_rand_gate_max_w < 0) phase2_flip_rand_gate_max_w = 0;
				  const char *env_phase2_flip_rand_gate_gates = getenv("IBEX_PHASE2_FLIP_RAND_GATE_GATES");
				  int phase2_flip_rand_gate_gates = env_phase2_flip_rand_gate_gates ? atoi(env_phase2_flip_rand_gate_gates) : 2;
				  if (phase2_flip_rand_gate_gates < 1) phase2_flip_rand_gate_gates = 1;
				  if (phase2_flip_rand_gate_gates > 8) phase2_flip_rand_gate_gates = 8;
				  // Pure-2bit architectural exploration knobs (must NOT add any per-VN extra state).
				  // 1) Stall-triggered escape: when a phase stalls, make w=2 boosting more decisive in retry phases.
				  const char *env_stall_w2_esc = getenv("IBEX_STALL_W2_ESC");
				  const bool stall_w2_esc = env_stall_w2_esc ? (atoi(env_stall_w2_esc) != 0) : (ibex_profile >= 2);
			  const char *env_stall_w2_esc_iters = getenv("IBEX_STALL_W2_ESC_ITERS");
			  int stall_w2_esc_iters = env_stall_w2_esc_iters ? atoi(env_stall_w2_esc_iters) : 8;
			  if (stall_w2_esc_iters < 1) stall_w2_esc_iters = 1;
				  const char *env_stall_w2_esc_min_iter = getenv("IBEX_STALL_W2_ESC_MIN_ITER");
				  int stall_w2_esc_min_iter =
				      env_stall_w2_esc_min_iter ? atoi(env_stall_w2_esc_min_iter) : ldpc_decoder_input.post_iteration;
				  if (stall_w2_esc_min_iter < 0) stall_w2_esc_min_iter = 0;
				  if (stall_w2_esc_min_iter > ldpc_decoder_input.iteration_limit)
				    stall_w2_esc_min_iter = ldpc_decoder_input.iteration_limit;
				  // 1b) PPBF-like probabilistic escape (retry-only, no per-VN state): boost w=1/2 with p(E).
				  const char *env_ppbf_esc = getenv("IBEX_PPBF_ESC");
				  const bool ppbf_esc = env_ppbf_esc ? (atoi(env_ppbf_esc) != 0) : (ibex_profile >= 2);
				  const char *env_ppbf_esc_iters = getenv("IBEX_PPBF_ESC_ITERS");
				  int ppbf_esc_iters = env_ppbf_esc_iters ? atoi(env_ppbf_esc_iters) : 8;
				  if (ppbf_esc_iters < 1) ppbf_esc_iters = 1;
				  const char *env_ppbf_esc_min_iter = getenv("IBEX_PPBF_ESC_MIN_ITER");
				  int ppbf_esc_min_iter =
				      env_ppbf_esc_min_iter ? atoi(env_ppbf_esc_min_iter) : ldpc_decoder_input.post_iteration;
				  if (ppbf_esc_min_iter < 0) ppbf_esc_min_iter = 0;
				  if (ppbf_esc_min_iter > ldpc_decoder_input.iteration_limit)
				    ppbf_esc_min_iter = ldpc_decoder_input.iteration_limit;
				  const int stall_count_w2_min_iter = std::min(stall_w2_esc_min_iter, ppbf_esc_min_iter);

				  // 1c) Stall-triggered "mode window" (global only; no per-VN state):
				  // When a retry phase stalls, open a short window to temporarily boost escape actions
				  // (e.g., PPBF-like p(E) and/or column-level escape acceptance).
				  const char *env_mode_win = getenv("IBEX_MODE_WIN");
				  const bool mode_win = env_mode_win ? (atoi(env_mode_win) != 0) : false;
				  const char *env_mode_win_len = getenv("IBEX_MODE_WIN_LEN");
				  int mode_win_len = env_mode_win_len ? atoi(env_mode_win_len) : 16;
				  if (mode_win_len < 1) mode_win_len = 1;
				  if (mode_win_len > ldpc_decoder_input.iteration_limit)
				    mode_win_len = ldpc_decoder_input.iteration_limit;
				  const char *env_mode_win_trig_iters = getenv("IBEX_MODE_WIN_TRIG_ITERS");
				  int mode_win_trig_iters =
				      env_mode_win_trig_iters ? atoi(env_mode_win_trig_iters)
				                              : std::min(stall_w2_esc_iters, ppbf_esc_iters);
				  if (mode_win_trig_iters < 1) mode_win_trig_iters = 1;
				  const char *env_mode_win_ppbf_boost = getenv("IBEX_MODE_WIN_PPBF_BOOST");
				  const bool mode_win_ppbf_boost =
				      env_mode_win_ppbf_boost ? (atoi(env_mode_win_ppbf_boost) != 0) : true;
				  const char *env_mode_win_col_esc_boost = getenv("IBEX_MODE_WIN_COL_ESC_BOOST");
				  const bool mode_win_col_esc_boost =
				      env_mode_win_col_esc_boost ? (atoi(env_mode_win_col_esc_boost) != 0) : true;

				  // 1d) PPBF p(E) annealing (global schedule; no per-VN state):
				  // Make p(E) hotter right after post_iteration, colder later to reduce tail mis-flips.
				  const char *env_ppbf_anneal = getenv("IBEX_PPBF_ESC_ANNEAL");
				  const bool ppbf_anneal = env_ppbf_anneal ? (atoi(env_ppbf_anneal) != 0) : false;
				  const char *env_ppbf_anneal_iters = getenv("IBEX_PPBF_ESC_ANNEAL_ITERS");
				  int ppbf_anneal_iters = env_ppbf_anneal_iters ? atoi(env_ppbf_anneal_iters) : 64;
				  if (ppbf_anneal_iters < 0) ppbf_anneal_iters = 0;
				  if (ppbf_anneal_iters > ldpc_decoder_input.iteration_limit)
				    ppbf_anneal_iters = ldpc_decoder_input.iteration_limit;

				  // 1e) NGDBF-lite discrete noise bump (LFSR-gated; no per-VN state):
				  // In stall windows, occasionally add +1 to the tiny energy proxy to increase escape probability.
				  const char *env_ngdbf_noise = getenv("IBEX_NGDBF_NOISE");
				  const bool ngdbf_noise = env_ngdbf_noise ? (atoi(env_ngdbf_noise) != 0) : false;
				  const char *env_ngdbf_noise_gates = getenv("IBEX_NGDBF_NOISE_GATES");
				  int ngdbf_noise_gates = env_ngdbf_noise_gates ? atoi(env_ngdbf_noise_gates) : 3; // ~1/8
				  if (ngdbf_noise_gates < 1) ngdbf_noise_gates = 1;
				  if (ngdbf_noise_gates > 8) ngdbf_noise_gates = 8;
				  const char *env_ngdbf_noise_min_w = getenv("IBEX_NGDBF_NOISE_MIN_W");
				  int ngdbf_noise_min_w = env_ngdbf_noise_min_w ? atoi(env_ngdbf_noise_min_w) : 2;
				  if (ngdbf_noise_min_w < 0) ngdbf_noise_min_w = 0;
				  const char *env_ngdbf_noise_only_stall = getenv("IBEX_NGDBF_NOISE_ONLY_STALL");
				  const bool ngdbf_noise_only_stall =
				      env_ngdbf_noise_only_stall ? (atoi(env_ngdbf_noise_only_stall) != 0) : true;
				  const char *env_ngdbf_noise_stall_iters = getenv("IBEX_NGDBF_NOISE_STALL_ITERS");
				  int ngdbf_noise_stall_iters =
				      env_ngdbf_noise_stall_iters ? atoi(env_ngdbf_noise_stall_iters)
				                                  : std::min(stall_w2_esc_iters, ppbf_esc_iters);
				  if (ngdbf_noise_stall_iters < 1) ngdbf_noise_stall_iters = 1;
				  // 2) Randomized layered schedule: rotate bit-scan start within each column (optional retry-only).
				  const char *env_rotate_k = getenv("IBEX_ROTATE_K");
				  const bool rotate_k = env_rotate_k ? (atoi(env_rotate_k) != 0) : (ibex_profile >= 1);
			  const char *env_rotate_k_phase1_only = getenv("IBEX_ROTATE_K_PHASE1_ONLY");
			  const bool rotate_k_phase1_only =
			      env_rotate_k_phase1_only ? (atoi(env_rotate_k_phase1_only) != 0) : (ibex_profile >= 1);
			  // 3) Tail stabilization: cap the number of toggles per column to avoid cascade mis-flips.
			  const char *env_max_toggles_per_col = getenv("IBEX_MAX_TOGGLES_PER_COL");
			  int max_toggles_per_col = env_max_toggles_per_col ? atoi(env_max_toggles_per_col) : 0;
			  if (max_toggles_per_col < 0) max_toggles_per_col = 0;
			  const char *env_max_toggles_tail_only = getenv("IBEX_MAX_TOGGLES_TAIL_ONLY");
			  const bool max_toggles_tail_only =
			      env_max_toggles_tail_only ? (atoi(env_max_toggles_tail_only) != 0) : true;
			  const char *env_restart_prng_skip = getenv("IBEX_RESTART_PRNG_SKIP");
			  const int restart_prng_skip = env_restart_prng_skip ? atoi(env_restart_prng_skip) : 73;
			  const char *env_restart_on_stall = getenv("IBEX_RESTART_ON_STALL");
			  const bool restart_on_stall =
			      env_restart_on_stall ? (atoi(env_restart_on_stall) != 0) : (ibex_profile >= 2);
		  const char *env_stall_iters = getenv("IBEX_STALL_ITERS");
		  int stall_iters = env_stall_iters ? atoi(env_stall_iters) : 32;
		  if (stall_iters < 1) stall_iters = 1;
		  if (stall_iters > ldpc_decoder_input.iteration_limit) stall_iters = ldpc_decoder_input.iteration_limit;
		  const char *env_stall_min_iter = getenv("IBEX_STALL_MIN_ITER");
		  int stall_min_iter = env_stall_min_iter ? atoi(env_stall_min_iter) : 200;
			  if (stall_min_iter < 0) stall_min_iter = 0;
			  if (stall_min_iter > ldpc_decoder_input.iteration_limit) stall_min_iter = ldpc_decoder_input.iteration_limit;
			  const char *env_restart_relax_post_gate = getenv("IBEX_RESTART_RELAX_POST_GATE");
			  const bool restart_relax_post_gate =
			      env_restart_relax_post_gate ? (atoi(env_restart_relax_post_gate) != 0) : false;
			  const char *env_phase2_relax_post_gate = getenv("IBEX_PHASE2_RELAX_POST_GATE");
			  const bool phase2_relax_post_gate =
			      env_phase2_relax_post_gate ? (atoi(env_phase2_relax_post_gate) != 0) : false;
			  const char *env_restart_split_budget = getenv("IBEX_RESTART_SPLIT_BUDGET");
			  const bool restart_split_budget =
			      env_restart_split_budget ? (atoi(env_restart_split_budget) != 0) : (ibex_profile >= 2);
		  const char *env_smooth_fail = getenv("IBEX_SMOOTH_FAIL");
		  const bool smooth_fail = (k_allow_extra_vn_state && env_smooth_fail) ? (atoi(env_smooth_fail) != 0) : false;
		  const char *env_smooth_win = getenv("IBEX_SMOOTH_WIN");
		  int smooth_win = env_smooth_win ? atoi(env_smooth_win) : 64;
		  if (smooth_win <= 0) smooth_win = 64;
		  if (smooth_win > ldpc_decoder_input.iteration_limit) smooth_win = ldpc_decoder_input.iteration_limit;

			  // Column-level global escape + backtracking (no per-VN state; optional small per-column buffers allowed).
			  const char *env_col_global_esc = getenv("IBEX_COL_GLOBAL_ESC");
			  const bool col_global_esc = env_col_global_esc ? (atoi(env_col_global_esc) != 0) : false;
			  const char *env_col_global_esc_min_phase = getenv("IBEX_COL_GLOBAL_ESC_MIN_PHASE");
			  int col_global_esc_min_phase = env_col_global_esc_min_phase ? atoi(env_col_global_esc_min_phase) : 1;
			  if (col_global_esc_min_phase < 0) col_global_esc_min_phase = 0;
			  if (col_global_esc_min_phase > 4) col_global_esc_min_phase = 4;
			  const char *env_col_esc_target_hot = getenv("IBEX_COL_ESC_TARGET_HOT");
			  const bool col_esc_target_hot = env_col_esc_target_hot ? (atoi(env_col_esc_target_hot) != 0) : false;
			  const char *env_col_esc_target_hot_min_phase = getenv("IBEX_COL_ESC_TARGET_HOT_MIN_PHASE");
			  int col_esc_target_hot_min_phase =
			      env_col_esc_target_hot_min_phase ? atoi(env_col_esc_target_hot_min_phase) : 2;
			  if (col_esc_target_hot_min_phase < 0) col_esc_target_hot_min_phase = 0;
			  if (col_esc_target_hot_min_phase > 4) col_esc_target_hot_min_phase = 4;
			  const char *env_col_esc_target_hot_min_iter = getenv("IBEX_COL_ESC_TARGET_HOT_MIN_ITER");
			  int col_esc_target_hot_min_iter =
			      env_col_esc_target_hot_min_iter ? atoi(env_col_esc_target_hot_min_iter) : ldpc_decoder_input.post_iteration;
			  if (col_esc_target_hot_min_iter < 0) col_esc_target_hot_min_iter = 0;
			  if (col_esc_target_hot_min_iter > ldpc_decoder_input.iteration_limit)
			    col_esc_target_hot_min_iter = ldpc_decoder_input.iteration_limit;
			  const char *env_col_global_esc_iters = getenv("IBEX_COL_GLOBAL_ESC_ITERS");
			  int col_global_esc_iters = env_col_global_esc_iters ? atoi(env_col_global_esc_iters) : 8;
			  if (col_global_esc_iters < 1) col_global_esc_iters = 1;
		  const char *env_col_global_esc_min_iter = getenv("IBEX_COL_GLOBAL_ESC_MIN_ITER");
		  int col_global_esc_min_iter =
		      env_col_global_esc_min_iter ? atoi(env_col_global_esc_min_iter) : ldpc_decoder_input.post_iteration;
		  if (col_global_esc_min_iter < 0) col_global_esc_min_iter = 0;
		  if (col_global_esc_min_iter > ldpc_decoder_input.iteration_limit)
		    col_global_esc_min_iter = ldpc_decoder_input.iteration_limit;
		  const char *env_col_global_esc_max_toggles = getenv("IBEX_COL_GLOBAL_ESC_MAX_TOGGLES");
		  int col_global_esc_max_toggles = env_col_global_esc_max_toggles ? atoi(env_col_global_esc_max_toggles) : 2;
		  if (col_global_esc_max_toggles < 1) col_global_esc_max_toggles = 1;
		  if (col_global_esc_max_toggles > 8) col_global_esc_max_toggles = 8;
		  const char *env_col_global_esc_soft_only = getenv("IBEX_COL_GLOBAL_ESC_SOFT_ONLY");
		  const bool col_global_esc_soft_only =
		      env_col_global_esc_soft_only ? (atoi(env_col_global_esc_soft_only) != 0) : true;
		  const char *env_col_backtrack = getenv("IBEX_COL_BACKTRACK");
		  const bool col_backtrack = env_col_backtrack ? (atoi(env_col_backtrack) != 0) : true;
		  const char *env_col_backtrack_slack = getenv("IBEX_COL_BACKTRACK_SLACK");
		  int col_backtrack_slack = env_col_backtrack_slack ? atoi(env_col_backtrack_slack) : 0;
		  if (col_backtrack_slack < 0) col_backtrack_slack = 0;
		  const char *env_col_backtrack_require_improve = getenv("IBEX_COL_BACKTRACK_REQUIRE_IMPROVE");
		  const bool col_backtrack_require_improve =
		      env_col_backtrack_require_improve ? (atoi(env_col_backtrack_require_improve) != 0) : true;

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
        printf("### LDPC DECODER C++ EARLY TERMINATION: %5d\n", syndrome_weight);
    }
  }
  if (ldpc_decoder_input.syndrome_cal_only) {
    finished = 1;
  }

	  likelihood_levels =
	      f_likelihood_levels(ldpc_decoder_input.nand_strobes, ldpc_decoder_parameters, syndrome_weight, h_matrix.rows);

	  // soft_data -> likelihood_level
	  std::vector<int16_t> vn_likelihood_init;
	  vn_likelihood_init.resize(h_matrix.cols * h_matrix.bits);
	  bool soft_data[2];
	  for (j = 0; j < h_matrix.cols; j++) {
	    for (k = 0; k < h_matrix.bits; k++) {
	      soft_data[0] = ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable;
	      soft_data[1] = ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2;
	      if ((soft_data[1] == 0) && (soft_data[0] == 0))
	        vn.c[j].b[k].likelihood = likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[0]]; // 3 verilog: 3
	      if ((soft_data[1] == 0) && (soft_data[0] == 1))
	        vn.c[j].b[k].likelihood = likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[1]]; // 2 verilog: 2
	      if ((soft_data[1] == 1) && (soft_data[0] == 0))
	        vn.c[j].b[k].likelihood = likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[2]]; // 1 verilog: 3
	      if ((soft_data[1] == 1) && (soft_data[0] == 1))
	        vn.c[j].b[k].likelihood = likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[3]]; // 0 verilog: 2
	      if (init_soft_bias && (VN_BITS <= 2) && (ldpc_decoder_input.soft_bits > 0)) {
	        const bool soft_unreliable_init =
	            soft_data[0] || ((ldpc_decoder_input.soft_bits > 1) && soft_data[1]);
	        vn.c[j].b[k].likelihood = soft_unreliable_init ? 1 : 0;
	      }
	      vn_likelihood_init[j * h_matrix.bits + k] = static_cast<int16_t>(vn.c[j].b[k].likelihood);

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
           ldpc_decoder_input.corrupted_codeword.c[71].b[405].bit_questionable2, vn.c[71].b[405].likelihood);
  }

	  // 4x5 statistic table init
	  for (j = 0; j <= 3; j++) {
	    for (k = 0; k <= 4; k++) {
	      hard_codeword.errors_at_level_and_weight[j][k] = 0;
      hard_codeword.correct_at_level_and_weight[j][k] = 0;
	    }
	  }

			  std::vector<int16_t> last_toggle_iter;
			  if (tabu1 || tabu_rev) {
			    last_toggle_iter.assign(h_matrix.cols * h_matrix.bits, static_cast<int16_t>(-2));
			  }

			  int prng_skip_steps = 0;

			  for (int phase = 0; (phase < restart_phases) && (finished == 0) && (give_up == 0); phase++) {
			    std::vector<uint8_t> w2_cand_mem;
			    if (w2_cand_strong) {
			      w2_cand_mem.assign(h_matrix.cols * h_matrix.bits, 0);
			    }
			    std::vector<int8_t> smooth_sum;
			    int phase_iter_limit = ldpc_decoder_input.iteration_limit;
			    if (restart_split_budget && (restart_phases > 1)) {
			      const int total_budget = ldpc_decoder_input.iteration_limit;
			      const int base_budget = total_budget / restart_phases;
			      const int rem_budget = total_budget % restart_phases;
			      phase_iter_limit = base_budget + ((phase < rem_budget) ? 1 : 0);
			      if (phase_iter_limit < 1)
			        phase_iter_limit = 1;
			    }
			    const int effective_smooth_win = std::min(smooth_win, phase_iter_limit);
			    const int smooth_start_iter =
			        smooth_fail ? std::max(0, phase_iter_limit - effective_smooth_win) : (1 << 30);
			    if (phase > 0) {
			      // Restart from the exact initial VN likelihoods and no flips.
			      for (j = 0; j < h_matrix.cols; j++) {
			        for (k = 0; k < h_matrix.bits; k++) {
			          const int vn_idx = j * h_matrix.bits + k;
			          vn.c[j].b[k].likelihood = vn_likelihood_init[vn_idx];
			          vn.c[j].b[k].flipped = 0;
			        }
			      }
			      cn = f_check_nodes(h_matrix, hard_codeword);
			      syndrome_weight = f_check_node_weight(h_matrix, cn);
			      syndrome_weight_delayed = syndrome_weight;
			      for (i = 0; i < 5; i++)
			        syndrome_weight_r[i] = syndrome_weight;
			      finished = (syndrome_weight == 0);
			      give_up = 0;
			      if (tabu1 || tabu_rev)
			        std::fill(last_toggle_iter.begin(), last_toggle_iter.end(), static_cast<int16_t>(-2));
			      if (w2_cand_strong && !w2_cand_mem.empty())
			        std::fill(w2_cand_mem.begin(), w2_cand_mem.end(), 0);
			      iteration = 0;
			      clock_cycles = (2 * h_matrix.cols) + 1;
			    }
			    prng_skip_steps = (restart_prng_skip > 0) ? (phase * restart_prng_skip) : 0;
			    int best_sw_phase = syndrome_weight;
			    int stall_count = 0;
			    int best_sw_phase_w2 = syndrome_weight;
			    int stall_count_w2 = 0;
			    int prev_iter_sw = syndrome_weight;
				    bool pushing_iter = true;
				    bool w2_boost_sw_delta_allow = true;
				    // Retry-phase diversity (no per-VN state): optionally use different aggr thresholds in phase>0.
				    int aggr_iter_hi_eff = aggr_iter_hi;
				    int aggr_iter_lo_eff = aggr_iter_lo;
				    int aggr_synd_th_eff = aggr_synd_th;
				    int aggr_strong_synd_th_eff = aggr_strong_synd_th;
				    if (phase > 0) {
				      if (phase1_aggr_iter_hi >= 0) aggr_iter_hi_eff = phase1_aggr_iter_hi;
				      if (phase1_aggr_iter_lo >= 0) aggr_iter_lo_eff = phase1_aggr_iter_lo;
				      if (phase1_aggr_synd_th >= 0) aggr_synd_th_eff = phase1_aggr_synd_th;
				      if (phase1_aggr_strong_synd_th >= 0) aggr_strong_synd_th_eff = phase1_aggr_strong_synd_th;
				    }
				    if (phase >= 2) {
				      if (phase2_aggr_iter_hi >= 0) aggr_iter_hi_eff = phase2_aggr_iter_hi;
				      if (phase2_aggr_iter_lo >= 0) aggr_iter_lo_eff = phase2_aggr_iter_lo;
				      if (phase2_aggr_synd_th >= 0) aggr_synd_th_eff = phase2_aggr_synd_th;
				      if (phase2_aggr_strong_synd_th >= 0) aggr_strong_synd_th_eff = phase2_aggr_strong_synd_th;
				    }
				    int w2_boost_to_phase_eff = 3;
				    bool w2_not_pushing_phase_eff = false;
				    if (phase > 0) {
				      w2_boost_to_phase_eff = phase1_w2_boost_to;
				      w2_not_pushing_phase_eff = phase1_w2_not_pushing;
				    }
				    if (phase >= 2) {
				      if (phase2_w2_boost_to >= 0)
				        w2_boost_to_phase_eff = phase2_w2_boost_to;
				      if (phase2_w2_not_pushing_override >= 0)
				        w2_not_pushing_phase_eff = (phase2_w2_not_pushing_override != 0);
				    }
				    const bool relax_post_gate_phase_eff =
				        (restart_relax_post_gate && (phase > 0)) || (phase2_relax_post_gate && (phase >= 2));
				    const bool flip_rand_gate_phase_eff = phase2_flip_rand_gate && (phase >= 2);
				    const int flip_rand_gate_max_w_eff = phase2_flip_rand_gate_max_w;
				    const int flip_rand_gate_gates_eff = phase2_flip_rand_gate_gates;
				    int col_esc_target_col = 0; // selected escape column for this iteration (updated at j==0)
				    int mode_win_rem = 0;
				    int mode_win_cooldown = 0;
				    int hot_col_prev = 0; // hottest column index from the previous iteration (global-only; no per-VN state)

				    while ((iteration < phase_iter_limit) && (finished == 0) && (give_up == 0)) {
				      const bool mode_win_trig =
				          mode_win && (phase > 0) && (iteration >= stall_count_w2_min_iter) &&
				          (iteration >= ldpc_decoder_input.post_iteration) &&
				          (stall_count_w2 >= mode_win_trig_iters);
			      if (mode_win_rem == 0) {
			        if (mode_win_cooldown > 0) {
			          mode_win_cooldown--;
			        } else if (mode_win_trig) {
			          mode_win_rem = std::min(mode_win_len, phase_iter_limit);
			          mode_win_cooldown = mode_win_len; // simple symmetric cooldown by default
			        }
				      }
				      const bool mode_win_active = mode_win && (mode_win_rem > 0);
				      int hot_col_next = hot_col_prev;
				      int hot_col_best_score = -1;
				    for (j = 0; j < h_matrix.cols; j++) {
			      clock_cycles++;
			      if ((iteration == (ldpc_decoder_input.post_iteration + 0)) && (j == 0)) {
	        for (i = 0; i < 256; i++)
	          prng_256.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
	        for (i = 0; i < 512; i++)
	          prng_512.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
	        for (i = 0; i < 512; i++)
	          prng_512.b[i] = (0x1fe0 >> (i % 16)) & 1; // to match verilog
	        if (prng_skip_steps > 0) {
	          for (int s = 0; s < prng_skip_steps; s++) {
	            prng_256 = f_256_bit_lfsr(prng_256);
	            prng_512 = f_512_bit_lfsr(prng_512);
	          }
	        }
		      } else if (iteration >= ldpc_decoder_input.post_iteration) {
		        prng_256 = f_256_bit_lfsr(prng_256);
		        prng_512 = f_512_bit_lfsr(prng_512);
		      }

	      if (col_global_esc && (j == 0) && (h_matrix.cols > 0)) {
	        const bool use_hot_target = col_esc_target_hot && (phase >= col_esc_target_hot_min_phase) &&
	                                    (iteration >= col_esc_target_hot_min_iter);
	        if (use_hot_target) {
	          col_esc_target_col = (h_matrix.cols > 0) ? (hot_col_prev % h_matrix.cols) : 0;
	        } else {
	          // Pick one column per iteration for escape to cap overhead. Derived from existing PRNG state.
	          int rnd = 0;
	          if (h_matrix.bits == 512) {
	            for (int bb = 0; bb < 8; bb++) {
	              const int idx = (iteration * 13 + bb * 37) & 511;
	              if (prng_512.b[idx])
	                rnd |= (1 << bb);
	            }
	          } else {
	            for (int bb = 0; bb < 8; bb++) {
	              const int idx = (iteration * 13 + bb * 37) & 255;
	              if (prng_256.b[idx])
	                rnd |= (1 << bb);
	            }
	          }
	          col_esc_target_col = (h_matrix.cols > 0) ? (rnd % h_matrix.cols) : 0;
	        }
	      }

      syndrome_weight_r[4] = syndrome_weight_r[3];
      syndrome_weight_r[3] = syndrome_weight_r[2];
      syndrome_weight_r[2] = syndrome_weight_r[1];
      syndrome_weight_r[1] = syndrome_weight_r[0];
      syndrome_weight_r[0] = syndrome_weight;
      if (iteration == 0)
        syndrome_weight_delayed = (j <= 3) ? syndrome_weight_r[0] : syndrome_weight_r[4];
      else
        syndrome_weight_delayed = syndrome_weight_r[4];
      const int prev_sw_col = syndrome_weight_r[3];
      if (push_dynamic && (push_mode == 1) && (j == 0)) {
        pushing_iter = (iteration == 0) ? true : (syndrome_weight >= prev_iter_sw);
        prev_iter_sw = syndrome_weight;
      }
      bool pushing = true;
      if (push_dynamic) {
        pushing = (push_mode == 1) ? pushing_iter : (syndrome_weight_delayed >= prev_sw_col);
      }
      w2_boost_sw_delta_allow = true;
      if (w2_boost_sw_delta_gate && (phase > 0) && (iteration >= w2_boost_sw_delta_min_iter)) {
        const int sw_delta_col = syndrome_weight_delayed - prev_sw_col;
        w2_boost_sw_delta_allow = (sw_delta_col >= 0) && (sw_delta_col <= w2_boost_sw_delta_hi);
      }

      if (VERBOSITY > 0)
        printf("### C++ ITERATION %4d, COLUMN %2d, SYNDROME WEIGHT: %4d DELAYED WEIGHT: %4d ###\n", iteration, j,
               syndrome_weight, syndrome_weight_delayed);

      // distribute different thr for post processing stage
      if (iteration >= ldpc_decoder_input.post_iteration) {
        hamming_weight_lt_post_thr = (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_post);
        hamming_weight_lt_circ_thr = (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_qc);
        post_trigger = hamming_weight_lt_circ_thr && ((iteration % 16) < ldpc_decoder_parameters.post_ratio) &&
                       ldpc_decoder_parameters.post_process_en;
        post_trigger2 = hamming_weight_lt_post_thr && ((iteration % 16) >= ldpc_decoder_parameters.post_ratio) &&
                        ldpc_decoder_parameters.post_process_en;

        if (VERBOSITY > 0) {
          printf("### ITERATION %4d POST PROCESSING   BEFORE UPDATE OF COLUMN: %2dWEIGHT: %4d DELAYED WEIGHT: %4d "
                 "POST PROCESS: %d THR_POST: %3d THR_QC: %3d POST_TRIGGER: %d POST_TRIGGER2: %d ###\n",
                 iteration, j, syndrome_weight, syndrome_weight_delayed, post_process,
                 ldpc_decoder_parameters.syndrome_weight_thr_post, ldpc_decoder_parameters.syndrome_weight_thr_qc,
                 post_trigger, post_trigger2);
        }
      } else {
        post_process = 0;
        post_trigger = 0;
        post_trigger2 = 0;
      }

      // be_aggressive is different with rtl and not used
      be_aggressive = (ldpc_decoder_input.soft_bits > 0) &&
                      (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr) && !post_trigger &&
                      !post_trigger2; // different with verilog, add two post trigger judge
      // be_aggressive = (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min <
      // ldpc_decoder_parameters.likelihood_thr) && (iteration < (ldpc_decoder_input.post_iteration << 1));
      // be_aggressive = (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min <
      // ldpc_decoder_parameters.likelihood_thr) && (iteration < (ldpc_decoder_input.post_iteration + 100));
	      // be_aggressive = (ldpc_decoder_input.soft_bits > 0) && (likelihood_levels.min <
	      // ldpc_decoder_parameters.likelihood_thr);
	      const bool rotate_k_eff =
	          rotate_k && (iteration >= ldpc_decoder_input.post_iteration) && (!rotate_k_phase1_only || (phase > 0));
	      int k_start = 0;
	      if (rotate_k_eff) {
	        if (h_matrix.bits == 512) {
	          for (int bb = 0; bb < 9; bb++) {
	            const int idx = (j * 13 + bb * 37) & 511;
	            if (prng_512.b[idx])
	              k_start |= (1 << bb);
	          }
	          k_start &= 511;
	        } else if (h_matrix.bits == 256) {
	          for (int bb = 0; bb < 8; bb++) {
	            const int idx = (j * 13 + bb * 37) & 255;
	            if (prng_256.b[idx])
	              k_start |= (1 << bb);
	          }
	          k_start &= 255;
	        } else {
	          // Fallback for non power-of-two circulants (not expected for IBEX 256/512 cases).
	          for (int bb = 0; bb < 16; bb++) {
	            const int idx = (j * 13 + bb * 37) & 255;
	            k_start = (k_start << 1) | (prng_256.b[idx] ? 1 : 0);
	          }
	          if (h_matrix.bits > 0)
	            k_start = k_start % h_matrix.bits;
	        }
	      }
		      const bool cap_toggles_eff =
		          (max_toggles_per_col > 0) && (!max_toggles_tail_only || (iteration >= ldpc_decoder_input.post_iteration));
		      int toggles_in_col = 0;
		      int col_weight_sum = 0;
		      for (int kk = 0; kk < h_matrix.bits; kk++) {
	        k = rotate_k_eff ? (((h_matrix.bits == 512) ? ((kk + k_start) & 511)
	                                                    : ((h_matrix.bits == 256) ? ((kk + k_start) & 255)
	                                                                              : ((kk + k_start) % h_matrix.bits))))
	                         : kk;
	        look = 0;
	        do_not_use_this_bit = 0;
	        do_not_use_this_bit |= ((h_matrix.extra_bits_of_parity > 0) && (j == (h_matrix.cols - h_matrix.rows)) &&
	                                (k >= h_matrix.extra_bits_of_parity));
        do_not_use_this_bit |= ((h_matrix.extra_bits_of_userdata > 0) && (j == (h_matrix.cols - h_matrix.rows - 1)) &&
                                (k >= h_matrix.extra_bits_of_userdata));
	        if ((VERBOSITY > 0) && do_not_use_this_bit)
	          printf("### DO NOT USE THIS BIT %2d %3d\n", j, k);
	        if (!do_not_use_this_bit) {
	          const int vn_idx = j * h_matrix.bits + k;
	          if (tabu1 && (iteration > 0) && (last_toggle_iter[vn_idx] == static_cast<int16_t>(iteration - 1))) {
	            // One-iteration tabu: skip updating nodes that just toggled last iteration to reduce ping-pong.
	            continue;
	          }
	          weight = 0;
		          for (i = 0; i < h_matrix.rows; i++) {
		            m = (k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits;
		            if (h_matrix.extra_bytes_of_parity == 0) {
	              if (h_matrix.occupied[i][j] && (cn.r[i].b[m] == 1))
	                weight++;
	            } else {
              if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)) && (cn.r[i].b[m] == 1))
                weight++;
              if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) && (cn.r[i].b[m] == 1) && h_matrix.mask[j][k])
                weight++;
              if (h_matrix.fade[i][j] && (cn.r[i].b[m] == 1) && !h_matrix.mask[j][k])
                weight++;
	            }
	          }
		          col_weight_sum += weight;
	          if (iteration == 0) {
	            if (ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_is_error)
	              hard_codeword.errors_at_level_and_weight[ldpc_decoder_input.corrupted_codeword.c[j].b[k].level][weight]++;
	            else
              hard_codeword
                  .correct_at_level_and_weight[ldpc_decoder_input.corrupted_codeword.c[j].b[k].level][weight]++;
          }
	          flipped_prev = vn.c[j].b[k].flipped;
	          int likelihood_prev = vn.c[j].b[k].likelihood;
	          look = (VERBOSITY > 0) && ((j == 33) && (k == 163));
	          const bool soft_unreliable =
	              (ldpc_decoder_input.soft_bits > 0) &&
	              (ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable ||
	               ((ldpc_decoder_input.soft_bits > 1) &&
	                ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2));

          if (iteration >= ldpc_decoder_input.post_iteration) {
            // normal disurbance
            prng_post_process = post_trigger &&
                                (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_qc) &&
                                ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
            // radical disturbance                    
            prng_post_process2 = post_trigger2 &&
                                 (syndrome_weight_delayed < ldpc_decoder_parameters.syndrome_weight_thr_post) &&
                                 ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
          } else {
            prng_post_process = 0;
            prng_post_process2 = 0;
          }
          if ((VERBOSITY > 0) && look) {
            printf("### C++ LOOK ITERATION %4d, SW: %4d, BIT: %2d %3d  ORIGINAL: %x CORRUPTED: %x    WEIGHT: %1d MIN: "
                   "%3d LIKELIHOOD: %3d POST_TRIGGER: %d %d POST: %d %d PRNG[%3d]: %d\n",
                   iteration, syndrome_weight_delayed, j, k, 0,
                   ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard, weight, likelihood_levels.min,
                   vn.c[j].b[k].likelihood, post_trigger, post_trigger2, prng_post_process, prng_post_process2, k,
                   prng_512.b[k]);
            f_print_s_512_bits(prng_512);
          }

          // Adjust weight if aggressive mode
	          // UP-GDBF inspired active-iteration gate for stochastic/escape actions (global only; no per-VN state).
	          int upgdbf_active_iter_eff = upgdbf_active_iter_base;
	          if ((phase > 0) && (upgdbf_active_phase_bonus > 0)) {
	            upgdbf_active_iter_eff = std::max(0, upgdbf_active_iter_eff - upgdbf_active_phase_bonus);
	          }
	          if (upgdbf_active_iter_eff > ldpc_decoder_input.iteration_limit)
	            upgdbf_active_iter_eff = ldpc_decoder_input.iteration_limit;
	          const bool upgdbf_stall_trigger =
	              upgdbf_stall_early && (phase > 0) && (iteration >= stall_count_w2_min_iter) &&
	              (stall_count_w2 >= std::min(stall_w2_esc_iters, ppbf_esc_iters));
	          const bool upgdbf_rand_enable = (iteration >= upgdbf_active_iter_eff) || upgdbf_stall_trigger;

	          bool aggr = ((VN_BITS <= 2) &&
	                       (((iteration >= aggr_iter_hi_eff) && (syndrome_weight < aggr_strong_synd_th_eff)) ||
	                        ((iteration >= aggr_iter_lo_eff) && (syndrome_weight < aggr_synd_th_eff)))) ||
	                      ((ldpc_decoder_input.soft_bits > 0) &&
	                       (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr &&
	                        ((VN_BITS <= 2) || !flipped_prev)));
	          int w = weight;
	          if (aggr) {
	            if (w1_stoch && !flipped_prev && (weight == 1) && pushing &&
	                (likelihood_prev >= (likelihood_levels.flip_thr - 1)) &&
	                (iteration >= ldpc_decoder_input.post_iteration) && upgdbf_rand_enable) {
	              // NGDBF-like noise injection in tail stage: stochastically promote w=1 to w=3 (attack only).
	              const int prng_idx = (h_matrix.bits == 512) ? ((k + 101) & 511) : ((k + 101) & 255);
	              const bool rand_boost = (h_matrix.bits == 512) ? prng_512.b[prng_idx] : prng_256.b[prng_idx];
	              if (rand_boost)
	                w = 3;
	            }
	            if (weight == 2) {
	              // Only boost attack; keep retract conservative to reduce oscillation.
	              if (!flipped_prev) {
	                if (!soft_guard || soft_unreliable) {
	                  const bool w2_soft_only_ok = !w2_boost_soft_only || soft_unreliable;
	                  if (w2_soft_only_ok) {
			                    const bool w2_strong_unflipped =
			                        (VN_BITS <= 2) && (likelihood_prev < (likelihood_levels.flip_thr - 1));
				                    const int w2_boost_to_eff = w2_boost_to_phase_eff;
			                    const bool stall_w2_esc_gate =
			                        mode_win ? mode_win_active : (stall_count_w2 >= stall_w2_esc_iters);
			                    const bool stall_w2_esc_active =
			                        stall_w2_esc && (phase > 0) && stall_w2_esc_gate &&
			                        (iteration >= stall_w2_esc_min_iter);
		                    const bool ppbf_esc_gate =
		                        mode_win ? mode_win_active : (stall_count_w2 >= ppbf_esc_iters);
		                    const bool ppbf_esc_active =
		                        ppbf_esc && upgdbf_rand_enable && (phase > 0) && ppbf_esc_gate &&
		                        (iteration >= std::max(ppbf_esc_min_iter, ldpc_decoder_input.post_iteration));
			                    bool w2_cand_allow = true;
			                    if (w2_cand_strong && push_dynamic && w2_strong_unflipped && !w2_cand_mem.empty()) {
			                      if (!pushing || (iteration < ldpc_decoder_input.post_iteration)) {
			                        w2_cand_mem[vn_idx] = 0;
	                      } else if (!w2_cand_mem[vn_idx]) {
	                        w2_cand_mem[vn_idx] = 1; // arm; require a second hit to boost
	                        w2_cand_allow = false;
	                      }
	                    } else if (w2_cand_strong && !w2_cand_mem.empty()) {
	                      // Clear stale cand as soon as the VN becomes weak/flipped/weight!=2.
	                      w2_cand_mem[vn_idx] = 0;
	                    }
	                    if (!w2_cand_allow) {
	                      // No boost on the first strong w=2 hit while stalled.
	                    } else
	                    if (w2_stoch && push_dynamic) {
	                      if (!upgdbf_rand_enable) {
	                        // Before active-iter enables randomness, keep w=2 unchanged (no stochastic boost).
	                      } else if (iteration >= ldpc_decoder_input.post_iteration) {
	                        if (pushing) {
	                          const int prng_idx1 = (h_matrix.bits == 512) ? ((k + 37) & 511) : ((k + 37) & 255);
	                          bool rand_boost =
	                              (h_matrix.bits == 512) ? prng_512.b[prng_idx1] : prng_256.b[prng_idx1];
	                          if (w2_stoch_xor) {
	                            const int prng_idx2 = (h_matrix.bits == 512) ? ((k + 173) & 511) : ((k + 173) & 255);
	                            const bool rand_boost2 =
	                                (h_matrix.bits == 512) ? prng_512.b[prng_idx2] : prng_256.b[prng_idx2];
	                            rand_boost ^= rand_boost2;
	                          }
	                          if (w2_boost_weak_only && w2_strong_unflipped) {
	                            // Strong unflipped VN: reduce w2-boost probability to avoid slowly pushing correct bits.
	                            const int prng_idx3 = (h_matrix.bits == 512) ? ((k + 211) & 511) : ((k + 211) & 255);
	                            const bool rand_gate3 =
	                                (h_matrix.bits == 512) ? prng_512.b[prng_idx3] : prng_256.b[prng_idx3];
	                            rand_boost &= rand_gate3; // ~1/4 boost probability (assuming independent)
	                          }
		                          if (w2_tail_guard && w2_strong_unflipped && !soft_unreliable) {
		                            // Tail-stage safety: when only a few checks remain unsatisfied, avoid w=2 boosting
		                            // strong+reliable VNs (a common source of rare mis-flips / error-floor packets).
		                            if (hamming_weight_lt_post_thr) {
		                              rand_boost = false;
		                            } else if (hamming_weight_lt_circ_thr) {
		                              const int prng_idx4 = (h_matrix.bits == 512) ? ((k + 233) & 511) : ((k + 233) & 255);
		                              const bool rand_gate4 =
		                                  (h_matrix.bits == 512) ? prng_512.b[prng_idx4] : prng_256.b[prng_idx4];
		                              rand_boost &= rand_gate4; // further reduce probability in tail
		                            }
		                          }
		                          if (stall_w2_esc_active) {
		                            // When stuck in a retry phase, make escaping more decisive for "bad" candidates.
		                            const bool cand_bad = (!w2_strong_unflipped) || soft_unreliable;
		                            if (cand_bad && (!w2_boost_weak_only || !w2_strong_unflipped))
		                              rand_boost = true;
		                          }
		                          if (rand_boost && w2_boost_sw_delta_allow)
		                            w = w2_boost_to_eff;
		                        } else {
		                          if (!w2_boost_only_when_pushing) {
		                            if (w2_boost_sw_delta_allow && (!w2_boost_weak_only || !w2_strong_unflipped))
		                              w = w2_boost_to_eff;
		                          } else if (stall_w2_esc_active) {
		                            // If a phase has stalled, allow a small amount of w=2 boost even when not pushing.
		                            const bool cand_bad = (!w2_strong_unflipped) || soft_unreliable;
		                            if (cand_bad && (!w2_boost_weak_only || !w2_strong_unflipped)) {
		                              const int prng_idx6 = (h_matrix.bits == 512) ? ((k + 101) & 511) : ((k + 101) & 255);
		                              const bool r6 = (h_matrix.bits == 512) ? prng_512.b[prng_idx6] : prng_256.b[prng_idx6];
		                              if (r6 && w2_boost_sw_delta_allow)
		                                w = w2_boost_to_eff;
		                            }
			                          } else if (w2_not_pushing_phase_eff) {
		                            // Retry-phase escape hatch: even if "only-when-pushing" is enabled, allow a very small
		                            // probability of w=2 boost when the trend is not pushing, to escape some trap sets.
		                            const bool cand_bad = (!w2_strong_unflipped) || soft_unreliable;
		                            if (cand_bad && (!w2_boost_weak_only || !w2_strong_unflipped)) {
		                              const int prng_idx4 = (h_matrix.bits == 512) ? ((k + 211) & 511) : ((k + 211) & 255);
		                              const int prng_idx5 = (h_matrix.bits == 512) ? ((k + 233) & 511) : ((k + 233) & 255);
		                              const bool r4 = (h_matrix.bits == 512) ? prng_512.b[prng_idx4] : prng_256.b[prng_idx4];
		                              const bool r5 = (h_matrix.bits == 512) ? prng_512.b[prng_idx5] : prng_256.b[prng_idx5];
		                              // ~1/4 (two independent gates) on top of the base rand_boost in the pushing path.
		                              if (r4 && r5 && w2_boost_sw_delta_allow)
		                                w = w2_boost_to_eff;
		                            }
		                          } else if (ppbf_esc_active) {
		                            // PPBF-like escape hatch: boost w=2 with probability p(E) based on a tiny energy proxy.
		                            const bool cand_bad = (!w2_strong_unflipped) || soft_unreliable;
		                            if (cand_bad && (!w2_boost_weak_only || !w2_strong_unflipped)) {
		                              int E = weight; // here weight==2
		                              if (soft_unreliable)
		                                E++;
		                              if (!pushing)
		                                E++;
		                              if (ngdbf_noise && (weight >= ngdbf_noise_min_w)) {
		                                const bool noise_stall_ok =
		                                    !ngdbf_noise_only_stall ||
		                                    ((iteration >= stall_count_w2_min_iter) &&
		                                     (stall_count_w2 >= ngdbf_noise_stall_iters));
		                                if (noise_stall_ok) {
		                                  bool bump = true;
		                                  for (int gg = 0; gg < ngdbf_noise_gates; gg++) {
		                                    const int idx = (h_matrix.bits == 512)
		                                                        ? ((k + 19 + gg * 37) & 511)
		                                                        : ((k + 19 + gg * 37) & 255);
		                                    const bool rb = (h_matrix.bits == 512) ? prng_512.b[idx] : prng_256.b[idx];
		                                    bump &= rb;
		                                  }
		                                  if (bump)
		                                    E++; // discrete +1 noise bump (NGDBF-lite)
		                                }
		                              }
		                              const int prng_idx7 = (h_matrix.bits == 512) ? ((k + 79) & 511) : ((k + 79) & 255);
		                              const int prng_idx8 = (h_matrix.bits == 512) ? ((k + 157) & 511) : ((k + 157) & 255);
		                              const int prng_idx9 = (h_matrix.bits == 512) ? ((k + 211) & 511) : ((k + 211) & 255);
		                              const bool r7 = (h_matrix.bits == 512) ? prng_512.b[prng_idx7] : prng_256.b[prng_idx7];
		                              const bool r8 = (h_matrix.bits == 512) ? prng_512.b[prng_idx8] : prng_256.b[prng_idx8];
		                              const bool r9 = (h_matrix.bits == 512) ? prng_512.b[prng_idx9] : prng_256.b[prng_idx9];
		                              const bool ppbf_hot =
		                                  (mode_win_active && mode_win_ppbf_boost && !hamming_weight_lt_post_thr) ||
		                                  (ppbf_anneal && !hamming_weight_lt_post_thr &&
		                                   (iteration < (ppbf_esc_min_iter + ppbf_anneal_iters)));
		                              const bool ppbf_cold = ppbf_anneal && !ppbf_hot;
		                              bool rand_boost = false;
		                              if (ppbf_hot) {
		                                if (E >= 4) {
		                                  rand_boost = r7; // ~1/2
		                                } else if (E == 3) {
		                                  rand_boost = r7; // ~1/2 (hotter)
		                                } else if (E == 2) {
		                                  rand_boost = r7 && r8 && r9; // ~1/8 (very small widening)
		                                }
		                              } else if (ppbf_cold) {
		                                if (E >= 4) {
		                                  rand_boost = r7 && r8; // ~1/4
		                                } else if (E == 3) {
		                                  rand_boost = r7 && r8 && r9; // ~1/8
		                                }
		                              } else {
		                                if (E >= 4) {
		                                  rand_boost = r7; // ~1/2
		                                } else if (E == 3) {
		                                  rand_boost = r7 && r8; // ~1/4
		                                }
		                              }
		                              if (rand_boost && w2_boost_sw_delta_allow)
		                                w = w2_boost_to_eff;
		                            }
		                          }
		                        }
		                      } else {
		                        // Before post_iteration, avoid stochastic/early w=2 boost to reduce chaotic mass-flips.
	                      }
	                    } else {
	                      if (w2_boost_sw_delta_allow && (!w2_boost_weak_only || !w2_strong_unflipped))
	                        w = w2_boost_to_eff;
	                    }
	                  }
	                }
	              }
	            } else if (!flipped_prev && (weight == 3)) {
	              w = 5; // 3->5 (attack)
	            } else if (!flipped_prev && (weight == 4)) {
	              w = 7; // 4->7 (attack)
            }
          }
          if ((VERBOSITY > 0) && look)
				            printf("C++ LOOK   ITERATION: %4d BEFORE LIKELIHOOD UPDATE SW: %4d CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x "
				                   "CORRUPTED: %x FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: %3d\n",
				                   iteration, syndrome_weight, j, k, 0, ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
				                   vn.c[j].b[k].flipped, weight, likelihood_levels.min, vn.c[j].b[k].likelihood);
	          const bool post_gate_mode = (mode_2bit >= 2) ? (!aggr || relax_post_gate_phase_eff) : true;
          const bool post_gate_push = post_only_when_pushing ? pushing : true;
          const bool post_gate_soft = soft_guard ? soft_unreliable : true;
          const bool post_gate = post_gate_mode && post_gate_push && post_gate_soft;
          const bool prng_post_process_eff = prng_post_process && post_gate;
          const bool prng_post_process2_eff = prng_post_process2 && post_gate;
	          vn.c[j].b[k].likelihood =
			              f_update_vn_post(vn.c[j].b[k].likelihood, w, likelihood_levels.min, likelihood_levels.max,
			                               prng_post_process_eff, prng_post_process2_eff, aggr, likelihood_levels.flip_thr, pushing);

          if ((VERBOSITY > 0) && look)
            printf("C++ LOOK   ITERATION: %4d AFTER LIKELIHOOD UPDATE SW: %4d CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x "
                   "CORRUPTED: %x FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: %3d\n",
                   iteration, syndrome_weight, j, k, 0, ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                   vn.c[j].b[k].flipped, weight, likelihood_levels.min, vn.c[j].b[k].likelihood);
	          vn.c[j].b[k].flipped = (vn.c[j].b[k].likelihood >= likelihood_levels.flip_thr);
	          if (toggle_strong && (VN_BITS <= 2) && (flipped_prev != vn.c[j].b[k].flipped) && (w >= 5)) {
	            // Implicit 1-step hysteresis without extra per-VN state:
	            // If a VN toggles due to strong evidence (w>=5), snap to a strong state to reduce oscillation.
	            vn.c[j].b[k].likelihood =
	                vn.c[j].b[k].flipped ? likelihood_levels.max : likelihood_levels.min;
	          }
	          if (tabu_rev && (iteration > 0) &&
	              (last_toggle_iter[vn_idx] == static_cast<int16_t>(iteration - 1)) &&
	              (flipped_prev != vn.c[j].b[k].flipped)) {
	            // TRGDBF-style 1-iteration tabu: allow likelihood update, but forbid immediate toggle back.
	            vn.c[j].b[k].likelihood = flipped_prev ? likelihood_levels.flip_thr : (likelihood_levels.flip_thr - 1);
	            vn.c[j].b[k].flipped = flipped_prev;
	          }

	          if ((VERBOSITY > 0) && look)
	            printf("C++ ITERATION: %4d AFTER FLIP CHECK  SW: %4d CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x CORRUPTED: %x "
	                   "FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: %3d THR: %3d\n",
	                   iteration, syndrome_weight, j, k, 0, ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
	                   vn.c[j].b[k].flipped, weight, likelihood_levels.min, vn.c[j].b[k].likelihood,
	                   likelihood_levels.flip_thr);

			          if (cap_toggles_eff && (flipped_prev != vn.c[j].b[k].flipped) &&
			              (toggles_in_col >= max_toggles_per_col)) {
			            // Tail-stage safety: suppress excessive toggles within a column and clamp to the boundary.
			            vn.c[j].b[k].likelihood = flipped_prev ? likelihood_levels.flip_thr : (likelihood_levels.flip_thr - 1);
			            vn.c[j].b[k].flipped = flipped_prev;
			          }

			          if (flip_rand_gate_phase_eff && hamming_weight_lt_circ_thr &&
			              (iteration >= ldpc_decoder_input.post_iteration) && (flipped_prev != vn.c[j].b[k].flipped) &&
			              (weight <= flip_rand_gate_max_w_eff)) {
			            // FM-PGDBF style probability gate: randomly suppress some "weak" flip events to break
			            // deterministic time cycles. Phase2-only to isolate side-effects; no per-VN state.
			            bool allow = false;
			            for (int gg = 0; gg < flip_rand_gate_gates_eff; gg++) {
			              const int idx = (h_matrix.bits == 512) ? ((k + 37 + gg * 73 + j * 11) & 511)
			                                                    : ((k + 37 + gg * 73 + j * 11) & 255);
			              const bool r = (h_matrix.bits == 512) ? prng_512.b[idx] : prng_256.b[idx];
			              allow |= r;
			            }
			            if (!allow) {
			              vn.c[j].b[k].likelihood = flipped_prev ? likelihood_levels.flip_thr : (likelihood_levels.flip_thr - 1);
			              vn.c[j].b[k].flipped = flipped_prev;
			            }
			          }

			          // Increment update syndrome if flipped state changed
			          if (flipped_prev != vn.c[j].b[k].flipped) {
			            toggles_in_col++;
			            if (tabu1 || tabu_rev)
		              last_toggle_iter[vn_idx] = static_cast<int16_t>(iteration);
		            for (i = 0; i < h_matrix.rows; i++) {
		              m = (k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits;
	              if (h_matrix.extra_bytes_of_parity == 0) {
	                if (h_matrix.occupied[i][j])
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
              } else {
                if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)))
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
                if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) && h_matrix.mask[j][k])
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
                if (h_matrix.fade[i][j] && !h_matrix.mask[j][k])
                  cn.r[i].b[m] = 1 - cn.r[i].b[m];
              }
            }
            if ((VERBOSITY > 0) && vn.c[j].b[k].flipped)
              printf("### C++ ITERATION %4d   FLIPPED BIT: %2d %3d ###\n", iteration, j, k);
            if ((VERBOSITY > 0) && !vn.c[j].b[k].flipped)
              printf("### C++ ITERATION %4d   UNFLIPPED BIT: %2d %3d ###\n", iteration, j, k);

            if (VERBOSITY > 0) {
              printf("### DECODER C++: AFTER TOGGLE, CODEWORD HAS CHECK NODES:\n");
              f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
            }
          }
	        }
	      }

	      if (col_esc_target_hot && (phase >= col_esc_target_hot_min_phase)) {
	        if (col_weight_sum > hot_col_best_score) {
	          hot_col_best_score = col_weight_sum;
	          hot_col_next = j;
	        }
	      }

	      // Column-level global escape + optional energy backtracking (no per-VN state; per-column tiny buffers only).
	      if (col_global_esc) {
	        int upgdbf_active_iter_eff2 = upgdbf_active_iter_base;
        if ((phase > 0) && (upgdbf_active_phase_bonus > 0)) {
          upgdbf_active_iter_eff2 = std::max(0, upgdbf_active_iter_eff2 - upgdbf_active_phase_bonus);
        }
        if (upgdbf_active_iter_eff2 > ldpc_decoder_input.iteration_limit)
          upgdbf_active_iter_eff2 = ldpc_decoder_input.iteration_limit;
        const bool upgdbf_stall_trigger2 =
            upgdbf_stall_early && (phase > 0) && (iteration >= stall_count_w2_min_iter) &&
            (stall_count_w2 >= std::min(stall_w2_esc_iters, ppbf_esc_iters));
        const bool upgdbf_rand_enable2 = (iteration >= upgdbf_active_iter_eff2) || upgdbf_stall_trigger2;

	        const bool col_esc_active =
	            col_global_esc && upgdbf_rand_enable2 && (phase >= col_global_esc_min_phase) && (j == col_esc_target_col) &&
	            (stall_count_w2 >= col_global_esc_iters) &&
	            hamming_weight_lt_circ_thr &&
	            (iteration >= std::max(col_global_esc_min_iter, ldpc_decoder_input.post_iteration));
        if (col_esc_active) {
          const int sw_before_escape = f_check_node_weight(h_matrix, cn);

          int max_w_in_col = -1;
          for (int kk2 = 0; kk2 < h_matrix.bits; kk2++) {
            const int k2 = kk2;
            bool dnu = false;
            dnu |= ((h_matrix.extra_bits_of_parity > 0) && (j == (h_matrix.cols - h_matrix.rows)) &&
                    (k2 >= h_matrix.extra_bits_of_parity));
            dnu |= ((h_matrix.extra_bits_of_userdata > 0) && (j == (h_matrix.cols - h_matrix.rows - 1)) &&
                    (k2 >= h_matrix.extra_bits_of_userdata));
            if (dnu)
              continue;

            int w2tmp = 0;
            for (int ii2 = 0; ii2 < h_matrix.rows; ii2++) {
              const int mm2 = (k2 + h_matrix.bits - h_matrix.element[ii2][j]) % h_matrix.bits;
              if (h_matrix.extra_bytes_of_parity == 0) {
                if (h_matrix.occupied[ii2][j] && (cn.r[ii2].b[mm2] == 1))
                  w2tmp++;
              } else {
                if (h_matrix.occupied[ii2][j] && (ii2 < (h_matrix.rows - 1)) && (cn.r[ii2].b[mm2] == 1))
                  w2tmp++;
                if (h_matrix.occupied[ii2][j] && (ii2 == (h_matrix.rows - 1)) && (cn.r[ii2].b[mm2] == 1) &&
                    h_matrix.mask[j][k2])
                  w2tmp++;
                if (h_matrix.fade[ii2][j] && (cn.r[ii2].b[mm2] == 1) && !h_matrix.mask[j][k2])
                  w2tmp++;
              }
            }
            if (w2tmp > max_w_in_col)
              max_w_in_col = w2tmp;
          }

          if (max_w_in_col >= 2) {
            const bool col_mode_win_hot = mode_win_active && mode_win_col_esc_boost;
            const int col_esc_max_toggles_eff =
                col_mode_win_hot ? std::min(8, col_global_esc_max_toggles + 1) : col_global_esc_max_toggles;
            int esc_k[8] = {};
            int16_t esc_like[8] = {};
            uint8_t esc_flip[8] = {};
            int esc_cnt = 0;

            for (int kk2 = 0; (kk2 < h_matrix.bits) && (esc_cnt < col_esc_max_toggles_eff); kk2++) {
              const int k2 = kk2;
              bool dnu = false;
              dnu |= ((h_matrix.extra_bits_of_parity > 0) && (j == (h_matrix.cols - h_matrix.rows)) &&
                      (k2 >= h_matrix.extra_bits_of_parity));
              dnu |= ((h_matrix.extra_bits_of_userdata > 0) && (j == (h_matrix.cols - h_matrix.rows - 1)) &&
                      (k2 >= h_matrix.extra_bits_of_userdata));
              if (dnu)
                continue;

              // Attack-only escape: only consider currently-unflipped bits (avoid oscillatory unflip).
              if (vn.c[j].b[k2].flipped)
                continue;
              const int weak_thr = std::max(0, likelihood_levels.flip_thr - 1);
              if (vn.c[j].b[k2].likelihood < weak_thr)
                continue;

              const bool soft0 = ldpc_decoder_input.corrupted_codeword.c[j].b[k2].bit_questionable;
              const bool soft1 = ldpc_decoder_input.corrupted_codeword.c[j].b[k2].bit_questionable2;
              const bool soft_unreliable2 =
                  soft0 || ((ldpc_decoder_input.soft_bits > 1) && soft1);
              if (col_global_esc_soft_only && !soft_unreliable2)
                continue;

              int w2tmp = 0;
              for (int ii2 = 0; ii2 < h_matrix.rows; ii2++) {
                const int mm2 = (k2 + h_matrix.bits - h_matrix.element[ii2][j]) % h_matrix.bits;
                if (h_matrix.extra_bytes_of_parity == 0) {
                  if (h_matrix.occupied[ii2][j] && (cn.r[ii2].b[mm2] == 1))
                    w2tmp++;
                } else {
                  if (h_matrix.occupied[ii2][j] && (ii2 < (h_matrix.rows - 1)) && (cn.r[ii2].b[mm2] == 1))
                    w2tmp++;
                  if (h_matrix.occupied[ii2][j] && (ii2 == (h_matrix.rows - 1)) && (cn.r[ii2].b[mm2] == 1) &&
                      h_matrix.mask[j][k2])
                    w2tmp++;
                  if (h_matrix.fade[ii2][j] && (cn.r[ii2].b[mm2] == 1) && !h_matrix.mask[j][k2])
                    w2tmp++;
                }
              }
              if (w2tmp != max_w_in_col)
                continue;

              if (cap_toggles_eff && (toggles_in_col >= max_toggles_per_col))
                break;

              // Energy-based probabilistic acceptance: p~1/2 for w>=3, p~1/4 for w==2.
              const int prng_idx7 = (h_matrix.bits == 512) ? ((k2 + 79) & 511) : ((k2 + 79) & 255);
              const int prng_idx8 = (h_matrix.bits == 512) ? ((k2 + 157) & 511) : ((k2 + 157) & 255);
              const int prng_idx9 = (h_matrix.bits == 512) ? ((k2 + 211) & 511) : ((k2 + 211) & 255);
              const bool r7 = (h_matrix.bits == 512) ? prng_512.b[prng_idx7] : prng_256.b[prng_idx7];
              const bool r8 = (h_matrix.bits == 512) ? prng_512.b[prng_idx8] : prng_256.b[prng_idx8];
              const bool r9 = (h_matrix.bits == 512) ? prng_512.b[prng_idx9] : prng_256.b[prng_idx9];
              bool accept = false;
              if (max_w_in_col >= 3) {
                accept = r7;
              } else if (max_w_in_col == 2) {
                accept = col_mode_win_hot ? (r7 && r8) : (r7 && r8 && r9); // ~1/4 vs ~1/8
              }
              if (!accept)
                continue;

              esc_k[esc_cnt] = k2;
              esc_like[esc_cnt] = static_cast<int16_t>(vn.c[j].b[k2].likelihood);
              esc_flip[esc_cnt] = static_cast<uint8_t>(vn.c[j].b[k2].flipped ? 1 : 0);
              esc_cnt++;

              // Force a single toggle by snapping likelihood across the threshold (pure 2bit, no extra state).
              const bool flipped_prev2 = vn.c[j].b[k2].flipped;
              vn.c[j].b[k2].likelihood = flipped_prev2 ? likelihood_levels.min : likelihood_levels.max;
              vn.c[j].b[k2].flipped = (vn.c[j].b[k2].likelihood >= likelihood_levels.flip_thr);

              if (flipped_prev2 != vn.c[j].b[k2].flipped) {
                toggles_in_col++;
                for (int ii2 = 0; ii2 < h_matrix.rows; ii2++) {
                  const int mm2 = (k2 + h_matrix.bits - h_matrix.element[ii2][j]) % h_matrix.bits;
                  if (h_matrix.extra_bytes_of_parity == 0) {
                    if (h_matrix.occupied[ii2][j])
                      cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                  } else {
                    if (h_matrix.occupied[ii2][j] && (ii2 < (h_matrix.rows - 1)))
                      cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                    if (h_matrix.occupied[ii2][j] && (ii2 == (h_matrix.rows - 1)) && h_matrix.mask[j][k2])
                      cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                    if (h_matrix.fade[ii2][j] && !h_matrix.mask[j][k2])
                      cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                  }
                }
              }
            }

            if (col_backtrack && (esc_cnt > 0)) {
              const int sw_after_escape = f_check_node_weight(h_matrix, cn);
              const bool need_backtrack =
                  col_backtrack_require_improve ? (sw_after_escape >= sw_before_escape)
                                                : (sw_after_escape > (sw_before_escape + col_backtrack_slack));
              if (need_backtrack) {
                // Backtrack: restore VN likelihoods + undo syndrome toggles for this column.
                for (int r = 0; r < esc_cnt; r++) {
                  const int k2 = esc_k[r];
                  const int16_t like_prev2 = esc_like[r];
                  const bool flip_prev2 = (esc_flip[r] != 0);
                  const bool flip_cur2 = vn.c[j].b[k2].flipped;
                  vn.c[j].b[k2].likelihood = like_prev2;
                  vn.c[j].b[k2].flipped = flip_prev2;
                  if (flip_cur2 != vn.c[j].b[k2].flipped) {
                    if (toggles_in_col > 0)
                      toggles_in_col--;
                    for (int ii2 = 0; ii2 < h_matrix.rows; ii2++) {
                      const int mm2 = (k2 + h_matrix.bits - h_matrix.element[ii2][j]) % h_matrix.bits;
                      if (h_matrix.extra_bytes_of_parity == 0) {
                        if (h_matrix.occupied[ii2][j])
                          cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                      } else {
                        if (h_matrix.occupied[ii2][j] && (ii2 < (h_matrix.rows - 1)))
                          cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                        if (h_matrix.occupied[ii2][j] && (ii2 == (h_matrix.rows - 1)) && h_matrix.mask[j][k2])
                          cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                        if (h_matrix.fade[ii2][j] && !h_matrix.mask[j][k2])
                          cn.r[ii2].b[mm2] = 1 - cn.r[ii2].b[mm2];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      // Shift check nodes for this column (for debug / alignment).
      for (i = 0; i < h_matrix.rows; i++) {
        if (h_matrix.occupied[i][j] || h_matrix.fade[i][j])
          for (k = 0; k < h_matrix.bits; k++)
            cn_shifted.r[i].b[k] = cn.r[i].b[(k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits];
      }
      if (VERBOSITY > 0) {
        printf("### DECODER C++: DECODING ITERATION %4d AT COLUMN %2d, CODEWORD HAS CHECK NODES:\n", iteration, j);
        f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
        printf("### DECODER C++: DECODING ITERATION %4d AT COLUMN %2d, CODEWORD HAS CHECK NODES SHIFTED:\n", iteration,
               j);
        f_print_check_nodes_shifted(cn, h_matrix, j);
      }

      syndrome_weight = f_check_node_weight(h_matrix, cn);
      if (syndrome_weight == 0) {
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

			    if (col_esc_target_hot && (phase >= col_esc_target_hot_min_phase)) {
			      hot_col_prev = hot_col_next;
			    }

				    if (smooth_fail && (finished == 0) && (iteration >= smooth_start_iter)) {
				      if (smooth_sum.empty())
				        smooth_sum.assign(h_matrix.cols * h_matrix.bits, 0);
			      for (int jj = 0; jj < h_matrix.cols; jj++) {
			        for (int kk = 0; kk < h_matrix.bits; kk++) {
			          const int vn_idx = jj * h_matrix.bits + kk;
			          const bool decoded_bit = ldpc_decoder_input.corrupted_codeword.c[jj].b[kk].bit_hard ^ vn.c[jj].b[kk].flipped;
			          smooth_sum[vn_idx] += decoded_bit ? -1 : 1; // 0->+1, 1->-1
			        }
			      }
			    }

					    if (mode_win_rem > 0)
					      mode_win_rem--;
					    clock_cycles++;
					    iteration++;
						    if ((finished == 0) && (give_up == 0) && (iteration >= stall_min_iter)) {
						      if (syndrome_weight < best_sw_phase) {
						        best_sw_phase = syndrome_weight;
						        stall_count = 0;
						      } else {
						        stall_count++;
						      }
						    }
							    if ((stall_w2_esc || ppbf_esc || col_global_esc || mode_win || ngdbf_noise) && (finished == 0) &&
							        (give_up == 0)) {
							      if (syndrome_weight < best_sw_phase_w2) {
							        best_sw_phase_w2 = syndrome_weight;
							        if (iteration >= stall_count_w2_min_iter)
							          stall_count_w2 = 0;
						      } else if (iteration >= stall_count_w2_min_iter) {
						        stall_count_w2++;
						      }
						    }
						    if (restart_on_stall && (finished == 0) && (give_up == 0) && ((phase + 1) < restart_phases) &&
						        (iteration >= stall_min_iter) && (stall_count >= stall_iters)) {
						      break; // early restart next phase
						    }
					  }

				    if (smooth_fail && (finished == 0) && !smooth_sum.empty()) {
				      // If the last phase fails, try a hard-decision majority vote over the tail window.
				      s_hard_codeword hard_codeword_smoothed = {};
			      for (int jj = 0; jj < h_matrix.cols; jj++) {
			        for (int kk = 0; kk < h_matrix.bits; kk++) {
			          const int vn_idx = jj * h_matrix.bits + kk;
			          const int sum = smooth_sum[vn_idx];
			          bool decoded_bit;
			          if (sum > 0)
			            decoded_bit = 0;
			          else if (sum < 0)
			            decoded_bit = 1;
			          else
			            decoded_bit = ldpc_decoder_input.corrupted_codeword.c[jj].b[kk].bit_hard ^ vn.c[jj].b[kk].flipped;
			          hard_codeword_smoothed.c[jj].b[kk] = decoded_bit;
			        }
			      }
			      s_check_nodes cn_smoothed = f_check_nodes(h_matrix, hard_codeword_smoothed);
			      const int sw_smoothed = f_check_node_weight(h_matrix, cn_smoothed);
			      if (sw_smoothed == 0) {
			        cn = cn_smoothed;
			        syndrome_weight = 0;
			        finished = 1;
			        for (int jj = 0; jj < h_matrix.cols; jj++) {
			          for (int kk = 0; kk < h_matrix.bits; kk++) {
			            vn.c[jj].b[kk].flipped =
			                ldpc_decoder_input.corrupted_codeword.c[jj].b[kk].bit_hard ^ hard_codeword_smoothed.c[jj].b[kk];
			          }
			        }
			      }
			    }
			  }

		  if (VERBOSITY > 0) {
		    printf("### AFTER %4d ITERATIONS, CODEWORD HAS CHECK NODES:\n", iteration);
	    f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
	  }

  ldpc_decoder_output.iterations = iteration;
  ldpc_decoder_output.clock_cycles = clock_cycles;
  ldpc_decoder_output.syndrome_weight_after = f_check_node_weight(h_matrix, cn);
  ldpc_decoder_output.failure = (ldpc_decoder_output.syndrome_weight_after != 0);
  ldpc_decoder_output.errors_in_userdata = 0;
  ldpc_decoder_output.errors_in_codeword = 0;
  for (j = 0; j < h_matrix.cols; j++) {
    for (k = 0; k < h_matrix.bits; k++) {
      ldpc_decoder_output.corrected_codeword.c[j].b[k] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard ^ vn.c[j].b[k].flipped;
      dec_do_blk[j * h_matrix.bits + k] = ldpc_decoder_output.corrected_codeword.c[j].b[k];
      if (vn.c[j].b[k].flipped)
        ldpc_decoder_output.errors_in_codeword++;
      if (vn.c[j].b[k].flipped && (j < (h_matrix.cols - h_matrix.rows)))
        ldpc_decoder_output.errors_in_userdata++;
      if ((VERBOSITY > 0) && (vn.c[j].b[k].flipped))
        printf("### C FLIPPED BIT %2d %3d ### %4d ERRORS IN CODEWORD, %4d ERRORS IN USERDATA\n", j, k,
               ldpc_decoder_output.errors_in_codeword, ldpc_decoder_output.errors_in_userdata);
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
    f_print_hard_codeword(ldpc_decoder_output.corrected_codeword, h_matrix.cols, h_matrix.bits);
  }
  cw_fail = ldpc_decoder_output.failure;
}
