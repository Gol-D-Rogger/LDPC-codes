#include <algorithm>
#include <cstdint>
#include <cstring>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>

#include "finite_lib.h"
#include "ldpc_codec.h"
#include "mod2convert.h"
#include "mod2dense.h"
#include "mod2sparse.h"
#include "vec_op.h"

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
          } else {
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

// LDPC code configuration
void ldpc_packet::ldpc_config(int m, int n, int sc) {

  // read parity check matrix
  if (sc == 256) {
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

    hm_m = m*sc;
    hm_n = n*sc;
    hm_k = hm_n - hm_m;

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

    pad_len = h_matrix.unused_bytes_of_userdata*8;

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

    int use_customized_matrix = 1; // 0 not change 1: M3's matrix
    bool print_customize_matrix = true;
    bool customize_config1 = (h_matrix.rows == 8) && (h_matrix.cols == 73);
    bool customize_config2 = (h_matrix.rows == 10) && (h_matrix.cols == 77);
    if ((customize_config1 || customize_config2) && (use_customized_matrix == 1)) {
      FILE *fp_occu, *fp_fade, *fp_elem;
      int file_n; // actual number of columns in the matrix file
      if (customize_config1)
      {
        // printf("read customized matrix for 8x73\n");
        // fp_occu = fopen("./src/matrice/8x73/occupied_matrix/LDPC_8x73ex512_w4_dense5_occupied_1_0.txt", "r");
        // fp_fade = fopen("./src/matrice/8x73/fade_matrix/LDPC_8x73ex512_w4_dense5_fade_1_0.txt", "r");
        // fp_elem = fopen("./src/matrice/8x73/matrix/LDPC_8x73ex512_w4_dense5_QC_H_1_0.txt", "r");
        file_n = 75;
        printf("read customized matrix for 8x%d, trim to 8x%d\n", file_n, h_matrix.cols);
        fp_occu = fopen("./src/matrice/8x75/occupied_matrix/LDPC_8x75ex512_w4_dense5_occupied_1.txt", "r");
        fp_fade = fopen("./src/matrice/8x75/fade_matrix/LDPC_8x75ex512_w4_dense5_fade_1.txt", "r");
        fp_elem = fopen("./src/matrice/8x75/matrix/LDPC_8x75ex512_w4_dense5_QC_H_1.txt", "r");
      }
      else
      {
        file_n = 77;
        printf("read customized matrix for 10x%d\n", file_n);
        fp_occu = fopen("./src/matrice/10x77/occupied_matrix/LDPC_10x77ex512_w4_dense5_occupied_1.txt", "r");
        fp_fade = fopen("./src/matrice/10x77/fade_matrix/LDPC_10x77ex512_w4_dense5_fade_1.txt", "r");
        fp_elem = fopen("./src/matrice/10x77/matrix/LDPC_10x77ex512_w4_dense5_QC_H_1.txt", "r");
      }
      if (file_n == h_matrix.cols) {
        // No trimming needed — direct read (original behavior)
        for (int i = 0; i < h_matrix.rows; i++)
          for (int j = 0; j < h_matrix.cols; j++) {
            fscanf(fp_occu, "%d", &h_matrix.occupied[i][j]);
            fscanf(fp_fade, "%d", &h_matrix.fade[i][j]);
            fscanf(fp_elem, "%d", &h_matrix.element[i][j]);
          }
      } else {
        // Trim: keep target_k payload cols from left, all M parity cols from right,
        // skip (file_n - h_matrix.cols) payload cols in between.
        // File layout:  [file_k payload] [M parity]   (file_k = file_n - M)
        // Target layout: [target_k payload] [M parity] (target_k = h_matrix.cols - M)
        int target_k = h_matrix.cols - h_matrix.rows;
        int col_offset = file_n - h_matrix.cols;
        int tmp_occu, tmp_fade, tmp_elem;
        for (int i = 0; i < h_matrix.rows; i++) {
          for (int j = 0; j < file_n; j++) {
            fscanf(fp_occu, "%d", &tmp_occu);
            fscanf(fp_fade, "%d", &tmp_fade);
            fscanf(fp_elem, "%d", &tmp_elem);
            int dst = -1;
            if (j < target_k)
              dst = j;                    // payload: direct map
            else if (j >= file_n - h_matrix.rows)
              dst = j - col_offset;       // parity: shift left by trimmed count
            if (dst >= 0) {
              h_matrix.occupied[i][dst] = tmp_occu;
              h_matrix.fade[i][dst] = tmp_fade;
              h_matrix.element[i][dst] = tmp_elem;
            }
          }
        }
      }
      fclose(fp_occu); fclose(fp_fade); fclose(fp_elem);
    } 
    // else if ((customize_config1 || customize_config2) && (use_customized_matrix == 2))
    // {
    //   printf("Read customized matrix generated by Chunhua\n");
    //   FILE *fp_elem;
    //   if (customize_config1)
    //     fp
    // }


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
    // if (VERBOSITY > 0)
    //f_print_h_matrix(h_matrix);
  }
} // ldpc_config


// LDPC decoder config
void ldpc_packet::ldpc_dec_config(int max_ldec_itr, float dec_alpha, int fin_mode,
                                  int fin_q_num, int fin_r_num, int fin_f_num) {
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

} // ldpc_dec_config


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
        unpacked_bit_location++;
      } else if (j == (h_matrix.cols - h_matrix.rows)) {
        if (k < (512 - skip_parity_bits_in_first_parity_column)) {
          tx_blk[bit_location++] = ldpc_encoder_output.c[j].b[k];
        }
        unpacked_bit_location++;
      } else if (j > (h_matrix.cols - h_matrix.rows)) {
        tx_blk[bit_location++] = ldpc_encoder_output.c[j].b[k];
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

void ldpc_packet::ldpc_decoder() {
  // add 0 padding
  
  vec_copy(det_blk, dec_di_blk, 0, 0, info_len);
  for (int i = 0; i < pad_len; i++)
    dec_di_blk[info_len + i] = max_llr_bin;
  vec_copy(det_blk, dec_di_blk, info_len, hm_k, h_matrix.extra_bits_of_parity);
  for (int i = 0; i < h_matrix.unused_bytes_of_parity * 8; i++)
    dec_di_blk[hm_k + h_matrix.extra_bits_of_parity + i] = max_llr_bin;
  vec_copy(det_blk, dec_di_blk, info_len + h_matrix.extra_bits_of_parity, hm_k + h_matrix.bits, (h_matrix.rows - 1) * h_matrix.bits);
  // `LAYER_G2` supports IBEX shortening mode where the first parity column is fractional
  // (`extra_bits_of_parity > 0`). For a full-parity codeword (`extra_bits_of_parity == 0`),
  // parity is contiguous `hm_m` bits and must start from the first parity column.
  // if (h_matrix.extra_bits_of_parity > 0) {
  //   vec_copy(det_blk, dec_di_blk, info_len, hm_k, h_matrix.extra_bits_of_parity);
  //   for (int i = 0; i < h_matrix.unused_bytes_of_parity * 8; i++)
  //     dec_di_blk[hm_k + h_matrix.extra_bits_of_parity + i] = max_llr_bin;
  //   vec_copy(det_blk, dec_di_blk, info_len + h_matrix.extra_bits_of_parity, hm_k + cir_sz, (bm_m - 1) * cir_sz);
  // } else {
  //   vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);
  // }
  
  // decoder core
  ldpc_dec_layer2(h_matrix);

  // NOTE: for IBEX (cir_sz==512), the transmitted block `det_blk` does not contain the padded tail bits in
  // the first parity column when `unused_bytes_of_parity > 0`. The decoder input `dec_di_blk` is always
  // sized to `hm_n` (full QC length), so we must explicitly insert those missing bits instead of blindly
  // copying `hm_m` bits (which would read past `det_blk`).
  // if ((cir_sz == 512) && (h_matrix.unused_bytes_of_parity > 0)) {
  //   const int skip_parity_bits_in_first_parity_column = h_matrix.unused_bytes_of_parity * 8;
  //   const int valid_parity_bits_in_first_parity_column = h_matrix.bits - skip_parity_bits_in_first_parity_column;
  //   const int parity_bits_in_det = blk_len - info_len;

  //   if (valid_parity_bits_in_first_parity_column > 0) {
  //     vec_copy(det_blk, dec_di_blk, info_len, hm_k, valid_parity_bits_in_first_parity_column);
  //   }

  //   for (int i = 0; i < skip_parity_bits_in_first_parity_column; i++)
  //     dec_di_blk[hm_k + valid_parity_bits_in_first_parity_column + i] = max_llr_bin;

  //   const int remaining_parity_bits = parity_bits_in_det - valid_parity_bits_in_first_parity_column;
  //   if (remaining_parity_bits > 0) {
  //     vec_copy(det_blk, dec_di_blk, info_len + valid_parity_bits_in_first_parity_column, hm_k + h_matrix.bits,
  //              remaining_parity_bits);
  //   }
  // } else {
  //   vec_copy(det_blk, dec_di_blk, info_len, hm_k, hm_m);
  // }

  // remove padding
  vec_copy(dec_do_blk, dec_blk, 0, 0, info_len);
  if (h_matrix.extra_bits_of_parity > 0)
  {
    vec_copy(dec_do_blk, dec_blk, hm_k, info_len, h_matrix.extra_bits_of_parity);
    vec_copy(dec_do_blk, dec_blk, hm_k + h_matrix.bits, info_len + h_matrix.extra_bits_of_parity, (h_matrix.rows - 1) * h_matrix.bits);
  } else {
    vec_copy(dec_do_blk, dec_blk, hm_k, info_len, hm_m);
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

void ldpc_packet::ldpc_dec_layer2(s_h_matrix h_matrix) 
{
  int col, layer_pre;
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
  int **e_pre;

  char *layer_synd;
  char *cn_dec_hd;
  char *vn_dec_hd;
  int hd_updated;
  int layer_synd_wt;
  int synd_pass_cnt = 0;
  int hd_stable_cnt = 0;

  // allocation
  dec_init = (char *)calloc(h_matrix.cols, sizeof(*dec_init));
  vec_set(dec_init, h_matrix.cols);

  cn_c_mem = (struct cn_msg **)calloc(h_matrix.rows, sizeof(*cn_c_mem));
  for (int i = 0; i < h_matrix.rows; i++)
    cn_c_mem[i] = (struct cn_msg *)calloc(h_matrix.bits, sizeof(*cn_c_mem[i]));
  cn_c_updt_cur = (struct cn_msg *)calloc(h_matrix.bits, sizeof(*cn_c_updt_cur));

  cn_q_mem = (float **)calloc(h_matrix.cols, sizeof(*cn_q_mem));
  for (int i = 0; i < h_matrix.cols; i++)
    cn_q_mem[i] = (float *)calloc(h_matrix.bits, sizeof(*cn_q_mem[i]));

  cn_r_new_pre  = (float *)calloc(h_matrix.bits, sizeof(*cn_r_new_pre));
  cn_app_pre    = (float *)calloc(h_matrix.bits, sizeof(*cn_app_pre));
  cn_app_cur    = (float *)calloc(h_matrix.bits, sizeof(*cn_app_cur));
  cn_q_sel_cur  = (float *)calloc(h_matrix.bits, sizeof(*cn_q_sel_cur));
  cn_r_old_cur  = (float *)calloc(h_matrix.bits, sizeof(*cn_r_old_cur));
  cn_q_updt_cur = (float *)calloc(h_matrix.bits, sizeof(*cn_q_updt_cur));

  cn_q_sign = (int **)calloc(5 * h_matrix.cols, sizeof(*cn_q_sign));
  for (int i = 0; i < 5 * h_matrix.cols; i++)
    cn_q_sign[i] = (int *)calloc(h_matrix.bits, sizeof(*cn_q_sign[i]));

  e_pre = (int **)calloc(h_matrix.rows, sizeof(*e_pre));
  for (int i = 0; i < h_matrix.rows; i++)
    e_pre[i] = (int *)calloc(h_matrix.cols, sizeof(*e_pre[i]));

  layer_synd = (char *)calloc(h_matrix.bits, sizeof(*layer_synd));
  vn_dec_hd  = (char *)calloc(h_matrix.bits, sizeof(*vn_dec_hd));
  cn_dec_hd  = (char *)calloc(h_matrix.bits, sizeof(*cn_dec_hd));

  // initialize decoder
  cw_fail = 1;
  cw_miscorr = 0;
  vec_copy(dec_di_blk, dec_do_blk, 0, 0, hm_n);

  for (int i = 0; i < h_matrix.rows; i++)
    for (int j = 0; j < h_matrix.cols; j++)
      e_pre[i][j] = -1;

  for (int i = 0; i < h_matrix.cols; i++) 
  {
    int tmp_pre = -1;
    int first_layer = -1;
    for (int j = 0; j < h_matrix.rows; j++)
    {
      if ((h_matrix.occupied[j][i]) || ((h_matrix.fade[j][i]) && h_matrix.extra_bytes_of_parity > 0)) 
      {
        e_pre[j][i] = tmp_pre;
        tmp_pre = j;
        if (first_layer < 0)
          first_layer = j;
      }
    }
    e_pre[first_layer][i] = tmp_pre;
  }

  for (int i = 0; i < h_matrix.cols; i++)
    for (int j = 0; j < h_matrix.bits; j++)
      cn_q_mem[i][j] = (float)llr_tbl[dec_di_blk[i * h_matrix.bits + j]];

  // iterative decoding
  for (int itr = 0; (itr < ldec_max_itr) && ((ldec_early_term_en == 0) || (cw_fail == 1)); itr++) {
    // Q sign mem index
    cir_cnt = 0;

    // layer decoding
    for (int layer = 0; layer < h_matrix.rows && ((ldec_early_term_en == 0) || (cw_fail == 1)); layer++) {
      // initilize HD mem
      hd_init = (vec_sum(dec_init, h_matrix.cols) != 0);

      // #ifdef _LDPC_DEBUG_DUMP
      //       printf("[LDPC DEBUG] Layer decoding @ iteration %d, layer %d ...\n", itr, layer);
      // #endif

      // init current layer C-MSG
      // C-MSG of previous iteration

      cn_c_sel_cur = cn_c_mem[layer];
      // C-MSG to be updt
      for (int i = 0; i < h_matrix.bits; i++) {
        cn_c_updt_cur[i].min1_val = 100000;
        cn_c_updt_cur[i].min2_val = 100000;
        cn_c_updt_cur[i].min1_pos = 0;
        cn_c_updt_cur[i].sign_tot = 1;
      }

      // Earlier termination init
      hd_updated = 0;
      vec_clr(layer_synd, h_matrix.bits);

      // Per circulant of the layer
      for (col=0; col < h_matrix.cols; col++) {
        if (((h_matrix.occupied[layer][col] == 0) && (h_matrix.fade[layer][col] == 0))
        || ((h_matrix.fade[layer][col] == 1)&& (h_matrix.extra_bytes_of_parity==0)))
          continue;


        cn_q_sel_pre = cn_q_mem[col];
        layer_pre = e_pre[layer][col];
        cn_c_sel_pre = cn_c_mem[layer_pre]; // read previous layer C msg

        // cal Rnew and APP
        for (int i = 0; i < h_matrix.bits; i++) {
          // Qmsg sign
          sign_tmp = (cn_q_sel_pre[i] >= 0) ? 1 : -1;

          // Rnew
          if (cn_c_sel_pre[i].min1_pos == col)
            cn_r_new_pre[i] = cn_c_sel_pre[i].min2_val * cn_c_sel_pre[i].sign_tot * sign_tmp;
          else
            cn_r_new_pre[i] = cn_c_sel_pre[i].min1_val * cn_c_sel_pre[i].sign_tot * sign_tmp;

          // APP in CN order of previous layer
          if (h_matrix.extra_bytes_of_parity == 0) {
            cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
          } else {
            if (h_matrix.occupied[layer_pre][col] && (layer_pre < (h_matrix.rows - 1)))
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else if (h_matrix.occupied[layer_pre][col] && (layer_pre == (h_matrix.rows - 1)) &&
                     h_matrix.mask[col][(i + h_matrix.element[layer_pre][col]) % h_matrix.bits])
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else if (h_matrix.fade[layer_pre][col] && !h_matrix.mask[col][(i + h_matrix.element[layer_pre][col]) % h_matrix.bits])
              cn_app_pre[i] = cn_r_new_pre[i] + cn_q_sel_pre[i];
            else
              cn_app_pre[i] = cn_q_sel_pre[i];
          }

          // Quantization
          if (finite_mode == 1) {
            cn_app_pre[i] =
                (float)Sat_Quan((double)cn_app_pre[i], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
          }
        }

        // APP shift
        // when decoder initilized, Q msg are in VN order
        if (dec_init[col] == 1) {
          shift_val1 = h_matrix.element[layer][col];
          shift_val2 = 0;
          dec_init[col] = 0;
        } else {
          shift_val1 = -1 * h_matrix.element[layer_pre][col] + h_matrix.element[layer][col];
          shift_val2 = -1 * h_matrix.element[layer_pre][col];
        }
        for (int i = 0; i < h_matrix.bits; i++) {
          cn_app_cur[i] = cn_app_pre[(i + shift_val1 + h_matrix.bits) % h_matrix.bits];
          vn_dec_hd[i]  = cn_app_pre[(i + shift_val2 + h_matrix.bits) % h_matrix.bits] >= 0 ? 0 : 1;
        }

        // CW converge check logic per circulant
        // 1. check if HD updated
        if (hd_updated == 0)
          if (vec_cmp(dec_do_blk, vn_dec_hd, col * h_matrix.bits, 0, h_matrix.bits) == 1)
            hd_updated = 1;

        vec_copy(vn_dec_hd, dec_do_blk, 0, col * h_matrix.bits, h_matrix.bits);

        // 2. accumulate syndrome
        vec_shift(vn_dec_hd, cn_dec_hd, h_matrix.bits, -1 * h_matrix.element[layer][col]);
        if (h_matrix.extra_bytes_of_parity == 0) {
          vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, h_matrix.bits);
        } else {
          if (h_matrix.occupied[layer][col] && (layer < (h_matrix.rows - 1))) {
            for (int i = 0; i < h_matrix.bits; i++) {
              if (!h_matrix.mask[col][(i + h_matrix.element[layer][col]) % h_matrix.bits])
                cn_dec_hd[i] = 0;
            }
          }
          if (h_matrix.fade[layer][col]) {
            for (int i = 0; i < h_matrix.bits; i++) {
              if (!h_matrix.mask[col][(i + h_matrix.element[layer][col]) % h_matrix.bits])
                cn_dec_hd[i] = 0;
            }
          }
          vec_mod2_add(cn_dec_hd, layer_synd, layer_synd, h_matrix.bits);
        }

        // calculate R_old and current Q, update current layer C and Q
        for (int i = 0; i < h_matrix.bits; i++) {
          if (cn_c_sel_cur[i].min1_pos == col)
            cn_r_old_cur[i] = cn_c_sel_cur[i].min2_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];
          else
            cn_r_old_cur[i] = cn_c_sel_cur[i].min1_val * cn_c_sel_cur[i].sign_tot * cn_q_sign[cir_cnt][i];

          // Q -= Rold
          if (h_matrix.extra_bytes_of_parity == 0) {
            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          } else {
            if (h_matrix.occupied[layer][col] && (layer == (h_matrix.rows - 1)) &&
                (!h_matrix.mask[col][(i + h_matrix.element[layer][col]) % h_matrix.bits])) {
              cn_r_old_cur[i] = 0;
            }
            if (h_matrix.fade[layer][col] && (h_matrix.mask[col][(i + h_matrix.element[layer][col]) % h_matrix.bits])) {
              cn_r_old_cur[i] = 0;
            }
            cn_q_updt_cur[i] = cn_app_cur[i] - cn_r_old_cur[i];
          }

          // Quantization
          if (finite_mode == 1) {
            cn_q_updt_cur[i] =
                (float)Sat_Quan((double)cn_q_updt_cur[i], finite_q_max, finite_q_min, finite_q_num, finite_f_num);
          }

          // update C
          // mask
          if (h_matrix.extra_bits_of_parity == 0) {
            sign_tmp = (cn_q_updt_cur[i] >= 0) ? 1 : -1;
            val_tmp = cn_q_updt_cur[i] * sign_tmp;
          } else {
            if (h_matrix.occupied[layer][col] && (layer == (h_matrix.rows - 1)) &&
                (!h_matrix.mask[col][(i + h_matrix.element[layer][col]) % h_matrix.bits])) {
              sign_tmp = 1;
              val_tmp = 100000;
            }
            if (h_matrix.fade[layer][col] && (h_matrix.mask[col][(i + h_matrix.element[layer][col]) % h_matrix.bits])) {
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
            cn_c_updt_cur[i].min1_pos = col;
          } else if (val_tmp < cn_c_updt_cur[i].min2_val) {
            cn_c_updt_cur[i].min2_val = val_tmp;
          }

          cn_q_sign[cir_cnt][i] = sign_tmp;
        }

        // update Q memory
        for (int i = 0; i < h_matrix.bits; i++) {
          cn_q_mem[col][i] = cn_q_updt_cur[i];
        }

        cir_cnt++;
      } // per circulant

      // update C_MSG per layer
      for (int i = 0; i < h_matrix.bits; i++) {
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
      layer_synd_wt = vec_sum(layer_synd, h_matrix.bits);
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
      if ((synd_pass_cnt >= h_matrix.rows) && (hd_stable_cnt >= h_matrix.rows - 1)) {
        cw_fail = 0;
        cnvg_itr = itr;
        cnvg_lyr = layer;
      }
    }
  }

  if ((cw_fail == 1) || (ldec_early_term_en == 0)) {
    cnvg_itr = ldec_max_itr - 1;
    cnvg_lyr = h_matrix.rows - 1;
  }

  // free all
  free(dec_init);
  for (int i = 0; i < h_matrix.rows; i++)
    free(cn_c_mem[i]);
  free(cn_c_mem);
  free(cn_c_updt_cur);
  for (int i = 0; i < h_matrix.cols; i++)
    free(cn_q_mem[i]);
  free(cn_q_mem);
  free(cn_r_new_pre);
  free(cn_app_pre);
  free(cn_app_cur);
  free(cn_q_sel_cur);
  free(cn_r_old_cur);
  free(cn_q_updt_cur);
  for (int i = 0; i < 5 * h_matrix.cols; i++)
    free(cn_q_sign[i]);
  free(cn_q_sign);
  for (int i = 0; i < h_matrix.rows; i++)
    free(e_pre[i]);
  free(e_pre);
  free(layer_synd);
  free(cn_dec_hd);
  free(vn_dec_hd);
} // ldpc_dec_layer2