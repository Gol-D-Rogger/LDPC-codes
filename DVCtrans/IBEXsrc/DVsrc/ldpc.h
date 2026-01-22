#define L   5
#define M   13
#define U   67
#define N   (M+U)
#define P   2048
#define SCALE 1
#define VN_BITS     8
#define VN_MIN      0
#define VN_MAX      255
#define VN_THR      248

#include <stdio.h>
#include <random>

struct s_h_matrix
{
    int rows;
    int cols;
    int bits;
    int column_weight;
    int bytes_of_userdata;
    int bytes_of_parity;
    int extra_bytes_of_parity;
    int extra_bytes_of_userdata;
    int unused_bytes_of_parity;
    int unused_bytes_of_userdata;
    int extra_bits_of_parity;
    int extra_bits_of_userdata;
    int min_rows;
    int max_rows;
    int delta[M];
    int element[M][N];
    int first_element[M];
    int last_element[M];
    int wraparound[M];
    int wrap_base[M];
    int wrap_num_deltas[M];
    int row_weight[M];
    int col_weight[N];
    bool occupied[M][N];
    bool fade[M][N];
    bool mask[N][P];
    bool parity_column[N];
    int bits_in_last_column;
};

struct s_h_matrix_x
{
    s_h_matrix x[4]; // extra parity columns 0 to 3
};

struct s_h_matrix_all
{
    s_h_matrix_x r[M-(L-1)]; // rows L to M
    int delta[M];
    int bits;
};

struct s_ldpc_matrix
{
    int rows[M];
    int cols[4][M];
    int delta[M];
    int element[4][M][M][N];
    int first_element[4][M][M];
    int last_element[M];
    int wraparound[4][M][M];
    int wrap_num_deltas[4][M][M];
    int row_weight[4][M][M];
    int col_weight[4][M][N];
    bool occupied[4][M][M][N];
    bool fade[4][M][M][N];
};

struct s_variable_node_bit
{
    short int likelihood;
    bool bit_hard;
    bool bit_questionable;
    bool bit_questionable2;
    short int level;
    bool flipped;
};

struct s_variable_node_column
{
    s_variable_node_bit b[P];
};

struct s_variable_nodes
{
    s_variable_node_column c[N];
};

struct s_check_node_row
{
    bool b[P];
};

struct s_check_nodes
{
    s_check_node_row r[M];
};

struct s_codeword_bit
{
    bool bit_hard;
    bool bit_questionable;
    bool bit_questionable2;
    short int level;
    short int syndrome_weight;
    bool bit_corrected;
    bool bit_is_error;
};

struct s_codeword_column
{
    s_codeword_bit b[P];
};

struct s_codeword
{
    s_codeword_column c[N];
};

struct s_hard_codeword_column
{
    bool b[P];
};

struct s_hard_codeword
{
    s_hard_codeword_column c[N];
    unsigned int errors_at_level_and_weight[4][5];
    unsigned int correct_at_level_and_weight[4][5];
    double probability_of_error_at_level_and_weight[4][5];
};

struct s_ldpc_decoder_parameters
{
    int post_process_en;
    int syndrome_weight_thr_qc;
    int syndrome_weight_thr_post;
    int likelihood_thr;
    int likelihood_init_coef_all[8][4];
    int likelihood_init_coef[4];
    int likelihood_init_fraction[3];
    int soft_bit_table[8];
    int post_ratio;
    int likelihood_map[8];
    bool early_terminate_dis;
    int early_terminate_thr[2][8];
    bool questionable_sense;
};

struct s_ldpc_decoder_input
{
    int nand_strobes;
    int soft_bits;
    int post_iteration;
    int iteration_limit;
    bool syndrome_cal_only;
    int errors_in_userdata;
    int errors;
    float rber;
    bool verbose;
    s_codeword corrupted_codeword;
};

struct s_ldpc_decoder_output
{
    int failure;
    int errors_in_userdata;
    int errors_in_codeword;
    int iterations;
    int clock_cycles;
    int syndrome_weight_before;
    int syndrome_weight_after;
    int early_termination;
    int col_cnt; // 0-based column index where decoding converged; -1 if not converged or no column processed
    s_hard_codeword corrected_codeword;
};

struct s_256_bits
{
    bool b[256];
};

struct s_512_bits
{
    bool b[512];
};

struct s_likelihood_levels
{
    int min;
    int max;
    int flip_thr;
    int weak;
    int strong;
    int level[4];
};

s_check_nodes cn;
s_h_matrix h_matrix;
s_ldpc_decoder_input ldpc_decoder_input;
s_ldpc_decoder_output ldpc_decoder_output;
s_ldpc_decoder_parameters ldpc_decoder_parameters;
s_ldpc_matrix ldpc_matrix;

void f_print_h_matrix(s_h_matrix h_matrix) {
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

// LDPC code configuration
s_h_matrix f_h_matrix(int bytes_of_userdata, int bytes_of_parity)
{
    int VERBOSITY = 0;
    int i;
    int j;
    int k;
    int bit;
    int n;
    bool m[M][N];
    int h_matrix_index[13];
    s_h_matrix h_matrix;

    h_matrix.bits = (bytes_of_userdata > 5000) ? 1024 : 512;
    h_matrix.rows = (h_matrix.bits == 1024) ? ((bytes_of_parity + 127) >> 7) : ((bytes_of_parity + 63) >> 6);
    h_matrix.cols = (h_matrix.bits == 1024) ? ((bytes_of_userdata + 127) >> 7) : ((bytes_of_userdata + 63) >> 6);
    h_matrix.cols += h_matrix.rows; // total columns = userdata + parity
    h_matrix.bytes_of_userdata = bytes_of_userdata;
    h_matrix.bytes_of_parity = bytes_of_parity;
    h_matrix.min_rows = L;
    h_matrix.max_rows = M;

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
    return h_matrix;
}

s_check_nodes f_check_nodes_clear(s_h_matrix h_matrix) {
  s_check_nodes cn;
  int i, j;
  for (i = 0; i < h_matrix.rows; i++) {
    for (j = 0; j < h_matrix.bits; j++) {
      cn.r[i].b[j] = 0;
    }
  }
  return cn;
}

void f_print_check_nodes(s_check_nodes data, int rows, int bits) {
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

void f_print_check_nodes_shifted(s_check_nodes data, s_h_matrix h_matrix, int column) {
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

int f_check_node_weight(s_h_matrix h_matrix, s_check_nodes cn) {
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

s_check_nodes f_check_nodes(s_h_matrix h_matrix, s_hard_codeword vn) {
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

s_256_bits f_256_bit_lfsr(s_256_bits data_in) {
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

s_512_bits f_512_bit_lfsr(s_512_bits data_in) {
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

void f_print_hard_codeword(s_hard_codeword data, int cols, int bits) {
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

void f_print_s_256_bits(s_256_bits s) {
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

void f_print_s_512_bits(s_512_bits s) {
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

// Function: IBEX Likelihood level Initial
s_likelihood_levels f_likelihood_levels(int strobes, s_ldpc_decoder_parameters ldpc_decoder_parameters,
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

int f_update_vn_post(int likelihood, int weight, int min_likelihood, int max_likelihood, bool post_process,
                                  bool post_process2, bool be_aggressive, int flip_threshold, bool pushing) {
  int likelihood_new;
  bool clamp_to_min;
  bool do_post_flipped;
  bool do_post_unflipped;
  bool do_aggr;
  bool flipped = (likelihood >= flip_threshold);

  int delta;
  if (VN_BITS <= 2)
    delta = (weight>>1) + ((weight&1) && pushing ? 1 : 0);
  else
    delta = weight;

  likelihood_new = flipped ? (likelihood - delta) : (likelihood + delta - 1); // when clamp!=1, same as strong
  // likelihood_new = flipped ? (likelihood - weight) : (likelihood + weight - 1); // when clamp!=1, same as strong

  bool flipped_new = (likelihood_new >= flip_threshold);
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

s_ldpc_decoder_output f_ldpc_decode(s_ldpc_decoder_input ldpc_decoder_input,
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

#ifdef _LDPC_DBG_DUMP
  static int dbg_made_dir = 0;
  if (!dbg_made_dir) { system("mkdir -p output >/dev/null 2>&1"); dbg_made_dir = 1; }
  FILE* dbg_sw_fp = fopen("./output/ldpc_dbg_sw.txt", "w");
  if (dbg_sw_fp)
    fprintf(dbg_sw_fp, "%6s %6s\n", "iter", "sw");
  FILE* dbg_tr_fp = NULL;
#endif

  s_check_nodes cn_shifted;
  s_likelihood_levels likelihood_levels;
  s_ldpc_decoder_output ldpc_decoder_output;
  ldpc_decoder_output.col_cnt = -1;
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

  int prev_sw = (iteration == 0) ? syndrome_weight_delayed : syndrome_weight_r[3];
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
        printf("### LDPC DECODER C++ EARLY TERMINATION: %5d\n", syndrome_weight);
    }
  }
  if (ldpc_decoder_input.syndrome_cal_only) {
    finished = 1;
  }

  likelihood_levels =
      f_likelihood_levels(ldpc_decoder_input.nand_strobes, ldpc_decoder_parameters, syndrome_weight, h_matrix.rows);

  // soft_data -> likelihood_level
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

  while ((iteration < ldpc_decoder_input.iteration_limit) && (finished == 0) && (give_up == 0)) {
#ifdef _LDPC_DBG_DUMP
    {
      char dbg_tr_path[128];
      snprintf(dbg_tr_path, sizeof(dbg_tr_path), "./output/ldpc_dbg_trace_iter%d.txt", iteration);
      dbg_tr_fp = fopen(dbg_tr_path, "w");
      if (dbg_tr_fp)
        fprintf(dbg_tr_fp, "%6s %6s %4s %3s %1s\n", "iter", "idx", "like", "wt", "f");
    }
#endif
    for (j = 0; j < h_matrix.cols; j++) {
      clock_cycles++;
      if ((iteration == (ldpc_decoder_input.post_iteration + 0)) && (j == 0)) {
        for (i = 0; i < 256; i++)
          prng_256.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
        for (i = 0; i < 512; i++)
          prng_512.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
        for (i = 0; i < 512; i++)
          prng_512.b[i] = (0x1fe0 >> (i % 16)) & 1; // to match verilog
      } else if (iteration >= ldpc_decoder_input.post_iteration) {
        prng_256 = f_256_bit_lfsr(prng_256);
        prng_512 = f_512_bit_lfsr(prng_512);
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
      for (k = 0; k < h_matrix.bits; k++) {
        look = 0;
        do_not_use_this_bit = 0;
        do_not_use_this_bit |= ((h_matrix.extra_bits_of_parity > 0) && (j == (h_matrix.cols - h_matrix.rows)) &&
                                (k >= h_matrix.extra_bits_of_parity));
        do_not_use_this_bit |= ((h_matrix.extra_bits_of_userdata > 0) && (j == (h_matrix.cols - h_matrix.rows - 1)) &&
                                (k >= h_matrix.extra_bits_of_userdata));
        if ((VERBOSITY > 0) && do_not_use_this_bit)
          printf("### DO NOT USE THIS BIT %2d %3d\n", j, k);
        if (do_not_use_this_bit) {
#ifdef _LDPC_DBG_DUMP
          if (dbg_tr_fp) {
            int idx = (j * h_matrix.bits) + k;
            fprintf(dbg_tr_fp, "%6d %6d %4d %3d %1d\n", iteration, idx, (int)vn.c[j].b[k].likelihood, -1,
                    (int)vn.c[j].b[k].flipped);
          }
#endif
        } else {
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
          bool aggr = (ldpc_decoder_input.soft_bits > 0) &&
                      (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr && !flipped_prev);
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
            printf("C++ LOOK   ITERATION: %4d BEFORE LIKELIHOOD UPDATE SW: %4d CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x "
                   "CORRUPTED: %x FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: %3d\n",
                   iteration, syndrome_weight, j, k, 0, ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                   vn.c[j].b[k].flipped, weight, likelihood_levels.min, vn.c[j].b[k].likelihood);
          vn.c[j].b[k].likelihood =
              f_update_vn_post(vn.c[j].b[k].likelihood, w, likelihood_levels.min, likelihood_levels.max,
                               prng_post_process, prng_post_process2, aggr, likelihood_levels.flip_thr, pushing);

          if ((VERBOSITY > 0) && look)
            printf("C++ LOOK   ITERATION: %4d AFTER LIKELIHOOD UPDATE SW: %4d CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x "
                   "CORRUPTED: %x FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: %3d\n",
                   iteration, syndrome_weight, j, k, 0, ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                   vn.c[j].b[k].flipped, weight, likelihood_levels.min, vn.c[j].b[k].likelihood);
	          vn.c[j].b[k].flipped = (vn.c[j].b[k].likelihood >= likelihood_levels.flip_thr);

#ifdef _LDPC_DBG_DUMP
	          if (dbg_tr_fp) {
	            int idx = (j * h_matrix.bits) + k;
	            fprintf(dbg_tr_fp, "%6d %6d %4d %3d %1d\n", iteration, idx, (int)vn.c[j].b[k].likelihood, weight,
	                    (int)vn.c[j].b[k].flipped);
	          }
#endif

	          if ((VERBOSITY > 0) && look)
	            printf("C++ ITERATION: %4d AFTER FLIP CHECK  SW: %4d CURRENT BIT FLIP: %2d %3d  ORIGINAL: %x CORRUPTED: %x "
	                   "FLIPPED: %x WEIGHT: %1d MIN: %3d LIKELIHOOD: %3d THR: %3d\n",
	                   iteration, syndrome_weight, j, k, 0, ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard,
                   vn.c[j].b[k].flipped, weight, likelihood_levels.min, vn.c[j].b[k].likelihood,
                   likelihood_levels.flip_thr);

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

      // 
      for (i = 0; i < h_matrix.rows; i++) {
        if (h_matrix.occupied[i][j] || h_matrix.fade[i][j])
          for (k = 0; k < h_matrix.bits; k++)
            cn_shifted.r[i].b[k] = cn.r[i].b[(k + h_matrix.bits - h_matrix.element[i][j]) % h_matrix.bits];
      }
      if (VERBOSITY > 0) {
        printf("### DECODER C++: DECODING ITERATION %4d AT COLUMN %2d, CODEWORD HAS CHECK NODES:\n", iteration, j);
        f_print_check_nodes(cn, h_matrix.rows, h_matrix.bits);
        printf("### DECODER C++: DECODING ITERATION %4d AT COLUMN %2d, CODEWORD HAS CHECK NODES SHIFTED:\n", iteration, j);
        f_print_check_nodes_shifted(cn, h_matrix, j);
      }

      syndrome_weight = f_check_node_weight(h_matrix, cn);
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
		    snprintf(dbg_cn_path, sizeof(dbg_cn_path), "./output/ldpc_dbg_cn_synd_iter%d.txt", iteration);
		    FILE* dbg_cn_fp = fopen(dbg_cn_path, "w");
		    if (dbg_cn_fp) {
		      fprintf(dbg_cn_fp, "%6s %8s\n", "idx", "syndrome");
		      for (int rr = 0; rr < h_matrix.rows; rr++) {
		        for (int kk = 0; kk < h_matrix.bits; kk++) {
		          const int idx = (rr * h_matrix.bits) + kk;
		          fprintf(dbg_cn_fp, "%6d %8d\n", idx, (int)cn.r[rr].b[kk]);
		        }
		      }
		      fclose(dbg_cn_fp);
		    }
#endif

#ifdef _LDPC_DBG_DUMP
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

#ifdef _LDPC_DBG_DUMP
	  if (dbg_sw_fp)
	    fclose(dbg_sw_fp);
	  if (dbg_tr_fp)
	    fclose(dbg_tr_fp);
#endif

	  return ldpc_decoder_output;
	}

s_hard_codeword f_ldpc_encode(s_hard_codeword ldpc_encoder_input, s_h_matrix h_matrix) {
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

  bool payload[U][512];
  bool parity[M][512];
  bool check_node[M][512];
  bool ldpc_matrix_occupied[M][N];
  bool ldpc_matrix_fade[M][N];
  unsigned int ldpc_matrix[M][N];
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

//////////////////////////// ERROR INJECTION FUNCTION ////////////////////////////
double normalCDF(double value)
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}

double normalInverseCDF(double value)
{
    double x = 1.0;
    double y = 0.0;
    for (int i = 0; i < 10000; i++)
    {
        y = normalCDF(x);
        x = x + (0.25 * (value - y));
    }
    return x;
}

struct s_rbers
{
    float rber[7];
};

struct s_col_bit
{
    int col;
    int bit;
};

struct s_row_bit
{
    int col;
    int bit;
};

struct s_error_pattern
{
    int errors;
    s_col_bit location[P];
};

struct s_distribution
{
    int errors_in_bin[8];
    int correct_in_bin[8];
};

double f_erf(double x)
{
    double y = 1.0 / ( 1.0 + 0.3275911 * x);
    return 1- (((((
        + 1.061405429  * y
        - 1.453152027) * y
        + 1.421413741) * y
        - 0.284496736) * y
        + 0.254829592) * y)
        * exp(-x * x);
}

// Returns the probability of x, given the distribution described by mu and sigma.
double f_pdf (double x, double mu, double sigma) 
{ 
  //Constants
  static const double pi = 3.14159265;
  return exp (-1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma *sqrt (2 *pi));
}
// Returns the probability of [-inf,x] of a gaussian distribution
double f_cdf (double x, double mu, double sigma) 
{
  return 0.5 * (1 + erf ((x - mu) / (sigma * sqrt (2.))));
}

double f_find_sigma (double x, double mu, double target_cdf) 
{
    int i = 0;
    double sigma = 1.0;
    double cdf = 0.0;
    double error = 1;
    while ((fabs(error) > 0.00000001) && (i < 10000)) {
        cdf = f_cdf (x, mu, sigma);
        error = target_cdf - cdf;
        sigma *= (1.0 - (1.0 * error));
        i++;
    }
    return sigma;
}

double f_find_x(double mu, double sigma, double target_cdf) 
{
    int i = 0;
    double x = mu;
    double cdf = 0.0;
    double error = 1;
    while ((fabs(error) > 0.00000001) && (i < 10000)) {
        cdf = f_cdf (x, mu, sigma);
        error = target_cdf - cdf;
        x += (0.1 * error);
        i++;
    }
    return x;
}

s_rbers f_find_error_probabilities(int strobes, float rber)
{
    float vt[7];
    float cdf[7];
    s_rbers rbers;
    float scale;
    float Eb_over_N0;
    float Eb_over_N0_dB;
    float distribution_center = -1.0;
    float rber_scaled;
    float target_cdf;
    float target_rber;
    float vt_scaled;
    float x[2];
    float mu[2];
    float sigma[2];
    rbers.rber[0] = rber;
    vt[0] = 0.0;
    mu[0] =-1.0;
    mu[1] = 1.0;
    if (strobes == 1)
    {
        vt[1] = 0.0;
        vt[2] = 0.0;
        vt[3] = 0.0;
        vt[4] = 0.0;
        vt[5] = 0.0;
        vt[6] = 0.0;
    }
    else
    {
        if (strobes == 3)
        {
            target_rber = 0.015;
            target_cdf = 1.0 - target_rber;
            sigma[0] = f_find_sigma(vt[0], mu[0], target_cdf);
            target_rber *= 0.10;
            target_cdf = 1.0 - target_rber;
            vt[2] = f_find_x(mu[0], sigma[0], target_cdf);
            vt[4] = vt[2];
            vt[6] = vt[2];
        }
        if (strobes == 5)
        {
            target_rber = 0.015;
            target_cdf = 1.0 - target_rber;
            sigma[0] = f_find_sigma(vt[0], mu[0], target_cdf);
            target_rber *= 0.05;
            target_cdf = 1.0 - target_rber;
            vt[4] = f_find_x(mu[0], sigma[0], target_cdf);
            vt[2] = 0.5 * vt[4];
            vt[6] = vt[4];
        }
        if (strobes == 7)
        {
            target_rber = 0.015;
            target_cdf = 1.0 - target_rber;
            sigma[0] = f_find_sigma(vt[0], mu[0], target_cdf);
            target_rber *= 0.03;
            target_cdf = 1.0 - target_rber;
            vt[6] = f_find_x(mu[0], sigma[0], target_cdf);
            vt[2] = (1.0 / 3.0) * vt[6];
            vt[4] = (2.0 / 3.0) * vt[6];
        }
        vt[1] = -vt[2];
        vt[3] = -vt[4];
        vt[5] = -vt[6];
        target_cdf = 1 - rber;
        sigma[0] = f_find_sigma(vt[0], mu[0], target_cdf);
        rbers.rber[0] = rber;
        for (int i = 1; i < 7; i++)
            rbers.rber[i] = 1.0 - f_cdf(vt[i], mu[0], sigma[0]);
    }
    return rbers;
}

struct s_error_injected_bit
{
    bool hard;
    bool soft0;
    bool soft1;
};

s_error_injected_bit f_inject_error(int strobes, s_rbers rbers, float random_number)
{
    s_error_injected_bit output_bit;
    output_bit.hard = (random_number < rbers.rber[0]);
    if (strobes == 3)
    {
        if ((rbers.rber[2] < random_number) && (random_number < rbers.rber[1]))
            output_bit.soft0 = 1;
        else
            output_bit.soft0 = 0;
        output_bit.soft1 = 0;
    }
    else if (strobes > 3)
    {
        if ((rbers.rber[2] < random_number) && (random_number < rbers.rber[1]))
        {
            output_bit.soft0 = 1;
            output_bit.soft1 = 1;
        }
        else if ((rbers.rber[4] < random_number) && (random_number < rbers.rber[3]))
        {
            output_bit.soft0 = 1;
            output_bit.soft1 = 0;
        }
        else if ((rbers.rber[6] < random_number) && (random_number < rbers.rber[5]))
        {
            output_bit.soft0 = 0;
            output_bit.soft1 = 1;
        }
        else
        {
            output_bit.soft0 = 0;
            output_bit.soft1 = 0;
        }
    }
    
    return output_bit;
}

void f_create_include_files()
{
    
}
