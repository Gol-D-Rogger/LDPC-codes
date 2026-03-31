#pragma once

#include <cstdint>

// LDPC basic parameters definition (add LDPC_ prefix to avoid naming conflicts)
#define LDPC_L 5
#define LDPC_M 13
#define LDPC_U 67
#define LDPC_N (LDPC_M + LDPC_U) // total number of columns
#define LDPC_P 512
#define LDPC_Pm1 511

#define DEFAULT_MATRIX 0

struct h_matrix {
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

  int delta[LDPC_M];
  short element[LDPC_M][LDPC_N];
  int first_element[LDPC_M];
  int last_element[LDPC_M];
  int wraparound[LDPC_M];
  int wrap_base[LDPC_M];
  int wrap_num_deltas[LDPC_M];
  int row_weight[LDPC_M];
  int col_weight[LDPC_N];

  bool occupied[LDPC_M][LDPC_N];
  bool fade[LDPC_M][LDPC_N];
  bool last_row_active_bits[LDPC_N][LDPC_P];
  bool parity_column[LDPC_N];
  int bits_in_last_column;

  struct range {
    unsigned short active_count;
    unsigned short offset;
  };

  unsigned char operational_h_matrix[LDPC_M][LDPC_N];
  range tile_ranges[2 * LDPC_N];

  h_matrix();
  h_matrix(const int bytes_of_userdata, const int bytes_of_parity,
           int matrix_sel = DEFAULT_MATRIX);

  void print(void);
  void setup_range_operational_matrix(void);
};
