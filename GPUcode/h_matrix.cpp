#include <algorithm>
#include <cstdio>

#include "h_matrix.h"

namespace {
void init_delta(int delta[LDPC_M], int matrix_sel) {
  static const int delta_m0[LDPC_M] = {0,  13, 19, 29,  41,  67, 73,
                                       79, 91, 97, 103, 111, 119};
  static const int delta_m1[LDPC_M] = {0,   293, 61,  479, 17,  173, 53,
                                       277, 307, 229, 67,  271, 307};
  const int *src = (matrix_sel == 1) ? delta_m1 : delta_m0;
  for (int i = 0; i < LDPC_M; ++i)
    delta[i] = src[i];
}
} // namespace

h_matrix::h_matrix() {
  rows = 0;
  cols = 0;
  bits = LDPC_P;
  column_weight = 0;
  bytes_of_userdata = 0;
  bytes_of_parity = 0;
  extra_bytes_of_parity = 0;
  extra_bytes_of_userdata = 0;
  unused_bytes_of_parity = 0;
  unused_bytes_of_userdata = 0;
  extra_bits_of_parity = 0;
  extra_bits_of_userdata = 0;
  min_rows = LDPC_L;
  max_rows = LDPC_M;
  bits_in_last_column = bits;

  for (int i = 0; i < LDPC_M; ++i) {
    delta[i] = 0;
    first_element[i] = -1;
    last_element[i] = -1;
    wraparound[i] = 0;
    wrap_base[i] = 0;
    wrap_num_deltas[i] = 0;
    row_weight[i] = 0;
    for (int j = 0; j < LDPC_N; ++j) {
      element[i][j] = -1;
      occupied[i][j] = false;
      fade[i][j] = false;
      operational_h_matrix[i][j] = 0;
    }
  }

  for (int j = 0; j < LDPC_N; ++j) {
    col_weight[j] = 0;
    parity_column[j] = false;
    for (int k = 0; k < LDPC_P; ++k)
      last_row_active_bits[j][k] = false;
  }

  for (int i = 0; i < 2 * LDPC_N; ++i) {
    tile_ranges[i].active_count = 0;
    tile_ranges[i].offset = 0;
  }
}

h_matrix::h_matrix(const int bytes_of_userdata, const int bytes_of_parity,
                   int matrix_sel)
    : h_matrix() {
  this->bytes_of_userdata = bytes_of_userdata;
  this->bytes_of_parity = bytes_of_parity;

  const int bytes_per_col = bits >> 3; // 64 bytes for LDPC_P=512
  rows = std::max(1, (bytes_of_parity + bytes_per_col - 1) / bytes_per_col);
  rows = std::min(rows, LDPC_M);

  int user_cols =
      std::max(1, (bytes_of_userdata + bytes_per_col - 1) / bytes_per_col);
  user_cols = std::min(user_cols, LDPC_U);
  cols = std::min(user_cols + rows, LDPC_N);

  unused_bytes_of_parity = rows * bytes_per_col - bytes_of_parity;
  extra_bytes_of_parity = (unused_bytes_of_parity == 0)
                              ? 0
                              : (bytes_per_col - unused_bytes_of_parity);
  unused_bytes_of_userdata = (cols - rows) * bytes_per_col - bytes_of_userdata;
  extra_bytes_of_userdata = (unused_bytes_of_userdata == 0)
                                ? 0
                                : (bytes_per_col - unused_bytes_of_userdata);
  extra_bits_of_parity = extra_bytes_of_parity << 3;
  extra_bits_of_userdata = extra_bytes_of_userdata << 3;
  bits_in_last_column =
      (extra_bits_of_parity == 0) ? bits : extra_bits_of_parity;

  for (int j = 0; j < cols; ++j)
    parity_column[j] = (j < rows);

  init_delta(delta, matrix_sel);

  // Build a conservative QC-like occupancy map from OCR intent:
  // - parity part has a guaranteed diagonal
  // - user-data part uses cyclic sparse occupancy
  const int user_cols_end = cols - rows;
  for (int i = 0; i < rows; ++i) {
    const int parity_col = user_cols_end + i;
    if (parity_col >= 0 && parity_col < cols)
      occupied[i][parity_col] = true;
    for (int j = 0; j < user_cols_end; ++j) {
      bool use_location = ((j + i) % LDPC_L) != 0;
      if (j == i % LDPC_L)
        use_location = false;
      if (j < LDPC_L && i == rows - 1)
        use_location = true;
      occupied[i][j] = use_location;
    }
  }

  // Fill cyclic shifts and weights.
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (occupied[i][j] || fade[i][j]) {
        const int d = (delta[i] == 0) ? 1 : delta[i];
        element[i][j] = static_cast<short>((j * d) % bits);
        if (first_element[i] < 0)
          first_element[i] = element[i][j];
        last_element[i] = element[i][j];
      } else {
        element[i][j] = -1;
      }
      if (occupied[i][j]) {
        row_weight[i]++;
        col_weight[j]++;
      }
    }
    if (first_element[i] < 0)
      first_element[i] = 0;
    if (last_element[i] < 0)
      last_element[i] = 0;
    wraparound[i] = (bits + first_element[i] - last_element[i]) % bits;
  }

  // Last-row active bits: parity tail may have partial valid bits.
  const int active_bits_of_parity =
      (extra_bits_of_parity == 0) ? bits : extra_bits_of_parity;
  for (int j = 0; j < cols; ++j) {
    if (!occupied[rows - 1][j])
      continue;
    for (int k = 0; k < bits; ++k) {
      const int bit = (k + bits - element[rows - 1][j]) % bits;
      if (bit < active_bits_of_parity)
        last_row_active_bits[j][k] = true;
    }
  }

  setup_range_operational_matrix();

  std::printf("[H_MATRIX] H-matrix configured (userdata=%dB parity=%dB rows=%d "
              "cols=%d bits=%d matrix_sel=%d)\n",
              this->bytes_of_userdata, this->bytes_of_parity, rows, cols, bits,
              matrix_sel);
}

void h_matrix::print(void) {
  for (int i = 0; i < rows; ++i) {
    std::printf("%02d | ", row_weight[i]);
    for (int j = 0; j < cols; ++j) {
      if (occupied[i][j] || fade[i][j])
        std::printf("%03X ", element[i][j] & 0x3FF);
      else
        std::printf("... ");
    }
    std::printf("\n");
  }
}

void h_matrix::setup_range_operational_matrix(void) {
  for (int i = 0; i < LDPC_M; ++i) {
    for (int j = 0; j < LDPC_N; ++j)
      operational_h_matrix[i][j] = 0;
  }

  tile_ranges[0] = {0, 0};      // empty tile
  tile_ranges[1] = {LDPC_P, 0}; // full tile
  int tile_ranges_idx = 2;

  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      if (!occupied[i][j] && !fade[i][j]) {
        operational_h_matrix[i][j] = 0;
        continue;
      }

      range curr_range = {0, 0};
      bool head_found = false;
      for (int k = 0; k < bits; ++k) {
        bool do_not_use_this_bit = false;
        do_not_use_this_bit =
            do_not_use_this_bit ||
            ((extra_bits_of_parity > 0) && (j == (cols - rows)) &&
             (k >= extra_bits_of_parity));
        if (do_not_use_this_bit)
          continue;

        bool bit_active = false;
        if (occupied[i][j] && (i < (rows - 1)))
          bit_active = true;
        else if (occupied[i][j] && (i == (rows - 1)) &&
                 last_row_active_bits[j][k])
          bit_active = true;
        else if (fade[i][j] && !last_row_active_bits[j][k])
          bit_active = true;

        if (bit_active) {
          curr_range.active_count++;
          if (!head_found) {
            head_found = true;
            curr_range.offset = static_cast<unsigned short>(k);
          }
        } else if (head_found) {
          head_found = false;
        }
      }

      if (curr_range.active_count == 0) {
        operational_h_matrix[i][j] = 0;
      } else if (curr_range.active_count >= bits) {
        operational_h_matrix[i][j] = 1;
      } else if (tile_ranges_idx < 2 * LDPC_N) {
        operational_h_matrix[i][j] =
            static_cast<unsigned char>(tile_ranges_idx);
        tile_ranges[tile_ranges_idx] = curr_range;
        tile_ranges_idx++;
      } else {
        operational_h_matrix[i][j] = 1;
      }
    }
  }
}
