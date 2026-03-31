#pragma once

#include "codeword.h"

#include "h_matrix.h"

struct check_nodes {

  static constexpr int WORD_COUNT = 8;

  // 8 words * 64 bits = 512

  static constexpr int WORD_SIZE = 64;

  static constexpr int MAX_ROW_COUNT = LDPC_M;

  uint64_t rows[MAX_ROW_COUNT][WORD_COUNT];

  check_nodes();

  void clear(void);

  void print(void);
};

struct variable_nodes {

  static constexpr int WORD_COUNT = 8;

  // 8 words * 64 bits = 512

  static constexpr int WORD_SIZE = 64;

  static constexpr int MAX_COL_COUNT = LDPC_N;

  uint64_t flipped[MAX_COL_COUNT][WORD_COUNT];

  uint8_t likelihood[MAX_COL_COUNT][WORD_COUNT * WORD_SIZE];

  variable_nodes();

  void clear();

  void print(void);
};

struct s_512_bits {

  bool b[512];
};

struct s_likelihood_levels {

  int min;

  int max;

  int flip_thr;

  int weak;

  int strong;

  int level[4];
};

class codec_support {

public:
  static const uint64_t UNIT = 1ULL;

  static void rotate_left_512(const uint64_t in[8], uint64_t out[8],
                              unsigned r);

  static void mask_range_512(uint64_t out[8], unsigned start, unsigned len);
  static check_nodes config_check_nodes(const h_matrix &h_matrix_ref,
                                        const codeword &vn);
  static uint64_t check_node_weight(const check_nodes &cn);
};
