#ifndef _IBEX_RTL_ENGINE_H
#define _IBEX_RTL_ENGINE_H

#include "ldpc_codec.h"

struct ibex_family_lut {
  bool loaded;
  int min_m;
  int max_m;
  int k_min;
  int k_max;
  int wrap_base[LDPC_MAX_ROWS];
  int wrap_base_deltas[LDPC_MAX_ROWS][LDPC_MAX_ROWS * 4];
  int wrap_base_deltas_len[LDPC_MAX_ROWS];
  bool has_first_mask_shift[LDPC_MAX_ROWS + 1];
  int first_mask_shift[LDPC_MAX_ROWS + 1][4];
  char source_path[512];
};

struct ibex_rtl_view {
  int shift[LDPC_MAX_ROWS][LDPC_MAX_COLS];
  bool mask[LDPC_MAX_COLS][LDPC_MAX_CIRC_BITS];
  int wrap_base[LDPC_MAX_ROWS];
  int wrap_base_delta[LDPC_MAX_ROWS];
  int lut_first_mask_shift[4];
  int derived_first_mask_shift;
  bool has_lut_first_mask_shift;
  bool has_derived_first_mask_shift;
};

const ibex_family_lut &get_ibex_family_lut();
void rotate_cn_row(s_check_node_row *row, int bits, int shift);
int first_active_shift_for_row(const s_h_matrix &h_matrix,
                                   const ibex_rtl_view &view, int row);
void build_rtl_view(const s_h_matrix &h_matrix,
                        const ibex_family_lut &lut,
                        ibex_rtl_view *view);
void config_wrap(const s_h_matrix &h_matrix,
                     const ibex_family_lut &lut,
                     ibex_rtl_view *view);

#endif
