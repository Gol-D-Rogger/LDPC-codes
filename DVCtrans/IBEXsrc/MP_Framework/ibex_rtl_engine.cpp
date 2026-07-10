#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "ibex_rtl_engine.h"

static int parse_sv_scalar_after_equals(const std::string &line,
                                            int default_value = 0) {
  const size_t eq = line.find('=');
  if (eq == std::string::npos)
    return default_value;
  size_t pos = eq + 1;
  while ((pos < line.size()) &&
         !std::isdigit(static_cast<unsigned char>(line[pos])) &&
         (line[pos] != '-')) {
    pos++;
  }
  if (pos >= line.size())
    return default_value;
  return std::atoi(line.c_str() + pos);
}

static std::vector<int> parse_sv_decimal_list(const std::string &line) {
  std::vector<int> values;
  size_t pos = 0;
  while ((pos = line.find("'d", pos)) != std::string::npos) {
    pos += 2;
    size_t end = pos;
    while ((end < line.size()) &&
           std::isdigit(static_cast<unsigned char>(line[end]))) {
      end++;
    }
    if (end > pos)
      values.push_back(std::atoi(line.substr(pos, end - pos).c_str()));
    pos = end;
  }
  return values;
}

static int parse_suffix_digits(const std::string &line,
                                   const char *prefix,
                                   const char *suffix,
                                   int default_value = -1) {
  const size_t begin = line.find(prefix);
  if (begin == std::string::npos)
    return default_value;
  const size_t digits_start = begin + std::strlen(prefix);
  const size_t digits_end = line.find(suffix, digits_start);
  if ((digits_end == std::string::npos) || (digits_end <= digits_start))
    return default_value;
  return std::atoi(line.substr(digits_start, digits_end - digits_start).c_str());
}

static bool try_load_family_lut_file(const char *path,
                                         ibex_family_lut *lut) {
  if (!path || !*path || !lut)
    return false;

  printf("[LDPC] Try family LUT SVH: %s\n", path);
  std::ifstream file(path);
  if (!file.good()) {
    printf("[LDPC WARN] Cannot open family LUT SVH: %s\n", path);
    return false;
  }

  std::memset(lut, 0, sizeof(*lut));
  for (int row = 0; row < LDPC_MAX_ROWS; row++) {
    for (int idx = 0; idx < (LDPC_MAX_ROWS * 4); idx++)
      lut->wrap_base_deltas[row][idx] = 0;
    lut->wrap_base_deltas_len[row] = 0;
  }

  std::string line;
  while (std::getline(file, line)) {
    if ((line.find("LDPC_MATRIX_MIN_M") != std::string::npos) ||
        (line.find("ldpc_matrix_min_m") != std::string::npos)) {
      lut->min_m = parse_sv_scalar_after_equals(line, 0);
    } else if ((line.find("LDPC_MATRIX_MAX_M") != std::string::npos) ||
               (line.find("ldpc_matrix_max_m") != std::string::npos)) {
      lut->max_m = parse_sv_scalar_after_equals(line, 0);
    } else if ((line.find("LDPC_MATRIX_LUT_K_MIN") != std::string::npos) ||
               (line.find("ldpc_matrix_k_min") != std::string::npos)) {
      lut->k_min = parse_sv_scalar_after_equals(line, 0);
    } else if ((line.find("LDPC_MATRIX_LUT_K_MAX") != std::string::npos) ||
               (line.find("ldpc_matrix_k_max") != std::string::npos)) {
      lut->k_max = parse_sv_scalar_after_equals(line, 0);
    } else if (line.find("ldpc_matrix_wrap_base_tmp_adj_ord") !=
               std::string::npos) {
      const std::vector<int> values = parse_sv_decimal_list(line);
      for (size_t i = 0; (i < values.size()) && (i < LDPC_MAX_ROWS); i++)
        lut->wrap_base[i] = values[i];
    } else if (line.find("ldpc_matrix_first_mask_shift") !=
               std::string::npos) {
      const int m = parse_suffix_digits(line, "ldpc_matrix_first_mask_shift",
                                            " = {");
      const std::vector<int> values = parse_sv_decimal_list(line);
      if ((m >= 0) && (m <= LDPC_MAX_ROWS)) {
        lut->has_first_mask_shift[m] = (values.size() >= 4);
        for (size_t i = 0; (i < values.size()) && (i < 4); i++)
          lut->first_mask_shift[m][i] = values[i];
      }
    } else if (line.find("ldpc_matrix_wrap_base_deltas") !=
               std::string::npos) {
      const int row = parse_suffix_digits(
          line, "ldpc_matrix_wrap_base_deltas", "_adj_ord");
      const std::vector<int> values = parse_sv_decimal_list(line);
      if ((row >= 0) && (row < LDPC_MAX_ROWS)) {
        lut->wrap_base_deltas_len[row] =
            std::min<int>(values.size(), LDPC_MAX_ROWS * 4);
        for (int i = 0; i < lut->wrap_base_deltas_len[row]; i++)
          lut->wrap_base_deltas[row][i] = values[i];
      }
    }
  }

  const bool ok = (lut->min_m > 0) && (lut->max_m >= lut->min_m) &&
                  (lut->k_max >= lut->k_min);
  if (!ok) {
    printf("[LDPC WARN] Invalid family LUT SVH format: %s (min_m=%d max_m=%d k_min=%d k_max=%d)\n",
           path, lut->min_m, lut->max_m, lut->k_min, lut->k_max);
    return false;
  }

  snprintf(lut->source_path, sizeof(lut->source_path), "%s", path);
  lut->loaded = true;
  printf("[LDPC] Loaded IBEX family LUT from %s (M=%d..%d, K=%d..%d)\n",
         path, lut->min_m, lut->max_m, lut->k_min, lut->k_max);
  return true;
}

const ibex_family_lut &get_ibex_family_lut() {
  static bool initialized = false;
  static ibex_family_lut lut;
  if (initialized)
    return lut;

  initialized = true;
  std::memset(&lut, 0, sizeof(lut));

  const char *default_paths[] = {
      "IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
      "../IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
      "../../../../IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
      "../../../../../IBEX/ibex_matrix/family_lut_output/ldpc_matrix_lut_fixed_wrap_base.svh",
  };
  const size_t default_paths_count =
      sizeof(default_paths) / sizeof(default_paths[0]);
  for (size_t i = 0; i < default_paths_count; ++i) {
    const char *path = default_paths[i];
    if (try_load_family_lut_file(path, &lut))
      return lut;
  }

  printf("[LDPC WARN] Failed to load IBEX family LUT sidecar; RTL-CN path will "
         "fall back to h_matrix wrap_base/wrap_num_deltas.\n");
  return lut;
}

static int lookup_wrap_base_delta(const ibex_family_lut &lut, int row,
                                      int m, int k) {
  if (!lut.loaded || (row < 0) || (row >= LDPC_MAX_ROWS) || (m < lut.min_m) ||
      (m > lut.max_m) || (k < lut.k_min) || (k > lut.k_max))
    return -1;
  const int num_k = lut.k_max - lut.k_min + 1;
  const int idx = (m - lut.min_m) * num_k + (k - lut.k_min);
  if ((idx < 0) || (idx >= lut.wrap_base_deltas_len[row]))
    return -1;
  return lut.wrap_base_deltas[row][idx];
}

static bool lookup_first_mask_shift(const ibex_family_lut &lut, int m,
                                        int k, int *value) {
  if (!value || !lut.loaded || (m < lut.min_m) || (m > lut.max_m) ||
      (k < lut.k_min) || (k > lut.k_max) || !lut.has_first_mask_shift[m])
    return false;
  *value = lut.first_mask_shift[m][k - lut.k_min];
  return true;
}

void rotate_cn_row(s_check_node_row *row, int bits, int shift) {
  if (!row || (bits <= 0))
    return;
  shift %= bits;
  if (shift < 0)
    shift += bits;
  if (shift == 0)
    return;

  bool rotated[LDPC_MAX_CIRC_BITS];
  for (int k = 0; k < bits; k++)
    rotated[k] = row->b[(k + bits - shift) % bits];
  for (int k = 0; k < bits; k++)
    row->b[k] = rotated[k];
}

int first_active_shift_for_row(const s_h_matrix &h_matrix,
                                   const ibex_rtl_view &view, int row) {
  for (int col = 0; col < h_matrix.cols; col++) {
    if ((row >= 0) && (row < h_matrix.rows) &&
        (h_matrix.occupied[row][col] || h_matrix.fade[row][col]) &&
        (view.shift[row][col] >= 0)) {
      return view.shift[row][col];
    }
  }
  return 0;
}

static const int WBD_SEARCH_MAX = 512;

static int search_wrap_base_delta(int last_element, int wrap_base,
                                  int delta, int bits) {
  if (delta == 0)
    return 0;
  for (int k = 0; k < WBD_SEARCH_MAX; k++) {
    if ((last_element + wrap_base + k * delta) % bits == wrap_base)
      return k;
  }
  return -1;
}

static bool env_enabled(const char *name) {
  const char *value = std::getenv(name);
  return value && *value && (std::strcmp(value, "0") != 0);
}

void build_rtl_view(const s_h_matrix &h_matrix,
                        const ibex_family_lut &lut,
                        ibex_rtl_view *view) {
  if (!view)
    return;

  std::memset(view, 0, sizeof(*view));
  for (int i = 0; i < LDPC_MAX_ROWS; i++) {
    for (int j = 0; j < LDPC_MAX_COLS; j++)
      view->shift[i][j] = -1;
  }

  const int m = h_matrix.rows;
  const int k_payload = h_matrix.cols - h_matrix.rows;
  const int last_row = h_matrix.rows - 1;
  const bool use_lut_shifts = lut.loaded;
  const bool force_rebuilt_mask = env_enabled("IBEX_RTL_FORCE_REBUILT_MASK");
  int first_element_last_row = -1;
  int first_parity_shift_last_row = -1;

  for (int row = 0; row < h_matrix.rows; row++) {
    const int wb = (use_lut_shifts && (row < LDPC_MAX_ROWS))
                       ? lut.wrap_base[row]
                       : h_matrix.wrap_base[row];
    view->wrap_base[row] = wb;

    int next_shift = (wb + h_matrix.delta[row]) % h_matrix.bits;
    int last_active_shift = -1;
    for (int col = 0; col < h_matrix.cols; col++) {
      if (h_matrix.occupied[row][col] || h_matrix.fade[row][col]) {
        const int rebuilt_shift = next_shift;
        view->shift[row][col] = use_lut_shifts
                                    ? rebuilt_shift
                                    : ((h_matrix.element[row][col] >= 0)
                                           ? h_matrix.element[row][col]
                                           : rebuilt_shift);
        last_active_shift = view->shift[row][col];
        if ((row == last_row) && (first_element_last_row < 0))
          first_element_last_row = view->shift[row][col];
        if ((row == last_row) && (col >= (h_matrix.cols - h_matrix.rows)) &&
            h_matrix.occupied[row][col] && (first_parity_shift_last_row < 0)) {
          first_parity_shift_last_row = view->shift[row][col];
        }
        next_shift = (next_shift + h_matrix.delta[row]) % h_matrix.bits;
      }
    }

    const int lut_wbd = use_lut_shifts
                            ? lookup_wrap_base_delta(lut, row, m, k_payload)
                            : -1;
    if (lut_wbd >= 0) {
      view->wrap_base_delta[row] = lut_wbd;
    } else {
      const int le = (last_active_shift >= 0) ? last_active_shift
                                              : h_matrix.last_element[row];
      const int searched = search_wrap_base_delta(le, wb, h_matrix.delta[row], h_matrix.bits);
      view->wrap_base_delta[row] =
          (searched >= 0) ? searched : h_matrix.wrap_num_deltas[row];
    }
  }

  if (h_matrix.extra_bits_of_parity > 0) {
    for (int col = 0; col < h_matrix.cols; col++) {
      if (!h_matrix.occupied[last_row][col] || (view->shift[last_row][col] < 0))
        continue;
      const int start_shift = view->shift[last_row][col];
      for (int bit = 0; bit < h_matrix.bits; bit++) {
        const int offset = (bit + h_matrix.bits - start_shift) % h_matrix.bits;
        if (offset < h_matrix.extra_bits_of_parity)
          view->mask[col][bit] = true;
      }
    }

    if (!force_rebuilt_mask) {
      for (int col = 0; col < h_matrix.cols; col++) {
        bool use_loaded_mask = false;
        for (int bit = 0; bit < h_matrix.bits; bit++) {
          if (h_matrix.mask[col][bit]) {
            use_loaded_mask = true;
            break;
          }
        }
        if (use_loaded_mask) {
          for (int bit = 0; bit < h_matrix.bits; bit++)
            view->mask[col][bit] = h_matrix.mask[col][bit];
        }
      }
    }
  }

  view->has_derived_first_mask_shift =
      (first_element_last_row >= 0) && (first_parity_shift_last_row >= 0);
  if (view->has_derived_first_mask_shift) {
    view->derived_first_mask_shift =
        (h_matrix.bits + first_element_last_row - first_parity_shift_last_row) %
        h_matrix.bits;
  }

  view->has_lut_first_mask_shift = false;
  for (int idx = 0; idx < 4; idx++) {
    int first_mask_shift = 0;
    if (lookup_first_mask_shift(lut, m, lut.k_min + idx,
                                    &first_mask_shift)) {
      view->lut_first_mask_shift[idx] = first_mask_shift;
      view->has_lut_first_mask_shift = true;
    }
  }

  if (view->has_lut_first_mask_shift && view->has_derived_first_mask_shift) {
    const int lut_idx = std::max(0, std::min(3, k_payload - lut.k_min));
    if (view->lut_first_mask_shift[lut_idx] != view->derived_first_mask_shift) {
      printf("[LDPC WARN] RTL-CN rebuilt first_mask_shift mismatch: M=%d K=%d "
             "derived=%d lut=%d\n",
             m, k_payload, view->derived_first_mask_shift,
             view->lut_first_mask_shift[lut_idx]);
    }
  }
}

void config_wrap(const s_h_matrix &h_matrix,
                     const ibex_family_lut &lut,
                     ibex_rtl_view *view) {
  (void)lut;
  if (!view)
    return;

  // Temporary alignment mode:
  // Do not rebuild matrix shifts from RTL/LUT chain.
  // Use the same matrix path as baseline decoder (`h_matrix` from ldpc_config)
  // so only decoder behavior differs in A/B comparison.
  std::memset(view, 0, sizeof(*view));
  for (int i = 0; i < LDPC_MAX_ROWS; i++) {
    for (int j = 0; j < LDPC_MAX_COLS; j++)
      view->shift[i][j] = -1;
  }

  for (int row = 0; row < h_matrix.rows; row++) {
    view->wrap_base[row] = h_matrix.wrap_base[row];
    view->wrap_base_delta[row] = h_matrix.wrap_num_deltas[row];
    for (int col = 0; col < h_matrix.cols; col++) {
      if (h_matrix.occupied[row][col] || h_matrix.fade[row][col]) {
        view->shift[row][col] = h_matrix.element[row][col];
      }
    }
  }

  for (int col = 0; col < h_matrix.cols; col++) {
    for (int bit = 0; bit < h_matrix.bits; bit++) {
      view->mask[col][bit] = h_matrix.mask[col][bit];
    }
  }

  view->has_lut_first_mask_shift = false;
  view->has_derived_first_mask_shift = false;
}

static void mkdir_if_missing_rtl(const char *path) {
  if (!path || !*path)
    return;
  struct stat st;
  if (stat(path, &st) == 0) {
    return;
  }
  if (mkdir(path, 0777) != 0 && errno != EEXIST)
    printf("[LDPC WARN] mkdir failed: %s (errno=%d)\n", path, errno);
}

static FILE *open_bf_ibex_sw_delta_dump_file_rtl() {
  mkdir_if_missing_rtl("./output");
  return fopen("./output/bf_ibex_sw_delta.dat", "w");
}

static FILE *open_bf_ibex_row_sw_dump_file_rtl() {
  mkdir_if_missing_rtl("./output");
  return fopen("./output/bf_ibex_row_sw.dat", "w");
}

static void log_bf_ibex_sw_rtl(FILE *fp, int iteration, const char *phase,
                               int col, int syndrome_weight) {
  if (!fp || !phase)
    return;
  fprintf(fp, "%8d %12s %6d %16d\n", iteration, phase, col, syndrome_weight);
}

static void log_bf_ibex_row_sw_header_rtl(FILE *fp, int rows) {
  if (!fp)
    return;
  fprintf(fp, "%8s %12s %6s %16s", "iter", "phase", "col", "syndrome_weight");
  for (int row = 0; row < rows; row++)
    fprintf(fp, " row%02d_sw", row);
  fputc('\n', fp);
}

static void log_bf_ibex_row_sw_rtl(FILE *fp, int iteration, const char *phase,
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
  for (int row = 0; row < h_matrix.rows; row++)
    fprintf(fp, " %8d", row_sw[row]);
  fputc('\n', fp);
}

#ifdef _LDPC_DBG_DUMP
static FILE *open_bf_ibex_rtl_cn_trace_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  mkdir_if_missing_rtl("./output");
  snprintf(path, sizeof(path), "./output/ldpc_dbg_rtl_cn_trace_cw%04u.txt",
           cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open RTL-CN trace dump: %s\n", path);
  return fp;
}

static FILE *open_bf_ibex_rtl_cn_like_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  mkdir_if_missing_rtl("./output");
  snprintf(path, sizeof(path),
           "./output/ldpc_dbg_bf_ibex_rtl_cn_like_cw%04u.txt", cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open RTL-CN likelihood dump: %s\n", path);
  return fp;
}

static FILE *open_bf_ibex_rtl_cn_post_ctrl_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  mkdir_if_missing_rtl("./output");
  snprintf(path, sizeof(path),
           "./output/ldpc_dbg_bf_ibex_rtl_cn_post_ctrl_cw%04u.txt",
           cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open RTL-CN post-ctrl dump: %s\n", path);
  return fp;
}

static FILE *open_bf_ibex_rtl_cn_post_eff_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  mkdir_if_missing_rtl("./output");
  snprintf(path, sizeof(path),
           "./output/ldpc_dbg_bf_ibex_rtl_cn_post_eff_cw%04u.txt", cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open RTL-CN post-eff dump: %s\n", path);
  return fp;
}

static FILE *open_bf_ibex_rtl_cn_post_trace_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  mkdir_if_missing_rtl("./output");
  snprintf(path, sizeof(path),
           "./output/ldpc_dbg_bf_ibex_rtl_cn_post_trace_cw%04u.txt",
           cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open RTL-CN post-trace dump: %s\n", path);
  return fp;
}

static FILE *open_bf_ibex_rtl_cn_unsat_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  mkdir_if_missing_rtl("./output");
  snprintf(path, sizeof(path),
           "./output/ldpc_dbg_bf_ibex_rtl_cn_unsat_cw%04u.txt", cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open RTL-CN unsat dump: %s\n", path);
  return fp;
}

static FILE *open_bf_ibex_rtl_cn_prng512_dump_file() {
  static unsigned int cw_seq = 0;
  char path[128];
  mkdir_if_missing_rtl("./output");
  snprintf(path, sizeof(path),
           "./output/ldpc_dbg_bf_ibex_rtl_cn_prng512_cw%04u.txt", cw_seq++);
  FILE *fp = fopen(path, "w");
  if (!fp)
    printf("[LDPC WARN] failed to open RTL-CN prng512 dump: %s\n", path);
  return fp;
}

static void dump_bf_ibex_rtl_cn_col_likelihood(FILE *fp,
                                               const s_variable_nodes &vn,
                                               int itr_cnt, int col_cnt,
                                               int cir_bits) {
  if (!fp || (col_cnt < 0) || (cir_bits <= 0))
    return;

  fprintf(fp, "itr_cnt:%03d col_cnt:%03d:", itr_cnt, col_cnt);
  for (int k = cir_bits - 1; k >= 0; k--)
    fprintf(fp, "%X", (unsigned int)vn.c[col_cnt].b[k].likelihood);
  fputc('\n', fp);
}

static void dump_bf_ibex_rtl_cn_col_unsat(FILE *fp,
                                          const unsigned char *unsat_cnt,
                                          int itr_cnt, int col_cnt,
                                          int cir_bits) {
  if (!fp || !unsat_cnt || (col_cnt < 0) || (cir_bits <= 0))
    return;

  fprintf(fp, "itr_cnt:%03d col_cnt:%03d:", itr_cnt, col_cnt);
  for (int k = cir_bits - 1; k >= 0; k--)
    fprintf(fp, "%1X", (unsigned int)unsat_cnt[k]);
  fputc('\n', fp);
}

static void log_bf_ibex_rtl_cn_post_ctrl(
    FILE *fp, int iteration, int col, int post_iter, int syndrome_weight,
    int syndrome_weight_delayed, int prev_sw, bool pushing, bool post_trigger,
    bool post_trigger2, int likelihood_min, int flip_thr, int thr_post,
    int thr_qc, int post_ratio) {
  if (!fp)
    return;
  fprintf(fp,
          "%8d %6d %9d %8d %8d %8d %8d %6d %6d %10d %8d %8d %8d %10d\n",
          iteration, col, post_iter, syndrome_weight, syndrome_weight_delayed,
          prev_sw, (int)pushing, (int)post_trigger, (int)post_trigger2,
          likelihood_min, flip_thr, thr_post, thr_qc, post_ratio);
}

static void log_bf_ibex_rtl_cn_post_eff(FILE *fp, int iteration, int col,
                                        int post1_cnt, int post2_cnt,
                                        int aggressive_cnt,
                                        int flipchg_cnt,
                                        int syndrome_weight_after_update) {
  if (!fp)
    return;
  fprintf(fp, "%8d %6d %8d %8d %10d %10d %16d\n", iteration, col, post1_cnt,
          post2_cnt, aggressive_cnt, flipchg_cnt,
          syndrome_weight_after_update);
}

static void dump_bf_ibex_rtl_cn_post_trace(
    FILE *fp, const s_variable_nodes &vn, int iteration, int col,
    int post_iteration, int syndrome_weight, int syndrome_weight_delayed,
    bool post_trigger, bool post_trigger2, int likelihood_min, int flip_thr,
    int thr_post, int thr_qc, int post_ratio, int post1_cnt, int post2_cnt,
    int aggressive_cnt, int flipchg_cnt, int syndrome_weight_after_update,
    int cir_bits) {
  if (!fp || (col < 0) || (cir_bits <= 0))
    return;

  fprintf(fp,
          "================================================================\n");
  fprintf(fp, "[RTL_CN POST] ITER=%04d COL=%02d POST_ITER=%d\n", iteration, col,
          post_iteration);
  fprintf(fp,
          "[CTRL] sw=%d sw_dly=%d trig=%d trig2=%d like_min=%d flip_thr=%d "
          "thr_post=%d thr_qc=%d post_ratio=%d\n",
          syndrome_weight, syndrome_weight_delayed, (int)post_trigger,
          (int)post_trigger2, likelihood_min, flip_thr, thr_post, thr_qc,
          post_ratio);
  fprintf(fp,
          "[EFF] post1=%d post2=%d aggr=%d flipchg=%d sw_after_update=%d\n",
          post1_cnt, post2_cnt, aggressive_cnt, flipchg_cnt,
          syndrome_weight_after_update);
  fprintf(fp, "[LIKE]\n");
  dump_bf_ibex_rtl_cn_col_likelihood(fp, vn, iteration, col, cir_bits);
  fputc('\n', fp);
  fflush(fp);
}

static char dbg_hex_digit(unsigned int val) {
  static const char kHex[] = "0123456789ABCDEF";
  return kHex[val & 0xF];
}

static void dump_bf_ibex_rtl_cn_prng512(FILE *fp, const s_512_bits &prng,
                                        int itr_cnt, int col_cnt) {
  if (!fp || (col_cnt < 0))
    return;

  fprintf(fp, "itr_cnt:%03d col_cnt:%03d:", itr_cnt, col_cnt);
  for (int k = 511; k >= 0; k -= 4) {
    unsigned int nibble = 0;
    nibble |= prng.b[k] ? 8u : 0u;
    nibble |= prng.b[k - 1] ? 4u : 0u;
    nibble |= prng.b[k - 2] ? 2u : 0u;
    nibble |= prng.b[k - 3] ? 1u : 0u;
    fputc(dbg_hex_digit(nibble), fp);
  }
  fputc('\n', fp);
}

static void dbg_fprint_bool_bits_hex(FILE *fp, const bool *bits, int nbits,
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
      fputc(dbg_hex_digit(nibble), fp);
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
    fputc(dbg_hex_digit(nibble), fp);
    hex_count++;
    if ((hex_count < total_hex) && ((hex_count % hex_wrap) == 0)) {
      fputc('\n', fp);
      if (prefix_cont)
        fputs(prefix_cont, fp);
    }
  }
}

static void dbg_fprint_shifted_bool_bits_hex(FILE *fp, const bool *bits,
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
      fputc(dbg_hex_digit(nibble), fp);
      hex_count++;
      if ((hex_count < total_hex) && ((hex_count % 1024) == 0))
        fputc('\n', fp);
      nibble = 0;
      nibble_bits = 0;
    }
  }

  if (nibble_bits != 0) {
    nibble <<= (4 - nibble_bits);
    fputc(dbg_hex_digit(nibble), fp);
    hex_count++;
    if ((hex_count < total_hex) && ((hex_count % 1024) == 0))
      fputc('\n', fp);
  }
}

static void dump_cn_synd_file_rtl(FILE *fp, const s_h_matrix &h_matrix,
                                  const s_check_nodes &cn) {
  if (!fp)
    return;

  const int syndrome_dump_bits =
      (h_matrix.bits >= 512) ? 512 : h_matrix.bits;
  const int chunk_bits = 128;
  for (int rr = 0; rr < h_matrix.rows; rr++) {
    for (int local_lo = 0; local_lo < syndrome_dump_bits;
         local_lo += chunk_bits) {
      const int len = std::min(chunk_bits, syndrome_dump_bits - local_lo);
      char prefix[32];
      snprintf(prefix, sizeof(prefix), "row%02d ", rr);
      dbg_fprint_bool_bits_hex(fp, cn.r[rr].b + local_lo, len, 1024, prefix,
                               NULL);
      fputc('\n', fp);
    }
  }
}

static void dump_cn_synd_canon_file_rtl(FILE *fp, const s_h_matrix &h_matrix,
                                        const s_check_nodes &cn,
                                        const int *accum_rot) {
  if (!fp || !accum_rot)
    return;

  const int syndrome_dump_bits =
      (h_matrix.bits >= 512) ? 512 : h_matrix.bits;
  const int chunk_bits = 128;
  for (int rr = 0; rr < h_matrix.rows; rr++) {
    const int canon_shift =
        ((accum_rot[rr] % h_matrix.bits) + h_matrix.bits) % h_matrix.bits;
    for (int local_lo = 0; local_lo < syndrome_dump_bits;
         local_lo += chunk_bits) {
      const int len = std::min(chunk_bits, syndrome_dump_bits - local_lo);
      char prefix[32];
      snprintf(prefix, sizeof(prefix), "row%02d a=%03d ", rr, canon_shift);
      fputs(prefix, fp);
      dbg_fprint_shifted_bool_bits_hex(fp, cn.r[rr].b, h_matrix.bits, local_lo,
                                       len, canon_shift, false);
      fputc('\n', fp);
    }
  }
}

static void dump_cn_synd_shifted_col_file_rtl(FILE *fp,
                                              const s_h_matrix &h_matrix,
                                              const ibex_rtl_view &rtl_view,
                                              const s_check_nodes &cn, int col,
                                              bool right_shift) {
  if (!fp)
    return;

  const int syndrome_dump_bits =
      (h_matrix.bits >= 512) ? 512 : h_matrix.bits;
  const int chunk_bits = 128;
  fprintf(fp, "col%02d\n", col);
  for (int rr = 0; rr < h_matrix.rows; rr++) {
    const int shift = rtl_view.shift[rr][col];
    for (int local_lo = 0; local_lo < syndrome_dump_bits;
         local_lo += chunk_bits) {
      const int len = std::min(chunk_bits, syndrome_dump_bits - local_lo);
      char prefix[32];
      snprintf(prefix, sizeof(prefix), "row%02d s=%03d ", rr, shift);
      fputs(prefix, fp);
      if (shift >= 0)
        dbg_fprint_shifted_bool_bits_hex(fp, cn.r[rr].b, h_matrix.bits,
                                         local_lo, len, shift, right_shift);
      else
        dbg_fprint_bool_bits_hex(fp, cn.r[rr].b + local_lo, len, 1024, NULL,
                                 NULL);
      fputc('\n', fp);
    }
  }
  fputc('\n', fp);
}

static void dump_cn_synd_shifted_file_rtl(FILE *fp,
                                          const s_h_matrix &h_matrix,
                                          const ibex_rtl_view &rtl_view,
                                          const s_check_nodes &cn,
                                          bool right_shift) {
  if (!fp)
    return;

  const int syndrome_dump_bits =
      (h_matrix.bits >= 512) ? 512 : h_matrix.bits;
  const int chunk_bits = 128;
  for (int col = 0; col < h_matrix.cols; col++) {
    fprintf(fp, "col%02d\n", col);
    for (int rr = 0; rr < h_matrix.rows; rr++) {
      const int shift = rtl_view.shift[rr][col];
      for (int local_lo = 0; local_lo < syndrome_dump_bits;
           local_lo += chunk_bits) {
        const int len = std::min(chunk_bits, syndrome_dump_bits - local_lo);
        char prefix[32];
        snprintf(prefix, sizeof(prefix), "row%02d s=%03d ", rr, shift);
        fputs(prefix, fp);
        if (shift >= 0)
          dbg_fprint_shifted_bool_bits_hex(fp, cn.r[rr].b, h_matrix.bits,
                                           local_lo, len, shift, right_shift);
        else
          dbg_fprint_bool_bits_hex(fp, cn.r[rr].b + local_lo, len, 1024, NULL,
                                   NULL);
        fputc('\n', fp);
      }
    }
    fputc('\n', fp);
  }
}

static void dump_cn_synd_snapshot_rtl(FILE *fp, const s_h_matrix &h_matrix,
                                      const ibex_rtl_view &rtl_view,
                                      const s_check_nodes &cn,
                                      const int *accum_rot,
                                      int iteration, int col,
                                      const char *phase, int syndrome_weight,
                                      bool include_shifted) {
  if (!fp || !phase)
    return;

  fprintf(fp,
          "================================================================\n");
  if (col >= 0)
    fprintf(fp, "[RTL_CN] ITER=%04d COL=%02d PHASE=%s SW=%d\n", iteration, col,
            phase, syndrome_weight);
  else
    fprintf(fp, "[RTL_CN] ITER=%04d PHASE=%s SW=%d\n", iteration, phase,
            syndrome_weight);

  fprintf(fp, "[RAW]\n");
  dump_cn_synd_file_rtl(fp, h_matrix, cn);

  if (accum_rot) {
    fprintf(fp, "[CANON_BY_ACCUM_ROT]\n");
    dump_cn_synd_canon_file_rtl(fp, h_matrix, cn, accum_rot);
  }

  if (include_shifted && (col >= 0)) {
    fprintf(fp, "[SEM_CUR_COL_RSHIFT COL=%02d]\n", col);
    dump_cn_synd_shifted_col_file_rtl(fp, h_matrix, rtl_view, cn, col, true);
    fprintf(fp, "[SEM_CUR_COL_LSHIFT COL=%02d]\n", col);
    dump_cn_synd_shifted_col_file_rtl(fp, h_matrix, rtl_view, cn, col, false);
  }

  fputc('\n', fp);
  fflush(fp);
}
#endif

void ldpc_packet::ldpc_dec_bf_ibex_rtl_cn(
    s_ldpc_decoder_input ldpc_decoder_input,
    s_ldpc_decoder_parameters ldpc_decoder_parameters, s_h_matrix h_matrix) {
  const bool experimental_rtl_cn =
      (getenv("IBEX_EXPERIMENTAL_RTL_CN") != NULL) &&
      (std::strcmp(getenv("IBEX_EXPERIMENTAL_RTL_CN"), "1") == 0);
  if (!experimental_rtl_cn) {
    ldpc_dec_bf_ibex(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);
    return;
  }

  const bool diag_ab = env_enabled("IBEX_RTL_CN_DIAG");
  const int VERBOSITY = 0;
  const int MAX_ERROR_COUNT = 4095;
  int iteration = 0;
  int clock_cycles = 0;
  int syndrome_weight = 0;
  int syndrome_weight_delayed = 0;
  bool finished = false;
  bool diag_mismatch_found = false;

  s_hard_codeword hard_codeword;
  s_variable_nodes vn;
  s_check_nodes cn;
  // s_check_nodes cn_delay1;
  // s_check_nodes cn_delay2;
  s_check_nodes cn_ref;
  s_likelihood_levels likelihood_levels;
  s_256_bits prng_256;
  s_512_bits prng_512;
  int prng_init[32];
  int syndrome_weight_r[5] = {0, 0, 0, 0, 0};
  int accum_rot[LDPC_MAX_ROWS];

  FILE *sw_delta_fp = open_bf_ibex_sw_delta_dump_file_rtl();
  FILE *row_sw_fp = open_bf_ibex_row_sw_dump_file_rtl();
  if (sw_delta_fp)
    fprintf(sw_delta_fp, "%8s %12s %6s %16s\n", "iter", "phase", "col",
            "syndrome_weight");
  log_bf_ibex_row_sw_header_rtl(row_sw_fp, h_matrix.rows);
#ifdef _LDPC_DBG_DUMP
  FILE *dbg_trace_fp = open_bf_ibex_rtl_cn_trace_dump_file();
  FILE *dbg_like_fp = open_bf_ibex_rtl_cn_like_dump_file();
  FILE *dbg_post_ctrl_fp = open_bf_ibex_rtl_cn_post_ctrl_dump_file();
  FILE *dbg_post_eff_fp = open_bf_ibex_rtl_cn_post_eff_dump_file();
  FILE *dbg_post_trace_fp = open_bf_ibex_rtl_cn_post_trace_dump_file();
  FILE *dbg_unsat_fp = open_bf_ibex_rtl_cn_unsat_dump_file();
  FILE *dbg_prng512_fp = open_bf_ibex_rtl_cn_prng512_dump_file();
  if (dbg_trace_fp)
    fprintf(dbg_trace_fp,
            "[RTL_CN TRACE] per-column and per-iteration syndrome dump\n");
  if (dbg_post_ctrl_fp)
    fprintf(dbg_post_ctrl_fp,
            "%8s %6s %9s %8s %8s %8s %8s %6s %6s %10s %8s %8s %8s %10s\n",
            "iter", "col", "post_iter", "sw", "sw_dly", "prev_sw",
            "pushing", "trig", "trig2", "like_min", "flip_thr",
            "thr_post", "thr_qc", "post_ratio");
  if (dbg_post_eff_fp)
    fprintf(dbg_post_eff_fp,
            "%8s %6s %8s %8s %10s %10s %16s\n", "iter", "col", "post1",
            "post2", "aggr", "flipchg", "sw_after_update");
  if (dbg_post_trace_fp)
    fprintf(dbg_post_trace_fp,
            "[RTL_CN POST TRACE] post-only per-column control/result/likelihood\n");
  if (dbg_unsat_fp)
    fprintf(dbg_unsat_fp,
            "[RTL_CN UNSAT] per-column per-bit unsat count, 1-digit hex per bit (MSB->LSB)\n");
  if (dbg_prng512_fp)
    fprintf(dbg_prng512_fp,
            "[RTL_CN PRNG512] per-column PRNG state after init/advance, 1-digit hex per 4 bits (MSB->LSB)\n");
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

  std::memset(&hard_codeword, 0, sizeof(hard_codeword));
  std::memset(&vn, 0, sizeof(vn));
  for (int j = 0; j < h_matrix.cols; j++) {
    for (int k = 0; k < h_matrix.bits; k++) {
      hard_codeword.c[j].b[k] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard;
      vn.c[j].b[k].bit_hard = hard_codeword.c[j].b[k];
      vn.c[j].b[k].flipped = 0;
    }
  }

  const ibex_family_lut &family_lut = get_ibex_family_lut();
  ibex_rtl_view rtl_view;
  build_rtl_view(h_matrix, family_lut, &rtl_view);

  if (family_lut.loaded) {
    for (int i = 0; i < h_matrix.rows; i++)
      for (int j = 0; j < h_matrix.cols; j++)
        if (rtl_view.shift[i][j] >= 0)
          h_matrix.element[i][j] = rtl_view.shift[i][j];
    if (h_matrix.extra_bits_of_parity > 0) {
      for (int j = 0; j < h_matrix.cols; j++)
        for (int k = 0; k < h_matrix.bits; k++)
          h_matrix.mask[j][k] = rtl_view.mask[j][k];
    }
  }

  cn = f_check_nodes(h_matrix, hard_codeword);
  if (diag_ab)
    cn_ref = cn;
  for (int i = 0; i < h_matrix.rows; i++) {
    const int first_shift =
        first_active_shift_for_row(h_matrix, rtl_view, i);
    rotate_cn_row(&cn.r[i], h_matrix.bits, first_shift);
    accum_rot[i] = first_shift;
  }
  // cn_delay1 = cn;
  // cn_delay2 = cn;
  if (diag_ab) {
    printf("[DIAG] Initial rotation applied. Verifying syndrome equivalence...\n");
    for (int i = 0; i < h_matrix.rows; i++) {
      for (int k = 0; k < h_matrix.bits; k++) {
        const int m = (k + h_matrix.bits - accum_rot[i]) % h_matrix.bits;
        if (cn.r[i].b[k] != cn_ref.r[i].b[m]) {
          printf("[DIAG] INIT ROT MISMATCH row=%d k=%d: rtl_cn=%d ref_at_m=%d (m=%d, rot=%d)\n",
                 i, k, (int)cn.r[i].b[k], (int)cn_ref.r[i].b[m], m, accum_rot[i]);
          diag_mismatch_found = true;
          break;
        }
      }
      if (diag_mismatch_found) break;
    }
    if (!diag_mismatch_found)
      printf("[DIAG] Initial rotation: OK (all rows match)\n");
  }
  syndrome_weight = f_check_node_weight(h_matrix, cn);
  clock_cycles = (2 * h_matrix.cols) + 1;
  if (syndrome_weight == 0)
    finished = true;
  ldpc_decoder_output.syndrome_weight_before = syndrome_weight;
  syndrome_weight_delayed = syndrome_weight;
  ldpc_decoder_output.early_termination = 0;
  ldpc_decoder_output.col_cnt = -1;

  if (ldpc_decoder_parameters.early_terminate_dis == 0) {
    int early_term_thr = 0;
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
      finished = true;
      ldpc_decoder_output.early_termination = 1;
    }
  }

  if (ldpc_decoder_input.syndrome_cal_only)
    finished = true;

  likelihood_levels = f_likelihood_levels(ldpc_decoder_input.nand_strobes,
                                          ldpc_decoder_parameters,
                                          syndrome_weight, h_matrix.rows);

  // RTL treats the initial syndrome calculation as a standalone iteration 0.
  // Prime the delayed syndrome-weight pipeline with that value so the first
  // column-updating pass starts at iteration 1 with a settled history.
  if (!finished) {
    for (int idx = 0; idx < 5; idx++)
      syndrome_weight_r[idx] = syndrome_weight;
    // RTL spends the first two clocks of iteration-1 col0 while still in
    // iteration 0, then the remaining two clocks at the start of iteration 1.
    for (int sw_clk = 0; sw_clk < 2; sw_clk++) {
      syndrome_weight_r[4] = syndrome_weight_r[3];
      syndrome_weight_r[3] = syndrome_weight_r[2];
      syndrome_weight_r[2] = syndrome_weight_r[1];
      syndrome_weight_r[1] = syndrome_weight_r[0];
      syndrome_weight_r[0] = syndrome_weight;
    }
    clock_cycles += 2;
    iteration = 1;
  }
  const int post_start_iteration =
      (ldpc_decoder_input.post_iteration < 1) ? 1
                                              : ldpc_decoder_input.post_iteration;
  // int sw_delay_clock_idx = 0;

  for (int j = 0; j < h_matrix.cols; j++) {
    for (int k = 0; k < h_matrix.bits; k++) {
      const bool pad_userdata =
          (j == (h_matrix.cols - h_matrix.rows - 1)) &&
          (k >= (h_matrix.bits - 8 * h_matrix.unused_bytes_of_userdata));
      const bool pad_parity =
          (j == (h_matrix.cols - h_matrix.rows)) &&
          (k >= (h_matrix.bits - 8 * h_matrix.unused_bytes_of_parity));
      if (pad_userdata || pad_parity)
        continue;

      const bool soft0 =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable;
      const bool soft1 =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2;
      if (!soft1 && !soft0)
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[0]];
      if (!soft1 && soft0)
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[1]];
      if (soft1 && !soft0)
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[2]];
      if (soft1 && soft0)
        vn.c[j].b[k].likelihood =
            likelihood_levels.level[ldpc_decoder_parameters.likelihood_map[3]];
    }
  }

  while ((iteration < ldpc_decoder_input.iteration_limit) && !finished) {
    // Initial syndrome calculation is iteration 0.  At each real iteration,
    // col0/col1 read this settled snapshot; col2 first sees col0's update.
    // Temporarily disable cn_delay for weight calculation debug.
    // cn_delay1 = cn;
    // cn_delay2 = cn;
    log_bf_ibex_sw_rtl(sw_delta_fp, iteration, "iter_pre", -1,
                       syndrome_weight);
    log_bf_ibex_row_sw_rtl(row_sw_fp, iteration, "iter_pre", -1, h_matrix, cn);
    for (int j = 0; j < h_matrix.cols; j++) {
      // Advance the syndrome-weight FIFO by RTL clock slots, not logical cols.
      const int last_col = h_matrix.cols - 1;
      const bool odd_col_count = ((h_matrix.cols & 1) != 0);
      const bool init_col0_hold = (iteration == 1) && (j == 0);
      const bool last_col_hold = (iteration >= 1) && (j == last_col);
      const bool steady_col0_hold =
          odd_col_count && (iteration >= 3) && (j == 0);
      const bool sample_before_shift = false;

      int sw_shift_cycles = 1;
      if (init_col0_hold || last_col_hold || steady_col0_hold)
        sw_shift_cycles = 2;

      int sw_pre_shift_cycles = sw_shift_cycles;
      int sw_post_shift_cycles = 0;
      if (last_col_hold) {
        sw_pre_shift_cycles = 1;
        sw_post_shift_cycles = sw_shift_cycles - 1;
      }

      clock_cycles += sw_shift_cycles;

      if ((iteration == post_start_iteration) && (j == 0)) {
        for (int i = 0; i < 256; i++)
          prng_256.b[i] = (prng_init[int(i / 16)] >> (i % 16)) & 1;
        for (int i = 0; i < 512; i++)
          prng_512.b[i] = (0x1fe0 >> (i % 16)) & 1;
      } else if (iteration >= post_start_iteration) {
        prng_256 = f_256_bit_lfsr(prng_256);
        prng_512 = f_512_bit_lfsr(prng_512);
      }
#ifdef _LDPC_DBG_DUMP
      if ((h_matrix.bits == 512) &&
          (iteration >= ldpc_decoder_input.post_iteration))
        dump_bf_ibex_rtl_cn_prng512(dbg_prng512_fp, prng_512, iteration, j);
#endif

      int sw_dly_sample = syndrome_weight_delayed;
      if (sample_before_shift)
        sw_dly_sample = syndrome_weight_r[4];

      const int sw_sample_clk =
          steady_col0_hold ? (sw_pre_shift_cycles - 1) : 0;
      for (int sw_clk = 0; sw_clk < sw_pre_shift_cycles; sw_clk++) {
        syndrome_weight_r[4] = syndrome_weight_r[3];
        syndrome_weight_r[3] = syndrome_weight_r[2];
        syndrome_weight_r[2] = syndrome_weight_r[1];
        syndrome_weight_r[1] = syndrome_weight_r[0];
        syndrome_weight_r[0] = syndrome_weight;

        // A held RTL column can occupy multiple clocks.  Use the first clock's
        // syndw for the column behavior while still advancing the FIFO through
        // all consumed clocks for subsequent columns.
        if (!sample_before_shift && (sw_clk == sw_sample_clk))
          sw_dly_sample = syndrome_weight_r[4];
      }
      syndrome_weight_delayed = sw_dly_sample;

      bool post_trigger = false;
      bool post_trigger2 = false;
      if (iteration >= ldpc_decoder_input.post_iteration) {
        const bool hamming_weight_lt_post_thr =
            syndrome_weight_delayed <
            ldpc_decoder_parameters.syndrome_weight_thr_post;
        const bool hamming_weight_lt_circ_thr =
            syndrome_weight_delayed <
            ldpc_decoder_parameters.syndrome_weight_thr_qc;
        post_trigger = hamming_weight_lt_circ_thr &&
                       ((iteration % 16) <
                        ldpc_decoder_parameters.post_ratio) &&
                       ldpc_decoder_parameters.post_process_en;
        post_trigger2 = hamming_weight_lt_post_thr &&
                        ((iteration % 16) >=
                         ldpc_decoder_parameters.post_ratio) &&
                        ldpc_decoder_parameters.post_process_en;
      }

      const int prev_sw =
          (iteration == 0) ? syndrome_weight_delayed : syndrome_weight_r[3];
      const bool pushing = (syndrome_weight_delayed >= prev_sw);

#ifdef _LDPC_DBG_DUMP
      if (iteration >= ldpc_decoder_input.post_iteration) {
        log_bf_ibex_rtl_cn_post_ctrl(
            dbg_post_ctrl_fp, iteration, j, ldpc_decoder_input.post_iteration,
            syndrome_weight, syndrome_weight_delayed, prev_sw, pushing,
            post_trigger, post_trigger2, likelihood_levels.min,
            likelihood_levels.flip_thr,
            ldpc_decoder_parameters.syndrome_weight_thr_post,
            ldpc_decoder_parameters.syndrome_weight_thr_qc,
            ldpc_decoder_parameters.post_ratio);
      }
#endif

      if (diag_ab && !diag_mismatch_found && (iteration == 0)) {
        for (int i = 0; i < h_matrix.rows; i++) {
          if (rtl_view.shift[i][j] < 0)
            continue;
          if (accum_rot[i] != rtl_view.shift[i][j]) {
            printf("[DIAG] ROT TRACK MISMATCH iter=%d col=%d row=%d accum=%d expected=%d (element=%d)\n",
                   iteration, j, i, accum_rot[i], rtl_view.shift[i][j],
                   h_matrix.element[i][j]);
            diag_mismatch_found = true;
          }
        }
      }

      int col_post1_cnt = 0;
      int col_post2_cnt = 0;
      int col_aggressive_cnt = 0;
      int col_flipchg_cnt = 0;
      unsigned char col_unsat_cnt[LDPC_MAX_CIRC_BITS] = {0};
      // const s_check_nodes &cn_for_weight = cn_delay2;
      const s_check_nodes &cn_for_weight = cn;
      for (int k = 0; k < h_matrix.bits; k++) {
        bool do_not_use_this_bit = false;
        do_not_use_this_bit |= ((h_matrix.extra_bits_of_parity > 0) &&
                                (j == (h_matrix.cols - h_matrix.rows)) &&
                                (k >= h_matrix.extra_bits_of_parity));
        do_not_use_this_bit |= ((h_matrix.extra_bits_of_userdata > 0) &&
                                (j == (h_matrix.cols - h_matrix.rows - 1)) &&
                                (k >= h_matrix.extra_bits_of_userdata));
        if (do_not_use_this_bit)
          continue;

        int weight = 0;
        for (int i = 0; i < h_matrix.rows; i++) {
          const int shift = rtl_view.shift[i][j];
          if (shift < 0)
            continue;
          if (h_matrix.extra_bytes_of_parity == 0) {
            if (h_matrix.occupied[i][j] && cn_for_weight.r[i].b[k])
              weight++;
          } else {
            if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)) &&
                cn_for_weight.r[i].b[k])
              weight++;
            if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
                cn_for_weight.r[i].b[k] && rtl_view.mask[j][k])
              weight++;
            if (h_matrix.fade[i][j] && cn_for_weight.r[i].b[k] &&
                !rtl_view.mask[j][k])
              weight++;
          }
        }
        col_unsat_cnt[k] = (unsigned char)weight;

        if (diag_ab && !diag_mismatch_found && (iteration == 0) && (k < 8)) {
          int ref_weight = 0;
          for (int i = 0; i < h_matrix.rows; i++) {
            if (rtl_view.shift[i][j] < 0)
              continue;
            const int m = (k + h_matrix.bits - h_matrix.element[i][j]) %
                          h_matrix.bits;
            if (h_matrix.extra_bytes_of_parity == 0) {
              if (h_matrix.occupied[i][j] && cn_ref.r[i].b[m])
                ref_weight++;
            } else {
              if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)) &&
                  cn_ref.r[i].b[m])
                ref_weight++;
              if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
                  cn_ref.r[i].b[m] && h_matrix.mask[j][k])
                ref_weight++;
              if (h_matrix.fade[i][j] && cn_ref.r[i].b[m] &&
                  !h_matrix.mask[j][k])
                ref_weight++;
            }
          }
          if (weight != ref_weight) {
            printf("[DIAG] WEIGHT MISMATCH iter=%d col=%d bit=%d: rtl=%d ref=%d\n",
                   iteration, j, k, weight, ref_weight);
            for (int i = 0; i < h_matrix.rows; i++) {
              if (rtl_view.shift[i][j] < 0) continue;
              const int m = (k + h_matrix.bits - h_matrix.element[i][j]) %
                            h_matrix.bits;
              printf("  row=%d occ=%d fade=%d shift=%d element=%d accum=%d "
                     "cn_rtl[%d]=%d cn_ref[%d]=%d\n",
                     i, (int)h_matrix.occupied[i][j], (int)h_matrix.fade[i][j],
                     rtl_view.shift[i][j], h_matrix.element[i][j], accum_rot[i],
                     k, (int)cn.r[i].b[k], m, (int)cn_ref.r[i].b[m]);
            }
            diag_mismatch_found = true;
          }
        }

        const bool flipped_prev = vn.c[j].b[k].flipped;
        bool prng_post_process = false;
        bool prng_post_process2 = false;
        if (iteration >= ldpc_decoder_input.post_iteration) {
          prng_post_process =
              post_trigger &&
              (syndrome_weight_delayed <
               ldpc_decoder_parameters.syndrome_weight_thr_qc) &&
              ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
          prng_post_process2 =
              post_trigger2 &&
              (syndrome_weight_delayed <
               ldpc_decoder_parameters.syndrome_weight_thr_post) &&
              ((h_matrix.bits == 512) ? prng_512.b[k] : prng_256.b[k]);
        }
        if (prng_post_process)
          col_post1_cnt++;
        if (prng_post_process2)
          col_post2_cnt++;

        const bool aggressive =
            (ldpc_decoder_input.soft_bits > 0) &&
            (likelihood_levels.min < ldpc_decoder_parameters.likelihood_thr &&
             !flipped_prev);
        if (aggressive)
          col_aggressive_cnt++;
        int adjusted_weight = weight;
        if (aggressive && (weight == 2))
          adjusted_weight = 3;
        if (aggressive && (weight == 3))
          adjusted_weight = 5;
        if (aggressive && (weight == 4))
          adjusted_weight = 7;

        vn.c[j].b[k].likelihood = f_update_vn_post(
            vn.c[j].b[k].likelihood, adjusted_weight, likelihood_levels.min,
            likelihood_levels.max, prng_post_process, prng_post_process2,
            aggressive, likelihood_levels.flip_thr, pushing);
        vn.c[j].b[k].flipped =
            (vn.c[j].b[k].likelihood >= likelihood_levels.flip_thr);

        if (flipped_prev != vn.c[j].b[k].flipped) {
          col_flipchg_cnt++;
          for (int i = 0; i < h_matrix.rows; i++) {
            const int shift = rtl_view.shift[i][j];
            if (shift < 0)
              continue;
            if (h_matrix.extra_bytes_of_parity == 0) {
              if (h_matrix.occupied[i][j])
                cn.r[i].b[k] = !cn.r[i].b[k];
            } else {
              if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)))
                cn.r[i].b[k] = !cn.r[i].b[k];
              if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
                  rtl_view.mask[j][k])
                cn.r[i].b[k] = !cn.r[i].b[k];
              if (h_matrix.fade[i][j] && !rtl_view.mask[j][k])
                cn.r[i].b[k] = !cn.r[i].b[k];
            }
          }
          if (diag_ab) {
            for (int i = 0; i < h_matrix.rows; i++) {
              if (rtl_view.shift[i][j] < 0)
                continue;
              const int m = (k + h_matrix.bits - h_matrix.element[i][j]) %
                            h_matrix.bits;
              if (h_matrix.extra_bytes_of_parity == 0) {
                if (h_matrix.occupied[i][j])
                  cn_ref.r[i].b[m] = !cn_ref.r[i].b[m];
              } else {
                if (h_matrix.occupied[i][j] && (i < (h_matrix.rows - 1)))
                  cn_ref.r[i].b[m] = !cn_ref.r[i].b[m];
                if (h_matrix.occupied[i][j] && (i == (h_matrix.rows - 1)) &&
                    h_matrix.mask[j][k])
                  cn_ref.r[i].b[m] = !cn_ref.r[i].b[m];
                if (h_matrix.fade[i][j] && !h_matrix.mask[j][k])
                  cn_ref.r[i].b[m] = !cn_ref.r[i].b[m];
              }
            }
          }
        }
      }

      const int syndrome_weight_after_update =
          f_check_node_weight(h_matrix, cn);

#ifdef _LDPC_DBG_DUMP
      if (iteration >= ldpc_decoder_input.post_iteration) {
        dump_bf_ibex_rtl_cn_col_unsat(dbg_unsat_fp, col_unsat_cnt, iteration,
                                      j, h_matrix.bits);
        dump_bf_ibex_rtl_cn_col_likelihood(dbg_like_fp, vn, iteration, j,
                                           h_matrix.bits);
      }
      if (iteration >= ldpc_decoder_input.post_iteration) {
        log_bf_ibex_rtl_cn_post_eff(dbg_post_eff_fp, iteration, j,
                                    col_post1_cnt, col_post2_cnt,
                                    col_aggressive_cnt, col_flipchg_cnt,
                                    syndrome_weight_after_update);
        dump_bf_ibex_rtl_cn_post_trace(
            dbg_post_trace_fp, vn, iteration, j,
            ldpc_decoder_input.post_iteration, syndrome_weight,
            syndrome_weight_delayed, post_trigger, post_trigger2,
            likelihood_levels.min, likelihood_levels.flip_thr,
            ldpc_decoder_parameters.syndrome_weight_thr_post,
            ldpc_decoder_parameters.syndrome_weight_thr_qc,
            ldpc_decoder_parameters.post_ratio, col_post1_cnt, col_post2_cnt,
            col_aggressive_cnt, col_flipchg_cnt, syndrome_weight_after_update,
            h_matrix.bits);
      }
      dump_cn_synd_snapshot_rtl(dbg_trace_fp, h_matrix, rtl_view, cn,
                                accum_rot, iteration, j, "post_update",
                                syndrome_weight_after_update, true);
#endif

      for (int i = 0; i < h_matrix.rows; i++) {
        if (rtl_view.shift[i][j] >= 0) {
          rotate_cn_row(&cn.r[i], h_matrix.bits, h_matrix.delta[i]);
          // rotate_cn_row(&cn_delay1.r[i], h_matrix.bits, h_matrix.delta[i]);
          // rotate_cn_row(&cn_delay2.r[i], h_matrix.bits, h_matrix.delta[i]);
          accum_rot[i] = (accum_rot[i] + h_matrix.delta[i]) % h_matrix.bits;
        }
      }

      syndrome_weight = f_check_node_weight(h_matrix, cn);
      for (int sw_clk = 0; sw_clk < sw_post_shift_cycles; sw_clk++) {
        syndrome_weight_r[4] = syndrome_weight_r[3];
        syndrome_weight_r[3] = syndrome_weight_r[2];
        syndrome_weight_r[2] = syndrome_weight_r[1];
        syndrome_weight_r[1] = syndrome_weight_r[0];
        syndrome_weight_r[0] = syndrome_weight;
      }
#ifdef _LDPC_DBG_DUMP
      dump_cn_synd_snapshot_rtl(dbg_trace_fp, h_matrix, rtl_view, cn, accum_rot,
                                iteration, j, "post_delta", syndrome_weight,
                                false);
#endif
      if (sw_delta_fp)
        log_bf_ibex_sw_rtl(sw_delta_fp, iteration, "col_post", j,
                           syndrome_weight);
      log_bf_ibex_row_sw_rtl(row_sw_fp, iteration, "col_post", j, h_matrix, cn);
      ldpc_decoder_output.col_cnt = j;
      // cn_delay2 = cn_delay1;
      // cn_delay1 = cn;
      if (syndrome_weight == 0) {
        finished = true;
        break;
      }
    }

#ifdef _LDPC_DBG_DUMP
    if (!finished) {
      for (int i = 0; i < h_matrix.rows; i++) {
        const int wrap_shift =
            (rtl_view.wrap_base[i] +
             (rtl_view.wrap_base_delta[i] * h_matrix.delta[i])) %
            h_matrix.bits;
        rotate_cn_row(&cn.r[i], h_matrix.bits, wrap_shift);
        accum_rot[i] = (accum_rot[i] + wrap_shift) % h_matrix.bits;
      }

      dump_cn_synd_snapshot_rtl(dbg_trace_fp, h_matrix, rtl_view, cn, accum_rot,
                                iteration, -1, "post_wrap", syndrome_weight,
                                false);

      if (diag_ab && !diag_mismatch_found && (iteration == 0)) {
        printf("[DIAG] Wrap verification after iter=%d:\n", iteration);
        for (int i = 0; i < h_matrix.rows; i++) {
          const int expected_fe =
              first_active_shift_for_row(h_matrix, rtl_view, i);
          const int wrap_shift =
              (rtl_view.wrap_base[i] +
               (rtl_view.wrap_base_delta[i] * h_matrix.delta[i])) %
              h_matrix.bits;
          if (accum_rot[i] != expected_fe) {
            printf("[DIAG] WRAP MISMATCH row=%d accum=%d expected_fe=%d "
                   "wb=%d wbd=%d delta=%d wrap_shift=%d\n",
                   i, accum_rot[i], expected_fe,
                   rtl_view.wrap_base[i], rtl_view.wrap_base_delta[i],
                   h_matrix.delta[i], wrap_shift);
            diag_mismatch_found = true;
          } else {
            printf("[DIAG]  row=%d OK (accum=%d == FE=%d) wb=%d wbd=%d delta=%d\n",
                   i, accum_rot[i], expected_fe,
                   rtl_view.wrap_base[i], rtl_view.wrap_base_delta[i],
                   h_matrix.delta[i]);
          }
        }
      }
    } else {
      dump_cn_synd_snapshot_rtl(dbg_trace_fp, h_matrix, rtl_view, cn, accum_rot,
                                iteration, -1, "iter_end", syndrome_weight,
                                false);
    }
#endif

    iteration++;
  }

#ifdef _LDPC_DBG_DUMP
  if (dbg_trace_fp)
    fclose(dbg_trace_fp);
  if (dbg_like_fp)
    fclose(dbg_like_fp);
  if (dbg_post_ctrl_fp)
    fclose(dbg_post_ctrl_fp);
  if (dbg_post_eff_fp)
    fclose(dbg_post_eff_fp);
  if (dbg_post_trace_fp)
    fclose(dbg_post_trace_fp);
  if (dbg_unsat_fp)
    fclose(dbg_unsat_fp);
  if (dbg_prng512_fp)
    fclose(dbg_prng512_fp);
#else
  if (sw_delta_fp)
    fclose(sw_delta_fp);
  if (row_sw_fp)
    fclose(row_sw_fp);
#endif

  ldpc_decoder_output.iterations = iteration;
  ldpc_decoder_output.clock_cycles = clock_cycles;
  ldpc_decoder_output.syndrome_weight_after = f_check_node_weight(h_matrix, cn);
  ldpc_decoder_output.failure =
      (ldpc_decoder_output.syndrome_weight_after != 0);
  ldpc_decoder_output.errors_in_userdata = 0;
  ldpc_decoder_output.errors_in_codeword = 0;
  for (int j = 0; j < h_matrix.cols; j++) {
    for (int k = 0; k < h_matrix.bits; k++) {
      ldpc_decoder_output.corrected_codeword.c[j].b[k] =
          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard ^
          vn.c[j].b[k].flipped;
      dec_do_blk[j * h_matrix.bits + k] =
          ldpc_decoder_output.corrected_codeword.c[j].b[k];
      if (vn.c[j].b[k].flipped)
        ldpc_decoder_output.errors_in_codeword++;
      if (vn.c[j].b[k].flipped && (j < (h_matrix.cols - h_matrix.rows)))
        ldpc_decoder_output.errors_in_userdata++;
    }
  }

  if (ldpc_decoder_output.early_termination) {
    ldpc_decoder_output.iterations = 0;
    ldpc_decoder_output.syndrome_weight_after = syndrome_weight;
    ldpc_decoder_output.failure = 1;
    ldpc_decoder_output.errors_in_userdata = 0;
    ldpc_decoder_output.errors_in_codeword = 0;
  } else if (ldpc_decoder_input.syndrome_cal_only) {
    ldpc_decoder_output.iterations = 0;
    ldpc_decoder_output.syndrome_weight_after = syndrome_weight;
    ldpc_decoder_output.failure = 0;
    ldpc_decoder_output.errors_in_userdata = 0;
    ldpc_decoder_output.errors_in_codeword = 0;
  } else if (ldpc_decoder_output.syndrome_weight_before == 0) {
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
    printf("### C++ RTL CN OUTPUT ###\n");
    f_print_hard_codeword(ldpc_decoder_output.corrected_codeword, h_matrix.cols,
                          h_matrix.bits);
  }

  init_synd_wt = ldpc_decoder_output.syndrome_weight_before;
  fina_synd_wt = ldpc_decoder_output.syndrome_weight_after;
  if (ldpc_decoder_output.failure)
    cnvg_itr = ldpc_decoder_input.iteration_limit - 1;
  else
    cnvg_itr = ldpc_decoder_output.iterations - 1;

  if (ldpc_decoder_output.failure || (ldpc_decoder_output.col_cnt < 0))
    cnvg_lyr = 0;
  else
    cnvg_lyr = ldpc_decoder_output.col_cnt;

  cw_fail = ldpc_decoder_output.failure;
}
