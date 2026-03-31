#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "svdpi.h"

extern "C" {
void ldpc_config(int h_m, int h_n, int h_sc, int h_st, int h_wt, int info_num,
                 int pad_bit, int ch_mode, int ch_para, int fdec_max_itr,
                 int ldec_max_itr, int llr0, int llr1, int alpha,
                 int finite_q_num, int finite_r_num, int sd_num,
                 svOpenArrayHandle v_ref_sv, int debug,
                 int sdlite_llr_config, int sdlite_llr0, int sdlite_llr1,
                 int sdlite_llr2, int sdlite_llr3, int dv_user_data_bytes,
                 int dv_parity_bytes, int post_iter, int max_iter,
                 int nand_strobes, int codeword_4k_8k,
                 int syndrome_weight_thr_qc, int syndrome_weight_thr_post,
                 int early_termination_dis, int post_process_en,
                 int ldpc_decoder_control_likelihood_0,
                 int ldpc_decoder_control_likelihood_1,
                 int ldpc_decoder_control_likelihood_2,
                 int ldpc_decoder_control_likelihood_3,
                 int ldpc_decoder_control_post, unsigned int ldpc_early_term_0,
                 unsigned int ldpc_early_term_1, unsigned int ldpc_early_term_2,
                 unsigned int ldpc_early_term_3, unsigned int ldpc_early_term_4,
                 unsigned int ldpc_early_term_5, unsigned int ldpc_early_term_6);
void ldpc_enc(svOpenArrayHandle usr_data_sv, svOpenArrayHandle enc_data_sv,
              int debug);
void ldpc_dec(svOpenArrayHandle det_data_sv, svOpenArrayHandle dec_data_sv,
              int *dec_unc_sv, int *init_synd_wt_sv, int *dec_cnvg_itr_sv,
              int *dec_cnvg_col_sv, int *fina_synd_wt_sv, int dec_mode,
              int sd_num, int debug, int h_n, int h_sc);
void ldpc_cleanup();
}

namespace {

struct Scenario {
  const char *name;
  int h_m_cfg;
  int h_n_cfg;
  int bm_n_out;
  int pad_bit;
  int user_bytes;
  int parity_bytes;
};

test_sv_array_t make_handle(void *data, int len, size_t elem_size) {
  test_sv_array_t handle;
  handle.data = data;
  handle.len = len;
  handle.elem_size = elem_size;
  return handle;
}

void fill_user_words(std::vector<int> &words, int info_bits) {
  for (int bit = 0; bit < info_bits; bit++) {
    const int value = ((bit / 11) ^ bit) & 1;
    const int word = bit / 32;
    const int bit_pos = 31 - (bit % 32);
    words[static_cast<size_t>(word)] |= (value << bit_pos);
  }
}

void run_scenario(const Scenario &scenario) {
  const int info_bits = scenario.user_bytes * 8;
  const int blk_bits = (scenario.user_bytes + scenario.parity_bytes) * 8;
  const int usr_words = (info_bits + 31) / 32;
  const int blk_words = (blk_bits + 31) / 32;
  const int out_words = (scenario.bm_n_out * 512) / 32;

  std::vector<int> vref(1, 0);
  std::vector<int> usr_data(static_cast<size_t>(usr_words), 0);
  std::vector<int> enc_data(static_cast<size_t>(blk_words), 0);
  std::vector<int> det_data(static_cast<size_t>(blk_words), 0);
  std::vector<int> dec_data(static_cast<size_t>(out_words), 0);

  fill_user_words(usr_data, info_bits);

  test_sv_array_t vref_h = make_handle(vref.data(), static_cast<int>(vref.size()),
                                       sizeof(vref[0]));
  test_sv_array_t usr_h = make_handle(usr_data.data(),
                                      static_cast<int>(usr_data.size()),
                                      sizeof(usr_data[0]));
  test_sv_array_t enc_h = make_handle(enc_data.data(),
                                      static_cast<int>(enc_data.size()),
                                      sizeof(enc_data[0]));
  test_sv_array_t det_h = make_handle(det_data.data(),
                                      static_cast<int>(det_data.size()),
                                      sizeof(det_data[0]));
  test_sv_array_t dec_h = make_handle(dec_data.data(),
                                      static_cast<int>(dec_data.size()),
                                      sizeof(dec_data[0]));

  int dec_unc = -1;
  int init_synd = -1;
  int cnvg_itr = -1;
  int cnvg_col = -1;
  int fina_synd = -1;

  std::cout << "[CMODEL TEST] " << scenario.name
            << " user_bytes=" << scenario.user_bytes
            << " parity_bytes=" << scenario.parity_bytes
            << " h_m=" << scenario.h_m_cfg
            << " h_n=" << scenario.h_n_cfg
            << " pad_bit=" << scenario.pad_bit << std::endl;

  ldpc_config(
      scenario.h_m_cfg, scenario.h_n_cfg, 512, 5, 4, info_bits,
      scenario.pad_bit, 0, 0, 24, 32, 30, -30, 625, 8, 6, 1, &vref_h,
      0,  // debug
      0, 0, 0, 0, 0,  // sdlite
      scenario.user_bytes, scenario.parity_bytes,  // dv bytes
      0, 48, 1, 0,    // post_iter, max_iter, nand_strobes, codeword_4k_8k
      0, 0, 0, 0,     // syndrome thr / early term / post process
      0, 0, 0, 0, 0,  // likelihood ctrl 0..3, post
      0u, 0u, 0u, 0u, 0u, 0u, 0u);  // early_term 0..6
  ldpc_enc(&usr_h, &enc_h, 0);

  std::memcpy(det_data.data(), enc_data.data(), blk_words * sizeof(int));
  ldpc_dec(&det_h, &dec_h, &dec_unc, &init_synd, &cnvg_itr, &cnvg_col,
           &fina_synd, 0, 1, 0, scenario.bm_n_out, 512);
  ldpc_cleanup();

  std::cout << "[CMODEL TEST] status"
            << " dec_unc=" << dec_unc
            << " init_synd=" << init_synd
            << " fina_synd=" << fina_synd
            << " cnvg_itr=" << cnvg_itr
            << " cnvg_col=" << cnvg_col << std::endl;
}

} // namespace

int main() {
  const Scenario scenarios[] = {
      {"parity320B", 5, 69, 70, 0, 4096, 320},
      {"parity321B", 5, 69, 70, 8, 4096, 321},
  };

  try {
    for (const Scenario &scenario : scenarios)
      run_scenario(scenario);
  } catch (const std::exception &ex) {
    std::cerr << "[CMODEL TEST ERROR] " << ex.what() << std::endl;
    ldpc_cleanup();
    return 1;
  }

  std::cout << "[CMODEL TEST] completed successfully" << std::endl;
  return 0;
}
