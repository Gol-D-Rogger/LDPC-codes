#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../ldpc_codec.h"

namespace {

enum class RunMode {
  kConfig,
  kEncode,
  kLayer,
  kBfIbex,
  kBfIbexRtlCn,
  kSyndromeAB,
  kAll,
};

struct Scenario {
  int m;
  int k;
  int n;
  int user_bytes;
  int parity_bytes;
};

const char *kDefaultMatrixRoot = "../../../../IBEX/ibex_matrix_flat_13rate";

RunMode parse_mode(const std::string &value) {
  if (value == "config")
    return RunMode::kConfig;
  if (value == "encode")
    return RunMode::kEncode;
  if (value == "layer")
    return RunMode::kLayer;
  if (value == "bf_ibex")
    return RunMode::kBfIbex;
  if (value == "bf_ibex_rtl_cn")
    return RunMode::kBfIbexRtlCn;
  if (value == "syndrome_ab")
    return RunMode::kSyndromeAB;
  if (value == "all")
    return RunMode::kAll;
  throw std::runtime_error("unknown mode: " + value);
}

void fill_payload(ldpc_packet &packet) {
  for (int i = 0; i < packet.info_len; i++)
    packet.usr_blk[i] = static_cast<char>(((i / 7) ^ i) & 1);
}

int count_bit_mismatch(const char *a, const char *b, int n_bits) {
  int mismatches = 0;
  for (int i = 0; i < n_bits; i++) {
    if ((a[i] & 1) != (b[i] & 1))
      mismatches++;
  }
  return mismatches;
}

Scenario make_scenario(int m, int k) {
  Scenario scenario;
  scenario.m = m;
  scenario.k = k;
  scenario.n = m + k;
  scenario.user_bytes = k * 64;
  scenario.parity_bytes = (m == 5) ? 320 : (m * 64 - 1);
  return scenario;
}

Scenario make_scenario_with_bytes(int m, int k, int user_bytes,
                                  int parity_bytes) {
  Scenario scenario = make_scenario(m, k);
  scenario.user_bytes = user_bytes;
  scenario.parity_bytes = parity_bytes;
  return scenario;
}

std::vector<Scenario> default_scenarios() {
  std::vector<Scenario> scenarios;
  for (int m = 5; m <= 17; m++) {
    for (int k = 64; k <= 67; k++)
      scenarios.push_back(make_scenario(m, k));
  }
  return scenarios;
}

void configure_packet(ldpc_packet &packet, const Scenario &scenario) {
  const int info_bits = scenario.user_bytes * 8;
  const int block_bits = (scenario.user_bytes + scenario.parity_bytes) * 8;
  const float vref[1] = {0.0f};

  packet.ch_config(info_bits, block_bits, CLEAN, 0.0f, 8);
  packet.ldpc_config(scenario.m, scenario.n, 512, 5, 4);
  packet.ldpc_pckt_alloc();
  packet.ch_llr_alloc(MANUAL, 1, const_cast<float *>(vref));
  packet.ch_llr_gen(1.875f, -1.875f, 7, 3);
  packet.ldpc_dec_config(24, 0, 32, 0.625f, 1, 8, 6, 4, 0, 0, 0, 0, 0);
  packet.ldpc_ibex_parameters(1, 48, 16, 0,
                              0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0);
  if (!packet.qc_bm || packet.total_cir <= 0)
    throw std::runtime_error("ldpc_config failed");
}

void prepare_clean_channel(ldpc_packet &packet) {
  packet.ch_transmit();
  packet.ch_detector();
  packet.ldpc_ibex_input(0, 50, 0, 1);
}

void inject_manual_error_pattern(ldpc_packet &packet) {
  const int error_indices[] = {7};
  for (int target_idx : error_indices) {
    if ((target_idx < 0) || (target_idx >= packet.blk_len))
      continue;

    int bit_index = 0;
    for (int j = 0; j < packet.h_matrix.cols; j++) {
      for (int k = 0; k < packet.h_matrix.bits; k++) {
        const bool pad_userdata =
            (j == (packet.h_matrix.cols - packet.h_matrix.rows - 1)) &&
            (k >=
             (packet.h_matrix.bits -
              8 * packet.h_matrix.unused_bytes_of_userdata));
        const bool pad_parity =
            (j == (packet.h_matrix.cols - packet.h_matrix.rows)) &&
            (k >=
             (packet.h_matrix.bits -
              8 * packet.h_matrix.unused_bytes_of_parity));
        if (pad_userdata || pad_parity)
          continue;

        if (bit_index == target_idx) {
          packet.ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard =
              packet.tx_blk[target_idx] ? 0 : 1;
          packet.ldpc_decoder_input.corrupted_codeword.c[j]
              .b[k]
              .bit_questionable = 1;
          packet.ldpc_decoder_input.corrupted_codeword.c[j]
              .b[k]
              .bit_questionable2 = 1;
          j = packet.h_matrix.cols;
          break;
        }
        bit_index++;
      }
    }
  }
  packet.raw_err_num = sizeof(error_indices) / sizeof(error_indices[0]);
  packet.ldpc_decoder_input.iteration_limit = 128;
}

void cleanup_packet(ldpc_packet &packet) {
  packet.ch_llr_clean();
  packet.ldpc_pckt_clean();
  packet.ldpc_clean();
}

void run_config_only(const Scenario &scenario) {
  ldpc_packet packet{};
  configure_packet(packet, scenario);
  cleanup_packet(packet);
}

void run_encode_only(const Scenario &scenario) {
  ldpc_packet packet{};
  configure_packet(packet, scenario);
  fill_payload(packet);
  packet.ldpc_ibex_encoder();
  cleanup_packet(packet);
}

void run_layer_decode(const Scenario &scenario) {
  ldpc_packet packet{};
  configure_packet(packet, scenario);
  fill_payload(packet);
  packet.ldpc_ibex_encoder();
  prepare_clean_channel(packet);
  packet.ldpc_decoder(LAYER_G2);

  const int mismatches = count_bit_mismatch(packet.dec_blk, packet.tx_blk, packet.blk_len);
  if (mismatches != 0) {
    std::ostringstream oss;
    oss << "layer decode mismatch: m=" << scenario.m
        << " k=" << scenario.k
        << " mismatches=" << mismatches;
    cleanup_packet(packet);
    throw std::runtime_error(oss.str());
  }

  cleanup_packet(packet);
}

void run_bf_ibex_decode(const Scenario &scenario) {
  ldpc_packet packet{};
  configure_packet(packet, scenario);
  fill_payload(packet);
  packet.ldpc_ibex_encoder();
  prepare_clean_channel(packet);
  packet.ldpc_decoder(BF_IBEX);

  const int mismatches = count_bit_mismatch(packet.dec_blk, packet.tx_blk, packet.blk_len);
  if (mismatches != 0) {
    std::ostringstream oss;
    oss << "bf_ibex decode mismatch: m=" << scenario.m
        << " k=" << scenario.k
        << " mismatches=" << mismatches;
    cleanup_packet(packet);
    throw std::runtime_error(oss.str());
  }

  cleanup_packet(packet);
}

void run_bf_ibex_rtl_cn_decode(const Scenario &scenario) {
  ldpc_packet packet{};
  configure_packet(packet, scenario);
  fill_payload(packet);
  packet.ldpc_ibex_encoder();
  prepare_clean_channel(packet);
  inject_manual_error_pattern(packet);
  packet.ldpc_decoder(BF_IBEX_RTL_CN);

  if (packet.init_synd_wt <= 0) {
    cleanup_packet(packet);
    throw std::runtime_error("bf_ibex_rtl_cn did not exercise non-zero syndrome path");
  }

  const int mismatches = count_bit_mismatch(packet.dec_blk, packet.tx_blk, packet.blk_len);
  if ((mismatches != 0) || (packet.fina_synd_wt != 0)) {
    std::ostringstream oss;
    oss << "bf_ibex_rtl_cn decode mismatch: m=" << scenario.m
        << " k=" << scenario.k
        << " mismatches=" << mismatches
        << " init_sw=" << packet.init_synd_wt
        << " fina_sw=" << packet.fina_synd_wt
        << " cnvg_itr=" << packet.cnvg_itr
        << " cnvg_lyr=" << packet.cnvg_lyr;
    cleanup_packet(packet);
    throw std::runtime_error(oss.str());
  }

  cleanup_packet(packet);
}

void run_syndrome_ab(const Scenario &scenario) {
  ldpc_packet pkt_std{};
  ldpc_packet pkt_rtl;
  configure_packet(pkt_std, scenario);
  configure_packet(pkt_rtl, scenario);
  fill_payload(pkt_std);
  fill_payload(pkt_rtl);
  pkt_std.ldpc_ibex_encoder();
  pkt_rtl.ldpc_ibex_encoder();
  prepare_clean_channel(pkt_std);
  prepare_clean_channel(pkt_rtl);
  inject_manual_error_pattern(pkt_std);
  inject_manual_error_pattern(pkt_rtl);

  pkt_std.ldpc_decoder(BF_IBEX);
  pkt_rtl.ldpc_decoder(BF_IBEX_RTL_CN);

  std::cout << "[AB] m=" << scenario.m << " k=" << scenario.k
            << "  STD: init_sw=" << pkt_std.init_synd_wt
            << " fina_sw=" << pkt_std.fina_synd_wt
            << " cnvg_itr=" << pkt_std.cnvg_itr
            << " cnvg_lyr=" << pkt_std.cnvg_lyr
            << "  RTL: init_sw=" << pkt_rtl.init_synd_wt
            << " fina_sw=" << pkt_rtl.fina_synd_wt
            << " cnvg_itr=" << pkt_rtl.cnvg_itr
            << " cnvg_lyr=" << pkt_rtl.cnvg_lyr
            << std::endl;

  const bool init_match = (pkt_std.init_synd_wt == pkt_rtl.init_synd_wt);
  const bool fina_match = (pkt_std.fina_synd_wt == pkt_rtl.fina_synd_wt);
  const bool itr_match = (pkt_std.cnvg_itr == pkt_rtl.cnvg_itr);
  const bool lyr_match = (pkt_std.cnvg_lyr == pkt_rtl.cnvg_lyr);

  if (!init_match || !fina_match || !itr_match || !lyr_match) {
    std::ostringstream oss;
    oss << "syndrome_ab mismatch: m=" << scenario.m
        << " k=" << scenario.k;
    if (!init_match)
      oss << " init_sw(" << pkt_std.init_synd_wt
          << "!=" << pkt_rtl.init_synd_wt << ")";
    if (!fina_match)
      oss << " fina_sw(" << pkt_std.fina_synd_wt
          << "!=" << pkt_rtl.fina_synd_wt << ")";
    if (!itr_match)
      oss << " cnvg_itr(" << pkt_std.cnvg_itr
          << "!=" << pkt_rtl.cnvg_itr << ")";
    if (!lyr_match)
      oss << " cnvg_lyr(" << pkt_std.cnvg_lyr
          << "!=" << pkt_rtl.cnvg_lyr << ")";
    cleanup_packet(pkt_std);
    cleanup_packet(pkt_rtl);
    throw std::runtime_error(oss.str());
  }

  cleanup_packet(pkt_std);
  cleanup_packet(pkt_rtl);
}

void run_scenario(const Scenario &scenario, RunMode mode) {
  std::cout << "[TEST] m=" << scenario.m
            << " k=" << scenario.k
            << " n=" << scenario.n
            << " user_bytes=" << scenario.user_bytes
            << " parity_bytes=" << scenario.parity_bytes
            << std::endl;

  switch (mode) {
  case RunMode::kConfig:
    run_config_only(scenario);
    break;
  case RunMode::kEncode:
    run_encode_only(scenario);
    break;
  case RunMode::kLayer:
    run_layer_decode(scenario);
    break;
  case RunMode::kBfIbex:
    run_bf_ibex_decode(scenario);
    break;
  case RunMode::kBfIbexRtlCn:
    run_bf_ibex_rtl_cn_decode(scenario);
    break;
  case RunMode::kSyndromeAB:
    run_syndrome_ab(scenario);
    break;
  case RunMode::kAll:
    run_config_only(scenario);
    run_encode_only(scenario);
    run_layer_decode(scenario);
    run_bf_ibex_decode(scenario);
    break;
  }
}

} // namespace

int main(int argc, char **argv) {
  RunMode mode = RunMode::kAll;
  std::vector<Scenario> scenarios = default_scenarios();
  int override_user_bytes = -1;
  int override_parity_bytes = -1;
  int single_m = -1;
  int single_k = -1;

  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if ((arg == "--mode") && (i + 1 < argc)) {
      mode = parse_mode(argv[++i]);
    } else if ((arg == "--matrix-root") && (i + 1 < argc)) {
      setenv("DVC_IBEX_MATRIX_ROOT", argv[++i], 1);
    } else if ((arg == "--single") && (i + 2 < argc)) {
      single_m = std::atoi(argv[++i]);
      single_k = std::atoi(argv[++i]);
    } else if ((arg == "--user-bytes") && (i + 1 < argc)) {
      override_user_bytes = std::atoi(argv[++i]);
    } else if ((arg == "--parity-bytes") && (i + 1 < argc)) {
      override_parity_bytes = std::atoi(argv[++i]);
    } else {
      std::cerr << "Usage: " << argv[0]
                << " [--mode config|encode|layer|bf_ibex|bf_ibex_rtl_cn|syndrome_ab|all]"
                << " [--matrix-root PATH]"
                << " [--single M K]"
                << " [--user-bytes N]"
                << " [--parity-bytes N]"
                << std::endl;
      return 2;
    }
  }

  if ((single_m >= 0) && (single_k >= 0)) {
    scenarios.clear();
    if ((override_user_bytes >= 0) || (override_parity_bytes >= 0)) {
      const Scenario base = make_scenario(single_m, single_k);
      scenarios.push_back(make_scenario_with_bytes(
          single_m, single_k,
          (override_user_bytes >= 0) ? override_user_bytes : base.user_bytes,
          (override_parity_bytes >= 0) ? override_parity_bytes
                                       : base.parity_bytes));
    } else {
      scenarios.push_back(make_scenario(single_m, single_k));
    }
  }

  if (std::getenv("DVC_IBEX_MATRIX_ROOT") == nullptr)
    setenv("DVC_IBEX_MATRIX_ROOT", kDefaultMatrixRoot, 1);

  try {
    for (const Scenario &scenario : scenarios)
      run_scenario(scenario, mode);
  } catch (const std::exception &ex) {
    std::cerr << "[TEST ERROR] " << ex.what() << std::endl;
    return 1;
  }

  std::cout << "[TEST] completed successfully" << std::endl;
  return 0;
}
