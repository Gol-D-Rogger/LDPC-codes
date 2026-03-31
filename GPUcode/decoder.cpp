#include <iomanip>
#include <iostream>

#include "decoder.h"

uint64_t logger::ITER_LIMIT = logger::MAX_ITER;

uint16_t decoder::prng_init[32] = {
    0x9365, 0x49f9, 0xda8f, 0xf7b2, 0x30ee, 0xef08, 0x1b73, 0x8c9a,
    0xc646, 0xb550, 0x2edb, 0x71cc, 0x5d27, 0xa8a1, 0x6214, 0x043d,
    0x9375, 0x89f9, 0xd98f, 0xf2b2, 0x31ee, 0xef18, 0x0b73, 0x4c9a,
    0xc946, 0xba50, 0x3edb, 0x71bc, 0x5327, 0xa8a1, 0x3214, 0x083d,
};

ldpc_decoder_input::ldpc_decoder_input() = default;

ldpc_decoder_input::ldpc_decoder_input(uint64_t nand_strobes,
                                       uint64_t soft_bits,
                                       uint64_t iteration_limit,
                                       uint64_t post_iteration) {
  this->nand_strobes = nand_strobes;
  this->soft_bits = soft_bits;
  this->iteration_limit = iteration_limit;
  this->post_iteration = post_iteration;
}

void ldpc_decoder_input::clear_cw(void) { corrupted_codeword.clear(); }

void ldpc_decoder_output::print_stats(void) {
  std::cout << "===== Codeword Decoding Stats =====\n";
  std::cout << "failure: " << failure << "\n";
  std::cout << "iterations: " << iterations << "\n";
  std::cout << "clock_cycles: " << clock_cycles << "\n";
  std::cout << "syndrome_weight_before: " << syndrome_weight_before << "\n";
  std::cout << "syndrome_weight_after: " << syndrome_weight_after << "\n";
  std::cout << "total_errors: " << total_errors << "\n";
  std::cout << "net_acc: " << net_acc << "\n";
}

void decoder_output_acc::print_stats(void) {
  std::cout << "===== Accumulated Decoding Stats =====\n";
  std::cout << "failure: " << failure << "\n";
  std::cout << "iterations: " << iterations << "\n";
  std::cout << "clock_cycles: " << clock_cycles << "\n";
  std::cout << "syndrome_weight_before: " << syndrome_weight_before << "\n";
  std::cout << "syndrome_weight_after: " << syndrome_weight_after << "\n";
  std::cout << "total_errors: " << total_errors << "\n";
  std::cout << "net_acc: " << net_acc << "\n";
}

decoder::decoder(h_matrix &h_matrix_ref, int nand_strobes, int post_ratio)
    : h_matrix_ref(h_matrix_ref) {
  ldpc_decoder_parameters.post_ratio = static_cast<uint16_t>(post_ratio);
  ldpc_decoder_parameters.likelihood_init_fraction[0] = 0;
  if (nand_strobes == 7) {
    ldpc_decoder_parameters.likelihood_init_fraction[1] = 8;
    ldpc_decoder_parameters.likelihood_init_fraction[2] = 4;
  } else {
    ldpc_decoder_parameters.likelihood_init_fraction[1] = 6;
    ldpc_decoder_parameters.likelihood_init_fraction[2] = 6;
  }
}

ldpc_decoder_output decoder::decode_planar(const ldpc_decoder_input &input) {
  ldpc_decoder_output out;
  out.net_acc = 1;

  check_nodes cn =
      config_check_nodes(h_matrix_ref, input.corrupted_codeword.hard);
  const uint64_t syndrome_weight = check_node_weight(cn);
  out.syndrome_weight_before = static_cast<uint32_t>(syndrome_weight);
  out.syndrome_weight_after = static_cast<uint32_t>(syndrome_weight);
  out.clock_cycles = static_cast<uint32_t>(h_matrix_ref.cols + 1);
  out.iterations = 0;
  out.total_errors = 0;
  out.failure = (syndrome_weight != 0) ? 1U : 0U;

  return out;
}

void decoder::update_vn_planar(uint64_t *vn_flipped, const uint8_t *weight_word,
                               unsigned char *vn_likelihood,
                               const s_likelihood_levels &likelihood_levels,
                               const bool *prng_512, bool post_trigger1,
                               bool post_trigger2, int bit_offset) {
  for (int i = 0; i < 8; ++i)
    vn_flipped[i] = 0ULL;

  for (int bit = 0; bit < 512; ++bit) {
    const int idx = bit_offset + bit;
    if (idx < 0)
      continue;

    const uint8_t w = weight_word[bit];
    unsigned char lk = vn_likelihood[idx];

    if (post_trigger1 || post_trigger2) {
      if (prng_512[bit]) {
        if (lk < static_cast<unsigned char>(likelihood_levels.level[3]))
          lk++;
      } else {
        if (lk > static_cast<unsigned char>(likelihood_levels.level[0]))
          lk--;
      }
    }

    if (w >= static_cast<uint8_t>(likelihood_levels.flip_thr)) {
      vn_flipped[bit >> 6] |= (UNIT << (bit & 63));
    }

    vn_likelihood[idx] = lk;
  }
}

s_likelihood_levels decoder::compute_likelihood_levels(
    int strobes, decoder::ldpc_decoder_params ldpc_decoder_parameters,
    int syndrome_weight, int rows) {
  s_likelihood_levels levels{};
  levels.min = 0;
  levels.max = (1 << ldpc_decoder_parameters.VN_BITS) - 1;
  levels.flip_thr = ldpc_decoder_parameters.likelihood_thr;
  levels.weak = levels.flip_thr / 2;
  levels.strong = levels.flip_thr;

  const int frac_idx = (strobes >= 7) ? 2 : 1;
  const int base = ldpc_decoder_parameters.likelihood_init_fraction[frac_idx];
  const int pressure = (rows > 0) ? (syndrome_weight * 16 / rows) : 0;

  levels.level[0] = 0;
  levels.level[1] = std::min(levels.max, base + pressure / 4);
  levels.level[2] = std::min(levels.max, base * 2 + pressure / 2);
  levels.level[3] = std::min(levels.max, base * 3 + pressure);
  return levels;
}

s_512_bits decoder::lfsr_512_bit(s_512_bits data_in) {
  s_512_bits data_out{};
  for (int i = 511; i > 0; --i)
    data_out.b[i] = data_in.b[i - 1];
  data_out.b[0] =
      data_in.b[511] ^ data_in.b[509] ^ data_in.b[502] ^ data_in.b[500];
  return data_out;
}

void logger::log_elapsed_time(
    const std::chrono::steady_clock::time_point &start_time,
    const std::chrono::steady_clock::time_point &end_time) {
  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                      end_time - start_time)
                      .count();
  std::cout << "[DECODER] elapsed: " << ms << " ms\n";
}

void logger::print_accumulated_stats(uint64_t accumulated_cw_count,
                                     ldpc_decoder_output &decode_stats) {
  if (accumulated_cw_count == 0) {
    std::cout << "[DECODER] no accumulated codewords\n";
    return;
  }
  std::cout << "[DECODER] CW: " << accumulated_cw_count
            << " failures: " << decode_stats.failure
            << " failure_rate: " << std::setprecision(8)
            << static_cast<double>(decode_stats.failure) /
                   static_cast<double>(accumulated_cw_count)
            << "\n";
}

void logger::print_accumulated_stats(uint64_t accumulated_cw_count,
                                     decoder_output_acc &decode_stats) {
  if (accumulated_cw_count == 0) {
    std::cout << "[DECODER] no accumulated codewords\n";
    return;
  }
  std::cout << "[DECODER] CW: " << accumulated_cw_count
            << " failures: " << decode_stats.failure
            << " failure_rate: " << std::setprecision(8)
            << static_cast<double>(decode_stats.failure) /
                   static_cast<double>(accumulated_cw_count)
            << "\n";
}

void logger::print_iter_stats(uint32_t *iter_info) {
  if (iter_info == nullptr)
    return;
  std::cout << "[DECODER] iteration histogram (first 16): ";
  for (int i = 0; i < 16; ++i) {
    std::cout << iter_info[i];
    if (i != 15)
      std::cout << ", ";
  }
  std::cout << "\n";
}
