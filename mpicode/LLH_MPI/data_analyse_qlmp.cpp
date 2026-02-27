#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <bitset>

using namespace std::chrono;
using namespace std;

#define NCOL 72
#include "ldpc.h"
#include <mpi.h>

struct s_decode_res {
  long int fail_cnt;
  long int total_iter;
  long int syndrome_weight_before;
  long int errors_in_codeword;
  long int actual_codewords;
};

struct s_input_parameter {
  int max_failures;
  int num_errors;
  int bytes_of_parity;
  int max_iter;
  int post_iter;
};

s_input_parameter input_parameter;
int *iter_stat;

void f_print_ldpc_decoder_output(s_ldpc_decoder_output s) {
  printf("        decode fail     =   %3d\n", s.failure);
  printf("        ecnt_userdata   =   %3d\n", s.errors_in_userdata);
  printf("        ecnt_codeword   =   %3d\n", s.errors_in_codeword);
  printf("        iteration       =   %3d\n", s.iterations);
  printf("        syndw_before    =   %3d\n", s.syndrome_weight_before);
  printf("        syndw_after     =   %3d\n", s.syndrome_weight_after);
  printf("        early_term      =   %3d\n", s.early_termination);
};

s_hard_codeword read_data_from_file(const char *filename, int cols,
                                    int bytes_of_userdata,
                                    int bytes_of_parity) {
  s_hard_codeword res;
  std::ifstream infile(filename);
  std::vector<uint8_t> data;

  if (!infile) {
    std::cerr << "Error: Unable to open file " << filename << std::endl;
    return res;
  }

  std::string token;
  while (infile >> token) {
    // convert hex string to uint32_t
    uint32_t value;
    std::stringstream ss(token);
    ss >> std::hex >> value;
    data.push_back(value % 256);
    value /= 256;
    data.push_back(value % 256);
    value /= 256;
    data.push_back(value % 256);
    value /= 256;
    data.push_back(value % 256);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    // printf("data[%d]=%x\n", i, data[i]);
  }

  int col_index = 0;
  int byte_index = 0;

  // user data
  while (byte_index < bytes_of_userdata) {
    uint8_t data_t = data[byte_index];
    int base = (byte_index % 64) * 8;
    for (int ii = 0; ii < 8; ii++) {
      res.c[col_index].b[base + ii] = data_t % 2;
      data_t /= 2;
    }
    byte_index++;
    if (byte_index % 64 == 0 && byte_index > 63)
      col_index++;
  }
  if (byte_index % 64 != 0)
    col_index++;

  // 1st parity
  int parity_byte_index = 0;
  while (byte_index < bytes_of_userdata + (bytes_of_parity % 64)) {
    uint8_t data_t = data[byte_index];
    int base = parity_byte_index * 8;
    for (int ii = 0; ii < 8; ii++) {
      res.c[col_index].b[base + ii] = data_t % 2;
      data_t /= 2;
    }
    byte_index++;
    parity_byte_index++;
  }

  if (bytes_of_parity % 64 != 0)
    col_index++;
  parity_byte_index = 0;
  while (byte_index < bytes_of_userdata + bytes_of_parity) {
    uint8_t data_t = data[byte_index];
    int base = (parity_byte_index % 64) * 8;
    for (int ii = 0; ii < 8; ii++) {
      res.c[col_index].b[base + ii] = data_t % 2;
      data_t /= 2;
    }
    byte_index++;
    parity_byte_index++;
    if (parity_byte_index % 64 == 0 && parity_byte_index > 63)
      col_index++;
  }

  return res;
}

s_decode_res
encode_errinj_decode(s_h_matrix h_matrix,
                     s_ldpc_decoder_parameters ldpc_decoder_parameters,
                     s_ldpc_decoder_input ldpc_decoder_input, int no_err_inj,
                     float hard_rber, long number_of_codewords) {
  s_ldpc_decoder_output ldpc_decoder_output;
  s_decode_res res{0, 0, 0, 0, 0};
  int extra_bits = h_matrix.extra_bits_of_userdata;
  float hard_rber_use;
  s_rbers rbers =
      f_find_error_probabilities(ldpc_decoder_input.nand_strobes, hard_rber);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    printf("rbers=%f    %f    %f    %f    %f    %f    %f\n", rbers.rber[0],
           rbers.rber[1], rbers.rber[2], rbers.rber[3], rbers.rber[4],
           rbers.rber[5], rbers.rber[6]);

  s_hard_codeword random_userdata = f_generate_random_userdata(
      h_matrix.rows, h_matrix.cols, h_matrix.bits, extra_bits);
  s_hard_codeword encoded_codeword = f_ldpc_encode(random_userdata, h_matrix);

  auto start = std::chrono::high_resolution_clock::now();
  for (long int num_cw = 0; num_cw < number_of_codewords; num_cw++) {
    for (int j = 0; j < h_matrix.cols; j++) {
      for (int k = 0; k < h_matrix.bits; k++) {
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_is_error = false;
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = encoded_codeword.c[j].b[k];
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 0;
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 0;
      }
    }
    if (no_err_inj == 0) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dist(0.0f, 1.0f);
      int error_bits = 0;
      int w03 = 0, w13 = 0, w02 = 0, w12 = 0, w01 = 0, w11 = 0, w00 = 0,
          w10 = 0;
      for (int j = 0; j < h_matrix.cols; j++) {
        for (int k = 0; k < h_matrix.bits; k++) {
          double vn;
          bool do_not_use_this_bit =
              ((h_matrix.extra_bits_of_parity > 0) &&
               (j == (h_matrix.cols - h_matrix.rows)) &&
               (k >= h_matrix.extra_bits_of_parity));
          do_not_use_this_bit |=
              ((h_matrix.extra_bits_of_userdata > 0) &&
               (j == (h_matrix.cols - h_matrix.rows - 1)) &&
               (k >= h_matrix.extra_bits_of_userdata));
          if (do_not_use_this_bit)
            vn = 1.00;
          else
            vn = dist(gen);

          if (ldpc_decoder_input.soft_bits == 1) {
            if ((rbers.rber[2] < vn) && (vn < rbers.rber[1])) {
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].level = 3;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 1;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 1;
            } else {
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].level = 0;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 0;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 0;
            }
          } else if (ldpc_decoder_input.soft_bits == 0) {
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 0;
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 0;
          } else if (ldpc_decoder_input.soft_bits == 2) {
            if ((rbers.rber[2] < vn) && (vn < rbers.rber[1])) {
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 1;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 1;

              if (vn < rbers.rber[0])
                w13++;
              else
                w03++;
            } else if ((rbers.rber[4] < vn) && (vn < rbers.rber[3])) {
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 0;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 1;

              if (vn < rbers.rber[0])
                w12++;
              else
                w02++;
            } else if ((rbers.rber[6] < vn) && (vn < rbers.rber[5])) {
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 1;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 0;

              if (vn < rbers.rber[0])
                w11++;
              else
                w01++;
            } else {
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 0;
              ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 0;

              if (vn < rbers.rber[0])
                w10++;
              else
                w00++;
            }
          }

          ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard =
              vn < rbers.rber[0] ? !encoded_codeword.c[j].b[k]
                                 : encoded_codeword.c[j].b[k];
        }
      }
    } else if (no_err_inj == -1) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dist(0.0f, 1.0f);

      static int once = false;
      const unsigned total_bits = h_matrix.cols * h_matrix.bits;
      const int active_bits = total_bits -
                              (512 - h_matrix.extra_bits_of_parity) -
                              (512 - h_matrix.extra_bits_of_userdata);
      int expected_errors = std::round(active_bits * rbers.rber[0]);
      int actual_errors = 0;

      if (!once) {
        std::cout << "active_bits: " << active_bits << "\n"
                  << "expected_errors: " << expected_errors << "\n";
        once = true;
      }

      while (actual_errors < expected_errors) {
        int candidate_pos = gen() % total_bits;
        int j = candidate_pos / h_matrix.bits;
        int k = candidate_pos % h_matrix.bits;

        // check if this bit should not be used for error injection
        bool do_not_use =
            ((h_matrix.extra_bits_of_parity > 0) &&
             (j == (h_matrix.cols - h_matrix.rows)) &&
             (k >= h_matrix.extra_bits_of_parity));
        do_not_use |= ((h_matrix.extra_bits_of_userdata > 0) &&
                       (j == (h_matrix.cols - h_matrix.rows - 1)) &&
                       (k >= h_matrix.extra_bits_of_userdata));

        // Skip if this bit is not allowed for error injection
        if (do_not_use ||
            ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_is_error)
          continue;

        // Mark the bit as error otherwise, and filp the hard value
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_is_error = 1;
        ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = !ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard;
        
        // hard cases don't need any soft info
        actual_errors++;
      }
    }

    ldpc_decoder_output =
        f_ldpc_decode(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);
    if (ldpc_decoder_output.failure)
      res.fail_cnt++;
    res.total_iter += ldpc_decoder_output.iterations;
    res.syndrome_weight_before += ldpc_decoder_output.syndrome_weight_before;
    res.errors_in_codeword += ldpc_decoder_output.errors_in_codeword;
    res.actual_codewords++;
    iter_stat[ldpc_decoder_output.iterations]++;
    long int total_fail_cnt = 0;
    if (num_cw % 1000 == 0) {
      MPI_Allreduce(&res.fail_cnt, &total_fail_cnt, 1, MPI_LONG, MPI_SUM,
                    MPI_COMM_WORLD);
      if (total_fail_cnt >= input_parameter.max_failures)
        return res;
      if (rank == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "Completed " << num_cw
                  << " codewords per process, and total " << total_fail_cnt
                  << " errors, Elapsed time: " << elapsed.count() << " ms\n";

        double time_per_codeword_ms = elapsed.count() / 1000;
        std::cout << "  ->Time per codeword (ms): " << std::fixed
                  << std::setprecision(6) << time_per_codeword_ms << "\n";
      }
    }
  }

  return res;
}

void test_custom_data() {
  s_ldpc_decoder_parameters ldpc_decoder_parameters;
  ldpc_decoder_parameters.post_process_en = 1;
  ldpc_decoder_parameters.syndrome_weight_thr_qc = 16;
  ldpc_decoder_parameters.syndrome_weight_thr_post = 16;
  ldpc_decoder_parameters.likelihood_thr = 128;
  int likelihood_init[8][4] = {{0, 4, 3, 3}, {0, 3, 3, 0}, {0, 3, 3, 2},
                               {0, 3, 2, 2}, {0, 2, 2, 2}, {0, 5, 4, 3},
                               {0, 2, 1, 1}, {0, 1, 1, 1}};
  // std::copy(&likelihood_init[0][0], &likelihood_init[0][0] + 8*4,
  // &ldpc_decoder_parameters.likelihood_init[0][0]);
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 4; j++)
      ldpc_decoder_parameters.likelihood_init_coef_all[i][j] =
          likelihood_init[i][j];

  ldpc_decoder_parameters.post_ratio = 12;
  ldpc_decoder_parameters.likelihood_map[0] = 0;
  ldpc_decoder_parameters.likelihood_map[1] = 1;
  ldpc_decoder_parameters.likelihood_map[2] = 2;
  ldpc_decoder_parameters.likelihood_map[3] = 3;

  ldpc_decoder_parameters.early_terminate_dis = 1;
  ldpc_decoder_parameters.questionable_sense = 1;
  int bytes_of_parity = 672, bytes_of_userdata = 4112;
  int cols = (bytes_of_parity / 64) + (bytes_of_userdata / 64);
  if (bytes_of_parity % 64 != 0)
    cols++;
  if (bytes_of_userdata % 64 != 0)
    cols++;

  s_h_matrix h_matrix = f_h_matrix(bytes_of_userdata, bytes_of_parity);
  s_hard_codeword codeword_in_hard = read_data_from_file(
      "hard_data_0", h_matrix.cols, bytes_of_userdata, bytes_of_parity);
  f_print_hard_codeword(codeword_in_hard, h_matrix.cols, h_matrix.bits);

  // hard decoding
  s_ldpc_decoder_input ldpc_decoder_input;

  ldpc_decoder_input.nand_strobes = 0;
  ldpc_decoder_input.soft_bits = 0;
  ldpc_decoder_input.post_iteration = 1024;
  ldpc_decoder_input.iteration_limit = 5048;
  ldpc_decoder_input.syndrome_cal_only = 0;
  s_codeword corrupted_codeword;
  for (int j = 0; j < h_matrix.cols; j++) {
    for (int k = 0; k < h_matrix.bits; k++) {
      ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard = codeword_in_hard.c[j].b[k];
      ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = 0;
      ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = 0;
    }
  }
  s_ldpc_decoder_output ldpc_decoder_output =
      f_ldpc_decode(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);
  printf("print decoder output\n");
  f_print_hard_codeword(ldpc_decoder_output.corrected_codeword, h_matrix.cols,
                        h_matrix.bits);
  f_print_ldpc_decoder_output(ldpc_decoder_output);

  // analyze soft bits quality
  printf("print soft data\n");
  s_hard_codeword codeword_in_soft = read_data_from_file(
      "soft_data_0527.txt", h_matrix.cols, bytes_of_userdata, bytes_of_parity);
  f_print_hard_codeword(codeword_in_soft, h_matrix.cols, h_matrix.bits);

  int error_cnt = 0;
  int cnt00 = 0, cnt01 = 0, cnt10 = 0, cnt11 = 0;
  int softbit_weak = 0, softbit_strong = 0;
  for (int j = 0; j < h_matrix.cols; j++) {
    for (int k = 0; k < h_matrix.bits; k++) {
      if (codeword_in_hard.c[j].b[k] !=
          ldpc_decoder_output.corrected_codeword.c[j].b[k]) {
        error_cnt++;
        if (codeword_in_soft.c[j].b[k] == 1)
          cnt11++;
        else
          cnt10++;
      } else {
        if (codeword_in_soft.c[j].b[k] == 1)
          cnt01++;
        else
          cnt00++;
      }
      if (codeword_in_soft.c[j].b[k] == 1)
        softbit_weak++;
      else
        softbit_strong++;
    }
  }

  printf(" Hard bit error cnt             : %d\n", error_cnt);
  printf(" Hard bit error, soft bit strong: %d\n", cnt10);
  printf(" Hard bit error, soft bit weak  : %d\n", cnt11);
  printf(" Hard bit noerr, soft bit strong: %d\n", cnt00);
  printf(" Hard bit noerr, soft bit weak  : %d\n", cnt01);
  printf(" Total soft bit weak cnt        : %d\n", softbit_weak);
  printf(" Total soft bit strong cnt      : %d\n", softbit_strong);
}

void test_custom_data_soft() {
  s_ldpc_decoder_parameters ldpc_decoder_parameters;
  ldpc_decoder_parameters.post_process_en = 1;
  ldpc_decoder_parameters.syndrome_weight_thr_qc = 16;
  ldpc_decoder_parameters.syndrome_weight_thr_post = 48;
  ldpc_decoder_parameters.likelihood_thr = 128;
  int likelihood_init[8][4] = {{0, 4, 3, 3}, {0, 3, 3, 0}, {0, 3, 3, 2},
                               {0, 3, 2, 2}, {0, 2, 2, 2}, {0, 5, 4, 3},
                               {0, 2, 1, 1}, {0, 1, 1, 1}};
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 4; j++)
      ldpc_decoder_parameters.likelihood_init_coef_all[i][j] =
          likelihood_init[i][j];

  ldpc_decoder_parameters.post_ratio = 12;
  ldpc_decoder_parameters.likelihood_map[0] = 0;
  ldpc_decoder_parameters.likelihood_map[1] = 1;
  ldpc_decoder_parameters.likelihood_map[2] = 2;
  ldpc_decoder_parameters.likelihood_map[3] = 3;

  ldpc_decoder_parameters.early_terminate_dis = 1;
  ldpc_decoder_parameters.questionable_sense = 1;
  int bytes_of_parity = 672, bytes_of_userdata = 4112;
  int cols = (bytes_of_parity / 64) + (bytes_of_userdata / 64);
  if (bytes_of_parity % 64 != 0)
    cols++;
  if (bytes_of_userdata % 64 != 0)
    cols++;

  s_h_matrix h_matrix = f_h_matrix(bytes_of_userdata, bytes_of_parity);
  s_hard_codeword codeword_in_hard = read_data_from_file(
      "hard_data_0", h_matrix.cols, bytes_of_userdata, bytes_of_parity);
  f_print_hard_codeword(codeword_in_hard, h_matrix.cols, h_matrix.bits);

  printf("print soft data\n");
  s_hard_codeword codeword_in_soft = read_data_from_file(
      "soft_data_0527.txt", h_matrix.cols, bytes_of_userdata, bytes_of_parity);
  f_print_hard_codeword(codeword_in_soft, h_matrix.cols, h_matrix.bits);

  // hard decoding
  s_ldpc_decoder_input ldpc_decoder_input;

  ldpc_decoder_input.nand_strobes = 3;
  ldpc_decoder_input.soft_bits = 1;
  ldpc_decoder_input.post_iteration = 1024;
  ldpc_decoder_input.iteration_limit = 4048;
  ldpc_decoder_input.syndrome_cal_only = 0;
  s_codeword corrupted_codeword;
  for (int j = 0; j < h_matrix.cols; j++) {
    for (int k = 0; k < h_matrix.bits; k++) {
      ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_hard =
          codeword_in_hard.c[j].b[k];
      ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable = codeword_in_soft.c[j].b[k];
      ldpc_decoder_input.corrupted_codeword.c[j].b[k].bit_questionable2 = codeword_in_soft.c[j].b[k];
    }
  }
  s_ldpc_decoder_output ldpc_decoder_output =
      f_ldpc_decode(ldpc_decoder_input, ldpc_decoder_parameters, h_matrix);
  printf("print decoder output\n");
  f_print_hard_codeword(ldpc_decoder_output.corrected_codeword, h_matrix.cols,
                        h_matrix.bits);
  f_print_ldpc_decoder_output(ldpc_decoder_output);
}

s_decode_res test_random_data(float rber, long number_of_codewords) {
  s_ldpc_decoder_parameters ldpc_decoder_parameters;
  ldpc_decoder_parameters.post_process_en = 1;
  ldpc_decoder_parameters.VN_BITS = 8;
  ldpc_decoder_parameters.syndrome_weight_thr_qc = 16;
  ldpc_decoder_parameters.syndrome_weight_thr_post = 48;
  ldpc_decoder_parameters.likelihood_thr = 128;
  s_decode_res res;
  int likelihood_init[8][4] = {{0, 7, 6, 5}, {0, 7, 5, 5}, {0, 6, 5, 4},
                               {0, 5, 5, 4}, {0, 4, 4, 3}, {0, 2, 2, 2},
                               {0, 4, 3, 2}, {0, 4, 3, 2}};
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 4; j++)
      ldpc_decoder_parameters.likelihood_init_coef_all[i][j] =
          likelihood_init[i][j];

  ldpc_decoder_parameters.post_ratio = 12;
  ldpc_decoder_parameters.likelihood_map[0] = 0;
  ldpc_decoder_parameters.likelihood_map[1] = 1;
  ldpc_decoder_parameters.likelihood_map[2] = 2;
  ldpc_decoder_parameters.likelihood_map[3] = 3;

  ldpc_decoder_parameters.early_terminate_dis = 1;
  ldpc_decoder_parameters.questionable_sense = 1;
  int bytes_of_parity = 448, bytes_of_userdata = 4096;
  s_ldpc_decoder_input ldpc_decoder_input;
  ldpc_decoder_input.nand_strobes = 7;

  if (ldpc_decoder_input.nand_strobes != 0)
    ldpc_decoder_parameters.VN_BITS = 8;
  ldpc_decoder_input.soft_bits = ldpc_decoder_input.nand_strobes == 0   ? 0
                                 : ldpc_decoder_input.nand_strobes == 3 ? 1
                                                                        : 2;

  ldpc_decoder_input.post_iteration = input_parameter.post_iter;
  ldpc_decoder_input.iteration_limit = input_parameter.max_iter;
  ldpc_decoder_input.syndrome_cal_only = 0;

  if (ldpc_decoder_input.nand_strobes == 5) {
    ldpc_decoder_parameters.likelihood_init_fraction[1] = 6;
    ldpc_decoder_parameters.likelihood_init_fraction[2] = 6;
  } else if (ldpc_decoder_input.nand_strobes == 7) {
    ldpc_decoder_parameters.likelihood_init_fraction[1] = 8;
    ldpc_decoder_parameters.likelihood_init_fraction[2] = 4;
  } else {
    ldpc_decoder_parameters.likelihood_init_fraction[1] =0;
    ldpc_decoder_parameters.likelihood_init_fraction[2] = 0;
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  s_h_matrix h_matrix = f_h_matrix(bytes_of_userdata, bytes_of_parity);
  if (rank == 0) {
    printf("Print decoding parameters ... \n");
    printf("    VN_BITS                       = %d \n", ldpc_decoder_parameters.VN_BITS);
    printf("    Nand strobes                  = %d \n", ldpc_decoder_input.nand_strobes);
    printf("    iter, post_iter               = %d %d \n", ldpc_decoder_input.iteration_limit, ldpc_decoder_input.post_iteration);
    printf("    user_data_bytes, parity_bytes = %d %d \n", bytes_of_userdata, bytes_of_parity);
    printf("    rber                          = %f \n", rber);
    printf("    post_ratio                    = %d \n", ldpc_decoder_parameters.post_ratio);
    printf("    numb_of_cw for each processor = %ld \n", number_of_codewords);
    int likelihood_use_index;
    if (h_matrix.rows == 5)
      likelihood_use_index = 0;
    else
      likelihood_use_index = h_matrix.rows - 6;

    printf("likelihood_init_coef = %d %d %d %d\n",
           ldpc_decoder_parameters.likelihood_init_coef_all[likelihood_use_index][0],
           ldpc_decoder_parameters.likelihood_init_coef_all[likelihood_use_index][1],
           ldpc_decoder_parameters.likelihood_init_coef_all[likelihood_use_index][2],
           ldpc_decoder_parameters.likelihood_init_coef_all[likelihood_use_index][3]);
    printf("Print decoding parameters end\n\n\n");
  }
  int no_err_inj = 0;
  res = encode_errinj_decode(h_matrix, ldpc_decoder_parameters,
                             ldpc_decoder_input, no_err_inj, rber,
                             number_of_codewords);

  return res;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // printf("MPI process d of d\n",rank,size);
  auto start = std::chrono::high_resolution_clock::now();

  std::time_t t = std::time(nullptr);
  char buffer[100];
  if (rank == 0) {
    std::strftime(buffer, sizeof(buffer),
                  "\n Current time: %Y-%m-%d %H:%M:%S", std::localtime(&t));
    std::cout << buffer << '\n';
  }
  float rber = 0.0108;
  long int total_codewords = 20; // number_of_codewords * size;

  input_parameter.num_errors = 0;
  input_parameter.max_iter = 1024;
  input_parameter.post_iter = 512;
  input_parameter.max_failures = 1000;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--max_failures")
      input_parameter.max_failures = atof(argv[i + 1]);
    else if (arg == "--num_codewords")
      total_codewords = atof(argv[i + 1]);
    else if (arg == "--rber")
      rber = atof(argv[i + 1]);
    else if (arg == "--num_errors")
      input_parameter.num_errors = atof(argv[i + 1]);
    else if (arg == "--bytes_of_parity")
      input_parameter.bytes_of_parity = atof(argv[i + 1]);
    // else printf("illegal parameter,will exit\n");exit(1);
  }

  iter_stat = new int[input_parameter.max_iter + 1];
  for (int i = 0; i < input_parameter.max_iter + 1; i++)
    iter_stat[i] = 0;
  long int number_of_codewords = total_codewords / size;
  if (total_codewords - number_of_codewords * size != 0 && rank == 0) {
    printf("total codes can not divede number of processes\n");
    exit(1);
  }

  s_decode_res res;

  long int total_fails = 0, total_iter = 0, total_syndrome_weight_before = 0,
           total_errors_in_codeword = 0, actual_codewords;
  float fail_rate;

  res = test_random_data(rber, number_of_codewords);
  if (rank == 0)
    printf("test finished, collecting result\n");
  MPI_Reduce(&res.fail_cnt, &total_fails, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&res.total_iter, &total_iter, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&res.syndrome_weight_before, &total_syndrome_weight_before, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&res.errors_in_codeword, &total_errors_in_codeword, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&res.actual_codewords, &actual_codewords, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  // MPI Reduce(iter_stat,iter _stat,input
  // _parameter.max_iter+1,MPI_INT,MPI_SUM,0,MPI _COMM_WORLD);

  if (rank == 0)
    MPI_Reduce(MPI_IN_PLACE, iter_stat, input_parameter.max_iter + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce(iter_stat, NULL, input_parameter.max_iter + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    fail_rate = total_fails * 1.0 / (actual_codewords * 1.0);
    printf("There are %ld of fails of %ld codewords, failure rate=%f\n",
           total_fails, actual_codewords, fail_rate);
    printf("Total iterations= %ld, average iterations=%f\n", total_iter,
           (float)total_iter * 1.0 / actual_codewords);
    printf("Total syndrome weight= %ld, average syndrome weight=%f\n",
           total_syndrome_weight_before,
           (float)total_syndrome_weight_before * 1.0 / actual_codewords);
    printf("Total errors in codeword= %ld, average errors in codeword=%f\n",
           total_errors_in_codeword,
           (float)total_errors_in_codeword * 1.0 / actual_codewords);
    printf("Print iteration statitics...\n");
    long long int sum = 0;
    for (int ii = 0; ii < input_parameter.max_iter + 1; ii++) {
      if (iter_stat[ii] != 0)
        printf("iter_stat[%4d] = %d\n", ii, iter_stat[ii]);
      sum += iter_stat[ii];
    }
    if (sum != total_codewords)
      printf("iter statitics wrong, actual=%lld, exp=%ld\n", sum,
             total_codewords);
    else
      printf("iter statitics correct, actual=%lld, exp=%ld\n", sum,
             total_codewords);
    printf("Print iteration statitics done!\n");
  }
  delete(iter_stat);

  MPI_Finalize();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsed = end - start;
  if (rank == 0)
    std::cout << "Elapsed time:" << elapsed.count() << "ms\n";
  t = std::time(nullptr);
  if (rank == 0) {
    std::strftime(buffer, sizeof(buffer),
                  "\n\n\n Current time: %Y-%m-%d %H:%M:%S", std::localtime(&t));
    std::cout << buffer << '\n';
  }
  return 0;
}
