#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include "ldpc_error_injector.h"

using namespace std;

ldpc_error_injector::ldpc_error_injector()
{
    int i;
}

ldpc_error_injector::~ldpc_error_injector()
{

}

extern "C"
void* ldpc_error_injector_new()
{
    return new ldpc_error_injector();
}

extern "C"
void ldpc_error_injector_inv (void* sv_inst,
                              svOpenArrayHandle sv_hardbit_data,
                              svOpenArrayHandle sv_softbit_data,
                              svOpenArrayHandle sv_softbit_data1,
                              const svOpenArrayHandle sv_hardbit_corrupted,
                              const svOpenArrayHandle sv_softbit0_corrupted,
                              const svOpenArrayHandle sv_softbit1_corrupted,
                              const svOpenArrayHandle sv_num_errors_in_codeword,
                              const svOpenArrayHandle sv_num_errors_in_userdata,
                              int user_data_bytes=0,
                              float rber=0.0,
                              int random_seed=0,
                              int nand_srobes=1)
{
    int soft_bits = 0;
    float random_number;
    s_rbers rbers;
    s_error_injected_bit error_injected_bit;
    int time_value = (unsigned int)time(NULL);
    int injected_errors_in_codeword = 0;
    int injected_errors_in_userdata = 0;
    const int userdata_bits = user_data_bytes * 8;

    // generate new dynamic arrays
    int size = svHigh(sv_hardbit_data, 1);
    unsigned char* hardbit_corrupted = new unsigned char[size + 1];
    unsigned char* softbit0_corrupted = new unsigned char[size + 1];
    unsigned char* softbit1_corrupted = new unsigned char[size + 1];

#ifdef _LDPC_DBG_ERR_INJ_IN_FIRST128
    const unsigned int seed = (random_seed == 0) ? (unsigned int)time_value : (unsigned int)random_seed;
    std::mt19937 gen(seed);
#ifndef LDPC_DBG_ERR_NUM
#define LDPC_DBG_ERR_NUM 100
#endif
    const int bit_len = size + 1;
    const int dbg_window_bits = (bit_len < 128) ? bit_len : 128;
    int dbg_err_num = LDPC_DBG_ERR_NUM;
    if (dbg_err_num < 0)
        dbg_err_num = 0;
    if (dbg_err_num > dbg_window_bits)
        dbg_err_num = dbg_window_bits;
    bool dbg_inject_mask[128] = {0};
    {
        int idx[128];
        for (int t = 0; t < dbg_window_bits; t++)
            idx[t] = t;
        for (int t = dbg_window_bits - 1; t > 0; t--) {
            std::uniform_int_distribution<int> pick(0, t);
            const int r = pick(gen);
            const int tmp = idx[t];
            idx[t] = idx[r];
            idx[r] = tmp;
        }
        for (int t = 0; t < dbg_err_num; t++)
            dbg_inject_mask[idx[t]] = true;
    }
#else
    std::random_device rd;
    std::mt19937 gen(rd());
#endif
    std::uniform_real_distribution<> dis(0.0, 1.0);
    if (random_seed == 0)
        srand(time_value);
    else
        srand(random_seed);

    rbers = f_find_error_probabilities(nand_srobes, rber);
    printf("rbers[0]=%f\n", rbers.rber[0]);
    printf("rbers[1]=%f\n", rbers.rber[1]);
    printf("rbers[2]=%f\n", rbers.rber[2]);

    printf("### LDPC ERROR INJECTOR C++: SIZE: %8d USER DATA BYTES: %4d\n", size, user_data_bytes);
    for (int i = 0; i <= size; i++)
    {
        random_number = dis(gen);
#ifdef _LDPC_DBG_ERR_INJ_IN_FIRST128
        if (i < dbg_window_bits) {
            const auto soft_sample = f_inject_error(nand_srobes, rbers, random_number);
            error_injected_bit.hard = dbg_inject_mask[i];
            error_injected_bit.soft0 = (nand_srobes > 1) ? soft_sample.soft0 : 0;
            error_injected_bit.soft1 = (nand_srobes > 3) ? soft_sample.soft1 : 0;
        } else {
            error_injected_bit.hard = 0;
            error_injected_bit.soft0 = 0;
            error_injected_bit.soft1 = 0;
        }
#else
        error_injected_bit = f_inject_error(nand_srobes, rbers, random_number);
        if (nand_srobes <= 1)
            error_injected_bit.soft0 = 0;
        if (nand_srobes <= 3)
            error_injected_bit.soft1 = 0;
#endif
        if (error_injected_bit.hard) {
            injected_errors_in_codeword++;
            if (i < userdata_bits)
                injected_errors_in_userdata++;
        }
        const unsigned char hard_in = *((unsigned char*)svGetArrElemPtr(sv_hardbit_data, i));
        const unsigned char soft0_in = (nand_srobes > 1) ? *((unsigned char*)svGetArrElemPtr(sv_softbit_data, i)) : 0;
        const unsigned char soft1_in = (nand_srobes > 3) ? *((unsigned char*)svGetArrElemPtr(sv_softbit_data1, i)) : 0;
        hardbit_corrupted[i] = error_injected_bit.hard ^ hard_in;
        softbit0_corrupted[i] = error_injected_bit.soft0 ^ soft0_in;
        softbit1_corrupted[i] = error_injected_bit.soft1 ^ soft1_in;
        unsigned char* hardbit_temp = (unsigned char*)svGetArrElemPtr(sv_hardbit_corrupted, i);
        *hardbit_temp = hardbit_corrupted[i];
        unsigned char* softbit0_temp = (unsigned char*)svGetArrElemPtr(sv_softbit0_corrupted, i);
        *softbit0_temp = softbit0_corrupted[i];
        unsigned char* softbit1_temp = (unsigned char*)svGetArrElemPtr(sv_softbit1_corrupted, i);
        *softbit1_temp = softbit1_corrupted[i];
    }

    int *c_num_errors_in_codeword = (int*)svGetArrElemPtr(sv_num_errors_in_codeword, 0);
    *c_num_errors_in_codeword = injected_errors_in_codeword;
    int *c_num_errors_in_userdata = (int*)svGetArrElemPtr(sv_num_errors_in_userdata, 0);
    *c_num_errors_in_userdata = injected_errors_in_userdata;
}                              
