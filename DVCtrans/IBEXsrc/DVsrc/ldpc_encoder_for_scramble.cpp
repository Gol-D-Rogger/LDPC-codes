#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "ldpc_encoder.h"
#include "ldpc_matrix_inverse.h"

int main()
{
    char encoded_codeword_file_name[128];
    char userdata_file_name[128];

    int user_data_bytes = 4176;
    int parity_bytes = 4608 - user_data_bytes;

    s_h_matrix h_matrix = f_h_matrix(user_data_bytes, parity_bytes);

    sprintf(userdata_file_name, "scrambler_output.txt");
    sprintf(encoded_codeword_file_name, "ldpc_encoded_output.txt");
    f_ldpc_encode_for_scrambler(userdata_file_name, user_data_bytes, parity_bytes, encoded_codeword_file_name);

    return 0;
}