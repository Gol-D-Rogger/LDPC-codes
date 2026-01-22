#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "ldpc_encoder.h"
#include "ldpc.h"

using namespace std;

ldpc_encoder::ldpc_encoder()
{
    int i;

    ROWS = 13;
    COLS = 80;
    BITS = 512;
    CWBITS = COLS*BITS;
}

ldpc_encoder::~ldpc_encoder()
{

}

void ldpc_encoder::encoder_inv(unsigned char *userBits, unsigned char *encodedBits, unsigned int user_data_bytes, unsigned int parity_bytes)
{
#include "ldpc_matrix_inverse.h"

    int DEBUG_MODE = 0;

    int i, j, k;
    int matrix_bits = 512;
    int matrix_rows = (parity_bytes + 63) >> 63;
    int matrix_cols = ((user_data_bytes + 63) >> 6) + matrix_rows;
    s_h_matrix h_matrix = f_h_matrix(user_data_bytes, parity_bytes);

    s_hard_codeword ldpc_encoder_input;
    int user_data_bits = user_data_bytes * 8;
    int parity_bits = parity_bytes * 8;
    int skip_parity_bits_in_first_parity_column = (h_matrix.cols * h_matrix.bits) - parity_bits;
    int last_parity_bit_in_first_parity_column = h_matrix.bits - skip_parity_bits_in_first_parity_column - 1;
    int first_parity = ((h_matrix.cols * h_matrix.bits) - parity_bits - 1);
    int bit_location = 0;
    int unpacked_bit_location = 0;

    for (j = 0; j < h_matrix.cols; j++)
    {
        for (i = 0; i < h_matrix.bits; i++)
        {
            if (bit_location < user_data_bits)
                ldpc_encoder_input.c[j].b[k] = userBits[bit_location++];
            else
                ldpc_encoder_input.c[j].b[k] = 0;
            if (DEBUG_MODE)
                ldpc_encoder_input.c[j].b[k] = 0;
        }
    }

    s_hard_codeword ldpc_encoder_output = f_ldpc_encode(ldpc_encoder_input, h_matrix);

    bit_location = 0;
    for (j = 0; j < h_matrix.cols; j++)
    {
        for (k = 0; k < h_matrix.bits; k++)
        {
            if (bit_location < user_data_bits)
            {
                encodedBits[bit_location++] = ldpc_encoder_output.c[j].b[k];
                unpacked_bit_location++;
            }
            else if (j == (h_matrix.cols - matrix_rows))
            {
                if (k < (512 - last_parity_bit_in_first_parity_column))
                    encodedBits[unpacked_bit_location++] = ldpc_encoder_output.c[j].b[k];
                unpacked_bit_location++;
            }
            else if (j > (h_matrix.cols - matrix_rows))
            {
                encodedBits[unpacked_bit_location++] = ldpc_encoder_output.c[j].b[k];
                unpacked_bit_location++;
            }
            else
            {
                unpacked_bit_location++;
            }
        }
    }
}

extern "C"
void* ldpc_encoder_new()
{
    return new ldpc_encoder();
}

extern "C"
void ldpc_encoder_delete(void* sv_inst) 
{
   ldpc_encoder* inst =  (ldpc_encoder*) sv_inst;
   inst->~ldpc_encoder();
}

extern "C"
void ldpc_encoder_encoder_inv(void* sv_inst, svOpenArrayHandle sv_userBits, svOpenArrayHandle sv_encodedBits, unsigned int user_data_bytes = 4096, unsigned int parity_bytes = 512)
{
    // Generate new dynamic array
    int size = svHigh(sv_userBits, 1);
    unsigned char* userBits = new unsigned char[size + 1];
    for (int i = 0; i <= size; i++) {
        userBits[i] = *((unsigned char*)svGetArrElemPtr(sv_userBits, i));
    }

    size = svHigh(sv_encodedBits, 1);
    unsigned char* encodedBits = new unsigned char[size + 1];
    ldpc_encoder* inst =  (ldpc_encoder*) sv_inst;
    inst->encoder_inv(userBits, encodedBits, user_data_bytes, parity_bytes);

    for (int i = 0; i < inst->CWBITS; i++) {
        unsigned char* encodedBits_tmp =  (unsigned char*)svGetArrElemPtr(sv_encodedBits, i);
        *encodedBits_tmp = encodedBits[i];
    }
}