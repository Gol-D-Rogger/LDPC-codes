class ldpc_encoder
{
public:

    int ROWS, COLS, BITS, CWBITS;

    int check_node[13][80];
    int variable_node[80][512];
    int syndrome[13][512];

public:
    ldpc_encoder();
    ~ldpc_encoder();

    void encoder_inv(unsigned char *userBits, unsigned char *encodeBits, unsigned int user_data_bytes, unsigned int parity_bytes);
    int check(int *decodedBits);
};
