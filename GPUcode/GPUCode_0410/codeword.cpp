#include <iostream>
#include <iomanip>
#include "codeword.h"

void codeword::clear(void) {
    for (auto &col : cols) {
        for (auto &elem : col)
            elem = 0;
    }

    return;
}

void codeword::print(void) {
    for (int idx = 0; idx < MAX_COL_COUNT; idx++) {
        std::cout << "Col " << std::dec << std::setw(2) << std::setfill('0') << idx++ << ": ";
        for (int i = 0; i < 8; i++)
            std::cout << std::uppercase << std::hex << std::setw(16) << std::setfill('0') << cols[idx][i] << " ";
        std::cout << "\n";
    }

    std::cout << std::endl;
    return;
}

decoder_input_cw::decoder_input_cw() {
    clear();
    return;
}

void decoder_input_cw::clear(void) {
    soft0.clear();
    soft1.clear();
    hard.clear();
}
