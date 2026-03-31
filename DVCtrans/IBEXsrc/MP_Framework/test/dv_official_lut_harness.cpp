#include <filesystem>
#include <iostream>

#include "../../DVsrc/ldpc.h"

int main() {
    std::filesystem::create_directories("output");
    std::cout << "[dv-official-lut] generating output/ldpc_matrix.vh and output/ldpc_matrix.h" << std::endl;
    f_create_include_files();
    std::cout << "[dv-official-lut] done" << std::endl;
    return 0;
}
