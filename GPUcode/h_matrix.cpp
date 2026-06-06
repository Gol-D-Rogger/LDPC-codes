#include <cstdio>
#include <cstdlib>

#include "h_matrix.h"

namespace {
void zero_h_matrix(h_matrix *matrix) {
    matrix->rows = 0;
    matrix->cols = 0;
    matrix->bits = LDPC_P;
    matrix->column_weight = 0;
    matrix->bytes_of_userdata = 0;
    matrix->bytes_of_parity = 0;
    matrix->extra_bytes_of_parity = 0;
    matrix->extra_bytes_of_userdata = 0;
    matrix->unused_bytes_of_parity = 0;
    matrix->unused_bytes_of_userdata = 0;
    matrix->extra_bits_of_parity = 0;
    matrix->extra_bits_of_userdata = 0;
    matrix->min_rows = LDPC_L;
    matrix->max_rows = LDPC_M;
    matrix->bits_in_last_column = 0;

    for (int i = 0; i < LDPC_M; i++) {
        matrix->delta[i] = 0;
        matrix->first_element[i] = 0;
        matrix->last_element[i] = 0;
        matrix->wraparound[i] = 0;
        matrix->wrap_base[i] = 0;
        matrix->wrap_num_deltas[i] = 0;
        matrix->row_weight[i] = 0;
        for (int j = 0; j < LDPC_N; j++) {
            matrix->element[i][j] = -1;
            matrix->occupied[i][j] = false;
            matrix->fade[i][j] = false;
            matrix->operational_h_matrix[i][j] = 0;
        }
    }

    for (int j = 0; j < LDPC_N; j++) {
        matrix->col_weight[j] = 0;
        matrix->parity_column[j] = false;
        for (int k = 0; k < LDPC_P; k++)
            matrix->last_row_active_bits[j][k] = false;
    }

    for (int i = 0; i < 2 * LDPC_N; i++)
        matrix->tile_ranges[i] = {0, 0};
}

bool base_use_location(int i, int j) {
    static const int avoid_row[50] = {
        0, 1, 2, 3, 4, 1, 3, 4, 2, 0, 4, 3, 1, 0, 2, 3, 1, 2, 0, 4, 0, 0, 4, 4, 2, 2, 3, 3, 1, 1, 3, 0, 3, 1, 4, 0, 4, 2, 2, 1, 1, 4, 3, 2, 0, 2, 4, 0, 1, 3,
    };

    if (j < 50)
        return i != avoid_row[j];

    int k = LDPC_N - 1 - j;
    if (j <= 59)
        return ((k + 4) % LDPC_L) != i;

    return ((k + 2) % LDPC_L) != i;
}
} // namespace

h_matrix::h_matrix() {
    zero_h_matrix(this);
    std::printf("*** Un-initialized H matrix constructed ***\n");
    return;
}

h_matrix::h_matrix(const int bytes_of_userdata, const int bytes_of_parity, int mx_cnfg) {
    zero_h_matrix(this);

    int VERBOSITY = 0;
    int i;
    int j;
    int k;
    int bit;
    int n;
    int h_matrix_index[LDPC_M];
    int matrix_sel = 0;

    this->bits = (bytes_of_userdata > 5000) ? 1024 : 512;
    this->rows = (this->bits == 1024) ? ((bytes_of_parity + 127) >> 7) : ((bytes_of_parity + 63) >> 6);
    this->cols = (this->bits == 1024) ? ((bytes_of_userdata + 127) >> 7) : ((bytes_of_userdata + 63) >> 6);
    this->cols += this->rows;

    this->bytes_of_userdata = bytes_of_userdata;
    this->bytes_of_parity = bytes_of_parity;
    this->min_rows = LDPC_L;
    this->max_rows = LDPC_M;

    this->unused_bytes_of_parity = (this->rows * (this->bits >> 3)) - this->bytes_of_parity;
    if (this->unused_bytes_of_parity != 0)
        this->extra_bytes_of_parity = (this->bits >> 3) - this->unused_bytes_of_parity;
    else
        this->extra_bytes_of_parity = 0;

    this->unused_bytes_of_userdata = ((this->cols - this->rows) * (this->bits >> 3)) - this->bytes_of_userdata;
    if (this->unused_bytes_of_userdata != 0)
        this->extra_bytes_of_userdata = (this->bits >> 3) - this->unused_bytes_of_userdata;
    else
        this->extra_bytes_of_userdata = 0;

    this->extra_bits_of_parity = this->extra_bytes_of_parity << 3;
    this->extra_bits_of_userdata = this->extra_bytes_of_userdata << 3;
    this->bits_in_last_column = this->extra_bits_of_userdata;

    if (VERBOSITY > 0)
        std::printf("#MATRIX: BITS: %4d ROWS: %2d COLS: %3d BYTES_OF_USERDATA: %5d BYTES_OF_PARITY: %4d EXTRA_BYTES_OF_USERDATA: %5d EXTRA_BYTES_OF_PARITY: %4d\n", this->bits, this->rows, this->cols, this->bytes_of_userdata, this->bytes_of_parity, this->extra_bytes_of_userdata, this->extra_bytes_of_parity);

    int rw[LDPC_M][LDPC_N];
    int rw_max[LDPC_M + 1];
    int rw_min[LDPC_M + 1];
    float rw_avg[LDPC_M + 1];
    bool ldpc_matrix_occupied[LDPC_M][LDPC_M][LDPC_N];
    bool ldpc_matrix_fade[LDPC_M][LDPC_M][LDPC_N];

    for (k = 0; k < LDPC_M; k++) {
        for (i = 0; i < LDPC_M; i++) {
            rw[k][i] = 0;
            for (j = 0; j < LDPC_N; j++) {
                ldpc_matrix_occupied[k][i][j] = 0;
                ldpc_matrix_fade[k][i][j] = 0;
            }
        }
    }

    for (i = LDPC_L - 1; i < LDPC_M; i++) {
        rw_avg[i] = (float)(4 * (LDPC_U + i + 1)) / (float)(i + 1);
        rw_min[i] = int(rw_avg[i]);
        rw_max[i] = int(rw_avg[i] + 0.999999);
        if (VERBOSITY > 0)
            std::printf("FOR MATRIX %2d, ROW_WEIGHT_AVG: %8.4f = %2d + %d/%d MIN:%2d MAX:%2d\n", i, rw_avg[i], int(rw_avg[i]), int(rw_avg[i] * (i + 1)) % (i + 1), (i + 1), rw_min[i], rw_max[i]);
    }

    srand(1);

    for (i = 0; i < LDPC_L; i++) {
        for (j = LDPC_N - LDPC_L; j < LDPC_N; j++) {
            k = LDPC_N - 1 - j;
            bool use_location = (k % LDPC_L) != i;
            if ((i == (LDPC_L - 1)) && (j == (LDPC_N - 1)))
                use_location = 0;

            if (use_location) {
                ldpc_matrix_occupied[LDPC_L - 1][i][j] = 1;
                rw[LDPC_L - 1][i]++;
            }
        }
    }

    for (i = 0; i < LDPC_L; i++) {
        for (j = 0; j < LDPC_U; j++) {
            bool use_location = base_use_location(i, j);
            if (use_location) {
                ldpc_matrix_occupied[LDPC_L - 1][i][j] = 1;
                rw[LDPC_L - 1][i]++;
            }
        }
    }

    if (VERBOSITY > 0) {
        for (i = 0; i < LDPC_M; i++) {
            k = LDPC_L - 1;
            std::printf("### MATRIX: %2d ROW:%2d WEIGHT: %2d OCCUPIED: ", k, i, rw[k][i]);
            for (j = 0; j < LDPC_N; j++)
                std::printf("%1x ", ldpc_matrix_occupied[k][i][j]);
            std::printf("\n");
        }
        std::printf("\n");
    }

    for (i = 0; i < LDPC_M; i++) {
        if (rw[LDPC_L - 1][i] > rw_max[LDPC_L - 1])
            rw_max[LDPC_L - 1] = rw[LDPC_L - 1][i];
    }

    for (k = LDPC_L; k < LDPC_M; k++) {
        for (i = 0; i < LDPC_M; i++) {
            rw[k][i] = 0;
            for (j = 0; j < LDPC_N; j++) {
                ldpc_matrix_occupied[k][i][j] = ldpc_matrix_occupied[k - 1][i][j];
                ldpc_matrix_fade[k][i][j] = ldpc_matrix_fade[k - 1][i][j];
                rw[k][i] += ldpc_matrix_occupied[k][i][j];
            }
        }

        for (i = 0; i < 5; i++) {
            ldpc_matrix_occupied[k][k - i][LDPC_N - k - 1] = 1;
            rw[k][k - i]++;
        }

        while (rw[k][k] < rw_min[k]) {
            i = (rand() % LDPC_U);
            if (ldpc_matrix_occupied[k][k][i] == 0) {
                ldpc_matrix_occupied[k][k][i] = 1;
                rw[k][k]++;
            }
        }

        for (j = 0; j < LDPC_N; j++) {
            if (ldpc_matrix_occupied[k][k][j]) {
                int rw_max_in_col = 0;
                for (i = 0; i < k; i++) {
                    if (ldpc_matrix_occupied[k][i][j] && (rw[k][i] > rw_max_in_col))
                        rw_max_in_col = rw[k][i];
                }

                bool found_it = 0;
                for (i = 0; i < k; i++) {
                    if (ldpc_matrix_occupied[k][i][j] && (rw[k][i] == rw_max_in_col) && !found_it) {
                        found_it = 1;
                        rw[k][i]--;
                        ldpc_matrix_occupied[k][i][j] = 0;
                        ldpc_matrix_fade[k][i][j] = 1;
                    }
                }
            }
        }

        if (VERBOSITY > 0) {
            for (i = 0; i < LDPC_M; i++) {
                std::printf("### MATRIX: %2d ROW: %2d WEIGHT: %2d OCCUPIED: ", k, i, rw[k][i]);
                for (j = 0; j < LDPC_N; j++)
                    std::printf("%1x ", ldpc_matrix_occupied[k][i][j]);
                std::printf("\n");
            }
            std::printf("\n");
        }
    }

    if (matrix_sel == 0) {
        this->delta[0] = 0;
        this->delta[1] = 13;
        this->delta[2] = 19;
        this->delta[3] = 29;
        this->delta[4] = 41;
        this->delta[5] = 67;
        this->delta[6] = 73;
        this->delta[7] = 79;
        this->delta[8] = 91;
        this->delta[9] = 97;
        this->delta[10] = 103;
        this->delta[11] = 111;
        this->delta[12] = 119;
    } else {
        this->delta[0] = 0;
        this->delta[1] = 293;
        this->delta[2] = 61;
        this->delta[3] = 479;
        this->delta[4] = 17;
        this->delta[5] = 173;
        this->delta[6] = 53;
        this->delta[7] = 277;
        this->delta[8] = 307;
        this->delta[9] = 229;
        this->delta[10] = 67;
        this->delta[11] = 271;
        this->delta[12] = 307;
    }

    h_matrix_index[0] = 0 + (3 * this->delta[0]);
    h_matrix_index[1] = 0 + (4 * this->delta[1]);
    h_matrix_index[2] = 0 + (4 * this->delta[2]);
    h_matrix_index[3] = 0 + (4 * this->delta[3]);
    h_matrix_index[4] = 0 + (3 * this->delta[4]);
    h_matrix_index[5] = 0;
    h_matrix_index[6] = 0;
    h_matrix_index[7] = 0;
    h_matrix_index[8] = 0;
    h_matrix_index[9] = 0;
    h_matrix_index[10] = 0;
    h_matrix_index[11] = 0;
    h_matrix_index[12] = 0;

    for (i = 0; i < this->rows; i++)
        this->last_element[i] = h_matrix_index[i];

    for (i = 0; i < this->rows; i++) {
        for (j = 0; j < this->cols; j++) {
            if (j < (this->cols - this->rows)) {
                this->fade[i][j] = ldpc_matrix_fade[this->rows - 1][i][j];
                this->occupied[i][j] = ldpc_matrix_occupied[this->rows - 1][i][j];
            } else {
                this->occupied[i][j] = ldpc_matrix_occupied[this->rows - 1][i][LDPC_N - this->cols + j];
                this->fade[i][j] = ldpc_matrix_fade[this->rows - 1][i][LDPC_N - this->cols + j];
            }
        }
    }

    for (i = 0; i < this->rows; i++)
        this->row_weight[i] = 0;

    for (j = 0; j < this->cols; j++)
        this->col_weight[j] = 0;

    for (j = 0; j < this->cols; j++) {
        n = this->cols - 1 - j;
        k = n;
        this->parity_column[j] = (j < this->rows);

        for (i = 0; i < this->rows; i++) {
            if (this->occupied[i][n] || this->fade[i][n]) {
                this->element[i][k] = h_matrix_index[i];
                this->first_element[i] = this->element[i][k];
                h_matrix_index[i] = (h_matrix_index[i] + this->bits - this->delta[i]) % this->bits;
                if (this->occupied[i][n]) {
                    this->row_weight[i]++;
                    this->col_weight[k]++;
                }
            } else {
                this->element[i][k] = -1;
            }
        }
    }

    for (i = 0; i < this->rows; i++)
        this->wraparound[i] = (this->bits + this->first_element[i] - this->last_element[i]) % this->bits;

    for (j = 0; j < LDPC_N; j++) {
        for (k = 0; k < LDPC_P; k++)
            this->last_row_active_bits[j][k] = 0;
    }

    int active_bits_of_parity = this->extra_bits_of_parity == 0 ? 512 : this->extra_bits_of_parity;
    for (j = 0; j < this->cols; j++) {
        if (this->occupied[this->rows - 1][j]) {
            for (k = 0; k < this->bits && k < LDPC_P; k++) {
                bit = (k + this->bits - this->element[this->rows - 1][j]) % this->bits;
                if (bit < active_bits_of_parity)
                    this->last_row_active_bits[j][k] = 1;
            }
        }
    }

    std::printf("[H_MATRIX] :: H-matrix with size %dB x %dB configured. [matrix_sel : %d]\n", bytes_of_userdata, bytes_of_parity, matrix_sel);

    FILE *fp_fade = NULL;
    FILE *fp_occu = NULL;
    FILE *fp_h = NULL;
    char pchk_file[100];
    char occupy_file[100];
    char fade_file[100];

    if (mx_cnfg == 0) {
        std::sprintf(pchk_file, "./matrice/matrix/LDPC_%dx%dex512_w4_dense5_QC_H.txt", this->rows, this->cols);
        std::sprintf(occupy_file, "./matrice/occupy/LDPC_%dx%dex512_w4_dense5_occupied.txt", this->rows, this->cols);
        std::sprintf(fade_file, "./matrice/fade/LDPC_%dx%dex512_w4_dense5_fade.txt", this->rows, this->cols);
    } else if (mx_cnfg == 1) {
        std::sprintf(pchk_file, "./rand_matrix/matrix/LDPC_%dx%dex512_w4_dense5_QC_H_30_0.txt", this->rows, this->cols);
        std::sprintf(occupy_file, "./rand_matrix/occupy/LDPC_%dx%dex512_w4_dense5_QC_H_occupied_30_0.txt", this->rows, this->cols);
        std::sprintf(fade_file, "./rand_matrix/fade/LDPC_%dx%dex512_w4_dense5_QC_H_fade_30_0.txt", this->rows, this->cols);
    }

    if ((mx_cnfg == 0) || (mx_cnfg == 1)) {
        std::printf("use customized matrix, mx_cnfg=%d\n", mx_cnfg);
        std::printf("fade_file= %s\n", fade_file);
        std::printf("occupy_file= %s\n", occupy_file);
        std::printf("element file= %s\n", pchk_file);

        fp_fade = std::fopen(fade_file, "r");
        fp_h = std::fopen(pchk_file, "r");
        fp_occu = std::fopen(occupy_file, "r");

        if (fp_fade && fp_h && fp_occu) {
            for (i = 0; i < this->rows; i++) {
                for (j = 0; j < this->cols; j++) {
                    int fade_val = 0;
                    int occupied_val = 0;
                    int element_val = -1;
                    std::fscanf(fp_fade, "%d", &fade_val);
                    std::fscanf(fp_occu, "%d", &occupied_val);
                    std::fscanf(fp_h, "%d", &element_val);
                    this->fade[i][j] = fade_val;
                    this->occupied[i][j] = occupied_val;
                    this->element[i][j] = element_val;
                }
            }
        }

        if (fp_fade)
            std::fclose(fp_fade);
        if (fp_h)
            std::fclose(fp_h);
        if (fp_occu)
            std::fclose(fp_occu);
    }

    this->setup_range_operational_matrix();
    this->print();
    std::printf("[H_MATRIX] :: Operational matrix with range info configured\n");
    return;
}

void h_matrix::print(void) {
    int i;
    int j;

    for (i = 0; i < this->rows; i++) {
        std::printf("%02d ", this->row_weight[i]);
        for (j = 0; j < this->cols; j++) {
            if (this->occupied[i][j] || this->fade[i][j])
                std::printf("%03X ", this->element[i][j]);
            else
                std::printf("...");
        }
        std::printf("\n");
    }

    std::printf("\n");

    for (i = 0; i < this->rows; i++) {
        std::printf("%02d,,", this->row_weight[i]);
        for (j = 0; j < this->cols; j++) {
            if (this->occupied[i][j])
                std::printf("x,");
            else if (this->fade[i][j])
                std::printf("y,");
            else
                std::printf(",");
        }
        std::printf("\n");
    }

    std::printf("\nOCCUPIED\n");
    for (i = 0; i < this->rows; i++) {
        std::printf("%02d,,", this->row_weight[i]);
        for (j = 0; j < this->cols; j++) {
            if (this->occupied[i][j])
                std::printf("1 ");
            else
                std::printf(".");
        }
        std::printf("\n");
    }

    std::printf("\nFADE\n");
    for (i = 0; i < this->rows; i++) {
        std::printf("%02d,,", this->row_weight[i]);
        for (j = 0; j < this->cols; j++) {
            if (this->fade[i][j])
                std::printf("1 ");
            else
                std::printf(".");
        }
        std::printf("\n");
    }

    std::printf("\n");
}

void h_matrix::setup_range_operational_matrix(void) {
    int i;
    int j;
    int k;
    bool do_not_use_this_bit = false;
    int tile_ranges_idx = 2;

    for (i = 0; i < LDPC_M; i++) {
        for (j = 0; j < LDPC_N; j++)
            operational_h_matrix[i][j] = 0;
    }

    for (j = 0; j < cols; j++) {
        for (i = 0; i < rows; i++) {
            if (!occupied[i][j] && !fade[i][j])
                continue;

            range curr_range = {0, 0};
            bool head_found = false;

            for (k = 0; k < bits && k < LDPC_P; k++) {
                do_not_use_this_bit = false;
                do_not_use_this_bit = do_not_use_this_bit || ((extra_bits_of_parity > 0) && (j == (cols - rows)) && (k >= extra_bits_of_parity));
                if (do_not_use_this_bit)
                    continue;

                bool bit_active = false;
                if (occupied[i][j] && (i < (rows - 1))) {
                    bit_active = true;
                } else if (occupied[i][j] && (i == (rows - 1)) && last_row_active_bits[j][k]) {
                    bit_active = true;
                } else if (fade[i][j] && !last_row_active_bits[j][k]) {
                    bit_active = true;
                }

                if (bit_active) {
                    curr_range.active_count++;
                    if (!head_found) {
                        head_found = true;
                        curr_range.offset = k;
                    }
                }

                if (head_found && !bit_active)
                    head_found = false;
            }

            if (curr_range.active_count == 0) {
                operational_h_matrix[i][j] = 0;
            } else if (curr_range.active_count == 512) {
                operational_h_matrix[i][j] = 1;
            } else {
                operational_h_matrix[i][j] = tile_ranges_idx;
                tile_ranges[tile_ranges_idx] = curr_range;
                tile_ranges_idx++;
            }
        }
    }

    range edge_case_null = {0, 0};
    tile_ranges[0] = edge_case_null;

    range edge_case_fill = {512, 0};
    tile_ranges[1] = edge_case_fill;

    return;
}
