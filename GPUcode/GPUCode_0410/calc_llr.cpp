// Approximate inverse Q-function
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

int Quantize(double x, int m, int k) {
    int t, t1;
    t = 1 << k;
    t = (int)(x * t + (x > 0 ? 0.5 : -0.5));
    t1 = (1 << m) - 1;
    t1 = t1 & t;
    return t1;
}

void Tru2ints(int x, int m, int *dx) {
    int t1 = x;
    int t2 = -1;
    t1 = (x >> (m - 1));
    if (t1 == 0)
        (*dx) = x;
    else {
        t1 = (1 << m) - 1;
        t2 = t2 ^ t1; // xor operation
        *dx = t2 | x;
    }
}

double Sat_Quan(double x, double Max_V, double Min_V, int m, int k) {
    double val;
    int q, p;
    val = x;
    val = (val > Max_V) ? Max_V : val;
    val = (val < Min_V) ? Min_V : val;
    q = Quantize(val, m, k);
    Tru2ints(q, m, &p);
    return (double)p / std::pow(2.0, k);
}

double inv_norm_cdf(double p) {
    assert(p > 0.0 && p < 1.0);
    // Coefficients
    const double a1 = -3.969683028665376e+01;
    const double a2 = 2.209460984245205e+02;
    const double a3 = -2.759285104469687e+02;
    const double a4 = 1.383577518672690e+02;
    const double a5 = -3.066479806614716e+01;
    const double a6 = 2.506628277459239e+00;
    const double b1 = -5.447609879822406e+01;
    const double b2 = 1.615858368580409e+02;
    const double b3 = -1.556989798598866e+02;
    const double b4 = 6.680131188771972e+01;
    const double b5 = -1.328068155288572e+01;
    const double c1 = -7.784894002430293e-03;
    const double c2 = -3.223964580411365e-01;
    const double c3 = -2.400758277161838e+00;
    const double c4 = -2.549732539343734e+00;
    const double c5 = 4.374664141464968e+00;
    const double c6 = 2.938163982698783e+00;
    const double d1 = 7.784695709041462e-03;
    const double d2 = 3.224671290700398e-01;
    const double d3 = 2.445134137142996e+00;
    const double d4 = 3.754408661907416e+00;
    const double plow = 0.02425;
    const double phigh = 1.0 - plow;
    double q, r;
    if (p < plow) {
        q = std::sqrt(-2.0 * std::log(p));
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }

    if (p > phigh) {
        q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }

    q = p - 0.5;
    r = q * q;
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / ((((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r) + 1.0);
}

double Qinv(double rber) {
    return inv_norm_cdf(1.0 - rber);
}

double rber_to_sigma(double rber) {
    double qinv = Qinv(rber);
    return 1.0 / qinv;
}

void ch_llr_alloc(int sd_bit, const float *vref_in, int llr_tot_bit, int llr_frac_bit, float hd0_llr, float hd1_llr, float rber, float llr_tbl_return[8], float &awgn_sigma) {
    int sd_num = sd_bit;
    int rd_num = sd_num;
    int bin_num = rd_num + 1;
    int max_llr_bin = 0;
    awgn_sigma = (float)rber_to_sigma((double)rber);

    if (sd_bit <= 1) {
        for (int i = 0; i < 4; i++)
            llr_tbl_return[i] = hd0_llr;
        for (int i = 4; i < 8; i++)
            llr_tbl_return[i] = hd1_llr;
        printf("final LLR table=%f %f %f %f %f %f %f %f\n", llr_tbl_return[0], llr_tbl_return[1], llr_tbl_return[2], llr_tbl_return[3], llr_tbl_return[4], llr_tbl_return[5], llr_tbl_return[6], llr_tbl_return[7]);
        return;
    }

    float *vref = (float *)std::calloc(rd_num, sizeof(*vref));
    for (int i = 0; i < rd_num; i++)
        vref[i] = vref_in[i];
    float *llr_tbl = (float *)std::calloc(bin_num, sizeof(*llr_tbl));
    char *bin_split = (char *)std::calloc(rd_num, sizeof(*bin_split));
    char *bin_asc_ord = (char *)std::calloc(rd_num + 1, sizeof(*bin_asc_ord));
    float *llr_asc_ord = (float *)std::calloc(rd_num + 1, sizeof(*llr_asc_ord));
    float *vref_asc_ord = (float *)std::calloc(rd_num, sizeof(*vref_asc_ord));
    char *sd_asc_ord = (char *)std::calloc(rd_num + 1, sizeof(*sd_asc_ord));
    int *bin_distr = (int *)std::calloc(bin_num, sizeof(*bin_distr));

    float awgn_sigma_sqaure = awgn_sigma * awgn_sigma;

    int llr_tot_num = llr_tot_bit;
    int llr_frac_num = llr_frac_bit;
    float llr_max = (std::pow(2.0f, llr_tot_num - 1) - 1.0f) / std::pow(2.0f, llr_frac_num);
    float llr_min = (1.0f - std::pow(2.0f, llr_tot_num - 1)) / std::pow(2.0f, llr_frac_num);
    // Step 1. Decide the bin to be splitted and determine bin indices in
    // ascending Vref order.
    for (int rd_indx = 0; rd_indx < rd_num; rd_indx++) 
    {
        if (rd_indx == 0)
        {
            bin_split[0] = 0;
            vref_asc_ord[0] = vref[0];
            bin_asc_ord[0] = 1;
            bin_asc_ord[1] = 0;
        }
        else if (vref[rd_indx] > vref_asc_ord[rd_indx - 1]) {
                bin_split[rd_indx] = bin_asc_ord[rd_indx];
                bin_asc_ord[rd_indx + 1] = bin_asc_ord[rd_indx];
                bin_asc_ord[rd_indx] = rd_indx + 1;
                vref_asc_ord[rd_indx] = vref[rd_indx];
        } 
        else {
            for (int i = 0; i < rd_indx; i++) {
                if (vref[rd_indx] < vref_asc_ord[i]) {
                    bin_split[rd_indx] = bin_asc_ord[i];

                    for (int j = rd_indx; j >= i; j--)
                        bin_asc_ord[j + 1] = bin_asc_ord[j];
                    bin_asc_ord[i] = rd_indx + 1;

                    for (int j = rd_indx - 1; j >= i; j--)
                        vref_asc_ord[j + 1] = vref_asc_ord[j];
                    vref_asc_ord[i] = vref[rd_indx];

                    break;
                }
            }
        }

        // Step 2. Determine LLR in ascending order of Vref.
        if (rd_num == 1) {
            llr_asc_ord[0] = hd1_llr;
            llr_asc_ord[1] = hd0_llr;
        } 
        else
        {
            llr_asc_ord[0] = (float)Sat_Quan((double)((-1.0f + vref_asc_ord[0]) / awgn_sigma_sqaure), llr_max, llr_min, llr_tot_num, llr_frac_num);
            llr_asc_ord[rd_num] = (float)Sat_Quan((double)((1.0f + vref_asc_ord[rd_num - 1]) / awgn_sigma_sqaure), llr_max, llr_min, llr_tot_num, llr_frac_num);

            for (int i = 1; i < rd_num; i++) 
                llr_asc_ord[i] = (float)Sat_Quan((double)((vref_asc_ord[i] + vref_asc_ord[i - 1]) / awgn_sigma_sqaure), llr_max, llr_min, llr_tot_num, llr_frac_num);
        }
    

        // Step 3. Determine soft data from NAND.
        int **nand_read;

        nand_read = (int **)std::calloc(rd_num + 1, sizeof(*nand_read));
        for (int i = 0; i <= rd_num; i++)
            nand_read[i] = (int *)std::calloc(rd_num, sizeof(*nand_read[i]));

        for (int i = 0; i < rd_num; i++) {
            for (int j = 0; j < rd_num; j++) {
                nand_read[i][j] = (vref_asc_ord[i] <= vref[j]) ? 1 : 0;
            }
        }

        for (int j = 0; j < rd_num; j++)
            nand_read[rd_num][j] = 0;
        for (int i = 0; i <= rd_num; i++) {
            sd_asc_ord[i] = nand_read[i][0];
            for (int j = 1; j < rd_num; j++)
                sd_asc_ord[i] = sd_asc_ord[i] * 2 + nand_read[i][j];
        }

        for (int i = 0; i <= rd_num; i++)
            std::free(nand_read[i]);
        std::free(nand_read);

        // Step 4. Determine LLR table.
        for (int i = 0; i < bin_num; i++) {
            for (int j = 0; j <= rd_num; j++) {
                if (bin_asc_ord[j] == i)
                    llr_tbl[i] = llr_asc_ord[j];
            }
        }

        max_llr_bin = bin_asc_ord[rd_num];
    }
    
    llr_tbl_return[0] = llr_tbl[0];
    llr_tbl_return[1] = llr_tbl[6];
    llr_tbl_return[2] = llr_tbl[4];
    llr_tbl_return[3] = llr_tbl[2];
    llr_tbl_return[4] = llr_tbl[7];
    llr_tbl_return[5] = llr_tbl[5];
    llr_tbl_return[6] = llr_tbl[3];
    llr_tbl_return[7] = llr_tbl[1];
  
    printf("final LLR table=%f %f %f %f %f %f %f %f\n", llr_tbl_return[0], llr_tbl_return[1], llr_tbl_return[2], llr_tbl_return[3], llr_tbl_return[4], llr_tbl_return[5], llr_tbl_return[6], llr_tbl_return[7]);
}
