#ifndef _GEN_LDPC_H
#define _GEN_LDPC_H

#include "mod2sparse.h"

typedef struct WD_pair
{
    int weight;
    int num;
} WD_pair;

typedef struct WD_vector
{
    int n;
    WD_pair *wd;
} WD_vector;

int gen_ldpc_matrix( mod2sparse *basePCH,   // place to store the base parity-check matrix 
                     mod2sparse *exPCH,     // place to store the expanded parity-check matrix
                     int rownum,            // number of rows in the base parity-check matrix
                     int colnum,            // number of columns in the base parity-check matrix
                     int ex_factor,         // expansion factor
                     int pad_bit,           // position of the padding bit in each row of the base matrix
                     int t_size,            // Identity matrix T size
                     WD_vector *row_dt,     // row weight distribution of the base matrix
                     WD_vector *col_dt,     // column weight distribution of the base matrix
                     int maximum_girth,     // maximum girth of the expanded matrix
                     int minimum_girth,     // minimum girth of the expanded matrix
                     int** mask_matrix );

// generate irregular LDPC codes files using back padding scheme
void gen_ldpc_files (char *codefile, char* maskfile, char *cyclefile, int filenum,
                     WD_vector *row_dt, WD_vector *col_dt,
                     int ex_factor, int pad_bit, int g,
                     int start_girth, int end_girth,
                     float expected_avg_girth);

// insert a sub matrix into a frame matrix
void insert_sub_matrix( mod2sparse *mainMatrix, 
                        int startRow,
                        int startCol,
                        int subDim,
                        int shift );

void insert_sub_matrix_row_mask(mod2sparse *mainMatrix, 
                                int startRow,
                                int startCol,
                                int mask_row_start,
                                int mask_row_end,
                                int subDim,
                                int shift );                        

void insert_sub_matrix_col_mask(mod2sparse *mainMatrix, 
                                int startRow,
                                int startCol,
                                int mask_col_start,
                                int mask_col_end,
                                int subDim,
                                int shift );
                                
// remove a sub matrix into a frame matrix
void remove_sub_matrix( mod2sparse *mainMatrix, 
                        int startRow,
                        int startCol,
                        int subDim,
                        int shift );                                

void remove_sub_matrix_row_mask(mod2sparse *mainMatrix, 
                                int startRow,
                                int startCol,
                                int mask_row_start,
                                int mask_row_end,
                                int subDim,
                                int shift );                
                                
void remove_sub_matrix_col_mask(mod2sparse *mainMatrix, 
                                int startRow,
                                int startCol,
                                int mask_col_start,
                                int mask_col_end,
                                int subDim,
                                int shift );
                                
// check the matrix PHI invertility
int check_PHI_inv( mod2sparse *H, int m, int n, int p, int g);
int check_PHI_inv_pad( mod2sparse *H, int m, int n, int p, int g, int pad_bit);                                
                     
// search cyc
int search_cyc4_mask(mod2sparse* basePCH, mod2entry* e, int** mask_matrix, int qc_size, int pad_bit);
int search_cyc6_mask(mod2sparse* basePCH, mod2entry* e, int** mask_matrix, int qc_size, int pad_bit);
int search_cyc4(mod2sparse* basePCH, mod2entry* e, int qc_size);
int search_cyc6(mod2sparse* basePCH, mod2entry* e, int qc_size);

int cyc4_cal_t(int** matrix, int** mask_matrix, int row, int col, int qc_size, int pad_bit);
int cyc6_cal_t(int** matrix, int** mask_matrix, int row, int col, int qc_size, int pad_bit);

// print mask matrix
void print_matrix_to_file(int** matrix, int row, int col, char* matrix_file);

#endif
