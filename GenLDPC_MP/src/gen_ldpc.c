#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alloc.h"
#include "bipartite_graph.h"
#include "intio.h"
#include "gen_ldpc.h"
#include "mod2convert.h"
#include "mod2sparse.h"
#include "open.h"
#include "queue.h"
#include "rcode.h"

int *mask_vec;

int gen_ldpc_matrix(mod2sparse *basePCH, mod2sparse *exPCH, int rownum, int colnum, int ex_factor, int pad_bit, int g, WD_vector *row_dt, WD_vector *col_dt, int maximum_girth, int minimum_girth, int** occupied_matrix, int** fade_matrix, mod2sparse *cropPCH0)
{
    int i, j, k, p, index, temp, try_cnt;
    int current_col, current_row;   
    int unoperated_col;
    int row_ones; // the number of 1's in all rows of base matrix
    int col_ones; // the number of 1's in all columns of base matrix
    int tot_rows, tot_cols;
    int count1, count2, ones_inserted, pseudo_weight, temp_counter;
    int tmp_girth, tmp_max_girth;

    int *row_cnt;
    int *col_cnt;
    int offset_pos;

    int max_rw=0;
    int mask_flag;
    int row_sel;

    tanner_graph *ex_g;
    mod2entry *e, *e_pre0, *e_pre1;
    int inserted;
    int shift0, shift1;
    int cycle_deg_constraint;
    int t_size = rownum - g;
    int last_element[17];

    // Matrix Parameter check:
    // 1. Check if the number of "1" in rows is equal to the number of "1" in columns
    // 2. Check if the number of rows(cols) in distribution vector is equal to the matrix

    row_ones = 0;
    col_ones = 0;
    tot_rows = 0;
    tot_cols = 0;
    for (i = 0; i < row_dt->n; i++)
    {
        row_ones += row_dt->wd[i].weight * row_dt->wd[i].num;
        tot_rows += row_dt->wd[i].num;
    }
    for (i = 0; i < col_dt->n; i++)
    {
        col_ones += col_dt->wd[i].weight * col_dt->wd[i].num;
        tot_cols += col_dt->wd[i].num;
    }

    if (row_ones != col_ones)
    {
        fprintf(stderr, "Error: The number of 1's in rows is not equal to the number of 1's in columns!\n");
        exit(1);
    }
    if (tot_rows != rownum)
    {
        fprintf(stderr, "Error: The number of rows in row distribution vector is not equal to the number of rows in base matrix!\n");
        exit(1);
    }
    if (tot_cols != colnum)
    {
        fprintf(stderr, "Error: The number of columns in column distribution vector is not equal to the number of columns in base matrix!\n");
        exit(1);
    }

    int pad_bit_num = pad_bit;
    int drop_len = pad_bit;
    int mask_len = ex_factor - drop_len;

    int max_mask_block_num = row_ones / tot_rows;

    int low_limit_mask = max_mask_block_num / (tot_rows-1);
    int up_limit_mask = low_limit_mask + 1;

    int row_mask_cnt[tot_rows];
    for (int i=0; i<tot_rows; i++)
        row_mask_cnt[i] = 0;

    int shift_delta[17];
    
    shift_delta[0]  =   0;
    shift_delta[1]  =   13;
    shift_delta[2]  =   19;
    shift_delta[3]  =   29;
    shift_delta[4]  =   41;
    shift_delta[5]  =   67;
    shift_delta[6]  =   73;
    shift_delta[7]  =   79;
    shift_delta[8]  =   91;
    shift_delta[9]  =   97;
    shift_delta[10] =   103;
    shift_delta[11] =   111;
    shift_delta[12] =   119;
    shift_delta[13] =   127;
    shift_delta[14] =   131;
    shift_delta[15] =   137;
    shift_delta[16] =   149;

    for (int i=0; i<rownum; i++)
    {
        for (int j=0; j<colnum; j++)
        {
            occupied_matrix[i][j] = 0;
            fade_matrix[i][j] = 0;
        }
    }

    // Initialize tanner graph for BFS in checking cycle and cycle degree
    ex_g = chk_alloc(1, sizeof(tanner_graph));
    ex_g->bit_node = chk_alloc(colnum * ex_factor, sizeof(node_entry));
    ex_g->check_node = chk_alloc(rownum * ex_factor, sizeof(node_entry));
    ex_g->matrix = exPCH;

    // assign node degree for each bit node, to be used in counting cycle degree
    offset_pos = 0;
    for (i = 0; i< col_dt->n; i++)
    {
        for (j = 0; j < col_dt->wd[i].num; j++)
        {
            for (k = 0; k < ex_factor; k++)
                if (col_dt->wd[i].weight != 2)
                    ex_g->bit_node[offset_pos * ex_factor + k].node_degree = col_dt->wd[i].weight;
                else
                    ex_g->bit_node[offset_pos * ex_factor + k].node_degree = 0;
            offset_pos++;
        }
    }

    for (i = 0; i<rownum*ex_factor; i++)
        ex_g->check_node[i].node_degree = 0;

    // Init counters of rows
    row_cnt = chk_alloc(rownum, sizeof(int));
    current_row = 0;

    for (i = 0; i < row_dt->n-1; i++)
    {
        for (j=0; j<row_dt->wd[i].num;)
        {
            // only 1 one in this row
            if (current_row < g - 1)
                row_cnt[current_row] = row_dt->wd[i].weight-4;
            else if (current_row == g - 1)
                row_cnt[current_row] = row_dt->wd[i].weight-3;
            else if (current_row < rownum - 1)
                row_cnt[current_row] = row_dt->wd[i].weight-1;
            else
                row_cnt[current_row] = row_dt->wd[i].weight;

            if (row_cnt[current_row]>max_rw)
                max_rw = row_cnt[current_row];
            
            current_row++;
            j++;
        }
    }

    // Insert G in matrix
    mod2sparse_clear(basePCH);
    for (i = 0; i < g; i++)
    {
        int use_location;
        int g_ind = 1;
        for (int j = colnum - g; j < colnum; j++)
        {
            k = colnum - 1 - j;
            use_location = (k % g) != i;
            if ((i == (g - 1)) && (j == (colnum - 1)))
                use_location = 0; // Make matrix invertible

            if (use_location)
            {
                occupied_matrix[i][j] = 1;
                e = mod2sparse_insert(basePCH, i, j);
                e->shift = base_g[i][j-colnum+g];
                insert_sub_matrix(exPCH, i*ex_factor, j * ex_factor, ex_factor, e->shift);
                g_ind++;
            }
        }
    }

    // insert E in matrix, except for the first column
    for (i=g; i<rownum; i++)
    {
        occupied_matrix[i][colnum-1-i] = 1;
        mod2sparse_insert(basePCH, i, colnum-1-i);
        if (i == rownum - 1)
            insert_sub_matrix_row_mask(exPCH, (rownum-1)*ex_factor, (colnum - rownum)*ex_factor, pad_bit, ex_factor-1, ex_factor, 0);
        else
            insert_sub_matrix(exPCH, i*ex_factor, (colnum-1-i)*ex_factor, ex_factor, 0);
    }

    // Init counters of cols
    col_cnt = chk_alloc(colnum, sizeof(int));
    current_col = 0;

    for (i = 0; i < col_dt->n; i++)
    {
        for (j = 0; j < col_dt->wd[i].num; j++)
        {
            if (current_col >= colnum - g)
                col_cnt[current_col] = 0; // already insert the G matrix
            else if (current_col > colnum - rownum && current_col < colnum - g)
                col_cnt[current_col] = col_dt->wd[i].weight-1;
            else
                col_cnt[current_col] = col_dt->wd[i].weight;
            current_col++;
        }
    }

    // Based on frame matrix, construct the LDPC matrix by inserting "1" to each
    // column. The number of "1" is equal to the weight of the column. The row 
    // position must be choosen from the given queue created previously.
    current_col = colnum - 1; // start form the last column
    ones_inserted = 0;
    tmp_max_girth = maximum_girth;
    temp_counter = 0;
    if (rownum == 5)
        row_sel = rownum + 1;
    else
        row_sel = rownum;

    for (i=0 ; i<colnum; i++)
    {
        pseudo_weight = col_cnt[current_col];
        mask_flag = mask_vec[current_col];

        for(k=0; k<pseudo_weight; k++)
        {
            inserted = 0;

            int is_drop_loc;
            if ((current_col <= colnum - rownum) && (mask_flag == 1))
            {
                is_drop_loc = 1;
                mask_flag = 0;
            }
            else
                is_drop_loc = 0;
            
            // select row
            temp = 0;
            for (index=0; index < row_sel - 1; index++)
            {
                if (row_cnt[index] > current_col) // row_weight > unoperated col
                {
                    current_row = index;
                    temp=1;
                    break;
                }
            }

            if (temp==0)
            {
                current_row = rand() % (row_sel-1);

                if (current_col >= colnum - rownum)
                {
                    while(row_cnt[current_row] == 0 || (current_row >= colnum-1-current_col)
                        || ((is_drop_loc==1) && (row_mask_cnt[current_row]>=up_limit_mask)))
                        current_row = rand() % (row_sel-1);
                }
                else
                {
                    while(row_cnt[current_row] == 0 || row_cnt[current_row] < max_rw)
                        current_row = rand() % (row_sel-1);
                    if ((is_drop_loc==1) && (row_mask_cnt[current_row]>=up_limit_mask))
                    {
                        is_drop_loc = 0;
                        mask_flag = 1;
                    }
                }
            }

            row_cnt[current_row]--;
            if(is_drop_loc==1)
                row_mask_cnt[current_row]++;

            for (tmp_girth = tmp_max_girth; tmp_girth >= minimum_girth; tmp_girth -= 2)
            {
                count2 = 0;
                while (!inserted)
                {
                    count1 = 0;
                    while (mod2sparse_find(basePCH, current_row, current_col)  || ((current_row < g)&&(current_col > (colnum-1-g))))
                    {
                        // if it has '1' in current position or in T matrix region
                        // re-find row
                        row_cnt[current_row]++;
                        if(is_drop_loc==1)
                            row_mask_cnt[current_row]--;

                        // select row
                        temp = 0;
                        for (index=0; index<row_sel-1; index++)
                        {
                            if (row_cnt[index] > current_col) // row_weight > unoperated col
                            {
                                temp=1;
                                current_row = index;
                                break;
                            }
                        }

                        if (temp==0)
                        {
                            current_row = rand() % (row_sel-1);

                            if (current_col >= colnum - rownum)
                            {
                                while(row_cnt[current_row] == 0 || (current_row >= colnum-1-current_col)
                                    || ((is_drop_loc==1) && (row_mask_cnt[current_row]>=up_limit_mask)))
                                    current_row = rand() % (row_sel-1);
                            }
                            else
                            {
                                while(row_cnt[current_row] == 0 || row_cnt[current_row] < max_rw)
                                    current_row = rand() % (row_sel-1);
                                if ((is_drop_loc==1) && (row_mask_cnt[current_row]>=up_limit_mask))
                                {
                                    is_drop_loc = 0;
                                    mask_flag = 1;
                                }
                            }
                        }

                        row_cnt[current_row]--;
                        if(is_drop_loc==1)
                            row_mask_cnt[current_row]++;

                        // repeat
                        count1++;
                        if (count1 > row_ones - t_size - ones_inserted - 18) // after trying all possible available row positions
                        {
                            free(row_cnt);
                            free(ex_g->bit_node);
                            free(ex_g->check_node);
                            free(ex_g);
                            return -1;
                        }
                    }

                    if (is_drop_loc==1 && current_col == colnum - rownum)
                    {
                        mod2sparse_insert(basePCH, current_row, current_col);
                        fade_matrix[current_row][current_col] = 1;
                    }
                    else if (is_drop_loc==1 && current_col != colnum - rownum)
                    {
                        mod2sparse_insert(basePCH, current_row, current_col);
                        mod2sparse_insert(basePCH, rownum-1, current_col);
                        fade_matrix[current_row][current_col] = 1;
                        occupied_matrix[rownum-1][current_col] = 1;
                    }
                    else
                    {
                        mod2sparse_insert(basePCH, current_row, current_col);
                        occupied_matrix[current_row][current_col] = 1;
                    }
                    count2++;

                    try_cnt = 0;
                    cycle_deg_constraint = 1;
                    // select shift
                    while ((try_cnt < ex_factor) && (cycle_deg_constraint == 1))
                    {
                        if (is_drop_loc==1 && current_col == colnum - rownum)
                        {
                            shift0 = rand() % ex_factor;
                            shift1 = 0;
                            insert_sub_matrix_col_mask(exPCH, current_row*ex_factor, current_col*ex_factor,
                                         0, drop_len-1, ex_factor, shift0);
                        }
                        else if (current_col == colnum -g - 1)
                        {
                            shift0 = 0;
                            insert_sub_matrix(exPCH, current_row*ex_factor, current_col*ex_factor,
                                         ex_factor, shift0);
                        }
                        else if (is_drop_loc==1 && current_col < colnum - rownum)
                        {
                            shift0 = rand() % ex_factor;
                            shift1 = rand() % ex_factor;
                            insert_sub_matrix_col_mask(exPCH, current_row*ex_factor, current_col*ex_factor,
                                         shift1%ex_factor, (pad_bit-1+shift1)%ex_factor, ex_factor, shift0);
                            insert_sub_matrix_row_mask(exPCH, (rownum-1)*ex_factor, current_col*ex_factor,
                                         pad_bit, ex_factor-1, ex_factor, shift1);
                        }
                        else
                        {
                            shift0 = rand() % ex_factor;
                            insert_sub_matrix(exPCH, current_row * ex_factor, current_col * ex_factor,
                                         ex_factor, shift0);
                        }

                        count2++;

                        cycle_deg_constraint = 0;
                        if(is_drop_loc==1)
                        {
                            for (p = 0; p<ex_factor; p++)
                            {
                                if((shift1 - shift0)%ex_factor <= (shift1 - shift0 + pad_bit - 1) % ex_factor)
                                {
                                    if ((p >= (shift1 - shift0)%ex_factor) && (p <= (shift1 - shift0 + pad_bit - 1) % ex_factor))
                                        continue;
                                }
                                else
                                {
                                    if ((0 <= p && p <= (shift1 - shift0 + pad_bit - 1) % ex_factor) ||
                                        ((shift1 - shift0)%ex_factor <= p && p <= ex_factor - 1))
                                        continue;
                                }
                                if(check_cycle_deg(ex_g, current_row*ex_factor+p, tmp_girth))
                                {
                                    cycle_deg_constraint = 1;
                                    break;
                                }
                            }
                            for (p = 0; p < pad_bit; p++)
                            {
                                if(check_cycle_deg(ex_g, (rownum-1)*ex_factor+p, tmp_girth))
                                {
                                    cycle_deg_constraint = 1;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            for (p = 0; p<ex_factor; p++)
                            {
                                if(check_cycle_deg(ex_g, current_row*ex_factor+p, tmp_girth))
                                {
                                    cycle_deg_constraint = 1;
                                    break;
                                }
                            }
                        }

                        if (cycle_deg_constraint)
                        {
                            if (is_drop_loc ==1 && current_col == colnum - rownum)
                            {
                                remove_sub_matrix_col_mask(exPCH, current_row*ex_factor, current_col*ex_factor, 0, drop_len-1, ex_factor, shift0);
                            }
                            else if(is_drop_loc==1 && current_col != colnum - rownum)
                            {
                                remove_sub_matrix_col_mask(exPCH, current_row*ex_factor, current_col*ex_factor, shift1%ex_factor, (pad_bit-1+shift1)%ex_factor, ex_factor, shift0);
                                remove_sub_matrix_row_mask(exPCH, (rownum-1)*ex_factor, current_col*ex_factor, pad_bit, ex_factor-1, ex_factor, shift1);
                            }
                            else
                            {
                                remove_sub_matrix(exPCH, current_row*ex_factor, current_col*ex_factor, ex_factor, shift0);
                            }
                        }
                        try_cnt++;
                    }
                    
                    if (cycle_deg_constraint)
                    {
                        if (is_drop_loc==1 && current_col == colnum - rownum)
                        {
                            fade_matrix[current_row][current_col] = 0;
                            mod2sparse_delete(basePCH, mod2sparse_find(basePCH, current_row, current_col));
                        }
                        else if(is_drop_loc==1 && current_col != colnum - rownum)
                        {
                            fade_matrix[current_row][current_col] = 0;
                            occupied_matrix[rownum-1][current_col] = 0;
                            mod2sparse_delete(basePCH, mod2sparse_find(basePCH, current_row, current_col));
                            mod2sparse_delete(basePCH, mod2sparse_find(basePCH, rownum-1, current_col));
                        }
                        else
                        {
                            occupied_matrix[current_row][current_col] = 0;
                            mod2sparse_delete(basePCH, mod2sparse_find(basePCH, current_row, current_col));
                        }

                        // change a row
                        row_cnt[current_row]++;
                        if (is_drop_loc==1)
                            row_mask_cnt[current_row]--;

                        // select row
                        temp = 0;
                        for (index=0; index<row_sel-1; index++)
                        {
                            if (row_cnt[index] > current_col) // row_weight > unoperated col
                            {
                                temp=1;
                                current_row = index;
                                break;
                            }
                        }

                        if (temp==0)
                        {
                            current_row = rand() % (row_sel-1);

                            if (current_col >= colnum - rownum)
                            {
                                while(row_cnt[current_row] == 0 || (current_row >= colnum-1-current_col)
                                    || ((is_drop_loc==1) && (row_mask_cnt[current_row]>=up_limit_mask)))
                                    current_row = rand() % (row_sel-1);
                            }
                            else
                            {
                                while(row_cnt[current_row] == 0 || row_cnt[current_row] < max_rw)
                                    current_row = rand() % (row_sel-1);
                                if ((is_drop_loc==1) && (row_mask_cnt[current_row]>=up_limit_mask))
                                {
                                    is_drop_loc = 0;
                                    mask_flag = 1;
                                }
                            }
                        }

                        row_cnt[current_row]--;
                        if (is_drop_loc==1)
                            row_mask_cnt[current_row]++;

                        if (count2 > row_ones - t_size - ones_inserted - 18) // after trying all possible available row positions
                        {
                            tmp_max_girth -= 2;
                            break;
                        }
                    }
                    else
                    {
                        inserted = 1;
                        if (is_drop_loc == 1 && current_col == colnum - rownum)
                        {
                            e = mod2sparse_find(basePCH, current_row, current_col);
                            e->shift = shift0; 
                        }
                        else if (current_col == colnum -g - 1)
                        {
                            e = mod2sparse_find(basePCH, current_row, current_col);
                            e->shift = shift0;
                        }
                        else if (is_drop_loc == 1 && current_col != colnum - rownum)
                        {
                            e = mod2sparse_find(basePCH, current_row, current_col);
                            e->shift = shift0; 
                            e = mod2sparse_find(basePCH, rownum-1, current_col);
                            e->shift = shift1;
                        }
                        else
                        {
                            e = mod2sparse_find(basePCH, current_row, current_col);
                            e->shift = shift0;
                        }
                    }
                }
            }

            if (!inserted)
            {
                free(row_cnt);
                free(ex_g->bit_node);
                free(ex_g->check_node);
                free(ex_g);
                return -2;
            }
            ones_inserted++;

            max_rw = 0;
            for (index=0; index<row_sel-1; index++)
            {
                if (mod2sparse_find(basePCH, index, current_col))
                    continue;
                if (row_cnt[index] > max_rw)
                    max_rw = row_cnt[index];
            }
        }

        max_rw = 0;
        for (index=0; index<rownum; index++)
        {
            if (row_cnt[index] > max_rw)
                max_rw = row_cnt[index];
        }
        /* first, insert G;
           second, insert B and E which don't include drop circulant;
           third, insert B and E which include drop circulant;
           last, insert A and C. 
        */

        current_col--;

        /*
        FILE *fp;
        fp = fopen("./output/LDPC_8x73x512_w6_dense6_QC_H_pad256.txt", "w");
        for (k=0; k < mod2sparse_rows(basePCH); k++)
        {
            int pre_col = -1;

            for (e = mod2sparse_first_in_row(basePCH, k);
                    !mod2sparse_at_end(e);
                    e = mod2sparse_next_in_row(e))
            {
                for (j = (pre_col+1); j<e->col; j++)
                    fprintf(fp, " -1");
                fprintf(fp, "%4d", e->shift);
                pre_col = e->col;
            }

            for (j=pre_col+1; j<mod2sparse_cols(basePCH); j++)
            fprintf(fp, " -1");
        fprintf(fp, "\n");
        }
        fclose(fp);
        */
    }

    /*
        if ((current_col==colnum-t_size) && (max_mask_block_num>0))
            current_col = colnum - rownum - 1;
        else if ((current_col==colnum-rownum-max_mask_block_num) && (max_mask_block_num>0))
            current_col = colnum - t_size - 1;
        else if ((current_col==colnum-rownum) && (max_mask_block_num>0))
            current_col = colnum - rownum - max_mask_block_num - 1;
        else
            current_col--;
    }
    */

    int cyc4_num_ibex, cyc6_num_ibex;
    cyc4_num_ibex = cyc4_cal_ibex(basePCH, fade_matrix, ex_factor, pad_bit);
    cyc6_num_ibex = cyc6_cal_ibex(basePCH, fade_matrix, ex_factor, pad_bit);

    int num_cyc4, num_cyc6;
    int num_cyc6_tmp;
    int shift_initial, shift_tmp;
    int free_cycle6_round = 6; // 消环次数
    for (int ii=0; ii<free_cycle6_round ; ii++)
    {
        for (int i=colnum-1; i>=0 ; i--)
        {
            for (e = mod2sparse_first_in_col(basePCH, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
            {
                if(e->col >= colnum-g-1)
                    continue;
                if ((e->col>colnum-rownum-1) && (e->row==colnum-e->col-1))
                    continue;
                if((fade_matrix[e->row][e->col] != 0) || (e->row==rownum-1))
                    continue;

                shift_initial = e->shift;
                shift_tmp = shift_initial;
                remove_sub_matrix(exPCH, e->row*ex_factor, e->col*ex_factor, ex_factor, shift_initial);
                num_cyc6_tmp = search_cyc6(basePCH, e, fade_matrix, rownum, ex_factor, pad_bit);

                for (k=1; k<ex_factor; k++)
                {
                    e->shift = (shift_initial + k) % ex_factor;
                    num_cyc4 = search_cyc4(basePCH, e, fade_matrix, rownum, ex_factor, pad_bit);
                    if (num_cyc4>0)
                        continue;
                    num_cyc6 = search_cyc6(basePCH, e, fade_matrix, rownum, ex_factor, pad_bit);
                    if (num_cyc6 < num_cyc6_tmp)
                    {
                        shift_tmp = e->shift;
                        num_cyc6_tmp = num_cyc6;
                    }
                }
                e->shift = shift_tmp;
                insert_sub_matrix(exPCH, e->row*ex_factor, e->col*ex_factor, ex_factor, shift_tmp);
            }
        }
    }

    // clipping matrix
    mod2sparse_copy(basePCH, cropPCH0);

    mod2entry *e_crop0;
    for (int i=0; i<rownum; i++)
    {
        for (int j=colnum-1; j>=0; j--)
        {
            if (!mod2sparse_find(basePCH, i, j))
                continue;
            e        = mod2sparse_find(basePCH, i, j);
            e_crop0 = mod2sparse_find(cropPCH0, i, j);
            if ((j >= colnum-rownum-6) && (j < colnum-rownum))
                mod2sparse_delete(cropPCH0, e_crop0);
            else
                e_crop0->shift = e->shift;
        }
    }
    
    free(row_cnt);
    free(col_cnt);
    free(ex_g->bit_node);
    free(ex_g->check_node);
    free(ex_g);

    return tmp_max_girth;
}

int gen4_ldpc_matrix(mod2sparse *basePCH, mod2sparse *exPCH, int rownum, int colnum, int ex_factor, int t_size, WD_vector *row_dt, WD_vector *col_dt, int maximum_girth, int minimum_girth)
{
    int i, j, k, p, index, temp, try_cnt;
    int current_col, current_row;   
    int row_ones; // the number of 1's in all rows of base matrix
    int col_ones; // the number of 1's in all columns of base matrix
    int tot_rows, tot_cols;
    int count1, count2, ones_inserted, pseudo_weight, temp_counter;
    int tmp_girth, tmp_max_girth;

    int *row_cnt;
    int offset_pos;

    tanner_graph *ex_g;
    mod2entry *e;
    int inserted;
    int shift;
    int cycle_deg_constraint;

    // Matrix Parameter check:
    // 1. Check if the number of "1" in rows is equal to the number of "1" in columns
    // 2. Check if the number of rows(cols) in distribution vector is equal to the matrix

    row_ones = 0;
    col_ones = 0;
    tot_rows = 0;
    tot_cols = 0;
    for (i = 0; i < row_dt->n; i++)
    {
        row_ones += row_dt->wd[i].weight * row_dt->wd[i].num;
        tot_rows += row_dt->wd[i].num;
    }
    for (i = 0; i < col_dt->n; i++)
    {
        col_ones += col_dt->wd[i].weight * col_dt->wd[i].num;
        tot_cols += col_dt->wd[i].num;
    }

    if (row_ones != col_ones)
    {
        fprintf(stderr, "Error: The number of 1's in rows is not equal to the number of 1's in columns!\n");
        exit(1);
    }
    if (tot_rows != rownum)
    {
        fprintf(stderr, "Error: The number of rows in row distribution vector is not equal to the number of rows in base matrix!\n");
        exit(1);
    }
    if (tot_cols != colnum)
    {
        fprintf(stderr, "Error: The number of columns in column distribution vector is not equal to the number of columns in base matrix!\n");
        exit(1);
    }

    // Initialize tanner graph for BFS in checking cycle and cycle degree
    ex_g = chk_alloc(1, sizeof(tanner_graph));
    ex_g->bit_node = chk_alloc(colnum * ex_factor, sizeof(node_entry));
    ex_g->check_node = chk_alloc(rownum * ex_factor, sizeof(node_entry));
    ex_g->matrix = exPCH;

    for (i=0; i < floor(col_dt->n/2); i++) {
        temp = col_dt->wd[i].weight;
        col_dt->wd[i].weight = col_dt->wd[col_dt->n - 1 - i].weight;
        col_dt->wd[col_dt->n - 1 - i].weight = temp;
        temp = col_dt->wd[i].num;
        col_dt->wd[i].num = col_dt->wd[col_dt->n - 1 - i].num;
        col_dt->wd[col_dt->n - 1 - i].num = temp;
    }

    // assign node degree for each bit node, to be used in counting cycle degree
    offset_pos = 0;
    for (i = 0; i< col_dt->n; i++) {
        for (j = 0; j < col_dt->wd[i].num; j++) {
            for (k = 0; k < ex_factor; k++)
                if (col_dt->wd[i].weight != 2)
                    ex_g->bit_node[offset_pos * ex_factor + k].node_degree = col_dt->wd[i].weight;
                else
                    ex_g->bit_node[offset_pos * ex_factor + k].node_degree = 0;
            offset_pos++;
        }
    }

    for (i = 0; i < rownum * ex_factor; i++)
        ex_g->check_node[i].node_degree = 0;

    // Init counters of rows
    row_cnt = chk_alloc(rownum, sizeof(int));
    current_row = 0;

    for (i = 0; i < row_dt->n; i++) {
        for (j=0; j<row_dt->wd[i].num; j++) {
            if (current_row < t_size)
                row_cnt[current_row] = row_dt->wd[i].weight - 1;
            else
                row_cnt[current_row] = row_dt->wd[i].weight;

            current_row++;
        }
    }

    // Insert T in matrix
    mod2sparse_clear(basePCH);
    for (i = 0; i < t_size; i++)
    {
        mod2sparse_insert(basePCH, i, colnum - t_size + i);
        shift = 0;
        insert_sub_matrix(exPCH, i*ex_factor, (colnum-t_size+i)*ex_factor, ex_factor, shift);
    }

    // Based on frame matrix, construct the LDPC matrix by inserting "1" to each
    // column. The number of "1" is equal to the weight of the column. The row 
    // position must be choosen from the given queue created previously.
    current_col = colnum - 1; // start form the last column
    ones_inserted = 0;
    tmp_max_girth = maximum_girth;
    temp_counter = 0;

    for (i= col_dt->n - 1 ; i>=0; i--) {
        for(j = 0; j < col_dt->wd[i].num; j++) {
            if (temp_counter < t_size)
                pseudo_weight = col_dt->wd[i].weight - 1; // already insert the T matrix
            else
                pseudo_weight = col_dt->wd[i].weight;

            for (k=0; k<pseudo_weight; k++) {

                inserted = 0;

                // select row
                temp = 0;
                for (index = 0; index < rownum; index++) {
                    if (row_cnt[index] > current_col) {
                        current_row = index;
                        temp=1;
                        break;
                    }
                }

                if (temp==0)
                {
                    current_row = rand() % rownum;
                    while(row_cnt[current_row] == 0)
                        current_row = rand() % rownum;
                }

                row_cnt[current_row]--;

                for (tmp_girth = tmp_max_girth; tmp_girth >= minimum_girth; tmp_girth -= 2) {
                    count2 = 0;
                    while (!inserted) {
                        count1 = 0;
                        while (mod2sparse_find(basePCH, current_row, current_col)  || ((current_row < t_size)&&(current_col > (colnum-1-t_size))))
                        {
                            // if it has '1' in current position or in T matrix region
                            // re-find row
                            row_cnt[current_row]++;

                            // select row
                            temp = 0;
                            for (index=0; index<rownum; index++) {
                                if (row_cnt[index] > current_col) {
                                    temp=1;
                                    current_row = index;
                                    break;
                                }
                            }

                            if (temp==0) {
                                current_row = rand() % rownum;

                                while(row_cnt[current_row] == 0)
                                    current_row = rand() % rownum;
                            }

                            row_cnt[current_row]--;

                            count1++;
                            if (count1 > row_ones - t_size - ones_inserted)
                            {
                                free(row_cnt);
                                free(ex_g->bit_node);
                                free(ex_g->check_node);
                                free(ex_g);
                                return -1;
                            }
                        }

                        mod2sparse_insert(basePCH, current_row, current_col);
                        count2++;

                        try_cnt = (current_row == 0) ? ex_factor-1 : 0;
                        cycle_deg_constraint = 1;

                        // select shift
                        while ((try_cnt < ex_factor) && (cycle_deg_constraint == 1)) {
                            if (current_row == 0)
                                shift = 0;
                            else
                                shift = rand() % ex_factor;
                                
                            insert_sub_matrix(exPCH, current_row*ex_factor, current_col*ex_factor, ex_factor, shift);
                            count2++;
        
                            cycle_deg_constraint = 0;

                            for (p = 0; p<ex_factor; p++)
                            {
                                if(check_cycle_deg(ex_g, current_row*ex_factor+p, tmp_girth))
                                {
                                    cycle_deg_constraint = 1;
                                    break;
                                }
                            }
             
                            if (cycle_deg_constraint)
                            {
                                remove_sub_matrix(exPCH, current_row*ex_factor, current_col*ex_factor, ex_factor, shift);
                            }
                            try_cnt++;
                        }
                        
                        if (cycle_deg_constraint)
                        {
                            mod2sparse_delete(basePCH, mod2sparse_find(basePCH, current_row, current_col));

                            // change a row
                            row_cnt[current_row]++;

                            // select row
                            temp = 0;
                            for (index = 0; index < rownum; index++)
                            {
                                if (row_cnt[index] > current_col) // row_weight > unoperated col
                                {
                                    temp = 1;
                                    current_row = index;
                                    break;
                                }
                            }

                            if (temp == 0)
                            {
                                current_row = rand() % rownum;

                                while(row_cnt[current_row] == 0)
                                    current_row = rand() % rownum;
                            }
                            row_cnt[current_row]--;

                            if (count2 > row_ones - t_size - ones_inserted)
                            {
                                tmp_max_girth -= 2;
                                break;
                            }
                        } else {
                            inserted = 1;
                            e = mod2sparse_find(basePCH, current_row, current_col);
                            e->shift = shift;
                        }
                    }
                }

                if (!inserted) {
                    free(row_cnt);
                    free(ex_g->bit_node);
                    free(ex_g->check_node);
                    free(ex_g);
                    return -2;
                }
                ones_inserted++;
            }

            current_col--;
            temp_counter++;
        }
    }

    free(row_cnt);
    free(ex_g->bit_node);
    free(ex_g->check_node);
    free(ex_g);

    return tmp_max_girth;
}


void insert_sub_matrix(mod2sparse *mainMatrix, int startRow, int startCol, int subDim, int shift)
{
    int i, rowPos, colPos;
    for (i=0; i<subDim; i++)
    {
        rowPos = startRow + i;
        colPos = startCol + (i + shift) % subDim;
        mod2sparse_insert(mainMatrix, rowPos, colPos);
    }
}

void insert_sub_matrix_row_mask(mod2sparse *mainMatrix, int startRow, int startCol, int mask_row_start, int mask_row_end, int subDim, int shift)
{
    int i, rowPos, colPos;
    for (i=0; i<subDim; i++)
    {
        if ((i>=mask_row_start) && (i<=mask_row_end))
            continue;

        rowPos = startRow + i;
        colPos = startCol + (i + shift) % subDim;
        mod2sparse_insert(mainMatrix, rowPos, colPos);
    }
}

void insert_sub_matrix_col_mask(mod2sparse *mainMatrix, int startRow, int startCol, int mask_col_start, int mask_col_end, int subDim, int shift)
{
    int i, rowPos, colPos;
    for (i=0; i<subDim; i++)
    {
        if (mask_col_start <= mask_col_end)
        {
            if ((i>=mask_col_start) && (i<=mask_col_end))
                continue;
        }
        else
        {
            if ((i>=0 && i<=mask_col_end) || (i>=mask_col_start && i<=subDim - 1))   
                continue;
        }

        colPos = startCol + i;
        rowPos = startRow + (i - shift + subDim) % subDim;
        mod2sparse_insert(mainMatrix, rowPos, colPos);
    }
}

void remove_sub_matrix(mod2sparse *mainMatrix, int startRow, int startCol, int subDim, int shift)
{
    int i, rowPos, colPos;
    for (i=0; i<subDim; i++)
    {
        rowPos = startRow + i;
        colPos = startCol + (i + shift) % subDim;
        mod2sparse_delete(mainMatrix, mod2sparse_find(mainMatrix, rowPos, colPos));
    }
}

void remove_sub_matrix_row_mask(mod2sparse *mainMatrix, int startRow, int startCol, int mask_row_start, int mask_row_end, int subDim, int shift)
{
    int i, rowPos, colPos;
    for (i=0; i<subDim; i++)
    {
        if ((i>=mask_row_start) && (i<=mask_row_end))
            continue;

        rowPos = startRow + i;
        colPos = startCol + (i + shift) % subDim;
        mod2sparse_delete(mainMatrix, mod2sparse_find(mainMatrix, rowPos, colPos));
    }
}

void remove_sub_matrix_col_mask(mod2sparse *mainMatrix, int startRow, int startCol, int mask_col_start, int mask_col_end, int subDim, int shift)
{
    int i, rowPos, colPos;
    for (i=0; i<subDim; i++)
    {
        if (mask_col_start <= mask_col_end)
        {
            if ((i>=mask_col_start) && (i<=mask_col_end))
                continue;
        }
        else
        {
            if ((i>=0 && i<=mask_col_end) || (i>=mask_col_start && i<=subDim - 1))   
                continue;
        }
    
        colPos = startCol + i; 
        rowPos = startRow + (i - shift + subDim) % subDim;
        mod2sparse_delete(mainMatrix, mod2sparse_find(mainMatrix, rowPos, colPos));
    }
}

void gen_ldpc_files(char *codefile, char *occupiedfile, char *fadefile, char *cyclefile, int filenum, WD_vector *row_dt, WD_vector *col_dt, int ex_factor, int pad_bit, 
                    int g, int start_girth, int end_girth, float expected_avg_girth)
{
    char file[256];
    int b_rows, b_cols;
    mod2sparse *baseMatrix, *exMatrix, *exMatrix_mask;
    mod2sparse *cropMatrix0;

    int i, j, file_counter, t_size;
    int gen, inv;
    FILE *temp_f;
    char postfix[5];

    tanner_graph *ex_g;
    float cur_avg_girth;
    mod2entry *e;
    int pre_col;
    int index;

    int** occupied_matrix;
    int** fade_matrix;
    int** H_matrix;
    int** H_matrix_tmp;
    int drop_len = pad_bit;
    int mask_len = ex_factor - pad_bit;
    int num_drop_block;

    int* row_weight_cnt;
    int* col_weight_cnt;
    int* row_drop_cnt;

    b_rows = 0;
    b_cols = 0;
    for (i = 0; i < row_dt->n; i++)
        b_rows += row_dt->wd[i].num;
    for (i = 0; i < col_dt->n; i++)
        b_cols += col_dt->wd[i].num;

    t_size = b_rows - g;

    baseMatrix = mod2sparse_allocate(b_rows, b_cols);
    cropMatrix0 = mod2sparse_allocate(b_rows, b_cols);
    exMatrix = mod2sparse_allocate(b_rows * ex_factor, b_cols * ex_factor);
    
    occupied_matrix = (int**)chk_alloc(b_rows, sizeof(*occupied_matrix));
    for (i=0; i<b_rows; i++)
        occupied_matrix[i] = (int*)chk_alloc(b_cols, sizeof(*occupied_matrix[i]));
    fade_matrix = (int**)chk_alloc(b_rows, sizeof(*fade_matrix));
    for (i=0; i<b_rows; i++)
        fade_matrix[i] = (int*)chk_alloc(b_cols, sizeof(*fade_matrix[i]));
    H_matrix = (int**)chk_alloc(b_rows, sizeof(*H_matrix));
    for (i=0; i<b_rows; i++)
        H_matrix[i] = (int*)chk_alloc(b_cols, sizeof(*H_matrix[i]));
    H_matrix_tmp = (int**)chk_alloc(b_rows, sizeof(*H_matrix_tmp));
    for (i=0; i<b_rows; i++)
        H_matrix_tmp[i] = (int*)chk_alloc(b_cols, sizeof(*H_matrix_tmp[i]));
    row_weight_cnt = (int*)chk_alloc(b_rows, sizeof(*row_weight_cnt));
    col_weight_cnt = (int*)chk_alloc(b_cols, sizeof(*col_weight_cnt));
    row_drop_cnt = (int*)chk_alloc(b_rows, sizeof(*row_drop_cnt));

    mask_vec = (int*)chk_alloc(b_cols-b_rows+1, sizeof(*mask_vec));
    for (i = 0; i<b_cols-b_rows; i++)
        mask_vec[i] = 0;
    mask_vec[b_cols-b_rows] = 1;
    i = 0;
    if (b_rows == 5)
        mask_vec[b_cols-b_rows] = 0;
    else
    {
        while (i < (int)(((g-1)*b_cols)/b_rows-1))
        {
            j = rand() % (b_cols-b_rows);
            if (mask_vec[j] == 0)
            {
                mask_vec[j] = 1;
                i++;
            }
        }
    }

    file_counter = 0;
    index = 0;
    printf("%d LDPC codes to be generated \n", filenum);
    while (file_counter < filenum)
    {
        if (index == 0)
            printf("Matrix %d construction ... \n", file_counter);

        printf("# Round %d ... \n", index+1);

        gen = gen_ldpc_matrix(baseMatrix, exMatrix, b_rows, b_cols,
                         ex_factor, drop_len, g, row_dt, col_dt, start_girth, end_girth, occupied_matrix, fade_matrix, cropMatrix0);
        
        if (gen >0)
            printf("Matrix is ready\n");
        else if (gen==-1)
            printf("All possible positions are tried, cannot find a position to insert a 1\n");
        else if (gen==-2)
            printf("All possible shifts are tried, cannot find a shift to insert a 1\n");
        
        if (gen>0)
        {
            printf("H Matrix is generateing ...\n");

            for (int i=0; i<b_rows; i++)
            {
                for (int j=0; j<b_cols; j++)
                {
                    H_matrix[i][j] = -1;
                    H_matrix_tmp[i][j]  =-1;
                }
            }

            for (int i=0; i<b_rows; i++)
            {
                for (e = mod2sparse_first_in_row(baseMatrix, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
                {
                    H_matrix[i][e->col] = e->shift;
                    if (occupied_matrix[i][e->col]==1)
                        H_matrix_tmp[i][e->col] = e->shift;
                }
            }
        }

        if (gen>0)
        {
            printf("PHI invertible checking ...\n");

            inv = check_PHI_inv_ibex(exMatrix, b_rows, b_cols, ex_factor,g*ex_factor);

            if (!inv)
                printf("PHI is not invertible!\n");
        }

        if ((gen>0) && (inv==1))
        {
            printf("| Girth checking ... \n");
            ex_g = chk_alloc(1, sizeof(tanner_graph));
            ex_g->bit_node = chk_alloc(mod2sparse_cols(exMatrix), sizeof(node_entry));
            ex_g->check_node = chk_alloc(mod2sparse_rows(exMatrix), sizeof(node_entry));
            ex_g->matrix = exMatrix;

            // Only save the expanded matrices with high average girth
            cur_avg_girth = avg_girth(ex_g);

            // col rand
            if (cur_avg_girth >= expected_avg_girth)
            {
                printf("# Matrix average girth %f \n", cur_avg_girth);
                
                // write matrix to a file
                printf("# Writing matrix to file ... \n");
                print_matrix_to_file(H_matrix, b_rows, b_cols, codefile);
                print_matrix_to_file(occupied_matrix, b_rows, b_cols, occupiedfile);
                print_matrix_to_file(fade_matrix, b_rows, b_cols, fadefile);

                printf("# Writing cycle number to file ...\n");
                FILE *fp_cyc;
                fp_cyc = fopen(cyclefile, "w");
                for (int i=0; i<b_rows; i++)
                {
                    row_weight_cnt[i] = 0;
                    row_drop_cnt[i] = 0;
                }
                for (int j=0; j<b_cols; j++)
                    col_weight_cnt[j] = 0;

                for (int i=0; i<b_rows; i++)
                    for (int j=0; j<b_cols; j++)
                    {
                        if (H_matrix[i][j] != -1)
                        {
                            row_weight_cnt[i]++;
                            col_weight_cnt[j]++;
                        }

                        if (fade_matrix[i][j] == 1)
                            row_drop_cnt[i]++;
                    }

                for (int i=0; i<b_rows; i++)
                    fprintf(fp_cyc, "rows%2d weight:%2d\n", i, row_weight_cnt[i]);
                fprintf(fp_cyc, "----------------------------------------------\n");
                for (int j=0; j<b_cols; j++)
                    fprintf(fp_cyc, "cols%3d weight:%2d\n", j, col_weight_cnt[j]);
                fprintf(fp_cyc, "----------------------------------------------\n");
                for (int i=0; i<b_rows; i++)
                    fprintf(fp_cyc, "rows%2d drop num = %3d\n", i, row_drop_cnt[i]);
                fprintf(fp_cyc, "----------------------------------------------\n");

                int cyc4_num = 0;
                int cyc6_num = 0;
                cyc4_num = cyc4_cal_ibex(baseMatrix, fade_matrix, ex_factor, pad_bit);
                cyc6_num = cyc6_cal_ibex(baseMatrix, fade_matrix, ex_factor, pad_bit);

                fprintf(fp_cyc, "row %2d, col %3d, cyc4 = %d \n", b_rows, b_cols, cyc4_num);
                fprintf(fp_cyc, "row %2d, col %3d, cyc6 = %d \n", b_rows, b_cols, cyc6_num);
                fprintf(fp_cyc, "----------------------------------------------\n");
                cyc4_num = cyc4_cal_ibex(cropMatrix0, fade_matrix, ex_factor, pad_bit);
                cyc6_num = cyc6_cal_ibex(cropMatrix0, fade_matrix, ex_factor, pad_bit);
                fprintf(fp_cyc, "row %2d, col %3d, cyc4 = %d \n", b_rows, b_cols-6, cyc4_num);
                fprintf(fp_cyc, "row %2d, col %3d, cyc6 = %d \n", b_rows, b_cols-6, cyc6_num);
                fprintf(fp_cyc, "----------------------------------------------\n");

                fclose(fp_cyc);

                file_counter++;
                index = 0;
            }
            else
            {
                index++;
            }

            free(ex_g->bit_node);
            free(ex_g->check_node);
            free(ex_g);
        }
        else
        {
            index++;
        }

        mod2sparse_free(baseMatrix);
        mod2sparse_free(exMatrix);

        baseMatrix = mod2sparse_allocate(b_rows, b_cols);
        exMatrix = mod2sparse_allocate(b_rows * ex_factor, b_cols * ex_factor);

        mod2sparse_free(cropMatrix0);
        cropMatrix0 = mod2sparse_allocate(b_rows, b_cols);
    }

    mod2sparse_free(baseMatrix);
    mod2sparse_free(exMatrix);
    mod2sparse_free(cropMatrix0);
    for (i=0; i<b_rows; i++)
        free(occupied_matrix[i]);
    free(occupied_matrix);
    for (i=0; i<b_rows; i++)
        free(fade_matrix[i]);
    free(fade_matrix);
    for (i=0; i<b_rows; i++)
        free(H_matrix[i]);
    free(H_matrix);
    for (i=0; i<b_rows; i++)
        free(H_matrix_tmp[i]);
    free(H_matrix_tmp);
    free(mask_vec);
}

void gen4_ldpc_files(char *codefile, int filenum, WD_vector *row_dt, WD_vector *col_dt, int ex_factor, int g, int start_girth, int end_girth, float expected_avg_girth)
{
    char file[60];
    int b_rows, b_cols;
    mod2sparse *baseMatrix, *exMatrix;

    int i, j, file_counter, t_size;
    int gen, inv;
    FILE *temp_f;
    char postfix[5];

    tanner_graph *ex_g;
    float cur_avg_girth;
    mod2entry *e;
    int pre_col;
    int index;

    b_rows = 0;
    b_cols = 0;
    for (i = 0; i < row_dt->n; i++)
        b_rows += row_dt->wd[i].num;
    for (i = 0; i < col_dt->n; i++)
        b_cols += col_dt->wd[i].num;

    t_size = b_rows - g;

    baseMatrix = mod2sparse_allocate(b_rows, b_cols);
    exMatrix = mod2sparse_allocate(b_rows * ex_factor, b_cols * ex_factor);

    file_counter = 0;
    index = 0;
    printf("%d LDPC codes to be generated \n", filenum);
    while (file_counter < filenum)
    {
        if (index == 0)
            printf("Matrix %d construction ... \n", file_counter);

        printf("# Round %d ... \n", index+1);

        gen = gen4_ldpc_matrix(baseMatrix, exMatrix, b_rows, b_cols,
                         ex_factor, t_size, row_dt, col_dt, start_girth, end_girth);
        
        if (gen > 0)
            printf("Matrix is ready\n");
        else if (gen == -1)
            printf("Construction failed due to placement! \n");
        else if (gen == -2)
            printf("Construction failed due to minimum girth! \n");

        if (gen > 0)
        {
            printf("PHI invertible checking ...\n");

            inv = check_PHI_inv(exMatrix, b_rows, b_cols, ex_factor,g * ex_factor);

            if (!inv)
                printf("PHI in not invertible\n");
        }

        if ((gen > 0) && (inv == 1))
        {
            printf("Girth checking ... \n");
            ex_g = chk_alloc(1, sizeof(tanner_graph));
            ex_g->bit_node = chk_alloc(mod2sparse_cols(exMatrix), sizeof(node_entry));
            ex_g->check_node = chk_alloc(mod2sparse_rows(exMatrix), sizeof(node_entry));
            ex_g->matrix = exMatrix;

            // Only save the expanded matrices with high average girth
            cur_avg_girth = avg_girth(ex_g);

            if (cur_avg_girth >= expected_avg_girth)
            {
                sprintf(postfix, "%d", file_counter + 1);
                sprintf(file, "%s", codefile);
                sprintf(file + strlen(codefile), "%s", postfix);

                printf("# Matrix average girth %f \n", cur_avg_girth);
                printf("# Writing matrix to file %s ...\n", file);

                temp_f = open_file_std(file, "wb");

                if (temp_f == NULL)
                {
                    fprintf(stderr, "Can't create files to store shift keys");
                    exit(1);
                }

                for (i = 0; i < b_rows; i++)
                {
                    pre_col = -1;

                    for (e = mod2sparse_first_in_row(baseMatrix, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
                    {
                        for (j = pre_col + 1; j < e->col; j++)
                            fprintf(temp_f, "  -1");

                        fprintf(temp_f, "%4d", e->shift);

                        pre_col = e->col;
                    }

                    for (j = (pre_col + 1); j < b_cols; j++)
                        fprintf(temp_f, "  -1");

                    fprintf(temp_f, "\n");
                }

                if (ferror(temp_f) || fclose(temp_f))
                {
                    fprintf(stderr, "Error writing to shift keys to file %s\n", file);
                    exit(1);
                }

                file_counter++;
                index = 0;
            }
            else
            {
                index++;
            }

            free(ex_g->bit_node);
            free(ex_g->check_node);
            free(ex_g);
        }
        else
        {
            index++;
        }

        mod2sparse_free(baseMatrix);
        mod2sparse_free(exMatrix);

        baseMatrix = mod2sparse_allocate(b_rows, b_cols);
        exMatrix = mod2sparse_allocate(b_rows * ex_factor, b_cols * ex_factor);
    }

    mod2sparse_free(baseMatrix);
    mod2sparse_free(exMatrix);
}

int check_PHI_inv(mod2sparse *H, int m, int n, int p, int g)
{
    mod2sparse *PHI, *E, *B, *D;
    mod2sparse *EB;
    mod2dense *PHI_D, *PHI_DI;
    int M, N, inv;

    M = mod2sparse_rows(H);
    N = mod2sparse_cols(H);

    PHI = mod2sparse_allocate(g, g);
    E = mod2sparse_allocate(g, (m*p-g));
    B = mod2sparse_allocate((m*p-g), g);
    D = mod2sparse_allocate(g, g);
    EB = mod2sparse_allocate(g, g);
    PHI_D = mod2dense_allocate(g, g);
    PHI_DI = mod2dense_allocate(g, g);

    // Extract the submatrix E, T, B, D from parity check matrix H
    sub_matrix(H, E, (m*p-g), M-1, (n-m)*p+g, N-1);
    sub_matrix(H, B, 0, (m*p-g)-1, (n-m)*p, (n-m)*p+g-1);
    sub_matrix(H, D, (m*p-g), M-1, (n-m)*p, (n-m)*p+g-1);

    // Calculate: PHI = EB + D
    mod2sparse_multiply(E, B, EB);
    mod2sparse_add(EB, D, PHI);

    mod2sparse_to_dense(PHI, PHI_D);

    inv = mod2dense_invert(PHI_D, PHI_DI);

    mod2sparse_free(PHI);
    mod2sparse_free(E);
    mod2sparse_free(B);
    mod2sparse_free(D);
    mod2sparse_free(EB);
    mod2dense_free(PHI_D);
    mod2dense_free(PHI_DI);

    return inv; // 1: invertible, 0: non-invertible
}

int check_PHI_inv_pad( mod2sparse *H, int m, int n, int p, int g, int pad_bit)
{
    mod2sparse *PHI, *E, *B, *D;
    mod2sparse *EB;
    mod2dense *PHI_D, *PHI_DI;
    int M, N, inv;
    int mask_len;

    mask_len = p - pad_bit;
    M = m*p;
    N = n*p;

    PHI = mod2sparse_allocate(g, g);
    E = mod2sparse_allocate(g, (m*p-g-mask_len));
    B = mod2sparse_allocate((m*p-g-mask_len), g);
    D = mod2sparse_allocate(g, g);
    EB = mod2sparse_allocate(g, g);
    PHI_D = mod2dense_allocate(g, g);
    PHI_DI = mod2dense_allocate(g, g);

    // Extract the submatrix E, T, B, D from parity check matrix H
    sub_matrix(H, E, (m*p-g), M-1, (n-m)*p+g, N-mask_len-1);
    sub_matrix(H, B, 0, (m*p-g)-mask_len-1, (n-m)*p, (n-m)*p+g-1);
    sub_matrix(H, D, (m*p-g), M-1, (n-m)*p, (n-m)*p+g-1);

    // Calculate: PHI = EB + D
    mod2sparse_multiply(E, B, EB);
    mod2sparse_add(EB, D, PHI);

    mod2sparse_to_dense(PHI, PHI_D);

    inv = mod2dense_invert(PHI_D, PHI_DI);

    mod2sparse_free(PHI);
    mod2sparse_free(E);
    mod2sparse_free(B);
    mod2sparse_free(D);
    mod2sparse_free(EB);
    mod2dense_free(PHI_D);
    mod2dense_free(PHI_DI);

    return inv; // 1: invertible, 0: non-invertible
}

int check_PHI_inv_ibex(mod2sparse *H, int m, int n, int p, int g)
{
    mod2sparse *PHI;
    mod2dense *PHI_D, *PHI_DI;
    int M, N, inv;
    int i, j;
    mod2entry *e;

    M = mod2sparse_rows(H);
    N = mod2sparse_cols(H);

    PHI = mod2sparse_allocate(g, g);
    PHI_D = mod2dense_allocate(g, g);
    PHI_DI = mod2dense_allocate(g, g);

    // Extract the submatrix G from parity check matrix H
    sub_matrix(H, PHI, 0, g-1, (N-g), N-1);

    // Calculate: PHI inv
    mod2sparse_to_dense(PHI, PHI_D);

    inv = mod2dense_invert(PHI_D, PHI_DI);

    mod2sparse_free(PHI);
    mod2dense_free(PHI_D);
    mod2dense_free(PHI_DI);

    return inv; // 1: invertible, 0: non-invertible
}

int search_cyc4_mask(mod2sparse* basePCH, mod2entry* e, int** mask_matrix, int qc_size, int pad_bit)
{
    int num4 = 0;

    mod2entry *e2, *e3, *e4;

    if ((pad_bit==qc_size) && (mask_matrix[e->row][e->col]==1))
        return 0;

    for (e2 = mod2sparse_first_in_col(basePCH, e->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
    {
        if ((pad_bit==qc_size) && (mask_matrix[e2->row][e2->col]==1))
            continue;
        if ((mask_matrix[e->row][e->col]==1 && mask_matrix[e2->row][e2->col]==2)
            || (mask_matrix[e->row][e->col]==2 && mask_matrix[e2->row][e2->col]==1))
            continue;
        if (e2->row == e->row)
            continue;

        e3 = mod2sparse_first_in_row(basePCH, e2->row);
        e4 = mod2sparse_first_in_row(basePCH, e->row);
        while (!mod2sparse_at_end(e3) && !mod2sparse_at_end(e4))
        {
            if (e3->col > e4->col)       
                e4 = mod2sparse_next_in_row(e4);
            else if (e3->col < e4->col)
                e3 = mod2sparse_next_in_row(e3);
            else
            {
                if (e3->col == e2->col){}
                else if ((pad_bit == qc_size)
                        && (mask_matrix[e3->row][e3->col]==1 || mask_matrix[e4->row][e4->col]==1)){}
                else if ((mask_matrix[e3->row][e3->col]==1 && mask_matrix[e4->row][e4->col]==2)
                        || (mask_matrix[e3->row][e3->col]==2 && mask_matrix[e4->row][e4->col]==1)){}
                else if ((e->shift - e2->shift + e3->shift - e4->shift) % qc_size == 0)
                    num4++;

                e3 = mod2sparse_next_in_row(e3);
                e4 = mod2sparse_next_in_row(e4);
            }
        }
    }
    return num4;
}


int search_cyc6_mask(mod2sparse* basePCH, mod2entry* e, int** mask_matrix, int qc_size, int pad_bit)
{
    int num6 = 0;

    mod2entry *e2, *e3, *e4, *e5, *e6;

    if ((pad_bit==qc_size) && (mask_matrix[e->row][e->col]==1))
        return 0;

    for (e2 = mod2sparse_first_in_col(basePCH, e->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
    {
        if (e2->row == e->row)
            continue;
        if ((pad_bit==qc_size) && (mask_matrix[e2->row][e2->col]==1))
            continue;
        if ((mask_matrix[e->row][e->col]==1 && mask_matrix[e2->row][e2->col]==2)
            || (mask_matrix[e->row][e->col]==2 && mask_matrix[e2->row][e2->col]==1))
            continue;

        for (e3 = mod2sparse_first_in_row(basePCH, e2->row); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
        {
            if (e3->col == e2->col)
                continue;
            if ((pad_bit==qc_size) && (mask_matrix[e3->row][e3->col]==1))
                continue;
            
            for (e6 = mod2sparse_first_in_row(basePCH, e->row); !mod2sparse_at_end(e6); e6 = mod2sparse_next_in_row(e6))
            {
                if ((pad_bit==qc_size) && (mask_matrix[e6->row][e6->col]==1))
                    continue;
                if (e6->col == e3->col || e6->col == e->col)
                    continue;

                e4 = mod2sparse_first_in_col(basePCH, e3->col);
                e5 = mod2sparse_first_in_col(basePCH, e6->col);
                while (!mod2sparse_at_end(e4) && !mod2sparse_at_end(e5))
                {
                    if (e4->row > e5->row)       
                        e5 = mod2sparse_next_in_col(e5);
                    else if (e4->row < e5->row)
                        e4 = mod2sparse_next_in_col(e4);
                    else
                    {
                        if (e4->row == e3->row || e5->row == e6->row){}
                        else if ((pad_bit == qc_size)
                                && (mask_matrix[e4->row][e4->col]==1 || mask_matrix[e5->row][e5->col]==1)){}
                        else if ((mask_matrix[e3->row][e3->col]==1 && mask_matrix[e4->row][e4->col]==2)
                                || (mask_matrix[e3->row][e3->col]==2 && mask_matrix[e4->row][e4->col]==1)){}
                        else if ((mask_matrix[e5->row][e5->col]==1 && mask_matrix[e6->row][e6->col]==2)
                                || (mask_matrix[e5->row][e5->col]==2 && mask_matrix[e6->row][e6->col]==1)){}
                        else if ((e->shift - e2->shift + e3->shift - e4->shift + e5->shift - e6->shift) % qc_size == 0)
                            num6++;

                        e4 = mod2sparse_next_in_col(e4);
                        e5 = mod2sparse_next_in_col(e5);
                    }
                }
            }
        }
    }
    return num6;
}

int search_cyc4(mod2sparse* basePCH, mod2entry* e, int** fade_matrix, int b_rows, int qc_size, int pad_bit)
{
    int num4 = 0;

    mod2entry *e2, *e3, *e4;

    if ((pad_bit==qc_size) && (fade_matrix[e->row][e->col]==1))
        return 0;

    for (e2 = mod2sparse_first_in_col(basePCH, e->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
    {
        if ((pad_bit == qc_size) && (fade_matrix[e2->row][e2->col]==1))
            continue;
        if ((fade_matrix[e->row][e->col]==1 && e2->row == b_rows-1)
            || (e->row == b_rows-1 && fade_matrix[e2->row][e2->col]==1))
            continue;
        if (e2->row == e->row)
            continue;

        e3 = mod2sparse_first_in_row(basePCH, e2->row);
        e4 = mod2sparse_first_in_row(basePCH, e->row);
        while (!mod2sparse_at_end(e3) && !mod2sparse_at_end(e4))
        {
            if (e3->col > e4->col)       
                e4 = mod2sparse_next_in_row(e4);
            else if (e3->col < e4->col)
                e3 = mod2sparse_next_in_row(e3);
            else
            {
                if (e3->col == e2->col){}
                else if ((pad_bit == qc_size)
                        && (fade_matrix[e3->row][e3->col]==1 || fade_matrix[e4->row][e4->col]==1)){}
                else if ((fade_matrix[e3->row][e3->col]==1 && e4->row == b_rows-1)
                        || (e3->row == b_rows-1 && fade_matrix[e4->row][e4->col]==1)){}
                else if ((e->shift - e2->shift + e3->shift - e4->shift) % qc_size == 0)
                    num4++;

                e3 = mod2sparse_next_in_row(e3);
                e4 = mod2sparse_next_in_row(e4);
            }
        }
    }
    return num4;
}

int search_cyc6(mod2sparse* basePCH, mod2entry* e, int** fade_matrix, int b_rows, int qc_size, int pad_bit)
{
    int num6 = 0;

    mod2entry *e2, *e3, *e4, *e5, *e6;

    if ((pad_bit==qc_size) && (fade_matrix[e->row][e->col]==1))
        return 0;

    for (e2 = mod2sparse_first_in_col(basePCH, e->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
    {
        if (e2->row == e->row)
            continue;
        if ((pad_bit==qc_size) && (fade_matrix[e2->row][e2->col]==1))
            continue;
        if ((fade_matrix[e->row][e->col]==1 && e2->row == b_rows-1)
            || (e->row == b_rows-1 && fade_matrix[e2->row][e2->col]==1))
            continue;

        for (e3 = mod2sparse_first_in_row(basePCH, e2->row); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
        {
            if (e3->col == e2->col)
                continue;
            if ((pad_bit==qc_size) && (fade_matrix[e3->row][e3->col]==1))
                continue;
            
            for (e6 = mod2sparse_first_in_row(basePCH, e->row); !mod2sparse_at_end(e6); e6 = mod2sparse_next_in_row(e6))
            {
                if ((pad_bit==qc_size) && (fade_matrix[e6->row][e6->col]==1))
                    continue;
                if (e6->col == e->col || e3->col == e6->col)
                    continue;

                e4 = mod2sparse_first_in_col(basePCH, e3->col);
                e5 = mod2sparse_first_in_col(basePCH, e6->col);
                while (!mod2sparse_at_end(e4) && !mod2sparse_at_end(e5))
                {
                    if (e4->row > e5->row)       
                        e5 = mod2sparse_next_in_col(e5);
                    else if (e4->row < e5->row)
                        e4 = mod2sparse_next_in_col(e4);
                    else
                    {
                        if (e4->row == e3->row || e5->row == e6->row){}
                        else if ((pad_bit == qc_size)
                                && (fade_matrix[e4->row][e4->col]==1 || fade_matrix[e5->row][e5->col]==1)){}
                        else if ((fade_matrix[e3->row][e3->col]==1 && e4->row == b_rows-1)
                                || (e3->row == b_rows-1 && fade_matrix[e4->row][e4->col]==1)){}
                        else if ((fade_matrix[e5->row][e5->col]==1 && e6->row == b_rows-1)
                                || (e5->row == b_rows-1 && fade_matrix[e6->row][e6->col]==1)){}
                        else if ((e->shift - e2->shift + e3->shift - e4->shift + e5->shift - e6->shift) % qc_size == 0)
                            num6++;

                        e4 = mod2sparse_next_in_col(e4);
                        e5 = mod2sparse_next_in_col(e5);
                    }
                }
            }
        }
    }
    return num6;
}

int cyc4_cal_t(int** matrix, int** mask_matrix, int row, int col, int qc_size, int pad_bit)
{
    int num4 = 0;
    int l0, l1, k0, k1;

    int drop_len = pad_bit;
    int mask_len = qc_size - pad_bit;
    int max_mask_len = drop_len>mask_len ? drop_len : mask_len;

    for (k0 = 0; k0<col; k0++)
    {
        for (l0 = 0; l0<row; l0++)
        {
            if ((matrix[l0][k0] == -1) || (pad_bit==qc_size) && (mask_matrix[l0][k0]==1))
                continue;
            for (l1=l0+1; l1<row; l1++)
            {
                if ((matrix[l1][k0] == -1) || (pad_bit==qc_size) && (mask_matrix[l1][k0]==1))
                    continue;
                if ((mask_matrix[l0][k0]==1 && mask_matrix[l1][k0]==2)
                    || (mask_matrix[l0][k0]==2 && mask_matrix[l1][k0]==1))
                    continue;

                for (k1 = k0+1; k1<col; k1++)
                {
                    if ((matrix[l0][k1] == -1) || (matrix[l1][k1] == -1)
                        || (pad_bit==qc_size) && (mask_matrix[l0][k1]==1)
                        || (pad_bit==qc_size) && (mask_matrix[l1][k1]==1))
                        continue;
                    if ((mask_matrix[l0][k1]==1 && mask_matrix[l1][k1]==2)
                        || (mask_matrix[l0][k1]==2 && mask_matrix[l1][k1]==1))
                        continue;

                    if ((matrix[l0][k0] - matrix[l1][k0] + matrix[l1][k1] - matrix[l0][k1]) % qc_size == 0)
                        num4++;
                }
            }
        }
    }
    return num4;
}

int cyc6_cal_t(int** matrix, int** mask_matrix, int row, int col, int qc_size, int pad_bit)
{
    int num6 = 0;
    int l0, l1, l2, k0, k1, k2;

    for (k0 = 0; k0<col; k0++)
    {
        for (l0 = 0; l0<row; l0++)
        {
            if ((matrix[l0][k0] == -1) || (pad_bit==qc_size) && (mask_matrix[l0][k0]==1))
                continue;

            for (l1=l0+1; l1<row; l1++)
            {
                if ((matrix[l1][k0] == -1) || (pad_bit==qc_size) && (mask_matrix[l1][k0]==1))
                    continue;
                if ((mask_matrix[l0][k0]==1 && mask_matrix[l1][k0]==2)
                    || (mask_matrix[l0][k0]==2 && mask_matrix[l1][k0]==1))
                    continue;

                for (k1 = k0+1; k1<col; k1++)
                {
                    if ((matrix[l1][k1] == -1) || (pad_bit==qc_size) && (mask_matrix[l1][k1]==1))
                        continue;

                    for (l2=0; l2<row; l2++)
                    {
                        if ((l2==l0) || (l2==l1))
                            continue;
                        if ((matrix[l2][k1] == -1) || (pad_bit==qc_size) && (mask_matrix[l2][k1]==1))
                            continue;
                        if ((mask_matrix[l1][k1]==1 && mask_matrix[l2][k1]==2)
                            || (mask_matrix[l1][k1]==2 && mask_matrix[l2][k1]==1))
                            continue;

                        for (k2=k0+1; k2<col; k2++)
                        {
                            if (k2==k1)
                                continue;
                            if ((matrix[l0][k2] == -1) || (matrix[l2][k2] == -1)
                                || (pad_bit==qc_size) && (mask_matrix[l0][k2]==1)
                                || (pad_bit==qc_size) && (mask_matrix[l2][k2]==1))
                                continue;
                            if ((mask_matrix[l0][k2]==1 && mask_matrix[l2][k2]==2)
                                || (mask_matrix[l0][k2]==2 && mask_matrix[l2][k2]==1))
                                continue;

                            if ((matrix[l0][k0] - matrix[l1][k0] + matrix[l1][k1] - matrix[l2][k1]
                                + matrix[l2][k2] - matrix[l0][k2]) % qc_size == 0)
                                num6++;
                        }
                    }
                }
            }
        }
    }
    return num6;
}

// print mask matrix
void print_matrix_to_file(int** matrix, int row, int col, char* matrix_file)
{
    FILE *savefile = fopen(matrix_file, "w");
    if (savefile == NULL)
    {
        printf("It is a invalid matrix savefile!!!\n");
        exit(-1);
    }

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            fprintf(savefile, "%4d", matrix[i][j]);
        }
        fprintf(savefile, "\n");
    }
    fclose(savefile);
}

void ldpc_rd_phck(char *pchk_file, mod2sparse *base_b, mod2sparse *base_e, mod2sparse *base_d)
{
    mod2entry *e;
    FILE *fp;
    int col_shift;

    fp = fopen(pchk_file, "r");

    if (fp == 0)
    {
        printf("[LDPC] Error: Wrong parity check matrix file name\n");
        exit(0);
    }

    for (int i=0; i<11; i++)
    {
        for (int j=0; j<76; j++)
        {
            fscanf(fp, "%d", &col_shift);

            if (col_shift > 0)
            {
                if (i>=0 && i<=3 && j>=65 && j<=70)
                {
                    e = mod2sparse_insert(base_b, i, j-65);
                    e->shift = col_shift;
                }
                if (i>=5 && j>=71 && j<=74)
                {
                    e = mod2sparse_insert(base_e, i-5, j-71);
                    e->shift = col_shift;
                }
                if (i>=5 && j>=65 && j<=70)
                {
                    e = mod2sparse_insert(base_d, i-5, j-65);
                    e->shift = col_shift;
                }
            }
        }
    }
    fclose(fp);
    printf("[LDPC] H matrix porting ready!\n");
}