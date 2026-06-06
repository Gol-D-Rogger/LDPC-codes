#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "bipartite_graph.h"
#include "queue.h"
#include "intio.h"
#include "gen_ldpc.h"

node_entry trace_back(node_entry u, int n)
{
    int i;
    for(i = 0; i < n; i++)
        u = *(u.parent);
    return u;
}

int count_cyc_deg(node_entry u, node_entry v, int n)
{
    int i, cycle_degree;
    cycle_degree = u.node_degree;
    for(i = 0; i < n; i++)
    {
        u = *(u.parent);
        cycle_degree += u.node_degree;
    }
    cycle_degree += v.node_degree;
    n-=2;
    for(i = 0; i < n; i++)
    {
        v = *(v.parent);
        cycle_degree += v.node_degree;
    }
    return cycle_degree;
}

int check_cycle(tanner_graph *g, int root_index, int length)
{
    int M, N, i;
    int current_index, tmp_index;
    node_entry *current_node, *tmp_node;
    node_entry p_node1, p_node2; // used for tracing back
    mod2entry *e;
    Queue Q;

    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    Q = CreateQueue(M + N);

    // Initialize check nodes
    for (i = 0; i < M; i++)
    {
        g->check_node[i].index = i;
        g->check_node[i].color = white;
        g->check_node[i].parent = NULL;
        g->check_node[i].depth = -1;
    }

    // Initialize variable nodes
    for (i = 0; i < N; i++)
    {
        g->bit_node[i].index = i + M;
        g->bit_node[i].color = white;
        g->bit_node[i].parent = NULL;
        g->bit_node[i].depth = -1;
    }

    g->check_node[root_index].color = gray;
    g->check_node[root_index].depth = 0;
    Enqueue(root_index, Q);

    while (!IsEmpty(Q))
    {
        current_index = FrontAndDequeue(Q);

        if (current_index < M) // current node is a check node
        {
            current_node = &(g->check_node[current_index]);
            e = mod2sparse_first_in_row(g->matrix, current_index);

            while (!mod2sparse_at_end(e))
            {
                tmp_index = mod2sparse_col(e);
                tmp_node = &(g->bit_node[tmp_index]);

                if (current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch (tmp_node->color)
                    {
                    case white:
                        tmp_node->color = gray;
                        tmp_node->parent = current_node;
                        tmp_node->depth = current_node->depth + 1;
                        Enqueue(tmp_node->index, Q);
                        break;

                    case gray:
                        if (2 * tmp_node->depth <= length)
                        {
                            p_node1 = trace_back(*current_node, current_node->depth - 1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                            if (p_node1.index != p_node2.index)
                            {
                                DisposeQueue(Q);
                                return 1; // find cycle with given length passing root
                            }
                        }
                        // else the cycle doesn't pass the root
                        break;

                    case black:
                        break;

                    default:
                        fprintf(stderr, "Invalid node color!\n");
                    }
                }

                e = mod2sparse_next_in_row(e);
            }
            current_node->color = black;
        }

        if (current_index >= M) // current node is a variable node
        {
            current_index -= M;
            current_node = &(g->bit_node[current_index]);

            // get the first check node connected to this bit node
            e = mod2sparse_first_in_col(g->matrix, current_index);

            while (!mod2sparse_at_end(e))
            {
                tmp_index = mod2sparse_row(e);
                tmp_node = &(g->check_node[tmp_index]);

                if (current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch (tmp_node->color)
                    {
                    case white:
                        tmp_node->color = gray;
                        tmp_node->parent = current_node;
                        tmp_node->depth = current_node->depth + 1;
                        Enqueue(tmp_node->index, Q);
                        break;

                    case gray:
                        if (2 * tmp_node->depth <= length)
                        {
                            p_node1 = trace_back(*current_node, current_node->depth - 1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                            if (p_node1.index != p_node2.index)
                            {
                                DisposeQueue(Q);
                                return 1;
                            }
                        }
                        // else the cycle doesn't pass the root
                        break;
                    }
                }

                e = mod2sparse_next_in_col(e);
            }
            current_node->color = black;
        }

        if (2 * current_node->depth > length)
        {
            DisposeQueue(Q);
            return 0; // no cycle with length <= n
        }
    }

    DisposeQueue(Q);
    return 0; // no cycle with length <= n
}

int check_cycle_deg(tanner_graph *g, int root_index, int length)
{
    int M, N, i;
    int current_index, tmp_index;
    int cyc_deg;
    node_entry *current_node, *tmp_node;
    node_entry p_node1, p_node2;
    mod2entry *e;
    Queue Q;

    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    Q = CreateQueue(M + N);

    // Initialize check nodes
    for (i = 0; i < M; i++)
    {
        g->check_node[i].index = i;
        g->check_node[i].color = white;
        g->check_node[i].parent = NULL;
        g->check_node[i].depth = -1;
    }

    // Initialize variable nodes
    for (i = 0; i < N; i++)
    {
        g->bit_node[i].index = i + M;
        g->bit_node[i].color = white;
        g->bit_node[i].parent = NULL;
        g->bit_node[i].depth = -1;
    }

    g->check_node[root_index].color = gray;
    g->check_node[root_index].depth = 0;
    Enqueue(root_index, Q);

    while (!IsEmpty(Q))
    {
        current_index = FrontAndDequeue(Q);

        if (current_index < M) // current node is a check node
        {
            current_node = &(g->check_node[current_index]);
            e = mod2sparse_first_in_row(g->matrix, current_index);

            while (!mod2sparse_at_end(e))
            {
                tmp_index = mod2sparse_col(e);
                tmp_node = &(g->bit_node[tmp_index]);

                if (current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch (tmp_node->color)
                    {
                    case white:
                        tmp_node->color = gray;
                        tmp_node->parent = current_node;
                        tmp_node->depth = current_node->depth + 1;
                        Enqueue(tmp_node->index, Q);
                        break;

                    case gray:
                        if (2 * tmp_node->depth <= length)
                        {
                            p_node1 = trace_back(*current_node, current_node->depth - 1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                            if (p_node1.index != p_node2.index)
                            {
                                DisposeQueue(Q);
                                return 1;
                            }
                        }
                        // check cycle degree constraint; only consider when cycle length is 8 or less
                        if ((2*tmp_node->depth == length+2) && (length <= 6))
                        {
                            p_node1 = trace_back(*current_node, current_node->depth - 1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                            if (p_node1.index != p_node2.index)
                            {
                                cyc_deg = count_cyc_deg(*tmp_node, *current_node, tmp_node->depth);
                                if (cyc_deg < 3 * tmp_node->depth)
                                {
                                    DisposeQueue(Q);
                                    return 1;
                                }
                            }
                        }
                        // else the cycle doesn't pass the root
                        break;

                    case black:
                        break;
                    default:
                        fprintf(stderr, "Invalid node color!\n");
                    }
                }

                e = mod2sparse_next_in_row(e);
            }
            current_node->color = black;
        }

        if (current_index >= M) // current node is a variable node
        {
            current_index -= M;
            current_node = &(g->bit_node[current_index]);

            // get the first check node connected to this bit node
            e = mod2sparse_first_in_col(g->matrix, current_index);

            while (!mod2sparse_at_end(e))
            {
                tmp_index = mod2sparse_row(e);
                tmp_node = &(g->check_node[tmp_index]);

                if (current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch (tmp_node->color)
                    {
                    case white:
                        tmp_node->color = gray;
                        tmp_node->parent = current_node;
                        tmp_node->depth = current_node->depth + 1;
                        Enqueue(tmp_node->index, Q);
                        break;

                    case gray:
                        if (2 * tmp_node->depth <= length)
                        {
                            p_node1 = trace_back(*current_node, current_node->depth - 1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                            if (p_node1.index != p_node2.index)
                            {
                                DisposeQueue(Q);
                                return 1;
                            }
                        }
                        // check cycle degree constraint; only consider when cycle length is 8 or less
                        if ((2*tmp_node->depth == length+2) && (length <= 6))
                        {
                            p_node1 = trace_back(*current_node, current_node->depth - 1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                            if (p_node1.index != p_node2.index)
                            {
                                cyc_deg = count_cyc_deg(*tmp_node, *current_node, tmp_node->depth);
                                if (cyc_deg < 3 * tmp_node->depth) // node with degree 2 can't in cycle with length>=6
                                {
                                    DisposeQueue(Q);
                                    return 1; // fail to meet cycle degree constraint
                                }
                            }
                        }
                        // else the cycle doesn't pass the root
                        break;
                    }
                }

                e = mod2sparse_next_in_col(e);
            }
            current_node->color = black;
        }
        if (2 * current_node->depth > length)
        {
            DisposeQueue(Q);
            return 0; // no cycle with length <= n
        }
    }
                
    DisposeQueue(Q);
    return 0; // no cycle with length <= n
}

/**************************************************************
The idea is:
1) Pick one node as root, count all the cycles(<MAX_CYCLE_LENGTH) passing
   through this node, using BFS to traverse the nodes(vetex) and edges.
2) Isolate this node by removing the edges connected to it. Go back to 1) until
    all nodes are isolated.
***************************************************************/

void count_cycles(tanner_graph *g, int *cycle)
{
    int M, N, i, j;
    int round;
    mod2entry *e;
    int tmp_col, tmp_row;
    tanner_graph g_cpy;

    g_cpy = *g;

    round = 0;
    M = mod2sparse_rows(g_cpy.matrix);
    N = mod2sparse_cols(g_cpy.matrix);

    for (i=0; i < MAX_CYCLE_LENGTH; i++)
        cycle[i] = 0;

    // Initialize check nodes
    for(i = 0; i < M; i++)
    {
        g_cpy.check_node[i].index = i;
        g_cpy.check_node[i].color = white;
        g_cpy.check_node[i].parent = NULL;
        g_cpy.check_node[i].depth = -1;
    }

    // Initialize variable nodes
    for(i = 0; i < N; i++)
    {
        g_cpy.bit_node[i].index = i+M;
        g_cpy.bit_node[i].color = white;
        g_cpy.bit_node[i].parent = NULL;
        g_cpy.bit_node[i].depth = -1;
    }

    // This loop counts cycles passing through check node; isolated it after it's done
    for(i = 0; i < M; i++)
    {
        if (round!=0)
        {
            // re-initialize check nodes
            for(j = i; j < M; j++)
            {
                g_cpy.check_node[j].color = white;
                g_cpy.check_node[j].parent = NULL;
                g_cpy.check_node[j].depth = -1;
            }

            // re-initialize variable nodes
            for(j = 0; j < N; j++)
            {
                g_cpy.bit_node[j].color = white;
                g_cpy.bit_node[j].parent = NULL;
                g_cpy.bit_node[j].depth = -1;
            }
        }
        round++;

        // BFS traverse ...
        BFS_traverse(g, &g_cpy.check_node[i], cycle);

        // isolate this check node by removing all the edges connected to it
        tmp_row = i;
        e = mod2sparse_first_in_row(g_cpy.matrix, tmp_row);
        while(!mod2sparse_at_end(e))
        {
            tmp_col = mod2sparse_col(e);
            e = mod2sparse_next_in_row(e);
            mod2sparse_delete(g_cpy.matrix, mod2sparse_find(g_cpy.matrix, tmp_row, tmp_col));
        }
    }
}

void enumerate_cycle(tanner_graph *g, cyc4_entry *cyc4_rec)
{
    int i, j;
    int M, N;
    int round;
    mod2entry *e;
    int tmp_col, tmp_row;
    int cyc4_cnt = 0;

    round = 0;
    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    // Initialize check nodes
    for(i = 0; i < M; i++)
    {
        g->check_node[i].index = i;
        g->check_node[i].color = white;
        g->check_node[i].parent = NULL;
        g->check_node[i].depth = -1;
    }

    // Initialize variable nodes
    for(i = 0; i < N; i++)
    {
        g->bit_node[i].index = i+M;
        g->bit_node[i].color = white;
        g->bit_node[i].parent = NULL;
        g->bit_node[i].depth = -1;
    }

    // This loop counts cycles passing through check node; isolated it after it's done
    for(i = 0; i < M; i++)
    {
        if(round != 0)
        {
            // re-initialize check nodes
            for(j = i; j < M; j++)
            {
                g->check_node[j].color = white;
                g->check_node[j].parent = NULL;
                g->check_node[j].depth = -1;
            }

            // re-initialize variable nodes
            for(j = 0; j < N; j++)
            {
                g->bit_node[j].color = white;
                g->bit_node[j].parent = NULL;
                g->bit_node[j].depth = -1;
            }
        }
        round++;

        // BFS traverse ...
        BFS_traverse_enum(g, &g->check_node[i], cyc4_rec, &cyc4_cnt);

        // isolate this check node by removing all the edges connected to it
        tmp_row = i;
        e = mod2sparse_first_in_row(g->matrix, tmp_row);
        while(!mod2sparse_at_end(e))
        {
            tmp_col = mod2sparse_col(e);
            e = mod2sparse_next_in_row(e);
            mod2sparse_delete(g->matrix, mod2sparse_find(g->matrix, tmp_row, tmp_col));
        }
    }
}

double BFS_traverse(tanner_graph *g, node_entry *root, int *cycle)
{
    double lambda = 0.5;
    double connectivity = 0.0;

    int M, N, i;
    int current_index, tmp_index, cyc_len;
    node_entry *current_node, *tmp_node;
    node_entry p_node1, p_node2;
    mod2entry *e;
    Queue Q;

    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    Q = CreateQueue(M+N);

    root->color = gray;
    root->depth = 0;
    Enqueue(root->index, Q);

    while( !IsEmpty(Q))
    {
        current_index = FrontAndDequeue(Q);

        if(current_index < M) // current node is a check node
        {
            current_node = &(g->check_node[current_index]);
            if (current_node->depth > MAX_CYCLE_LENGTH/2)
            {
                DisposeQueue(Q);
                return 0;
            }

            // get the first bit node connected to this check node
            e = mod2sparse_first_in_row(g->matrix, current_index);

            while(!mod2sparse_at_end(e))
            {
                tmp_index = mod2sparse_col(e);
                tmp_node = &(g->bit_node[tmp_index]);

                if(current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch(tmp_node->color)
                    {
                        case white:
                            tmp_node->color = gray;
                            tmp_node->parent = current_node;
                            tmp_node->depth = current_node->depth + 1;
                            Enqueue(tmp_node->index, Q);
                            break;

                        case gray:
                            p_node1 = trace_back(*current_node, current_node->depth-1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth-1);
                            if(p_node1.index != p_node2.index)
                            {
                                cyc_len = 2 * tmp_node->depth;
                                cycle[cyc_len]++;
                            }
                            // else the cycle does not pass the root
                            break;

                        case black:
                            break;
                        default:
                            fprintf(stderr, "Invalid node color!\n");
                    }
                }
                e = mod2sparse_next_in_row(e);
            }
            current_node->color = black;
        }
        
        if (current_index >= M) // current node is a variable node
        {
            current_index -= M;
            current_node = &(g->bit_node[current_index]);

            if (current_node->depth > MAX_CYCLE_LENGTH/2)
            {
                DisposeQueue(Q);
                return 0;
            }

            // get the first check node connected to this bit node
            e = mod2sparse_first_in_col(g->matrix, current_index);

            while (!mod2sparse_at_end(e))
            {
                // retrieve the index of the check node
                tmp_index = mod2sparse_row(e);
                tmp_node = &(g->check_node[tmp_index]);

                if (current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch (tmp_node->color)
                    {
                    case white:
                        tmp_node->color = gray;
                        tmp_node->parent = current_node;
                        tmp_node->depth = current_node->depth + 1;
                        connectivity += (pow(lambda, tmp_node->depth));
                        Enqueue(tmp_node->index, Q);
                        break;

                    case gray:
                        p_node1 = trace_back(*current_node, current_node->depth - 1);
                        p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                        if (p_node1.index != p_node2.index)
                        {
                            cyc_len = 2 * tmp_node->depth;
                            cycle[cyc_len]++;
                        }
                        // else the cycle does not pass the root
                        break;

                    case black:
                        break;
                    default:
                        fprintf(stderr, "Invalid node color!\n");
                    }
                }

                // get next check node connected to this bit node
                e = mod2sparse_next_in_col(e);
            }
            current_node->color = black;
        }
    }
    DisposeQueue(Q);
    return connectivity;
}

int BFS_traverse_enum(tanner_graph *g, node_entry *root, cyc4_entry *cyc4_rec, int *cyc4_cnt)
{
    int M, N;
    int current_index, tmp_index;
    node_entry *current_node, *tmp_node;
    node_entry p_node1, p_node2;
    mod2entry *e;
    Queue Q;

    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    Q = CreateQueue(M + N);

    root->color = gray;
    root->depth = 0;
    Enqueue(root->index, Q);

    while (!IsEmpty(Q))
    {
        current_index = FrontAndDequeue(Q);

        if (current_index < M) // current node is a check node
        {
            current_node = &(g->check_node[current_index]);
            if (current_node->depth > MAX_CYCLE_LENGTH / 2)
            {
                DisposeQueue(Q);
                return 0;
            }

            // get the first bit node connected to this check node
            e = mod2sparse_first_in_row(g->matrix, current_index);

            while (!mod2sparse_at_end(e))
            {
                tmp_index = mod2sparse_col(e);
                tmp_node = &(g->bit_node[tmp_index]);

                if (current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch (tmp_node->color)
                    {
                    case white:
                        tmp_node->color = gray;
                        tmp_node->parent = current_node;
                        tmp_node->depth = current_node->depth + 1;
                        Enqueue(tmp_node->index, Q);
                        break;

                    case gray:
                    
                        break;

                    case black:
                        break;

                    default:
                        fprintf(stderr, "Invalid node color!\n");
                    }
                }
                e = mod2sparse_next_in_row(e);
            }
            current_node->color = black;
        }

        if (current_index >= M) // current node is a variable node
        {
            current_index -= M;
            current_node = &(g->bit_node[current_index]);

            if (current_node->depth > MAX_CYCLE_LENGTH / 2)
            {
                DisposeQueue(Q);
                return 0;
            }

            // get the first check node connected to this bit node
            e = mod2sparse_first_in_col(g->matrix, current_index);

            while (!mod2sparse_at_end(e))
            {
                // retrieve the index of the check node
                tmp_index = mod2sparse_row(e);
                tmp_node = &(g->check_node[tmp_index]);

                if (current_node->parent == NULL || tmp_index != current_node->parent->index)
                {
                    switch (tmp_node->color)
                    {
                    case white:
                        tmp_node->color = gray;
                        tmp_node->parent = current_node;
                        tmp_node->depth = current_node->depth + 1;
                        Enqueue(tmp_node->index, Q);
                        break;

                    case gray:
                        if (current_node->depth == 1)
                        {
                            p_node1 = trace_back(*current_node, current_node->depth - 1);
                            p_node2 = trace_back(*tmp_node, tmp_node->depth - 1);
                            if (p_node1.index != p_node2.index)
                            {
                                cyc4_rec[*cyc4_cnt].r1 = root->index;
                                cyc4_rec[*cyc4_cnt].r2 = tmp_node->index;
                                cyc4_rec[*cyc4_cnt].c1 = p_node2.index - M;
                                cyc4_rec[*cyc4_cnt].c2 = p_node1.index - M;

                                (*cyc4_cnt)++;
                            }
                        }
                        break;

                    case black:
                        break;

                    default:
                        fprintf(stderr, "Invalid node color!\n");
                    }
                }

                // get next check node connected to this bit node
                e = mod2sparse_next_in_col(e);
            }
            current_node->color = black;
        }
    }
    DisposeQueue(Q);
    return 1;
}

void DFS_count_cycles(tanner_graph *g,  // Tanner graph associated with the parity-check matrix
                      int *cycle)       // place to store the cycles counter
{
    int M, N, i;
    int wh, gr, bl; // for debug

    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    for (i = 0; i < MAX_CYCLE_LENGTH; i++)
        cycle[i] = 0;

    // Initialize check nodes
    for (i = 0; i < M; i++)
    {
        g->check_node[i].index = i;
        g->check_node[i].color = white;
        g->check_node[i].parent = NULL;
        g->check_node[i].depth = -1;
    }

    // Initialize variable nodes
    for (i = 0; i < N; i++)
    {
        g->bit_node[i].index = i + M;   // let the index of bit node following check node index, so
                                        // all the nodes in graph have a unique identified index.
        g->bit_node[i].color = white;
        g->bit_node[i].parent = NULL;
        g->bit_node[i].depth = -1;
    }

    // start from root node (just use check node as root, because every bit node will
    // be connected at least to 1 check node, if it's not isolated)
    for (i = 0; i < M; i++)
    {
        if (g->check_node[i].color == white)
        {
            g->check_node[i].depth = 0;
            g->check_node[i].color = gray;
            DFS_visit(g, &g->check_node[i], cycle);
        }
    }

    // for debug
    wh = gr = bl = 0;
    for (i = 0; i < M; i++)
    {
        switch (g->check_node[i].color)
        {
        case white:
            wh++;
            break;
        case gray:
            gr++;
            break;
        case black:
            bl++;
            break;
        default:
            fprintf(stderr, "There is node with invalid color!\n");
        }
    }

    for (i = 0; i < N; i++)
    {
        switch (g->bit_node[i].color)
        {
        case white:
            wh++;
            break;
        case gray:
            gr++;
            break;
        case black:
            bl++;
            break;
        default:
            fprintf(stderr, "There is node with invalid color!\n");
        }
    }

    fprintf(stderr, "white node = %d\n, gray node = %d\n, black node = %d\n", wh, gr, bl);
}

void DFS_visit(tanner_graph *g, node_entry *u, int *cycle)
{
    mod2entry *e;
    int tmp_col, tmp_row;
    node_entry *tmp_bit_node, *tmp_check_node;
    int cycle_length;
    int M, N;

    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    if (u->index < M && u->index >= 0) // u is a check node
    {
        // get the first bit node connected to this check node
        e = mod2sparse_first_in_row(g->matrix, u->index);

        while (!mod2sparse_at_end(e))
        {
            tmp_col = mod2sparse_col(e);
            tmp_bit_node = &(g->bit_node[tmp_col]);

            if (u->parent == NULL || tmp_col != (u->parent->index - M))
            {
                switch (tmp_bit_node->color)
                {
                case white:
                    tmp_bit_node->color = gray;
                    tmp_bit_node->parent = u;
                    tmp_bit_node->depth = u->depth + 1;
                    DFS_visit(g, tmp_bit_node, cycle);
                    break;

                case gray:
                    cycle_length = u->depth - tmp_bit_node->depth + 1;
                    if (cycle_length < 0)
                        fprintf(stderr, "cycle length cannot be negative!\n");
                    if (cycle_length < MAX_CYCLE_LENGTH)
                        cycle[cycle_length]++;
                    break;

                case black:
                    break;
                default:
                    fprintf(stderr, "Invalid node color!\n");
                }
            }
            e = mod2sparse_next_in_row(e);
        }
        u->color = black;
        return;
    }

    else if (u->index >= M && u->index < M + N) // u is a bit node
    {
        // get the first check node connected to this bit node
        e = mod2sparse_first_in_col(g->matrix, u->index - M);

        while (!mod2sparse_at_end(e))
        {
            tmp_row = mod2sparse_row(e);
            tmp_check_node = &(g->check_node[tmp_row]);

            if (u->parent == NULL || tmp_row != u->parent->index)
            {
                switch (tmp_check_node->color)
                {
                case white:
                    tmp_check_node->color = gray;
                    tmp_check_node->depth = u->depth + 1;
                    tmp_check_node->parent = u;
                    DFS_visit(g, tmp_check_node, cycle);
                    break;

                case gray:
                    cycle_length = u->depth - tmp_check_node->depth + 1;
                    if (cycle_length < 0)
                        fprintf(stderr, "cycle length cannot be negative!\n");
                    if (cycle_length < MAX_CYCLE_LENGTH)
                        cycle[cycle_length]++;
                    break;

                case black:
                    // according to the property of DFS, there should not be any back edge to black node
                    break;
                default:
                    fprintf(stderr, "Invalid node color!\n");
                }
            }
            e = mod2sparse_next_in_col(e);
        }
        u->color = black;
        return;
    }

    else
        fprintf(stderr, "Node index out of range!\n");

    return;
}

void check_phck_cycles(char *file, char *log, int *cycle)
{
    mod2sparse *readMatrix;
    tanner_graph *g;
    int i;
    FILE *f, *l;
    char logline[256];

    f = fopen(file, "rb");
    if (f == NULL)
    {
        fprintf(stderr, "Cannot open file %s to read!\n", file);
        exit(1);
    }
    if (intio_read(f) != ('P'<<8)+0x80)
    {
        fprintf(stderr, "File %s is not in the correct format!\n", file);
        exit(1);
    }

    readMatrix = mod2sparse_read(f);
    mod2sparse_print(stdout, readMatrix);

    g = chk_alloc(1, sizeof(tanner_graph));
    g->bit_node = chk_alloc(mod2sparse_cols(readMatrix), sizeof(node_entry));
    g->check_node = chk_alloc(mod2sparse_rows(readMatrix), sizeof(node_entry));
    g->matrix = readMatrix;

    if ((l = fopen(log, "a")) == NULL)
        printf("Fail to write %s\n", log);
    count_cycles(g, cycle);
    for (i=0; i<30; i++)
        printf("Cycle[%d] = %d\n", i, cycle[i]);

    fputs(file, l);
    fputs("\n", l);
    for (i = 0; i < 30; i++)
        sprintf(logline + 4*i, "%4d", i);
    fputs(logline, l);
    fputs("\n", l);
    for (i = 0; i < 30; i++)
        sprintf(logline + 4*i, "%4d", cycle[i]);
    fputs(logline, l);
    fputs("\n", l);
    
    free(g->bit_node);
    free(g->check_node);
    free(g);
    free(readMatrix);
    fclose(f);
    fclose(l);
}

double cycle_effect(int *cycle, int n)
{
    int i;
    double ce, alpha;
    alpha = 1/sqrt(30);
    ce = 0.0;

    for(i = 0; i < n; i++)
        ce += cycle[i] * pow(alpha, (i-4));

    return ce;
}

float avg_girth(tanner_graph *g)
{
    int i, j, M, N;
    float g_sum, n_sum;
    int girth[7];

    M = mod2sparse_rows(g->matrix);
    N = mod2sparse_cols(g->matrix);

    for (i=0; i<7; i++)
        girth[i] = 0;

    for (j=0; j<M; j++)
    {
        if (check_cycle(g, j, 4))
            girth[0]++;
        else
        {
            if(check_cycle(g, j, 6))
                girth[1]++;
            else
            {
                if (check_cycle(g, j, 8))
                    girth[2]++;
                else
                {
                    if (check_cycle(g, j, 10))
                        girth[3]++;
                    else
                    {
                        if (check_cycle(g, j, 12))
                            girth[4]++;
                        else
                        {
                            if (check_cycle(g, j, 14))
                                girth[5]++;
                            else
                            {
                                if (check_cycle(g, j, 16))
                                girth[6]++;
                            }
                        }
                    }
                }
            }
        }
    }

    g_sum = 0.0;
    n_sum = 0.0;
    for (i=0; i<7; i++)
    {
        g_sum += girth[i] * (2*i + 4);
        n_sum += girth[i];
    }

    return g_sum/n_sum;
}


int count_4cyc(mod2sparse *matrix)
{
    mod2entry *e1, *e2, *e3, *e4;
    int i, j, row_weight;
    int cnt = 0;

    for (i=0; i < matrix->n_rows; i++)
    {
        row_weight = 0;
        for (e1 = mod2sparse_first_in_row(matrix, i); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_row(e1))
        {   
            row_weight++;
            for (e2 = mod2sparse_first_in_col(matrix, e1->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if (e2->row == i)
                    continue;
                for (e3 = mod2sparse_first_in_row(matrix, e2->row); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
                {
                    if (e3->col == e1->col || e3->col > e1->col)
                        continue;
                    for (e4 = mod2sparse_first_in_col(matrix, e3->col); !mod2sparse_at_end(e4); e4 = mod2sparse_next_in_col(e4))
                    {
                        if (e4->row == i)
                            cnt++;
                    }
                }
            }
        }

        for (j=0; j<row_weight; j++)
        {
            e1 = mod2sparse_first_in_row(matrix, i);
            mod2sparse_delete(matrix, e1);
        }
    }

    return cnt;
}


int enumerate_4cyc(mod2sparse *matrix, cyc4_entry *cyc4_rec)
{
    mod2entry *e1, *e2, *e3, *e4;
    int i, j, row_weight;
    int cnt = 0;

    for (i=0; i < matrix->n_rows; i++)
    {
        row_weight = 0;
        for (e1 = mod2sparse_first_in_row(matrix, i); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_row(e1))
        {   
            row_weight++;
            for (e2 = mod2sparse_first_in_col(matrix, e1->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if (e2->row == i)
                    continue;
                for (e3 = mod2sparse_first_in_row(matrix, e2->row); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
                {
                    if (e3->col == e1->col || e3->col > e1->col)
                        continue;
                    for (e4 = mod2sparse_first_in_col(matrix, e3->col); !mod2sparse_at_end(e4); e4 = mod2sparse_next_in_col(e4))
                    {
                        if (e4->row == i)
                        {
                            cyc4_rec[cnt].r1 = i;
                            cyc4_rec[cnt].r2 = e2->row;
                            cyc4_rec[cnt].c1 = e1->col;
                            cyc4_rec[cnt].c2 = e3->col;
                            cnt++;
                        }
                    }
                }
            }
        }

        for (j=0; j<row_weight; j++)
        {
            e1 = mod2sparse_first_in_row(matrix, i);
            mod2sparse_delete(matrix, e1);
        }
    }

    return cnt;
}

void remove_4cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num)
{
    int i, j;
    int diff1, diff2, form_4cyc;
    mod2entry *current_entry, *e;

    for (i=0; i<cyc4_num; i++)
    {
        current_entry = mod2sparse_find(basePCH, cyc4_rec[i].r1, cyc4_rec[i].c1);
        cyc4_rec[i].s1 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc4_rec[i].r1, cyc4_rec[i].c2);
        cyc4_rec[i].s2 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc4_rec[i].r2, cyc4_rec[i].c1);
        cyc4_rec[i].s3 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc4_rec[i].r2, cyc4_rec[i].c2);
        cyc4_rec[i].s4 = &current_entry->shift;
    }

    for (i = 0; i < cyc4_num; i++)
    {
        diff1 = *cyc4_rec[i].s2 - *cyc4_rec[i].s1;
        diff2 = *cyc4_rec[i].s4 - *cyc4_rec[i].s3;
        if (diff1 < 0)
            diff1 += p;
        if (diff2 < 0)
            diff2 += p;
        if (diff1%p == diff2%p) // this 4-cycle will be removed by removing one entry
            form_4cyc = 1;
        else
            form_4cyc = 0;

        while (form_4cyc)
        {
            switch(rand() % 4)
            {
                case 0:
                    if (*cyc4_rec[i].s1 != p)
                    *cyc4_rec[i].s1 = rand() % p;
                    break;
                case 1:
                    if (*cyc4_rec[i].s2 != p)
                    *cyc4_rec[i].s2 = rand() % p;
                    break;
                case 2:
                    if (*cyc4_rec[i].s3 != p)
                    *cyc4_rec[i].s3 = rand() % p;
                    break;
                case 3:
                    if (*cyc4_rec[i].s4 != p)
                    *cyc4_rec[i].s4 = rand() % p;
                    break;
            }

            diff1 = *cyc4_rec[i].s2 - *cyc4_rec[i].s1;
            diff2 = *cyc4_rec[i].s4 - *cyc4_rec[i].s3;
            if (diff1 < 0)
                diff1 += p;
            if (diff2 < 0)
                diff2 += p;
            if (diff1%p == diff2%p) 
                continue;

            for(j=0; j<i; j++)
            {
                diff1 = *cyc4_rec[j].s2 - *cyc4_rec[j].s1;
                diff2 = *cyc4_rec[j].s4 - *cyc4_rec[j].s3;
                if (diff1 < 0)
                    diff1 += p;
                if (diff2 < 0)
                    diff2 += p;
                if (diff1%p == diff2%p)
                    break;
            }
            if (j == i)
                form_4cyc = 0;
        }
    }
}


int count_6cyc(mod2sparse *matrix)
{
    mod2entry *e1, *e2, *e3, *e4, *e5, *e6;
    int i, j, row_weight;
    int cnt = 0;

    for (i=0; i < matrix->n_rows; i++)
    {
        row_weight = 0;
        for (e1 = mod2sparse_first_in_row(matrix, i); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_row(e1))
        {   
            row_weight++;
            for (e2 = mod2sparse_first_in_col(matrix, e1->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if (e2->row == i)
                    continue;
                for (e3 = mod2sparse_first_in_row(matrix, e2->row); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
                {
                    if (e3->col == e1->col)
                        continue;
                    for (e4 = mod2sparse_first_in_col(matrix, e3->col); !mod2sparse_at_end(e4); e4 = mod2sparse_next_in_col(e4))
                    {
                        if (e4->row == e2->row || e4->row == i)
                            continue;
                        for (e5 = mod2sparse_first_in_row(matrix, e4->row); !mod2sparse_at_end(e5); e5 = mod2sparse_next_in_row(e5))
                        {
                            if (e5->col == e3->col || e5->col == e1->col)
                                continue;
                            for (e6 = mod2sparse_first_in_col(matrix, e5->col); !mod2sparse_at_end(e6); e6 = mod2sparse_next_in_col(e6))
                            {
                                if ((e6->row == i) && (e5->col > e1->col))
                                    cnt++;
                            }
                        }
                    }
                }
            }
        }

        for (j=0; j<row_weight; j++)
        {
            e1 = mod2sparse_first_in_row(matrix, i);
            mod2sparse_delete(matrix, e1);
        }
    }
    return cnt;
}

int enumerate_6cyc(mod2sparse *matrix, cyc6_entry *cyc6_rec)
{
    mod2entry *e1, *e2, *e3, *e4, *e5, *e6;
    int i, j, row_weight;
    int cnt = 0;

    for (i=0; i < matrix->n_rows; i++)
    {
        row_weight = 0;
        for (e1 = mod2sparse_first_in_row(matrix, i); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_row(e1))
        {   
            row_weight++;
            for (e2 = mod2sparse_first_in_col(matrix, e1->col); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if (e2->row == i)
                    continue;
                for (e3 = mod2sparse_first_in_row(matrix, e2->row); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
                {
                    if (e3->col == e1->col)
                        continue;
                    for (e4 = mod2sparse_first_in_col(matrix, e3->col); !mod2sparse_at_end(e4); e4 = mod2sparse_next_in_col(e4))
                    {
                        if (e4->row == e2->row || e4->row == i)
                            continue;
                        for (e5 = mod2sparse_first_in_row(matrix, e4->row); !mod2sparse_at_end(e5); e5 = mod2sparse_next_in_row(e5))
                        {
                            if (e5->col == e3->col || e5->col == e1->col)
                                continue;
                            for (e6 = mod2sparse_first_in_col(matrix, e5->col); !mod2sparse_at_end(e6); e6 = mod2sparse_next_in_col(e6))
                            {
                                if ((e6->row == i) && (e5->col > e1->col))
                                {
                                    cyc6_rec[cnt].r1 = i;
                                    cyc6_rec[cnt].r2 = e2->row;
                                    cyc6_rec[cnt].r3 = e4->row;
                                    cyc6_rec[cnt].c1 = e1->col;
                                    cyc6_rec[cnt].c2 = e3->col;
                                    cyc6_rec[cnt].c3 = e5->col;
                                    cnt++;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (j=0; j<row_weight; j++)
        {
            e1 = mod2sparse_first_in_row(matrix, i);
            mod2sparse_delete(matrix, e1);
        }
    }
    return cnt;
}

void remove_6cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num, cyc6_entry *cyc6_rec, int cyc6_num)
{
    int i, j;
    int diff1, diff2, diff3, diffc1, diffc2, diffc3, form_4cyc, form_6cyc;
    int s1, s2, s3, s4, s5, s6;
    mod2entry *current_entry, *e;
    int cnt = 0;
    int cnt1 = 0;

    form_4cyc = 0;

    for (i = 0; i < cyc6_num; i++)
    {
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r1, cyc6_rec[i].c1);
        cyc6_rec[i].s1 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r1, cyc6_rec[i].c3);
        cyc6_rec[i].s2 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r2, __min(cyc6_rec[i].c1, cyc6_rec[i].c2));
        cyc6_rec[i].s3 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r2, __max(cyc6_rec[i].c1, cyc6_rec[i].c2));
        cyc6_rec[i].s4 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r3, __min(cyc6_rec[i].c2, cyc6_rec[i].c3));
        cyc6_rec[i].s5 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r3, __max(cyc6_rec[i].c2, cyc6_rec[i].c3));
        cyc6_rec[i].s6 = &current_entry->shift;

        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r1, cyc6_rec[i].c1);
        cyc6_rec[i].ss1 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r2, cyc6_rec[i].c1);
        cyc6_rec[i].ss2 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, __min(cyc6_rec[i].r2, cyc6_rec[i].r3), cyc6_rec[i].c2);
        cyc6_rec[i].ss3 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, __max(cyc6_rec[i].r2, cyc6_rec[i].r3), cyc6_rec[i].c2);
        cyc6_rec[i].ss4 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, __min(cyc6_rec[i].r1, cyc6_rec[i].r3), cyc6_rec[i].c3);
        cyc6_rec[i].ss5 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, __max(cyc6_rec[i].r1, cyc6_rec[i].r3), cyc6_rec[i].c3);
        cyc6_rec[i].ss6 = &current_entry->shift;
    }

    for (i=0; i<cyc6_num; i++)
    {
        diff1 = *cyc6_rec[i].s2 - *cyc6_rec[i].s1;
        diff2 = *cyc6_rec[i].s4 - *cyc6_rec[i].s3;
        diff3 = *cyc6_rec[i].s6 - *cyc6_rec[i].s5;

        if (diff1 < 0)
            diff1 += p;
        if (diff2 < 0)
            diff2 += p;
        if (diff3 < 0)
            diff3 += p;

        diffc1 = *cyc6_rec[i].ss2 - *cyc6_rec[i].ss1;
        diffc2 = *cyc6_rec[i].ss4 - *cyc6_rec[i].ss3;
        diffc3 = *cyc6_rec[i].ss6 - *cyc6_rec[i].ss5;

        if (diffc1 < 0)
            diffc1 += p;
        if (diffc2 < 0)
            diffc2 += p;
        if (diffc3 < 0)
            diffc3 += p;

        if ( (diff1%p == diff2%p && diff1%p == diff3%p && diff2%p == diff3%p)
            || (diffc1%p == diffc2%p && diffc1%p == diffc3%p && diffc2%p == diffc3%p) 
            || (diffc1%p + diffc2%p == diffc3%p )
            || (diffc1%p + diffc3%p == diffc2%p) 
            || (diffc2%p + diffc3%p == diffc1%p)
            || (abs(diffc1%p - diffc2%p) == diffc3%p)    
            || (abs(diffc1%p - diffc3%p) == diffc2%p)    
            || (abs(diffc2%p - diffc3%p) == diffc1%p)
            || (diff1%p + diff2%p == diff3%p)
            || (diff1%p + diff3%p == diff2%p)
            || (diff2%p + diff3%p == diff1%p)
            || (abs(diff1%p - diff2%p) == diff3%p)
            || (abs(diff1%p - diff3%p) == diff2%p)
            || (abs(diff2%p - diff3%p) == diff1%p)
            || (diffc1%p + diffc2%p + diffc3%p == p)
            || (diff1%p + diff2%p + diff3%p == p)
            || (cyc6_rec[i].s1 == cyc6_rec[i].s2)
            || (cyc6_rec[i].s3 == cyc6_rec[i].s4)
            || (cyc6_rec[i].s5 == cyc6_rec[i].s6)
            || ((*cyc6_rec[i].s1 + *cyc6_rec[i].s2 + *cyc6_rec[i].s3 + *cyc6_rec[i].s4 + *cyc6_rec[i].s5 + *cyc6_rec[i].s6) % p == 0)
            || (cyc6_rec[i].s1 == cyc6_rec[i].s6)
            || (cyc6_rec[i].s2 == cyc6_rec[i].s3)
            || (cyc6_rec[i].s4 == cyc6_rec[i].s5)
        )
            form_6cyc = 1;
        else
            form_6cyc = 0;

        s1 = *cyc6_rec[i].s1;
        s2 = *cyc6_rec[i].s2;
        s3 = *cyc6_rec[i].s3;
        s4 = *cyc6_rec[i].s4;
        s5 = *cyc6_rec[i].s5;
        s6 = *cyc6_rec[i].s6;

        cnt = 0;
        while (form_6cyc || form_4cyc)
        {
            switch(rand() % 6)
            {
                case 0:
                    if (*cyc6_rec[i].s1 != p)
                    *cyc6_rec[i].s1 = rand() % p;
                    break;
                case 1:
                    if (*cyc6_rec[i].s2 != p)
                    *cyc6_rec[i].s2 = rand() % p;
                    break;
                case 2:
                    if (*cyc6_rec[i].s3 != p)
                    *cyc6_rec[i].s3 = rand() % p;
                    break;
                case 3:
                    if (*cyc6_rec[i].s4 != p)
                    *cyc6_rec[i].s4 = rand() % p;
                    break;
                case 4:
                    if (*cyc6_rec[i].s5 != p)
                    *cyc6_rec[i].s5 = rand() % p;
                    break;
                case 5:
                    if (*cyc6_rec[i].s6 != p)
                    *cyc6_rec[i].s6 = rand() % p;
                    break;
            }
            
            diff1 = *cyc6_rec[i].s2 - *cyc6_rec[i].s1;
            diff2 = *cyc6_rec[i].s4 - *cyc6_rec[i].s3;
            diff3 = *cyc6_rec[i].s6 - *cyc6_rec[i].s5;

            if (diff1 < 0)
                diff1 += p;
            if (diff2 < 0)
                diff2 += p;
            if (diff3 < 0)
                diff3 += p;

            diffc1 = *cyc6_rec[i].ss2 - *cyc6_rec[i].ss1;
            diffc2 = *cyc6_rec[i].ss4 - *cyc6_rec[i].ss3;
            diffc3 = *cyc6_rec[i].ss6 - *cyc6_rec[i].ss5;

            if (diffc1 < 0)
                diffc1 += p;
            if (diffc2 < 0)
                diffc2 += p;
            if (diffc3 < 0)
                diffc3 += p;

            if ( (diff1%p == diff2%p && diff1%p == diff3%p && diff2%p == diff3%p)
                || (diffc1%p == diffc2%p && diffc1%p == diffc3%p && diffc2%p == diffc3%p) 
                || (diffc1%p + diffc2%p == diffc3%p )
                || (diffc1%p + diffc3%p == diffc2%p) 
                || (diffc2%p + diffc3%p == diffc1%p)
                || (abs(diffc1%p - diffc2%p) == diffc3%p)    
                || (abs(diffc1%p - diffc3%p) == diffc2%p)    
                || (abs(diffc2%p - diffc3%p) == diffc1%p)
                || (diff1%p + diff2%p == diff3%p)
                || (diff1%p + diff3%p == diff2%p)
                || (diff2%p + diff3%p == diff1%p)
                || (abs(diff1%p - diff2%p) == diff3%p)
                || (abs(diff1%p - diff3%p) == diff2%p)
                || (abs(diff2%p - diff3%p) == diff1%p)
                || (diffc1%p + diffc2%p + diffc3%p == p)
                || (diff1%p + diff2%p + diff3%p == p)
                || (cyc6_rec[i].s1 == cyc6_rec[i].s2)
                || (cyc6_rec[i].s3 == cyc6_rec[i].s4)
                || (cyc6_rec[i].s5 == cyc6_rec[i].s6)
                || ((*cyc6_rec[i].s1 + *cyc6_rec[i].s2 + *cyc6_rec[i].s3 + *cyc6_rec[i].s4 + *cyc6_rec[i].s5 + *cyc6_rec[i].s6) % p == 0)
                || (cyc6_rec[i].s1 == cyc6_rec[i].s6)
                || (cyc6_rec[i].s2 == cyc6_rec[i].s3)
                || (cyc6_rec[i].s4 == cyc6_rec[i].s5)
            )
                continue;

            for(j=0; j<i; j++)
            {
                diff1 = *cyc6_rec[i].s2 - *cyc6_rec[i].s1;
                diff2 = *cyc6_rec[i].s4 - *cyc6_rec[i].s3;
                diff3 = *cyc6_rec[i].s6 - *cyc6_rec[i].s5;
                if (diff1 < 0)
                    diff1 += p;
                if (diff2 < 0)
                    diff2 += p;
                if (diff3 < 0)
                    diff3 += p;

                diffc1 = *cyc6_rec[i].ss2 - *cyc6_rec[i].ss1;
                diffc2 = *cyc6_rec[i].ss4 - *cyc6_rec[i].ss3;
                diffc3 = *cyc6_rec[i].ss6 - *cyc6_rec[i].ss5;

                if (diffc1 < 0)
                    diffc1 += p;
                if (diffc2 < 0)
                    diffc2 += p;
                if (diffc3 < 0)
                    diffc3 += p;

                if ( (diff1%p == diff2%p && diff1%p == diff3%p && diff2%p == diff3%p)
                    || (diffc1%p == diffc2%p && diffc1%p == diffc3%p && diffc2%p == diffc3%p) 
                    || (diffc1%p + diffc2%p == diffc3%p )
                    || (diffc1%p + diffc3%p == diffc2%p) 
                    || (diffc2%p + diffc3%p == diffc1%p)
                    || (abs(diffc1%p - diffc2%p) == diffc3%p)    
                    || (abs(diffc1%p - diffc3%p) == diffc2%p)    
                    || (abs(diffc2%p - diffc3%p) == diffc1%p)
                    || (diff1%p + diff2%p == diff3%p)
                    || (diff1%p + diff3%p == diff2%p)
                    || (diff2%p + diff3%p == diff1%p)
                    || (abs(diff1%p - diff2%p) == diff3%p)
                    || (abs(diff1%p - diff3%p) == diff2%p)
                    || (abs(diff2%p - diff3%p) == diff1%p)
                    || (diffc1%p + diffc2%p + diffc3%p == p)
                    || (diff1%p + diff2%p + diff3%p == p)
                    || (cyc6_rec[i].s1 == cyc6_rec[i].s2)
                    || (cyc6_rec[i].s3 == cyc6_rec[i].s4)
                    || (cyc6_rec[i].s5 == cyc6_rec[i].s6)
                    || ((*cyc6_rec[i].s1 + *cyc6_rec[i].s2 + *cyc6_rec[i].s3 + *cyc6_rec[i].s4 + *cyc6_rec[i].s5 + *cyc6_rec[i].s6) % p == 0)
                    || (cyc6_rec[i].s1 == cyc6_rec[i].s6)
                    || (cyc6_rec[i].s2 == cyc6_rec[i].s3)
                    || (cyc6_rec[i].s4 == cyc6_rec[i].s5)
                )
                    break;
            }
            if (j == i)
                form_6cyc = 0;

            for (j=0; j<cyc4_num; j++)
            {
                diff1 = *cyc4_rec[j].s2 - *cyc4_rec[j].s1;
                diff2 = *cyc4_rec[j].s4 - *cyc4_rec[j].s3;

                if (diff1 < 0)
                    diff1 += p;
                if (diff2 < 0)
                    diff2 += p;
                if (diff1%p == diff2%p)
                {
                    form_4cyc = 1;
                    break;
                }
            }
            if (j == cyc4_num)
                form_4cyc = 0;

            cnt++;
            if (cnt > 20)
            {
                cnt1++;
                *cyc6_rec[i].s1 = s1;
                *cyc6_rec[i].s2 = s2;
                *cyc6_rec[i].s3 = s3;
                *cyc6_rec[i].s4 = s4;
                *cyc6_rec[i].s5 = s5;
                *cyc6_rec[i].s6 = s6;
                break;
            }
        }
    }
}

int reduce_6cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num)
{
    int i, j;
    int diff1, diff2, form_4cyc;
    mod2sparse *exPCH_cpy;
    mod2entry *e;
    int cyc6_num = 0x0FFFFFFF;
    int ones_num = 0;
    int try_num;
    int *shift_cpy, **shift;
    int cnt;

    // get the number of 1's in base matrix
    for (i=0; i<basePCH->n_rows; i++)
        for (e = mod2sparse_first_in_row(basePCH, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
            ones_num++;

    // back up shift value
    shift_cpy = chk_alloc(ones_num, sizeof(int));
    shift = chk_alloc(ones_num, sizeof(int *));
    j = 0;
    for (i=0; i<basePCH->n_rows; i++)
        for (e = mod2sparse_first_in_row(basePCH, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
        {
            shift_cpy[j] = e->shift;
            shift[j] = &e->shift;
            j++;
        }

    try_num = 0;
    while (try_num < 1000)
    {
        for (i=0; i<ones_num; i++)
            shift_cpy[i] = *shift[i];
        for (i=0; i<10; i++)
            *shift[rand() % ones_num] = rand() % p;

        remove_4cyc(basePCH, exPCH, p, cyc4_rec, cyc4_num);

        exPCH_cpy = mod2sparse_allocate(exPCH->n_rows, exPCH->n_cols);
        for (i=0; i<basePCH->n_rows; i++)
            for (e = mod2sparse_first_in_row(basePCH, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
                insert_sub_matrix(exPCH_cpy, p * e->row, p * e->col, p, e->shift);

        cnt = count_6cyc(exPCH_cpy);
        if (cnt < cyc6_num)
        {
            try_num = 0;
            cyc6_num = cnt;
            printf("6-cycle number = %d\n", cyc6_num);
        }
        else
        {
            try_num++;
            for (i=0; i<ones_num; i++)
                *shift[i] = shift_cpy[i];
        }

        mod2sparse_free(exPCH_cpy);
        if (cyc6_num == 0)
            break;
    }

    return cyc6_num;
}


void new_remove_6cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num, cyc6_entry *cyc6_rec, int cyc6_num)
{
    int i, j;
    int diff1, diff2, form_4cyc, form_6cyc;
    int s1, s2, s3, s4, s5, s6;
    mod2entry *current_entry;
    int cnt = 0;
    int cnt1 = 0;

    tanner_graph *g;
    mod2sparse *exPCH_cpy;

    // Initialize tanner graph for BFS
    g = chk_alloc(1, sizeof(tanner_graph));
    g->bit_node = chk_alloc( exPCH->n_cols, sizeof(node_entry));
    g->check_node = chk_alloc( exPCH->n_rows, sizeof(node_entry));

    form_4cyc = 0;
    for (i = 0; i < cyc6_num; i++)
    {
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r1, cyc6_rec[i].c1);
        cyc6_rec[i].s1 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r1, cyc6_rec[i].c3);
        cyc6_rec[i].s2 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r2, __min(cyc6_rec[i].c1, cyc6_rec[i].c2));
        cyc6_rec[i].s3 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r2, __max(cyc6_rec[i].c1, cyc6_rec[i].c2));
        cyc6_rec[i].s4 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r3, __min(cyc6_rec[i].c2, cyc6_rec[i].c3));
        cyc6_rec[i].s5 = &current_entry->shift;
        current_entry = mod2sparse_find(basePCH, cyc6_rec[i].r3, __max(cyc6_rec[i].c2, cyc6_rec[i].c3));
        cyc6_rec[i].s6 = &current_entry->shift;
    }

    for (i=0; i<cyc6_num; i++)
    {
        printf("\n %d", i);

        exPCH_cpy = mod2sparse_allocate(exPCH->n_rows, exPCH->n_cols);
        insert_sub_matrix(exPCH_cpy, p * cyc6_rec[i].r1, p * cyc6_rec[i].c1, p, *cyc6_rec[i].s1);
        insert_sub_matrix(exPCH_cpy, p * cyc6_rec[i].r1, p * cyc6_rec[i].c3, p, *cyc6_rec[i].s2);
        insert_sub_matrix(exPCH_cpy, p * cyc6_rec[i].r2, p * __min(cyc6_rec[i].c1, cyc6_rec[i].c2), p, *cyc6_rec[i].s3);
        insert_sub_matrix(exPCH_cpy, p * cyc6_rec[i].r2, p * __max(cyc6_rec[i].c1, cyc6_rec[i].c2), p, *cyc6_rec[i].s4);
        insert_sub_matrix(exPCH_cpy, p * cyc6_rec[i].r3, p * __min(cyc6_rec[i].c2, cyc6_rec[i].c3), p, *cyc6_rec[i].s5);
        insert_sub_matrix(exPCH_cpy, p * cyc6_rec[i].r3, p * __max(cyc6_rec[i].c2, cyc6_rec[i].c3), p, *cyc6_rec[i].s6);
        g->matrix = exPCH_cpy;

        if (check_cycle(g, p * cyc6_rec[i].r1, 6))
        {
            form_6cyc = 1;
            cnt1++;
        }
        else
            form_6cyc = 0;
        mod2sparse_free(exPCH_cpy);
    }

    printf("\n There are %d 6-cycles in base matrix to form %d 6cycles in expansion matrix\n", cnt1, cnt1*64);

    free(g->bit_node);
    free(g->check_node);
    free(g);
}


int count_weight4(mod2sparse *mat)
{
    int i1, i2, i3, i4, j;
    char *u; // Pointer to unpacked vector to multiply, length is the number of cols of mat
    char *v; // Pointer to unpacked vector to store the result, length is the number of rows of mat
    int cnt = 0;
    int zero_flag;

    u = chk_alloc(mat->n_cols, sizeof(char));
    v = chk_alloc(mat->n_rows, sizeof(char));

    for (i1 = 0; i1 < mat->n_cols; i1++)
        for (i2 = i1 + 1; i2 < mat->n_cols; i2++)
            for (i3 = i2 + 1; i3 < mat->n_cols; i3++)
                for (i4 = i3 + 1; i4 < mat->n_cols; i4++)
                {
                    for (j = 0; j < mat->n_cols; j++)
                        u[j] = 0;
                    u[i1] = 1;
                    u[i2] = 1;
                    u[i3] = 1;
                    u[i4] = 1;
                    mod2sparse_mulvec(mat, u, v);
                    zero_flag = 0;
                    for (j = 0; j < mat->n_rows; j++)
                    {
                        if(v[j]==1)
                            zero_flag = 1;
                    }
                    if(!zero_flag)
                        cnt++;
                }

    free(u);
    free(v);
    return cnt;
}

int enumerate_weight4(mod2sparse *mat, weight4_entry *weight4_rec)
{
    int i1, i2, i3, i4, j;
    char *u; // Pointer to unpacked vector to multiply, length is the number of cols of mat
    char *v; // Pointer to unpacked vector to store the result, length is the number of rows of mat
    int cnt = 0;
    int zero_flag;

    u = chk_alloc(mat->n_cols, sizeof(char));
    v = chk_alloc(mat->n_rows, sizeof(char));

    for (i1 = 0; i1 < mat->n_cols; i1++)
        for (i2 = i1 + 1; i2 < mat->n_cols; i2++)
            for (i3 = i2 + 1; i3 < mat->n_cols; i3++)
                for (i4 = i3 + 1; i4 < mat->n_cols; i4++)
                {
                    for (j = 0; j < mat->n_cols; j++)
                        u[j] = 0;
                    u[i1] = 1;
                    u[i2] = 1;
                    u[i3] = 1;
                    u[i4] = 1;
                    mod2sparse_mulvec(mat, u, v);
                    zero_flag = 0;
                    for (j = 0; j < mat->n_rows; j++)
                    {
                        if(v[j]==1)
                            zero_flag = 1;
                    }
                    if(!zero_flag)
                    {
                        weight4_rec[cnt].c[0] = i1;
                        weight4_rec[cnt].c[1] = i2;
                        weight4_rec[cnt].c[2] = i3;
                        weight4_rec[cnt].c[3] = i4;
                        cnt++;
                    }
                }

    free(u);
    free(v);
    return cnt;
}

int chk_expand_weight4(mod2sparse *mat, int p, weight4_entry w4_rec)
{
    mod2sparse *baseMat, *exMat;
    mod2entry *e, *e1;
    int row_weight;

    int i1, i2, i3, i4, c1, c2, c3, c4, j1, j2, k;
    char *u; // Pointer to unpacked vector to multiply, length is the number of cols of mat
    char *v; // Pointer to unpacked vector to store the result, length is the number of rows of mat
    int zero_flag = 1;
    int cnt;
    
    u = chk_alloc(4*p, sizeof(char));
    v = chk_alloc(mat->n_rows*p, sizeof(char));
    baseMat = mod2sparse_allocate(mat->n_rows, 4);
    exMat = mod2sparse_allocate(mat->n_rows*p, 4*p);

    for (k = 0; k < 4; k++)
    {
        for (e = mod2sparse_first_in_col(mat, w4_rec.c[k]); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
        {
            insert_sub_matrix(exMat, e->row * p, k*p, p, e->shift);
            mod2sparse_insert(baseMat, e->row, k);
            e1 = mod2sparse_find(baseMat, e->row, k);
            e1->shift = e->shift;
        }
    }

    for (e = mod2sparse_first_in_col(baseMat, 0); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
    {
        row_weight = mod2sparse_count_row(baseMat, e->row);
        if (row_weight == 2)
            break;
    }

    if (row_weight == 2)
    {
        i1 = e->shift;
        e = mod2sparse_next_in_row(e);

        i2 = -1;
        i3 = -1;
        i4 = -1;

        switch (e->col)
        {
            case 1:
                if (e->shift == p)
                    i2 = 0;
                else
                    i2 = e->shift;
                break;
            case 2:
                if (e->shift == p)
                    i3 = 0;
                else
                    i3 = e->shift;
                break;
            case 3:
                if (e->shift == p)
                    i4 = 0;
                else
                    i4 = e->shift;
                break;
            default:
                fprintf(stderr, "Not in valid column!\n");
                exit(1);
        }

        if (e->shift < 0 || e->shift > p)
        {
            fprintf(stderr, "Shift %d out of bound!\n", e->shift);
            exit(1);
        }

        for (j1 = 0; j1 < p; j1++)
        {
            for (j2 = 0; j2 < p; j2++)
            {
                for (k = 0; k < 4*p; k++)
                    u[k] = 0;
                u[i1] = 1;

                if (i2 != -1)
                {
                    u[p+i2] = 1;
                    u[2*p+j1] = 1;
                    u[3*p+j2] = 1;
                }
                else if (i3 != -1)
                {
                    u[p+j1] = 1;
                    u[2*p+i3] = 1;
                    u[3*p+j2] = 1;
                }
                else if (i4 != -1)
                {
                    u[p+j1] = 1;
                    u[2*p+j2] = 1;
                    u[3*p+i4] = 1;
                }

                zero_flag = 0;
                mod2sparse_mulvec(exMat, u, v);
                for (k = 0; k < exMat->n_rows; k++)
                {
                    if(v[k]==1)
                    {
                        zero_flag = 1;
                        break;
                    }
                }
                if(!zero_flag)
                {
                    printf("Find a weight-4 codeword in expansion matrix!\n");
                    free(u);
                    free(v);
                    mod2sparse_free(baseMat);
                    mod2sparse_free(exMat);
                    return 1;
                }
            }
        }
    }
    else if (row_weight == 4)
    {
        i1 = 0;
        for (i2 = 0; i2 < p; i2++)
        {
            for (i3 = 0; i3 < p; i3++)
            {
                for (i4 = 0; i4 < p; i4++)
                {
                    for (k = 0; k < 4*p; k++)
                        u[k] = 0;
                    u[i1] = 1;
                    u[p+i2] = 1;
                    u[2*p+i3] = 1;
                    u[3*p+i4] = 1;

                    mod2sparse_mulvec(exMat, u, v);
                    zero_flag = 0;
                    for (k = 0; k < mat->n_rows * p; k++)
                    {
                        if(v[k]==1)
                            zero_flag = 1;
                    }
                    if(!zero_flag)
                    {
                        printf("Find a weight-4 codeword in expansion matrix!\n");
                        free(u);
                        free(v);
                        mod2sparse_free(baseMat);
                        mod2sparse_free(exMat);
                        return 1;
                    }
                }
            }
        }
    }
    else
    {
        fprintf(stderr, "Row weight %d not supported!\n", row_weight);
        exit(1);
    }

    free(u);
    free(v);
    mod2sparse_free(baseMat);
    mod2sparse_free(exMat);
    return 0;
}


void expand_matrix(mod2sparse *exMatrix, mod2sparse *baseMatrix, int exFactor)
{
    mod2entry *e;
    int baseRows, baseCols;
    int rowPos, colPos;
    int i, j;

    baseRows = mod2sparse_rows(baseMatrix);
    baseCols = mod2sparse_cols(baseMatrix);

    for (i = 0; i < baseRows; i++)
    {
        j = 0;
        for (e = mod2sparse_first_in_row(baseMatrix, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
        {
            rowPos = e->row * exFactor;
            colPos = e->col * exFactor;
            insert_sub_matrix(exMatrix, rowPos, colPos, exFactor, e->shift);
        }
    }
}

int cyc4_cal(mod2sparse *basePCH, int qc_size)
{
    int num4 = 0;
    int row_num = mod2sparse_rows(basePCH);
    int col_num = mod2sparse_cols(basePCH);
    mod2entry *e1, *e2, *e3, *e4;

    for (int k0 = 0; k0 < col_num; k0++)
        for (e1 = mod2sparse_first_in_col(basePCH, k0); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_col(e1))
            for (e2 = mod2sparse_next_in_col(e1); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                e3 = mod2sparse_next_in_row(e1);
                e4 = mod2sparse_next_in_row(e2);

                while(!mod2sparse_at_end(e3) && !mod2sparse_at_end(e4))
                {
                    if (e3->col > e4->col)
                        e4 = mod2sparse_next_in_row(e4);
                    else if (e3->col < e4->col)
                        e3 = mod2sparse_next_in_row(e3);
                    else
                    {
                        if ((e1->shift - e2->shift + e4->shift - e3->shift ) % qc_size == 0)
                            num4++;
                        e3 = mod2sparse_next_in_row(e3);
                        e4 = mod2sparse_next_in_row(e4);
                    }
                }
            }


    return num4;
}

int cyc4_cal_mask(mod2sparse *basePCH, int **mask_matrix, int qc_size, int pad_bit)
{
    int num4 = 0;
    int row_num = mod2sparse_rows(basePCH);
    int col_num = mod2sparse_cols(basePCH);
    mod2entry *e1, *e2, *e3, *e4;

    for (int k0 = 0; k0 < col_num; k0++)
        for (e1 = mod2sparse_first_in_col(basePCH, k0); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_col(e1))
        {
            if ((pad_bit==qc_size) && (mask_matrix[e1->row][e1->col]==1))
                continue;

            for (e2 = mod2sparse_next_in_col(e1); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if ((pad_bit==qc_size) && (mask_matrix[e2->row][e2->col]==1))
                    continue;
                if ((mask_matrix[e1->row][e1->col]==1 && mask_matrix[e2->row][e2->col]==2)
                    || (mask_matrix[e1->row][e1->col]==2 && mask_matrix[e2->row][e2->col]==1))
                    continue;

                e3 = mod2sparse_next_in_row(e1);
                e4 = mod2sparse_next_in_row(e2);

                while(!mod2sparse_at_end(e3) && !mod2sparse_at_end(e4))
                {
                    if (e3->col > e4->col)
                        e4 = mod2sparse_next_in_row(e4);
                    else if (e3->col < e4->col)
                        e3 = mod2sparse_next_in_row(e3);
                    else
                    {
                        if ((pad_bit==qc_size) && (mask_matrix[e3->row][e3->col]==1) || mask_matrix[e4->row][e4->col]==1){}
                        else if ((mask_matrix[e3->row][e3->col]==1 && mask_matrix[e4->row][e4->col]==2)
                            || (mask_matrix[e3->row][e3->col]==2 && mask_matrix[e4->row][e4->col]==1)){}
                        else if ((e1->shift - e2->shift + e4->shift - e3->shift ) % qc_size == 0)
                            num4++;
                        e3 = mod2sparse_next_in_row(e3);
                        e4 = mod2sparse_next_in_row(e4);
                    }
                }
            }
        }

    return num4;
}

int cyc4_cal_ibex(mod2sparse *basePCH, int **fade_matrix, int qc_size, int pad_bit)
{
    int num4 = 0;
    int row_num = mod2sparse_rows(basePCH);
    int col_num = mod2sparse_cols(basePCH);
    mod2entry *e1, *e2, *e3, *e4;

    for (int k0 = 0; k0 < col_num; k0++)
    {
        for (e1 = mod2sparse_first_in_col(basePCH, k0); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_col(e1))
        {
            if ((pad_bit==qc_size) && (fade_matrix[e1->row][e1->col]==1))
                continue;

            for (e2 = mod2sparse_next_in_col(e1); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if ((pad_bit==qc_size) && (fade_matrix[e2->row][e2->col]==1))
                    continue;
                if (fade_matrix[e1->row][e1->col]==1 && e2->row == row_num-1)
                    continue;

                e3 = mod2sparse_next_in_row(e1);
                e4 = mod2sparse_next_in_row(e2);

                while(!mod2sparse_at_end(e3) && !mod2sparse_at_end(e4))
                {
                    if (e3->col > e4->col)
                        e4 = mod2sparse_next_in_row(e4);
                    else if (e3->col < e4->col)
                        e3 = mod2sparse_next_in_row(e3);
                    else
                    {
                        if ((pad_bit==qc_size) && (fade_matrix[e3->row][e3->col]==1 || fade_matrix[e4->row][e4->col]==1)) {}
                        else if (fade_matrix[e3->row][e3->col]==1 && e4->row==row_num-1) {}
                        else if ((e1->shift - e2->shift + e4->shift - e3->shift ) % qc_size == 0)
                        {
                            num4++;
                            printf("col%3d | col%3d \n", e1->col, e3->col);
                            printf("row%2d: %3d | %3d\n", e1->row, e1->shift, e3->shift);
                            printf("row%2d: %3d | %3d\n", e2->row, e2->shift, e4->shift);
                            printf("**************************\n");
                        }
                        e3 = mod2sparse_next_in_row(e3);
                        e4 = mod2sparse_next_in_row(e4);
                    }
                }
            }
        }
    }


    return num4;
}

int cyc6_cal(mod2sparse *basePCH, int qc_size)
{
    int num6 = 0;
    int row_num = mod2sparse_rows(basePCH);
    int col_num = mod2sparse_cols(basePCH);
    mod2entry *e1, *e2, *e3, *e4, *e5, *e6;

    for (int k0 = 0; k0 < col_num; k0++)
        for (e1 = mod2sparse_first_in_col(basePCH, k0); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_col(e1))
            for (e2 = mod2sparse_next_in_col(e1); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
                for (e6 = mod2sparse_next_in_row(e1); !mod2sparse_at_end(e6); e6 = mod2sparse_next_in_row(e6))
                    for (e3 = mod2sparse_next_in_row(e2); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
                    {
                        if(e6->col == e3->col)
                            continue;

                        e4 = mod2sparse_first_in_col(basePCH, e3->col);
                        e5 = mod2sparse_first_in_col(basePCH, e6->col);

                        while(!mod2sparse_at_end(e4) && !mod2sparse_at_end(e5))
                        {
                            if (e4->row < e5->row)
                                e4 = mod2sparse_next_in_col(e4);
                            else if (e5->row < e4->row)
                                e5 = mod2sparse_next_in_col(e5);
                            else
                            {
                                if (((e1->shift - e2->shift + e3->shift - e4->shift + e5->shift - e6->shift ) % qc_size == 0)
                                    && e4->row!=e1->row && e4->row!=e2->row)
                                    num6++;
                                e4 = mod2sparse_next_in_col(e4);
                                e5 = mod2sparse_next_in_col(e5);
                            }
                        }
                    }

    return num6;
}

int cyc6_cal_mask(mod2sparse *basePCH, int **mask_matrix, int qc_size, int pad_bit)
{
    int num6 = 0;
    int row_num = mod2sparse_rows(basePCH);
    int col_num = mod2sparse_cols(basePCH);
    mod2entry *e1, *e2, *e3, *e4, *e5, *e6;

    for (int k0 = 0; k0 < col_num; k0++)
    {
        for (e1 = mod2sparse_first_in_col(basePCH, k0); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_col(e1))
        {
            if ((pad_bit==qc_size) && (mask_matrix[e1->row][e1->col]==1))
                continue;

            for (e2 = mod2sparse_next_in_col(e1); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if ((pad_bit==qc_size) && (mask_matrix[e2->row][e2->col]==1))
                    continue;
                if ((mask_matrix[e1->row][e1->col]==1 && mask_matrix[e2->row][e2->col]==2)
                    || mask_matrix[e1->row][e1->col]==2 && mask_matrix[e2->row][e2->col]==1)
                    continue;

                for (e6 = mod2sparse_next_in_row(e1); !mod2sparse_at_end(e6); e6 = mod2sparse_next_in_row(e6))
                {
                    if ((pad_bit==qc_size) && (mask_matrix[e6->row][e6->col]==1))
                        continue;

                    for (e3 = mod2sparse_next_in_row(e2); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
                    {
                        if ((pad_bit==qc_size) && (mask_matrix[e3->row][e3->col]==1))
                            continue;
                        if (e6->col == e3->col)
                            continue;

                        e4 = mod2sparse_first_in_col(basePCH, e3->col);
                        e5 = mod2sparse_first_in_col(basePCH, e6->col);
                        while (!mod2sparse_at_end(e4) && !mod2sparse_at_end(e5))
                        {
                            if (e4->row < e5->row)
                                e4 = mod2sparse_next_in_col(e4);
                            else if (e5->row < e4->row)
                                e5 = mod2sparse_next_in_col(e5);
                            else
                            {
                                if ((pad_bit==qc_size) && (mask_matrix[e4->row][e4->col]==1) || (mask_matrix[e5->row][e5->col]==1)){}
                                else if ((mask_matrix[e3->row][e3->col]==1) && (mask_matrix[e4->row][e4->col]==2)
                                    || (mask_matrix[e3->row][e3->col]==2) && (mask_matrix[e4->row][e4->col]==1)){}
                                else if ((mask_matrix[e5->row][e5->col]==1) && (mask_matrix[e6->row][e6->col]==2)
                                    || (mask_matrix[e5->row][e5->col]==2) && (mask_matrix[e6->row][e6->col]==1)){}
                                else if (((e1->shift - e2->shift + e3->shift - e4->shift + e5->shift - e6->shift ) % qc_size == 0)
                                        && e4->row!=e1->row && e4->row!=e2->row)
                                    num6++;
                                
                                e4 = mod2sparse_next_in_col(e4);
                                e5 = mod2sparse_next_in_col(e5);
                            }
                        }
                    }
                }
            }
        }
    }

    return num6;
}

int cyc6_cal_ibex(mod2sparse *basePCH, int **fade_matrix, int qc_size, int pad_bit)
{
    int num6 = 0;
    int row_num = mod2sparse_rows(basePCH);
    int col_num = mod2sparse_cols(basePCH);
    mod2entry *e1, *e2, *e3, *e4, *e5, *e6;

    for (int k0 = 0; k0 < col_num; k0++)
    {
        for (e1 = mod2sparse_first_in_col(basePCH, k0); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_col(e1))
        {
            if ((pad_bit==qc_size) && (fade_matrix[e1->row][e1->col]==1))
                continue;

            for (e2 = mod2sparse_next_in_col(e1); !mod2sparse_at_end(e2); e2 = mod2sparse_next_in_col(e2))
            {
                if ((pad_bit==qc_size) && (fade_matrix[e2->row][e2->col]==1))
                    continue;
                if (((fade_matrix[e1->row][e1->col]==1) && (e2->row == row_num - 1))
                    || ((e1->row == row_num - 1) && (fade_matrix[e2->row][e2->col]==1)))
                    continue;

                for (e6 = mod2sparse_next_in_row(e1); !mod2sparse_at_end(e6); e6 = mod2sparse_next_in_row(e6))
                {
                    if ((pad_bit==qc_size) && (fade_matrix[e6->row][e6->col]==1))
                        continue;

                    for (e3 = mod2sparse_next_in_row(e2); !mod2sparse_at_end(e3); e3 = mod2sparse_next_in_row(e3))
                    {
                        if ((pad_bit==qc_size) && (fade_matrix[e3->row][e3->col]==1))
                            continue;
                        if (e6->col == e3->col)
                            continue;

                        e4 = mod2sparse_first_in_col(basePCH, e3->col);
                        e5 = mod2sparse_first_in_col(basePCH, e6->col);
                        while (!mod2sparse_at_end(e4) && !mod2sparse_at_end(e5))
                        {
                            if (e4->row < e5->row)
                                e4 = mod2sparse_next_in_col(e4);
                            else if (e5->row < e4->row)
                                e5 = mod2sparse_next_in_col(e5);
                            else
                            {
                                if ((pad_bit==qc_size) && (fade_matrix[e4->row][e4->col]==1) || (fade_matrix[e5->row][e5->col]==1)){}
                                else if ((fade_matrix[e3->row][e3->col]==1 && e4->row == row_num - 1)
                                    || e3->row == row_num - 1 && fade_matrix[e4->row][e4->col]==1) {}
                                else if ((fade_matrix[e5->row][e5->col]==1 && e6->row == row_num - 1)
                                    || e5->row == row_num - 1 && fade_matrix[e6->row][e6->col]==1) {}
                                else if (((e1->shift - e2->shift + e3->shift - e4->shift + e5->shift - e6->shift ) % qc_size == 0)
                                        && e4->row!=e1->row && e4->row!=e2->row)
                                    num6++;
                                
                                e4 = mod2sparse_next_in_col(e4);
                                e5 = mod2sparse_next_in_col(e5);
                            }
                        }
                    }
                }
            }
        }
    }

    return num6;
}