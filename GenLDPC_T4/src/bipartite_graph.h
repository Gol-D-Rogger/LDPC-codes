#ifndef _BIPARTITE_GRAPH_H
#define _BIPARTITE_GRAPH_H

#include "mod2sparse.h"

#define max_depth 100 // maximum girth = 2*max_depth
#define max_girth 200

#define MAX_CYCLE_LENGTH 8

#define __max(a,b) (((a) > (b)) ? (a) : (b))
#define __min(a,b) (((a) < (b)) ? (a) : (b))

inline double in_pow(double a, int b)
{
    double tmp = 1;
    for (int i = 0; i < b; i++)
        tmp *= a;
    return tmp;
}

typedef struct node_entry
{
    int index; 
    int depth;
    int node_degree;
    struct node_entry *parent;
    enum {white, gray, black, isolated} color;
} node_entry;

typedef struct tanner_graph
{
    int n_nodes;
    node_entry *check_node; // array of check nodes
    node_entry *bit_node;   // array of variable nodes
    mod2sparse *matrix;
} tanner_graph;

// structure to record nodes forming 4-cycle and shift values in base matrix
typedef struct cyc4_entry
{
    // for each 4-cycle there are two check node: r1, r2, two bits node: c1, c2,
    // and 4 shift values: s1, s2, s3, s4
    int r1, r2, c1, c2; // cycle formed by (r1,c1), (r1,c2), (r2,c1), (r2,c2)
    int *s1; // correspond to row r1 and col c1 position
    int *s2; // correspond to row r1 and col c2 position
    int *s3; // correspond to row r2 and col c1 position
    int *s4; // correspond to row r2 and col c2 position
} cyc4_entry;

// structure to record nodes forming 6-cycle and shift values in base matrix
typedef struct cyc6_entry
{
    // for each 6-cycle there are three check node: r1, r2, r3, three bits node: c1, c2, c3,
    // and 6 shift values: s1, s2, s3, s4, s5, s6
    int r1, r2, r3, c1, c2, c3;
    int *s1, *s2, *s3, *s4, *s5, *s6;
    int *ss1, *ss2, *ss3, *ss4, *ss5, *ss6;
} cyc6_entry;

// structure to record the sets of 4 dependant colomns
typedef struct weight4_entry
{
    int c[4];
} weight4_entry;

// Using Breadth-first search from a given root node, check if there is
// any cycle with given length pass the root
int check_cycle(tanner_graph *g, int root_index, int length);

// Check cycle degree in addition to cycle
int check_cycle_deg(tanner_graph *g, int root_index, int length);

// Trace back n level along the Breadth-first search tree
node_entry trace_back(node_entry u, int n);

// Using BFS to count cycles
void count_cycles(tanner_graph *g, int *cycle);

// Enumerate cycles, and store all 4 cycles in cyc4_rec
void enumerate_cycle(tanner_graph *g, cyc4_entry *cyc4_rec);

// Count cycle degree
int count_cycle_deg(node_entry u, node_entry v, int n);

// sub-rountine used for count_cycles
double BFS_traverse(tanner_graph *g, node_entry *root, int *cycle);

// sub-rountine used for eunmerate cycles
int BFS_traverse_enum(tanner_graph *g, node_entry *root, cyc4_entry *cyc4_rec, int *cyc_cnt);

// This is incorrect way to count cycles. 
// I keep it here for reference as Depth-first search
void DFS_count_cycles(tanner_graph *g, int *cycle);

// Sub-routine used for DFS_count_cycles
void DFS_visit(tanner_graph *g, node_entry *u, int *cycle);

// Check the cycles distribution in a parity check matrix
void check_phck_cycles(char *file, char *log, int *cycle);

// calculate the cycle effect
double cycle_effect(int *cycle, int n);

float avg_girth(tanner_graph *g);

int count_4cyc(mod2sparse *matrix);
int enumerate_4cyc(mod2sparse *matrix, cyc4_entry *cyc4_rec);
void remove_4cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num);

int count_6cyc(mod2sparse *matrix);
int enumerate_6cyc(mod2sparse *matrix, cyc6_entry *cyc6_rec);
void remove_6cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num, cyc6_entry *cyc6_rec, int cyc6_num);

int reduce_6cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num);
void new_remove_6cyc(mod2sparse *basePCH, mod2sparse *exPCH, int p, cyc4_entry *cyc4_rec, int cyc4_num, cyc6_entry *cyc6_rec, int cyc6_num);

int count_weight4(mod2sparse *basePCH);

void expand_matrix(mod2sparse *exMatrix, mod2sparse *baseMatrix, int exFactor);
int chk_expand_weight4(mod2sparse *mat, int p, weight4_entry w4_rec);

int cyc4_cal(mod2sparse* basePCH, int qc_size);
int cyc4_cal_mask(mod2sparse* basePCH, int** mask_matrix, int qc_size, int pad_bit);
int cyc4_cal_ibex(mod2sparse* basePCH, int** fade_matrix, int qc_size, int pad_bit);
int cyc6_cal(mod2sparse* basePCH, int qc_size);
int cyc6_cal_mask(mod2sparse* basePCH, int** mask_matrix, int qc_size, int pad_bit);
int cyc6_cal_ibex(mod2sparse* basePCH, int** fade_matrix, int qc_size, int pad_bit);

#endif
