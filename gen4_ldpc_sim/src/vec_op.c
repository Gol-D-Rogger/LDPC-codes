#include "vec_op.h"

// Clear a vector (set all elements to 0)
void vec_clr(char *vec, int len) {
    for (int i = 0; i < len; i++) {
        vec[i] = 0;
    }
}

// Clear an integer vector (set all elements to 0)
void vec_clr_int(int *vec, int len) {
    memset(vec, 0, len * sizeof(int));
}

// Set all elements of a vector to 1
void vec_set(char *vec, int len) {
    memset(vec, 1, len * sizeof(char));
}

// Copy a portion of one vector to another
void vec_copy(char *vs, char *vd, int vs_si, int vd_si, int d) {
    for (int i = 0; i < d; i++) {
        vd[vd_si + i] = vs[vs_si + i];
    }
}

// Compare two vectors
// Returns <0 if vec1 < vec2, 0 if equal, >0 if vec1 > vec2 (same as memcmp)
int vec_cmp(char *vs, char *vd, int vs_si, int vd_si, int d) {
    int mismatch = 0;

    for (int i = 0; i < d; i++) {
        if (vd[vd_si + i] != vs[vs_si + i]) {
            mismatch = 1;
            break;
        }
    }
    return mismatch;
}

// Add two vectors element-wise
void vec_add(char *vec1, char *vec2, char *result, int len) {
    for (int i = 0; i < len; i++) {
        result[i] = vec1[i] + vec2[i];
    }
}

// Increment a vector by another vector element-wise
void vec_incr(char *vec1, char *vec2, int len) {
    for (int i = 0; i < len; i++) {
        vec1[i] += vec2[i];
    }
}

// Increment an integer vector by a char vector element-wise
void vec_incr_int(int *vec1, char *vec2, int len) {
    for (int i = 0; i < len; i++) {
        vec1[i] += vec2[i];
    }
}

// Perform element-wise modulo-2 addition
void vec_mod2_add(char *va, char *vb, char *vs, int len) {
    for (int i = 0; i < len; i++) {
        vs[i] = (va[i] + vb[i]) % 2;
    }
}

// Shift a vector by one position
void vec_shift(char *vs, char *vd, int d, int shift) {
    for (int i = 0; i < d; i++)
    {
        vd[(d+i+shift)%d] = vs[i];
    }
}

// Sum all elements in a vector
int vec_sum(char *vec, int len) {
    int sum = 0;
    for (int i = 0; i < len; i++) {
        sum += vec[i];
    }
    return sum;
}

// Print a vector to a file
void vec_print(char *pblk, int pblk_len, FILE *stream_w) {
    int tmp;

    for (int i=0; i<pblk_len; i = i+32)
    {
        for (int j=0; j<8; j++)
        {
            tmp = 0;
            for (int k=0; k<4; k++)
            {
                tmp = tmp *2 + pblk[i+j*4+k];
            }

            fprintf(stream_w, "%X", tmp);
        }
        fprintf(stream_w,"\n");
    }
}

// Find maximum value in a vector
int vec_max(char *vec, int len) {
    if (len <= 0) return 0;
    
    int max_val = vec[0];
    for (int i = 1; i < len; i++) {
        if (vec[i] > max_val) {
            max_val = vec[i];
        }
    }
    return max_val;
}

// Find all positions where vector equals target value
// Returns number of positions found
int vec_find(char *vec, int len, int target_val, int *positions) {
    int count = 0;
    for (int i = 0; i < len; i++) {
        if (vec[i] == target_val) {
            positions[count] = i;
            count++;
        }
    }
    return count;
}

void vec_copy_fill(char *vs, int vs_sz, char* vd, int vs_si, int vd_si, int d, char fill_val)
{
    for (int i = 0; i < d; i++)
    {
        if (vs_si + i < vs_sz)
            vd[vd_si + i] = vs[vs_si + i];
        else
            vd[vd_si + i] = fill_val;
    }
}

void vec_mask(char *vec, int start, int mask_len, char mask_val)
{
    for (int i = 0; i < mask_len; i++)
    {
        vec[start + i] = mask_val;
    }
}