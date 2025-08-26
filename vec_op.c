#include "vec_op.h"

// Clear a vector (set all elements to 0)
void vec_clr(char *vec, int len) {
    memset(vec, 0, len * sizeof(char));
}

// Set all elements of a vector to 1
void vec_set(char *vec, int len) {
    memset(vec, 1, len * sizeof(char));
}

// Copy a portion of one vector to another
void vec_copy(char *dest, char *src, int dest_offset, int src_offset, int len) {
    memmove(dest + dest_offset, src + src_offset, len * sizeof(char));
}

// Compare two vectors
// Returns <0 if vec1 < vec2, 0 if equal, >0 if vec1 > vec2 (same as memcmp)
int vec_cmp(char *vec1, char *vec2, int offset1, int offset2, int len) {
    return memcmp(vec1 + offset1, vec2 + offset2, len * sizeof(char));
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

// Perform element-wise modulo-2 addition
void vec_mod2_add(char *result, char *vec1, char *vec2, int len) {
    for (int i = 0; i < len; i++) {
        result[i] = (vec1[i] + vec2[i]) % 2;
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