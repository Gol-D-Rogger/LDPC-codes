// include vector operations header file
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rand.h"

void vec_clr(char*, int);
void vec_clr_int(int*, int);
void vec_set(char*, int);
void vec_copy(char*, char*, int, int, int);
int vec_cmp(char*, char*, int, int, int);
void vec_add(char*, char*, char*, int);
void vec_incr(char*, char *, int);
void vec_incr_int(int*, char*, int);
void vec_mod2_add(char*, char*, char*, int);
void vec_shift(char *, char *, int , int);
int vec_sum(char*, int);
void vec_print(char *, int , FILE *);
int vec_max(char*, int);
int vec_find(char*, int, int, int*);

void vec_copy_fill(char*, int, char*, int, int, int, char);
void vec_mask(char*, int, int, char);