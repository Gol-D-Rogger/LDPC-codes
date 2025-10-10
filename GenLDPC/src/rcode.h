#include "mod2sparse.h"
#include "mod2dense.h"

extern mod2sparse *H;	/* Parity check matrix */

extern int M;    /* Number of rows in H */
extern int N;    /* Number of columns in H */

extern char type; /*Type of generator matrix representation */
extern int *cols; /*Order of columns in generator matrix*/

extern mod2sparse *L, *U; /* sparse LU decomposition, if type=='s' */
extern int *rows; /*Order of rows in generator matrix*/

extern mod2dense *G; /*Dense or mixed representation of generator matrix, if type =='d', or type=='M'*/

void read_pchk (char *);
void read_gen (char *, int,int);