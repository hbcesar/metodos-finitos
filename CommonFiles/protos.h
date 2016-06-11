#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "heads.h"

#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))

/*----------------------------------------------------------------------------
 * MATRIX HEADER FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern int      COMPARE_array            (const void * a, const void * b);
extern void     MATRIX_readCSR           (MAT* A, char* p);
extern double   MATRIX_aij               (MAT* A, int i, int j);
extern void     MATRIX_printCSR          (MAT* A);
extern void     MATRIX_printFULL         (MAT* A);
extern long int MATRIX_envelope          (MAT* A);
extern long int MATRIX_bandwidth         (MAT* A);
extern void     MATRIX_clean             (MAT* A);
extern void     MATRIX_matvec            (MAT* A, double* x, double* b);
extern void     MATRIX_forward           (MAT* L, double* b, double* y);
extern void     MATRIX_backward          (MAT* U, double* y, double* x);
extern void     MATRIX_permutation       (MAT* A, int* p);
extern void     MATRIX_writeCSR          (MAT* A, double* f, int* s, int nP, int bandwidth);
extern void     MATRIX_printLU           (MAT* A, MAT* L, MAT* U);

/*----------------------------------------------------------------------------
 * GRAPH FUNCTIONS PROTOTYPE IN CSR FORMAT
 *--------------------------------------------------------------------------*/
extern int      COMPARE_degr_ASC         (const void * a, const void * b);
extern int      COMPARE_dist_degr_DES    (const void * a, const void * b);
extern int      COMPARE_dist_degr_ASC    (const void * a, const void * b);
extern int      GRAPH_degree             (MAT* A, int x);
extern int*     GRAPH_adjacent           (MAT* A, int x);
extern int*     GRAPH_bfs                (MAT* A, int x, int* dist);
extern int*     GRAPH_bfs_RCM            (MAT* A, int x, int* dist);
extern int      GRAPH_LS_depth           (int* LS, int n);
extern int      GRAPH_LS_width           (int* LS, int n);
extern LIST*    GRAPH_LS_last_level      (MAT* A, int* LS, int n);
extern int*     GRAPH_LS_peripheral      (MAT* A, int *node_s, int* node_e);

/*----------------------------------------------------------------------------
 * LINKED LIST FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern LIST*    LIST_insert_IF_NOT_EXIST (LIST* L, int x);
extern LIST*    LIST_remove              (LIST* L, int x);
extern LIST*    LIST_remove_first        (LIST* L);
extern void     LIST_print               (LIST* L);
extern int      LIST_first               (LIST* L);
extern void     LIST_destroy             (LIST* L);

/*----------------------------------------------------------------------------
 * PRECONDITIONERS FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern void     SPARMAT_setup            (SparMAT* mat, int n);
extern void     SPARILU_setup            (SparILU* lu, int n);
extern void     SPARILU_row              (SparILU* lu, int nrow);
extern void     CSRto_SPARMAT            (MAT* A, SparMAT* mat);
extern void     SPARILU_toCSR            (SparILU* lu, MAT* L, MAT* U);
extern int      LEVEL_OF_FILL            (SparMAT* csmat, SparILU* lu, int p);
extern void     ILUP                     (SparMAT* csmat, SparILU* lu, int p);
extern void     SPARMAT_clean            (SparMAT* mat);
extern void     SPARILU_clean            (SparILU* lu);
extern void     SPARILU_print            (SparILU* lu);
extern void     QSPLIT                   (double *a, int *ind, int n, int Ncut);
extern void     ILUT                     (SparMAT* csmat, SparILU* lu, int lfil, double tol);

/*----------------------------------------------------------------------------
 * REORDERING FUNCTIONS PROTOTYPE
 *--------------------------------------------------------------------------*/
extern void     REORDERING_RCM_opt       (MAT* A, int** p);
extern void     REORDERING_RCM           (MAT* A, int** p);
