#ifndef GMRES_H
#define GMRES_H

#include "../CommonFiles/protos.h"
#include "../CommonFiles/Vector/Vector.h"

typedef struct Solution {

    /* the output vector */
    Vector x;

    /* how many iterations */
    unsigned int iterations;

} Solution;

/* destroy a solution struct */
void delete_solution(Solution s);

/* multiply a given CSR matrix to any Vector */
void matrix_vector_multiply_CSR(MAT* A, Vector b, Vector result);

/* undo the vector permutation */
Vector rearrange_solution(Vector b, int* p);

/* solve LUx = b */
void lu_solver(MAT *L, MAT *U, Vector b, Vector result);

/* the GMRES solver */
Solution gmres_solver(MAT *A, Vector b, double tol, unsigned int kmax, unsigned int lmax);

/* the GMRES main function */
Solution gmres_lu(MAT *A, MAT *L, MAT *U, Vector b, double tol, unsigned int kmax, unsigned int lmax);

#endif
