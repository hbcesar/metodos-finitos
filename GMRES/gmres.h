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

/* multiply a given CSR matrix to any Vector */
void matrix_vector_multiply_CSR(MAT* A, Vector b, Vector result);

/* the GMRES solver */
Solution gmres_solver(MAT *A, Vector b, double tol, unsigned int kmax, unsigned int lmax);

/* the GMRES main function */
Solution gmres(MAT *A, MAT *L, MAT *U, Vector b, double tol, unsigned int kmax, unsigned int lmax);

#endif
