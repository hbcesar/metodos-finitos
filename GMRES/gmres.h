#ifndef GMRES_H
#define GMRES_H

#include "../CommonFiles/protos.h"
#include "../CommonFiles/Vector/Vector.h"

typedef struct Solution {

    /* the output vector */
    Vector x;

    /* how many iterations */
    unsigned int iterations;

};

/* multiply a given CSR matrix to any Vector */
void matrix_vector_multiply_CSR(MAT* A, Vector b, Vector result);

/* the GMRES solver */
Solution gmres_solver(Mat *A, Vector b, double tol, double kmax, double lmax);

/* the GMRES main function */
Solution gmres(Mat *A, Mat *L, Mat *U, Vector b, double tol, unsigned int kmax, unsigned int nmax);

#endif
