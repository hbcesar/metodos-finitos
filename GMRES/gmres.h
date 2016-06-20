#ifndef GMRES_H
#define GMRES_H

#define true 1
#define false 0

#include "../CommonFiles/protos.h"
#include "../CommonFiles/Vector/Vector.h"

typedef struct GMRES_Parameters
{
    /* how many iterations */
    unsigned int lmax;

    /* how many vectors in Krylov space */
    /* each element here is a different krylov set */
    /* e.g kmax = [20, 50, 100] */
    unsigned int kmax[3];

    /* how many krylov spaces */
    unsigned int krilov_space_index;

    /* the tolerance */
    double tol;

    /* precondiotinig flag */
    /* precondiotioner = [false, true] */
    int preconditioner;

    /* fill-in values */
    /* fill_in = [0, 2, 10] */
    unsigned int fill_in[3];

    /* the current fill_in */
    unsigned int ilu;

    /* reordering flag */
    /* [false, true] */
    int reordering;

    /* the current solving matrix name */
    char* name;

} GMRES_Parameters, *GMRES_ParametersPtr;

/* build the base case parameters */
GMRES_ParametersPtr get_base_parameters();

typedef struct Solution
{

    /* the output vector */
    Vector x;

    /* how many iterations */
    unsigned int iterations;

    /* spent time */
    double time;

    /* vector with rho infos */
    Vector rhos;

    /* the rhos history vector size */
    unsigned int rho_last_index;

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
