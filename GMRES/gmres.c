#include <stdlib.h>
#include <stdio.h>
#include "gmres.h"

#include "../CommonFiles/heads.h"

//Multiplicacao Matriz x Vetor em CSR, onde MAT* a é a matriz CSR e Vector b é o vetor.
void matrix_vector_multiply_CSR(MAT* A, Vector b, Vector result)
{
    if (A->m != b.size || b.size != result.size)
    {
        printf("\nWrong dimensions! The A matrix collumns and the vectors b and result sizes must all be the same!\n");
        exit(-30);
    }

    int i, j, l_begin, l_end, col;
    int n = A->n;

    /* syntactic sugar */
    double *AA = A->AA;
    int *IA = A->IA;
    int *JA = A->JA;
    double *r = result.v;
    double *v = b.v;

    for(i = 0; i < n; i++)
    {
        r[i] = 0;

        l_begin = IA[i];
        l_end = IA[i+1];

        for (j = l_begin; j < l_end; j++)
        {
            col = JA[j];
            r[i] += AA[col] * v[col];
        }

    }

}
/* the GMRES solver */
Solution gmres_solver(MAT *A, Vector b, double tol, unsigned int kmax, unsigned int lmax) {

    /* the common helpers */
    unsigned int i, iplus1, j, jplus1, k, kmax1, iter, iter_gmres;
    unsigned int n = b.size;
    double rho, r, tmp;
    double *hv, *uv, *prev_v, *next_v;

    kmax1 = kmax + 1;

    /* syntatic sugar */
    /* access the vector b */
    double *bv = b.v;

    /* get the epson value */
    double epson = tol * EuclideanNorm(b);

    /* the solution */
    Solution sol;

    /* build the temporary x vector */
    /* this vector will be constantly */
    /* updated until we find the solution */
    sol.x = BuildVector(n);
    /* the direct access  */
    double *x0 = sol.x.v;

    /* allocate the c and s array */
    double *c = (double*) malloc(kmax*sizeof(double));
    double *s = (double*) malloc(kmax*sizeof(double));

    /* allocate the y array */
    double *y = (double*) malloc(kmax1*sizeof(double));

    /* allocate the error array, starts with zero value */
    double *e = (double*) calloc(kmax1, sizeof(double));

    /* build the u vector array */
    Vector u[kmax];

    /* build the h vector array */
    Vector h[kmax];

    /* allocate each u and h vector inside the array'*/
    for (i = 0; i < kmax; ++i)
    {
        u[i] = BuildVector(n);
        h[i] = BuildVector(kmax1);
    }

    /* get the residual direct access */
    double *rv = u[0].v;

    /* the GMRES main outside loop */
    /* we are going to break this loop */
    /* if we get a tiny little error */
    for (iter = 0; iter < lmax; ++iter)
    {
        /* reset the i counter */
        i = 0;

        /* first let's find the residual/error vector */
        /* start */

        matrix_vector_multiply_CSR(A, sol.x, u[i]);

        /* let's remember: r = b - Ax */
        /* but we got now just the Ax, so ... */
        for (j = 0; j < n; ++j) {
            /* subtract = b - Ax */
            rv[j] = bv[j] - rv[j];
        }

        /* end */

        /* we need the euclidean norm (i.e the scalar error value) */
        tmp = EuclideanNorm(u[i]);

        /* let's normalize the residual vector */
        ScaleVector(u[i], 1.0/tmp);

        /* update the rho value */
        rho = e[i] = tmp;

        /* the internal loop, for restart purpose */
        while (rho > epson && i < kmax)
        {
            /* set the iplus value */
            iplus1 = i+1;

            /* get the next direction vector */
            matrix_vector_multiply_CSR(A, u[i], u[iplus1]);

            /* Gram-Schmidt process */
            /* GRAM-SCHMIDT begin */
            /* update the h vector */
            hv = h[i].v;

            /* get the next vector direct access */
            next_v = u[iplus1].v;

            for (j = 0; j <= i; j++)
            {
                /* get the prev vector direct access */
                prev_v = u[j].v;

                tmp = InnerProduct(u[iplus1], u[j]);
                hv[j] = tmp;

                for (k = 0; k < n; k++)
                {
                    next_v[k] -= tmp * prev_v[k];
                }
            }

            /* get the euclidean norm of the last direction vector */
            tmp = EuclideanNorm(u[iplus1]);

            /* update the next h vector */
            hv[i+1] = tmp;

            /* normalize the direction vector */
            ScaleVector(u[iplus1], 1.0/tmp);

            /* GRAM-SCHMIDT end */

            /* QR algorithm */
            for (j = 0, jplus1 = 1; j < i; ++j, ++jplus1)
            {
                hv[j] = c[j]*hv[j] + s[j]*hv[jplus1];
                hv[jplus1] = -s[j]*hv[j] + c[j]*hv[jplus1];
            }

            /* update the residual value */
            r = sqrt(hv[i]*hv[i] + hv[iplus1]*hv[iplus1]);

            c[i] = hv[i]/r;

            s[i] = hv[iplus1]/r;

            hv[i] = r;

            hv[iplus1] = 0.0;

            e[iplus1] = -s[i]*e[i];

            e[i] *= c[i];

            rho = fabs(e[iplus1]);

            /* update the i counter */
            i += 1;
        }

        /* get the iteration counter */
        iter_gmres = i - 1;

        /* update the y value */
        y[iter_gmres] = e[iter_gmres]/h[iter_gmres].v[iter_gmres];

        for (i = iter_gmres - 1; 0 <= i; --i)
        {
            /* update the h vector direct access */
            hv = h[i].v;

            for (j = i + 1; j < iter_gmres + 1; ++j)
            {
                e[i] -= hv[j] * y[j];
            }

            y[i] = e[i]/hv[i];
        }

        for (i = 0; i < n; ++i)
        {
            /* get the current u vector */
            uv = u[i].v;

            for (j = 0; j < kmax + 1; ++j)
            {
                /* update the solution vector */
                x0[i] += uv[j] * y[j];
            }
        }

        /* reset the error */
        for (i = 0; i < kmax1; ++i)
        {
            e[i] = 0.0;
        }

        /* reset the h matrix */
        for (i = 0; i < kmax; ++i)
        {
            hv = h[i].v;
            for (j = 0; j < kmax1; j++)
            {
                hv[j] = 0.0;
            }
        }

    }

    /* update the iteration counter*/
    sol.iterations = iter;

    return sol;
}

/* the GMRES solver */
Solution gmres(MAT *A, MAT *L, MAT *U, Vector b, double tol, unsigned int kmax, unsigned int lmax)
{
    return gmres_solver(A, b, tol, kmax, lmax);
}
