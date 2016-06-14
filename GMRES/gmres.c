#include <stdlib.h>
#include <stdio.h>
#include "gmres.h"

#include "../CommonFiles/heads.h"

/* destroy a solution struct */
void delete_solution(Solution s)
{
    DeleteVector(s.x);
}

//Multiplicacao Matriz x Vetor em CSR, onde MAT* a é a matriz CSR e Vector b é o vetor.
void matrix_vector_multiply_CSR(MAT* A, Vector b, Vector result)
{
    if (A->m != b.size || b.size != result.size)
    {
        printf("\nWrong dimensions! The A matrix collumns and the vectors B and X dimension must all have the same size!\n");
        exit(-30);
    }

    int i, j, l_begin, l_end;
    int n = b.size;

    /* syntactic sugar */
    double *AA = A->AA;
    int *IA = A->IA;
    int *JA = A->JA;
    double *r = result.v;
    double *bv = b.v;

    /* show the matrix */
    for(i = 0; i < n; i++)
    {
        r[i] = 0.0;

        l_begin = IA[i];
        l_end = IA[i+1];

        for (j = l_begin; j < l_end; j++)
        {
            r[i] += AA[j] * bv[JA[j]];
        }
    }

    return;
}

/* solve LUx = b */
void lu_solver(MAT *L, MAT *U, Vector b, Vector result)
{
    if (L->n != U->n || L->n != b.size || b.size != result.size)
    {
        printf("\nWrong dimensions! The L/U matrix collumns and the vectors B and X dimension must all have the same size!\n");
        exit(-31);
    }

    /* helpers */
    int i, j, l_begin, l_end;
    /* get the system dimension */
    unsigned int n = b.size;

    /* syntactic sugar */
    double *AA = L->AA;
    int *IA = L->IA;
    int *JA = L->JA;
    double *D = U->D;

    /* get direct access to b and result vectors */
    double *bv = b.v;
    double *rv = result.v;

    /* temp vector */
    Vector p = BuildVector(n);
    /* the p vector direct access */
    double *pv = p.v;

    /* solve Lp= b */
    for (i = 0; i < n; ++i)
    {
        /* copy the current value */
        pv[i] = bv[i];

        /* get the line elements boundary */
        l_begin = IA[i];
        l_end = IA[i+1];

        for (j = l_begin; j < l_end; ++j)
        {
            /* update the current value */
            pv[i] -= AA[j] * pv[JA[j]];
        }

        /* L->D is the diagonal with all values equals to 1.0 */
        /* so we don't need to apply the usual division' */

    }

    /* change the matrixes pointers */
    AA = U->AA;
    JA = U->JA;
    IA = U->IA;

    /* solve Ux = p */
    for (i = n - 1; 0 <= i; --i)
    {
        /* copy the current value */
        rv[i] = pv[i];

        /* get the line elements boundary */
        l_begin = IA[i];
        l_end = IA[i+1];

        for (j = l_begin; j < l_end; ++j)
        {
            /* update the current value */
            rv[i] -= AA[j] * rv[JA[j]];
        }

        /* the D is the U matrix diagonal */
        rv[i] /= D[i];

    }

    /* remove the temp vector */
    DeleteVector(p);

    return;
}

// Reordena vetor b' para obter b
Vector rearrange_solution(Vector b, int* p)
{
    int i;
    unsigned int size = b.size;
    Vector vB = BuildVectorWithValue(size, 1);

    //acucar sintatico
    double* b_linha = b.v;
    double* b_original = vB.v;

    for (i = 0; i < size; i++)
    {
        b_original[i] = b_linha[p[i]];
    }

    return vB;
}

/* the GMRES solver */
Solution gmres_solver(MAT *A, Vector b, double tol, unsigned int kmax, unsigned int lmax)
{

    /* the common helpers */
    int i, iplus1, j, jplus1, k, kmax1, iter, iter_gmres;
    unsigned int n = b.size;
    double rho, r, tmp;
    double *hv, *prev_v, *next_v;
    double hji, hj1i;

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
    double *c = (double*) calloc(kmax, sizeof(double));
    double *s = (double*) calloc(kmax, sizeof(double));

    /* allocate the y array */
    double *y = (double*) calloc(kmax1, sizeof(double));

    /* allocate the error array, starts with zero value */
    double *e = (double*) malloc(kmax1*sizeof(double));

    /* build the u vector array */
    Vector u[kmax1];

    /* build the h vector array */
    Vector h[kmax1];

    /* allocate each u and h vector inside the array'*/
    for (i = 0; i < kmax1; ++i)
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
        /* first let's find the residual/error vector */
        /* start */

        matrix_vector_multiply_CSR(A, sol.x, u[0]);

        /* let's remember: r = b - Ax */
        /* but we got now just the Ax, so ... */
        for (j = 0; j < n; ++j) {
            /* subtract = b - Ax */
            rv[j] = bv[j] - rv[j];
        }

        /* end */

        /* we need the euclidean norm (i.e the scalar error value) */
        rho = EuclideanNorm(u[0]);

        if (rho != 0.0) {

            /* let's normalize the residual vector */
            ScaleVector(u[0], 1.0/rho);

        } else {

            /* we got a trivial solution */
            i += +1;
            for (i = 0; i < n; i++)
            {
                x0[i] = 0.0;
            }

            break;
        }

        /* update the rho value */
        e[0] = rho;
        /* reset the error */
        for (j = 1; j < kmax1; ++j)
        {
            e[j] = 0.0;
        }

        /* reset the h matrix */
        for (j = 0; j < kmax1; ++j)
        {
            hv = h[j].v;
            for (k = 0; k < kmax1; ++k)
            {
                hv[k] = 0.0;
            }
        }

        /* reset the i counter */
        i = 0;

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

            for (j = 0; j < iplus1; ++j)
            {
                /* get the prev vector direct access */
                prev_v = u[j].v;

                tmp = InnerProduct(u[j], u[iplus1]);
                hv[j] = tmp;

                /* update the current direction vector */
                for (k = 0; k < n; k++)
                {
                    next_v[k] -= tmp * prev_v[k];
                }
            }

            /* get the euclidean norm of the last direction vector */
            tmp = EuclideanNorm(u[iplus1]);

            /* update the next h vector */
            hv[iplus1] = tmp;

            if (0 == tmp)
            {
                kmax = i;
                i += 1;
                break;
            }

            /* normalize the direction vector */
            ScaleVector(u[iplus1], 1.0/tmp);

            /* GRAM-SCHMIDT end */

            /* QR algorithm */
            if (0 < i)
            {
                for (j = 0, jplus1 = 1; j < i; ++j, ++jplus1)
                {
                    hji = c[j]*hv[j] + s[j]*hv[jplus1];
                    hj1i = -s[j]*hv[j] + c[j]*hv[jplus1];
                    hv[j] = hji;
                    hv[jplus1] = hj1i;
                }

            }

            /* update the residual value */
            r = sqrt(hv[i]*hv[i] + hv[iplus1]*hv[iplus1]);

            /* update the cosine value */
            c[i] = hv[i]/r;

            /* update the sine value */
            s[i] = hv[iplus1]/r;

            /* update the current position inside the h vector */
            hv[i] = r;

            /* update the next position inside the h vector */
            hv[iplus1] = 0.0;

            /* rotate the error  */
            e[iplus1] = -s[i]*e[i];
            e[i] *= c[i];

            /* update the rho value, the current error */
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

            for (j = iter_gmres; j > i; --j)
            {
                e[i] -= h[j].v[i] * y[j];
            }
            y[i] = e[i]/hv[i];

        }

        iter_gmres += 1;
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < iter_gmres; ++k)
            {
                /* update the solution vector */
                /* it's so tricky! */
                x0[j] += y[k] * u[k].v[j];
            }

        }

        if (rho < epson)
        {
            break;
        }

    }

    /* update the iteration counter*/
    sol.iterations = iter + 1;

    /* remove the auxiliary arrays */
    free(c);
    free(s);
    free(y);
    free(e);

    /* remove the h and u vectors */
    for (i = 0; i < kmax1; ++i)
    {
        DeleteVector(u[i]);
        DeleteVector(h[i]);
    }

    return sol;
}

/* the GMRES solver */
Solution gmres_lu(MAT *A, MAT *L, MAT *U, Vector b, double tol, unsigned int kmax, unsigned int lmax)
{

    /* the common helpers */
    int i, iplus1, j, jplus1, k, kmax1, iter, iter_gmres;
    unsigned int n = b.size;
    double rho, r, tmp;
    double *hv, *prev_v, *next_v;
    double hji, hj1i;

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
    double *c = (double*) calloc(kmax, sizeof(double));
    double *s = (double*) calloc(kmax, sizeof(double));

    /* allocate the y array */
    double *y = (double*) calloc(kmax1, sizeof(double));

    /* allocate the error array, starts with zero value */
    double *e = (double*) malloc(kmax1*sizeof(double));

    /* build the u vector array */
    Vector u[kmax1];

    /* build the h vector array */
    Vector h[kmax1];

    /* allocate each u and h vector inside the array'*/
    for (i = 0; i < kmax1; ++i)
    {
        u[i] = BuildVector(n);
        h[i] = BuildVector(kmax1);
    }

    /* build the auxiliary vector */
    Vector aux = BuildVector(n);

    /* get the auxiliary vector direct access */
    double *av = aux.v;

    /* the GMRES main outside loop */
    /* we are going to break this loop */
    /* if we get a tiny little error */
    for (iter = 0; iter < lmax; ++iter)
    {
        /* first let's find the residual/error vector */
        /* start */
        matrix_vector_multiply_CSR(A, sol.x, aux);

        /* let's remember: r = b - Ax */
        /* we have r' = M'(b - Ax) */
        /* but we got now just the Ax, so ... */
        for (j = 0; j < n; ++j) {
            /* subtract = b - Ax */
            av[j] = bv[j] - av[j];
        }

        /* solve Mr' = r */
        /* let's save the result in the first direction vector' */
        lu_solver(L, U, aux, u[0]);

        /* end */

        /* we need the euclidean norm (i.e the scalar error value) */
        rho = EuclideanNorm(u[0]);

        if (rho != 0.0) {

            /* let's normalize the residual vector */
            ScaleVector(u[0], 1.0/rho);

        } else {

            /* we got a trivial solution */
            i += +1;
            for (i = 0; i < n; i++)
            {
                x0[i] = 0.0;
            }

            break;
        }

        /* update the rho value */
        e[0] = rho;
        /* reset the error */
        for (j = 1; j < kmax1; ++j)
        {
            e[j] = 0.0;
        }

        /* reset the h matrix */
        for (j = 0; j < kmax1; ++j)
        {
            hv = h[j].v;
            for (k = 0; k < kmax1; ++k)
            {
                hv[k] = 0.0;
            }
        }

        /* reset the i counter */
        i = 0;

        /* the internal loop, for restart purpose */
        while (rho > epson && i < kmax)
        {
            /* set the iplus value */
            iplus1 = i+1;

            /* get the next direction vector */
            matrix_vector_multiply_CSR(A, u[i], aux);

            /* get the next direction vector */
            lu_solver(L, U, aux, u[iplus1]);

            /* Gram-Schmidt process */
            /* GRAM-SCHMIDT begin */
            /* update the h vector */
            hv = h[i].v;

            /* get the next vector direct access */
            next_v = u[iplus1].v;

            for (j = 0; j < iplus1; ++j)
            {
                /* get the prev vector direct access */
                prev_v = u[j].v;

                tmp = InnerProduct(u[j], u[iplus1]);
                hv[j] = tmp;

                /* update the current direction vector */
                for (k = 0; k < n; k++)
                {
                    next_v[k] -= tmp * prev_v[k];
                }
            }

            /* get the euclidean norm of the last direction vector */
            tmp = EuclideanNorm(u[iplus1]);

            /* update the next h vector */
            hv[iplus1] = tmp;

            if (0 == tmp)
            {
                kmax = i;
                i += 1;
                break;
            }

            /* normalize the direction vector */
            ScaleVector(u[iplus1], 1.0/tmp);

            /* GRAM-SCHMIDT end */

            /* QR algorithm */
            if (0 < i)
            {
                for (j = 0, jplus1 = 1; j < i; ++j, ++jplus1)
                {
                    hji = c[j]*hv[j] + s[j]*hv[jplus1];
                    hj1i = -s[j]*hv[j] + c[j]*hv[jplus1];
                    hv[j] = hji;
                    hv[jplus1] = hj1i;
                }

            }

            /* update the residual value */
            r = sqrt(hv[i]*hv[i] + hv[iplus1]*hv[iplus1]);

            /* update the cosine value */
            c[i] = hv[i]/r;

            /* update the sine value */
            s[i] = hv[iplus1]/r;

            /* update the current position inside the h vector */
            hv[i] = r;

            /* update the next position inside the h vector */
            hv[iplus1] = 0.0;

            /* rotate the error  */
            e[iplus1] = -s[i]*e[i];
            e[i] *= c[i];

            /* update the rho value, the current error */
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

            for (j = iter_gmres; j > i; --j)
            {
                e[i] -= h[j].v[i] * y[j];
            }
            y[i] = e[i]/hv[i];

        }

        iter_gmres += 1;
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < iter_gmres; ++k)
            {
                /* update the solution vector */
                /* it's so tricky! */
                x0[j] += y[k] * u[k].v[j];
            }

        }

        if (rho < epson)
        {
            break;
        }

    }

    /* update the iteration counter*/
    sol.iterations = iter + 1;

    /* remove the auxiliary arrays */
    free(c);
    free(s);
    free(y);
    free(e);

    /* remove the h and u vectors */
    for (i = 0; i < kmax1; ++i)
    {
        DeleteVector(u[i]);
        DeleteVector(h[i]);
    }

    /* remove the auxiliary vector */
    DeleteVector(aux);

    return sol;
}
