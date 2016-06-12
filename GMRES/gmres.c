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
Solution gmres_solver(Mat *A, Vector b, double tol, double kmax, double lmax) {

    /* the common helpers */
    unsigned int i, iplus1, j, jplus1, k, iter;
    unsigned int n = b.size;
    double rho, r, tmp;

    /* syntatic sugar */
    /* access the vector b */
    double bv = b.v;

    /* get the epson value */
    double epson = tol * EuclideanNorm(b);

    /* the solution */
    Solution sol;

    /* build the temporary x vector /*
    /* this vector will be constantly */
    /* updated until we find the solution */
    sol.x = BuildVector(n);
    /* the direct access  */
    double x0 = sol.x.v;

    k = kmax + 1;

    /* allocate the c and s array */
    double *c = (double*) malloc(k*sizeof(double));
    double *s = (double*) malloc(k*sizeof(double));

    /* allocate the error array */
    double *e = (double*) malloc(k*sizeof(double));

    /* build the u vector array */
    Vector u[kmax];
    /* build the h vector array */
    Vector h[kmax];
    /* allocate each u and h vector inside the array'*/
    for (i = 0; i < k; ++i)
    {
        u[i] = BuildVector(n);
        h[i] = BuildVector(k);
    }

    /* get the residual direct access */
    double *rv = u[0].v;

    /* temp vector */
    Vector tmp = BuildVector(n);

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
            /* start */
            /* update the h vector */
            double *hv = h[i].v;

            /* get the next vector direct access */
            double *next_v = u[iplus1].v;

            for (j = 0; j <= i; j++)
            {
                /* get the prev vector direct access */
                double *prev_v = u[j].v;

                tmp = InnerProduct(u[iplus1], u[j]);
                hv[j] = tmp;

                for (k = 0; k < n; k++)
                {
                    next_v[k] -= tmp*prev_v[k];
                }
            }

            tmp = EuclideanNorm(u[iplus1]);
            hv[i+1] = tmp;
            ScaleVector(u[iplus1], 1.0/tmp);

            /* end */

            /* QR algorithm */
            for (j = 0, jplus1 = 1; j < i; ++j, ++jplus1)
            {
                hv[j] = c[j]*hv[j] + s[j]*hv[jplus1];
                hv[jplus1] = -s[j]*hv[j] + c[j]*hv[jplus1];
            }

            r = sqrt(hv[i]*hv[i] + hv[iplus1]*hv[iplus1]);

            c[i] = hv[i]/r;

            s[i] = hv[iplus1]/r;

            hv[i] = r;

            hv[iplus1] = 0.0;

            e[iplus1] = -s[i]*e[i];

            e[i] = c[i]*e[i];

            rho = fabs(e[iplus1]);

            /* update the i counter */
            i += 1;
        }

        i -= 1;
        for (j = i - 1; j >= 0; --j)
        {
            // y[j] = e[j] -
        }

    }

    /* update the iteration counter*/
    sol.iterations = iter;

    return sol;
}

/* the GMRES solver */
Solution gmres(Mat *A, Mat *L, Mat *U, Vector b, double tol, unsigned int kmax, unsigned int nmax)
{
    int v_size = b.size;

    Solution s;
    s.iterations = 0;
    s.x = BuildVector(v_size);

    /* syntactic sugar */
    double n = A->n;

    /* the temporary vector, this vector is constantly updated */
    /* until we reach the desired solution */
    double x0 = (double*) calloc(n, sizeof(double));

    /*  */
    double y = (double*) calloc(kmax, sizeof(double));

    /* the residual vector */
    Vector r = BuildVector(n);

    Mat V;
    V=eye(n,kmax);
    H=zeros(kmax+1,kmax);
    nitertotal = 0;


    while( l <= lmax )
    {

        r=b-A*x0;
        beta=norm(r);
        V(:,1)=r/beta;

        ebar = eye(kmax+1,1)*beta;
        epson = tol * norm(b);

        niter = 0;

        for i=1:kmax
            niter++;
            w=A*V(:,i);
            for j=1:i
                H(j,i)=V(:,j)'*w;
                w=w-H(j,i)*V(:,j);
            end
            H(i+1,i)=norm(w);
            if H(i+1,i)==0
                kmax=i;
                break
            end
            V(:,i+1)=w/H(i+1,i);

        % ...Aplica as rotações de Givens
        % ...Algoritmo QR

            if( i > 1)
                for j=1:i-1
                    hji = c(j)*H(j,i) + s(j)*H(j+1,i);
                    hj1i = -s(j)*H(j,i) + c(j)*H(j+1,i);
                    H(j,i) = hji;
                    H(j+1,i) = hj1i;
                end
            end

            r = H(i,i)*H(i,i) + H(i+1,i)*H(i+1,i);
            r = sqrt(r);

            c(i) = H(i,i)/r;
            s(i) = H(i+1,i)/r;
            H(i,i) = r;
            H(i+1,i) = 0.0;

            ebar(i+1) = -s(i)*ebar(i);
            ebar(i) = c(i)*ebar(i);

            nitertotal = nitertotal+1;
            if( abs(ebar(i+1)) < epson )
                break;
            end
        end % Fim da iteração GMRES

        residuo = abs(ebar(niter+1));

% ... Resolve o sistema liner Hy=e

        y(niter) = ebar(niter)/H(niter,niter);
        i = niter-1;
        while(i >= 1)
            j = niter;

            while( j > i )

                ebar(i) -= H(i,j)*y(j);
                j--;

            end;

            y(i) = ebar(i)/H(i,i);
            i--;

        end

        x = x0+V(:,1:niter)*y(1:niter);

        if( residuo < epson )
            break;
        else
            x0 = x;
        end
        l++;

    }

    printf("(GMRES) Ciclo %d: iteracao = %d\n",l,niter);

    return s;
}
