#include <stdlib.h>
#include <stdio.h>
#include "gmres.h"

#include "../CommonFiles/heads.h"

//Multiplicacao Matriz x Vetor em CSR, onde MAT* a é a matriz CSR e Vector b é o vetor.
Vector matrix_vector_multiply_CSR(MAT* A, Vector b)
{
    if (A->m != b.size)
    {
        printf("\nWrong dimensions! The A matrix collumns and the b vector size should be the same!\n");
        exit(-30);
    }

    int i, j, l_begin, l_end, col;
    int n = A->n;
    Vector result = BuildVector(b.size);

    // syntactic sugar
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

    return result;

}
