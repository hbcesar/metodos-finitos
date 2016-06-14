#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "./CommonFiles/protos.h"
#include "CommonFiles/Vector/Vector.h"
#include "GMRES/gmres.h"

double get_time ()
{
    struct timeval tv; gettimeofday(&tv, 0);
    return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}


int main (int argc, char* argv[])
{
    double time;

    if (argc != 2)
    {
        printf("\n Erro! Sem arquivo da matriz (.mtx)");
        printf("\n Modo de usar: ./program <nome_da_matriz> Saindo... [main]\n\n");
        return 0;
    }

    MAT *A = (MAT*) malloc (sizeof(MAT));
    MATRIX_readCSR (A,argv[1]);

    /*---------------------------------------------*/
    /*----MULTIPLICANDO A MATRIZ POR UM VETOR------*/
    Vector x = BuildVectorWithValue(A->m, 1);
    Vector result = BuildVector(A->m);
    matrix_vector_multiply_CSR(A, x, result);
    printf("\nThe x vector: ");
    ShowVector(x);
    printf("\nResult: ");
    ShowVector(result);
    DeleteVector(x);
    DeleteVector(result);

    /*---------------------------------------------*/
    /*---------------GMRES SOLVER------------------*/
    printf("\nGMRES solver\n");
    Vector b = BuildVectorWithRandomValues(A->m, 10);
    printf("\nThe b vector:");
    ShowVector(b);

    Solution sol;
    sol = gmres_solver(A, b, 0.001, 5, 2000);
    printf("\nThe x vector:");
    ShowVector(sol.x);
    printf("\nGMRES Iterations: %d\n", sol.iterations);

    delete_solution(sol);
    DeleteVector(b);

    /*---------------------------------------------*/
    /*---COMO USAR O ALGORITMO ILUP----------------*/
    /*---------------------------------------------*/
    MAT *L = (MAT*) malloc(sizeof(MAT));						// Alocando matriz L
    MAT *U = (MAT*) malloc(sizeof(MAT));						// Alocando matriz U

    SparMAT* mat = (SparMAT*) malloc(sizeof(SparMAT));				// Alocando estruturas para o ILU(p)
    SparILU* lu  = (SparILU*) malloc(sizeof(SparILU));

    printf("\n  [ CALCULANDO PRECONDICIONADOR ILU ]\n");
    /*---START TIME---------------> */ time = get_time();
    CSRto_SPARMAT (A,mat);								// Convertendo CSR para estrutura especial
    ILUP          (mat,lu,2);							// Algoritmo ILU(p)
    SPARILU_toCSR (lu,L,U);								// Convertendo estrutura especial para CSR
    /*---FINAL TIME---------------> */ time = (get_time() - time)/100.0;
    printf("  - Tempo total              : %.6f sec\n", time);

    SPARILU_clean (lu);								// Liberando memória da estrutura lu
    SPARMAT_clean (mat);								// Liberando memória da estrutura mat

    /* L contém a parte estritamente inferior de M / L->D contém a diagonal = 1.0 */
    /* U contém a parte estritamente superior de M / U->D contém a diagonal       */
    MATRIX_printLU (A,L,U);

    /*---------------------------------------------*/
    /*---COMO USAR O REORDENAMENTO RCM-------------*/
    /*---------------------------------------------*/
    int *p;									// Vetor de permutação
    int  bandwidth;

    bandwidth = (int) MATRIX_bandwidth(A);					// Calcula Largura de Banda da matriz original
    printf("\n  [ REORDENANDO com RCM ]\n");
    printf("  - Largura de Banda inicial : %d\n", bandwidth);

    /*---START TIME---------------> */ time = get_time();
    REORDERING_RCM_opt(A,&p);						// Aplica o reordenamento RCM na matriz A
    MATRIX_permutation(A,p); 						// Aplica a permutação em A para trocar linhas e colunas
    /*---FINAL TIME---------------> */ time = (get_time() - time)/100.0;

    bandwidth = (int) MATRIX_bandwidth(A);					// Calcula Largura de Banda da matriz reordenada
    printf("  - Largura de Banda final   : %d\n", bandwidth);
    printf("  - Tempo total              : %.6f sec\n\n", time);


    //JOSIAS, esse todo foi eu que fiz tá? Nao foi a profs nao.
    /*********************** TODO *******************************
    *
    * 
    *	Função para solucionar o sistema linear
    *
    *	Entrada:
    *		Matriz CSR: MAT *A
    *		Matriz Precondicionamento: MAT *L e MAT *U
    *		Vetor de Termos Independentes:
    * 		Tolerância: double tol
    *		Número máximo de vetores na base de Krylov:
    *		Número máximo de Iteraçoes:
    *		(caso o vetor de termos independentes não esteja permutado,  é necessário enviar também o vetor de permutações)
    *
    *
    *	Saída: Solução do Sistema Linear e o Número total de iterações realizadas pelo GMRES.
    *
    *************************************************************/


    free(p);
    MATRIX_clean(A);
    MATRIX_clean(L);
    MATRIX_clean(U);

    return 0;
}
