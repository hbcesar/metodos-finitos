#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "./CommonFiles/protos.h"
#include "./CommonFiles/Vector/Vector.h"
#include "./GMRES/gmres.h"

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
    /*---COMO USAR O REORDENAMENTO RCM-------------*/
    /*---------------------------------------------*/

    // Vetor de permutação
    int *p;
    int  bandwidth;

    // Calcula Largura de Banda da matriz original
    bandwidth = (int) MATRIX_bandwidth(A);

    printf("\n  [ REORDENANDO com RCM ]\n");
    printf("  - Largura de Banda inicial : %d\n", bandwidth);

    /*---START TIME---------------> */ time = get_time();
    // Aplica o reordenamento RCM na matriz A
    REORDERING_RCM_opt(A,&p);

    // Aplica a permutação em A para trocar linhas e colunas
    MATRIX_permutation(A,p);
    /*---FINAL TIME---------------> */ time = (get_time() - time)/100.0;

    // Calcula Largura de Banda da matriz reordenada
    bandwidth = (int) MATRIX_bandwidth(A);
    printf("  - Largura de Banda final   : %d\n", bandwidth);
    printf("  - Tempo total              : %.6f sec\n\n", time);

    /*---------------------------------------------*/
    /*---COMO USAR O ALGORITMO ILUP----------------*/
    /*---------------------------------------------*/
    // Alocando matriz L
    MAT *L = (MAT*) malloc(sizeof(MAT));

    // Alocando matriz U
    MAT *U = (MAT*) malloc(sizeof(MAT));

    // Alocando estruturas para o ILU(p)
    SparMAT* mat = (SparMAT*) malloc(sizeof(SparMAT));
    SparILU* lu  = (SparILU*) malloc(sizeof(SparILU));

    printf("\n  [ CALCULANDO PRECONDICIONADOR ILU ]\n");

    /*---START TIME---------------> */ time = get_time();
    // Convertendo CSR para estrutura especial
    CSRto_SPARMAT(A,mat);

    // Algoritmo ILU(p)
    ILUP(mat,lu,2);

    // Convertendo estrutura especial para CSR
    SPARILU_toCSR(lu,L,U);

    /*---FINAL TIME---------------> */
    time = (get_time() - time)/100.0;
    printf("  - Tempo total              : %.6f sec\n", time);

    // Liberando memória da estrutura lu
    SPARILU_clean(lu);

    // Liberando memória da estrutura mat
    SPARMAT_clean(mat);

    /* L contém a parte estritamente inferior de M / L->D contém a diagonal = 1.0 */
    /* U contém a parte estritamente superior de M / U->D contém a diagonal       */
    // MATRIX_printLU(A,L,U);

    /*---------------------------------------------*/
    /*---------------GMRES SOLVER------------------*/
    printf("\nGMRES solver\n");
//
//     /* create a vector with all values equal to 1.0 */
    Vector ones = BuildVectorWithValue(A->m, 1.0);

    /* the b independent vector */
    Vector b = BuildVector(A->m);

    /* get the b values - see 4.2 - Leitura de Matrizes */
    matrix_vector_multiply_CSR(A, ones, b);

    /* the solution must provide de x vector and how many iterations */
    /* you can change the Solution struct definition in order to obtain more info from the solver'' */
    Solution sol;

    // caso do gmres sem pré-condicionamento
    //sol = gmres_solver(A, b, 1e-08, 20, 300);

    // caso do gmres com pré-condicionamento
    sol = gmres_lu(A, L, U, b, 10e-08, 20, 600);

    // desfazer a permutação
    Vector x = rearrange_solution(sol.x, p);

    // printf("\nThe x vector:");
    // ShowVector(x);
    // printf("\nGMRES Iterations: %d\n", sol.iterations);

    /* TODO
     * Realizar todos os testes, com as variações definidas no trabalho
     * Salvar as informações apropriadas em um arquivo, se possível já com sintaxe do octave
     */

    delete_solution(sol);
    DeleteVector(b);
    DeleteVector(ones);
    DeleteVector(x);

    free(p);
    MATRIX_clean(A);
    MATRIX_clean(L);
    MATRIX_clean(U);

    return 0;
}
