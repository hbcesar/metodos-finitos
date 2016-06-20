#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "./CommonFiles/protos.h"
#include "./CommonFiles/Vector/Vector.h"
#include "./GMRES/gmres.h"



double get_time ()
{
    struct timeval tv; gettimeofday(&tv, 0);
    return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}

MAT* read_matrix(char *filename)
{
    MAT *A = (MAT*) malloc (sizeof(MAT));
    MATRIX_readCSR (A, filename);

    return A;
}

FILE* open_test_file(char *file)
{
    FILE *f = fopen(file, "a");

    if(f == NULL){
        printf("Could not create the output file. Does the enclosing folder exist?\n...Aborting...\n");
        exit(1);
    }

    return f;
}

void close_test_file(FILE *f)
{
    /* set the plot commands */
    fprintf(f, "plot(rhox_0_0_0_20, rhos_0_0_0_20, rhox_0_1_0_20, rhos_0_1_0_20, rhox_0_1_2_20, rhos_0_1_2_20);\n");
    fprintf(f, "grid minor;\n");
    fprintf(f, "figure;\n");
    fprintf(f, "plot(solx_0_0_0_20, solution_0_0_0_20);\n");
    fprintf(f, "end;\n");

    /* TODO */
    /*
     * We need to place the right plot commands in the octave file
     *
     */

    fclose(f);
}

void printSolution(Solution s, FILE *f, GMRES_ParametersPtr par){

    // unsigned int iterations = s.iterations;
    // double time = s.time;
    double* sol = s.x.v;
    double* r0 = s.rhos.v;

    int i = 0;

    /* helpers */
    /* just to avoid many indirect access */
    unsigned int kmax = par->kmax[par->krilov_space_index];
    unsigned int fill_in = par->fill_in[par->ilu];
    int reordering = par->reordering;
    int preconditioner = par->preconditioner;

    /* set the solution vector */
    /* start */
    if (preconditioner)
    {
        fprintf(f, "iterations_%d_1_%d_%u = %u\n", reordering, fill_in, kmax, s.iterations);

        /* save the current solution required time */
        fprintf(f, "time_%d_1_%d_%u = %lf\n", reordering, fill_in, kmax, s.time);

        /* save the current error x coordinate */
        fprintf(f, "\nrhox_%d_1_%d_%u = 0:%d;\n", reordering, fill_in, kmax, s.iterations - 1);

        /* save the error history */
        fprintf(f, "\nrhos_%d_1_%d_%u = [ ", reordering, fill_in, kmax);

        /* set the rho history vector */
        /* start */
        for(i = 0; i < s.rho_last_index; i++)
        {
            fprintf(f, " %.4lf", log(fabs(r0[i])));
        }

        fprintf(f, " ];\n");
        /* end */

        /* save the current solution x coordinate range */
        fprintf(f, "\nsolx_%d_1_%d_%u = 0:%d;\n", reordering, fill_in, kmax, s.x.size - 1);

        /* save the current solution vector */
        fprintf(f, "\nsolution_%d_1_%d_%u = [ ", reordering, fill_in, kmax);
    }
    else
    {
        /* save t*he current solution statistics */
        fprintf(f, "iterations_%d_0_%u = %u\n", reordering, kmax, s.iterations);

        /* save the current solution required time */
        fprintf(f, "time_%d_0_%u = %lf\n", reordering, kmax, s.time);

        /* save the current error x coordinate */
        fprintf(f, "\nrhox_%d_0_0_%u = 0:%d;\n", reordering, kmax, s.iterations - 1);

        fprintf(f, "\nrhos_%d_0_0_%d = [ ", reordering, kmax);

        /* set the rho history vector */
        /* start */
        for(i = 0; i < s.rho_last_index; i++)
        {
            fprintf(f, " %.4lf", log(fabs(r0[i])));
        }

        fprintf(f, " ];\n");
        /* end */

        /* save the current solution x coordinate range */
        fprintf(f, "\nsolx_%d_0_0_%u = 0:%d;\n ", reordering, kmax, s.x.size - 1);

        /* save the current solution vector */
        fprintf(f, "\nsolution_%d_0_0_%u = [ ", reordering, kmax);
    }

    for(i = 0; i < s.x.size; i++){
        fprintf(f, " %.4lf", sol[i]);

        if (i % 128 == 0)
        {
            fprintf(f, " ... \n");
        }
    }

    fprintf(f, " ];\n\n");
    /* end */

    return;
}

void printTime(GMRES_ParametersPtr par, Solution s){
  //Save the table with each matrix infos in a txt file
  char str[80];
  strcpy (str,"testes/");
  strcat (str,par->name);
  strcat (str,"_log.txt");
  printf("salvando em: %s\n", str);

  FILE *log = fopen( str, "a" );

  unsigned int kmax = par->kmax[par->krilov_space_index];
  unsigned int fill_in = par->fill_in[par->ilu];
  int reordering = par->reordering;
  int preconditioner = par->preconditioner;


  fprintf(log, "Condicionamento: %d\n", preconditioner);
  fprintf(log, "Reordenamento: %d\n", reordering);
  fprintf(log, "ILU Fill in: %d\n", fill_in);
  fprintf(log, "K: %d\n", kmax);
  fprintf(log, "Iterações: %d\n", s.iterations);
  fprintf(log, "Tempo total: %.4lfs\n\n\n", s.time);

  fclose(log);
}


void solver(MAT *A, FILE *output, GMRES_ParametersPtr params)
{
    /* some helpers */
    double time;
    double total_time;

    // Vetor de permutação
    int *p;

    printf("\n  [ GMRES Solver ] \n");
    if (params->reordering)
    {
        /*---------------------------------------------*/
        /*---COMO USAR O REORDENAMENTO RCM-------------*/
        /*---------------------------------------------*/
        //Get the time before the algorithm begins
        total_time = get_time();

        // int  bandwidth;

        // Calcula Largura de Banda da matriz original
        // bandwidth = (int) MATRIX_bandwidth(A);

        // printf("\n  [ REORDENANDO com RCM ]\n");
        // printf("  - Largura de Banda inicial : %d\n", bandwidth);

        /*---START TIME---------------> */ time = get_time();
        // Aplica o reordenamento RCM na matriz A
        REORDERING_RCM_opt(A, &p);

        // Aplica a permutação em A para trocar linhas e colunas
        MATRIX_permutation(A, p);
        /*---FINAL TIME---------------> */ time = (get_time() - time)/100.0;

        // Calcula Largura de Banda da matriz reordenada
        // bandwidth = (int) MATRIX_bandwidth(A);
        // printf("  - Largura de Banda final   : %d\n", bandwidth);
        // printf("  - Tempo total              : %.6f sec\n\n", time);

        printf("  [ Solving with reordering ]\n");
    }
    else
    {
        printf("  [ Solving without reordering ]\n");
    }

    /* declaring the preconditioning variables */
    MAT *L, *U;
    SparMAT *mat;
    SparILU *lu;

    if (params->preconditioner)
    {
        /*---------------------------------------------*/
        /*---COMO USAR O ALGORITMO ILUP----------------*/
        /*---------------------------------------------*/
        // Alocando matriz L
        L = (MAT*) malloc(sizeof(MAT));

        // Alocando matriz U
        U = (MAT*) malloc(sizeof(MAT));

        // Alocando estruturas para o ILU(p)
        mat = (SparMAT*) malloc(sizeof(SparMAT));
        lu = (SparILU*) malloc(sizeof(SparILU));

        // printf("\n  [ CALCULANDO PRECONDICIONADOR ILU ]\n");

        /*---START TIME---------------> */ time = get_time();
        // Convertendo CSR para estrutura especial
        CSRto_SPARMAT(A,mat);

        // Algoritmo ILU(p)
        ILUP(mat, lu, params->fill_in[params->ilu]);

        // Convertendo estrutura especial para CSR
        SPARILU_toCSR(lu,L,U);

        // <------FINAL TIME--------------->
        time = (get_time() - time)/100.0;
        // printf("  - Tempo total              : %.6f sec\n", time);

        // Liberando memória da estrutura lu
        SPARILU_clean(lu);

        // Liberando memória da estrutura mat
        SPARMAT_clean(mat);

        /* L contém a parte estritamente inferior de M / L->D contém a diagonal = 1.0 */
        /* U contém a parte estritamente superior de M / U->D contém a diagonal       */
        // MATRIX_printLU(A,L,U);
    }

    /*---------------------------------------------*/
    /*---------------GMRES SOLVER------------------*/

    /* create a vector with all values equal to 1.0 */
    Vector ones = BuildVectorWithValue(A->m, 1.0);

    /* the b independent vector */
    Vector b = BuildVector(A->m);

    /* get the b values - see 4.2 - Leitura de Matrizes */
    matrix_vector_multiply_CSR(A, ones, b);

    /* the solution must provide de x vector and how many iterations */
    /* you can change the Solution struct definition in order to obtain more info from the solver'' */
    Solution sol;

    /* build the raw solution vector */
    if(params->preconditioner)
    {
        printf("  [ Solving with preconditioning ] -");
        sol = gmres_lu(A, L, U, b, params->tol, params->kmax[params->krilov_space_index], params->lmax);

        /* remove the preconditioners MATs */
        MATRIX_clean(L);
        MATRIX_clean(U);

        printf(" Done\n");

    }
    else
    {
        printf("  [ Solving without preconditioning ] -");
        sol = gmres_solver(A, b, params->tol, params->kmax[params->krilov_space_index], params->lmax);

        printf(" Done\n");
    }

    /* Get total time at the end of the algorithm */
    total_time = (get_time() - total_time)/100.0;
    sol.time = total_time;

    if (params->reordering)
    {
        // desfazer a permutação
        Vector x = rearrange_solution(sol.x, p);

        /* copy the x vector back to the solution, just */
        SwapVectors(sol.x, x);

        /* remover o array de permutação */
        free(p);

        /* remove the x vector */
        DeleteVector(x);

    }


    /* Print the processing time as well as the iterations number */
    /*  printTime(input, sol); */

    /* Print the solution in a file */
    printSolution(sol, output, params);
    printTime(params, sol);

    /* remove the solution internal vectors */
    delete_solution(sol);

    /* delete the allocated vectors */
    DeleteVector(b);
    DeleteVector(ones);

    return;
}

int main (int argc, char* argv[])
{
    /* create a general parameter struct */
    GMRES_ParametersPtr params = get_base_parameters();

    /* matrix and tests filenames */
    char *dubcova[3] = {"matrizes/Dubcova2.mtx", "testes/Dubcova2/Dubcova2.m", "Dubcova2" };
    char *rail[3] = {"matrizes/rail_5177.mtx", "testes/rail_5177/rail_5177.m", "rail_5177" };
    char *aft01[3] = {"matrizes/aft01.mtx", "testes/aft01/aft01.m", "aft01" };
    char *fem3d[3] = {"matrizes/FEM_3D_thermal1.mtx", "testes/FEM_3D_thermal1/FEM_3D_thermal1.m", "FEM_3D_thermal1" };


    /* build the main matrix array */
    char **matrices[4] = {dubcova, rail, aft01, fem3d};

    /*******************************************************/
    /******************* EXPERIMENTS ***********************/
    /*******************************************************/
    /* helpers */
    int i, j, k, l, m;



    /* iterate over the matrices array */
    /* TODO INSERT ALL MATRICES INSIDE THE char **matrices[] array */
    for (m = 0; m < 4; ++m)
    {
        /* open the current output test file */
        FILE *f = open_test_file(matrices[m][1]);

        /* set the main function */
        fprintf(f, "function %s()\n", matrices[m][2]);
        printf("\nMatrix: %s\n", matrices[m][2]);
        strcpy(params->name, matrices[m][2]);

        /* iterate over the parameters values */
        for (i = 0; i < 3; ++i)
        {
            /* set the current space index */
            params->krilov_space_index = i;

            /* iterate over the reordering options */
            for (j = 0; j < 2; ++j)
            {
                /* set the reordering flag */
                params->reordering = j;

                /* iterate over the preconditioning options */
                for (k = 0; k < 2; ++k)
                {
                    /* set the preconditioner index */
                    params->preconditioner = k;

                    if (false != params->preconditioner)
                    {
                        /* iterate over the ilu fill in values */
                        for (l = 0; l < 2; ++l)
                        {
                            /* set the fill in index*/
                            params->ilu = l;

                            /* open the current matrix files */
                            MAT *A = (MAT*) malloc(sizeof(MAT));
                            MATRIX_readCSR(A, matrices[m][0]);

                            /* solve with the current params */
                            /* matrix[m][0] is the matrix filename */
                            /* matrix[m][1] is the matrix test filename */
                            solver(A, f, params);

                            /* clear the current matrix */
                            MATRIX_clean(A);
                        }
                    }
                    else
                    {
                        /* open the current matrix files */
                        MAT *A = (MAT*) malloc(sizeof(MAT));
                        MATRIX_readCSR(A, matrices[m][0]);

                        /* solve with the current params */
                        /* matrix[m][0] is the current matrix filename */
                        /* matrix[m][1] is the current matrix test filename */
                        solver(A, f, params);

                        /* clear the current matrix */
                        MATRIX_clean(A);
                    }
                }

            }
        }

        /* close the test file */
        close_test_file(f);
    }

    printf("\n\nDone\n\n");

    return 0;
}
