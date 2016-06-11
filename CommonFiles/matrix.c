/*----------------------------------------------------------------------------
 * MATRIX FUNCTIONS: IN CSR FORMAT
 *--------------------------------------------------------------------------*/
#include "protos.h"


int COMPARE_array (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr3 <  ((ARRAY*)b)->arr3) return -1;
	if (((ARRAY*)a)->arr3 >  ((ARRAY*)b)->arr3) return  1;
	if (((ARRAY*)a)->arr3 == ((ARRAY*)b)->arr3)
	{
		if (((ARRAY*)a)->arr2 < ((ARRAY*)b)->arr2) return -1;
		if (((ARRAY*)a)->arr2 > ((ARRAY*)b)->arr2) return  1;
	}
	return 0;
}

/*----------------------------------------------------------------------------
 * Read matrix from file in MM format to a CSR structure
 *--------------------------------------------------------------------------*/
void MATRIX_readCSR (MAT* A, char* p)
{
	int M, N, nz;   
	int i, j, row, colunm, elem = 0;
	double value;
	char line[1025];
	FILE* f;

	if ((f = fopen(p, "r")) == NULL) 
		exit(1);

	do 
	{
		if (fgets(line,1025,f) == NULL) 
			exit(0);
	}while (line[0] == '%');
   
	sscanf(line,"%d %d %d", &N, &M, &nz);

	/* reseve memory for matrices */
	A->m   = M;
	A->n   = N;
	A->nz  = nz;
    
	A->AA  = (double *) calloc(nz, sizeof (double));
	A->D   = (double *) calloc(N,  sizeof (double));
	A->JA  = (int    *) calloc(nz, sizeof (int));
	A->IA  = (int    *) calloc(N+1,sizeof (int));
	
	for (i = 0; i < nz; ++i)
	{
		fscanf(f, "%d %d %lf\n", &colunm, &row, &value);
		A->AA[i]   = value;
		A->JA[i]   = colunm - 1;
		elem      += 1;
		A->IA[row] = elem;
	}
	
	/* Adjusting IA array */
	for (i = 1; i < N + 1; ++i)
		if (A->IA[i] == 0) 
			A->IA[i] = A->IA[i-1];
	
	/* Diagonal */
	if (M == N) /* square matrix */
	{
		for (i = 0; i < A->n; ++i)
		{
			int k1 = A->IA[i];
			int k2 = A->IA[i+1]-1;
			for (j = k1; j <= k2; ++j)
				if (i == A->JA[j]) 
					A->D[i] = A->AA[j];
		}
	}
	fclose(f);
}

/*----------------------------------------------------------------------------
 * Compute the matrix bandwidth
 *--------------------------------------------------------------------------*/
long int MATRIX_bandwidth (MAT* A)
{
	int i, min;
	int n = A->n;
	unsigned long int bandwidth=0;
	
	for (i = 0; i < n; ++i)
	{
		min = (i - A->JA[A->IA[i]]);
		if (bandwidth < min) bandwidth = min; 
	}
	
	return bandwidth;
}

/*----------------------------------------------------------------------------
 * Perform the operation P*A*P' in CSR format
 *--------------------------------------------------------------------------*/
void MATRIX_permutation (MAT* A, int* p)
{

	int i, j, k;
	int n   = A->n;
	int nz  = A->nz;  

	MAT*  B     = malloc(      sizeof(MAT));
	      B->AA = malloc( nz  *sizeof(double));
	      B->D  = malloc( n   *sizeof(double));
	      B->JA = malloc( nz  *sizeof(int));
	      B->IA = malloc((n+1)*sizeof(int));
  
	B->IA[0] = 0;
	B->n     = n;
	B->m     = n;
	B->nz    = nz;

	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = A->IA[p[i]]; j <= A->IA[p[i]+1] - 1; ++j)
		{
			B->AA[k] = A->AA[j];
			B->JA[k] = A->JA[j];
			      k  = k + 1;
		}
		B->IA[i+1] = k;    
	}
	
	ARRAY* a = malloc (nz*sizeof(ARRAY));
	int*   q = malloc (n *sizeof(int));
	
	for (i = 0; i < n; ++i) 
		q[p[i]] = i; 
	
	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = B->IA[i]; j <= B->IA[i+1] - 1; ++j)
		{
			a[k].arr1 = B->AA[j];
			a[k].arr2 = q[B->JA[j]];
			a[k].arr3 = i;
				k = k + 1;
		}
		A->IA[i+1] = k;    
	}
	
	qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	
	for (i = 0; i < nz; ++i)
	{
		A->AA[i] = a[i].arr1;
		A->JA[i] = a[i].arr2;
	}
	
	free(a);
	free(q);
	MATRIX_clean(B);
}

/*----------------------------------------------------------------------------
 * Print L and U matrices
 *--------------------------------------------------------------------------*/
void MATRIX_printLU (MAT* A, MAT* L, MAT* U)
{
	int i, j, k, c;
	int n = L->n;
	printf ("\n\n");
	for (i = 0; i < n; ++i)
	{
		for (k = 0; k < n; ++k)
		{
		        c = 0;
			for (j = L->IA[i]; j < L->IA[i+1]; ++j)
				if (k == L->JA[j]){ printf ("%6.2lf ",L->AA[j]); c = 1; }
			if  (i == k) { printf ("%6.2lf ",U->D[i]); c = 1; }
 			for (j = U->IA[i]; j < U->IA[i+1]; ++j)
				if (k == U->JA[j]){ printf ("%6.2lf ",U->AA[j]); c = 1; }
			if (c == 0) printf ("%6.2lf ",0.0);
		}
		printf ("\n");
	}
	printf ("\n");
}

/*----------------------------------------------------------------------------
 * Free up memory allocated for MAT struct
 *--------------------------------------------------------------------------*/
void MATRIX_clean (MAT* A)
{
	free(A->AA);
	free(A->JA);
	free(A->IA);
	free(A->D);
	free(A);  
}

