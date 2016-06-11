//Multiplicacao Matriz x Vetor em CSR, onde MAT* a é a matriz CSR e Mat* b é o vetor.
double* multiplicacaoCSR(MAT* a, double* b)
{
	int i, j, linha_inicio, linha_fim;
	double *p = (double*) malloc (sizeof(double));

	for(i = 0; i < n; i++)
	{
		p[i] = 0;
		
		linha_inicio = a->IA[i];
		linha_fim = a->IA[i];


		for (j = linha_inicio; j < linha_fim; j++)
		{
			p[i] = p[i] + a->AA[j] * b(a->JA[j]);
		}

	}

	return p;

}