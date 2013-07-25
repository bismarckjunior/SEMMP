#include "memory.h"

int *iVector(int size)
{
	int *v;
   
	v = (int*)calloc(size, sizeof(int));
	if (!v) eprintf("allocation failure in ivector()");
	
	return v;
}
/*****************************************************************************/

double *rVector(int size)
{
	double *v;
   
	v = (double*)calloc(size, sizeof(double));
	if (!v) eprintf("allocation failure in rvector()");
	
	return v;
}
/*****************************************************************************/

double **rMatrix(int rows, int cols)
{
	int i;
	double **m;
	
	m = (double **)calloc(rows, sizeof(double *));
	if (!m) eprintf("allocation failure 1 in rmatrix():");
	
	for (i = 0; i < rows; i++) {
		m[i] = (double *)calloc(cols, sizeof(double));
		if (!m[i]) eprintf("allocation failure 2 in rmatrix():");
	}
	
	return m;
}
/*****************************************************************************/

int **iMatrix(int rows, int cols)
{
	int i;
	int **m;
	
	m = (int **)calloc(rows, sizeof(int *));
	if (!m) eprintf("allocation failure 1 in imatrix():");
	
	for (i = 0; i < rows; i++) {
		m[i] = (int *)calloc(cols, sizeof(int));
		if (!m[i]) eprintf("allocation failure 2 in imatrix():");
	}
	
	return m;
}
/*****************************************************************************/

void freeiVector(int *v)
{
	free(v);

	return;
}
/*****************************************************************************/

void freerVector(double *v)
{
	free(v);

	return; 
}
/*****************************************************************************/

void freeiMatrix(int **m, int rows, int cols)
{
	int i;

	for(i = 0; i < rows; i++){
		free(m[i]);
	}
	free(m);

	return;
}
/*****************************************************************************/

void freerMatrix(double **m, int rows, int cols)
{
	int i;

	for(i = 0; i < rows; i++){
		free(m[i]);
	}
	free(m);

	return;
}
/*****************************************************************************/

