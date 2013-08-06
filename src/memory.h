#ifndef MEMORY_H
#define MEMORY_H

#include <malloc.h>
#include "eprintf.h"

/******* memory allocation utilities ********/
int *iVector(int);
int **iMatrix(int, int);
void freeiVector(int *);
void freeiMatrix(int **, int, int);
double *rVector(int);
double **rMatrix(int, int);
void freerVector(double *);
void freerMatrix(double **, int, int);
/********************************************/


template <class T>
int*** array3D(int , int , int );



///template <class T>
void freeArray(int***, int rows, int cols, int lays);

#endif