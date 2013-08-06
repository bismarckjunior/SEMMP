#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "block.h"
#include "eprintf.h"


/* boundary conditions */
#define PRESSURE_GRADIENT_ESPECIFIED  1 
#define PRESSURE_ESPECIFIED	        2  //double definition (well.h)
#define NORTH			'N'
#define EAST			'E'
#define WEST			'W'
#define SOUTH			'S'
#define ABOVE			'A'
#define BELOW			'B'

/* bools */
#define TRUE			1   
#define FALSE			0

typedef struct{
	int row;		///< first block row or unique block 
	int col;		///< first block column or unique block 
	int lay;		///< first block layer or unique block
	int rowf;		///< last block row
	int colf;		///< last block column
	int layf;		///< last block layer
	int type;		///< type of boundary condition
	char side;      ///< side: NORTH, EAST, SOUTH, WEST
	double value;	///< value of boundary condition
} Boundary;


void setBoundaryConditions(int *j, int *Map, double *Ax, double *b, 
									 double N, double E, double W, double S, 
									 double A, double B,
									 Boundary* boundary, Block *block);
#endif