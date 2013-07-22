#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "boundary.h"

typedef struct{
	int row;		//row of first block or unique block 
	int col;		//col of first block or unique block 
	int rowf;		//row of last block 
	int colf;		//col of last block
	int type;		//type of boundary condition
	char side;      //side: NORTH, EAST, SOUTH, WEST
	double value;	//value of boundary condition
} Boundary;

void setBoundaryConditions(int*, int*, double*, double*, double, 
						   double, double, double, Boundary*, Block*);
#endif