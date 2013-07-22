#ifndef GRID_H
#define GRID_H

#include "block.h"

void setMatrixStructure(Parameters *, Block *, int *, int *, double *,int *, 
						int *, double *, int *);

void setTransmissibilityMatrix(Parameters *, Block *, Well *, double **, 
							   double *, double *, double *, double *, int *, 
							   Boundary*);

void setInitialPressure(Parameters *, Block *, Well *, double **, double *, 
						double *, double *);

void setCartesianTransmissibilities(Parameters *, Block *, Boundary *); 

void setCylindricalTransmissibilities(Parameters *, Block *, Boundary *, 
									  Well *);

#endif