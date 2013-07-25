#ifndef GRID_H
#define GRID_H

#include <math.h>
#include "block.h"
#include "well.h"
#include "iosemmp.h"
#include "eprintf.h"
#include "parameters.h"
#include "umfpack.h"
#include "boundary.h"
#include "memory.h"

#define EPSILON			(double)1e-8



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

void setWells(Parameters *, Block *, int **, Well *);

double normDelta(double *, double *, int );

#endif