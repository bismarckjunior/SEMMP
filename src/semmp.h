#include "block.h"
#include "well.h"
#include "boundary.h"
#include "parameters.h"
#include "iosemmp.h"
#include "grid.h"


/* indexes for fluid properties tables */
#define PRESS			0
#define FVF				1
#define INV_FVF			1	
#define VISC			2
#define INV_FVFVISC		2
#define GAMMA			3
#define NCOLPROPS		4

/*max column*/
#define MAXCOLDISPLAY	5   
#define MAXCOLREPORT	100  

/* strings length */
#define SIDENAME		6		 //Side name: North, South, East, West
#define LENGTHSN		20		 //section name length
#define LENGTHFN		70		 //file name length
#define LENGTHSP		25		 //section properties length
#define LINESIZE		1313	 //0.0.2 max line size in report (MAXCOLREPORT)
#define LENBCTYPE		30		 //boundary condition type
#define LENWELLTYPE		20		 //see: /* well types */ 
#define LENBLOCKTYPE	9		 //block type: ACTIVE, INACTIVE

/* extensions for output files*/
#define REPORTEXT		".txt"   //report extension
#define OUTPRESSUREEXT	".vtk"   //output pressure extension

/* default values */
#define NONP			-1		 //do not print Np
#define NOWELL			-1		 //block does not have a well
#define ACTIVEBLOCK		+1		 //active block
#define INVALIDBLOCK	 0		 //invalid enter
#define INACTIVEBLOCK	-1		 //inactive block
#define NOBOUNDARY		-1		 //block does not have boundary condition
#define UNDEF       "**UNDEF**"  //file name not found
//#define CALCULATEG      -1	     //flag to calculate G            
#define ISBONDARY       -3		 //flag to verify if face is bondary
#define NOTRANSMISSIBILITY -1

/* boundary conditions */
#define PRESSURE_GRADIENT_ESPECIFIED  1 
//#define PRESSURE_ESPECIFIED	        2
#define NORTH			'N'
#define EAST			'E'
#define WEST			'W'
#define SOUTH			'S'

/* bools */
#define TRUE			1   
#define FALSE			0

/* well types */
#define RATE_ESPECIFIED		1    //well type
#define PRESSURE_ESPECIFIED	2	 //well type
#define DING				0.25 //Ding's coupling model

/* constants */
#define BETAC			(double)1.127
#define ALPHAC			(double)5.614583
#define GAMMAC			(double)0.21584e-3
#define GRAV			(double)32.174
#define EPSILON			(double)1e-8
#define PI				(double)3.14159265358979

#define BARREL_TO_CF	(double)5.614584


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

double normDelta(double *, double *, int);