#ifndef IOSEMMP_H
#define IOSEMMP_H


#include <stdio.h>
#include <string>
#include <sstream>
#include "eprintf.h"
#include "well.h"
#include "block.h"
#include "parameters.h"
#include "boundary.h"
#include "memory.h"
#include <direct.h> //mkdir function
extern "C" {
#include <iniparser.h> 
} 



/* indexes for fluid properties tables */
#define PRESS			0
#define FVF				1
#define INV_FVF			1	
#define VISC			2
#define INV_FVFVISC		2
#define GAMMA			3
#define NCOLPROPS		4


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

#define ALPHAC			(double)5.614583


/* strings length */
#define SIDENAME		6		 //Side name: North, South, East, West
#define LENGTHSN		20		 //section name length
#define LENGTHFN		70		 //file name length
#define LENGTHSP		25		 //section properties length
#define LINESIZE		1313	 //0.0.2 max line size in report (MAXCOLREPORT)
#define LENBCTYPE		30		 //boundary condition type
#define LENWELLTYPE		20		 //see: /* well types */ 
#define LENBLOCKTYPE	9		 //block type: ACTIVE, INACTIVE

/*max column*/
#define MAXCOLDISPLAY	5   
#define MAXCOLREPORT	100  

/* extensions for output files*/
#define REPORTEXT		".txt"   //report extension
#define OUTPRESSUREEXT	".vtk"   //output pressure extension



/*ISEMMP**********************************************************************/

//char *createKey(std::string k, int num, std::string secName);


/*Gets table value
 * only for tables with equal step size in the independent variable: x*/
double getTableVal(double **, int, double, int, double, int);


/** 
 * Read table from a file.
 * @param filein filename
 * @param nrow number of rows
 * @param ncol number of columms
 * @return table table read
 */
template <class T>
T **readTableFile(char filein[], int nrow, int ncol);

/*Reads a list of parameters from a inifile*/
void readIniFile(int, char *[], Parameters *);

/*Reads outputFile and sets data*/
Out *readAndSetOuts(Parameters*, Block*, int**, char*, double*);

/*Reads main out and sets data*/
void readAndSetMainOut(Parameters*, double*);

/*Reads well parameters file and sets data*/
Well *readAndSetWellParameters(Parameters *);

/*Reads and sets geometry*/
Block *readAndSetGeometry(Parameters *, int **);

/*Reads the fluid properties from a file (pressure, Formation Volume Factor, 
 *viscosity, especific gravity) and stores p, 1/FVF, 1/(FVF*mu) and gamma */
double **readFluidProperties(Parameters *);

/*Reads boundary conditions file and sets data*/
Boundary *readAndSetBoundaryConditions(Parameters*, int**, Block*);

/*Modifies block to ACTIVE or INACTIVE*/
void modifyGeometry(int**, Parameters*);

/*Modifies transmissibilities from a file*/
void modifyTransmissibilities(int **, Block *, Parameters *);


/*****************************************************************************/

/*OSEMMP**********************************************************************/
/**
 * Writes outs to display or report
 *
 * @param f file
 * @return void
 */
void writeOuts(FILE*, int, Block*, Well *, Out*, double, double*, double*);

/*Writes table header*/
void writeTableHeader(FILE*, Parameters*, Out*, int, double*);

/*Prints report header and display*/
void header(FILE *, Parameters *);

/*Writes a file with all pressure blocks*/
void writePressureFile(Parameters *, int, Block *, double *);

/*****************************************************************************/

#endif