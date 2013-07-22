#ifndef IOSEMMP_H
#define IOSEMMP_H

#include "well.h"
#include "block.h"
#include "parameters.h"


/*ISEMMP**********************************************************************/
/*Gets table value
 * only for tables with equal step size in the independent variable: x*/
double getTableVal(double **, int, double, int, double, int);

/*Read table of double from a file*/
void rTableFile(char [], double **, int, int);

/*Read table of integer from a file*/
void iTableFile(char [], int **, int, int);

/*Reads a list of parameters from a inifile*/
void readIniFile(int, char *[], Parameters *);

/*Reads outputFile and sets data*/
Out *readAndSetOuts(Parameters*, Block*, int**, char*, double*);

/*Reads main out and sets data*/
void *readAndSetMainOut(Parameters*, double*);

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