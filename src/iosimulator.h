#ifndef IOSIMULATOR_H
#define IOSIMULATOR_H

#include <string>
//#include "iniparser.h"
//#include "block.h"
//#include "parameters.h"
//#include "well.h"
//#include "boundarycondition.h"

class out {
	bool p;			//bool for display pressure
	bool qw;		//bool for display flow rate
	bool np;		//bool for cumulative oil production	
	bool pwf;		//bool for display bottom hole pressure
};

/*Read table from a file*/
template <class T>
T readTableFile(std::string filename, int row, int col);

/*Read table of integer from a file*/
void iTableFile(char [], int **, int, int);

/*Gets table value
 * only for tables with equal step size in the independent variable: x*/
double getTableVal(double **, int, double, int, double, int);

/*Reads a list of parameters from a inifile*/
void readIniFile(int, char *[], Parameters *);

/*Reads and sets geometry*/
Block *readAndSetGeometry(Parameters *, int **);

/*Reads well parameters file and sets data*/
Well *readAndSetWellParameters(Parameters *);

/*Reads main out and sets data*/
void *readAndSetMainOut(Parameters*, double*);

/*Reads outputFile and sets data*/
Out *readAndSetOuts(Parameters*, Block*, int**, char*, double*);

/*Reads boundary conditions file and sets data*/
Boundary *readAndSetBoundaryConditions(Parameters*, int**, Block*);

/*Prints header of report and display*/
void header(FILE *, Parameters *);

/*Writes table header: Np, P, Pwf, Qw*/
void writeTableHeader(FILE*, Parameters*, Out*, int, double*);

/*Write outs to display or report*/
void writeOuts(FILE*, int, Block*, Well *, Out*, double, double*, double*);

/*Writes a file with all pressure blocks*/
void writePressureFile(Parameters *, int, Block *, double *);

/*Modifies block to ACTIVE or INACTIVE*/
void modifyGeometry(int**, Parameters*);

/*Modifies transmissibilities from a file*/
void modifyTransmissibilities(int **, Block *, Parameters *);

/*Reads the fluid properties from a file (pressure, Formation Volume Factor, 
 *viscosity, especific gravity) and stores p, 1/FVF, 1/(FVF*mu) and gamma */
double **readFluidProperties(Parameters *);

#endif 