#ifndef PARAMETERS_H
#define PARAMETERS_H

#define LENGTHFN	70		//file name length

typedef struct {
	double iniTime;			//initial time 
	double cf;				//fluid compressibility [1/psi]
	double cphi;			//rock compressibility [1/psi]
	double dpprops;			//step size in pressure [psia]
	double p0;				//reference pressure [psia]
	double z0;				//detph for reference pressure [ft]
	double invBSC;			//multiplicative inverse of BSC 
	double rhoSC;			//density in standart conditions [lbm/ft³]
	double dt;				//time step
	double dtMultiplier;	//time step multiplier
	
	int dtLogScale;			//bool for log scale in time
	int ncol;				//number of columns
	int nrow;				//number of rows
	int nSteps;				//number of steps
	int outPressureSteps;   //steps for print pressures in outPressureFile
	int nBlocks;			//number of blocks
	int nfprop;				//number of lines in fluid table
	int nwells;				//number of wells
	int nMatrix;			//number of elements in main matrix
	int displaySteps;		//steps for display on terminal
	int reportSteps;		//steps for print report of blocks (in reportFile)
	int nDisplayBlocks;		//number of blocks to display on terminal
	int nReportBlocks;		//number of blocks to the report (in outputFile)
	int nBoundary;          //number of boundary contitions (in bcFile)
	int isCylindrical;		//bool for cylindrical coordenates

	char reportFile[LENGTHFN];			//report file
	char outPressureFile[LENGTHFN];		//out pressure file
	char projectDir[LENGTHFN];			//project path
	char projectName[LENGTHFN];			//project name
	char iniFile[LENGTHFN];				//inifile
	char geometryFile[LENGTHFN];		//geometry file
	char modifiedBlocksFile[LENGTHFN]; 
	char modTransmissibilityFile[LENGTHFN]; 
	char porosityFile[LENGTHFN];		//porosity file
	char kxFile[LENGTHFN];				//kx file
	char kyFile[LENGTHFN];				//ky file
	char thicknessFile[LENGTHFN];		//thickness file
	char zTopFile[LENGTHFN];			//block top depth file
	char dxFile[LENGTHFN];				//dx file
	char dyFile[LENGTHFN];				//dy file
	char fluidPropFile[LENGTHFN];		//fluid properties file
	char wellsFile[LENGTHFN];			//wells file
	char outputFile[LENGTHFN];			//output file
	char bcFile[LENGTHFN];				//boundary conditions file
} Parameters;

#endif