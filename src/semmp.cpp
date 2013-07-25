/* ------------------------------------------------------------------------- */
/*                                                                           */
/*        SINGLE PHASE - COMPRESSIBLE AND SLIGHTLY-COMPRESSIBLE FLOW         */
/*     CARTESIAN (variable depth/variable thickness) - FINITE DIFFERENCE     */
/*                           BLOCK CENTERED GRID                             */
/*                           SIMULATOR Ver. 0.0.2                            */
/*                                                                           */
/*                           capico@lenep.uenf.br                            */
/*                         bismarckgomes@google.com                          */
/*                                                                           */
/*                                                                           */
/* ------------------------------------------------------------------------- */
/* REFERENCES:                                                               */
/*                                                                           */
/* Ertekin T., Abou-Kassem J. H., King G. R., (2001), Basic Applied Reservoir*/
/* Simulation, Society of Petroleum Engineers, Richardson, TX.               */
/*                                                                           */
/* Aziz K., Settari A., (1979), Petroleum Reservoir Simulation, Applied      */ 
/* Science Publishers, London.                                               */
/* ------------------------------------------------------------------------- */
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
//#include <malloc.h>
//#include <limits.h>
//#include <float.h>
//#include "umfpack.h"
//#include "iniparser.h"
//#include "eprintf.h"
#include "semmp.h"

/* ------------------------------------------------------------------------- */
/*                                   MAIN                                    */
/* ------------------------------------------------------------------------- */

int main (int argc, char *argv[]) 
{
	int i, step = 0, addSteps = 0;
	double **fprops; 
	double *pressureN, *pressure, *pressureNp1, *ptmp;
	double t, dt, delta = 0.0;
	double cumulativeProd = NONP;
	
	/* vetors for triplets format for sparse matrix */
	int *Ti, *Tj;
	double *Tx;
	
	/* vetors for compressed column format for sparse matrix */
	int *Ai, *Ap;
	double *Ax;

	double *b;
	int *Map;
	double *null = (double *)NULL ;
	double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
	void *Symbolic, *Numeric ;
	int status;
	int **geom; /* temporary array for geometry */
	Block *grid;
	Well *wells;
	Parameters par;
	Out *report;
	Out *display;
	FILE *reportFile;
	Boundary *boundary;
	/*************************************************************************/

	readIniFile(argc, argv, &par);
	
	header(stdout, &par); 
	
	geom  = iMatrix(par.nrow, par.ncol);

	wells = readAndSetWellParameters(&par);

	grid  = readAndSetGeometry(&par, geom);

	boundary = readAndSetBoundaryConditions(&par, geom, grid); 

	if (!par.isCylindrical) setWells(&par, grid, geom, wells);

	readAndSetMainOut(&par, &cumulativeProd);

	report =  readAndSetOuts(&par, grid, geom, "report", &cumulativeProd);

	display = readAndSetOuts(&par, grid, geom, "display", &cumulativeProd);

	freeiMatrix(geom, par.nrow, par.ncol);

	if (par.isCylindrical) //0.0.3
		setCylindricalTransmissibilities(&par, grid, boundary, wells);
	else
		setCartesianTransmissibilities(&par, grid, boundary); 
	
	fprops = readFluidProperties(&par);

	pressureN   = rVector(par.nBlocks); /* pressure at t */
	pressure    = rVector(par.nBlocks); /* pressure for iterative solution */
	pressureNp1 = rVector(par.nBlocks); /* pressure at t + delta_t */

	setInitialPressure(&par, grid, wells, fprops, pressureN, pressure, 
		pressureNp1);

	/* arrays for matrix structure, triplet and compressed column forms */
	Ap  = iVector(par.nBlocks + 1);
	Ai  = iVector(par.nMatrix);
	Ax  = rVector(par.nMatrix);
	b   = rVector(par.nBlocks);
	Ti  = iVector(par.nMatrix);
	Tj  = iVector(par.nMatrix);
	Tx  = rVector(par.nMatrix);
	Map = iVector(par.nMatrix);
	
	setMatrixStructure(&par, grid, Ti, Tj, Tx, Ap, Ai, Ax, Map);
	freeiVector(Ti);
	freeiVector(Tj);
	freerVector(Tx);
	
	umfpack_di_defaults(Control);
	Control[UMFPACK_PRL] = 5.0;
    status = umfpack_di_symbolic(par.nBlocks, par.nBlocks, Ap, Ai, Ax, 
		&Symbolic, Control, Info);
	
	setTransmissibilityMatrix(&par, grid, wells, fprops, pressureN,
							  pressure, Ax, b, Map, boundary);
	printf("\nPress Enter to continue ...\n");
	getchar();

	/**************  Table Header  **************/
	writeTableHeader(stdout, &par, display, par.nDisplayBlocks, 
					 &cumulativeProd);
	if (par.reportSteps != 0 || cumulativeProd != NONP){
		reportFile = fopen(par.reportFile, "w");
		header(reportFile, &par);
		writeTableHeader(reportFile, &par, report, par.nReportBlocks, &cumulativeProd);
	}
	/********************************************/
	
	/************* print inicial pressure**********************/
	if(par.outPressureSteps != 0)
		writePressureFile(&par, step, grid, pressureN);
	
	if(par.reportSteps != 0 || cumulativeProd != NONP)
		writeOuts(reportFile, par.nReportBlocks, grid, wells, report, 
			par.iniTime, pressureN, &cumulativeProd);

	if(par.displaySteps != 0 || cumulativeProd != NONP)
		writeOuts(stdout, par.nDisplayBlocks, grid, wells, display, 
			par.iniTime, pressureN, &cumulativeProd); 
	/*******************************************************/


	/*****************************************************************/
	/* ************************* TIME LOOP ************************* */
	
	dt = par.dt;
	par.dt = par.iniTime + par.dt;
	t = par.iniTime;
	while(step < par.nSteps){
		if (par.dtLogScale){
			par.dt = t*(pow(10.0, dt)-1);
		}
		
		do{
			ptmp = pressure;
			pressure = pressureNp1;
			pressureNp1 = ptmp;
			
			setTransmissibilityMatrix(&par, grid, wells, fprops, pressureN,
				pressure, Ax, b, Map, boundary);
			status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, 
				Control, Info);
			status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, pressureNp1, b, 
				Numeric, Control, Info);
			umfpack_di_free_numeric(&Numeric);

			delta = normDelta(pressure, pressureNp1, par.nBlocks);
		}
		while(delta > EPSILON);
		
		/*IMB
		if (par.nwells == 0){
			IMB = 0;
			cumulativeProd = 1;
			for (i=0; i<par.nBlocks; i++){
				IMB += grid[i].Vac*grid[i].phi*(
					getTableVal(fprops, PRESS, pressureNp1[i], INV_FVFVISC, par.dpprops, par.nfprop)
					-getTableVal(fprops, PRESS, pressureN[i], INV_FVFVISC, par.dpprops, par.nfprop));
			}
		}//*/

		ptmp = pressureN;
		pressureN = pressureNp1;
		pressureNp1 = ptmp;

		if (step==0 && !par.dtLogScale) par.dt = dt;

		t += par.dt;

		par.dt *= par.dtMultiplier;

		step++;
		
		/******* Cumulative Production ***********/
		if (cumulativeProd != NONP){
			for(i=0; i < par.nwells; i++)
				if(wells[i].qw < 0){
					wells[i].np += -wells[i].qw*par.dt;
					cumulativeProd += -wells[i].qw*par.dt; 
				}
		}
		/*****************************************/
		
		/************* outputs *********************************/
		if(par.outPressureSteps != 0 && step%par.outPressureSteps == 0)
			writePressureFile(&par, step, grid, pressureN);
		
		if(par.reportSteps != 0 && step%par.reportSteps == 0 
												|| cumulativeProd != NONP)
			writeOuts(reportFile, par.nReportBlocks, grid, wells, report, t, 
				pressureN, &cumulativeProd); //&IMB

		if(par.displaySteps != 0 && step%par.displaySteps == 0
												|| cumulativeProd != NONP)
			writeOuts(stdout, par.nDisplayBlocks, grid, wells, display, t, 
				pressureN, &cumulativeProd); //&IMB
		/*******************************************************/

		if(step == par.nSteps){
			printf("\a\n Steps finished. How much more steps ? ");
			scanf("%d", &addSteps);
			par.nSteps += addSteps;
			if (addSteps != 0)
				writeTableHeader(stdout, &par, display, par.nDisplayBlocks,
								 &cumulativeProd);
		}

	}

	if (par.reportSteps) fclose(reportFile);

	umfpack_di_free_symbolic(&Symbolic);
	free(grid);
	free(fprops);
	freerVector(pressureN);
	freerVector(pressure);
	freerVector(pressureNp1);
	freeiVector(Ap);
	freeiVector(Ai);
	freerVector(Ax);
	freerVector(b);
	freeiVector(Map);
       
	return 0;
}
/* ------------------------------------------------------------------------- */
/*                              END  MAIN                                    */
/* ------------------------------------------------------------------------- */









