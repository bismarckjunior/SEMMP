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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <limits.h>
#include <float.h>
#include "umfpack.h"
#include "iniparser.h"
#include "eprintf.h"
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
	FILE *reportFile, *f;
	Boundary *boundary;
	/*************************************************************************/

	readIniFile(argc, argv, &par);
	
	header(stdout, &par); 
	
	geom  = iMatrix(par.nrow, par.ncol, par.nlay);

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

int *iVector(int size)
{
	int *v;
   
	v = (int*)calloc(size, sizeof(int));
	if (!v) eprintf("allocation failure in ivector()");
	
	return v;
}
/*****************************************************************************/

double *rVector(int size)
{
	double *v;
   
	v = (double*)calloc(size, sizeof(double));
	if (!v) eprintf("allocation failure in rvector()");
	
	return v;
}
/*****************************************************************************/

double **rMatrix(int rows, int cols)
{
	int i;
	double **m;
	
	m = (double **)calloc(rows, sizeof(double *));
	if (!m) eprintf("allocation failure 1 in rmatrix():");
	
	for (i = 0; i < rows; i++) {
		m[i] = (double *)calloc(cols, sizeof(double));
		if (!m[i]) eprintf("allocation failure 2 in rmatrix():");
	}
	
	return m;
}
/*****************************************************************************/

int **iMatrix(int rows, int cols)
{
	int i;
	int **m;
	
	m = (int **)calloc(rows, sizeof(int *));
	if (!m) eprintf("allocation failure 1 in imatrix():");
	
	for (i = 0; i < rows; i++) {
		m[i] = (int *)calloc(cols, sizeof(int));
		if (!m[i]) eprintf("allocation failure 2 in imatrix():");
	}
	
	return m;
}
/*****************************************************************************/



void freeiVector(int *v)
{
	free(v);

	return;
}
/*****************************************************************************/

void freerVector(double *v)
{
	free(v);

	return;
}
/*****************************************************************************/

void freeiMatrix(int **m, int rows, int cols)
{
	int i;

	for(i = 0; i < rows; i++){
		free(m[i]);
	}
	free(m);

	return;
}
/*****************************************************************************/

void freerMatrix(double **m, int rows, int cols)
{
	int i;

	for(i = 0; i < rows; i++){
		free(m[i]);
	}
	free(m);

	return;
}
/*****************************************************************************/

void rTableFile(char filein[], double **table, int nrow, int ncol)
{
	float tmp;
	int col, row;
	FILE *mf;
		
	if ( (mf = fopen(filein, "r")) == NULL ){
		eprintf("opening table file %s:", filein);
	}
	
	for (row = 0; row < nrow; row++) {
		for (col = 0; col < ncol; col++) {
			
			fscanf(mf, "%f", &tmp);
			table[row][col] = (double)tmp;
		}
	}
	fclose(mf);
	
	return;
}
/*****************************************************************************/

void iTableFile(char filein[], int **table, int nrow, int ncol)
{
	float tmp;
	int col, row;
	FILE *mf;
		
	if ( (mf = fopen(filein, "r")) == NULL ){
		eprintf("opening table file:");
	}
	
	for (row = 0; row < nrow; row++) {
		for (col = 0; col < ncol; col++) {
			fscanf(mf, "%f", &tmp);
			table[row][col] = (int)tmp;
		}
	}
	fclose(mf);
	
	return;
}
/*****************************************************************************/

Block *readAndSetGeometry(Parameters *par, int **geom)
{
	int nrow, ncol, nlay, col, row, lay, localIndex, aux, i = 0, j = 0, k = 0;
	int isGeoFile, isPhiFile, rowf, colf, line = 0;
	int isKxFile, isKyFile, isDxFile, isDyFile, isDzFile, isZtopFile;
	double phi, kx, ky, dx, dy, dz, ztop; 
	double **porositytmp;
	double **dytmp, **kxtmp, **kytmp, **thicknesstmp, **zToptmp, **dxtmp;
	char dirFileName[LENGTHFN];
	Block *grid;
	FILE *mf;	
	
	setprogname("set geometry");
	
	nrow = par->nrow;
	ncol = par->ncol;
	nlay = par->nlay;
	par->nBlocks = 0;	

	isGeoFile = (strcmp(par->geometryFile, "full") != 0); 

	if (isGeoFile) {
		sprintf(dirFileName, "%s%s", par->projectDir, par->geometryFile);
		if ( (mf = fopen(dirFileName, "r")) == NULL ){ 
			eprintf("opening geometry file:");
		}
	}

	/* temporary arrays for blocks properties */
	if (isPhiFile = !atof(par->porosityFile)){
		porositytmp = rMatrix(nrow, ncol); 
		sprintf(dirFileName, "%s%s", par->projectDir, par->porosityFile);
		rTableFile(dirFileName, porositytmp, nrow, ncol);
	}
	else {
		phi = atof(par->porosityFile); 
	}
	if (isKxFile = !atof(par->kxFile)){
		kxtmp = rMatrix(nrow, ncol);
		sprintf(dirFileName, "%s%s", par->projectDir, par->kxFile);
		rTableFile(dirFileName, kxtmp, nrow, ncol);
	}
	else {
		kx = atof(par->kxFile);
	}
	if (isKyFile = !atof(par->kyFile)){
		kytmp = rMatrix(nrow, ncol);
		sprintf(dirFileName, "%s%s", par->projectDir, par->kyFile);
		rTableFile(dirFileName, kytmp, nrow, ncol);
	}
	else {
		ky = atof(par->kyFile);
	}
	if (isDzFile = !atof(par->thicknessFile)){
		thicknesstmp = rMatrix(nrow, ncol);
		sprintf(dirFileName, "%s%s", par->projectDir, par->thicknessFile);
		rTableFile(dirFileName, thicknesstmp, nrow, ncol);
	}
	else {
		dz = atof(par->thicknessFile);
	}
	if (isZtopFile = !atof(par->zTopFile)){
		zToptmp = rMatrix(nrow, ncol);
		sprintf(dirFileName, "%s%s", par->projectDir, par->zTopFile);
		rTableFile(dirFileName, zToptmp, nrow, ncol);
	}
	else {
		ztop = atof(par->zTopFile);
	}
	if (isDxFile = !atof(par->dxFile)){
		dxtmp = rMatrix(1, ncol);
		sprintf(dirFileName, "%s%s", par->projectDir, par->dxFile);
		rTableFile(dirFileName, dxtmp, 1, ncol);
	}
	else {
		dx = atof(par->dxFile);
	}
	if (isDyFile = !atof(par->dyFile)){
		dytmp = rMatrix(nrow, 1);
		sprintf(dirFileName, "%s%s", par->projectDir, par->dyFile);
		rTableFile(dirFileName, dytmp, nrow, 1);
	}
	else {
		dy = atof(par->dyFile);
	}
	/**************************************************/		 

	if (strcmp(par->modifiedBlocksFile, UNDEF) != 0){
		modifyGeometry(geom, par);
		setprogname("set geometry"); 
	}

	aux = 1;
	for (row = 0; row < nrow; row++) {
		for (col = 0; col < ncol; col++) {
		
			if (isGeoFile) 
				fscanf(mf, "%d", &aux);
			if(aux > 0 && geom[row][col] != INACTIVEBLOCK) { 
				++(par->nBlocks);
				geom[row][col] = (par->nBlocks);
			}
		}
	}
	if (isGeoFile)	fclose(mf);

	/* number of non zero entries on the the matrix for the system of linear */
	/* equations, T*p = b */  
	par->nMatrix = par->nBlocks*5;

	/* array for the geometric properties structures */
	grid = (Block *)calloc(par->nBlocks, sizeof(Block));
	if (!grid) eprintf("allocation failure for grid blocks:");
	
	
	/* set blocks data */
	for (row = 0; row < nrow; row++) {
		for (col = 0; col < ncol; col++) {
			if(geom[row][col] > 0) { 
				localIndex = (geom[row][col] - 1);

				grid[localIndex].row = row;
				grid[localIndex].col = col;
				grid[localIndex].dx = isDxFile ? dxtmp[0][col] : dx;
				grid[localIndex].dy = isDyFile ? dytmp[row][0] : dy;
				grid[localIndex].dz = isDzFile ? thicknesstmp[row][col] : dz;
				grid[localIndex].z  = (isZtopFile ? zToptmp[row][col] : ztop)
					+ 0.5*(isDzFile ? thicknesstmp[row][col] : dz);
				grid[localIndex].Vac = 
					(grid[localIndex].dx) * 
					(grid[localIndex].dy) * 
					(grid[localIndex].dz) / ALPHAC;
				grid[localIndex].phi = isPhiFile ? porositytmp[row][col] : phi;
				grid[localIndex].kx = (isKxFile ? kxtmp[row][col] : kx)/1000.0; 
				grid[localIndex].ky = (isKyFile ? kytmp[row][col] : ky)/1000.0; 
				
				grid[localIndex].H = localIndex;
				grid[localIndex].E = 
					geom[(row+nrow)%nrow][(col+ncol+1)%ncol] - 1;
				grid[localIndex].W = 
					geom[(row+nrow)%nrow][(col+ncol-1)%ncol] - 1;
				grid[localIndex].N = 
					geom[(row+nrow+1)%nrow][(col+ncol)%ncol] - 1;
				grid[localIndex].S = 
					geom[(row+nrow-1)%nrow][(col+ncol)%ncol] - 1;

				grid[localIndex].isWellBlock = NOWELL; 

				grid[localIndex].isBoundary = NOBOUNDARY;
				grid[localIndex].nBoundary = 0;

			}
		}
	}

	if(strcmp(par->modTransmissibilityFile, UNDEF)!=0){
		modifyTransmissibilities(geom, grid, par);
		setprogname("set geometry");
	}
	

	if(isKxFile)	freerMatrix(kxtmp, nrow, ncol);
	if(isKyFile)	freerMatrix(kytmp, nrow, ncol);
	if(isDzFile)	freerMatrix(thicknesstmp, nrow, ncol);
	if(isDxFile)	freerMatrix(dxtmp, 1, ncol);
	if(isDyFile)	freerMatrix(dytmp, nrow, 1);
	if(isPhiFile)	freerMatrix(porositytmp, nrow, ncol);
	if(isZtopFile)	freerMatrix(zToptmp, nrow, ncol);
	
	unsetprogname();
		
	return grid;
	
}
/*****************************************************************************/

void setCartesianTransmissibilities(Parameters *par, Block *grid, Boundary *boundary)
{
	int i, j, bcIndex;
	int bool_E, bool_W, bool_N, bool_S;
	double dxi, dyi, dzi, kxi, kyi;
	double dxip, dzip, kxip, dxim, dzim, kxim; 
	double dyjp, dzjp, kyjp, dyjm, dzjm, kyjm;

	setprogname("set transmissibilities");

	/****** constant part of the transmissibilities ******/
	for(i = 0; i < par->nBlocks; i++){

		/*********** for nonflow boundary conditions *********/
		if(grid[i].E < 0){ 
			if (grid[i].E != ISBONDARY)
				grid[i].Gxph = NOTRANSMISSIBILITY;
			grid[i].E = grid[i].H;			
			par->nMatrix--;
		}
		
		if(grid[i].W < 0){
			if (grid[i].W != ISBONDARY)
				grid[i].Gxmh = NOTRANSMISSIBILITY;
			grid[i].W = grid[i].H;
			par->nMatrix--;
		}
		
		if(grid[i].N < 0){
			if (grid[i].N != ISBONDARY)
				grid[i].Gyph = NOTRANSMISSIBILITY;
			grid[i].N = grid[i].H;
			par->nMatrix--;
		}
		
		if(grid[i].S < 0){
			if (grid[i].S != ISBONDARY)
				grid[i].Gymh = NOTRANSMISSIBILITY;
			grid[i].S = grid[i].H;
			par->nMatrix--;
		}
		/*****************************************************/

		dxi = grid[i].dx;
		dyi = grid[i].dy;
		dzi = grid[i].dz;
		kxi = grid[i].kx;
		kyi = grid[i].ky;

		bool_E = bool_W = bool_N = bool_S = FALSE;

		for(j=0; j < grid[i].nBoundary; j++){
			bcIndex = grid[i].boundary[j];
			bool_E = boundary[bcIndex].side == EAST  ? TRUE : bool_E;
			bool_W = boundary[bcIndex].side == WEST  ? TRUE : bool_W;
			bool_N = boundary[bcIndex].side == NORTH ? TRUE : bool_N;
			bool_S = boundary[bcIndex].side == SOUTH ? TRUE : bool_S;
		}
		
		if (grid[i].Gxph != NOTRANSMISSIBILITY || bool_E){
			dxip = grid[grid[i].E].dx;
			dzip = grid[grid[i].E].dz;
			kxip = grid[grid[i].E].kx;

			grid[i].Gxph = (2.0*BETAC*dyi*dzi*dyi*dzip*kxi*kxip) 
				/ (kxi*dyi*dzi*dxip + kxip*dyi*dzip*dxi);
		}
		else{
			grid[i].Gxph = 0.0;
		}

		
		if (grid[i].Gxmh != NOTRANSMISSIBILITY || bool_W){
			dxim = grid[grid[i].W].dx;
			dzim = grid[grid[i].W].dz;
			kxim = grid[grid[i].W].kx;

			grid[i].Gxmh = (2.0*BETAC*dyi*dzi*dyi*dzim*kxi*kxim) 
				/ (kxi*dyi*dzi*dxim + kxim*dyi*dzim*dxi);
		}
		else{
			grid[i].Gxmh = 0.0;
		}
		
		if (grid[i].Gyph != NOTRANSMISSIBILITY || bool_N){
			dyjp = grid[grid[i].N].dy;
			dzjp = grid[grid[i].N].dz;
			kyjp = grid[grid[i].N].ky;

			grid[i].Gyph = (2.0*BETAC*dxi*dzi*dxi*dzjp*kyi*kyjp) 
				/ (kyi*dxi*dzi*dyjp + kyjp*dxi*dzjp*dyi);
		}
		else{
			grid[i].Gyph = 0.0;
		}
		
		if (grid[i].Gymh != NOTRANSMISSIBILITY || bool_S){
				dyjm = grid[grid[i].S].dy;
				dzjm = grid[grid[i].S].dz;
				kyjm = grid[grid[i].S].ky; 

				grid[i].Gymh = (2.0*BETAC*dxi*dzi*dxi*dzjm*kyi*kyjm)
					/(kyi*dxi*dzi*dyjm + kyjm*dxi*dzjm*dyi);
		}
		else{
			grid[i].Gymh = 0.0;
		}

		if(grid[i].Gxph < 0.0 || grid[i].Gxmh < 0.0 || grid[i].Gyph < 0.0
			|| grid[i].Gymh < 0.0){
				weprintf("negative geometrical factor G in block (%d, %d)",
					grid[i].row, grid[i].col);
		}
	}
	/*************************************************************************/

	unsetprogname();

	return;
}
/*****************************************************************************/

void setCylindricalTransmissibilities(Parameters *par, Block *grid, Boundary 
								 *boundary, Well *wells)
{
	int i, j, n, bcIndex;
	int bool_E, bool_W, bool_N, bool_S;
	double dri, ri, rim, rip, rimh, riph, kri, krip, krim;
	double doj, dojm, dojp, koj, kojp, kojm,dz,dzip,dzim,dzjp,dzjm;
	double **drtmp, *r, *rh;
	char drFileName[LENGTHFN];

	setprogname("set cylindrical transmissibilities");

	/*
	 _ _ _ _ _______       _______       _______ _ _ _ _ 
	|       |       | ... |       | ... |       |       |
	| r[0]  |  r[1] | ... |  r[i] | ... |  r[n] | r[n+1]|
	|_ _ _ _|_______| ... |_______| ... |_______|_ _ _ _|
	        |             |       |             | 
	      rh[0]        rh[i-1]  rh[i]         rh[n]   
	*/

	/* half-block radius *///0.0.3
	if (!atof(par->dxFile)){
		drtmp = rMatrix(1, par->ncol);
		sprintf(drFileName, "%s%s", par->projectDir, par->dxFile);
		rTableFile(drFileName, drtmp, 1, par->ncol);
	}
	else printf("Error!");

	n = par->ncol-1; //Depende se tem uma coluna de zeros!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	r  = rVector(n + 2);	
	rh = rVector(n + 1); 
	rh[0] = wells[0].rw;
	for (i=0; i < n; i++){
		rh[i+1] = rh[i] + drtmp[0][i];
		r[i+1]  = sqrt(rh[i]*rh[i+1]);		
	}

	/* boundary-block radius *///r12 = (r2 - r1) / log(r2/r1)
	r[0] = rh[0];
	while(abs(r[0]-r[1]*exp((r[0]-r[1])/rh[0])) > EPSILON){
		r[0] = r[1]*exp((r[0]-r[1])/rh[0]);
	}
	r[0] = r[1]*exp((r[0]-r[1])/rh[0]);
	r[n+1] = rh[n];
	while(abs(r[n+1]-r[n]-rh[n]*log(r[n+1]/r[n])) > EPSILON){
		r[n+1] = r[n] + rh[n]*log(r[n+1]/r[n]);
	}
	r[n+1] = r[n] + rh[n]*log(r[n+1]/r[n]);
	/****** constant part of the transmissibilities ******/
	for(i = 0; i < par->nBlocks; i++){

		/*********** for nonflow boundary conditions *********/
		if(grid[i].E < 0){ 
			grid[i].E = grid[i].H;
		}
		
		if(grid[i].W < 0){
			grid[i].W = grid[i].H;
		}
		
		if(grid[i].N < 0){
			grid[i].N = grid[i].H;
		}
		
		if(grid[i].S < 0){
			grid[i].S = grid[i].H;
		}
		/*****************************************************/
		
		bool_E = bool_W = FALSE;

		for(j=0; j < grid[i].nBoundary; j++){
			bcIndex = grid[i].boundary[j];
			bool_E = boundary[bcIndex].side == EAST ? TRUE : bool_E;
			bool_W = boundary[bcIndex].side == WEST ? TRUE : bool_W;
		}

		n   = grid[i].col+1;
		ri  = r[n];
		dz  = grid[i].dz;
		dri = rh[n]-rh[n-1];
		doj = grid[i].dy;
		kri	= grid[i].kx;
		koj = grid[i].ky;
		riph = rh[n];
		rimh = rh[n-1];
		
		if (grid[i].Gxph < 0 ){
			if( grid[i].E != grid[i].H || bool_E == TRUE){ 
				rip = r[n+1];
				dzip = grid[grid[i].E].dz;
				krip = grid[grid[i].E].kx;
				
				grid[i].Gxph =	  BETAC*doj*(kri*krip*dz*dzip) / 
						 ( krip*dzip*log(riph/ri) + kri*dz*log(rip/riph) );
			}
			else{
				grid[i].Gxph = 0.0;
			}
		}
		
		if (grid[i].Gxmh < 0 ){
			if( grid[i].W != grid[i].H || bool_W == TRUE){ 
				rim = r[n-1]; 
				dzim = grid[grid[i].E].dz;
				krim = grid[grid[i].W].kx;
				
				grid[i].Gxmh = -BETAC*doj*(kri*krim*dz*dzim) / 
						 ( krim*dzim*log(rimh/ri) + kri*dz*log(rim/rimh) );
			}
			else{
				grid[i].Gxmh = 0.0;
			}
		}
		
		if (grid[i].Gyph < 0 ){
			if( grid[i].N != grid[i].H ){ 
				dojp = grid[grid[i].N].dy;
				dzjp = grid[grid[i].N].dz;
				kojp = grid[grid[i].N].ky;
				grid[i].Gyph = BETAC*dri/ri*(2*doj*dojp)/(doj*kojp*dzjp + dojp*koj*dz);
			}
			else{
				grid[i].Gyph = 0.0;
			}
		}
		
		if (grid[i].Gymh < 0 ){
			if( grid[i].S != grid[i].H ){
				dojm = grid[grid[i].S].dy;
				dzjm = grid[grid[i].S].dz;
				kojm = grid[grid[i].S].ky;
				grid[i].Gymh = BETAC*dri/ri*(2*doj*dojm)/(doj*kojm*dzjm + dojm*koj*dz);
			}
			else{
				grid[i].Gymh = 0.0;
			}
		}
		
		if(grid[i].Gxph < 0.0 || grid[i].Gxmh < 0.0 || grid[i].Gyph < 0.0
			|| grid[i].Gymh < 0.0){
				weprintf("negative geometrical factor G in block (%d, %d)",
					grid[i].row, grid[i].col);
		}
	}
	/*************************************************************************/

	unsetprogname();

	return;
}
/*****************************************************************************/

double **readFluidProperties(Parameters *par)
/* read the fluid properties from a file (pressure, Formation Volume Factor, 
/* viscosity, especific gravity) and stores p, 1/FVF, 1/(FVF*mu) and gamma */
{
	int i;
	char dirFileName[LENGTHFN];
	double **fluidProps;
	
	setprogname("read fluid properties");

	par->invBSC = 1.0; 

	/* for pressure, formation-volume-factor, viscosity, especific gravity */
	fluidProps = rMatrix(par->nfprop, NCOLPROPS);
	sprintf(dirFileName, "%s%s", par->projectDir, par->fluidPropFile);
	rTableFile(dirFileName, fluidProps, par->nfprop, NCOLPROPS);

	for(i = 0; i < par->nfprop; i++){
		fluidProps[i][INV_FVF] = 1.0 / fluidProps[i][FVF];
		fluidProps[i][INV_FVFVISC] = 
			fluidProps[i][INV_FVF] / fluidProps[i][VISC];
	}

	unsetprogname();

	return fluidProps;
}
/*****************************************************************************/

void setInitialPressure(Parameters *par, Block *grid, Well *wells, 
						double **fprops, double *p, double *pp, double *ppp)
{
	int i = 0, j = 0;
	double delta, dP, totalVolume = 0.0, Jw;
	
	setprogname("set initial pressure");

	for(i = 0; i < par->nBlocks; i++){
		p[i] = par->p0 
			+ (grid[i].z - par->z0) 
			* getTableVal(fprops, PRESS, par->p0, GAMMA, par->dpprops, 
			par->nfprop);

		totalVolume += grid[i].Vac*grid[i].phi;
	}

	do{
		for(i = 0; i < par->nBlocks; i++){			
			pp[i] = p[i];
			
			p[i] = par->p0 
				+ (grid[i].z - par->z0) 
				* getTableVal(fprops, PRESS, pp[i], GAMMA, par->dpprops,
				par->nfprop);
		}
		delta = normDelta(p, pp, par->nBlocks);
	}while(delta > EPSILON);

	if (par->cf == 0.0){
		for(i = 0; i < par->nBlocks; i++) ppp[i] = pp[i] = p[i]=7750;
		unsetprogname();
		return;
	}

	/************* gross estimative of the first pressure change *************/
	for(i = 0; i < par->nBlocks; i++){
		ppp[i] = pp[i] = p[i]; /* for use with slightly compressible model */
		
		if(grid[i].isWellBlock != NOWELL){
			if(wells[grid[i].isWellBlock].type == RATE_ESPECIFIED){				
				dP = (wells[grid[i].isWellBlock].qw *(par->iniTime+par->dt))
					/ (par->cf * totalVolume);
				ppp[grid[i].N] = 
					ppp[grid[i].S] =
					ppp[grid[i].E] =
					ppp[grid[i].W] =
					ppp[grid[i].H] += dP;
			}
			if(wells[grid[i].isWellBlock].type == PRESSURE_ESPECIFIED){
				Jw = (wells[grid[i].isWellBlock].gw) 
					* getTableVal(fprops, PRESS, p[i], INV_FVFVISC, par->dpprops,
					par->nfprop);
				dP = (-Jw*(p[i] - wells[grid[i].isWellBlock].pwf))
					/ (par->cf * totalVolume);
				ppp[grid[i].N] = 
					ppp[grid[i].S] =
					ppp[grid[i].E] =
					ppp[grid[i].W] =
					ppp[grid[i].H] += dP;
			}
		}
	}
	/*************************************************************************/

	unsetprogname();

	return;
}
/*****************************************************************************/

double normDelta(double *a, double *b, int n)
{
	int i;
	double value = 0.0, aux;

	for(i = 0; i < n; i++){
		aux = (a[i] - b[i]);
		value += aux*aux;
	}
	/*value = sqrt(value / n);*/
	value = sqrt(value);

	return value;
}
/*****************************************************************************/

double getTableVal(double **table, int colx, double x, int colfx, double delta,
				   int rowmax)
/* only for tables with equal step size in the independent variable, x */
{
	int j;
	double fx;

	j = (int)((x - table[0][colx])/delta);

	if(j < 0 || j > (rowmax - 2))
		eprintf("pressure row out of range in properties file (%f)",x);

	fx = table[j][colfx] + 
		(x - table[j][colx]) * (table[j+1][colfx] - table[j][colfx]) 
		/ (table[j+1][colx] - table[j][colx]);

	return fx;
}
/*****************************************************************************/

void setMatrixStructure(Parameters *par, Block *grid, int *Ti, int *Tj, 
						double *Tx, int *Ap, int *Ai, double *Ax, int *Map)
{
	int i, j = 0, status;
	
	for(i = 0; i < par->nBlocks; i++){
		Ti[j] = i;
		Tj[j] = i;
		Tx[j] = i;
		j++;

		if(grid[i].E != grid[i].H){
			Ti[j] = i;
			Tj[j] = grid[i].E;
			Tx[j] = grid[i].E;
			j++;
		}
		if(grid[i].W != grid[i].H){
			Ti[j] = i;
			Tj[j] = grid[i].W;
			Tx[j] = grid[i].W;
			j++;
		}
		if(grid[i].N != grid[i].H){
			Ti[j] = i;
			Tj[j] = grid[i].N;
			Tx[j] = grid[i].N;
			j++;
		}
		if(grid[i].S != grid[i].H){
			Ti[j] = i;
			Tj[j] = grid[i].S;
			Tx[j] = grid[i].S;
			j++;
		}
	}

	status = umfpack_di_triplet_to_col(par->nBlocks, par->nBlocks, 
		par->nMatrix, Ti, Tj, Tx, Ap, Ai, Ax, Map);

	return;
}
/*****************************************************************************/

void setTransmissibilityMatrix(Parameters *par, Block *grid, Well *wells,
							   double **fprops, double *pN, double *p, 
							   double *Ax, double *b, int *Map, 
							   Boundary *boundary)
{
 	int i, j = 0;
	double Gamma_dt, E, EG, W, WG, N, NG, S, SG, HG, QG, Jw;
	
	for(i = 0; i < par->nBlocks; i++){
		
		if (par->cf != 0.0){
		/* For slightly compressible flow */
		Gamma_dt = (grid[i].Vac/par->dt)*(grid[i].phi*par->invBSC*par->cf);
		}
		else{
		/* For compressible and slightly compressible flow (see setInitialPressure)*/
		 Gamma_dt = (p[i] - pN[i]) != 0.0 ?
			(grid[i].Vac/par->dt)*grid[i].phi
			*(getTableVal(fprops, PRESS, p[i], INV_FVF, par->dpprops, par->nfprop)
			- getTableVal(fprops, PRESS, pN[i], INV_FVF, par->dpprops, par->nfprop))
			/(p[i] - pN[i]) : 0.0;//*/
		}
		
		if(grid[i].Gxph > 0.0){ 
			E = grid[i].Gxph*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].E]), 
				INV_FVFVISC, par->dpprops, par->nfprop);
			EG = E*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].E]), GAMMA, 
				par->dpprops, par->nfprop);
		}
		else{
			E  = 0.0;
			EG = 0.0;
		}
		
		if(grid[i].Gxmh > 0.0){ 
			W = grid[i].Gxmh*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].W]), 
				INV_FVFVISC, par->dpprops, par->nfprop);
			WG = W*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].W]), GAMMA, 
				par->dpprops, par->nfprop);
		}
		else{
			W  = 0.0;
			WG = 0.0;
		}

		if(grid[i].Gyph > 0.0){
			N = grid[i].Gyph*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].N]), 
				INV_FVFVISC, par->dpprops, par->nfprop);
			NG = N*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].N]), GAMMA, 
				par->dpprops, par->nfprop);
		}
		else{
			N  = 0.0;
			NG = 0.0;
		}
		
		if(grid[i].Gymh > 0.0){
			S = grid[i].Gymh*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].S]), 
				INV_FVFVISC, par->dpprops, par->nfprop);
			SG = S*getTableVal(fprops, PRESS, 0.5*(p[i]+p[grid[i].S]), GAMMA, 
				par->dpprops, par->nfprop);
		}
		else{
			S  = 0.0;
			SG = 0.0;
		}

		HG = -(EG + WG + NG + SG);
		QG =  EG*grid[grid[i].E].z
			+ WG*grid[grid[i].W].z
			+ NG*grid[grid[i].N].z
			+ SG*grid[grid[i].S].z
			+ HG*grid[i].z;
		
		Ax[Map[j]] = -(Gamma_dt + E + W + N + S);
		b[i] = -(Gamma_dt*pN[i] - QG);

		/* wells */
		if(grid[i].isWellBlock != NOWELL){
			if(wells[grid[i].isWellBlock].type == RATE_ESPECIFIED){
				Jw = (wells[grid[i].isWellBlock].gw) 
					* getTableVal(fprops, PRESS, p[i], INV_FVFVISC, 
					par->dpprops, par->nfprop);

				wells[grid[i].isWellBlock].pwf = 
					(wells[grid[i].isWellBlock].qw/Jw) + p[i];
				
				b[i] -= wells[grid[i].isWellBlock].qw;
			}
			if(wells[grid[i].isWellBlock].type == PRESSURE_ESPECIFIED){
				Jw = (wells[grid[i].isWellBlock].gw) 
					* getTableVal(fprops, PRESS, p[i], INV_FVFVISC, 
					par->dpprops, par->nfprop);
				
				wells[grid[i].isWellBlock].qw = 
					-Jw*(p[i] - wells[grid[i].isWellBlock].pwf);
				
				Ax[Map[j]] -= Jw;
				b[i] -= Jw*wells[grid[i].isWellBlock].pwf;
			}
		}
		
		if(grid[i].nBoundary == 0){ //NOBOUNDARY 

			j++;  // diagonal (equation) 

			if(E > 0.0){
				Ax[Map[j]] = E;
				j++;	
			}
			if(W > 0.0){
				Ax[Map[j]] = W;
				j++;
			}
			if(N > 0.0){
				Ax[Map[j]] = N;
				j++;
			}
			if(S > 0.0){
				Ax[Map[j]] = S;
				j++;
			}
		}
		else 
			setBoundaryConditions(&j, Map, Ax, &b[i], N, E, W, S, boundary, 
								  &grid[i]);
	}

	return;
}
/*****************************************************************************/

void readIniFile(int argc, char *argv[], Parameters *par)
/* this function reads a list of parameters from a inifile and set 
the simulation parameters structure */
{
	dictionary *ini;
	char *str, *p;
	char iniFile[LENGTHFN];
	int pos;
	
	setprogname("read ini file");
	
	/********* reading project name *********/
	if (argc < 2) {
		printf("<inifile> not specified...using default\n\n");
		strcpy(par->projectName, "Projects/fiveSpot");
		 
	} 
	else
		strcpy(par->projectName, argv[1]);

	sprintf(par->projectDir, "%s/", par->projectName);

	//if the project is inside other folder
	p = strtok(par->projectName, "/");
	while (p!=NULL){
		strcpy(par->projectName, p);
		p = strtok(NULL, "/");
	}
	
	sprintf(par->iniFile, "%s%s.ini", par->projectDir, par->projectName);
	
	ini = iniparser_load(par->iniFile);
	if(ini == NULL) {
		eprintf("cannot parse ini file:");
	}
	/****************************************/

	/******* reading data in ini file *******/
    par->iniTime = iniparser_getdouble(ini, "control:initialtime", -1);
	par->dt      = iniparser_getdouble(ini, "control:dt", -1);
	par->dtMultiplier = iniparser_getdouble(ini, "control:dtMultiplier", 1);
	par->dtLogScale = iniparser_getboolean(ini, "control:dtLogScale", FALSE);
	par->isCylindrical = iniparser_getboolean(ini, 
		"reservoir description:cylindrical", FALSE); 
	par->nSteps  = iniparser_getint(ini, "control:maxsteps", -1);
	par->ncol = iniparser_getint(ini, "reservoir description:ncol", -1);
	par->nrow = iniparser_getint(ini, "reservoir description:nrow", -1);
	par->nlay = iniparser_getint(ini, "reservoir description:nlay", 1);
	par->cf   = iniparser_getdouble(ini, "reservoir description:cf", -1);
	par->cphi = iniparser_getdouble(ini, "reservoir description:cr", -1);
	par->p0   = iniparser_getdouble(ini, "reservoir description:refpres", -1);
	par->z0   = iniparser_getdouble(ini, "reservoir description:refdepth", -1);
	par->rhoSC = iniparser_getdouble(ini, "fluid properties:rho", -1);	
	par->nfprop = iniparser_getint(ini, "fluid properties:fproprows", -1);
	par->dpprops = iniparser_getdouble(ini, "fluid properties:dp", -1);

	str = iniparser_getstring(ini, "control:outputfile", NULL);
	strcpy(par->outputFile, str); 
	
	str = iniparser_getstring(ini, "reservoir description:geo", NULL);
	strcpy(par->geometryFile, str);
	
	str = iniparser_getstring(ini, "reservoir description:modified", UNDEF);
	strcpy(par->modifiedBlocksFile, str);

	str = iniparser_getstring(ini, "reservoir description:modTransmissibility", UNDEF);
	strcpy(par->modTransmissibilityFile, str); 
	
	str = iniparser_getstring(ini, "reservoir description:phi", NULL);
	strcpy(par->porosityFile, str);

	str = iniparser_getstring(ini, "reservoir description:kx", NULL);
	strcpy(par->kxFile, str);

	str = iniparser_getstring(ini, "reservoir description:ky", NULL);
	strcpy(par->kyFile, str);

	str = iniparser_getstring(ini, "reservoir description:dz", NULL);
	strcpy(par->thicknessFile, str);

	str = iniparser_getstring(ini, "reservoir description:ztop", NULL);
	strcpy(par->zTopFile, str);

	str = iniparser_getstring(ini, "reservoir description:dx", NULL);
	strcpy(par->dxFile, str);

	str = iniparser_getstring(ini, "reservoir description:dy", NULL);
	strcpy(par->dyFile, str);

	str = iniparser_getstring(ini, "control:bcFile", UNDEF);
	strcpy(par->bcFile, str);
	
	str = iniparser_getstring(ini, "fluid properties:propsfile", NULL);
	strcpy(par->fluidPropFile, str);
	
	str = iniparser_getstring(ini, "wells:wellsfile", UNDEF); 
	strcpy(par->wellsFile, str);
	
	iniparser_freedict(ini);
	unsetprogname();

	return;
}
/*****************************************************************************/

Well *readAndSetWellParameters(Parameters *par)
{
	Well *wells;
	dictionary *ini;
	int i, nsec;
	char SectionName[LENGTHSN];
	char num[3];
	char str[LENGTHSP], type[LENWELLTYPE], dirFileName[LENGTHFN];
	
	if(strcmp(par->wellsFile, UNDEF) == 0){
		par->nwells = 0;
		return 0;
	}

	setprogname("well parameters");

	sprintf(dirFileName, "%s%s", par->projectDir, par->wellsFile);
	ini = iniparser_load(dirFileName);
	if(ini == NULL) 
		eprintf("cannot parse wellfile:");
	

	/**** Number of sections ****/
	nsec = 0;			// Number of sections in file
	i = 0;
	while(i < iniparser_getnsec(ini)){
		sprintf(str, "Well %i", nsec+1); 
		if (iniparser_find_entry(ini, str)) nsec++;
		i++;
	}
	/****************************/

	par->nwells = par->isCylindrical ? 1 : 
		iniparser_getint(ini, "main:nwells", 0); 

	if (nsec > par->nwells){
		weprintf("inconsistencies in the wells configuration file");
		weprintf("the last %i wells are not being considered", 
			   nsec-par->nwells);
	}
	else if(nsec < par->nwells ){
		weprintf("inconsistencies in the wells configuration file");
		weprintf("there are only %i well sections in order. adopting nwells = %i",
			   nsec, par->nwells = nsec);
	}

	wells = (Well *)calloc(par->nwells, sizeof(Well));
	if (!wells) eprintf("allocation failure for wells struct:");

	for(i = 0; i < par->nwells; i++){
		itoa(i+1, num, 10);  
		sprintf(SectionName, "well %s:", num);
		
		sprintf(str, "%s%s", SectionName, "qw");
		wells[i].qw = iniparser_getdouble(ini, str, 0);

		sprintf(str, "%s%s", SectionName, "pwf");
		wells[i].pwf = iniparser_getdouble(ini, str, 0);
	
		sprintf(str, "%s%s", SectionName, "type");
		strcpy(type, _strupr(iniparser_getstring(ini, str, NULL)));  

		sprintf(str, "%s%s", SectionName, "row");
		wells[i].row = iniparser_getint(ini, str, -1);

		sprintf(str, "%s%s", SectionName, "col");
		wells[i].col = iniparser_getint(ini, str, -1);
		
		sprintf(str, "%s%s", SectionName, "dx");
		wells[i].dx = iniparser_getdouble(ini, str, 0); 

		sprintf(str, "%s%s", SectionName, "dy");
		wells[i].dy = iniparser_getdouble(ini, str, 0); 

		sprintf(str, "%s%s", SectionName, "rw");
		wells[i].rw = iniparser_getdouble(ini, str, -1);

		sprintf(str, "%s%s", SectionName, "h");
		wells[i].h  = iniparser_getdouble(ini, str, -1);

		sprintf(str, "%s%s", SectionName, "s");
		wells[i].s  = iniparser_getdouble(ini, str, 0);

		wells[i].np = 0; 

		if(wells[i].row < 0 || wells[i].row >= par->nrow){ 
			eprintf("row in [Well %i] outside of reservoir geometry", i+1); 
		}

		if(wells[i].col < 0 || wells[i].col >= par->ncol){
			eprintf("col in [Well %i] outside of reservoir geometry", i+1);
		}
		
		if (strcmp(type, "RATE_ESPECIFIED") == 0 ||	strcmp(type, "RATE") == 0)
			wells[i].type = RATE_ESPECIFIED;
		else if (strcmp(type, "PRESSURE_ESPECIFIED") == 0 ||
			     strcmp(type, "PRESSURE") == 0)
			wells[i].type = PRESSURE_ESPECIFIED;			
		else
			eprintf("type in [Well %i] is invalid", i+1);
	}
	
	iniparser_freedict(ini);

	unsetprogname();
	
	return wells;
}
/*****************************************************************************/

void writePressureFile(Parameters *par, int step, Block *grid, double *pressure)
{
	FILE *fp;
	char extension[] = OUTPRESSUREEXT, fileName[64];
	int col, row, s;
	
	sprintf(fileName, "%s_%.6d%s", par->outPressureFile, step, extension);
	
	if ( (fp = fopen(fileName, "w")) == NULL) {
		eprintf("opening output file");
	}
	
	s = 0;
	while(s < par->nBlocks) {
		for (row = 0; row < par->nrow; row++) {
			for (col = 0; col < par->ncol; col++) {
				
				if (grid[s].col == col && 
					grid[s].row == row) {
						fprintf(fp, "%g ", pressure[s]);
						s++;
					}
					else {
						fprintf(fp, "%d ", 0);
					}
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);	

	return;
}
/*****************************************************************************/

void setWells(Parameters *par, Block *grid, int **geom, Well *wells)
{
	int i, j, localIndex;
	double req, kH, offsetx, offsety;

	setprogname("well parameters");

	/****** wells ******/
	for(i = 0; i < par->nwells; i++){
		localIndex = (geom[wells[i].row][wells[i].col] - 1);
		grid[localIndex].isWellBlock = i;
		
		offsetx = fabs(wells[i].dx)/grid[localIndex].dx;
		offsety = fabs(wells[i].dy)/grid[localIndex].dy; 

		/* tolerance 1% - interior well block */
		if(offsetx < 0.01 && offsety < 0.01
			&& grid[localIndex].E != grid[localIndex].H
			&& grid[localIndex].W != grid[localIndex].H
			&& grid[localIndex].N != grid[localIndex].H
			&& grid[localIndex].S != grid[localIndex].H){

			req = reqPeaceman(grid, localIndex);
			//req = 0.14*sqrt(pow(grid[localIndex].dx,2.0)+pow(grid[localIndex].dy,2.0));
			//req = reqDing(grid, localIndex, req);
		}
		else if(grid[localIndex].N != grid[localIndex].S 
			&&  grid[localIndex].E != grid[localIndex].W) {
			
			req = reqAbouKassemAziz(grid, localIndex);
		}
		else{
			weprintf("\n\nWARNING????? Well model??\n\n");
		}

		kH = sqrt(grid[localIndex].kx * grid[localIndex].ky);

		wells[i].gw = (2.0*PI*BETAC*kH*wells[i].h) 
			/ (log(req/wells[i].rw) + wells[i].s);
	}

	unsetprogname();
	
	return;
}
/*****************************************************************************/

double reqPeaceman(Block *grid, int localIndex)
/* Peaceman D. W., Interpretation of Well-Block pressures in numerical 
reservoir simulation with nonsquare grid blocks and anisotropic permeability,
SPEJ (June 1978), pp. 531-43. */
{
	double req;
	
	req = 0.28*sqrt( sqrt(grid[localIndex].ky/grid[localIndex].kx) 
		* grid[localIndex].dx * grid[localIndex].dx 
		+ sqrt(grid[localIndex].kx/grid[localIndex].ky) 
		* grid[localIndex].dy * grid[localIndex].dy )
		/ ( pow(grid[localIndex].kx/grid[localIndex].ky, 0.25) 
		+ pow(grid[localIndex].ky/grid[localIndex].kx, 0.25) );

	return req;
}
/*****************************************************************************/

double reqDing(Block *grid, int localIndex, double req)
/* Yu Ding, Gérard Renard, Luce Weill, Representation of Wells in Numerical 
Reservoir Simulation, SPE (February 1998), pp. 18-23. */
{
	double dxi, dyi, dxip, dxim, dyjp, dyjm;
	double alphaN, alphaS, alphaE, alphaW;

	//double req = rw;
	
	//if(grid[localIndex].isWellBlock == RATE_ESPECIFIED){
		dxi = grid[localIndex].dx;
		dyi = grid[localIndex].dy;
		
		dxip = grid[grid[localIndex].E].dx;
		dxim = grid[grid[localIndex].W].dx;
		
		dyjp = grid[grid[localIndex].N].dy;
		dyjm = grid[grid[localIndex].S].dy;

		alphaE = (dxi + dxip)*atan(dyi/dxi)/(dyi*log((dxi + dxip)/(2.0*req)));
		alphaW = (dxi + dxim)*atan(dyi/dxi)/(dyi*log((dxi + dxim)/(2.0*req)));
		alphaN = (dyi + dyjp)*atan(dxi/dyi)/(dxi*log((dyi + dyjp)/(2.0*req)));
		alphaS = (dyi + dyjm)*atan(dxi/dyi)/(dxi*log((dyi + dyjm)/(2.0*req)));

		// on this block
		grid[localIndex].Gxph *= alphaE; // EAST			
		grid[localIndex].Gxmh *= alphaW; // WEST
		grid[localIndex].Gyph *= alphaN; // NORTH
		grid[localIndex].Gymh *= alphaS; // SOUTH

		// on its neighborhood
		grid[grid[localIndex].E].Gxmh *= alphaE; // EAST			
		grid[grid[localIndex].W].Gxph *= alphaW; // WEST
		grid[grid[localIndex].N].Gymh *= alphaN; // NORTH
		grid[grid[localIndex].S].Gyph *= alphaS; // SOUTH
	//}
	return req;
}
/*****************************************************************************/

double reqAbouKassemAziz(Block *grid, int localIndex)
/* Abou-Kassem J. H. and Aziz K., Analytical well models for reservoir 
simulation, SPEJ (Aug 1985), pp. 573-79. */ 
{
	double req;
	
	req = 0.0;

	return req;
}
/*****************************************************************************/

void writeOuts(FILE* local_write, int nBlocksOut, Block* grid, Well *wells, 
			   Out* outs, double time, double *pressure, double *cumulativeProd)
{
	int i, j;
	
	fprintf(local_write, "\n %.3E ", time);

	if(*cumulativeProd != NONP) 
		fprintf(local_write, "%12.3E ", *cumulativeProd);
	
	for(i = 0; i < nBlocksOut; i++){	
		if(outs[i].p)
			fprintf(local_write, "%12.3E ", pressure[outs[i].localIndex]);

		if(outs[i].pwf)
			fprintf(local_write, "%12.3E ", 
				wells[grid[outs[i].localIndex].isWellBlock].pwf);

		if(outs[i].qw)
			fprintf(local_write, "%12.3E ", 
				wells[grid[outs[i].localIndex].isWellBlock].qw);
		
		if(outs[i].np) 
			fprintf(local_write, "%12.3E ", 
				wells[grid[outs[i].localIndex].isWellBlock].np);
	}
	fflush(local_write);
}
/*****************************************************************************/

void header(FILE *file, Parameters *par)
{
	fprintf(file,"\
   =========================================================================  \n\n\
                    SINGLE PHASE - SLIGHTLY-COMPRESSIBLE-FLOW                 \n\
       CARTESIAN (variable depth/variable thickness) - FINITE DIFFERENCE      \n\
                              BLOCK CENTERED GRID                             \n\
                              SIMULATOR Ver. 0.0.3                            \n\n\
                              capico@lenep.uenf.br                            \n\
                           bismarckjunior@outlook.com                         \n\n\
   =========================================================================  \n\n\
     Project Name: %s.ini ", par->projectName);
	
	fprintf(file, "\n\n\
  nrow:                [%d]\n\
  ncol:                [%d]\n\
  nlay:                [%d]\n\
  initialTime:         [%g]\n\
  maxSteps:            [%d]\n\
  dt:                  [%g]\n\
  cf:                  [%g]\n\
  cphi:                [%g]\n\
  geometry:            [%s]\n\
  modGeom:             [%s]\n\
  modTransmissibility: [%s]\n\
  porosity:            [%s]\n\
  kx:                  [%s]\n\
  ky:                  [%s]\n\
  thickness:           [%s]\n\
  zTop:                [%s]\n\
  dx:                  [%s]\n\
  dy:                  [%s]\n\
  rhoSC:               [%g]\n\
  fproprows:           [%d]\n\
  dpprops:             [%g]\n\
  props file:          [%s]\n\
  p0:                  [%g]\n\
  z0:                  [%g]\n\
  wells file:          [%s]\n\
  output file:         [%s]\n\
  boundary file:       [%s]\n\n\
  ----------------------------------------------------------------------------\n", 
		par->nrow, par->ncol, par->nlay, par->iniTime, par->nSteps, par->dt, par->cf, 
		par->cphi, par->geometryFile, par->modifiedBlocksFile, 
		par->modTransmissibilityFile,par->porosityFile, 
		par->kxFile, par->kyFile, par->thicknessFile, par->zTopFile, 
		par->dxFile, par->dyFile, par->rhoSC, par->nfprop, par->dpprops, 
		par->fluidPropFile, par->p0, par->z0, par->wellsFile, 
		par->outputFile, par->bcFile);
}
/*****************************************************************************/

void writeTableHeader(FILE *local_out, Parameters *par, Out *outs, int n_outs,
					  double *cumulativeProd)
{
	int i;
	char line[LINESIZE];
	char types[LINESIZE];
	char blocks[LINESIZE];
	char str[LENGTHFN];
	
	sprintf(line,  "-----------");
	sprintf(types, "  Time[Day]  ");
	sprintf(blocks,"             ");
	
	if (*cumulativeProd != NONP){ 
		strcat(line,  "-------------");
		strcat(types, "  Np [bbl]   ");
		sprintf(str,  "             ");
		strcat(blocks,str);
	}

	for(i = 0; i < n_outs; i++){
		if(outs[i].p){
			strcat(line,  "-------------");
			strcat(types, "  P [psi]    ");
			sprintf(str, " (%03i,%03i)   ", outs[i].col, outs[i].row); 
			strcat(blocks,str);
		}
		if(outs[i].pwf){
			strcat(line,  "-------------");
			strcat(types, " Pwf [psi]   ");
			sprintf(str, " (%03i,%03i)   ", outs[i].col, outs[i].row);
			strcat(blocks, str);
		}
		if(outs[i].qw){
			strcat(line,  "-------------");
			strcat(types, " Qw[bbl/d]   ");
			sprintf(str, " (%03i,%03i)   ", outs[i].col, outs[i].row);
			strcat(blocks, str);
		}
		if(outs[i].np){ 
			strcat(line,  "-------------");
			strcat(types, "  Np [bbl]   ");
			sprintf(str, " (%03i,%03i)   ", outs[i].col, outs[i].row);
			strcat(blocks, str);
		}
	}
	if (n_outs != 0)
		fprintf(local_out, "\n%s\n%s\n%s\n%s", line, types, blocks, line);
	else
		fprintf(local_out, "\n%s\n%s\n%s", line, types, line);
	fflush(local_out);
}
/*****************************************************************************/

Out *readAndSetOuts(Parameters *par, Block *grid, int **geom, char *type, 
					double *cumulativeProd)
{
	dictionary *ini;
	Out *outs;
	int i, nCol, nOutSec, nOutEsp, maxcol;
	char SectionName[LENGTHSN], str[LENGTHSP], dirFileName [LENGTHFN], num[3];
		
	sprintf(str, "%s output", type);
	setprogname(str);
	
	sprintf(dirFileName, "%s%s", par->projectDir, par->outputFile);
	ini = iniparser_load(dirFileName);
	if(ini == NULL) 
		eprintf("cannot parse output:");
		
	/**** Number of sections ****/
	nOutSec = 0;			// Number of sections in file
	i = 0;
	while(i < iniparser_getnsec(ini)){
		sprintf(str, "%s %i", type, nOutSec+1); 
		if (iniparser_find_entry(ini, str)) nOutSec++;
		i++;
	}
	/****************************/

	// Number of sections especified
	sprintf(str, "main:n%sBlocks", type);
	nOutEsp = iniparser_getint(ini, str, 0); 

	if (nOutEsp < nOutSec){
		weprintf("inconsistencies in the ouput configuration file");
		weprintf("the last %i output are not being considered", 
			   nOutSec-nOutEsp);
	}
	else if(nOutEsp > nOutSec){
		weprintf("inconsistencies in the output configuration file");
		weprintf("there are only %i output sections in order. \
			\t\tadopting n%sBlocks = %i",
			   nOutSec, type, nOutEsp = nOutSec);
	}

	if(strcmp(type, "report") == 0)
		par->nReportBlocks = nOutEsp;
	else if (strcmp(type, "display") == 0)
		par->nDisplayBlocks = nOutEsp;

	outs = (Out *)calloc(nOutEsp, sizeof(Out));
	if (!outs) eprintf("allocation failure for out struct:");
	
	nCol=0; 
	
	maxcol = strcmp(type, "display") == 0 ? MAXCOLDISPLAY : 
		     strcmp(type, "report")  == 0 ? MAXCOLREPORT  : 0;
	
	if ( *cumulativeProd != NONP ) nCol++;  

	for(i = 0; i < nOutEsp; i++){
		itoa(i+1, num, 10); 
		sprintf(SectionName, "%s %s:", type, num);
		
		sprintf(str, "%s%s", SectionName, "row");
		outs[i].row = iniparser_getint(ini, str, -1);
	
		sprintf(str, "%s%s", SectionName, "col");
		outs[i].col = iniparser_getint(ini, str, -1);

		if(outs[i].row <0 || outs[i].row >= par->nrow) 
			eprintf("row in [%s %i] outside of reservoir geometry!", type, i+1); 
			
		if(outs[i].col <0 || outs[i].col >= par->ncol)
			eprintf("col in [%s %i] outside of reservoir geometry!", type, i+1);

		outs[i].localIndex = geom[outs[i].row][outs[i].col]-1;

		sprintf(str, "%s%s", SectionName, "p");
		outs[i].p = iniparser_getboolean(ini, str, FALSE);
		if (outs[i].p){
			nCol++;
			if(nCol>maxcol) outs[i].p = FALSE; 
		}

		sprintf(str, "%s%s", SectionName, "pwf");
		outs[i].pwf = grid[outs[i].localIndex].isWellBlock == NOWELL ?
						 FALSE : iniparser_getboolean(ini, str, FALSE);
		if (outs[i].pwf){
			nCol++;
			if(nCol>maxcol) outs[i].pwf = FALSE; 
		}

		sprintf(str, "%s%s", SectionName, "qw");
		outs[i].qw = grid[outs[i].localIndex].isWellBlock == NOWELL ?
						 FALSE : iniparser_getboolean(ini, str, FALSE);
		if (outs[i].qw){
			nCol++;
			if(nCol>maxcol) outs[i].qw = FALSE;
		}
		
		sprintf(str, "%s%s", SectionName, "np");
		outs[i].np = grid[outs[i].localIndex].isWellBlock == NOWELL ?
						 FALSE : iniparser_getboolean(ini, str, FALSE);
		if (outs[i].np){
			nCol++;
			if(nCol>maxcol) outs[i].np = FALSE;
		}		

		if(nCol>maxcol) {
			if(strcmp(type, "report") == 0)
				par->nReportBlocks = i+1;
			else if (strcmp(type, "display") == 0)
				par->nDisplayBlocks = i+1;
			break;
		}
	}
	
	if(nCol>maxcol){
		weprintf("inconsistencies in the ouput configuration file");
		weprintf("The maximum number of column for %s is %i. \n\r\
                %sing only %i columns.", type, maxcol, type, maxcol); 
	}
	iniparser_freedict(ini);
	unsetprogname();
	return outs;
}


/*****************************************************************************/

void *readAndSetMainOut(Parameters *par, double *cumulativeProd)
{
	dictionary *ini;
	char dirFileName[LENGTHFN];
	char *reportFile, *outPressureFile;
	
	setprogname("output data");

	sprintf(dirFileName, "%s%s", par->projectDir, par->outputFile);
	ini = iniparser_load(dirFileName);
	if(ini == NULL) 
		eprintf("cannot parse output:");
	
	par->displaySteps = iniparser_getint(ini, "main:displaySteps", 0);
	par->reportSteps  = iniparser_getint(ini, "main:reportSteps", 0);
	par->outPressureSteps = iniparser_getint(ini, "main:outPressureSteps", 0);
	
	if (TRUE==iniparser_getboolean(ini, "main:np", FALSE))
		*cumulativeProd = 0;
	
	reportFile = iniparser_getstring(ini, "main:reportFile", "report");
	outPressureFile = iniparser_getstring(ini, "main:outPressureFile", 
		"outPressure");
	
	if(*cumulativeProd != NONP || par->reportSteps!=0 || 
		par->outPressureSteps != 0){
		sprintf(dirFileName, "%s%s", par->projectDir, "Outs");
		mkdir(dirFileName);
	}
	else
		strcpy(dirFileName, par->projectDir);
		
	sprintf(par->reportFile, "%s/%s%s", dirFileName, reportFile, REPORTEXT); 
	sprintf(par->outPressureFile, "%s/%s", dirFileName, outPressureFile); 

	iniparser_freedict(ini);
	unsetprogname();
}
/*****************************************************************************/

void modifyGeometry(int **geom, Parameters *par){
	FILE *mf;
	char direction, dirFileName[LENGTHFN], type[LENBLOCKTYPE];
	int filler, row, col, rowi=0, coli=0, rowf=0, colf=0, line=0;

	setprogname("modify geometry");

	sprintf(dirFileName, "%s%s", par->projectDir, par->modifiedBlocksFile);
	if ( (mf = fopen(dirFileName, "r")) == NULL )
		eprintf("opening modified blocks file:");
				
	while( fscanf(mf,"\n%c ", &direction)>0 ){
		line++;

		switch(direction){
			case 'H':
			case 'h':
				if(fscanf(mf,"%d %d %d %s", &row, &coli, &colf, &type)==4){
					filler=strcmp(_strupr(type),"ACTIVE") == 0 ? ACTIVEBLOCK:
						   strcmp(_strupr(type),"INACTIVE")==0 ? INACTIVEBLOCK: 
						   INVALIDBLOCK;
					if ( coli < 0 || coli >= par->ncol || row < 0 ||
						 colf < 0 || colf >= par->ncol || row >= par->nrow )
						 break;

					for ( col=coli; col<=colf; col+= (colf>coli)-(colf<coli) )
						geom[row][col] = filler;
				}
				else
					eprintf("line %d in \"%s\" could not be understood", line, 
							par->modifiedBlocksFile);
				col = colf;
				break;

			case 'V':
			case 'v':
				if(fscanf(mf,"%d %d %d %s", &col, &rowi, &rowf, &type)==4){
					filler=strcmp(_strupr(type),"ACTIVE") == 0 ? ACTIVEBLOCK:
						   strcmp(_strupr(type),"INACTIVE")==0 ? INACTIVEBLOCK: 
						   INVALIDBLOCK;
					if ( rowi < 0 || rowi >= par->nrow || col < 0 ||
						 rowf < 0 || rowf >= par->nrow || col >= par->ncol )
						 break;

					for ( row=rowi; row<=rowf; row+=(rowf>rowi)-(rowf<rowi) )
						geom[row][col] = filler;
				}
				else
					eprintf("line %d in \"%s\" could not be understood", line, 
							par->modifiedBlocksFile);
				row = rowf;
				break;

			case 'B':
			case 'b':
				if(fscanf(mf,"%d %d %s", &row, &col, &type)==3){
					filler=strcmp(_strupr(type),"ACTIVE") == 0 ? ACTIVEBLOCK:
						   strcmp(_strupr(type),"INACTIVE")==0 ? INACTIVEBLOCK: 
						   INVALIDBLOCK;
					if ( row<0 || row>=par->nrow || col<0 || col>=par->ncol)
						 break;
					geom[row][col] = filler;
				}
				else
					eprintf("line %d in \"%s\" could not be understood", line, 
							par->modifiedBlocksFile);
				break;

			default:
				eprintf("line %d in \"%s\" could not be understood", line, 
						par->modifiedBlocksFile);
				break;
		}

		/*********** Exceptions **********/
		if( rowi < 0 || rowi >= par->nrow || row < 0 ||  
			rowf < 0 || rowf >= par->nrow || row >= par->nrow)
			eprintf("row in line %d in \"%s\" is outside of reservoir geometry", 
					line, par->modifiedBlocksFile);
		if( coli < 0 || coli >= par->ncol || col < 0 || 
			colf < 0 || colf >= par->ncol || col >= par->ncol)
			eprintf("col in line %d in \"%s\" is outside of reservoir geometry", 
					line, par->modifiedBlocksFile);
		/********************************/
	}
	fclose(mf);
	unsetprogname();
}

/*****************************************************************************/

void modifyTransmissibilities(int **geom, Block *grid, Parameters *par){
	FILE *mf;
	char direction, dirFileName[LENGTHFN], side[LENBLOCKTYPE];
	int localIndex, row, col, rowi=0, coli=0, rowf=0, colf=0, line=0;
	float filler;

	setprogname("modify transmissibility");
	
	sprintf(dirFileName, "%s%s", par->projectDir, 
		par->modTransmissibilityFile);

	if ( (mf = fopen(dirFileName, "r")) == NULL )
		eprintf("opening modified transmissibilities file:");
				
	while( fscanf(mf,"\n%c ", &direction)>0 ){
		line++;

		switch(direction){ 
			case 'H':
			case 'h':
				if(fscanf(mf,"%d %d %d %s", &row, &coli, &colf, &side)==4){
					if ( coli < 0 || coli >= par->ncol || row < 0 ||
						 colf < 0 || colf >= par->ncol || row >= par->nrow )
						 break;

					switch(_strupr(side)[0]){
						case NORTH:
							for ( col=coli; col<=colf; col+= (colf>coli)-(colf<coli) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gymh = NOTRANSMISSIBILITY; 
								grid[grid[localIndex].N].Gyph = NOTRANSMISSIBILITY; 
							}
							break;
						case SOUTH:
							for ( col=coli; col<=colf; col+= (colf>coli)-(colf<coli) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gyph = NOTRANSMISSIBILITY; 
								grid[grid[localIndex].S].Gymh = NOTRANSMISSIBILITY; 
							}
							break;
						case WEST:
							for ( col=coli; col<=colf; col+= (colf>coli)-(colf<coli) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gxmh = NOTRANSMISSIBILITY; 
								grid[grid[localIndex].W].Gxph = NOTRANSMISSIBILITY; 
							}
							break;
						case EAST:
							for ( col=coli; col<=colf; col+= (colf>coli)-(colf<coli) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gxph = NOTRANSMISSIBILITY; 
								grid[grid[localIndex].E].Gxmh = NOTRANSMISSIBILITY; 
							}
							break;
					}					
				}
				else
					eprintf("line %d in \"%s\" could not be understood", line, 
							par->modTransmissibilityFile);
				col = colf;
				break;

			case 'V':
			case 'v':
				if(fscanf(mf,"%d %d %d %s", &col, &rowi, &rowf, &side)==4){
					if ( rowi < 0 || rowi >= par->nrow || col < 0 ||
						 rowf < 0 || rowf >= par->nrow || col >= par->ncol )
						 break;

					switch(_strupr(side)[0]){
						case NORTH:
							for ( row=rowi; row<=rowf; row+=(rowf>rowi)-(rowf<rowi) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gymh = NOTRANSMISSIBILITY; 
								grid[grid[localIndex].N].Gyph = NOTRANSMISSIBILITY; 
							}
							break;
						case SOUTH:
							for ( row=rowi; row<=rowf; row+=(rowf>rowi)-(rowf<rowi) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gyph = NOTRANSMISSIBILITY; 
								grid[grid[localIndex].S].Gymh = NOTRANSMISSIBILITY; 
							}
							break;
						case WEST:
							for ( row=rowi; row<=rowf; row+=(rowf>rowi)-(rowf<rowi) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gxmh = NOTRANSMISSIBILITY; 
								grid[grid[localIndex].W].Gxph = NOTRANSMISSIBILITY; 
							}
							break;
						case EAST:
							for ( row=rowi; row<=rowf; row+=(rowf>rowi)-(rowf<rowi) ){
								localIndex = geom[row][col]-1;
								grid[localIndex].Gxph = NOTRANSMISSIBILITY;
								grid[grid[localIndex].E].Gxmh = NOTRANSMISSIBILITY;
							}
							break;
					}	
				}
				else
					eprintf("line %d in \"%s\" could not be understood", line, 
							par->modTransmissibilityFile);
				row = rowf;
				break;

			case 'B':
			case 'b':
				if(fscanf(mf,"%d %d %s", &row, &col, &side)==3){
					if ( row<0 || row>=par->nrow || col<0 || col>=par->ncol)
						 break;
					localIndex = geom[row][col]-1;

					switch(_strupr(side)[0]){
						case NORTH:
							grid[localIndex].Gymh = NOTRANSMISSIBILITY;
							grid[grid[localIndex].N].Gyph = NOTRANSMISSIBILITY;
							break;
						case SOUTH:
							grid[localIndex].Gyph = NOTRANSMISSIBILITY;
							grid[grid[localIndex].S].Gymh = NOTRANSMISSIBILITY;
							break;
						case WEST:
							grid[localIndex].Gxmh = NOTRANSMISSIBILITY;
							grid[grid[localIndex].W].Gxph = NOTRANSMISSIBILITY;
							break;
						case EAST:
							grid[localIndex].Gxph = NOTRANSMISSIBILITY;
							grid[grid[localIndex].E].Gxmh = NOTRANSMISSIBILITY;
							break;
					}

				}
				else
					eprintf("line %d in \"%s\" could not be understood", line, 
							par->modTransmissibilityFile);
				break;

			default:
				eprintf("line %d in \"%s\" could not be understood", line, 
						par->modTransmissibilityFile);
				break;
		}

		/*********** Exceptions **********/
		if( rowi < 0 || rowi >= par->nrow || row < 0 ||  
			rowf < 0 || rowf >= par->nrow || row >= par->nrow)
			eprintf("row in line %d in \"%s\" is outside of reservoir geometry", 
				line, par->modTransmissibilityFile);
		if( coli < 0 || coli >= par->ncol || col < 0 || 
			colf < 0 || colf >= par->ncol || col >= par->ncol)
			eprintf("col in line %d in \"%s\" is outside of reservoir geometry", 
				line, par->modTransmissibilityFile);
		/********************************/
	}
	fclose(mf);
	unsetprogname();
}

/*****************************************************************************/


Boundary *readAndSetBoundaryConditions(Parameters *par, int **geom, Block *grid)
{
	Boundary *boundary;
	dictionary *ini;
	int i, j, nsec, localIndex;
	int add_col, add_row, drow, dcol, nblocks;
	char num[3], SectionName[LENGTHSN];
	char side[SIDENAME], str[LENGTHSP], type[LENBCTYPE], dirFileName[LENGTHFN];
	
	setprogname("boundary conditions");

	if(strcmp(par->bcFile, UNDEF) == 0){
		par->nBoundary = 0;
		return 0;
	}

	sprintf(dirFileName, "%s%s", par->projectDir, par->bcFile);
	ini = iniparser_load(dirFileName);
	if(ini == NULL) 
		eprintf("cannot parse boundaryfile:");
	
	/**** Number of sections ****/
	nsec = 0;			// Number of sections in file
	i = 0;
	while(i < iniparser_getnsec(ini)){
		sprintf(str, "Boundary %i", nsec+1); 
		if (iniparser_find_entry(ini, str)) nsec++;
		i++;
	}
	/****************************/

	par->nBoundary = iniparser_getint(ini, "main:nboundary", -1);
	
	if (nsec > par->nBoundary){
		weprintf("inconsistencies in the boundary configuration file");
		weprintf("the last %i boundary are not being considered", 
			   nsec-par->nBoundary);
	}
	else if(nsec < par->nBoundary ){
		weprintf("inconsistencies in the boundary configuration file");
		weprintf("there are only %i boundary sections. adopting nboundary = %i",
			   nsec, par->nBoundary = nsec);
	}

	boundary = (Boundary *)calloc(par->nBoundary, sizeof(Boundary));
	if (!boundary) eprintf("allocation failure for boundary struct:");

	for(i = 0; i < par->nBoundary; i++){
		itoa(i+1, num, 10); 
		sprintf(SectionName, "Boundary %s:", num);
		
		sprintf(str, "%s%s", SectionName, "row");
		boundary[i].row = iniparser_getint(ini, str, -1);

		sprintf(str, "%s%s", SectionName, "col");
		boundary[i].col = iniparser_getint(ini, str, -1);

		sprintf(str, "%s%s", SectionName, "rowf");
		boundary[i].rowf = iniparser_getint(ini, str, boundary[i].row);

		sprintf(str, "%s%s", SectionName, "colf");
		boundary[i].colf = iniparser_getint(ini, str, boundary[i].col);
		
		sprintf(str, "%s%s", SectionName, "type");
		strcpy(type, _strupr(iniparser_getstring(ini, str, NULL)));
		
		sprintf(str, "%s%s", SectionName, "side");
		strcpy(side, _strupr(iniparser_getstring(ini, str, NULL)));
		
		sprintf(str, "%s%s", SectionName, "value");
		boundary[i].value = iniparser_getdouble(ini, str, NOBOUNDARY);

		if(boundary[i].row < 0 || boundary[i].row >= par->nrow){
			eprintf("row in [Boundary %i] outside of reservoir geometry", i+1); 
		}
		if(boundary[i].col < 0 || boundary[i].col >= par->ncol){
			eprintf("col in [Boundary %i] outside of reservoir geometry", i+1); 
		}
		if(boundary[i].rowf < 0 || boundary[i].rowf >= par->nrow){
			eprintf("rowf in [Boundary %i] outside of reservoir geometry", i+1); 
		}
		if(boundary[i].colf < 0 || boundary[i].colf >= par->ncol){
			eprintf("colf in [Boundary %i] outside of reservoir geometry", i+1); 
		}		

		if (strcmp(side, "E")  == 0  || strcmp(side, "EAST")  == 0)
			boundary[i].side = EAST;
		else if (strcmp(side, "W")  == 0  || strcmp(side, "WEST")  == 0)
				boundary[i].side = WEST;
		else if (strcmp(side, "N")  == 0  || strcmp(side, "NORTH") == 0)
				boundary[i].side = NORTH;
		else if (strcmp(side, "S")  == 0  || strcmp(side, "SOUTH") == 0)
				boundary[i].side = SOUTH;
		else 
			eprintf("side in [Boundary %i] is invalid", i+1); 
		
		if (strcmp(type, "PRESSURE_GRADIENT_ESPECIFIED") == 0 ||
			strcmp(type, "GRADIENT_ESPECIFIED") == 0 || 
			strcmp(type, "GRADIENT") == 0)
			boundary[i].type = PRESSURE_GRADIENT_ESPECIFIED;
		else if (strcmp(type, "PRESSURE_ESPECIFIED") == 0 ||
			     strcmp(type, "PRESSURE") == 0)
			 boundary[i].type = PRESSURE_ESPECIFIED;			
		else
			eprintf("type in [Boundary %i] is invalid", i+1); //0.0.2

		//Setting
		add_row = boundary[i].rowf > boundary[i].row  ? 1 : 
				  boundary[i].rowf == boundary[i].row ? 0 : -1;
		add_col = boundary[i].colf > boundary[i].col  ? 1 : 
				  boundary[i].colf == boundary[i].col ? 0 : -1;

		drow = abs(boundary[i].rowf - boundary[i].row);
		dcol = abs(boundary[i].colf - boundary[i].col);

		if ( drow!=0 && dcol!=0 && drow!=dcol ) {
			eprintf("[Boundary %i] are not vertical, horizontal or diagonal", i+1);
		}

		nblocks = drow==dcol ? drow : abs(drow-dcol); 

		for(j = 0; j <= nblocks ; j++){
			localIndex = 
				geom[boundary[i].row+add_row*j][boundary[i].col+add_col*j] - 1;
			
			grid[localIndex].boundary = 
				(int*) realloc(grid[localIndex].boundary, 
							   (++grid[localIndex].nBoundary)*sizeof(int));
			grid[localIndex].boundary[grid[localIndex].nBoundary-1] = i;

			switch(boundary[i].side){
				case NORTH:
					if (grid[localIndex].N != ISBONDARY){
						grid[grid[localIndex].N].S = ISBONDARY;
						grid[localIndex].N = ISBONDARY;
					}
					break;
				case SOUTH:
					if (grid[localIndex].S != ISBONDARY){
						grid[grid[localIndex].S].N = ISBONDARY;
						grid[localIndex].S = ISBONDARY;
					}
					break;
				case WEST:
					if (grid[localIndex].W != ISBONDARY){
						grid[grid[localIndex].W].E = ISBONDARY;
						grid[localIndex].W = ISBONDARY;
					}
					break;
				case EAST:
					if (grid[localIndex].E != ISBONDARY){
						grid[grid[localIndex].E].W = ISBONDARY;
						grid[localIndex].E = ISBONDARY;
					}
					break;
					
			}//*/
		}
	}

	iniparser_freedict(ini);

	unsetprogname();
	
	return boundary;
}
/*****************************************************************************/

void setBoundaryConditions(int *j, int *Map, double *Ax, double *b, 
									 double N, double E, double W, double S, 
									 Boundary* boundary, Block *block)
{
	int bool_E, bool_W, bool_N, bool_S, i, k, m, tmp;
	int diag_index = *j;

	bool_E = bool_W = bool_N = bool_S = TRUE;
	
	setprogname("boundary conditions");

	for(k=block->nBoundary-1; k >= 0; k--){

		i = block->boundary[k];
		
		if (boundary[i].type == PRESSURE_ESPECIFIED){
			if(boundary[i].side == EAST && bool_E == TRUE){
				Ax[Map[diag_index]] -= E;
				*b -= 2*E*boundary[i].value;
				bool_E = FALSE;
			}
			else if(boundary[i].side == WEST && bool_W == TRUE){
				Ax[Map[diag_index]] -= W;
				*b -= 2*W*boundary[i].value;
				bool_W = FALSE;
			}
			else if(boundary[i].side == NORTH && bool_N == TRUE){
				Ax[Map[diag_index]] -= N;
				*b -= 2*N*boundary[i].value;
				bool_N = FALSE;
			}
			else if(boundary[i].side == SOUTH && bool_S == TRUE){
				Ax[Map[diag_index]] -= S;
				*b -= 2*S*boundary[i].value;
				bool_S = FALSE;
			}
			else{
				weprintf(
"there are more than one boundary conditions to the same\
                         block. Forgeting [Boundary %i] for block (%d,%d).", 
					i, block->row, block->col);	
				for(m=i; m < block->nBoundary-1; m++){
					tmp = block->boundary[m+1];
					block->boundary[m+1] = block->boundary[m];
					block->boundary[m] = tmp;
				}
				block->nBoundary--;
			}
		}
	
		else if (boundary[i].type == PRESSURE_GRADIENT_ESPECIFIED){
			if(boundary[i].side == EAST && bool_E == TRUE){
				Ax[Map[diag_index]] += E;
				*b += block->dx*E*boundary[i].value;
				bool_E = FALSE;
			}
			else if(boundary[i].side == WEST && bool_W == TRUE){
				Ax[Map[diag_index]] += W;
				*b -= block->dx*W*boundary[i].value;
				bool_W = FALSE;
			}
			else if(boundary[i].side == NORTH && bool_N == TRUE){
				Ax[Map[diag_index]] += N;
				*b -= block->dy*N*boundary[i].value;
				bool_N = FALSE;
			}
			else if(boundary[i].side == SOUTH && bool_S == TRUE){
				Ax[Map[diag_index]] += S;
				*b -= block->dy*S*boundary[i].value;
				bool_S = FALSE;
			}
			else{
				weprintf(
"there are more than one boundary conditions to the same\
                         block. Forgeting [Boundary %i] for block (%d,%d).", 
				  i, block->row, block->col);	
				for(m=i; m < block->nBoundary-1; m++){
					tmp = block->boundary[m+1];
					block->boundary[m+1] = block->boundary[m];
					block->boundary[m] = tmp;
				}
				block->nBoundary--;
			}
		}
	}

	(*j)++;
	if(E > 0.0 && bool_E == TRUE){
		Ax[Map[*j]] = E;
		(*j)++;
	}
	if(W > 0.0 && bool_W == TRUE){
		Ax[Map[*j]] = W;
		(*j)++;
	}
	if(N > 0.0 && bool_N == TRUE){
		Ax[Map[*j]] = N;
		(*j)++;
	}
	if(S > 0.0 && bool_S == TRUE){
		Ax[Map[*j]] = S;
		(*j)++;
	}

	unsetprogname();
}