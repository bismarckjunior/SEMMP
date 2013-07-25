#include "grid.h"


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

void setCylindricalTransmissibilities(Parameters *par, Block *grid, Boundary 
								 *boundary, Well *wells)
{
	int i, j, n, bcIndex;
	int bool_E, bool_W;
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
		//readTableFile(drFileName, drtmp, 1, par->ncol);
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


void setWells(Parameters *par, Block *grid, int **geom, Well *wells)
{
	int i, localIndex;
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

			req = wells[i].reqPeaceman(grid[localIndex]);
			//req = 0.14*sqrt(pow(grid[localIndex].dx,2.0)+pow(grid[localIndex].dy,2.0));
			//req = reqDing(grid, localIndex, req);
		}
		else if(grid[localIndex].N != grid[localIndex].S 
			&&  grid[localIndex].E != grid[localIndex].W) {
			
			req = wells[i].reqAbouKassemAziz(grid, localIndex);
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
