#include "iosemmp.h"




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

template <class T>
T** readTableFile(char filein[], int nrow, int ncol)
{
	float tmp;
	int col, row;
	FILE *mf;
	T **table = rMatrix(nrow, ncol);
		
	if ( (mf = fopen(filein, "r")) == NULL ){
		eprintf("opening table file %s:", filein);
	}
	
	for (row = 0; row < nrow; row++) {
		for (col = 0; col < ncol; col++) {
			
			fscanf(mf, "%f", &tmp);
			table[row][col] = (T)tmp;
		}
	}
	fclose(mf);
	
	return table;
}
/*****************************************************************************/

double** readFluidProperties(Parameters* par)
{
	int i;
	char dirFileName[LENGTHFN];
	double** fluidProps;
	
	setprogname("read fluid properties");

	par->invBSC = 1.0; 

	/* for pressure, formation-volume-factor, viscosity, especific gravity */
	sprintf(dirFileName, "%s%s", par->projectDir, par->fluidPropFile);
	fluidProps = readTableFile<double>(dirFileName, par->nfprop, NCOLPROPS);

	for(i = 0; i < par->nfprop; i++){
		fluidProps[i][INV_FVF] = 1.0 / fluidProps[i][FVF];
		fluidProps[i][INV_FVFVISC] = 
			fluidProps[i][INV_FVF] / fluidProps[i][VISC];
	}

	unsetprogname();

	return fluidProps;
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

void readAndSetMainOut(Parameters *par, double *cumulativeProd)
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

Block *readAndSetGeometry(Parameters *par, int **geom)
{
	int nrow, ncol, col, row, localIndex, aux, i = 0, j = 0;
	int isGeoFile, isPhiFile, line = 0;
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
		sprintf(dirFileName, "%s%s", par->projectDir, par->porosityFile);
		porositytmp = readTableFile<double>(dirFileName, nrow, ncol);
	}
	else {
		phi = atof(par->porosityFile); 
	}
	if (isKxFile = !atof(par->kxFile)){
		sprintf(dirFileName, "%s%s", par->projectDir, par->kxFile);
		kxtmp = readTableFile<double>(dirFileName, nrow, ncol);
	}
	else {
		kx = atof(par->kxFile);
	}
	if (isKyFile = !atof(par->kyFile)){
		sprintf(dirFileName, "%s%s", par->projectDir, par->kyFile);
		kytmp = readTableFile<double>(dirFileName, nrow, ncol);
	}
	else {
		ky = atof(par->kyFile);
	}
	if (isDzFile = !atof(par->thicknessFile)){
		sprintf(dirFileName, "%s%s", par->projectDir, par->thicknessFile);
		thicknesstmp = readTableFile<double>(dirFileName, nrow, ncol);
	}
	else {
		dz = atof(par->thicknessFile);
	}
	if (isZtopFile = !atof(par->zTopFile)){
		sprintf(dirFileName, "%s%s", par->projectDir, par->zTopFile);
		zToptmp = readTableFile<double>(dirFileName, nrow, ncol);
	}
	else {
		ztop = atof(par->zTopFile);
	}
	if (isDxFile = !atof(par->dxFile)){
		sprintf(dirFileName, "%s%s", par->projectDir, par->dxFile);
		dxtmp = readTableFile<double>(dirFileName, 1, ncol);
	}
	else {
		dx = atof(par->dxFile);
	}
	if (isDyFile = !atof(par->dyFile)){
		sprintf(dirFileName, "%s%s", par->projectDir, par->dyFile);
		dytmp = readTableFile<double>(dirFileName, nrow, 1);
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
		
		wells[i] = Well(ini, i+1);

		if(wells[i].row < 0 || wells[i].row >= par->nrow){ 
			eprintf("row in [Well %i] outside of reservoir geometry", i+1); 
		}

		if(wells[i].col < 0 || wells[i].col >= par->ncol){
			eprintf("col in [Well %i] outside of reservoir geometry", i+1);
		}
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
		par->nrow, par->ncol, par->iniTime, par->nSteps, par->dt, par->cf, 
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
