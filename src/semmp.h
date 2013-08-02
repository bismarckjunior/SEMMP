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

typedef struct {
	double iniTime;
	double cf;
	double cphi;
	double dpprops;
	double p0;
	double z0;
	double invBSC;
	double rhoSC;
	double dt;
	double dtMultiplier;
	
	int dtLogScale;
	int ncol;
	int nrow;
	int nlay;
	int nSteps;
	int outPressureSteps;   // steps for print pressures in outPressureFile
	int nBlocks;
	int nfprop;
	int nwells;
	int nMatrix;
	int displaySteps;		// steps for display on terminal
	int reportSteps;		// steps for print report of blocks (in reportFile)
	int nDisplayBlocks;		// number of blocks to display on terminal
	int nReportBlocks;		// number of blocks to the report (in outputFile)
	int nBoundary;          // number of boundary contitions (in bcFile)
	int isCylindrical;		// bool for cylindrical coordenates

	char reportFile[LENGTHFN];  
	char outPressureFile[LENGTHFN];
	char projectDir[LENGTHFN];
	char projectName[LENGTHFN];  //0.0.2
	char iniFile[LENGTHFN];
	char geometryFile[LENGTHFN];
	char modifiedBlocksFile[LENGTHFN]; 
	char modTransmissibilityFile[LENGTHFN]; //0.0.3
	char porosityFile[LENGTHFN];
	char kxFile[LENGTHFN];
	char kyFile[LENGTHFN];
	char thicknessFile[LENGTHFN];
	char zTopFile[LENGTHFN];
	char dxFile[LENGTHFN];
	char dyFile[LENGTHFN];
	char fluidPropFile[LENGTHFN];
	char wellsFile[LENGTHFN];   
	char outputFile[LENGTHFN];  
	char bcFile[LENGTHFN];      //boundary conditions file
} Parameters;

typedef struct {
	/* constant part of the transmissibilities */
	double Gxph; /* Gx_i+1/2 */ 
	double Gxmh; /* Gx_i-1/2 */
	double Gyph; /* Gy_j+1/2 */
	double Gymh; /* Gy_j-1/2 */
	
	double dx;
	double dy;
	double dz;
	double z;
	double Vac;  //V_b*alpha_c
	double phi;
	double kx;
	double ky;

	int isWellBlock;

	int isBoundary;
	int* boundary;
	int nBoundary;

	int row;
	int col;
	
	int H; /* Here */
	int E; /* East */
	int W; /* West */
	int N; /* North */
	int S; /* South */
	int A; /* Above */
	int B; /* Below */
} Block;

typedef struct {
	int		row;
	int		col;	// column
	int		lay;	// layer
	int		type;	// pressure or flow rated specified
	double	dx;		// x distance to the block's SW corner
	double	dy;		// y distance to the block's SW corner
	double	rw;		// wellbore radius
	double	s;		// skin factor
	double	h;		// height
	double	gw;		// constant part of well index
	double	qw;		// flow rate
	double	pwf;	// bottom hole pressure
	double  np;		// cumulative oil production
} Well;

typedef struct {
	int row;
	int col;		// column
	int lay;		// layer
	int p;			// bool for display pressure
	int qw;			// bool for display flow rate
	int np;			// cumulative oil production	
	int pwf;		// bool for display bottom hole pressure
	int localIndex;	// block index in grid
} Out;

typedef struct{
	int row;		//row of first block or unique block 
	int col;		//col of first block or unique block 
	int lay;		//layer of first block or unique block 
	int rowf;		//row of last block 
	int colf;		//column of last block
	int layf;		//layer of last block
	int type;		//type of boundary condition
	char side;      //side: NORTH, EAST, SOUTH, WEST
	double value;	//value of boundary condition
} Boundary;


/********************************************/
/******* memory allocation utilities ********/
int *iVector(int);
int **iMatrix(int, int);
int ***iMatrix3D(int, int, int);
void freeiVector(int *);
void freeiMatrix(int **, int, int);
double *rVector(int);
double **rMatrix(int, int);
void freerVector(double *);
void freerMatrix(double **, int, int);
/********************************************/

/********************************************/
/******* file input/output utilities ********/
void rTableFile(char [], double **, int, int);
void iTableFile(char [], int **, int, int);
void readIniFile(int, char *[], Parameters *);
Block *readAndSetGeometry(Parameters *, int **);
double **readFluidProperties(Parameters *);
/********************************************/

void setInitialPressure(Parameters *, Block *, Well *, double **, double *, 
						double *, double *);

double normDelta(double *, double *, int);

double getTableVal(double **, int, double, int, double, int);

void setMatrixStructure(Parameters *, Block *, int *, int *, double *,int *, 
						int *, double *, int *);

void setTransmissibilityMatrix(Parameters *, Block *, Well *, double **, 
							   double *, double *, double *, double *, int *, 
							   Boundary*);

void writePressureFile(Parameters *, int, Block *, double *);

void setCartesianTransmissibilities(Parameters *, Block *, Boundary *); 

void setCylindricalTransmissibilities(Parameters *, Block *, Boundary *, 
									  Well *);

Well *readAndSetWellParameters(Parameters *);

void setWells(Parameters *, Block *, int **, Well *);

double reqAbouKassemAziz(Block *, int);

double reqPeaceman(Block *, int);

double reqDing(Block*, int, double);

void modifyGeometry(int**, Parameters*);

void modifyTransmissibilities(int **, Block *, Parameters *);

Boundary *readAndSetBoundaryConditions(Parameters*, int**, Block*);

void setBoundaryConditions(int*, int*, double*, double*, double, 
						   double, double, double, Boundary*, Block*);

Out *readAndSetOuts(Parameters*, Block*, int**, char*, double*);

void *readAndSetMainOut(Parameters*, double*);

void writeOuts(FILE*, int, Block*, Well *, Out*, double, double*, double*);

void writeTableHeader(FILE*, Parameters*, Out*, int, double*);

void header(FILE *, Parameters *);