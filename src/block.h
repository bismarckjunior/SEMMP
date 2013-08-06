#ifndef BLOCK_H
#define BLOCK_H


typedef struct {
	int row;		///< block row
	int col;		///< block column 
	int lay;		///< block layer 
	int p;			///< bool for display pressure
	int qw;			///< bool for display flow rate
	int np;			///< cumulative oil production	
	int pwf;		///< bool for display bottom hole pressure
	int localIndex;	///< block index in grid
} Out;


typedef struct {
	/*block indexes*/
	int H;				///< Here
	int E;				///< East
	int W;				///< West
	int N;				///< North
	int S;				///< South
	int A;				///< Above
	int B;				///< Below
	
	/*well block*/
	int isWellBlock;	///< is well block?

	/*boundaries conditions*/
	int isBoundary;		///< is boundary?
	int* boundary;		///< indexes to Boundary vector
	int nBoundary;		///< number of boundaries

	/*block properties*/
	int row;			///< block row
	int col;			///< block column
	int lay;			///< block layer
	double dx;			///< control volume dimension along x direction [ft]
	double dy;			///< control volume dimension along y direction [ft]
	double dz;			///< control volume dimension along z direction [ft]
	double z;			///< block center depth [ft]
	double Vac;			///< bulk volume/ALPHAC
	double phi;			///< porosity
	double kx;			///< permeability in the x direction [mD]
	double ky;			///< permeability in the y direction [mD]
	double kz;			///< permeability in the z direction [mD]

	/* constant part of the transmissibilities */
	double Gxph;		///< Gx_(i+1/2)  
	double Gxmh;		///< Gx_(i-1/2) 
	double Gyph;		///< Gy_(j+1/2) 
	double Gymh;		///< Gy_(j-1/2)
	double Gzph;		///< Gz_(k+1/2) 
	double Gzmh;		///< Gz_(k-1/2)
	
} Block;

#endif