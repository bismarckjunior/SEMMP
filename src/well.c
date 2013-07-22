#include "well.h"

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

