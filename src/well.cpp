#include "well.h"

char *key(std::string k, int num)
{
	std::stringstream key_tmp;
	char key[LENGTHSP];

	key_tmp << "well "<< num << ":" << k;
	strcpy(key, key_tmp.str().c_str());

	return key;
}
/*****************************************************************************/ 

Well::Well(dictionary *ini, int num)
{
	row = iniparser_getint(ini, key("row", num), -1);

	col = iniparser_getint(ini, key("col", num), -1);

	dx = iniparser_getdouble(ini, key("dx", num), 0); 

	dy = iniparser_getdouble(ini, key("dy", num), 0); 

	rw = iniparser_getdouble(ini, key("rw", num), -1);

	h  = iniparser_getdouble(ini, key("h", num), -1);

	s  = iniparser_getdouble(ini, key("s", num), 0);

	pwf = iniparser_getdouble(ini, key("pwf", num), 0);

	qw = iniparser_getdouble(ini, key("qw", num), 0);

	std::string type_ (iniparser_getstring(ini, key("type", num), NULL));  

	np = 0; 

	if (type_ == "RATE_ESPECIFIED" == 0 ||	type_ == "RATE")
		type = RATE_ESPECIFIED;
	else if (type_ == "PRESSURE_ESPECIFIED" || type_ == "PRESSURE")
		type = PRESSURE_ESPECIFIED;			
	else
		eprintf("type in [Well %i] is invalid", num);
}
/*****************************************************************************/ 

double Well::reqPeaceman(Block block)
{
	double req;
	
	req = 0.28*sqrt( sqrt(block.ky/block.kx) * block.dx * block.dx 
		+ sqrt(block.kx/block.ky) * block.dy * block.dy )
		/ ( pow(block.kx/block.ky, 0.25) + pow(block.ky/block.kx, 0.25) );

	return req;
}
/*****************************************************************************/

double Well::reqDing(Block *grid, int localIndex)
{
	double dxi, dyi, dxip, dxim, dyjp, dyjm;
	double alphaN, alphaS, alphaE, alphaW;

	double req = reqPeaceman(grid[localIndex]);

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

double Well::reqAbouKassemAziz(Block *grid, int localIndex)
{
	double req;
	
	req = 0.0;

	return req;
}
/*****************************************************************************/