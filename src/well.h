#ifndef WELL_H
#define WELL_H

#include <math.h>
#include <string>
#include <sstream>
#include "eprintf.h"
#include "block.h"
#include "parameters.h"
extern "C" {
#include <iniparser.h> 
} 


/* well types */
#define RATE_ESPECIFIED		1    
#define PRESSURE_ESPECIFIED	2	
#define DING				0.25 //Ding's constant for req

#define PI				(double)3.14159265358979
#define BETAC			(double)1.127

#define LENGTHSP		25		 //section properties length


class Well {
public:
	int		row;	// row block
	int		col;	// col block
	int		type;	// pressure or flow rated specified
	double	dx;		// x offset from center
	double	dy;		// y offset from center
	double	rw;		// wellbore radius
	double	s;		// skin factor
	double	h;		// height (open to flow)
	double	gw;		// constant part of well index
	double	qw;		// flow rate
	double	pwf;	// bottom hole pressure
	double  np;		// cumulative oil production


	/** 
	 * Constructor 
	 *
	 * @param ini dictionary to well file
	 * @param num well number section
	 */
	Well(dictionary* ini, int num);

	/**
	 * Set Well
	 *
	 * @param block wellblock
	 */
	void setWell(Block block);

	/** 
	 * Abou Kassem & Aziz's well model
	 *
	 * Abou-Kassem J. H. and Aziz K., Analytical well models for reservoir 
	 * simulation, SPEJ (Aug 1985), pp. 573-79. 
	 *
	 * @param grid grid
	 * @param localIndex index of well in grid
	 * @return equivalent wellblock radius
	 */
	double reqAbouKassemAziz(Block *grid, int localIndex);

	/** 
	 * Peaceman's well model
	 *
	 * Peaceman D. W., Interpretation of Well-Block pressures in numerical 
	 * reservoir simulation with nonsquare grid blocks and anisotropic 
	 * permeability,SPEJ (June 1978), pp. 531-43. 
	 *
	 * @param block wellblock
	 * @return equivalent wellblock radius
	 */
	double reqPeaceman(Block block);

	/** 
	 * Ding et al's well model
	 *
	 * Yu Ding, Gérard Renard, Luce Weill, Representation of Wells in Numerical 
	 * Reservoir Simulation, SPE (February 1998), pp. 18-23.
	 *
	 * @param grid grid
	 * @param localIndex index of well in grid
	 * @return equivalent wellblock radius
	 */
	double reqDing(Block *grid, int localIndex);
};


char *key(std::string k, int num);


#endif