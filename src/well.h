#ifndef WELL_H
#define WELL_H

#include "parameters.h"

/* well types */
//#define RATE_ESPECIFIED		1    
//#define PRESSURE_ESPECIFIED	2	
//#define DING				0.25 //Ding's constant for req

//#define PI				(double)3.14159265358979
//#define BETAC			(double)1.127

typedef struct {
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
} Well ;

/*Sets well model*/
//Well(Block *);

/*Sets constant part of well index*/
//void setGw();

/*Peaceman's well model*/
double reqPeaceman(Block *, int);

/*Ding et al's well model*/
double reqDing(Block*, int, double);

/*Abou Kassem & Aziz's well model*/
double reqAbouKassemAziz(Block *, int);

Well *readAndSetWellParameters(Parameters *);

void setWells(Parameters *, Block *, int **, Well *);

#endif