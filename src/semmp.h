
/* well types */
//#define RATE_ESPECIFIED		1    //well type
//#define PRESSURE_ESPECIFIED	2	 //well type
//#define DING				0.25 //Ding's coupling model

/* constants */
//#define BETAC			(double)1.127
//#define ALPHAC			(double)5.614583
#define GAMMAC			(double)0.21584e-3
#define GRAV			(double)32.174
//#define EPSILON			(double)1e-8
//#define PI				(double)3.14159265358979

#define BARREL_TO_CF	(double)5.614584


#include "block.h"
#include "well.h"
#include "boundary.h"
#include "parameters.h"
#include "iosemmp.h"
#include "grid.h"
#include "memory.h"
#include "umfpack.h"
#include "eprintf.h"






