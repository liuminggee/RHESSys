/*--------------------------------------------------------------*/
/* 								*/
/*		compute_potential_rain__interception		*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_potential_rain_interception  - computes amount 	*/
/*		rain that can be intercepted by the canopy.	*/
/*								*/
/*	SYNOPSIS						*/
/*	compute_potential_rain_interception(   			*/
/*			int	,				*/
/*			double	,				*/
/*			struct	canopy_strata_object	*);	*/
/*								*/
/*	returns:						*/
/*	potential_interception (m) - amount of rain that can be */
/*		intercepted by the canopy.			*/
/*								*/
/*	OPTIONS							*/
/*	rain (m) - amount of rain on stratum for whole day	*/
/*	canopy_strata_object - state of canopy strata		*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	Computes the amount of rain that is interceptible	*/
/*	given the rain during the day and the current canopy	*/
/*	start of day storage.  Note that this interceptible	*/
/*	rain may be later dripped or evaporated if the code 	*/
/*	choses in compute_rain_stored.				*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

double	compute_potential_rain_interception(
											int	verbose_flag,
											double	rain,
											struct	canopy_strata_object	*stratum)
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration				*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
    double	potential_interception = 0;
	double	interception_coef;
	/*--------------------------------------------------------------*/
	/*	Compute amount potentially intercepted.			*/
	/*	m = m2PlANT / m2ground *  ( (kg  / m2 * day * m2PLANT )	*/
	/*		* ( 1 m3 H20 / 1000 kg H20 )			*/
	/*	limit incoming rain by gap_fraction			*/
	/*--------------------------------------------------------------*/
	interception_coef = 1-stratum[0].gap_fraction;

    if (stratum[0].defaults[0][0].epc.veg_type != NON_VEG) {
        //07132023LML confine maximum intercepted water to 4 mm according to https://en.wikipedia.org/wiki/Interception_(water)
        double maxint = min(0.004,//	stratum[0].epv.proj_pai
                            stratum[0].epv.proj_pai_when_red //NREN 20180804 in update_phenology, if there is no beetle attack,the proj_pai_when_red = proj_pai
                             * stratum[0].defaults[0][0].specific_rain_capacity);
        potential_interception = min(interception_coef * rain
                                     ,maxint - stratum[0].rain_stored);
    } //else {
        //10052023LML should be zero
        //potential_interception = min(rain, (
        //	stratum[0].defaults[0][0].specific_rain_capacity
        //	- stratum[0].rain_stored));
    //}

    //07132023LML confine maximum intercepted water to 4 mm according to https://en.wikipedia.org/wiki/Interception_(water)
	potential_interception = max(potential_interception, 0.0);

    //printf(" rain:%lf rain_stored:%lf potential_interception:%lf proj_lai:%lf pai:%lf rain_capacity:%lf\n"
    //       ,rain*1000
    //       ,stratum[0].rain_stored*1000
    //       ,potential_interception*1000
    //       ,stratum[0].epv.proj_lai
    //       ,stratum[0].epv.proj_pai_when_red
    //       ,stratum[0].defaults[0][0].specific_rain_capacity*1000);

	return( potential_interception );
} /*end compute_potential_rain_interception */
