/*--------------------------------------------------------------*/
/* 								*/
/*		compute_rain_stored				*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_rain_stored -  rain intercepted by stratum and 	*/
/*		updates rain flaaing through.			*/
/*								*/
/*	SYNOPSIS						*/
/*	double	compute_rain_stored(				*/
/*			double	*,				*/
/*			struct	canopy_strata_object	*);	*/
/*								*/
/*	returns:						*/
/*	rain_storage (m) - amount of rain left in canopy	*/
/*								*/
/*	also updates the following state variables:		*/
/*								*/
/*	*rain (m) - rain falling on next stratum		*/
/*	stratum[0].potential_evaporation (m/day) - amount	*/
/*		of water that can be still evaporated that day	*/
/*	OPTIONS							*/
/*	*rain (m) - amount of rain on stratum for whole day	*/
/*	canopy_strata_object - state of canopy strata		*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	Updates amount of rain stored after using supplied 	*/
/*	daytime potential evaporation,  rainfall 			*/
/*	to evaporate and input respectively to the storage.	*/
/*	Evaporation from the storage is performed first.	*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

double	compute_rain_stored(
							int	verbose_flag,
							double	*rain,
							struct	canopy_strata_object	*stratum)
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration				*/
	/*--------------------------------------------------------------*/
	double	compute_potential_rain_interception(
		int	,
		double	,
		struct	canopy_strata_object *);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	double	potential_interception;
	double	potential_evaporation;
	double	rain_storage;
	double	storage_evaporated;
	double	potential_interception_evaporated;
	double	throughfall;
	/*--------------------------------------------------------------*/
	/*	Transfer potential evaporation to a working variable.	*/
	/*--------------------------------------------------------------*/
	potential_evaporation = stratum[0].potential_evaporation;
	/*--------------------------------------------------------------*/
	/*	Compute amount potentially intercepted.			*/
	/*	m = m2PlANT / m2ground *  ( (kg  / m2 * day * m2PLANT )	*/
	/*		* ( 1 m3 H20 / 1000 kg H20 )			*/
	/*--------------------------------------------------------------*/
	potential_interception = compute_potential_rain_interception(
		verbose_flag,
		*rain,
		stratum );
	if( verbose_flag >2)
		printf("%8.6f %8.6f ",*rain, stratum[0].rain_stored);
	if( verbose_flag >2)
		printf("%8.6f ",potential_interception);
	if (verbose_flag == -5) {
		printf("\n          RAIN STORED:rain=%lf rainstor=%lf potint=%lf potevap=%lf",
			   *rain,
			   stratum[0].rain_stored,
			   potential_interception,
			   potential_evaporation);
	}

	/*--------------------------------------------------------------*/
	/*	Add amount not potentially intercepted to throughfall.	*/
	/*	m = m - m 						*/
	/*--------------------------------------------------------------*/
	throughfall  = *rain - potential_interception;
	/*--------------------------------------------------------------*/
	/*	Compute amount of storage evaporated.			*/
	/*	m = m							*/
	/*--------------------------------------------------------------*/
	storage_evaporated =	min(potential_evaporation,stratum[0].rain_stored);


    //printf("potential_evaporation=%.1f rain_stored=%.1f storage_evaporated=%.1f\n",
    //       potential_evaporation*1000,
    //       stratum[0].rain_stored*1000,
    //       storage_evaporated*1000);


	/*--------------------------------------------------------------*/
	/*	Update amount of rain in storage after evaporation.	*/
	/*	m = m							*/
	/*--------------------------------------------------------------*/
	stratum[0].rain_stored  -= storage_evaporated;
	/*--------------------------------------------------------------*/
	/*	Update potentail evaporation.				*/
	/*	m = m							*/
	/*--------------------------------------------------------------*/
	potential_evaporation -= storage_evaporated;
	/*--------------------------------------------------------------*/
	/*	Compute amount of potential interception evaporated.	*/
	/* 	m = m							*/
	/*--------------------------------------------------------------*/
	potential_interception_evaporated  =
		min( potential_evaporation, potential_interception);
	/*--------------------------------------------------------------*/
	/*	Update the potential interception after evaporation.	*/
	/* 	m -= m							*/
	/*--------------------------------------------------------------*/
	potential_interception -= potential_interception_evaporated;
	/*--------------------------------------------------------------*/
	/*	Compute amount of evaporation that happened.		*/
	/* 	m = m + m						*/
	/*--------------------------------------------------------------*/
	stratum[0].evaporation =  storage_evaporated +
		potential_interception_evaporated;

    //printf("potential_interception_evaporated=%.1f evaporation=%.1f\n",
    //       potential_interception_evaporated*1000,
    //       stratum[0].evaporation*1000);

	/*--------------------------------------------------------------*/
	/*	Adjust the amount of remaining potential evaporation	*/
	/*--------------------------------------------------------------*/
	stratum[0].potential_evaporation -= stratum[0].evaporation;
	if( verbose_flag >2)
		printf("%8.6f ",stratum[0].evaporation);
	if( verbose_flag >2)
		printf("%8.6f ",potential_interception_evaporated);
	if( verbose_flag >2)
		printf("%8.6f ",storage_evaporated);
	/*--------------------------------------------------------------*/
	/*  Update rain storage.										*/
	/*	m = m							*/
	/*--------------------------------------------------------------*/
    double max_rain_storage = stratum[0].epv.all_pai_when_red //NREN 20180804 in update_phenology, if there is not beetle attack, the all_api_b == all_api
            * stratum[0].defaults[0][0].specific_rain_capacity;
    double extra_rain_storage = max(stratum[0].rain_stored + potential_interception - max_rain_storage,0);
    rain_storage = stratum[0].rain_stored + potential_interception - extra_rain_storage;
	if( verbose_flag > 2)
		printf("%8.6f ",rain_storage);
	/*--------------------------------------------------------------*/
	/*	Update rain throughfall.									*/
	/*	m += m							*/
	/*--------------------------------------------------------------*/
    //throughfall += max(potential_interception
    //	- (rain_storage - stratum[0].rain_stored),0);
    throughfall += extra_rain_storage;
	if( verbose_flag > 2)
		printf("%8.6f ",throughfall);

    if (verbose_flag == -5) {
		printf("                    RAIN STORED END:rainstor=%lf throughfall=%lf oldrainstor=%lf evap=%lf potevap=%lf",
			   rain_storage,
			   throughfall,
			   stratum[0].rain_stored,
			   stratum[0].evaporation,
			   stratum[0].potential_evaporation);
    }

	*rain = throughfall;

	return( rain_storage);
} /*end compute_rain_stored*/
