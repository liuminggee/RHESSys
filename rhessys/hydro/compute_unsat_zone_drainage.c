/*--------------------------------------------------------------*/
/* 								*/
/*		compute_unsat_zone_drainage			*/
/*								*/
/*	NAME							*/
/*	compute_unsat_zone_drainage - estimates vertical 	*/
/*		drainage from the unsat to sat zone.		*/
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	compute_unsat_zone_drainage(				*/
/*				int	,			*/	
/*				double	,			*/
/*				double	,			*/
/*				double	,			*/
/*				double	)			*/
/*								*/
/*	returns:						*/
/*	unsat_zone_drainage - (m water) drainage from unsat to	*/
/*			sat zone.				*/
/*								*/
/*	OPTIONS							*/
/*	int verbose_flag 					*/
/*	double	m - Ksat decay parameter			*/
/*	double	z - (m) depth to the water table		*/
/*	double Ksat_0 - (m/day) sat. hydraulic conductivity	*/
/*				at the surface.			*/
/*	double	potential_drainage - (m water) difference	*/
/*		between current unsat_zone_storage and current	*/
/*		field capacity.					*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	This routine is designed to estimate unsaturated zone	*/
/*	drainage from soils with "TOPMODEL" properties under	*/
/*	the assumptions that:					*/
/*								*/
/*	i) a field capacity for the unsat zone is known.	*/
/*	ii) the drainage will bethe minimum of drainage to	*/ 
/*		field capacity or the amount resulting from	*/
/*		the current Ksat at the water table.		*/
/*								*/
/*	iii) the change in sat_deficit due to drainage does	*/
/*		not significantly affect field capacity.	*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	To avoid confilcts with modified versions of TOPMODEL	*/
/*	we call a Ksat(z) curve.				*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double	compute_unsat_zone_drainage(
									int	verbose_flag,
									int	curve,
									double	p2,
									double	S,
									double	m,
									double	z,
									double	Ksat_0,
									double	potential_drainage )
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration									*/
	/*--------------------------------------------------------------*/
	double	Ksat_z_curve(
		int,
		double,
		double,
		double);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	double	Ksat;
	double	Ksat1;
	double	Ksat2;
	double	unsat_zone_drainage;
	/*--------------------------------------------------------------*/
	/*	Compute Ksat 		(m water/day)			*/
	/*								*/
	/*	This undersestimates drainage since Ksat will increase	*/
	/*	as drainage progresses.					*/
	/*--------------------------------------------------------------*/
	Ksat1  = Ksat_z_curve(
		verbose_flag,
		m,
		z,
		Ksat_0);

	if (curve == 1)
        Ksat2 = Ksat_0 * pow(S,(2./p2+3.));
	else
        Ksat2 = Ksat_0 * sqrt(S) * pow(1. - pow(1.-pow(S,1./p2),p2),2.);

	Ksat = min(Ksat1,Ksat2);

	/*--------------------------------------------------------------*/
	/*	Compute unsat zone drainage.				*/
	/*--------------------------------------------------------------*/
	unsat_zone_drainage = max(min( potential_drainage,	Ksat ),0);

	return(unsat_zone_drainage);
} /*compute_unsat_zone_drainage*/

double	compute_unsat_zone_drainage_patch(
                                    int	verbose_flag,
                                    struct patch_object *patch,
                                    float num_per_day,
                                    int rz,
                                    int topmodel,
                                    int potential_drainage_use_rootzone)
{
    /*--------------------------------------------------------------*/
    /*	Local function declaration									*/
    /*--------------------------------------------------------------*/
    double	Ksat_z_curve(
        int,
        double,
        double,
        double);
    /*--------------------------------------------------------------*/
    /*	Local variable definition.									*/
    /*--------------------------------------------------------------*/
    double	Ksat;
    double	Ksat1;
    double	Ksat2;
    double	unsat_zone_drainage;
    /*--------------------------------------------------------------*/
    /*	Compute Ksat 		(m water/day)			*/
    /*								*/
    /*	This undersestimates drainage since Ksat will increase	*/
    /*	as drainage progresses.					*/
    /*--------------------------------------------------------------*/
    struct soil_default *psoildef = patch->soil_defaults[0];
    int	curve = psoildef[0].theta_psi_curve;
    double	p2 = psoildef[0].pore_size_index;
    double	m = psoildef[0].mz_v;
    double	S;
    double	z;
    double	potential_drainage;
    double	Ksat_0;

    if (rz == 1) {
        S = patch->rootzone.S;
        z = patch->rootzone.depth;
        potential_drainage = patch->rz_storage - patch->rootzone.field_capacity;
    } else {
        S = patch->S;
        z = patch->sat_deficit_z;
        if (potential_drainage_use_rootzone)
            potential_drainage = patch->rz_storage - patch->rootzone.field_capacity;
        else potential_drainage = patch->unsat_storage - patch->field_capacity;
    }
    if (topmodel == 1) Ksat_0 = psoildef[0].Ksat_0 / num_per_day / 2.;
    else Ksat_0 = psoildef[0].Ksat_0_v / num_per_day / 2.;


    Ksat1  = Ksat_z_curve(
        verbose_flag,
        m,
        z,
        Ksat_0);

    if (curve == 1)
        Ksat2 = Ksat_0 * pow(S,(2./p2+3.));
    else
        Ksat2 = Ksat_0 * sqrt(S) * pow(1. - pow(1.-pow(S,1./p2),p2),2.);

    Ksat = min(Ksat1,Ksat2);

    /*--------------------------------------------------------------*/
    /*	Compute unsat zone drainage.				*/
    /*--------------------------------------------------------------*/
    unsat_zone_drainage = max(min( potential_drainage,	Ksat ),0);

    return(unsat_zone_drainage);
} /*compute_unsat_zone_drainage*/
