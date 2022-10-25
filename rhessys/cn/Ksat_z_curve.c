/*--------------------------------------------------------------*/
/*								*/
/*		Ksat_z_curve					*/
/*								*/
/*	Ksat_z_curve						*/
/*								*/
/*	NAME							*/
/*	Ksat_z_curve						*/
/*								*/
/*	SYNOPSIS						*/
/*	double	Ksat_z_curve(					*/
/*			int,					*/
/*			double, 	 			*/
/*			double )				*/
/*								*/
/*	returns: 						*/
/*	Ksat - (m/day) saturated hydraulic conductivity		*/
/*								*/
/*	OPTIONS							*/
/*	int verbose_flag					*/
/*	double m - TOPMODEL Ksat decay parameter		*/	
/*	double z - depth from surface				*/	
/*	Ksat_0 - (m/day) - saturated hydraulic conductivity 	*/
/*			at the surface 				*/ 
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	Uses a negative exponential decay of conductivity with 	*/
/*	depth.  						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

double	Ksat_z_curve(	
					 int	verbose_flag,
					 double	m,
					 double	z,
					 double	Ksat_0 )
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration									*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Specify curve						*/
	/*--------------------------------------------------------------*/
    if ( m > ZERO )
        return Ksat_0 * exp( -1. * z / m );
    else
        return Ksat_0;
} /*end Ksat_z_curve*/
