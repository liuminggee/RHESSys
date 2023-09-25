/*--------------------------------------------------------------*/
/* 								*/
/*		compute_ra_overstory.c		*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_ra_overstory - aerodynamic cond	*/
/*		for ovsretory layers in patch.			*/
/*								*/
/*	SYNOPSIS						*/
/*	double	compute_ra_overstory(		*/
/*						int,		*/
/*						double,		*/
/*						double,		*/
/*						double,		*/
/*						double,		*/
/*						*double);	*/
/*								*/
/*	returns:						*/
/*	ra - aerodynamic resistance of stratum			*/
/*								*/
/*	OPTIONS							*/
/*	int - verbose flag					*/
/*	u - wind speed measured in representative stratum at 	*/
/*			height z above ground surface. (m/s)	*/
/*	cn - wind attenuation coefficient			*/
/*	z - screen height					*/
/*	h - average height of canopy stratum (m)		*/
/*	z_u - average height of next lower layer (m)		*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	Computes aerodynamic resistance of what should be	*/
/*	the top canopy layer in a patch.  However, it may be	*/
/*	sufficient to consider any layer not close to the 	*/
/*	surface as having a log decay profile of windspeed 	*/
/*	which would qualify it for this algorithm.		*/
/*								*/
/*	No attempt has been made to add stability corrections	*/
/*		or to make use of knowledge of other strata or	*/
/*		or to take into account snowpack on/over strata	*/
/*								*/
/*	Reference: Heddeland, I and Lettenmaier, D. (1995)	*/
/*		   "Hydrological Modelling of Boreal Forest	*/
/*			Ecosystems" Water Resource Series #149,	*/
/*			Dept. of Civil Engineering, U of	*/
/*			Washington (pg. 38)			*/
/*								*/
/*	d and zo formula Taken from Xuewen Wang's version 	*/
/*	of rhessys C code.  					*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	A fatal error condition arises if the screen heigh	*/
/*	is lower than the zero plane height.  To avoid this the	*/
/*	calling routine should use data acquired at a higher	*/
/*	screen height or adjusted to a higher screen height.	*/ 
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"

double	compute_ra_overstory(
							 int	verbose_flag,
							 double  cn,
							 double	*u,
							 double	z,
							 double	h_o,
							 double  h_u,
							 double  *ga)
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration									*/
	/*--------------------------------------------------------------*/
	double compute_toc_wind(int,
		double,
		double,
		double);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	double	d_o, d_u;
	double	zo_o, zo_u;
	double	ra;
	double	ra_u;
	double 	u_o;

	/*--------------------------------------------------------------*/
	/* comput the zero plane displacement d (m)			*/
	/*								*/
	/*	Equation supplied by Xuewen Wang rhessys C code.	*/
	/*--------------------------------------------------------------*/
/*
	d_o =  pow(10.0, (0.979 * log10(h_o+0.001) - 0.154));
	d_u =  max(0.01*h_o, pow(10.0, (0.979 * log10(h_u+0.001) - 0.154)));
*/
    //09212023LML in case very low canopy, assume height is 0.01 m
    if (close_enough(h_o,0)) {
        h_o = 0.01;
        h_u = 0.5 * h_o;
    }

	d_o = 0.7 * h_o;
	d_u = 0.7 * h_u;

	/*--------------------------------------------------------------*/
	/*	Compute the roughness length zo (m)			*/
	/*								*/
	/*	Equation supplied by Xuewen Wang rhessys C code.	*/
	/*	Richard Fernandes:  This should change with snow pack .	*/
	/*--------------------------------------------------------------*/
/*
	zo_o = pow(10.0, (0.997 * log10(h_o+0.001) - 0.883));
	zo_u = max(0.01*h_o, pow(10.0, (0.997 * log10(h_u+0.001) - 0.883)));
*/
	zo_o = 0.1 * h_o;
	zo_u = 0.1 * h_u;


    //12232022LML with : R.H. Silversides (1978) Forest and Airport Wind Speeds, Atmosphere-Ocean,
    //16:3, 293-299, DOI: 10.1080/07055900.1978.9649036
    //calclulate wind speed above 2 meter above canopy height
    double z0a = 0.01; //(m) roughness of airport or bareground
    double za = z; //(m) screen hight of measurement
    double lamda = zo_o / (h_o - d_o);
    double z_f = 2.0;
    double uf_ua_ratio = pow(zo_o / z0a,0.07) * log(1. / lamda + z_f / zo_o)
                         / log(za / z0a);
    //printf("height:%lf uf_ua_ratio:%lf wind_10m:%lf ", h_o, uf_ua_ratio, *u);

    *u *= uf_ua_ratio;                                                          //(m/s) windspeed 2 meters above canopy height
    z = h_o + 2;    //the window speed already set at 2 meters above canopy height

    //printf(" wind_2m_above_canopy:%lf new_screen_hight:%lf\n", *u, z);

	/* Threshold code with ref height always above canopy height - see Heddeland p37*/
    //12232022LML commented out
    //if ( z < (h_o + 2)){
    //	z = h_o + 2;
    //}

	/*--------------------------------------------------------------*/
	/*	Compute the resistance to momentum transfer from a sourc*/
	/*	at the reference height but still in the the canopy	*/
	/*	stratum. ga (m/s)					*/
	/*--------------------------------------------------------------*/
	ra = pow( log( (z-d_o)/zo_o ) / 0.41, 2.0 ) / *u;
	ra_u = ra + log ( (z-d_o)/zo_o ) * h_o * exp(-1*cn)
		* ( exp(-1*cn*(d_u+zo_u)/h_o) - exp(-1*cn*(d_o+zo_o)/h_o))
		/ ( *u * 0.41 * 0.41 * cn * (h_o - d_o));

    //printf("*u:%lf\n",*u);


	/*--------------------------------------------------------------*/
	/*      Compute windspeed at top of logarithmic profile.        */
	/*      A log wind profile is used to compute toc wind          */
	/*      speed than an exponential profile is used to get        */
	/*      the within canopy attenuation.                          */
	/*      We assume the log profile starts at 0.1*h_o             */
	/*--------------------------------------------------------------*/
	u_o = compute_toc_wind( 0, *u, z, h_o) * exp(cn*(0.1-1));
	/*--------------------------------------------------------------*/
	/*	Compute exponential decay  through the canopy			*/
	/*--------------------------------------------------------------*/
	*u = max((u_o * exp(cn * (max(h_u, 0.1*h_o)/h_o - 1))), 0.0);
	/*--------------------------------------------------------------*/
	/*	if, this canopy below extends to 0.1*ho of the surface	*/
	/* 	include a logarithmic profile componenet of the near	*/
	/*	surface resistance					*/
	/*--------------------------------------------------------------*/
	if (h_u <= 0.1*h_o)
		ra_u += pow(log( (0.1*h_o)/ max(zo_u, 0.01 )),2.0) / (*u * 0.41*0.41);
	/*--------------------------------------------------------------*/
	/*	update conductance below this patch			*/
	/*--------------------------------------------------------------*/
	*ga = 1/ra_u;
	return(ra);
} /*compute_ra_overstory*/
