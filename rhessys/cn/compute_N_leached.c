/*--------------------------------------------------------------*/
/*                                                              */
/*		compute_N_leached				*/
/*                                                              */
/*  NAME                                                        */
/*		compute_N_leached				*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  void compute_N_leached(int					*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	,		*/
/*					double	);		*/
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*								*/
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"


double	compute_N_leached(int verbose_flag,
            double total_nitrate[],
			double Qout, 
			double s1, 
			double s2, 
            //double m,
            //double gamma,
			double n_0, 
			double p, 
            double N_decay_rate[],
			double z2_N, 
			double z2_water, 
            double N_absorption_rate[],
            double nleached[]
                          /*,
            double *transmissivity*/)
			
	{ 
	/*------------------------------------------------------*/ 
	/*	Local Function Declarations.						*/ 
	/*------------------------------------------------------*/
    	double  compute_delta_water(
                int,
                double,
                double,
                double,
                double,
                double);


	double  compute_z_final(
		int,
		double,
		double,
		double,
		double,
		double);
	double  compute_N_absorbed(
		int,
		double,
		double,
		double,
		double,
		double);
		
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
    //int didx_bot, didx_top;
    //double *nleached = (double*) malloc(LEACH_ELEMENT_counts * sizeof(double));
    double navail, nabsorbed;
    //double theta, sat_deficit;
    //double Q, Qtotal;
	double z1, z2;
	double	available_water,septic_depth;

    for (int i = 0; i < LEACH_ELEMENT_counts; i++) {
        nleached[i] = 0.0;
    }

	/*------------------------------------------------------*/
	/* nitrate export only occurs when Qout > 0.0		*/ 
	/*------------------------------------------------------*/
	if (Qout > ZERO) {
      if (s1 < 0.0) s1 = 0.0;
      if (s2 < s1) s2 = s1;
	/*------------------------------------------------------*/
	/*	first look at the case of return flow		*/
	/*	for return flow we must estimate the sat_deficit */
	/*	that would account for the flow			*/
	/*	(assuming all water leaves, so Qout/theta here */
	/*	is 1)						*/
	/*------------------------------------------------------*/
      if (close_enough(s1, 0.0) && close_enough(s2, 0.0)) {
		z2 = -1.0 * p * log (1 - (Qout) / (p * n_0));
		z1 = 0.0;
        for (int i = 0; i < LEACH_ELEMENT_counts; i++) {
          if (N_decay_rate[i] > ZERO) {
            navail = total_nitrate[i]
                / (1.0 - exp(-1.0 * N_decay_rate[i] * z2_N) )
                * (exp(-1.0 * N_decay_rate[i] * z1)
                - exp(-1.0 * N_decay_rate[i] * (z2)));
          } else {
            navail = total_nitrate[i] * (z2-z1)/z2_N;
          }
          nabsorbed = compute_N_absorbed(verbose_flag,
						z1,
						z2,
                        N_absorption_rate[i],
						p,
						n_0); 
	/*------------------------------------------------------*/
	/* in return flow Qout/theta = 1 so			*/
	/*------------------------------------------------------*/
          if (nabsorbed <= navail) {
            nleached[i] = navail - nabsorbed;
          }
        } //i
      } else {
	/*------------------------------------------------------*/
	/*	now for regular subsurface flow			*/
	/*	integrate through the saturated zone		*/
	/*------------------------------------------------------*/
        z2 = compute_z_final(
			verbose_flag,
			n_0,
			p,
			z2_water,		
			0.0,
			-s2);
        z1 = compute_z_final(
			verbose_flag,
			n_0,
			p,
			z2_water,		
			0.0,
			-s1);
        available_water = compute_delta_water(
            verbose_flag,
            n_0,p,z2_water,
            z2,
            z1);
        for (int i = 0; i < LEACH_ELEMENT_counts; i++) {
          if (N_decay_rate[i] > 0.0) {
            navail = total_nitrate[i]
                     / (1.0 - exp(-1.0 * N_decay_rate[i] * z2_N) )
                     * (exp(-1.0 * N_decay_rate[i] * z1)
                     - exp(-1.0 * N_decay_rate[i] * (z2)));

          }	else {
            septic_depth = -1.0*N_decay_rate[i];
			if (z1 > septic_depth)
                navail = 0.0;
			else
                navail = total_nitrate[i] * (z2-z1)/(z2_N - septic_depth);
          }
	
	/*------------------------------------------------------*/
	/* N-leached is mass flux of soluble nitrate	*/
	/* i.e n_avail / theta * outflow			*/
	/*------------------------------------------------------*/
          nabsorbed = compute_N_absorbed(verbose_flag,
			z1,
			z2,
            N_absorption_rate[i],
			p,
			n_0); 

          if (nabsorbed > navail) navail = 0;
          else navail -= nabsorbed;
		
          if (available_water > ZERO)
            nleached[i] = navail * min(1.,Qout/available_water);
	    }
      }
    }
    //return(nleached);
    return(0);
} /* end compute_N_leached */

