/*--------------------------------------------------------------*/
/* 								*/
/*		compute_potential_decomp					*/
/*								*/
/*								*/
/*	NAME							*/
/*	compute_potential_decomp -  					*/
/*		performs decomposition and updates soil/litter	*/
/*		carbon and nitrogen stores			*/
/*								*/
/*	SYNOPSIS						*/
/*	int compute_potential_decomp(					*/
/*			double,					*/
/*			double,					*/
/*			double,					*/
/*			double,					*/
/*			struct	soil_c_object	*		*/
/*			struct	soil_n_object	*		*/
/*			struct	litter_c_object	*		*/
/*			struct	litter_n_object	*		*/
/*			struct	cdayflux_patch_object *		*/
/*			struct	ndayflux_patch_object *		*/
/*				)				*/
/*								*/
/*	returns:						*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*		modified from Peter Thornton (1998)		*/
/*			dynamic - 1d-bgc ver4.0			*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

int compute_potential_decomp(double tsoil, double maxpsi,
							 double minpsi,
							 double theta,
							 double std,
							 struct  soil_c_object   *cs_soil,
							 struct  soil_n_object   *ns_soil,
							 struct  litter_c_object *cs_litr,
							 struct  litter_n_object *ns_litr,
							 struct cdayflux_patch_struct *cdf,
							 struct ndayflux_patch_struct *ndf,
							 struct patch_object *patch
							 )
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/

	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int ok;
	double rate_scalar, t_scalar, w_scalar;
	double a,b,c,d;
	double tk, thetai;
	double rfl1s1, rfl2s2,rfl4s3,rfs1s2,rfs2s3,rfs3s4;
	double kl1_base,kl2_base,kl4_base,ks1_base,ks2_base,ks3_base,ks4_base;
	double kl1,kl2,kl4,ks1,ks2,ks3,ks4;
	double cn_l1,cn_l2,cn_l4,cn_s1,cn_s2,cn_s3,cn_s4;
	double plitr1c_loss, plitr3c_loss, plitr2c_loss, plitr4c_loss;
	double psoil1c_loss, psoil2c_loss, psoil3c_loss, psoil4c_loss;
	double pmnf_l1s1,pmnf_l2s2,pmnf_l3l2, pmnf_l4s3,pmnf_s1s2,pmnf_s2s3,pmnf_s3s4,pmnf_s4;
	double potential_immob,mineralized;
	double weight1, weight2, theta1, theta2;
	double p_lignin, rate_landclim, rate_landclim_daily, psi, psi_max, psi_min, w_scalar_bgc;
	int nlimit, i;
	double et, w_scalar2, p_l1, p_l2, rate_scalar2, litr_decomp;

	p_lignin = 0.0, rate_landclim = 0.0, rate_landclim_daily = 0.0, psi = 0.0, psi_max = 0.0, psi_min = 0.0, w_scalar_bgc = 0.0;
	t_scalar = 0.0, w_scalar = 0.0; p_l1=0.0; p_l2=0.0, rate_scalar2 = 0.0;
	#define NUM_NORMAL  10 	/* resolution of normal distribution */
	double NORMAL[10]= {0,0,0.253,0.524,0.842,1.283,-0.253,-0.524,-0.842,-1.283};

    if ((cs_litr->litr1c + cs_litr->litr2c + cs_litr->litr3c + cs_litr->litr4c) > ZERO)
    {
    p_l1 = (cs_litr->litr1c / (cs_litr->litr1c + cs_litr->litr2c + cs_litr->litr3c + cs_litr->litr4c))*100;
    p_l2 = (cs_litr->litr2c / (cs_litr->litr1c + cs_litr->litr2c + cs_litr->litr3c + cs_litr->litr4c))*100;
    p_lignin = (cs_litr->litr4c / (cs_litr->litr1c + cs_litr->litr2c + cs_litr->litr3c + cs_litr->litr4c))*100;//0-100 this is yearly
    }
    else {
        p_l1 = 0.0;
        p_l2 = 0.0;
        p_lignin = 0.0;

    }

	ok = 0;
	/* calculate the rate constant scalar for soil temperature,
	assuming that the base rate constants are assigned for non-moisture
	limiting conditions at 25 C. The function used here is taken from
	Lloyd, J., and J.A. Taylor, 1994. On the temperature dependence of
	soil respiration. Functional Ecology, 8:315-323.
	This equation is a modification of their eqn. 11, changing the base
	temperature from 10 C to 25 C, since most of the microcosm studies
	used to get the base decomp rates were controlled at 25 C. */
	if (tsoil < -10.0){
		/* no decomp processes for tsoil < -10.0 C */
		t_scalar = 0.0;
	}
	else{
		tk = tsoil + 273.15;
		t_scalar = exp(308.56*((1.0/71.02)-(1.0/(tk-227.13))));
	}
	cs_litr->t_scalar = t_scalar;
	/* calculate the rate constant scalar for soil water content.
	use same empirical function as control on nitrification from NGas model
	but set parameters to reduce sensitivity to water stress
	(Parton et al, 1996 Global Biogeochemical cycles, 10:3, 401-412 ) */
    if (patch[0].soil_defaults[0][0].decom_model == 1) // 1 is default  2  is FireBGC 3 is landclim and
    {
	a=0.68; b=2.5; c=0.0012; d=2.84;
	w_scalar = 0.0;
	if (std > ZERO) {
		for (i=0; i<NUM_NORMAL; i++) {
			thetai = theta + NORMAL[i]*std;
			thetai = min(1.0, thetai);
			thetai = max(0.1, thetai);
			w_scalar  += (pow( ((thetai -b) / (a-b)), d*(b-a)/(a-c))
					* pow( ((thetai-c)/ (a-c)), d)) * 1.0/NUM_NORMAL;
			}
	}
	else {
		if ((theta <= ZERO) || (theta > 1.0))
			theta = 1.0;
		if (theta  > c) {
			w_scalar  = pow( ((theta -b) / (a-b)), d*(b-a)/(a-c))
			* pow( ((theta-c)/ (a-c)), d);
				}
		else
			w_scalar = 0.000001;
	}


	if (w_scalar > ZERO)
		w_scalar = sqrt(w_scalar);//here is the reason cause mismatch
	else
		w_scalar = ZERO;

	rate_scalar = w_scalar * t_scalar;
	/* assign output variables */
	//cs_litr->t_scalar = t_scalar;
	cs_litr->w_scalar = w_scalar;
	}//line 106 if

	else if(patch[0].soil_defaults[0][0].decom_model == 3) // 3 landclim
    {
    /*---------------below are different decomposition models, FireBGC and landClim*/

    et = patch[0].acc_year.et_decom_mean;
    double et_year = et*365*1000;

    //printf("\n [et (mm) %lf] \n", et_year);
    if(p_lignin > ZERO)
    {
        //rate_landclim = (-1.31369/365 + 0.05350* et_year + 0.18472 * (et_year/p_lignin))*0.01; //linear change -1.31369/365, this is year line 244, et is daily *365 to scale up
        rate_landclim = (-1.31369 + 0.05350* et_year + 0.18472 * (et_year/p_lignin))*0.01;
        //rate_landclim = sqrt(rate_landclim);
        //rate_landclim_daily = rate_landclim*100;
        rate_landclim_daily = 1.0 - pow((1.0 - rate_landclim), 0.002739726); //convert from yearly to daily plot_fig6_cc line 266 0.002739726 = 1/365 don't use 1/365
       // printf("\n [rate_landclim_daily %lf] \n", rate_landclim_daily);
    }
    else {rate_landclim_daily = 0.000001;}

    if (rate_landclim_daily < ZERO || isinf(rate_landclim_daily) || isnan(rate_landclim_daily))
        rate_landclim_daily = 0.000001; // give something to start with

    // printf("\n decomposition landclim [rate_landclim_daily %lf], [et_year(mm) %lf] [p_lignin %lf]\n", rate_landclim_daily, et_year, p_lignin);
        /*-----------------------must track the water scalar for compute_cwd_decay.c --*/
    a=0.68; b=2.5; c=0.0012; d=2.84;
	w_scalar2 = 0.0;
    if ((theta <= ZERO) || (theta > 1.0))
			theta = 1.0;
    if (theta  > c) {
			w_scalar2  = pow( ((theta -b) / (a-b)), d*(b-a)/(a-c))
			* pow( ((theta-c)/ (a-c)), d);
				}
		else
			w_scalar2 = 0.000001;

	if (w_scalar2 > ZERO)
		w_scalar2 = sqrt(w_scalar2);//here is the reason cause mismatch
	else
		w_scalar2 = ZERO;

	//rate_scalar = w_scalar * t_scalar;
	/* assign output variables */
	//cs_litr->t_scalar = t_scalar;
	cs_litr->w_scalar = w_scalar2;

    } // end if line 140


    else if (patch[0].soil_defaults[0][0].decom_model == 2)
    {
    // for FireBGCv2 model Keane et al 2011 and Biome-BGC (running & Coughlan 1988)
    // moisture/water scalar
    //psi = patch[0].psi;
    psi = patch[0].canopy_strata[0][0].epv.psi;
    //psi_max = -1.0*patch[0].soil_defaults[0][0].psi_max; //0.01
    //psi_min = -10;
    psi_max = -0.5; //psi_open
    psi_min = -2; //psi_close

    w_scalar_bgc = log(psi_min/psi)/log(psi_min/psi_max);
    if(w_scalar_bgc < ZERO || isinf(w_scalar_bgc) || isnan(w_scalar_bgc)){
        w_scalar_bgc = 0.000001;
    }

    //cs_litr->t_scalar = t_scalar; //these two scalar is used in other cwd_decay too so must track,but for landclim need to consider carefully
	cs_litr->w_scalar = w_scalar_bgc;

	rate_scalar = w_scalar_bgc * t_scalar; // use same temperature scalar as rhessys
    } //end if line 158
	/*---------------above different decomposition models-------*/
	//use rate_scalar
	// stratum[0].epv.psi
	//double psi;            /* (MPa) water potential of soil and leaves */
	//litr1c is labile; litr2c is cellulose; litr3c is shield cellulose; litr4c is lignan
	/*---------------------------------------------------------*/
	/* calculate compartment C:N ratios */
	if ((cs_litr->litr1c > ZERO) && (ns_litr->litr1n > ZERO ))	cn_l1 = cs_litr->litr1c/ns_litr->litr1n;
		else cn_l1 = LIVELAB_CN;
	if ((cs_litr->litr2c > ZERO) && (ns_litr->litr2n > ZERO))	cn_l2 = cs_litr->litr2c/ns_litr->litr2n;
		else cn_l2 = CEL_CN;
	if ((cs_litr->litr4c > ZERO) && (ns_litr->litr4n > ZERO))	cn_l4 = cs_litr->litr4c/ns_litr->litr4n;
		else cn_l4 = LIG_CN;
	cn_s1 = SOIL1_CN;
	cn_s2 = SOIL2_CN;
	cn_s3 = SOIL3_CN;
	cn_s4 = SOIL4_CN;
	/* respiration fractions for fluxes between compartments */
	rfl1s1 = 0.39;
	rfl2s2 = 0.55;
	rfl4s3 = 0.29;
	rfs1s2 = 0.28;
	rfs2s3 = 0.46;
	rfs3s4 = 0.55;

	/* calculate the corrected rate constants from the rate scalar and their
	base values. All rate constants are (1/day) */
    if (patch[0].soil_defaults[0][0].decom_model == 1 || patch[0].soil_defaults[0][0].decom_model == 2)
    {
	kl1_base = 0.7;      /* labile litter pool */
	kl2_base = 0.07;     /* cellulose litter pool */
	kl4_base = 0.014;     /* lignin litter pool */
	ks1_base = 0.07;    /* fast microbial recycling pool */
	ks2_base = 0.014;    /* medium microbial recycling pool */
	ks3_base = 0.0014;   /* slow microbial recycling pool */
	ks4_base = 0.0001;   /* recalcitrant SOM (humus) pool */

	kl1 = kl1_base * rate_scalar;
	kl2 = kl2_base * rate_scalar;
	kl4 = kl4_base * rate_scalar;
	ks1 = ks1_base * rate_scalar;
	ks2 = ks2_base * rate_scalar;
	ks3 = ks3_base * rate_scalar;
	ks4 = ks4_base * rate_scalar;
	//save these scalar
	cs_litr->rate_scalar = rate_scalar;//only for output
	cs_litr->litr_decomp = kl1*p_l1*0.01 + kl2*p_l2*0.01 + kl4*p_lignin*0.01; // only for ouput not used in model this is for same as R code, p unit is 0-100 not 0-1
    }
    else if (patch[0].soil_defaults[0][0].decom_model == 3)

    {

        kl1_base = 0.7;      /* labile litter pool */
        kl2_base = 0.07;     /* cellulose litter pool */
        kl4_base = 0.014;     /* lignin litter pool */
        ks1_base = 0.07;    /* fast microbial recycling pool */
        ks2_base = 0.014;    /* medium microbial recycling pool */
        ks3_base = 0.0014;   /* slow microbial recycling pool */
        ks4_base = 0.0001;   /* recalcitrant SOM (humus) pool */
        rate_scalar = cs_litr->w_scalar * cs_litr->t_scalar; // this is rhessyd default
        // how to assign the daily decomposition rate of landclim to different pool
        // another solution is scale up based rate_landclim_daily/(kl1_base*rate_scalar + kl2_base*rate_scalar + kl4_base*rate_scalar)
      kl1 = rate_scalar * kl1_base;
      kl2 = rate_scalar * kl2_base;
      kl4 = rate_scalar * kl4_base;
      //printf("\n [kl1 %lf] \n", kl1);
      litr_decomp = kl1*p_l1*0.01 + kl2*p_l2*0.01 + kl4*p_lignin*0.01;
      cs_litr->litr_decomp = litr_decomp;
      if(litr_decomp > ZERO)
      {
        rate_scalar2 = rate_landclim_daily / litr_decomp;
      }
      else {
        rate_scalar2 = 1;
      }
     // for soil decomp use the original model

      kl1 = rate_scalar * kl1_base*rate_scalar2;
      kl2 = rate_scalar * kl2_base*rate_scalar2;
      kl4 = rate_scalar * kl4_base*rate_scalar2;

     // soil use default
      ks1 = rate_scalar * ks1_base;// *rate_scalar2;
      ks2 = rate_scalar * ks2_base;// *rate_scalar2;
      ks3 = rate_scalar * ks3_base;// *rate_scalar2;
      ks4 = rate_scalar * ks4_base;// *rate_scalar2;

      cs_litr->rate_landclim_year = rate_landclim;
      cs_litr->rate_landclim_daily = rate_landclim_daily;
      //cs_litr->litr_decomp = kl1*p_l1*0.01 + kl2*p_l2*0.01 + kl4*p_lignin*0.01; // this is for output the same same as R code, p unit is 0-100 not 0-1
    //cs_litr->litr_decomp = kl1+ kl2 + kl4;
    }

	/* initialize the potential loss and mineral N flux variables */
	plitr1c_loss = plitr2c_loss = plitr3c_loss = plitr4c_loss = 0.0;
	psoil1c_loss = psoil2c_loss = psoil3c_loss = psoil4c_loss = 0.0;
	pmnf_l1s1 = pmnf_l2s2 = pmnf_l3l2 = pmnf_l4s3 = 0.0;
	pmnf_s1s2 = pmnf_s2s3 = pmnf_s3s4 = pmnf_s4 = 0.0;

	/* calculate the non-nitrogen limited fluxes between litter and
	soil compartments. These will be ammended for N limitation if it turns
	out the potential gross immobilization is greater than potential gross
	mineralization. */
	/* 1. labile litter to fast microbial recycling pool */
	if ((cs_litr->litr1c > ZERO) && (ns_litr->litr1n > ZERO)) {
		plitr1c_loss = kl1 * cs_litr->litr1c;
		pmnf_l1s1 = (plitr1c_loss * (1.0 - rfl1s1 - (cn_s1/cn_l1)))/cn_s1;
	}
	/* 2. cellulose litter to medium microbial recycling pool */
	if ((ns_litr->litr2n > ZERO) && (cs_litr->litr2c > ZERO)) {
		plitr2c_loss = kl2 * cs_litr->litr2c;
		pmnf_l2s2 = (plitr2c_loss * (1.0 - rfl2s2 - (cn_s2/cn_l2)))/cn_s2;
	}
	/* 2b. shield cellulose litter to goes to cellulose litter pool */
	/* respiration fractions not available to assume the same as for lignan (the "shield") */
	if ((ns_litr->litr3n > ZERO) && (cs_litr->litr3c > ZERO)) {
		plitr3c_loss = kl4 * cs_litr->litr3c;
		pmnf_l3l2 = (plitr3c_loss * (1.0 - rfl4s3 - (cn_l2/cn_l2)))/cn_l2;
	}
	/* 3. lignin litter to slow microbial recycling pool */
	if ((ns_litr->litr4n > ZERO) && (cs_litr->litr4c > ZERO)) {
		plitr4c_loss = kl4 * cs_litr->litr4c;
		pmnf_l4s3 = (plitr4c_loss * (1.0 - rfl4s3 - (cn_s3/cn_l4)))/cn_s3;
	}
	/* 4. fast microbial recycling pool to medium microbial recycling pool */
	if ((ns_soil->soil1n > ZERO) && (cs_soil->soil1c > ZERO)) {
		psoil1c_loss = ks1 * cs_soil->soil1c;
		pmnf_s1s2 = (psoil1c_loss * (1.0 - rfs1s2 - (cn_s2/cn_s1)))/cn_s2;
	}
	/* 5. medium microbial recycling pool to slow microbial recycling pool */
	if ((ns_soil->soil2n > ZERO) && (cs_soil->soil2c > ZERO)) {
		psoil2c_loss = ks2 * cs_soil->soil2c;
		pmnf_s2s3 = (psoil2c_loss * (1.0 - rfs2s3 - (cn_s3/cn_s2)))/cn_s3;
	}
	/* 6. slow microbial recycling pool to recalcitrant SOM pool */
	if ((ns_soil->soil3n > ZERO) && (cs_soil->soil3c > ZERO)) {
		psoil3c_loss = ks3 * cs_soil->soil3c;
		pmnf_s3s4 = (psoil3c_loss * (1.0 - rfs3s4 - (cn_s4/cn_s3)))/cn_s4;
	}
	/* 7. mineralization of recalcitrant SOM */
	if ((ns_soil->soil4n > ZERO) && (cs_soil->soil4c > ZERO)) {
		psoil4c_loss = ks4 * cs_soil->soil4c;
		pmnf_s4 = -psoil4c_loss/cn_s4;
	}

	/* determine if there is sufficient mineral N to support potential
	immobilization. Immobilization fluxes are positive, mineralization fluxes
	are negative */
	nlimit = 0;
	potential_immob = 0.0;
	mineralized = 0.0;
	if (pmnf_l1s1 > 0.0) potential_immob += pmnf_l1s1;
	else mineralized += -pmnf_l1s1;

	if (pmnf_l2s2 > 0.0) potential_immob += pmnf_l2s2;
	else mineralized += -pmnf_l2s2;

	if (pmnf_l3l2 > 0.0) potential_immob += pmnf_l3l2;
	else mineralized += -pmnf_l3l2;

	if (pmnf_l4s3 > 0.0) potential_immob += pmnf_l4s3;
	else mineralized += -pmnf_l4s3;

	if (pmnf_s1s2 > 0.0) potential_immob += pmnf_s1s2;
	else mineralized += -pmnf_s1s2;

	if (pmnf_s2s3 > 0.0) potential_immob += pmnf_s2s3;
	else mineralized += -pmnf_s2s3;

	if (pmnf_s3s4 > 0.0) potential_immob += pmnf_s3s4;
	else mineralized += -pmnf_s3s4;
	mineralized += -pmnf_s4;

	/* save the potential fluxes until plant demand has been assessed,
	to allow competition between immobilization fluxes and plant growth
	demands */
	ndf->mineralized = mineralized;
	ndf->potential_immob = potential_immob;
	cdf->plitr1c_loss = plitr1c_loss;
	ndf->pmnf_l1s1 = pmnf_l1s1;
	cdf->plitr2c_loss = plitr2c_loss;
	ndf->pmnf_l2s2 = pmnf_l2s2;
	cdf->plitr3c_loss = plitr3c_loss;
	ndf->pmnf_l3l2 = pmnf_l3l2;
	cdf->plitr4c_loss = plitr4c_loss;
	ndf->pmnf_l4s3 = pmnf_l4s3;
	cdf->psoil1c_loss = psoil1c_loss;
	ndf->pmnf_s1s2 = pmnf_s1s2;
	cdf->psoil2c_loss = psoil2c_loss;
	ndf->pmnf_s2s3 = pmnf_s2s3;
	cdf->psoil3c_loss = psoil3c_loss;
	ndf->pmnf_s3s4 = pmnf_s3s4;
	cdf->psoil4c_loss = psoil4c_loss;
	ndf->pmnf_s4 = pmnf_s4;
	cdf->kl4 = kl4;

	return(ok);
} /* end compute_potential_decomp.c */

