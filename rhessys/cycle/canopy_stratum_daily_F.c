/*--------------------------------------------------------------*/
/* 								*/
/*			canopy_stratum_daily_F			*/
/*								*/
/*	NAME							*/
/*	canopy_statrum_daily_F 					*/
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	void canopy_stratum_daily_F 				*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "rhessys.h"
#include "phys_constants.h"
#include "functions.h"



void	canopy_stratum_daily_F(
							   struct	world_object		      *world,
							   struct	basin_object		      *basin,
							   struct	hillslope_object	    *hillslope,
							   struct	zone_object		        *zone,
							   struct	patch_object		      *patch,
							   struct layer_object		      *layer,
							   struct canopy_strata_object 	*stratum,
                          struct canopy_strata_object  *shadow_strata,
							   struct command_line_object	  *command_line,
							   struct	tec_entry		          *event,
							   struct date 			            current_date)
{
	/*--------------------------------------------------------------*/
	/*	Local function declaration				*/
	/*--------------------------------------------------------------*/
	long	julday(struct date calendar_date);

	double	compute_diffuse_radiative_fluxes(
		int,
		double *,
		double *,
		double,
		double,
		double,
		double,
		double,
		double);

	double	compute_diffuse_radiative_PAR_fluxes(
		int	,
		double	*,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	);

	double	compute_direct_radiative_fluxes(
		int	,
		double	*,
		double  *,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	);


	int	compute_maint_resp(	double,
		double,
		struct  cstate_struct *,
		struct  nstate_struct *,
		struct epconst_struct *,
		struct  metvar_struct *,
		struct cdayflux_struct *);

	double	 compute_snow_stored(
		int	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	*,
		double	*,
		struct	canopy_strata_object *);

	double	 compute_rain_stored(
		int	,
		double	*,
		struct	canopy_strata_object *);

	double  compute_ra_overstory(
		int     ,
		double  ,
		double *,
		double  ,
		double  ,
		double	,
		double *);

	double  compute_ra_understory(
		int     ,
		double  ,
		double *,
		double  ,
		double  ,
		double  *);

	double  compute_ra_surface(
		int     ,
		double  ,
		double *,
		double  ,
		double  ,
		double  *);

	double	compute_vascular_stratum_conductance(
		int	,
		int	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double, int,
		struct canopy_strata_object *,
		struct patch_object *);

	double	compute_nonvascular_stratum_conductance(
		int	,
		double	,
		double	,
		double	,
		double	,
		double	);

	double	compute_surface_heat_flux(
		int	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	);

	double	penman_monteith(
		int	,
		double	,
		double	,
		double	,
		double	,
		double	,
		double	,
		int	);

	int	compute_farq_psn(
		struct psnin_struct * ,
		struct psnout_struct * ,
		int);

	double	compute_potential_N_uptake_Dickenson(
		struct	epconst_struct,
		struct	epvar_struct *,
		struct cstate_struct *cs,
		struct nstate_struct *ns,
		struct cdayflux_struct *);

	double	compute_potential_N_uptake_Waring(
		struct	epconst_struct,
		struct	epvar_struct *,
		struct cstate_struct *cs,
		struct nstate_struct *ns,
		struct cdayflux_struct *);

	double	compute_potential_N_uptake_combined(
		struct	epconst_struct,
		struct	epvar_struct *,
		struct cstate_struct *cs,
		struct nstate_struct *ns,
		struct cdayflux_struct *);

	double	compute_potential_N_uptake(
		struct	epconst_struct,
		struct	epvar_struct *,
		struct cstate_struct *cs,
		struct nstate_struct *ns,
		struct cdayflux_struct *);

	void	update_mortality(
		struct epconst_struct,
		struct cstate_struct *,
		struct cdayflux_struct *,
		struct cdayflux_patch_struct *,
		struct nstate_struct *,
		struct ndayflux_struct *,
		struct ndayflux_patch_struct *,
		struct litter_c_object *,
		struct litter_n_object *,
		int,
		struct mortality_struct);

	struct 	canopy_strata_object *construct_empty_shadow_strata(
		struct command_line_object *,
		FILE	*,
		struct	patch_object *,
		int     num_world_base_stations,
		struct  base_station_object **world_base_stations,
		struct	default_object	*defaults);


  void	update_shadow_strata(
	  struct	world_object		      *world,
	  struct    patch_object              *patch,
	  struct canopy_strata_object 	*stratum,
      struct canopy_strata_object   *shadow_strata,
	  struct command_line_object	  *command_line,
	  struct date 			            current_date);


   void compute_snag_decay(
					  struct cstate_struct *cs,
					  struct nstate_struct *ns,
					  struct cdayflux_patch_struct *cdf, //daily carbon flux
					  struct ndayflux_patch_struct *ndf,  // daily nitrogen flux
					  int year_delay,
					  int half_life,
					  struct date current_day); // current day that

   void compute_redneedle_decay(
                      struct epconst_struct epc,
					  //double cover_fraction,
					  struct cstate_struct *cs,
					  struct nstate_struct *ns,
					  struct litter_c_object *cs_litr,
					  struct litter_n_object *ns_litr,
					  struct cdayflux_patch_struct *cdf, //daily carbon flux
					  struct ndayflux_patch_struct *ndf,  // daily nitrogen flux
					  int year_delay,
					  int half_life,
					  struct date current_day);

    void update_deadroot_decay(
					  struct epconst_struct epc,
					  struct cstate_struct *cs,
					  struct nstate_struct *ns,
					  struct litter_c_object *cs_litr,
					  struct litter_n_object *ns_litr,
					  struct cdayflux_patch_struct *cdf, //daily carbon flux
					  struct ndayflux_patch_struct *ndf,  // daily nitrogen flux
					  int half_life);


	/*--------------------------------------------------------------*/
	/*  Local variable definition.                                  */
	/*--------------------------------------------------------------*/
	double tmid;
	double  assim_sunlit;
	double  assim_shade;
	double	dC13_sunlit, dC13_shade;
	double  APAR_direct_sunlit;
	double	dry_evaporation;
	double	ga, gasnow;
	double	Kdown_direct;
	double  Kup_direct;
	double	Kdown_diffuse;
	double  Kup_diffuse;
	double	PAR_diffuse;
	double	PAR_direct;
	double  total_incoming_PAR;
	double	perc_sunlit;
	double	potential_evaporation_rate = 0.0;
	double	potential_evaporation_rate_night = 0.0;
	double	potential_evaporation_rate_day = 0.0;
	double	potential_rainy_evaporation_rate;
	double	potential_rainy_evaporation_rate_night;
	double	potential_rainy_evaporation_rate_day;
	double	rainy_evaporation;
	double	rain_throughfall;
	double	NO3_throughfall;
	double	NO3_stored;
	double	rnet_evap;
	double  rnet_evap_night;
	double  rnet_evap_day;
	double	rnet_trans, rnet_trans_sunlit, rnet_trans_shade;
	double	snow_throughfall;
	double	transpiration;
	double	transpiration_rate;
	double	transpiration_rate_sunlit;
	double	transpiration_rate_shade;
	double	potential_transpiration;
	double	potential_transpiration_rate;
	double	potential_transpiration_rate_sunlit;
	double	potential_transpiration_rate_shade;
	double  wind, windcan, windsnow, ustar;
	double leafcloss_perc;

	double m_APAR_sunlit;
	double m_tavg_sunlit;
	double m_LWP_sunlit;
	double m_CO2_sunlit;
	double m_tmin_sunlit;
	double m_vpd_sunlit;
	double m_APAR_shade;
	double m_tavg_shade;
	double m_LWP_shade;
	double m_CO2_shade;
	double m_tmin_shade;
	double m_vpd_shade;

	double	fraction_direct_K_used;
	double	fraction_diffuse_K_used;
	double	fraction_direct_APAR_used;
	double	fraction_diffuse_APAR_used;
	double  fraction_L_used;
	double	fraction_L_used_night = 0.0;
	double  fraction_L_used_day = 0.0;
	double  fraction_surfheat_used;
	double  lhvap;
	double	Rnet_used, Rnet_canopy;
	double	APAR_used;
	double K_reflectance, PAR_reflectance;

	double dum;
	double max_snow_albedo_increase, wetfrac;
	double deltaT;

     int ok=0; //for beetle outbreak
     double time_diff;
     double leaf_time_diff;
     double year_delay, leaf_year_delay;
     struct date impact_date;
     struct  dated_sequence	clim_event;
     struct  dated_sequence2 clim_event1;
     struct  dated_sequence2 clim_event2;
     int inx;
     double overstory_burn, redneedlec_burn, delay_redneedlec_burn, snagc_burn, delay_snagc_burn, redneedlen_burn, delay_redneedlen_burn, snagn_burn, delay_snagn_burn;
//     int num_records;

	struct	psnin_struct	psnin;
	struct	psnout_struct	psnout;
	struct mortality_struct mort;
	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8d -444.1 ",julday(current_date)-2449000);
	if ( command_line[0].verbose_flag > 1 )
		printf("%8.2f %8.2f %8.2f %8.2f ",
		patch[0].Kdown_direct,
		patch[0].Kdown_diffuse,
		patch[0].PAR_direct/1000,
		patch[0].PAR_diffuse/1000);

  /*--------------------------------------------------------------*/
	/*	Initialize stratum variables.				*/
	/*--------------------------------------------------------------*/
	stratum[0].Kstar_diffuse = 0.0;
	stratum[0].APAR_diffuse = 0.0;
	stratum[0].Kstar_direct = 0.0;
	stratum[0].APAR_direct = 0.0;
	stratum[0].potential_evaporation = 0.0;
	 m_APAR_sunlit=0.0;
	 m_tavg_sunlit=0.0;
	 m_LWP_sunlit=0.0;
	 m_CO2_sunlit=0.0;
	 m_tmin_sunlit=0.0;
	 m_vpd_sunlit=0.0;
	 m_APAR_shade=0.0;
	 m_tavg_shade=0.0;
	 m_LWP_shade=0.0;
	 m_CO2_shade=0.0;
	 m_tmin_shade=0.0;
	 m_vpd_shade=0.0;
	stratum[0].mult_conductance.vpd = 0.0;
	stratum[0].mult_conductance.tmin = 0.0;
	stratum[0].mult_conductance.tavg = 0.0;
	stratum[0].mult_conductance.CO2 = 0.0;
	stratum[0].mult_conductance.LWP = 0.0;
	stratum[0].mult_conductance.APAR = 0.0;
	dum = 0.0;
	wetfrac = 0.0;
	rnet_evap = 0.0;
	rnet_evap_night = 0.0;
	rnet_evap_day = 0.0;

	stratum[0].canopy_drip = 0.0;
    stratum[0].assim_sunlit = 0.0;
    stratum[0].assim_shade = 0.0;
    stratum[0].trans_sunlit = 0.0;
    stratum[0].trans_shade = 0.0;

	K_reflectance = 0.0;
	PAR_reflectance = 0.0;
	fraction_direct_K_used = 0.0;
	fraction_diffuse_K_used = 0.0;
	fraction_direct_APAR_used = 0.0;
	fraction_diffuse_APAR_used = 0.0;
	fraction_L_used = 0.0;
	fraction_surfheat_used = 0.0;
	Rnet_used = 0.0;
	Rnet_canopy = 0.0;
	APAR_used = 0.0;

	max_snow_albedo_increase = 0.2;

	stratum[0].evaporation = 0.0;

	deltaT = 0.0;

	/*--------------------------------------------------------------*/
	/*	Initialize temporary variables for transmitted fluxes.	*/
	/*--------------------------------------------------------------*/
	Kdown_diffuse = patch[0].Kdown_diffuse ;
	Kup_diffuse = patch[0].Kup_diffuse;
	PAR_diffuse = patch[0].PAR_diffuse;
	Kdown_direct = patch[0].Kdown_direct;
	Kup_direct = patch[0].Kup_direct;
	PAR_direct = patch[0].PAR_direct;
	rain_throughfall = patch[0].rain_throughfall;
	snow_throughfall = patch[0].snow_throughfall;
	ga = patch[0].ga;
	gasnow = patch[0].gasnow;
	wind  = patch[0].wind;
	windcan  = patch[0].wind;
	windsnow  = patch[0].windsnow;
	ustar  = patch[0].ustar;
	transpiration = 0.0;
	potential_transpiration = 0.0;
	rainy_evaporation = 0;
	dry_evaporation = 0;
	total_incoming_PAR = PAR_diffuse + PAR_direct;
	NO3_stored=0;
    NO3_throughfall=patch[0].NO3_throughfall;                                   //06072022LML

	ustar = patch[0].ustar;

	lhvap = (2.5023e6 - 2430.54 * zone[0].metv.tday)/1000.0; /* KJ/kg H2O */

	double daylength = zone[0].metv.dayl;
	double nightlength = SECONDS_PER_DAY - daylength;
	double day_proportion = daylength / SECONDS_PER_DAY;
	double night_proportion = 1 - day_proportion;

	double rain_duration_day = zone[0].rain_duration * day_proportion;
	double rain_duration_night = zone[0].rain_duration - rain_duration_day;

	/* Lstar calcs are done in patch daily F AFTER this routine, so using yesterday's	*/
	/* patch total canopy Lstar and scaling back to this stratum by cover fraction.		*/
	/*stratum[0].Lstar = patch[0].Lstar_canopy / stratum[0].cover_fraction;*/

	/*--------------------------------------------------------------*/
	/*	perform plant grazing losses				*/
	/*	(if grow flag is on)					*/
	/*--------------------------------------------------------------*/
	if ((stratum[0].cs.leafc > ZERO) && (patch[0].grazing_Closs > ZERO) && command_line[0].grow_flag == 1 ) {
			leafcloss_perc  = patch[0].grazing_Closs * (stratum[0].ns.leafn/stratum[0].cs.leafc)
						 / patch[0].grazing_mean_nc / stratum[0].cs.leafc;
			mort.mort_leafc  = leafcloss_perc;
			mort.mort_cpool = 0.0;
			mort.mort_deadleafc = 0.0;
			mort.mort_livestemc = 0.0;
			mort.mort_deadstemc = 0.0;
			mort.mort_livecrootc = 0.0;
			mort.mort_deadcrootc = 0.0;
			mort.mort_frootc = 0.0;
			update_mortality(stratum[0].defaults[0][0].epc,
				&(stratum[0].cs),
				&(stratum[0].cdf),
				&(patch[0].cdf),
				&(stratum[0].ns),
				&(stratum[0].ndf),
				&(patch[0].ndf),
				&(patch[0].litter_cs),
				&(patch[0].litter_ns),
				2,
				mort);
			printf("\n completed %lf from %lf", leafcloss_perc, stratum[0].cs.leafc);

	}

	/*	NEW CHECK TO SEE IF STRATUM IS VEGETATED, OTHERWISE SKIP	*/
	/*	TO END TO UPDATE RADIATION BY COVER FRACTION				*/

	if ((stratum[0].defaults[0][0].epc.veg_type != NON_VEG) && (stratum[0].epv.proj_pai > ZERO) && (stratum[0].cover_fraction > ZERO)) {

		if (stratum[0].snow_stored > ZERO) {
			/* REFLECT TOO HIGH WITH SNOW IN CANOPY SO REDUCING BY WET FRACTION */
			/* Wet fraction from Deardorff 1978) */
			//wetfrac = pow(stratum[0].snow_stored/(stratum[0].epv.all_pai * stratum[0].defaults[0][0].specific_snow_capacity), 0.6667); //20181101
			wetfrac = pow(stratum[0].snow_stored/(stratum[0].epv.all_pai_when_red * stratum[0].defaults[0][0].specific_snow_capacity), 0.6667); //NREN chang all_pai_all to all_pai_when_red
			wetfrac = min(wetfrac, 1.0);
			/* From Kuusinen 2012 boreal pine forest snow albedo measurments. Max albedo increase due to snow */
			/* interception ~ 0.3 and no change with interception < 50%. Assuming linear scaling between  */
			/* 50% and 100% wet fraction. */
			if (wetfrac < 0.5) {
				K_reflectance = stratum[0].defaults[0][0].K_reflectance;
				PAR_reflectance = stratum[0].defaults[0][0].PAR_reflectance;
				}
			else {
				K_reflectance = stratum[0].defaults[0][0].K_reflectance + 2.0 * (max_snow_albedo_increase * (wetfrac - 0.5));
				PAR_reflectance = stratum[0].defaults[0][0].PAR_reflectance + 2.0 * (max_snow_albedo_increase * (wetfrac - 0.5));
				}
			/*K_reflectance = snow_albedo * wetfrac + stratum[0].defaults[0][0].K_reflectance * (1.0 - wetfrac);
			PAR_reflectance = snow_albedo * wetfrac + stratum[0].defaults[0][0].PAR_reflectance * (1.0 - wetfrac);*/
		}
		else {
			K_reflectance = stratum[0].defaults[0][0].K_reflectance;
			PAR_reflectance = stratum[0].defaults[0][0].PAR_reflectance;
		}

		if ( command_line[0].verbose_flag == -5 ){
			printf("\n     STRATUM DAILY ALBEDO: K_refl=%lf PAR_refl=%lf snow_stor=%lf wetfrac=%lf snowcap=%lf",
				   K_reflectance,
				   PAR_reflectance,
				   stratum[0].snow_stored,
				   wetfrac,
				   (stratum[0].epv.all_pai_when_red * stratum[0].defaults[0][0].specific_snow_capacity));
		}


	/*--------------------------------------------------------------*/
	/*  Intercept diffuse radiation.                                */
	/*  We assume that the zone slope == patch slope.               */
	/*  We also assume that radiation reflected into the upper      */
	/*      hemisphere is lost.                                     */
	/*  We do not make adjustements for chaanging gap fraction over */
	/*      the diurnal cycle - using a suitable mean gap fraction  */
	/*      instead.                                                */
	/*  We do take into account the patch level horizon which will  */
	/*      allow simulation of clear-cuts/clearings with low sza's */
	/*  Note that for zone level temperature and radiation          */
	/*      computation we took into account only zone level horizon*/
	/*      since we figured that climate above the zone was well   */
	/*      mixed.                                                  */
	/*--------------------------------------------------------------*/
	if ( command_line[0].verbose_flag > 2 ){
		printf("\n%4d %4d %4d -444.4 \n",
			current_date.day, current_date.month, current_date.year);
		printf("\n %f %f %f %f %f %f %f",
			Kdown_diffuse,
			Kdown_direct,
			-1*stratum[0].defaults[0][0].epc.ext_coef,
			stratum[0].gap_fraction,
			stratum[0].epv.proj_pai_when_red,
			basin[0].theta_noon,
			K_reflectance);
		printf("\n%4d %4d %4d -444.4b \n",
			current_date.day, current_date.month, current_date.year);
	}
	stratum[0].Kstar_diffuse = compute_diffuse_radiative_fluxes(
		command_line[0].verbose_flag,
		&(Kdown_diffuse),
		&(Kup_diffuse),
		Kdown_direct,
		-1*stratum[0].defaults[0][0].epc.ext_coef,
		stratum[0].gap_fraction,
		//stratum[0].epv.proj_pai,
		stratum[0].epv.proj_pai_when_red, //NREN 20180804 consider the red phase after beetle attack
		basin[0].theta_noon,
		K_reflectance);
	if ( command_line[0].verbose_flag > 2 )
		printf("\n%8d -444.5 ",julday(current_date)-2449000);
	stratum[0].APAR_diffuse = compute_diffuse_radiative_PAR_fluxes(
		command_line[0].verbose_flag,
		&(PAR_diffuse),
		PAR_direct,
		-1*stratum[0].defaults[0][0].epc.ext_coef,
		stratum[0].gap_fraction,
		stratum[0].epv.proj_lai, // DONOT change here lai 20181101
		//stratum[0].epv.proj_lai_when_red, //psn using the real LAI
		basin[0].theta_noon,
		PAR_reflectance);
	/*--------------------------------------------------------------*/
	/*  Intercept direct radiation.                                 */
	/*      hard to measure for each strata.  We could use top      */
	/*      of canopy albedo but this integrates the effect of the  */
	/*      entire canopy.  Furthermore, it requires in general     */
	/*      knowledge of the canopy BRDF - which we want to avoid.  */
	/*      Instead we assume a certain reflectance and             */
	/*      transmittance for each strata's canopy elements.        */
	/*  We assume that the zone slope == patch slope.               */
	/*  We also assume that radiation reflected into the upper      */
	/*      hemisphere is lost.                                     */
	/*  We do not make adjustements for chaanging gap fraction over */
	/*      the diurnal cycle - using a suitable mean gap fraction  */
	/*      instead.                                                */
	/*  We do take into account the patch level horizon which will  */
	/*      allow simulation of clear-cuts/clearings with low sza's */
	/*  Note that for zone level temperature and radiation          */
	/*      computation we took into account only zone level horizon*/
	/*      since we figured that climate above the zone was well   */
	/*      mixed.                                                  */
	/*--------------------------------------------------------------*/
	if ( command_line[0].verbose_flag > 2 )
		printf("\n%8d -444.2 ",julday(current_date)-2449000);
	stratum[0].Kstar_direct = compute_direct_radiative_fluxes(
		command_line[0].verbose_flag,
		&(Kdown_direct),
		&(Kup_direct),
		-1*stratum[0].defaults[0][0].epc.ext_coef,
		stratum[0].gap_fraction,
		//stratum[0].epv.proj_pai,
		stratum[0].epv.proj_pai_when_red, //NREN 20180804
		basin[0].theta_noon,
		K_reflectance,
		K_reflectance);

	if ( command_line[0].verbose_flag >2  )
		printf("\n%d %d %d  -444.3 ",
		current_date.year, current_date.month, current_date.day);
	stratum[0].APAR_direct = compute_direct_radiative_fluxes(
		command_line[0].verbose_flag,
		&(PAR_direct),
		&(dum),
		-1*stratum[0].defaults[0][0].epc.ext_coef,
		stratum[0].gap_fraction,
		stratum[0].epv.proj_pai, //DNOT Change here 20181101 this is for psn
		//stratum[0].epv.proj_pai_when_red, //NREN 20180804 20181030 Ning Ren
		basin[0].theta_noon,
		PAR_reflectance,
		PAR_reflectance);


	/*--------------------------------------------------------------*/
	/*  Calculate net longwave.										*/
	/*--------------------------------------------------------------*/
	if (command_line[0].evap_use_longwave_flag) {
		compute_Lstar_canopy(command_line[0].verbose_flag, //not related to LAI or pai NREN 20180804
				stratum[0].Kstar_direct + stratum[0].Kstar_diffuse,
				stratum[0].snow_stored,
				zone,
				patch,
				stratum);
	}


	if ( stratum[0].Kstar_direct < -1 ) {
			printf("CANOPY_START ID=%d: pai=%lf snowstor=%lf APARused=%lf APARdir=%lf APAR=%lf Rnet_used=%lf Kstardir=%lf Kstar=%lf Lstar=%lf \n",
				   stratum[0].ID,
				   stratum[0].epv.proj_pai_when_red, //NREN
				   stratum[0].snow_stored,
				   APAR_used,
				   stratum[0].APAR_direct,
				   (stratum[0].APAR_direct+stratum[0].APAR_diffuse)/86.4,
				   Rnet_used/86.4,
				   stratum[0].Kstar_direct/86.4,
				   (stratum[0].Kstar_direct+stratum[0].Kstar_diffuse)/86.4,
				   stratum[0].Lstar/86.4);
		}


	/*--------------------------------------------------------------*/
	/*	record fraction of PAR absorption 			*/
	/*	and seasonal low water content				*/
	/*	(for use in dynamic version to limit allocation)	*/
	/*--------------------------------------------------------------*/
	stratum[0].epv.min_vwc = min((patch[0].unsat_storage / patch[0].sat_deficit),
		stratum[0].epv.min_vwc);
	if ( command_line[0].verbose_flag >2  )
		printf("\n fparabs %f min_vwc %f",stratum[0].epv.max_fparabs,
		stratum[0].epv.min_vwc);
	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8d -444.6 ",julday(current_date)-2449000);

	/*--------------------------------------------------------------*/
	/*	Compute conductance aerodynamic.			*/
	/*--------------------------------------------------------------*/
	if ( patch[0].snowpack.height <= stratum[0].epv.height ){
		/*--------------------------------------------------------------*/
		/*		Compute by brute force if usrat_overu given.	*/
		/*--------------------------------------------------------------*/
        if ( close_enough(stratum[0].defaults[0][0].ustar_overu, -999.9 )){
			if ( command_line[0].verbose_flag == -5 ){
				printf("\n     STRATUM DAILY F: PAI0=%lf tmin=%lf tmax=%lf tavg=%lf RH=%lf edewpt=%lf vpd=%lf Kdir=%lf Kdif=%lf Lstar=%lf Lstarpch=%lf\n          K_refl=%lf PAR_refl=%lf snow_stor=%lf wetfrac=%lf snowcap=%lf",
				    stratum[0].epv.proj_pai_when_red,
					zone[0].metv.tmin,
					zone[0].metv.tmax,
					zone[0].metv.tavg,
					zone[0].relative_humidity,
					zone[0].e_dewpoint,
					zone[0].metv.vpd,
					stratum[0].Kstar_direct/86.4,
					stratum[0].Kstar_diffuse/86.4,
					stratum[0].Lstar/86.4,
					patch[0].Lstar_canopy/86.4,
					K_reflectance,
					PAR_reflectance,
					   stratum[0].snow_stored,
					   wetfrac,
					   (stratum[0].epv.all_pai_when_red * stratum[0].defaults[0][0].specific_snow_capacity));
			}
			/*--------------------------------------------------------------*/
			/*		Highest layer in patch.				*/
			/*--------------------------------------------------------------*/
            if ( close_enough(stratum[0].epv.height, patch[0].layers[0].height )){
				stratum[0].ga = 1.0 / compute_ra_overstory(
					command_line[0].verbose_flag,
					stratum[0].defaults[0][0].wind_attenuation_coeff,
					&(wind),
					zone[0].base_stations[0][0].screen_height,
					stratum[0].epv.height,
					layer[0].base,
					&(ga));
			}

			/*--------------------------------------------------------------*/
			/*		Layer is not the highest and is >0.1highest ht.	*/
			/*--------------------------------------------------------------*/
			else if (stratum[0].epv.height > (patch[0].layers[0].height * 0.1) ){
				if ( command_line[0].verbose_flag == -5 ){
					printf("\nCASE2");
				}
				windcan = wind;
				 stratum[0].ga = 1.0 / compute_ra_understory(
					command_line[0].verbose_flag,
				//	stratum[0].defaults[0][0].wind_attenuation_coeff * stratum[0].epv.proj_pai,
					stratum[0].defaults[0][0].wind_attenuation_coeff * stratum[0].epv.proj_pai_when_red, //NREN 20180804
					&(wind),
					stratum[0].epv.height,
					layer[0].base,
					&(ga));
				windsnow = wind;
				gasnow = ga;
			}
			/*--------------------------------------------------------------*/
			/*		Layer is <0.1highest ht. in height.		*/
			/*--------------------------------------------------------------*/
			else{
				if ( command_line[0].verbose_flag == -5 ){
					printf("\nCASE3");
				}
				windcan = wind;
				stratum[0].ga = 1.0 /
					compute_ra_surface(
					command_line[0].verbose_flag,
				//	stratum[0].defaults[0][0].wind_attenuation_coeff * stratum[0].epv.proj_pai,
					stratum[0].defaults[0][0].wind_attenuation_coeff * stratum[0].epv.proj_pai_when_red,  //NREN 20180804
					&(wind),
                    stratum[0].epv.height,
                    layer[0].base,
					&(ga));
				windsnow = wind;
				gasnow = ga;
			}
		}
		else{

			stratum[0].ga =
				pow(zone[0].wind * stratum[0].defaults[0][0].ustar_overu,2.0);
			ga = stratum[0].ga;
		}
	}
	else{
		stratum[0].ga = 0.0;
	}

	/*--------------------------------------------------------------*/
	/*	Factor in stability correction.				*/
	/*--------------------------------------------------------------*/
	stratum[0].ga = stratum[0].ga * patch[0].stability_correction ;
	if ( command_line[0].verbose_flag > 1 )
		printf("%8.6f ",stratum[0].ga);
	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8d -444.8 ",julday(current_date)-2449000);

	/*--------------------------------------------------------------*/
	/*	Determine non-vascular condductance to evaporation.		*/
	/*	This conductance represnets the inverse of an additional 	*/
	/*	resistance to vapour flux from the stratum rian storage		*/
	/*	surface over and above aerodynamic resistances that also	*/
	/*	affect turbulent heat transfer.  						 	*/
	/*	                            								*/
	/*	A linear relationship is currently assumed with the amount	*/
	/*	of relative rain stored - or with the patch unsat zone storage  */
	/*	relative to water equiv depth top water table if roots are present */
	/*	.  Parameters for the relationship	*/
	/*	are supplied via the stratum default file.					*/
	/*--------------------------------------------------------------*/
	if ( stratum[0].epv.all_pai > 0 ){
		stratum[0].gsurf = compute_nonvascular_stratum_conductance(
			command_line[0].verbose_flag,
			stratum[0].rain_stored+rain_throughfall,
		//	stratum[0].epv.all_pai      //20181101
			stratum[0].epv.all_pai_when_red //NREN 20180804
			* stratum[0].defaults[0][0].specific_rain_capacity,
			stratum[0].defaults[0][0].epc.gl_c,
			stratum[0].defaults[0][0].gsurf_slope,
			stratum[0].defaults[0][0].gsurf_intercept);
	}
	else{
		stratum[0].gsurf = compute_nonvascular_stratum_conductance(
			command_line[0].verbose_flag,
			patch[0].unsat_storage+rain_throughfall,
			patch[0].sat_deficit,
			stratum[0].defaults[0][0].epc.gl_c,
			stratum[0].defaults[0][0].gsurf_slope,
			stratum[0].defaults[0][0].gsurf_intercept);
	}

	if ( command_line[0].verbose_flag > 1 )
		printf("%8.6f ",stratum[0].gs);
	if ( command_line[0].verbose_flag > 1 )
		printf("\n%4d %4d %4d -444.9 \n ",
		current_date.day, current_date.month, current_date.year);
	if ( command_line[0].verbose_flag > 1 )
		printf("%8.4f %8.4f %8.4f %8.4f \n ",
		stratum[0].Kstar_direct,stratum[0].Kstar_diffuse,
		stratum[0].APAR_direct,stratum[0].APAR_diffuse);
	if ( command_line[0].verbose_flag > 1 ) {
		printf("\n%4d %4d %4d -444.9.1 \n ",current_date.day, current_date.month,
			current_date.year);
		printf(" %f\n", stratum[0].snow_stored);
	}
	/*--------------------------------------------------------------*/
	/*	Update snow storage ( this also updates the patch level	*/
	/*	snow throughfall and the stratum level Kstar )	 	*/
	/*	Snow on the canopy is first sublimated before any ET	*/
	/*	can take place.  This may seem like we are hobbling	*/
	/*	the energy budget since the sublimation reduces the	*/
	/*	energy used for stomatal conductance etc. however we	*/
	/*	suggest that when it snows it is likely ET is small.	*/
	/*--------------------------------------------------------------*/

	/* Lundberg 1994 reduce conductance for snow vs. rain by factor of 10 */
	stratum[0].snow_stored = compute_snow_stored(
		command_line[0].verbose_flag,
		zone[0].metv.tavg,
		zone[0].e_dewpoint,
		stratum[0].ga/10.0,
		zone[0].metv.pa,
		1000.0,
		&(snow_throughfall),
		&(rain_throughfall),
		stratum);

	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8.4f %f ",stratum[0].snow_stored, stratum[0].sublimation);
	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8d -444.10 ",julday(current_date)-2449000);
	if ( command_line[0].verbose_flag > 1 )
		printf("%8.4f %8.4f %8.4f %8.4f ",stratum[0].Kstar_direct,
		stratum[0].Kstar_diffuse,stratum[0].APAR_direct,
		stratum[0].APAR_diffuse);


	if (stratum[0].epv.proj_lai > ZERO)
		perc_sunlit = (stratum[0].epv.proj_lai_sunlit)/(stratum[0].epv.proj_lai);//LAI no need to change for EVP not considering dead lwaf 20180804
	else
		perc_sunlit = 1.0;

	if (stratum[0].epv.proj_lai_shade > ZERO && zone[0].metv.dayl > ZERO)
		stratum[0].ppfd_shade = (stratum[0].APAR_diffuse * (1-perc_sunlit)) /
                    zone[0].metv.dayl / max(1.0,stratum[0].epv.proj_lai_shade); //11022022LML add max(1.0)
    else
		stratum[0].ppfd_shade = 0.0;

	if (stratum[0].epv.proj_lai_sunlit > ZERO && zone[0].metv.dayl > ZERO)
		stratum[0].ppfd_sunlit = (stratum[0].APAR_direct + stratum[0].APAR_diffuse * perc_sunlit) /
                    zone[0].metv.dayl / max(1.0,stratum[0].epv.proj_lai_sunlit); //11022022LML add max(1.0)
	else
		stratum[0].ppfd_sunlit = 0.0;


    //printf("month:%d day:%d StrataID:%d APAR_direct(W/m2):%lf APAR_diffuse:%lf LAI_sunlit:%lf\n",
    //       (int)current_date.month,(int)current_date.day,stratum[0].ID
    //        ,stratum[0].APAR_direct*0.21739/zone[0].metv.dayl,stratum[0].APAR_diffuse*0.21739/zone[0].metv.dayl,
    //        stratum[0].epv.proj_lai_sunlit);

	/*--------------------------------------------------------------*/
	/*	Compute stratum conductance	m/s			*/
	/*								*/
	/*								*/
	/*	Note that the APAR supplied reduces to the rate of	*/
	/*	APAR (kJ/sec) absorbed during the daytime.  We do this	*/
	/*	since:							*/
	/*								*/
	/*	gs is only meaningful in the day			*/
	/*	gs responds to the instantaneous rather than daily	*/
	/*		averaged APAR.  				*/
	/*								*/
	/*	Also note that we will not void a shorwtwave energy	*/
	/*	budget since stratum conductance is used to scale 	*/
	/*	potnetial_evaporation; which we make sure is limited	*/
	/*	by the available energy. 				*/
	/*								*/
	/*	Note: PAR values are converted fro umol photon/m2*day to*/
	/*		umol photon/m2*sec				*/
	/*--------------------------------------------------------------*/
	/*
	if ( command_line[0].verbose_flag > 1 )
		printf("%f ",stratum[0].APAR_direct/zone[0].metv.dayl+
		stratum[0].APAR_diffuse/zone[0].metv.dayl);		*/
	/*--------------------------------------------------------------*/
	/*								*/
	/*	if there is still snow sitting on the canopy gs should be zero */
	/*								*/
	/*--------------------------------------------------------------*/

	if (stratum[0].snow_stored < ZERO)  {
	stratum[0].gs_sunlit = compute_vascular_stratum_conductance(
		command_line[0].verbose_flag,
		stratum[0].defaults[0][0].epc.psi_curve,
		stratum[0].defaults[0][0].epc.ppfd_coef,
		stratum[0].defaults[0][0].epc.gl_c,
		stratum[0].defaults[0][0].lai_stomatal_fraction,
		stratum[0].defaults[0][0].epc.psi_open,
		stratum[0].defaults[0][0].epc.psi_close,
		stratum[0].defaults[0][0].epc.psi_threshold,
		stratum[0].defaults[0][0].epc.psi_slp,
		stratum[0].defaults[0][0].epc.psi_intercpt,
		stratum[0].defaults[0][0].epc.gl_smax,
		stratum[0].defaults[0][0].epc.topt,
		stratum[0].defaults[0][0].epc.tcoef,
		stratum[0].defaults[0][0].epc.tmax,
		stratum[0].defaults[0][0].epc.vpd_open,
		stratum[0].defaults[0][0].epc.vpd_close,
		stratum[0].ppfd_sunlit,
		stratum[0].epv.proj_lai_sunlit,
		stratum[0].epv.psi,
		zone[0].metv.tsoil,
		zone[0].metv.tday,
		zone[0].metv.vpd, zone[0].CO2,
		stratum[0].defaults[0][0].epc.coef_CO2,
		stratum[0].ID,
		stratum, patch);

	m_APAR_sunlit = stratum[0].mult_conductance.APAR;
	m_tavg_sunlit = stratum[0].mult_conductance.tavg;
	m_LWP_sunlit = stratum[0].mult_conductance.LWP;
	m_CO2_sunlit = stratum[0].mult_conductance.CO2;
	m_tmin_sunlit = stratum[0].mult_conductance.tmin;
	m_vpd_sunlit = stratum[0].mult_conductance.vpd;

	/*
	stratum[0].potential_gs_sunlit = compute_vascular_stratum_conductance(
		command_line[0].verbose_flag,
		stratum[0].defaults[0][0].epc.psi_curve,
		stratum[0].defaults[0][0].epc.ppfd_coef,
		stratum[0].defaults[0][0].epc.gl_c,
		stratum[0].defaults[0][0].lai_stomatal_fraction,
		stratum[0].defaults[0][0].epc.psi_open,
		stratum[0].defaults[0][0].epc.psi_close,
		stratum[0].defaults[0][0].epc.psi_threshold,
		stratum[0].defaults[0][0].epc.psi_slp,
		stratum[0].defaults[0][0].epc.psi_intercpt,
		stratum[0].defaults[0][0].epc.gl_smax,
		stratum[0].defaults[0][0].epc.topt,
		stratum[0].defaults[0][0].epc.tcoef,
		stratum[0].defaults[0][0].epc.tmax,
		stratum[0].defaults[0][0].epc.vpd_open,
		stratum[0].defaults[0][0].epc.vpd_close,
		stratum[0].ppfd_sunlit,
		stratum[0].epv.proj_lai_sunlit,
		9999.0,
		zone[0].metv.tsoil,
		stratum[0].defaults[0][0].epc.topt,
		stratum[0].defaults[0][0].epc.vpd_open-1.0, zone[0].CO2,
		stratum[0].defaults[0][0].epc.coef_CO2,
		stratum[0].ID,
		stratum, patch);
	*/

	/* subsituting a simpler max conductance value */
	stratum[0].potential_gs_sunlit = stratum[0].defaults[0][0].epc.gl_smax *
			stratum[0].defaults[0][0].lai_stomatal_fraction * stratum[0].epv.proj_lai_sunlit;


	stratum[0].gs_shade = compute_vascular_stratum_conductance( // evp for vascular so no need to change lai 20181030
		command_line[0].verbose_flag,
		stratum[0].defaults[0][0].epc.psi_curve,
		stratum[0].defaults[0][0].epc.ppfd_coef,
		stratum[0].defaults[0][0].epc.gl_c,
		stratum[0].defaults[0][0].lai_stomatal_fraction,
		stratum[0].defaults[0][0].epc.psi_open,
		stratum[0].defaults[0][0].epc.psi_close,
		stratum[0].defaults[0][0].epc.psi_threshold,
		stratum[0].defaults[0][0].epc.psi_slp,
		stratum[0].defaults[0][0].epc.psi_intercpt,
		stratum[0].defaults[0][0].epc.gl_smax,
		stratum[0].defaults[0][0].epc.topt,
		stratum[0].defaults[0][0].epc.tcoef,
		stratum[0].defaults[0][0].epc.tmax,
		stratum[0].defaults[0][0].epc.vpd_open,
		stratum[0].defaults[0][0].epc.vpd_close,
		stratum[0].ppfd_shade,
		stratum[0].epv.proj_lai_shade,
		stratum[0].epv.psi,
		zone[0].metv.tsoil,
		zone[0].metv.tday,
		zone[0].metv.vpd, zone[0].CO2,
		stratum[0].defaults[0][0].epc.coef_CO2,
		stratum[0].ID,
		stratum, patch);

	m_APAR_shade = stratum[0].mult_conductance.APAR;
	m_tavg_shade = stratum[0].mult_conductance.tavg;
	m_LWP_shade = stratum[0].mult_conductance.LWP;
	m_CO2_shade = stratum[0].mult_conductance.CO2;
	m_tmin_shade = stratum[0].mult_conductance.tmin;
	m_vpd_shade = stratum[0].mult_conductance.vpd;



	/*
	stratum[0].potential_gs_shade = compute_vascular_stratum_conductance(
		command_line[0].verbose_flag,
		stratum[0].defaults[0][0].epc.psi_curve,
		stratum[0].defaults[0][0].epc.ppfd_coef,
		stratum[0].defaults[0][0].epc.gl_c,
		stratum[0].defaults[0][0].lai_stomatal_fraction,
		stratum[0].defaults[0][0].epc.psi_open,
		stratum[0].defaults[0][0].epc.psi_close,
		stratum[0].defaults[0][0].epc.psi_threshold,
		stratum[0].defaults[0][0].epc.psi_slp,
		stratum[0].defaults[0][0].epc.psi_intercpt,
		stratum[0].defaults[0][0].epc.gl_smax,
		stratum[0].defaults[0][0].epc.topt,
		stratum[0].defaults[0][0].epc.tcoef,
		stratum[0].defaults[0][0].epc.tmax,
		stratum[0].defaults[0][0].epc.vpd_open,
		stratum[0].defaults[0][0].epc.vpd_close,
		stratum[0].ppfd_shade,
		stratum[0].epv.proj_lai_shade,
		9999.0,
		zone[0].metv.tsoil,
		stratum[0].defaults[0][0].epc.topt,
		stratum[0].defaults[0][0].epc.vpd_open-1.0, zone[0].CO2,
		stratum[0].defaults[0][0].epc.coef_CO2,
		stratum[0].ID,
		stratum, patch);
	*/

	/* subsituting a simpler max conductance value */
	stratum[0].potential_gs_shade = stratum[0].defaults[0][0].epc.gl_smax *
			stratum[0].defaults[0][0].lai_stomatal_fraction * stratum[0].epv.proj_lai_shade;



	/* keep track of conductance multipliers actually used an indication of stress */
	stratum[0].mult_conductance.APAR = (m_APAR_shade*stratum[0].epv.proj_lai_shade +
		m_APAR_sunlit * stratum[0].epv.proj_lai_sunlit)
		/ (stratum[0].epv.proj_lai_shade+stratum[0].epv.proj_lai_sunlit);

	stratum[0].mult_conductance.tavg = (m_tavg_shade*stratum[0].epv.proj_lai_shade +
		m_tavg_sunlit * stratum[0].epv.proj_lai_sunlit)
		/ (stratum[0].epv.proj_lai_shade+stratum[0].epv.proj_lai_sunlit);

	stratum[0].mult_conductance.LWP = (m_LWP_shade*stratum[0].epv.proj_lai_shade +
		m_LWP_sunlit * stratum[0].epv.proj_lai_sunlit)
		/ (stratum[0].epv.proj_lai_shade+stratum[0].epv.proj_lai_sunlit);

	stratum[0].mult_conductance.CO2 = (m_CO2_shade*stratum[0].epv.proj_lai_shade +
		m_CO2_sunlit * stratum[0].epv.proj_lai_sunlit)
		/ (stratum[0].epv.proj_lai_shade+stratum[0].epv.proj_lai_sunlit);

	stratum[0].mult_conductance.tmin = (m_tmin_shade*stratum[0].epv.proj_lai_shade +
		m_tmin_sunlit * stratum[0].epv.proj_lai_sunlit)
		/ (stratum[0].epv.proj_lai_shade+stratum[0].epv.proj_lai_sunlit);

	stratum[0].mult_conductance.vpd = (m_vpd_shade*stratum[0].epv.proj_lai_shade +
		m_vpd_sunlit * stratum[0].epv.proj_lai_sunlit)
		/ (stratum[0].epv.proj_lai_shade+stratum[0].epv.proj_lai_sunlit);


	}
	else {
		stratum[0].gs_sunlit = 0.0;
		stratum[0].gs_shade = 0.0;
		stratum[0].potential_gs_sunlit = 0.0;
		stratum[0].potential_gs_shade = 0.0;
		}


	stratum[0].gs = stratum[0].gs_sunlit + stratum[0].gs_shade;

	/*--------------------------------------------------------------*/
	/*	Determine heat flux between stratum and surface.			*/
	/*	representative of the surface temperature.  Obviously   */
	/*	for strata NOT at the surface one would have 0 heat cap */
	/*	kJ/day							*/
	/*	We had the current input of rain as well.		*/
	/*--------------------------------------------------------------*/
	if (stratum[0].epv.height == 0) {
		stratum[0].surface_heat_flux =
			-1 * compute_surface_heat_flux(
			command_line[0].verbose_flag,
			stratum[0].snow_stored,
			stratum[0].rain_stored+patch[0].rain_throughfall,
			//stratum[0].epv.all_pai *
			stratum[0].epv.all_pai_when_red* //NREN 20180804 20181101
			stratum[0].defaults[0][0].specific_rain_capacity,
			zone[0].metv.tmin,
			zone[0].metv.tnightmax,
			zone[0].metv.tsoil,
			patch[0].soil_defaults[0][0].deltaz,
			stratum[0].defaults[0][0].min_heat_capacity,
			stratum[0].defaults[0][0].max_heat_capacity);

	}

	double surface_heat_flux_day = day_proportion * patch[0].surface_heat_flux;
	double surface_heat_flux_night = patch[0].surface_heat_flux - surface_heat_flux_day;

	if ( command_line[0].verbose_flag > 2 )
		printf("\n%8d -444.11 ",julday(current_date)-2449000);
	/*--------------------------------------------------------------*/
	/*	Compute evaporation and transpiration RATES (m/s)	*/
	/*	for daylight period .					*/
	/*	The rainy rate assumes a vpd of 10Pa.			*/
	/*	Note if surface heat flux makes evap negative we	*/
	/*	have condensation.  					*/
	/*--------------------------------------------------------------*/
	double rnet_evap_night = 1000 * (stratum[0].Lstar_night + surface_heat_flux_night) / nightlength;
	double rnet_evap_day = 1000 * (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar_day + surface_heat_flux_day) / daylength;

	if (rnet_evap_night < ZERO) rnet_evap_night = 0.0;
	if (rnet_evap_day < ZERO) rnet_evap_night = 0.0;
	rnet_evap = rnet_evap_night + rnet_evap_day;


	if ( command_line[0].verbose_flag == -5 ){
		printf("\n          BEFORE EVAP: rnet_evap=%lf Kstar=%lf Lstar=%lf surfheat=%lf",
				rnet_evap,
			   (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse)/86.4,
			   stratum[0].Lstar/86.4,
			   stratum[0].surface_heat_flux/86.4);
		}


	/*--------------------------------------------------------------*/
	/*	Estimate potential evap rates.				*/
	/*--------------------------------------------------------------*/
	if ((stratum[0].gsurf > ZERO) && (stratum[0].ga > ZERO) && (rnet_evap > ZERO)) {
		potential_evaporation_rate_night = penman_monteith(
				command_line[0].verbose_flag,
				zone[0].metv.tnight,
				zone[0].metv.pa,
				zone[0].metv.vpd_night,
				rnet_evap_night,
				1/stratum[0].gsurf,
				1/stratum[0].ga,
				2);
		potential_evaporation_rate_night = max(0.0, potential_evaporation_rate_night);

		potential_evaporation_rate_day = penman_monteith(
				command_line[0].verbose_flag,
				zone[0].metv.tday,
				zone[0].metv.pa,
				zone[0].metv.vpd_day,
				rnet_evap_day,
				1/stratum[0].gsurf,
				1/stratum[0].ga,
				2);
		potential_evaporation_rate_day = max(0.0, potential_evaporation_rate_day);

		// Daily weighted average rate
		potential_evaporation_rate = (night_proportion * potential_evaporation_rate_night) +
				(day_proportion * potential_evaporation_rate_day);

		potential_rainy_evaporation_rate_night = penman_monteith(
				command_line[0].verbose_flag,
				zone[0].metv.tnight,
				zone[0].metv.pa,
				0,
				rnet_evap_night,
				1/stratum[0].gsurf,
				1/stratum[0].ga,
				2);
		potential_rainy_evaporation_rate_night = max(0.0, potential_rainy_evaporation_rate_night);

		potential_rainy_evaporation_rate_day = penman_monteith(
				command_line[0].verbose_flag,
				zone[0].metv.tday,
				zone[0].metv.pa,
				0,
				rnet_evap_day,
				1/stratum[0].gsurf,
				1/stratum[0].ga,
				2);
		potential_rainy_evaporation_rate_day = max(0.0, potential_rainy_evaporation_rate_day);

		// Daily weighted average rate
		potential_rainy_evaporation_rate = (night_proportion * potential_rainy_evaporation_rate_night) +
				(day_proportion * potential_rainy_evaporation_rate_day);

	}
	 else {
		potential_rainy_evaporation_rate_day = 0.0;
		potential_rainy_evaporation_rate_night  = 0.0;
		potential_evaporation_rate_day=0.0;
		potential_evaporation_rate_night=0.0;
	}

	if ( command_line[0].verbose_flag > 2  )
		printf("\n%8d -444.12 ",julday(current_date)-2449000);
	/*--------------------------------------------------------------*/
	/*	Transpiration rate.					*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Reduce rnet for transpiration by ratio of lai to pai	*/
	/*--------------------------------------------------------------*/

	if (zone[0].metv.dayl > ZERO) {
		rnet_trans_sunlit = 1000 * ((stratum[0].Kstar_direct + (perc_sunlit)*stratum[0].Kstar_diffuse)
				+ perc_sunlit * (stratum[0].Lstar_day + surface_heat_flux_day) ) / daylength
						* ( stratum[0].epv.proj_lai / stratum[0].epv.proj_pai );

		rnet_trans_shade = 1000 * (( (1.0-perc_sunlit)*stratum[0].Kstar_diffuse)
				+ (1.0-perc_sunlit) * (stratum[0].Lstar_day + surface_heat_flux_day) ) / daylength
						* ( stratum[0].epv.proj_lai / stratum[0].epv.proj_pai );
	}
	else {
		rnet_trans_sunlit = 0.0;
		rnet_trans_shade = 0.0;
	}

	rnet_trans_shade = max(rnet_trans_shade, 0.0);
	rnet_trans_sunlit = max(rnet_trans_sunlit, 0.0);

    //if ((patch[0].ID == 64301) && (rnet_trans_sunlit > 0)) {
    //  printf("rnet_trans_sunlit:%lf rnet_trans_shade:%lf proj_lai:%lf",rnet_trans_sunlit*1000
    //         ,rnet_trans_shade*1000
    //         ,stratum[0].epv.proj_lai);
    //}

	if ( command_line[0].verbose_flag == -5 ){
		printf("\n          BEFORE TRANS: rnet_trans_sun=%lf rnet_trans_shade=%lf p_sunlit=%lf Kstar=%lf Lstar=%lf surfheat=%lf",
			   rnet_trans_sunlit,
			   rnet_trans_shade,
			   perc_sunlit,
			   (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse)/86.4,
			   stratum[0].Lstar/86.4,
			   stratum[0].surface_heat_flux/86.4);
	}


	if ( (rnet_trans_sunlit > ZERO ) &&
		(stratum[0].defaults[0][0].lai_stomatal_fraction > ZERO ) &&
		(stratum[0].gs_sunlit > ZERO) && ( stratum[0].ga > ZERO) ){
		transpiration_rate_sunlit = penman_monteith(
			command_line[0].verbose_flag,
			zone[0].metv.tday,
			zone[0].metv.pa,
			zone[0].metv.vpd_day,
			rnet_trans_sunlit,
			1/stratum[0].gs_sunlit,
			1/stratum[0].ga,
			2) ;
		potential_transpiration_rate_sunlit = penman_monteith(
			command_line[0].verbose_flag,
			zone[0].metv.tday,
			zone[0].metv.pa,
			zone[0].metv.vpd_day,
			rnet_trans_sunlit,
			1/stratum[0].potential_gs_sunlit,
			1/stratum[0].ga,
			2) ;
		}
	else{
		transpiration_rate_sunlit = 0.0;
		potential_transpiration_rate_sunlit = 0.0;
	}
	if ( (rnet_trans_shade > ZERO ) &&
		(stratum[0].defaults[0][0].lai_stomatal_fraction > ZERO ) &&
		(stratum[0].gs_shade > ZERO) && ( stratum[0].ga > ZERO) ){
		transpiration_rate_shade = penman_monteith(
			command_line[0].verbose_flag,
			zone[0].metv.tday,
			zone[0].metv.pa,
			zone[0].metv.vpd_day,
			rnet_trans_shade,
			1/stratum[0].gs_shade,
			1/stratum[0].ga,
			2) ;
		potential_transpiration_rate_shade = penman_monteith(
			command_line[0].verbose_flag,
			zone[0].metv.tday,
			zone[0].metv.pa,
			zone[0].metv.vpd_day,
			rnet_trans_shade,
			1/stratum[0].potential_gs_shade,
			1/stratum[0].ga,
			2) ;
	}
	else{
		transpiration_rate_shade = 0.0;
		potential_transpiration_rate_shade = 0.0;
	}

	transpiration_rate = transpiration_rate_sunlit +  transpiration_rate_shade;
	potential_transpiration_rate = potential_transpiration_rate_sunlit +  potential_transpiration_rate_shade;
    stratum[0].trans_sunlit = transpiration_rate_sunlit;
    stratum[0].trans_shade = transpiration_rate_shade;

	/*--------------------------------------------------------------*/
	/*	Compute potential evaporation of stratum. 		*/
	/*	Weighted by rain and non rain periods of the daytime	*/
	/*	m/day = m/s * (sec/day)					*/
	/*								*/
	/*	Note that Kstar is converted from Kj/m2*day to W/m2	*/
	/*--------------------------------------------------------------*/
	double potential_evaporation_night = (potential_evaporation_rate_night
					* (nightlength - rain_duration_night) )
					+ (potential_rainy_evaporation_rate_night
					* rain_duration_night);
	double potential_evaporation_day = (potential_evaporation_rate_day
						* (daylength - rain_duration_day) )
						+ (potential_rainy_evaporation_rate_day
						* rain_duration_day);
	stratum[0].potential_evaporation = potential_evaporation_night + potential_evaporation_day;

    //if (stratum[0].potential_evaporation > 0.01)
    //printf(" potential_evaporation:%lf potential_evaporation_night:%.1f potential_evaporation_day:%.1f daylength(hr):%.1f rain_duration_day:%.1f nightlength:%.1f rain_duration_night:%f potential_evaporation_rate_night:%.1f potential_rainy_evaporation_rate_night:%.1f potential_evaporation_rate_day:%.1f potential_rainy_evaporation_rate_day:%.1f gsurf:%lf gs_sunlit:%lf gs_shade:%lf gs:%lf ga:%lf rsurf:%lf rs:%lf ra:%lf rs_sun:%lf rs_shade:%lf rs_sun_pot:%lf rs_shade_pot:%lf\n"
    //        ,stratum[0].potential_evaporation * 1000
    //        ,potential_evaporation_night*1000
    //        ,potential_evaporation_day*1000
    //        ,daylength/3600
    //        ,rain_duration_day/3600
    //        ,nightlength/3600
    //        ,rain_duration_night/3600
    //        ,potential_evaporation_rate_night * 86400 * 1000
    //        ,potential_rainy_evaporation_rate_night * 86400 * 1000
    //        ,potential_evaporation_rate_day * 86400 * 1000
    //        ,potential_rainy_evaporation_rate_day * 86400 * 1000
    //        ,stratum[0].gsurf
    //        ,stratum[0].gs_sunlit
    //        ,stratum[0].gs_shade
    //        ,stratum[0].gs
    //        ,stratum[0].ga
    //        ,1./(1e-12+stratum[0].gsurf)
    //        ,1./(1e-12+stratum[0].gs)
    //        ,1./(1e-12+stratum[0].ga)
    //        ,1./(1e-12+stratum[0].gs_sunlit)
    //        ,1./(1e-12+stratum[0].gs_shade)
    //        ,1./(1e-12+stratum[0].potential_gs_sunlit)
    //        ,1./(1e-12+stratum[0].potential_gs_shade));

	/*--------------------------------------------------------------*/
	/*	If this stratum is on the surface and has a non-zero 	*/
	/*	rootin gdepth we assume it first takes water from soil. */
	/*	In this case we move some of the potential		*/
	/*	evaporation to porential transpiration.  The criteria 	*/
	/*	used to determine how much is based on the 		*/
	/*	estimated cap rise  during the day.		*/
	/*--------------------------------------------------------------*/
	if  ( (stratum[0].epv.height == 0 ) && ( stratum[0].rootzone.depth > 0 )
		&& (stratum[0].epv.proj_lai > ZERO) ){
		if ( stratum[0].potential_evaporation > ZERO ){
			transpiration = min(stratum[0].potential_evaporation,
				patch[0].cap_rise);
			potential_transpiration = min(stratum[0].potential_evaporation,
				patch[0].cap_rise);
			stratum[0].potential_evaporation -=	transpiration;
		}
		else{
			transpiration = 0.0;
			potential_transpiration = 0.0;
		}
	}
	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8d -444.13 ",julday(current_date)-2449000);
	if ( command_line[0].verbose_flag > 1 )
		printf("%8.6f %8.6f %8.6f ",potential_evaporation_rate*1000,
		transpiration_rate*1000, stratum[0].potential_evaporation*1000 );
	if ( command_line[0].verbose_flag > 2 )
		printf("\n%4d %4d %4d -444.14 ",current_date.day, current_date.month,
		current_date.year);
	/*--------------------------------------------------------------*/
	/*	Update rain storage ( this also updates the patch level	*/
	/*	rain_throughfall and stratum[0].potential_evaporation	*/
	/*--------------------------------------------------------------*/

	if ( command_line[0].verbose_flag == -5 ){
		printf("\n          BEFORE RAINSTORED: evap=%lf trans=%lf potevaprate=%lf potrainevaprate=%lf potevap=%lf DAYL?=%lf DRD?=%lf MULT1?=%lf MULT2?=%lf",
				   stratum[0].evaporation,
				   stratum[0].transpiration_unsat_zone + stratum[0].transpiration_sat_zone,
				   potential_evaporation_rate*1000.0,
				   potential_rainy_evaporation_rate*1000.0,
				   stratum[0].potential_evaporation,
				zone[0].metv.dayl,
			   zone[0].rain_duration,
			   (zone[0].metv.dayl - (zone[0].rain_duration * zone[0].metv.dayl/86400) ),
			   (zone[0].rain_duration * zone[0].metv.dayl/86400)  );
		}


	stratum[0].rain_stored  = compute_rain_stored(
		command_line[0].verbose_flag,
		&(rain_throughfall),
		stratum);

    //if (patch[0].ID == 8066)
    //    printf("stratum[0].ID:%d cf:%.2f rain_stored(mm):%lf\n",stratum[0].ID,stratum[0].cover_fraction,stratum[0].rain_stored);

	if (stratum[0].rain_stored > 0){
	    NO3_stored = (stratum[0].rain_stored + stratum[0].snow_stored)
	      	/ (stratum[0].rain_stored + stratum[0].snow_stored + rain_throughfall + snow_throughfall)
		* (stratum[0].NO3_stored + patch[0].NO3_throughfall);
	    NO3_throughfall = (rain_throughfall + snow_throughfall)
		/ (stratum[0].rain_stored + stratum[0].snow_stored + rain_throughfall + snow_throughfall)
		* (stratum[0].NO3_stored + patch[0].NO3_throughfall);

	}
	else{
	    if (rain_throughfall > 0){
		NO3_stored = (stratum[0].rain_stored + stratum[0].snow_stored)
	      	/ (stratum[0].rain_stored + stratum[0].snow_stored + rain_throughfall + snow_throughfall)
		* (stratum[0].NO3_stored + patch[0].NO3_throughfall);

		NO3_throughfall =  (rain_throughfall + snow_throughfall)
		/ (stratum[0].rain_stored + stratum[0].snow_stored + rain_throughfall + snow_throughfall)
		* (stratum[0].NO3_stored + patch[0].NO3_throughfall);

	    }
	    else{
                NO3_stored = stratum[0].NO3_stored + patch[0].NO3_throughfall;
		NO3_throughfall = 0;


	    }
	}

	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8d -444.15 ",julday(current_date)-2449000);
	if ( command_line[0].verbose_flag > 1 )
		printf("%8.4f %8.4f %8.4f %8.4f ",stratum[0].Kstar_direct,
		stratum[0].Kstar_diffuse,stratum[0].APAR_direct,
		stratum[0].APAR_diffuse);
	/*--------------------------------------------------------------
	 *	Separate the evaporation into rainy and dry
	 *	assuming rainy happens as much as it can.
	 *
	 *  Note:  For transpiration, we only care about
	 *  evaporation during daytime hours, not the entire 24-hour period.
	 *  However, we are now accounting for daytime and nighttime
	 *  potential evaporation rates separately, but stratum->evaporation
	 *  is the entire 24-hour flux.   This is not consistent, and may be
	 *  problematic.
	 *--------------------------------------------------------------*/
	if ( stratum[0].evaporation > ZERO ){
		rainy_evaporation =  min((zone[0].rain_duration * zone[0].metv.dayl/86400) *
			potential_rainy_evaporation_rate,
			stratum[0].evaporation );
		dry_evaporation = stratum[0].evaporation - rainy_evaporation;
	}
	else{
		rainy_evaporation = 0;
		dry_evaporation = stratum[0].evaporation ;
	}

	/*--------------------------------------------------------------*/
	/*	Compute Canopy transpiration.	m/day			*/
	/*								*/
	/*	We assume the potential_evaporation has been 		*/
	/*	reduced by the amount need to evaporate the storage 	*/
	/*	and incident precipitation.  				*/
	/*								*/
	/*	we  assume that part of the day				*/
	/*	which had a dry canopy which was transpiring.		*/
	/*								*/
	/*	The advantage of doing it via potentia evaporation 	*/
	/*	include:						*/
	/*								*/
	/*	1. the potential evaporation is computed once using	*/
	/*	the vpd and Rnet and ra of the WHOLE day.  We do not	*/
	/*	reduce Rnet after we evaporate water stored but		*/
	/*	reduce the potential evap. directly.			*/
	/*								*/
	/*	2. more importantly, the leaf conductance is computed	*/
	/*	on the full value average daytime Rnet and the 		*/
	/*	Rnet used in transpiration is the daily averaged value	*/
	/*	which properly reflects the balance in demand between	*/
	/*	Rnet and vpd.						*/
	/*								*/
	/*	we assume no transpiration happens during evaporation.	*/
	/*	or rain hours.						*/
	/*								*/
	/*	All of this is only done for vascular strata.		*/
	/*--------------------------------------------------------------*/
	if  ((stratum[0].defaults[0][0].lai_stomatal_fraction > ZERO) && (stratum[0].epv.proj_lai > ZERO)){
		if ( stratum[0].potential_evaporation > ZERO ){
            double t_trans = max(0.0,(zone[0].metv.dayl - (zone[0].rain_duration * zone[0].metv.dayl/86400) -
                    dry_evaporation / potential_evaporation_rate));
            transpiration  = transpiration_rate * t_trans;
            potential_transpiration  = potential_transpiration_rate * t_trans;

            //if (potential_transpiration_rate*86400*1000 > 20)
            //    printf("trans:%lf p_trans:%f dry_evap:%f p_evap_rate:%f p_trans_rate:%f dry_evap:%f frain:%f ftrans:%f\n"
            //           ,transpiration*1000
            //           ,potential_transpiration*1000
            //           ,dry_evaporation*1000
            //           ,(potential_evaporation_rate*86400*1000)
            //           ,(potential_transpiration_rate*86400*1000)
            //           ,dry_evaporation*1000
            //           ,zone[0].rain_duration/86400
            //           ,t_trans/zone[0].metv.dayl);


		}
		else{
			transpiration  = 0.0;
			potential_transpiration = 0.0;
		}
	}
	transpiration = max(transpiration, 0.0); // MCK, test
	potential_transpiration = max(potential_transpiration, 0.0);

	stratum[0].PET = potential_transpiration;

	/*--------------------------------------------------------------*/
	/*	Separate the transpiration into unsat and sat zone	*/
	/*	transpiration demands based on the rooting depth and	*/
	/*	the patch depth to water table.				*/
	/*								*/
	/*	This transpiration may actually be reduced if there is	*/
	/*	not enough water i nthe respective zone.  But that 	*/
	/*	event should indicate that maybe the leaf water 	*/
	/*	potential control on gs or the soil moisture control on */
	/*	gsurf (for bryophytes) is not well calibrated.		*/
	/*--------------------------------------------------------------*/
	if ( stratum[0].rootzone.depth > ZERO ){
		stratum[0].transpiration_sat_zone = transpiration
			* max(0, 1 - ( patch[0].sat_deficit_z
			/ stratum[0].rootzone.depth ) );
		stratum[0].transpiration_unsat_zone = transpiration
			- stratum[0].transpiration_sat_zone;
	}
	else{
		if ( patch[0].sat_deficit_z > ZERO ){
			stratum[0].transpiration_unsat_zone = transpiration;
			stratum[0].transpiration_sat_zone = 0;
		}
		else{
			stratum[0].transpiration_unsat_zone = 0;
			stratum[0].transpiration_sat_zone = transpiration;
		}
	}

		/*--------------------------------------------------------------*/
		/* Remove energy used for ET from Kstar and APAR */
		/*--------------------------------------------------------------*/

		Rnet_used = (stratum[0].evaporation + stratum[0].transpiration_unsat_zone
						+ stratum[0].transpiration_sat_zone) * lhvap * 1000;
		Rnet_canopy = (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse
					    + stratum[0].Lstar + stratum[0].surface_heat_flux)
						- Rnet_used;

		/* Remove energy for tracking purposes, although not accurate since we are */
		/* transferring this energy into sensible heat. */
		/* ASSUMES ONLY NEGATIVE FLUX CAN BE LSTAR... SURFACE HEAT?? */
		if (Rnet_used > 0.0) {
			if ( (stratum[0].Lstar > 0) && (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar + stratum[0].surface_heat_flux > 0) ){
				fraction_L_used = stratum[0].Lstar
					/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar + stratum[0].surface_heat_flux);
				if (stratum[0].Lstar_night > 0.0) {
					fraction_L_used_night = stratum[0].Lstar_night
							/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar + stratum[0].surface_heat_flux);
				}
				if (stratum[0].Lstar_day > 0.0) {
					fraction_L_used_day = stratum[0].Lstar_day
							/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar + stratum[0].surface_heat_flux);
				}
				fraction_direct_K_used = stratum[0].Kstar_direct
					/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar + stratum[0].surface_heat_flux);
				fraction_diffuse_K_used = stratum[0].Kstar_diffuse
					/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar + stratum[0].surface_heat_flux);
				fraction_surfheat_used = stratum[0].surface_heat_flux
					/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].Lstar + stratum[0].surface_heat_flux);
			}
			else if ( stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].surface_heat_flux > 0.0 ) {
				fraction_direct_K_used = stratum[0].Kstar_direct
					/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].surface_heat_flux);
				fraction_diffuse_K_used = stratum[0].Kstar_diffuse
					/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].surface_heat_flux);
				fraction_surfheat_used = stratum[0].surface_heat_flux
					/ (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse + stratum[0].surface_heat_flux);
			}
			if ( stratum[0].APAR_direct > 0 ){
				fraction_direct_APAR_used = 1.0
					/ ( 1 + stratum[0].APAR_diffuse/stratum[0].APAR_direct);
			}
			else
				fraction_direct_APAR_used = 0.0;
			if (stratum[0].APAR_diffuse > 0 ){
				fraction_diffuse_APAR_used = 1.0
					/ ( 1 + stratum[0].APAR_direct/stratum[0].APAR_diffuse);
			}
			else
				fraction_diffuse_APAR_used = 0.0;

			APAR_used = ( ( (fraction_direct_K_used * Rnet_used) + (fraction_diffuse_K_used * Rnet_used) )
						 / ( stratum[0].Kstar_direct + stratum[0].Kstar_diffuse ) )
						* (stratum[0].APAR_direct + stratum[0].APAR_diffuse);
			stratum[0].Lstar -= Rnet_used * fraction_L_used;
			stratum[0].Lstar_night -= Rnet_used * fraction_L_used_night;
			stratum[0].Lstar_day -= Rnet_used * fraction_L_used_day;
			stratum[0].Kstar_direct -= Rnet_used * fraction_direct_K_used;
			stratum[0].Kstar_diffuse -= Rnet_used * fraction_diffuse_K_used;
			stratum[0].surface_heat_flux -= Rnet_used * fraction_surfheat_used;
			stratum[0].APAR_direct -= APAR_used * fraction_direct_APAR_used;
			stratum[0].APAR_diffuse -= APAR_used * fraction_diffuse_APAR_used;

			/* Zero out negative APARs. K and Lstars ok since we are converting to temp change in canopy daily F */
			stratum[0].APAR_direct = max(stratum[0].APAR_direct,0.0);
			stratum[0].APAR_diffuse = max(stratum[0].APAR_diffuse,0.0);

		} /* end if Rnet_used > 0 */

		if ( stratum[0].APAR_direct < -1 ) {
			printf("CANOPY: APARused=%lf APARdir=%lf APAR=%lf Rnet_used=%lf Kstardir=%lf Kstar=%lf Lstar=%lf \n",
				   APAR_used,
				   stratum[0].APAR_direct,
				   (stratum[0].APAR_direct+stratum[0].APAR_diffuse)/86.4,
				   Rnet_used/86.4,
				   stratum[0].Kstar_direct/86.4,
				   (stratum[0].Kstar_direct+stratum[0].Kstar_diffuse)/86.4,
				   stratum[0].Lstar/86.4);
		}


		/* Now assume remaining energy is balanced by sensible heat exchange. */
		/* Reverse calculate Tcanopy. */
		/* This should only affect canopy longwave for snow and surface fluxes. */
		deltaT = (Rnet_canopy * 1000 / 86400) /* energy in J/m2/s */
					/ (1.292 - ( 0.00428 * zone[0].metv.tavg ))   /* density of air in kg/m3 */
					/ CP    /* heat capacity of air in J/kg/C */
					* (1/stratum[0].ga);   /* aerodynamic resistance in s/m */

		/* Empirical reduction in deltaT to account for feedbacks from changes in other */
		/* components of the energy balance (mostly longwave). Tested with a variety of */
		/* daily climate values and solving iteratively for steady-state canopy temp    */
		/* and values range from 0.75-0.95 with a mean of 0.85. */
		deltaT *= 0.85;
		/*deltaT = 0.0;*/


        if ( command_line[0].verbose_flag == -5 ){
		printf("\n          ENERGY REDUCTION: evap=%lf trans=%lf potevaprate=%lf potrainevaprate=%lf potevap=%lf totalrad=%lf Rnet=%lf Rnet_used=%lf tavg=%lf\n          Kdir=%lf Kdif=%lf Lstar=%lf surfheat=%lf\n          Lstarpch=%lf K_refl=%lf PAR_refl=%lf snowstor1000=%lf rainstor1000=%lf wetfrac=%lf snowcap=%lf raindur=%lf rnet_post=%lf deltaT=%lf",
			   stratum[0].evaporation,
			   stratum[0].transpiration_unsat_zone + stratum[0].transpiration_sat_zone,
			   potential_evaporation_rate,
			   potential_rainy_evaporation_rate,
			   stratum[0].potential_evaporation,
			   (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse
					+ stratum[0].Lstar + stratum[0].surface_heat_flux)/86.4,
			   rnet_evap,
			   Rnet_used/86.4,
			   zone[0].metv.tavg,
			   stratum[0].Kstar_direct/86.4,
			   stratum[0].Kstar_diffuse/86.4,
			   stratum[0].Lstar/86.4,
			   stratum[0].surface_heat_flux/86.4,
			   patch[0].Lstar_canopy/86.4,
			   K_reflectance,
			   PAR_reflectance,
			   stratum[0].snow_stored*1000,
			   stratum[0].rain_stored*1000,
			   wetfrac,
			   (stratum[0].epv.all_pai * stratum[0].defaults[0][0].specific_snow_capacity),
			   zone[0].rain_duration * zone[0].metv.dayl/86400,
			   Rnet_canopy,
			   deltaT);
        }




		/*--------------------------------------------------------------*/

	if ( command_line[0].verbose_flag > 1 )
		printf("\n%8d -444.16 ",julday(current_date)-2449000);
	if ( command_line[0].verbose_flag > 1 )
		printf("%8.6f ", transpiration);
	/*--------------------------------------------------------------*/
	/*	Do respiration and photosynthesis only for plants	*/
	/*--------------------------------------------------------------*/
	if (stratum[0].defaults[0][0].lai_stomatal_fraction > ZERO ) {
		/*--------------------------------------------------------------*/
		/*	perform maintenance respiration				*/
		/*	only fluxes are computed here; stores are updated later	*/
		/*	in canopy_stratum_daily_growth				*/
		/*--------------------------------------------------------------*/
		if (compute_maint_resp(
			stratum[0].defaults[0][0].mrc.q10,
			stratum[0].defaults[0][0].mrc.per_N,
			&(stratum[0].cs),
			&(stratum[0].ns),
			&(stratum[0].defaults[0][0].epc),
			&(zone[0].metv),
			&(stratum[0].cdf)	)){
			fprintf(stderr,"Error in compute_maint_resp() from bbgc.c... Exiting\n");
			exit(EXIT_FAILURE);
		}
		if ((stratum[0].epv.all_lai > ZERO) && (stratum[0].snow_stored < ZERO))  {
			/*--------------------------------------------------------------*/
			/*	convert maintenance resp from kg/m2*day to umol/m2*s	*/
			/*	only fluxes are computed here; stores are updated later	*/
			/*	in canopy_stratum_daily_growth				*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*	Set up input array for stratum assimilation		*/
			/*--------------------------------------------------------------*/

			/*--------------------------------------------------------------*/
			/* this all needs to be repeated for sunlit and for shaded 	*/
			/*	and for potential and actual psn`			*/
			/* note that at present potential only takes into account	*/
			/*	water (i.e not nitrogen) limitations			*/
			/*--------------------------------------------------------------*/
			/* potential sunlit psn						*/
			/*--------------------------------------------------------------*/
			if (stratum[0].defaults[0][0].epc.veg_type == C4GRASS)
				psnin.c3 = 0;
			else
				psnin.c3 = 1;
			if (zone[0].metv.dayl > ZERO)
				psnin.Rd = stratum[0].cdf.leaf_day_mr /
				(stratum[0].epv.proj_lai
				* zone[0].metv.dayl*12.011e-9);
			else
				psnin.Rd = 0.0;

			psnin.pa = zone[0].metv.pa;
			psnin.co2 = zone[0].CO2;
			psnin.flnr = stratum[0].defaults[0][0].epc.flnr;
			psnin.t = zone[0].metv.tday;
			psnin.irad = stratum[0].ppfd_sunlit;
			if ((stratum[0].cs.leafc > ZERO) && (stratum[0].epv.proj_sla_sunlit > ZERO))
				psnin.lnc = stratum[0].ns.leafn / (stratum[0].cs.leafc * 1.0)
					/ stratum[0].epv.proj_sla_sunlit;
			else
				psnin.lnc = 0.0;
			/*--------------------------------------------------------------*/
			/* note multiply by 1000; accounted for in compute_farq_psn by a /1000 */
			/*	this is done for numerical precision			*/
			/*--------------------------------------------------------------*/
			psnin.g = max(stratum[0].defaults[0][0].epc.gl_smax ,
				stratum[0].defaults[0][0].epc.gl_c ) *
				stratum[0].defaults[0][0].lai_stomatal_fraction * 1000 / 1.6;
			if ((psnin.lnc > 0.0) && (psnin.irad > 0.0))
				compute_farq_psn(&psnin, &psnout, 1);
			else
				psnout.A = 0.0;
			assim_sunlit = psnout.A;

            //printf("A:%f g:%f lc:%f irad:%f pa:%f t:%f\n",psnout.A,psnin.g,
            //       psnin.lnc,psnin.irad,psnin.pa,psnin.t);

            double assim_sunlit_pot = assim_sunlit;


			/*--------------------------------------------------------------*/
			/* potential shade psn						*/
			/*--------------------------------------------------------------*/
			if (zone[0].metv.dayl > ZERO)
				psnin.Rd = stratum[0].cdf.leaf_day_mr  /
				(stratum[0].epv.proj_lai
				* zone[0].metv.dayl*12.011e-9);
			else
				psnin.Rd = 0.0;

			psnin.irad = stratum[0].ppfd_shade;
			if ((stratum[0].cs.leafc > ZERO) && (stratum[0].epv.proj_sla_shade > ZERO))
				psnin.lnc = stratum[0].ns.leafn / (stratum[0].cs.leafc * 1.0)
					/ stratum[0].epv.proj_sla_shade ;
			else
				psnin.lnc = 0.0;

			if ((psnin.lnc > 0.0) && (psnin.irad > 0.0))
				compute_farq_psn(&psnin, &psnout, 1);
			else
				psnout.A = 0.0;
			assim_shade = psnout.A;
            double assim_shade_pot = assim_shade;

			/*--------------------------------------------------------------*/
			/* total potential psn						*/
			/*	We add the daytime leaf mr to the assim since it has	*/
			/*	been subtracted in both the farq and accounted for in	*/
			/*	total mr.						*/
			/*--------------------------------------------------------------*/
			stratum[0].cdf.potential_psn_to_cpool = (assim_sunlit*stratum[0].epv.proj_lai_sunlit
					+ assim_shade * stratum[0].epv.proj_lai_shade)
					*zone[0].metv.dayl*12.011e-9 + stratum[0].cdf.leaf_day_mr;

			/*--------------------------------------------------------------*/
			/* actual sunlit psn						*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*	convert water vapour to co2 conductance.		*/
			/*--------------------------------------------------------------*/
			if (zone[0].metv.dayl > ZERO)
				psnin.Rd = stratum[0].cdf.leaf_day_mr /
				(stratum[0].epv.proj_lai
				* zone[0].metv.dayl*12.011e-9);
			else
				psnin.Rd = 0.0;
			/*--------------------------------------------------------------*/
			/* note multiply by 1000; accounted for in compute_farq_psn by a /1000 */
			/*	this is done for numerical precision			*/
			/*--------------------------------------------------------------*/
			if (stratum[0].epv.proj_lai_sunlit > ZERO)
				psnin.g = stratum[0].gs_sunlit * 1000 / 1.6 / stratum[0].epv.proj_lai_sunlit;
			else
				psnin.g = 0.0;
			psnin.irad = stratum[0].ppfd_sunlit;
			if ((stratum[0].cs.leafc > ZERO) && (stratum[0].epv.proj_lai_sunlit > ZERO))
				psnin.lnc = stratum[0].ns.leafn / (stratum[0].cs.leafc * 1.0)
					/ stratum[0].epv.proj_sla_sunlit;
			else
				psnin.lnc = 0.0;
			if ((psnin.lnc > 0.0) && (psnin.irad > 0.0)) {
				if ( compute_farq_psn(&psnin, &psnout, 1)){
					fprintf(stderr,
						"FATAL ERROR: in canopy_stratum_daily_F error in farquhar");
					exit(EXIT_FAILURE);
				}
			}
			else
				psnout.A = 0.0;
			assim_sunlit = psnout.A;
			dC13_sunlit = psnout.dC13;
            stratum[0].assim_sunlit = assim_sunlit;

            //if (assim_sunlit <= 0.0 && stratum[0].cs.leafc > 0) {
            //    printf("\ny:%d m:%d d:%d assim_sunlit:%lf leafc:%lf leafn:%lf lai_sun:%lf\n"
            //           ,current_date.year
            //           ,current_date.month
            //           ,current_date.day
            //           ,assim_sunlit
            //           ,stratum[0].cs.leafc*1000
            //           ,stratum[0].ns.leafn*1000
            //           ,stratum[0].epv.proj_lai_sunlit
            //           );
            //}

			/*--------------------------------------------------------------*/
			/* actual shade psn						*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*	convert water vapour to co2 conductance.		*/
			/*--------------------------------------------------------------*/
			if (zone[0].metv.dayl > ZERO)
			psnin.Rd = stratum[0].cdf.leaf_day_mr  /
				(stratum[0].epv.proj_lai
				* zone[0].metv.dayl*12.011e-9);
			else
				psnin.Rd = 0.0;
			/*--------------------------------------------------------------*/
			/* note multiply by 1000; accounted for in compute_farq_psn by a /1000 */
			/*	this is done for numerical precision			*/
			/*--------------------------------------------------------------*/
			if (stratum[0].epv.proj_lai_shade > ZERO)
				psnin.g = stratum[0].gs_shade *1000 / 1.6 / stratum[0].epv.proj_lai_shade;
			else
				psnin.g = 0.0;
			psnin.irad = stratum[0].ppfd_shade;
			if ((stratum[0].cs.leafc > ZERO) && (stratum[0].epv.proj_lai_shade > ZERO))
				psnin.lnc = stratum[0].ns.leafn / (stratum[0].cs.leafc * 1.0)
					/ stratum[0].epv.proj_sla_shade;
			else
				psnin.lnc = 0.0;
			if ((psnin.lnc > 0.0) && (psnin.irad > 0.0)) {
				if ( compute_farq_psn(&psnin, &psnout, 1)){
					fprintf(stderr,
						"FATAL ERROR: in canopy_stratum_daily_F error in farquhar");
					exit(EXIT_FAILURE);
				}
			}
			else
				psnout.A = 0.0;
			assim_shade = psnout.A;
			dC13_shade = psnout.dC13;
            stratum[0].assim_shade = assim_shade;
			/*--------------------------------------------------------------*/
			/* total actual psn						*/
			/*--------------------------------------------------------------*/
			stratum[0].cdf.psn_to_cpool = (assim_sunlit*stratum[0].epv.proj_lai_sunlit
					+ assim_shade * stratum[0].epv.proj_lai_shade)
					*zone[0].metv.dayl*12.011e-9 + stratum[0].cdf.leaf_day_mr;


            //if (stratum[0].ID == 813381 && assim_sunlit_pot <= 0.){
            //    printf("mon:%d\tday:%d\tID:%d\tpsn_to_cpool(gC/m2/day):%lf\tassim_sunlit_pot:%lf\tassim_sunlit:%lf\tlai_sunlit:%lf\tppfd_sunlit:%llf\tassim_shade_pot:%lf\tassim_shade:%lf\tlai_shade:%lf\tppfd_shade:%lf\tdayl:%f\tday_mr:%lf\tleafc:%lf\tleafn:%lf\tproj_sla_sunlit:%lf\tproj_sla_shade:%lf\tgs_sunlit:%lf\tgs_shade:%lf\n",
            //            current_date.month,current_date.day,stratum[0].ID,stratum[0].cdf.psn_to_cpool*1000.0,assim_sunlit_pot,assim_sunlit,stratum[0].epv.proj_lai_sunlit,stratum[0].ppfd_sunlit
            //            ,assim_shade_pot,assim_shade,stratum[0].epv.proj_lai_shade,stratum[0].ppfd_shade,zone[0].metv.dayl,stratum[0].cdf.leaf_day_mr
            //            ,stratum[0].cs.leafc,stratum[0].ns.leafn,stratum[0].epv.proj_sla_sunlit,stratum[0].epv.proj_sla_shade
            //            ,stratum[0].gs_sunlit,stratum[0].gs_shade);
            //}



			if ((assim_sunlit + assim_shade) > ZERO)
				stratum[0].dC13 = (assim_sunlit * dC13_sunlit + assim_shade * dC13_shade)/(assim_sunlit+assim_shade);
			else
				stratum[0].dC13 = 0.0;
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
		} /* end if LAI > O  ** snow stored < 0*/

		else {
			stratum[0].cdf.psn_to_cpool = 0.0;
			stratum[0].cdf.potential_psn_to_cpool = 0.0;
		}


	} /* end if stomatal_fraction > 0 */
	else {
		stratum[0].cdf.psn_to_cpool = 0.0;
		stratum[0].cdf.potential_psn_to_cpool = 0.0;
	}

	/*--------------------------------------------------------------*/
	/*	perform growth related computations (if grow flag is on */
	/*--------------------------------------------------------------*/
	if ((command_line[0].grow_flag > 0) && (stratum[0].defaults[0][0].epc.veg_type != NON_VEG)) {
		/*--------------------------------------------------------------*/
		/*	compute N uptake from the soil 				*/
		/*--------------------------------------------------------------*/

		switch(stratum[0].defaults[0][0].epc.allocation_flag) {

		case CONSTANT: /* constant allocation */
			stratum[0].ndf.potential_N_uptake =	compute_potential_N_uptake(
				stratum[0].defaults[0][0].epc,
				&(stratum[0].epv),
				&(stratum[0].cs),
				&(stratum[0].ns),
				&(stratum[0].cdf));
			break;
		case DICKENSON:
			stratum[0].ndf.potential_N_uptake =compute_potential_N_uptake_Dickenson(
				stratum[0].defaults[0][0].epc,
				&(stratum[0].epv),
				&(stratum[0].cs),
				&(stratum[0].ns),
				&(stratum[0].cdf));
			break;
		case WARING:
			stratum[0].ndf.potential_N_uptake =compute_potential_N_uptake_Waring(
				stratum[0].defaults[0][0].epc,
				&(stratum[0].epv),
				&(stratum[0].cs),
				&(stratum[0].ns),
				&(stratum[0].cdf));
			break;
		case COMBINED:
			stratum[0].ndf.potential_N_uptake =compute_potential_N_uptake_combined(
				stratum[0].defaults[0][0].epc,
				&(stratum[0].epv),
				&(stratum[0].cs),
				&(stratum[0].ns),
				&(stratum[0].cdf));
			break;
		} /* end switch */
	}

	}
	/*--------------------------------------------------------------*/
	/*	Increment the transmitted fluxes from this patch layer	*/
	/*	by weighting the fluxes in this stratum by its cover	*/
	/*	fraction - we have check cover fractions sum to 1 in 	*/
	/*	a layer when constructing the patch.			*/
	/*--------------------------------------------------------------*/

    //10102023LML: there is a problem for layer structure of strata.
    //The cover_fraction should not be stratum fraction over the patch,
    //but just be the fraction over this layer
    double stratum_cover_fraction_layer = stratum[0].cover_fraction;
    if (layer->count > 1) {
        double total_cover = 0.;
        for (int s=0 ; s<layer->count; s++) {
            total_cover += patch[0].canopy_strata[layer->strata[s]]->cover_fraction;
        }
        if (!close_enough(total_cover,0))
            stratum_cover_fraction_layer = (1. - layer->null_cover) * stratum_cover_fraction_layer / total_cover;
        else stratum_cover_fraction_layer = 0;
    }
    //for ( stratum=0 ; stratum<patch[0].layers[layer].count; stratum++ ){


    patch[0].Kdown_direct_final += Kdown_direct * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].Kup_direct_final += Kup_direct * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].PAR_direct_final += PAR_direct * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].Kdown_diffuse_final += Kdown_diffuse * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].Kup_diffuse_final += Kup_diffuse * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].PAR_diffuse_final += PAR_diffuse * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].Ldown_final += patch[0].Ldown * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].Kstar_canopy_final += (stratum[0].Kstar_direct + stratum[0].Kstar_diffuse) * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].LE_canopy_final += Rnet_used * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
	patch[0].rain_throughfall_final += rain_throughfall
        * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
	patch[0].snow_throughfall_final += snow_throughfall
        * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
	patch[0].NO3_throughfall_final += NO3_throughfall
        * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
	stratum[0].NO3_stored = NO3_stored;

    //if (patch[0].ID == 7918)
    //  printf("canopy_stratum_daily_F.c strata:%d rain_throughfall_final:%lf\n"
    //        ,stratum[0].ID,patch[0].rain_throughfall_final*1000);

	// // calculate et just like in patch, but also by cover fraction, a variable we include MCK, patch - level variable multiplied by stratum cover fraction, just like those variables above

	patch[0].ga_final = ga;
	patch[0].gasnow_final = gasnow;
    patch[0].wind_final += wind * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].windsnow_final += windsnow * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].ustar_final += ustar * stratum_cover_fraction_layer; //stratum[0].cover_fraction;
    patch[0].T_canopy_final += (zone[0].metv.tavg + deltaT) * stratum_cover_fraction_layer; //stratum[0].cover_fraction;

	/* track variables for fire spread */
	if (command_line[0].firespread_flag == 1) {
//		printf("***********\n###########\nCurrent stratum height %lf\n***************\n#################\n",stratum[0].epv.height);
		if(stratum[0].epv.height<=patch[0].soil_defaults[0][0].understory_height_thresh)
		{
//			printf("Found an understory stratum!\n");
            //patch[0].fire.understory_et = (patch[0].fire_defaults[0][0].ndays_average*patch[0].fire.understory_et  +
            //(stratum[0].transpiration_sat_zone + stratum[0].transpiration_unsat_zone))/(patch[0].fire_defaults[0][0].ndays_average + 1);
            enqueue_FIFO_Queue(&patch[0].fire.Q_understory_et,(stratum[0].transpiration_sat_zone + stratum[0].transpiration_unsat_zone));

		//	+ stratum[0].evaporation))/ // MCK: look at et without the evaporation part

            //patch[0].fire.understory_pet = (patch[0].fire_defaults[0][0].ndays_average*patch[0].fire.understory_pet  +
            //(stratum[0].PET))/
            //(patch[0].fire_defaults[0][0].ndays_average + 1);
            enqueue_FIFO_Queue(&patch[0].fire.Q_understory_pet,stratum[0].PET);

//		printf("Debugging understory et and pet,stratum PET, stratum transpiration sat and unsat, stratum evaporation: %lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",patch[0].fire.understory_et,patch[0].fire.understory_pet, stratum[0].PET,stratum[0].transpiration_sat_zone, stratum[0].transpiration_unsat_zone,stratum[0].evaporation);

		}
	}

	if ( command_line[0].verbose_flag == -5 ){
	printf("\n          STRATUM DAILY END Kdir=%lf Kdif=%lf Kupdir=%lf Kupdif=%lf Lstar=%lf Lstarpch=%lf Rnet_used=%lf Kstar_can=%lf rnetevap=%lf dayl=%lf Tcan=%lf \n          ??? topt=%lf tcoef=%lf tmax=%lf",
		   Kdown_direct/86.4,
		   Kdown_diffuse/86.4,
		   patch[0].Kup_direct_final/86.4,
		   patch[0].Kup_diffuse_final/86.4,
		   stratum[0].Lstar/86.4,
		   patch[0].Lstar_canopy/86.4,
		   Rnet_used/86.4,
		   patch[0].Kstar_canopy_final/86.4,
		   rnet_evap,
		   zone[0].metv.dayl,
		   patch[0].T_canopy_final,
		   stratum[0].defaults[0][0].epc.topt,
		   stratum[0].defaults[0][0].epc.tcoef,
		   stratum[0].defaults[0][0].epc.tmax);
	}

	/*------------------------------------------------------------------------*/
	/*	If spinup option is set, update the shadow stratum until the targets  */
	/*	have been met                                                       	*/
	/*------------------------------------------------------------------------*/
	if(command_line[0].vegspinup_flag > 0){
    update_shadow_strata(world, patch, stratum, shadow_strata, command_line, current_date);
  }

   /***********************************************************************/
   /* if beetle attack first update the snag pool after the delay years */
   /**********************************************************************/

      	if (patch[0].base_stations != NULL && stratum[0].snag_sequence.seq !=NULL) {
		inx = patch[0].base_stations[0][0].dated_input[0].beetle_attack.inx; //this inx is the last index of climate event
       // num_records = patch[0].base_stations[0][0].dated_input[0].beetle_attack.num_records;

       // clim_event = patch[0].base_stations[0][0].dated_input[0].beetle_attack.seq[0];

        if (inx >-999) {
        /* update the snag pool from burning before you update it, only happens if firespread command line is on and there is snag pool*/
        /* due to the snagc is not the total snag, when snagc =0 the delay_snag is >0, so the condition snagc>0, make the fire can not burn the delay_snag*/
        /* so it is better to use the Snagc >= 0 condition for next step calculations, and just to make sure there is no negative values in flux out */
        /* if you put snagc and delay_snag combine with redneedle and delay_redneedle together condition may cause not burn once redneedle is zero */
        /* if this solution not working, maybe build two if conditions sepeartedly for snag and delay snag then using conditions and snag >=0 delay_snag>=0 */
            if (command_line[0].firespread_flag == 1 && (stratum[0].cs.snagc + stratum[0].cs.delay_snagc) >= ZERO && (stratum[0].cs.redneedlec + stratum[0].cs.delay_redneedlec) >= ZERO && patch[0].overstory_burn > ZERO){

                //burn the carbon pool
               //printf("\n the overstory burn is %lf \n", patch[0].overstory_burn);
                overstory_burn = patch[0].overstory_burn;
                redneedlec_burn = stratum[0].cs.redneedlec * overstory_burn;
                delay_redneedlec_burn = stratum[0].cs.delay_redneedlec * overstory_burn;

                snagc_burn = stratum[0].cs.snagc * overstory_burn;
                delay_snagc_burn = stratum[0].cs.delay_snagc * overstory_burn;

                //burn the redneedle and snag carbon

                stratum[0].cs.redneedlec = max(0,(stratum[0].cs.redneedlec - redneedlec_burn));
                stratum[0].cs.delay_redneedlec = max(0,(stratum[0].cs.delay_redneedlec - delay_redneedlec_burn)); //fix bug NREN 20190923

                stratum[0].cs.snagc = max(0, (stratum[0].cs.snagc - snagc_burn));
                stratum[0].cs.delay_snagc = max(0, (stratum[0].cs.delay_snagc - delay_snagc_burn));
                //burn the nitrogen pool
                redneedlen_burn = stratum[0].ns.redneedlen * overstory_burn;
                delay_redneedlen_burn = stratum[0].ns.delay_redneedlen * overstory_burn;

                snagn_burn = stratum[0].ns.snagn * overstory_burn;
                delay_snagn_burn = stratum[0].ns.delay_snagn * overstory_burn;

                //burn the redneedle and snag nitrogen

                stratum[0].ns.redneedlen = max(0, (stratum[0].ns.redneedlen - redneedlen_burn));
                stratum[0].ns.delay_redneedlen = max(0, (stratum[0].ns.delay_redneedlen - delay_redneedlen_burn));

                stratum[0].ns.snagn = max(0, (stratum[0].ns.snagn - snagn_burn));
                stratum[0].ns.delay_snagn = max(0, (stratum[0].ns.delay_snagn - delay_snagn_burn));


            }

        clim_event = patch[0].base_stations[0][0].dated_input[0].beetle_attack.seq[0];
            /*****************************************************************/
            /* if the beetle attack, update the dead leaf( red needle pool) */
            /****************************************************************/
           leaf_year_delay = world[0].defaults[0].beetle[0].leaf_year_delay;
           year_delay = world[0].defaults[0].beetle[0].year_delay;  // this is the snag delay time

        // this is to make sure if you have multiple attack, the current date is only updating the red needle and snag that corresponding the correct output break sequency
        // so here instead of hard coded as 2.5 years or 180 days, we need to figure out a better way that
        if (julday(clim_event.edate) + 0.8*leaf_year_delay*365 < julday(current_date)) { // what does this mean? to make sure the attack happens half year before the attack

       for (inx=0; inx < world[0].defaults[0].beetle[0].num_snag_sequence; inx=inx+24){ // here 300 is hard coded here, should figure out some other method.

        clim_event1 = stratum[0].redneedle_sequence.seq[inx];
        clim_event2 = stratum[0].snag_sequence.seq[inx];

        if (clim_event2.Cvalue > ZERO && clim_event2.Cvalue<100 && command_line[0].firespread_flag ==1 &&  patch[0].overstory_burn > ZERO)
            {
                //printf("\n the overstory burn is %lf \n", patch[0].overstory_burn);
                overstory_burn = patch[0].overstory_burn;
                redneedlec_burn = clim_event1.Cvalue * overstory_burn;
                snagc_burn = clim_event2.Cvalue * overstory_burn;

                redneedlen_burn = clim_event1.Nvalue * overstory_burn;
                snagn_burn = clim_event2.Nvalue * overstory_burn;

                stratum[0].redneedle_sequence.seq[inx].Cvalue -= redneedlec_burn;
                stratum[0].redneedle_sequence.seq[inx].Nvalue -= redneedlen_burn;
                stratum[0].snag_sequence.seq[inx].Cvalue -= snagc_burn;
                stratum[0].snag_sequence.seq[inx].Nvalue -= snagn_burn;
            }

        // if Cvalue >=0 then will updating understory too (OR updating pools that is zero),but due to the understory is zero snag pool
        // so it is not an issue
        if ((clim_event2.edate.year != 0) &&(clim_event2.Cvalue > ZERO)&&(clim_event2.Cvalue<100) && (julday(clim_event2.edate) + (int)(leaf_year_delay*365.256) == julday(current_date)))
                  {

                   if (inx ==0) {
                      if(stratum[0].ID == 7788) {
                   printf(" \n updating the leaf pool1 the inx is %d \n", inx); } // end if at stratum.ID
                    stratum[0].cs.redneedlec = stratum[0].redneedle_sequence.seq[inx].Cvalue;// there is getting the wrong index of value or call stratum first they call the daily so no value?
                    stratum[0].ns.redneedlen = stratum[0].redneedle_sequence.seq[inx].Nvalue;
                    stratum[0].cs.delay_redneedlec -= stratum[0].redneedle_sequence.seq[inx].Cvalue; //move out the delayed snag pool to ready to decay pool (decaying pool )
                    stratum[0].ns.delay_redneedlen -= stratum[0].redneedle_sequence.seq[inx].Nvalue; //NREN 20180728 is it possible this value <0?


                    }
                    else if (inx >0) {
                        if (stratum[0].ID == 7788) {
                    printf(" updating the leaf pool 2, the inx is %d \n", inx);}
                    stratum[0].cs.redneedlec += stratum[0].redneedle_sequence.seq[inx].Cvalue;
                    stratum[0].ns.redneedlen += stratum[0].redneedle_sequence.seq[inx].Nvalue;
                    stratum[0].cs.delay_redneedlec -= stratum[0].redneedle_sequence.seq[inx].Cvalue; //move out the delayed snag pool to ready to decay pool (decaying pool )
                    stratum[0].ns.delay_redneedlen -= stratum[0].redneedle_sequence.seq[inx].Nvalue; //NREN 20180728 is it possible this value <0?

                    }
                        } // if leaf_year_delay

        //how to make sure the model only update the overstory, not updating the understory since this snag is coresponding with stratum, if not over, the c&N is zero so still fine
        if ((clim_event2.edate.year != 0 ) && (clim_event2.Cvalue > ZERO)&&(clim_event2.Cvalue<100) &&(julday(clim_event2.edate) + (int)(year_delay*365.256) == julday(current_date)))
                   {

              if (inx ==0) {
                  if (stratum[0].ID == 7788) {
                    printf(" updating the snag pool the inx is %d\n", inx); }
                    stratum[0].cs.snagc = stratum[0].snag_sequence.seq[inx].Cvalue;// there is getting the wrong index of value or call stratum first they call the daily so no value?
                    stratum[0].ns.snagn = stratum[0].snag_sequence.seq[inx].Nvalue; // here need to change
                    stratum[0].cs.delay_snagc -= stratum[0].snag_sequence.seq[inx].Cvalue; //move the delayed(Wating ) snag pool to ready to decay (decaying pool)
                    stratum[0].ns.delay_snagn -= stratum[0].snag_sequence.seq[inx].Nvalue;  //NREN 20180728

                    }
              else if (inx >0) {
                    if (stratum[0].ID == 7788) {
                printf("updating the snag pool2 the inx is %d\n", inx); }//why??                                                       }
                stratum[0].cs.snagc += stratum[0].snag_sequence.seq[inx].Cvalue;
                stratum[0].ns.snagn += stratum[0].snag_sequence.seq[inx].Nvalue;
                stratum[0].cs.delay_snagc -= stratum[0].snag_sequence.seq[inx].Cvalue; //move the delayed(Wating ) snag pool to ready to decay (decaying pool)
                stratum[0].ns.delay_snagn -= stratum[0].snag_sequence.seq[inx].Nvalue;  //NREN 20180728

                             /*first move the snag _sequences to the dacay pool */
                }
                break;
                }

            } // enf for
                }

        } //end if (inx>-999)

        } //end 1959

        if (stratum[0].cs.snagc > ZERO ) {
                 compute_snag_decay(
                       &(stratum[0].cs),
                       &(stratum[0].ns),
                       &(patch[0].cdf), // here can change the patch[0].cdf to patch[0]
                       &(patch[0].ndf),
                       world[0].defaults[0].beetle[0].year_delay,
                       world[0].defaults[0].beetle[0].half_life,
                       current_date);
        }

        if (stratum[0].cs.redneedlec > ZERO) {

                    compute_redneedle_decay(
                    stratum[0].defaults[0][0].epc,
                    &(stratum[0].cs),
                    &(stratum[0].ns),
                    &(patch[0].litter_cs),
                    &(patch[0].litter_ns),
                    &(patch[0].cdf),
                    &(patch[0].ndf),
                    world[0].defaults[0].beetle[0].leaf_year_delay,
                    world[0].defaults[0].beetle[0].leaf_half_life,
                    current_date );
        }
        //add the dead root decay casued by beetle attack NREN 20190911
        if (stratum[0].cs.dead_rootc_beetle >0 && world[0].defaults[0].beetle[0].root_alive == 2){

                    update_deadroot_decay(
                    stratum[0].defaults[0][0].epc,
                    &(stratum[0].cs),
                    &(stratum[0].ns),
                    &(patch[0].litter_cs),
                    &(patch[0].litter_ns),
                    &(patch[0].cdf),
                    &(patch[0].ndf),
                    world[0].defaults[0].beetle[0].deadroot_half_life);

        }



	/*--------------------------------------------------------------*/
	/*      update accumlator variables                             */
	/*--------------------------------------------------------------*/
	if((command_line[0].output_flags.monthly == 1)&&(command_line[0].c != NULL)){
		stratum[0].acc_month.psn += stratum[0].cdf.psn_to_cpool - stratum[0].cdf.total_mr;
		stratum[0].acc_month.lwp += stratum[0].epv.psi;
		stratum[0].acc_month.leafc += stratum[0].cs.leafc;
		stratum[0].acc_month.rootc += stratum[0].cs.frootc+stratum[0].cs.live_crootc+stratum[0].cs.dead_crootc;
		stratum[0].acc_month.stemc += stratum[0].cs.live_stemc+stratum[0].cs.dead_stemc;
		stratum[0].acc_month.length += 1;
	}
	/* output fire effect */
	if ((command_line[0].output_flags.yearly == 1) && (command_line[0].c != NULL || command_line[0].f != NULL)){
        //06012022LML stratum[0].acc_year.psn += stratum[0].cdf.psn_to_cpool - stratum[0].cdf.total_mr;
		stratum[0].acc_year.lwp += stratum[0].epv.psi;
        if (close_enough(stratum[0].acc_year.minNSC, -999))
			stratum[0].acc_year.minNSC = stratum[0].cs.cpool;
		else
			stratum[0].acc_year.minNSC = min(stratum[0].cs.cpool, stratum[0].acc_year.minNSC);

#ifdef JMG_MORE_YEARLY_OUTPUT
        //12162022LML stratum[0].acc_year.LAIx += (stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc) * stratum[0].defaults[0][0].epc.proj_sla;
        stratum[0].acc_year.LAIx += stratum[0].cs.leafc * stratum[0].defaults[0][0].epc.proj_sla;
        stratum[0].acc_year.height += stratum[0].epv.height;
        stratum[0].acc_year.AGBc += ( // above ground biomass carbon weight (kgC/m2)
                stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc + // leafc
                stratum[0].cs.live_stemc + stratum[0].cs.dead_stemc + stratum[0].cs.livestemc_store + stratum[0].cs.deadstemc_store + stratum[0].cs.livestemc_transfer + stratum[0].cs.deadstemc_transfer + // stemc
                stratum[0].cs.cwdc + stratum[0].cs.cpool); // remaining biomass
        stratum[0].acc_year.BGBc += ( // below ground biomass carbon weight (kgC/m2), root carbon mainly
                stratum[0].cs.frootc + stratum[0].cs.live_crootc + stratum[0].cs.dead_crootc); // live croot + live froot + dead croot
        //stratum[0].acc_year.gpp += stratum[0].cs.net_psn - (stratum[0].cdf.total_mr + stratum[0].cdf.total_gr); // check on this
        //stratum[0].acc_year.gpp += stratum[0].cdf.psn_to_cpool;
        //stratum[0].acc_year.resp += stratum[0].cdf.total_mr + stratum[0].cdf.total_gr;
        //stratum[0].acc_year.resp += stratum[0].cdf.total_mr + stratum[0].cdf.total_gr;
        //stratum[0].acc_year.npp += stratum[0].cs.net_psn;
        //stratum[0].acc_year.npp += stratum[0].cdf.psn_to_cpool - stratum[0].cdf.total_mr - stratum[0].cdf.total_gr;
        stratum[0].acc_year.rootdepth += stratum[0].rootzone.depth;        
#endif

        stratum[0].acc_year.length += 1;
	}
	return;
} /*end canopy_stratum_daily_F.c*/
