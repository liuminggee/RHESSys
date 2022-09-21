/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_patch						*/
/*																*/
/*	output_yearly_patch - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_yearly_patch - outputs current contents of a patch.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_patch( int basinID, int hillID, int zoneID,								*/
/*					struct	patch_object	*patch,				*/
/*					struct	date	date,  						*/
/*					FILE 	*outfile)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs spatial structure according to commandline			*/
/*	specifications to specific files							*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	output_yearly_patch(
				int basinID, int hillID, int zoneID,
				struct	patch_object	*patch,
				struct	date	current_date,
				FILE *outfile)
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
#ifdef JMG_MORE_YEARLY_OUTPUT
    struct	canopy_strata_object *strata;
    strata = patch[0].canopy_strata[(patch[0].layers[0].strata[0])];

    double plantc, plantn, litrc, litrn, soilc, soiln, aleafc, aleafn, afrootc, afrootn, awoodc, awoodn;
    plantc = plantn = litrc= litrn = soilc = soiln = aleafc = aleafn = afrootc = afrootn = awoodc = awoodn = 0.0;

    litrc += patch[0].litter_cs.litr1c + patch[0].litter_cs.litr2c + patch[0].litter_cs.litr3c + patch[0].litter_cs.litr4c;
    litrn += patch[0].litter_ns.litr1n + patch[0].litter_ns.litr2n + patch[0].litter_ns.litr3n + patch[0].litter_ns.litr4n;
    soilc += patch[0].soil_cs.soil1c + patch[0].soil_cs.soil2c + patch[0].soil_cs.soil3c + patch[0].soil_cs.soil4c;
    soiln += patch[0].soil_ns.soil1n + patch[0].soil_ns.soil2n + patch[0].soil_ns.soil3n + patch[0].soil_ns.soil4n;

    aleafc += strata->cover_fraction * (strata->cs.leafc
        + strata->cs.leafc_store + strata->cs.leafc_transfer );

    aleafn += strata->cover_fraction * (strata->ns.leafn
        + strata->ns.leafn_store + strata->ns.leafn_transfer );

    afrootc += strata->cover_fraction
        * (strata->cs.frootc + strata->cs.frootc_store
        + strata->cs.frootc_transfer);

    afrootn += strata->cover_fraction
        * (strata->ns.frootn + strata->ns.frootn_store
        + strata->ns.frootn_transfer);
    awoodc += strata->cover_fraction * (strata->cs.live_crootc
        + strata->cs.live_stemc + strata->cs.dead_crootc
        + strata->cs.dead_stemc + strata->cs.livecrootc_store
        + strata->cs.livestemc_store + strata->cs.deadcrootc_store
        + strata->cs.deadstemc_store + strata->cs.livecrootc_transfer
        + strata->cs.livestemc_transfer + strata->cs.deadcrootc_transfer
        + strata->cs.deadstemc_transfer
        + strata->cs.cwdc + strata->cs.cpool);
    awoodn += strata->cover_fraction * (strata->ns.live_crootn
        + strata->ns.live_stemn + strata->ns.dead_crootn
        + strata->ns.dead_stemn + strata->ns.livecrootn_store
        + strata->ns.livestemn_store + strata->ns.deadcrootn_store
        + strata->ns.deadstemn_store + strata->ns.livecrootn_transfer
        + strata->ns.livestemn_transfer + strata->ns.deadcrootn_transfer
        + strata->ns.deadstemn_transfer
        + strata->ns.cwdn + strata->ns.npool + strata->ns.retransn);

    plantc += aleafc + afrootc + awoodc;
    plantn += aleafn + afrootn + awoodn;
#endif

	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/

	if (patch[0].acc_year.length > 0) 
		patch[0].acc_year.theta /= patch[0].acc_year.length;
#ifdef JMG_MORE_YEARLY_OUTPUT
        patch[0].acc_year.rz_storage /= patch[0].acc_year.length; // JMG09082022
        patch[0].acc_year.unsat_storage /= patch[0].acc_year.length; // JMG09082022
#endif
	if (patch[0].acc_year.recharge > ZERO)
		patch[0].acc_year.recharge_wyd /= patch[0].acc_year.recharge;

    fprintf(outfile,
#ifndef JMG_MORE_YEARLY_OUTPUT
            "%lf %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
#else
            "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
#endif
			current_date.year,
			basinID,
			hillID,
			zoneID,
			patch[0].ID,
#ifndef JMG_MORE_YEARLY_OUTPUT
			patch[0].acc_year.num_threshold,
			patch[0].acc_year.peaksweday,
			patch[0].acc_year.meltday,
			patch[0].acc_year.leach * 1000.0,
			patch[0].acc_year.denitrif * 1000.0,
			patch[0].acc_year.DOC_loss * 1000.0,
			patch[0].acc_year.DON_loss * 1000.0,
			patch[0].acc_year.psn * 1000.0,
			patch[0].acc_year.trans * 1000.0,
			patch[0].acc_year.et * 1000.0,
			patch[0].acc_year.lai,
			patch[0].acc_year.nitrif * 1000.0,
			patch[0].acc_year.mineralized * 1000.0,
			patch[0].acc_year.uptake * 1000.0,
			patch[0].acc_year.theta,
			patch[0].acc_year.sm_deficit*1000.0,
			patch[0].acc_year.snowpack * 1000.0,
			patch[0].acc_year.maxtrans * 1000.0,
			patch[0].acc_year.maxpet * 1000.0,
			patch[0].acc_year.streamflow * 1000.0,
			patch[0].acc_year.Qin_total * 1000.0,
			patch[0].acc_year.Qout_total * 1000.0,
			patch[0].acc_year.rec_wyd,
			patch[0].acc_year.rec_pet_wyd,
			patch[0].acc_year.ndays_sat, patch[0].acc_year.ndays_sat70, 
			patch[0].acc_year.midsm_wyd,
            patch[0].z,
            patch[0].acc_year.TPET*1000.0,
            patch[0].acc_year.PET*1000.0,
            patch[0].acc_year.PE*1000.0,
            patch[0].acc_year.pcp*1000.0, patch[0].acc_year.burn,
			patch[0].acc_year.snowin*1000.0, patch[0].acc_year.potential_recharge*1000.0,
			patch[0].acc_year.recharge*1000.0, patch[0].acc_year.potential_recharge_wyd,
            patch[0].acc_year.recharge_wyd
#else
            patch[0].acc_year.pcp*1000.0,
            patch[0].acc_year.streamflow * 1000.0, // streamflow (mm/yr)
            patch[0].acc_year.baseflow * 1000.0, // baseflow (mm/yr) // JMG09082022
            patch[0].acc_year.returnflow * 1000.0, // returnflow (mm/yr)
            patch[0].acc_year.rz_storage * 1000.0, // avg daily rootzone storage (mm) // JMG09082022
            patch[0].acc_year.unsat_storage * 1000.0, // avg daily unsat storage (mm) // JMG09082022
            patch[0].acc_year.gw_drainage * 1000.0, // avg daily gw storage (mm) // JMG09082022
            patch[0].acc_year.overland_flow * 1000.0, // JMG09122022
            patch[0].acc_year.et * 1000.0, // evapotranspiration
            patch[0].acc_year.ndays_sat, // number days at saturation (sat_def - unsat_stor - rz_stor <= 0)
            patch[0].acc_year.ndays_sat70, // number of days where rz moisture is > 70% saturation (rootzone.S > 0.7)
            plantc,
            plantn,
            litrc,
            litrn,
            soilc,
            soiln,
            patch[0].acc_year.lai, // JMG09082022 shows max LAI on the year (not the average)
            patch[0].acc_year.psn * 1000.0,
            patch[0].acc_year.uptake * 1000.0,
            patch[0].acc_year.leach * 1000.0,
            patch[0].acc_year.DON_loss * 1000.0,
            patch[0].acc_year.denitrif * 1000.0,
            patch[0].acc_year.nitrif * 1000.0,
            patch[0].acc_year.mineralized * 1000.0
#endif
            );


	/*--------------------------------------------------------------*/
	/*      reset accumulator variables                             */
	/*--------------------------------------------------------------*/
	patch[0].acc_year.num_threshold = 0;
	patch[0].acc_year.peaksweday = 0;
	patch[0].acc_year.meltday = 0;
	patch[0].acc_year.leach = 0.0;
	patch[0].acc_year.denitrif = 0.0;
	patch[0].acc_year.nitrif = 0.0;
	patch[0].acc_year.DOC_loss = 0.0;
	patch[0].acc_year.DON_loss = 0.0;
	patch[0].acc_year.psn = 0.0;
	patch[0].acc_year.et = 0.0;
	patch[0].acc_year.trans = 0.0;
    patch[0].acc_year.TPET = 0.0;
    patch[0].acc_year.PET = 0.0;
    patch[0].acc_year.PE = 0.0;
	patch[0].acc_year.lai = 0.0;
	patch[0].acc_year.mineralized = 0.0;
	patch[0].acc_year.uptake = 0.0;
	patch[0].acc_year.theta = 0.0;
	patch[0].acc_year.sm_deficit = 0.0;
	patch[0].acc_year.snowpack = 0.0;
	patch[0].acc_year.length = 0;
	patch[0].acc_year.wyd = 0;
	patch[0].acc_year.rec_wyd = 0;
	patch[0].acc_year.day7trans = 0;
	patch[0].acc_year.maxtrans = 0;
	patch[0].acc_year.midsm_wyd = 0;
	patch[0].acc_year.ndays_sat = 0;
	patch[0].acc_year.ndays_sat70 = 0;
	patch[0].acc_year.Qin_total = 0.0;
	patch[0].acc_year.Qout_total = 0.0;
	patch[0].acc_year.pcp = 0.0;
	patch[0].acc_year.streamflow = 0.0;    
#ifdef JMG_MORE_YEARLY_OUTPUT
    patch[0].acc_year.baseflow = 0.0; // JMG09082022
    patch[0].acc_year.returnflow = 0.0; // JMG09082022
    patch[0].acc_year.rz_storage = 0.0; // JMG09082022
    patch[0].acc_year.unsat_storage = 0.0; // JMG09082022
    patch[0].acc_year.gw_drainage = 0.0; // JMG09082022
    patch[0].acc_year.overland_flow = 0.0; // JMG09122022
#endif
	patch[0].acc_year.burn = 0.0;
	patch[0].acc_year.snowin = 0.0;
	patch[0].acc_year.potential_recharge = 0.0;
	patch[0].acc_year.potential_recharge_wyd = 0.0;
	patch[0].acc_year.recharge = 0.0;
	patch[0].acc_year.recharge_wyd = 0.0;

	return;

} /*end output_csv_yearly_patch*/

