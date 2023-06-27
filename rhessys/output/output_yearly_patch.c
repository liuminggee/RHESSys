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
                FILE *outfile
#ifdef JMG_TRACKING
                ,struct simtime *simtime
#endif
        )
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/

	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/

	if (patch[0].acc_year.length > 0) 
		patch[0].acc_year.theta /= patch[0].acc_year.length;

	if (patch[0].acc_year.recharge > ZERO)
		patch[0].acc_year.recharge_wyd /= patch[0].acc_year.recharge;

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic[] = "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#else
    char out_basic[] = "%d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#endif

#ifdef JMG_TRACKING
    char out_format[1000] = "%d %d %d ";
    strcat(out_format,out_basic);
#else
    char out_format[1000] = "";
    strcat(out_format, out_basic);
#endif

    fprintf(outfile, out_format,

#ifdef JMG_TRACKING
            simtime->sday,
            simtime->smth,
            simtime->syr,
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
            patch[0].acc_year.ndays_sat,
            patch[0].acc_year.ndays_sat70,
			patch[0].acc_year.midsm_wyd,
            patch[0].area,
            patch[0].acc_year.TPET*1000.0,
            patch[0].acc_year.PET*1000.0,
            patch[0].acc_year.PE*1000.0,
            patch[0].acc_year.pcp*1000.0,
            patch[0].acc_year.burn,
            patch[0].acc_year.snowin*1000.0,
            patch[0].acc_year.potential_recharge*1000.0,
            patch[0].acc_year.recharge*1000.0,
            patch[0].acc_year.potential_recharge_wyd,
            patch[0].acc_year.recharge_wyd,
            patch[0].acc_year.p_gw_drainage*1000,
            (patch[0].acc_year.pcp-patch[0].acc_year.et-patch[0].acc_year.streamflow-patch[0].acc_year.p_gw_drainage)*1000.0
#else
            patch[0].acc_year.pcp*1000.0,
            patch[0].acc_year.streamflow * 1000.0, // streamflow (mm/yr)
            patch[0].acc_year.baseflow * 1000.0, // baseflow (mm/yr) // JMG09082022
            patch[0].acc_year.returnflow * 1000.0, // returnflow (mm/yr)
            patch[0].acc_year.rz_storage * 1000.0 / patch[0].acc_year.length, // avg daily rootzone storage (mm) // JMG09082022
            patch[0].acc_year.unsat_storage * 1000.0 / patch[0].acc_year.length, // avg daily unsat storage (mm) // JMG09082022
            patch[0].acc_year.gw_drainage * 1000.0, // avg daily gw storage (mm) // JMG09082022
            patch[0].acc_year.overland_flow * 1000.0, // JMG09122022
            patch[0].acc_year.et * 1000.0, // evapotranspiration
            patch[0].acc_year.ndays_sat, // number days at saturation (sat_def - unsat_stor - rz_stor <= 0)
            patch[0].acc_year.ndays_sat70, // number of days where rz moisture is > 70% saturation (rootzone.S > 0.7)
            //plantc,
            patch[0].acc_year.plantc / patch[0].acc_year.length,
            //plantn,
            patch[0].acc_year.plantn / patch[0].acc_year.length,
            //litrc,
            patch[0].acc_year.litrc / patch[0].acc_year.length,
            //litrn,
            patch[0].acc_year.litrn / patch[0].acc_year.length,
            //soilc,
            patch[0].acc_year.soilc / patch[0].acc_year.length,
            //soiln,
            patch[0].acc_year.soiln / patch[0].acc_year.length,
            patch[0].acc_year.lai, // JMG09082022 shows max LAI on the year (not the average)
            patch[0].acc_year.psn * 1000.0,
            patch[0].acc_year.uptake * 1000.0,
            patch[0].acc_year.leach * 1000.0,
            patch[0].acc_year.DON_loss * 1000.0,
            patch[0].acc_year.denitrif * 1000.0,
            patch[0].acc_year.nitrif * 1000.0,
            patch[0].acc_year.mineralized * 1000.0,
            patch[0].acc_year.sat_deficit * 1000.0,
            patch[0].acc_year.sat_deficit_z * 1000.0
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
    patch[0].acc_year.p_gw_drainage = 0.0;
#ifdef JMG_MORE_YEARLY_OUTPUT    
    patch[0].acc_year.baseflow = 0.0; // JMG09082022
    patch[0].acc_year.returnflow = 0.0; // JMG09082022
    patch[0].acc_year.rz_storage = 0.0; // JMG09082022
    patch[0].acc_year.unsat_storage = 0.0; // JMG09082022
    patch[0].acc_year.gw_drainage = 0.0; // JMG09082022
    patch[0].acc_year.overland_flow = 0.0; // JMG09122022
    patch[0].acc_year.pcp = 0.0;
    patch[0].acc_year.pch_base_flow = 0.0; // JMG09122022

    patch[0].acc_year.n_deposition = 0.0; // JMG10112022
    patch[0].acc_year.plantc = 0.0; // JMG10112022
    patch[0].acc_year.plantn = 0.0; // JMG10112022
    patch[0].acc_year.AGBc = 0.0; // JMG10112022
    patch[0].acc_year.BGBc = 0.0; // JMG10112022
    patch[0].acc_year.soilc = 0.0; // JMG10112022
    patch[0].acc_year.soiln = 0.0; // JMG10112022
    patch[0].acc_year.litrc = 0.0; // JMG10112022
    patch[0].acc_year.litrn = 0.0; // JMG10112022
#endif
	patch[0].acc_year.burn = 0.0;
	patch[0].acc_year.snowin = 0.0;
	patch[0].acc_year.potential_recharge = 0.0;
	patch[0].acc_year.potential_recharge_wyd = 0.0;
	patch[0].acc_year.recharge = 0.0;
	patch[0].acc_year.recharge_wyd = 0.0;

	return;

} /*end output_csv_yearly_patch*/

