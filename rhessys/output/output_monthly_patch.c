/*--------------------------------------------------------------*/
/* 																*/
/*					output_monthly_patch						*/
/*																*/
/*	output_monthly_patch - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_monthly_patch - outputs current contents of a patch.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_monthly_patch(										*/
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

void	output_monthly_patch(
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
	int check;
	if (patch[0].acc_month.length == 0) patch[0].acc_month.length = 1;

	if (patch[0].acc_month.leach > 0.0)
		patch[0].acc_month.leach = log(patch[0].acc_month.leach*1000.0*1000.0);

    char out_basic[] = "%d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %8.3f %f %f %f %f %f %f \n";

#ifdef JMG_TRACKING
    char out_format[1000] = "%d %d %d ";
    strcat(out_format,out_basic);
#else
    char out_format[1000] = "";
    strcat(out_format, out_basic);
#endif
		
    check = fprintf(outfile, out_format,

#ifdef JMG_TRACKING
        simtime->sday,
        simtime->smth,
        simtime->syr,
#endif

		current_date.month,
		current_date.year,
		basinID,
		hillID,
		zoneID,
		patch[0].ID,
		patch[0].acc_month.leach * 1000.0,
		patch[0].acc_month.denitrif * 1000.0,
		patch[0].acc_month.sm_deficit  / patch[0].acc_month.length * 1000.0,
		patch[0].acc_month.et * 1000.0,
		patch[0].acc_month.psn * 1000.0,
		patch[0].acc_month.DOC_loss * 1000.0,
		patch[0].acc_month.DON_loss * 1000.0,
		patch[0].acc_month.lai,
		patch[0].acc_month.nitrif * 1000.0,
		patch[0].acc_month.mineralized * 1000.0,
		patch[0].acc_month.uptake * 1000.0,
		patch[0].acc_month.theta / patch[0].acc_month.length * 100.0,
		patch[0].acc_month.snowpack * 1000.0 ,
		patch[0].area,
		patch[0].soil_ns.nitrate+patch[0].surface_NO3,
		patch[0].soil_ns.sminn

		);

	if (check <= 0) {
		fprintf(stdout,
			"\nWARNING: output error has occured in output_monthly_patch");
	}
	/*--------------------------------------------------------------*/
	/*	reset accumulator variables				*/
	/*--------------------------------------------------------------*/
	patch[0].acc_month.sm_deficit = 0.0;
	patch[0].acc_month.et = 0.0;
	patch[0].acc_month.psn = 0.0;
	patch[0].acc_month.lai = 0.0;
	patch[0].acc_month.length = 0;
	patch[0].acc_month.leach= 0.0;
	patch[0].acc_month.DOC_loss = 0.0;
	patch[0].acc_month.DON_loss = 0.0;
	patch[0].acc_month.denitrif= 0.0;
	patch[0].acc_month.nitrif= 0.0;
	patch[0].acc_month.mineralized = 0.0;
	patch[0].acc_month.uptake = 0.0;
	patch[0].acc_month.theta = 0.0;
	patch[0].acc_month.snowpack = 0.0;
	return;
} /*end output_monthly_patch*/
