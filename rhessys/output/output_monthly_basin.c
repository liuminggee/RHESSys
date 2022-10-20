/*--------------------------------------------------------------*/
/* 																*/
/*					output_monthly_basin						*/
/*																*/
/*	output_monthly_basin - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_monthly_basin - outputs current contents of a basin.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_monthly_basin(										*/
/*					struct	basin_object	*basin,				*/
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

void	output_monthly_basin(
							 struct	basin_object	*basin,
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
	if (basin[0].acc_month.length == 0) basin[0].acc_month.length = 1;

	//if (basin->route_list->num_patches > 0)
	//	basin[0].acc_month.length /= (basin->route_list->num_patches);
  int patchCount = 0;
  for( int i = 0; i < basin[0].num_hillslopes; i++ ) {
    struct hillslope_object *hillslope = basin[0].hillslopes[i];
    if (hillslope->route_list != NULL) { //N REN 2019/04/07 2019/03/31 solve no lateral flow segment default bug
    patchCount += hillslope->route_list->num_patches;
    }
  }
  if( patchCount == 0 ) patchCount = 1;
  basin[0].acc_month.length /= patchCount;

#ifdef LIU_TRACKING_BASIN_LITTERC
  struct cdayflux_patch_struct *cdf = &basin[0].acc_month.cdf;
  double leafc_to_litrc = cdf->leafc_to_litr1c + cdf->cwdc_to_litr2c + cdf->cwdc_to_litr3c + cdf->cwdc_to_litr4c;
  double frootc_to_litrc = cdf->frootc_to_litr1c + cdf->frootc_to_litr2c + cdf->frootc_to_litr3c + cdf->frootc_to_litr4c;
  double cwdc_to_litrc = cdf->cwdc_to_litr2c + cdf->cwdc_to_litr3c + cdf->cwdc_to_litr4c;
  double stemc_to_litrc = cdf->stemc_to_litr1c;
  double mort_to_litrc = cdf->m_leafc_to_litr1c + cdf->m_leafc_to_litr2c + cdf->m_leafc_to_litr3c + cdf->m_leafc_to_litr4c
                         + cdf->m_frootc_to_litr1c + cdf->m_frootc_to_litr2c + cdf->m_frootc_to_litr3c + cdf->m_frootc_to_litr4c
                         + cdf->m_leafc_store_to_litr1c + cdf->m_frootc_store_to_litr1c + cdf->m_livestemc_store_to_litr1c + cdf->m_deadstemc_store_to_litr1c + cdf->m_livecrootc_store_to_litr1c + cdf->m_deadcrootc_store_to_litr1c
                         + cdf->m_leafc_transfer_to_litr1c + cdf->m_frootc_transfer_to_litr1c + cdf->m_livestemc_transfer_to_litr1c + cdf->m_deadstemc_transfer_to_litr1c + cdf->m_livecrootc_transfer_to_litr1c + cdf->m_deadcrootc_transfer_to_litr1c
                         + cdf->m_gresp_store_to_litr1c + cdf->m_gresp_transfer_to_litr1c;
  double do_litrc_loss = cdf->do_litr1c_loss + cdf->do_litr2c_loss + cdf->do_litr3c_loss + cdf->do_litr4c_loss;
  double m_litrc_to_atmos = cdf->m_litr1c_to_atmos + cdf->m_litr2c_to_atmos + cdf->m_litr3c_to_atmos + cdf->m_litr4c_to_atmos;
  double litrc_to_atmos = cdf->litterc_to_atmos;
  double litrc_to_soilc = cdf->litterc_to_soilc;
#endif

#ifdef LIU_TRACKING_BASIN_LITTERC
    char out_basic[] = "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#else
    char out_basic[] = "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#endif

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
		basin[0].ID,
		basin[0].acc_month.streamflow * 1000.0,
		basin[0].acc_month.stream_NO3 * 1000.0,
		basin[0].acc_month.denitrif * 1000.0,
		basin[0].acc_month.DOC_loss * 1000.0,
		basin[0].acc_month.DON_loss * 1000.0,
		basin[0].acc_month.et * 1000.0,
		basin[0].acc_month.psn * 1000.0,
		basin[0].acc_month.lai/basin[0].acc_month.length ,
		basin[0].acc_month.nitrif * 1000.0,
		basin[0].acc_month.mineralized * 1000.0,
		basin[0].acc_month.uptake * 1000.0
#ifdef LIU_TRACKING_BASIN_LITTERC
        ,
        leafc_to_litrc,
        frootc_to_litrc,
        cwdc_to_litrc,
        stemc_to_litrc,
        mort_to_litrc,
        do_litrc_loss,
        m_litrc_to_atmos,
        litrc_to_atmos,
        litrc_to_soilc
#endif
		);
	if (check <= 0) {
		fprintf(stdout,
			"\nWARNING: output error has occured in output_monthly_basin");
	}
	/*--------------------------------------------------------------*/
	/*	reset accumulator variables				*/
	/*--------------------------------------------------------------*/
	basin[0].acc_month.streamflow = 0.0;
	basin[0].acc_month.stream_NO3 = 0.0;
	basin[0].acc_month.et = 0.0;
	basin[0].acc_month.psn = 0.0;
	basin[0].acc_month.length = 0;
	basin[0].acc_month.DOC_loss = 0.0;
	basin[0].acc_month.DON_loss = 0.0;
	basin[0].acc_month.denitrif= 0.0;
	basin[0].acc_month.lai = 0.0;
	basin[0].acc_month.nitrif = 0.0;
	basin[0].acc_month.mineralized = 0.0;
	basin[0].acc_month.uptake = 0.0;
#ifdef LIU_TRACKING_BASIN_LITTERC
    memset(&basin[0].acc_month.cdf, 0, sizeof(struct cdayflux_patch_struct));
#endif
	return;
} /*end output_monthly_basin*/
