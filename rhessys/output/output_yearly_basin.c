/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_basin						*/
/*																*/
/*	output_yearly_basin - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_yearly_basin - outputs current contents of a basin.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_basin(										*/
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

void	output_yearly_basin(
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


//	if (basin->route_list->num_patches > 0)
//		basin[0].acc_year.length /= basin->route_list->num_patches;
  int patchCount = 0;
  for( int i = 0; i < basin[0].num_hillslopes; i++ ) {
    struct hillslope_object *hillslope = basin[0].hillslopes[i];
    if (hillslope->route_list != NULL) { //N REN 2019/04/07 2019/03/31 solve no lateral flow segment default bug
    patchCount += hillslope->route_list->num_patches;
    }
  }
  if( patchCount == 0 ) patchCount = 1;
  basin[0].acc_year.length /= patchCount;

    if (basin[0].acc_year.length == 0) basin[0].acc_year.length = 1;

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic[] = "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#else
    char out_basic[] = "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf\n";
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

        current_date.year,
        basin[0].ID,
#ifndef JMG_MORE_YEARLY_OUTPUT
        basin[0].acc_year.streamflow * 1000.0,
        basin[0].acc_year.stream_NO3 * 1000.0,
        basin[0].acc_year.denitrif * 1000.0, // gN/m2
        basin[0].acc_year.DOC_loss * 1000.0,
        basin[0].acc_year.DON_loss * 1000.0,
        basin[0].acc_year.et * 1000.0,
        basin[0].acc_year.psn * 1000.0,
        basin[0].acc_year.lai/ basin[0].acc_year.length,
        basin[0].acc_year.nitrif * 1000.0,
        basin[0].acc_year.mineralized * 1000.0,
        basin[0].acc_year.uptake * 1000.0,
        basin[0].acc_year.num_threshold,
        basin[0].acc_year.TPET * 1000.0,
        basin[0].acc_year.PET * 1000.0,
        basin[0].acc_year.PE * 1000.0
#else
            basin[0].acc_year.pcp * 1000.0,
            basin[0].acc_year.et * 1000.0,
            basin[0].acc_year.TPET * 1000.0,
            basin[0].acc_year.PET * 1000.0,
            basin[0].acc_year.PE * 1000.0,
            basin[0].acc_year.streamflow * 1000.0,
            basin[0].acc_year.baseflow * 1000.0,
            basin[0].acc_year.gw_drainage * 1000.0,
            basin[0].acc_year.rz_storage * 1000.0 / basin[0].acc_year.length,
            basin[0].acc_year.unsat_storage * 1000.0 / basin[0].acc_year.length,
            basin[0].acc_year.hill_gw_storage * 1000.0 / basin[0].acc_year.length,
            basin[0].acc_year.n_deposition * 1000,
            basin[0].acc_year.denitrif * 1000.0 / basin[0].acc_year.length,
            basin[0].acc_year.soilc / basin[0].acc_year.length,
            basin[0].acc_year.soiln / basin[0].acc_year.length,
            basin[0].acc_year.litrc / basin[0].acc_year.length,
            basin[0].acc_year.litrn / basin[0].acc_year.length,
            basin[0].acc_year.plantc / basin[0].acc_year.length,
            basin[0].acc_year.plantn / basin[0].acc_year.length,
            basin[0].acc_year.AGBc / basin[0].acc_year.length,
            basin[0].acc_year.BGBc / basin[0].acc_year.length,
            basin[0].acc_year.lai / basin[0].acc_year.length,
            basin[0].acc_year.sat_deficit * 1000.0/ basin[0].acc_year.length,
            basin[0].acc_year.sat_deficit_z * 1000.0 / basin[0].acc_year.length
#endif
        );
    if (check <= 0) {
        fprintf(stdout,
            "\nWARNING: output error has occured in output_yearly_basin");
    }
    /*--------------------------------------------------------------*/
    /*	reset accumulator variables				*/
    /*--------------------------------------------------------------*/
    basin[0].acc_year.streamflow = 0.0;
    basin[0].acc_year.stream_NO3 = 0.0;
    basin[0].acc_year.et = 0.0;
    basin[0].acc_year.TPET = 0.0;
    basin[0].acc_year.PET = 0.0;
    basin[0].acc_year.PE = 0.0;
    basin[0].acc_year.psn = 0.0;
    basin[0].acc_year.lai = 0.0;
    basin[0].acc_year.length = 0;
    basin[0].acc_year.DOC_loss = 0.0;
    basin[0].acc_year.DON_loss = 0.0;
    basin[0].acc_year.denitrif= 0.0;
    basin[0].acc_year.nitrif= 0.0;
    basin[0].acc_year.num_threshold = 0;
    basin[0].acc_year.mineralized = 0.0;
    basin[0].acc_year.uptake = 0.0;

#ifdef JMG_MORE_YEARLY_OUTPUT
    basin[0].acc_year.pcp = 0.0; // JMG10112022
    basin[0].acc_year.baseflow = 0.0; // JMG10112022
    basin[0].acc_year.hill_base_flow = 0.0; // JMG10112022
    basin[0].acc_year.gw_drainage = 0.0; // JMG10112022
    basin[0].acc_year.rz_storage = 0.0; // JMG10112022
    basin[0].acc_year.unsat_storage = 0.0; // JMG10112022
    basin[0].acc_year.hill_gw_storage = 0.0; // JMG10112022

    basin[0].acc_year.n_deposition = 0.0; // JMG10112022
    basin[0].acc_year.soilc = 0.0; // JMG10112022
    basin[0].acc_year.soiln = 0.0; // JMG10112022
    basin[0].acc_year.litrc = 0.0; // JMG10112022
    basin[0].acc_year.litrn = 0.0; // JMG10112022
    basin[0].acc_year.plantc = 0.0; // JMG10112022
    basin[0].acc_year.plantn = 0.0; // JMG10112022
    basin[0].acc_year.AGBc = 0.0; // JMG10112022
    basin[0].acc_year.BGBc = 0.0; // JMG10112022

    basin[0].acc_year.sat_deficit = 0.0;
    basin[0].acc_year.sat_deficit_z = 0.0;
#endif

    return;
} /*end output_yearly_basin*/
