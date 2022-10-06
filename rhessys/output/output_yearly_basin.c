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
                             FILE *outfile)
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

struct	patch_object  *patch;
struct	zone_object	*zone;
struct hillslope_object *hillslope;

int h,z,p,c;
h = z = p = c = 0;

double andep, asoilc, asoiln, alitrc, alitrn, aplantc, aplantn, aAGBc, aBGBc, hill_area, aarea;
andep = asoilc = asoiln = alitrc = alitrn = aplantc = aplantn = aAGBc = aBGBc = hill_area = aarea = 0.0;

    for (h=0; h < basin[0].num_hillslopes; h++){
        hillslope = basin[0].hillslopes[h];
        hill_area = 0.0;
        for (z=0; z< hillslope[0].num_zones; z++){
            zone = hillslope[0].zones[z];
            for (p=0; p< zone[0].num_patches; p++){
                patch = zone[0].patches[p];

                andep += patch[0].acc_year.n_deposition * patch[0].area;
                asoilc += patch[0].acc_year.soilc/patch[0].acc_year.length * patch[0].area;
                asoiln += patch[0].acc_year.soiln/patch[0].acc_year.length * patch[0].area;
                alitrc += patch[0].acc_year.litrc/patch[0].acc_year.length * patch[0].area;
                alitrn += patch[0].acc_year.litrn/patch[0].acc_year.length * patch[0].area;
                aplantc += patch[0].acc_year.plantc/patch[0].acc_year.length * patch[0].area;
                aplantn += patch[0].acc_year.plantn/patch[0].acc_year.length * patch[0].area;
                aAGBc += patch[0].acc_year.AGBc/patch[0].acc_year.length * patch[0].area;
                aBGBc += patch[0].acc_year.BGBc/patch[0].acc_year.length * patch[0].area;

                // TODO: Move this to the correct location in the correct file; these values could be used in another output script (yearly_growth_patch_output.c for example)
                // As of now, if growth flag is not on, then this script is never called and these values are never reset
                /*--------------------------------------------------------------*/
                /*	reset accumulator variables				*/
                /*--------------------------------------------------------------*/
                patch[0].acc_year.n_deposition = 0.0; // JMG09272022
                patch[0].acc_year.soilc = 0.0; // JMG09272022
                patch[0].acc_year.soiln = 0.0; // JMG09272022
                patch[0].acc_year.litrc = 0.0; // JMG09272022
                patch[0].acc_year.litrn = 0.0; // JMG09272022
                patch[0].acc_year.plantc = 0.0; // JMG09272022
                patch[0].acc_year.plantn = 0.0; // JMG09272022
                patch[0].acc_year.AGBc = 0.0; // JMG09272022
                patch[0].acc_year.BGBc = 0.0; // JMG09272022

                aarea +=  patch[0].area;
                hill_area += patch[0].area;
            }
        }
    }

    andep /= aarea;
    asoilc /= aarea;
    asoiln /= aarea;
    alitrc /= aarea;
    alitrn /= aarea;
    aplantc /= aarea;
    aplantn /= aarea;
    aAGBc /= aarea;
    aBGBc /= aarea;
#endif

    check = fprintf(outfile,
#ifndef JMG_MORE_YEARLY_OUTPUT
        "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf\n",
#else
        "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
#endif
        current_date.year,
        basin[0].ID,
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
#ifdef JMG_MORE_YEARLY_OUTPUT
        ,andep*1000.0, // gN/m2
        asoilc,
        asoiln,
        alitrc,
        alitrn,
        aplantc,
        aplantn,
        aAGBc,
        aBGBc
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
    basin[0].acc_year.lai = 0.0;
    basin[0].acc_year.mineralized = 0.0;
    basin[0].acc_year.uptake = 0.0;
    return;
} /*end output_yearly_basin*/
