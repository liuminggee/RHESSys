/*--------------------------------------------------------------*/
/*                                                              */
/*		find_patch_in_hillslope					*/
/*                                                              */
/*  NAME                                                        */
/*		find_patch_in_hillslope					*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/*  find_patch_in_hillslope( struct basin_object *basin)   	*/
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/* finds a patch in a hillslope. */
/*                                                              */
/*                                                              */
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"

struct patch_object *find_patch_in_hillslope( int patch_ID, int zone_ID, 
								struct hillslope_object *hillslope)
{
	/*------------------------------------------------------*/
	/*	Local Function Definition. 							*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	int i;
	int fnd;
	struct zone_object *zone;
	struct patch_object *patch;

	/*--------------------------------------------------------------*/
	/*	find zones						*/
	/*--------------------------------------------------------------*/
	i = 0;
	fnd = 0;
	while ( (fnd == 0) && (i >= 0) && (i < hillslope[0].num_zones)) {
		if (hillslope[0].zones[i][0].ID == zone_ID) {
			zone = hillslope[0].zones[i];
			fnd = 1;
		}
		else {
			i += 1;
		}
	}
	if (fnd == 0) {
		fprintf(stderr,
			"FATAL ERROR: Could not find zone %d in hillslope %d in find_patch_in_hillslope\n",zone_ID, hillslope[0].ID);
        //exit(EXIT_FAILURE);
        return 0; //09222023LML
	}
	/*--------------------------------------------------------------*/
	/*	find patches						*/
	/*--------------------------------------------------------------*/
	i = 0;
	fnd = 0;
	while ( (fnd == 0) && (i >= 0) && (i < zone[0].num_patches)) {
		if (zone[0].patches[i][0].ID == patch_ID) {
			patch = zone[0].patches[i];
			fnd = 1;
		}
		else {
			i += 1;
		}
	}
	if (fnd == 0) {
		fprintf(stderr,
			"FATAL ERROR: Could not find patch %d in zone %d hill %d\n",
			patch_ID, zone_ID, hillslope[0].ID);
        //exit(EXIT_FAILURE);
        return 0; //09222023LML
	}
	return(patch);
}
