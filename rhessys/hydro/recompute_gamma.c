
/*--------------------------------------------------------------*/
/* 								*/
/*		recompute_gamma			*/
/*								*/
/*	NAME							*/
/*	recompute_gamma - recomputes total and individual 	*/
/*	gamma values for subsurface routing to replace		*/
/*	tographic gradients by water table gradients		*/ 
/*								*/
/*								*/
/*	SYNOPSIS						*/
/*	recompute_gamma(					*/
/*			 struct  patch_object *,		*/
/*				double)				*/
/*								*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*								*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "rhessys.h"
#include "phys_constants.h"

double recompute_gamma( struct patch_object *patch,
			 double total_gamma)
{
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/ 
	/*--------------------------------------------------------------*/ 
    double adjustment = 0.0;
    double water_table_z1, water_table_z2;
	/*--------------------------------------------------------------*/ 
	/*	for now, if water table is above the surface we		*/
	/*	we set saturation deficit to zero, since return flow	*/
	/*	is modelled separately and should not be taken into	*/
	/*	account in modelling surface gradients			*/
	/*--------------------------------------------------------------*/ 
    //z1 = patch[0].z;
	if (patch[0].sat_deficit_z > ZERO)	
        water_table_z1	 = (patch[0].z - patch[0].sat_deficit_z);
	else
        water_table_z1 = patch[0].z;
    int numnb = patch[0].innundation_list[0].num_neighbours;
    if (numnb > 0)
        for (int i =0; i < numnb; i++) {
            struct  neighbour_object *nb = &patch[0].innundation_list[0].neighbours[i];
            struct  patch_object *tpatch = nb->patch;
            if (tpatch[0].sat_deficit_z > 0)
                water_table_z2	 = (tpatch[0].z - tpatch[0].sat_deficit_z);
			else
                water_table_z2 = tpatch[0].z;
            if (fabs(patch[0].z-tpatch[0].z) > ZERO) {
                adjustment += max(((water_table_z1 - water_table_z2) / (patch[0].z - tpatch[0].z) *
                    nb->gamma),0.0);
            }
		}
	else
		adjustment = 1.0;
    return(adjustment * total_gamma);
} /*recompute_gamma*/
