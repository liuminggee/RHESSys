/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_growth_patch						*/
/*																*/
/*	output_yearly_growth_patch - creates output_growth files objects.*/
/*																*/
/*	NAME														*/
/*	output_yearly_growth_patch - output_growths current contents of a patch.*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_growth_patch(int basinID, int hillID, int zoneID,	*/
/*					struct	patch_object	*patch,				*/
/*					struct	date	date,  						*/
/*					FILE 	*outfile)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	output_growths spatial structure according to commandline			*/
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

void	output_yearly_growth_patch(
                int basinID, int hillID, int zoneID,
                struct	patch_object	*patch,
                struct	date	current_date,
                FILE *outfile
#ifdef JMG_TRACKING
                ,struct simtime *simtime
#endif
        )
{
    /*--------------------------------------------------------------*/
    /*	Local function definition.									*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*	Local variable definition.									*/
    /*--------------------------------------------------------------*/

    double aleafc = 0;
    double aleafn = 0;
    double afrootc = 0;
    double afrootn = 0;
    double awoodc = 0;
    double awoodn = 0;
    struct	canopy_strata_object 	*strata;

    for (int layer=0 ; layer<patch[0].num_layers; layer++ ){
        for (int c=0 ; c<patch[0].layers[layer].count; c++ ){
            strata = patch[0].canopy_strata[(patch[0].layers[layer].strata[c])];
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
        }
    }

#ifdef JMG_MORE_YEARLY_OUTPUT
    double denitrif;
    denitrif = 0.0;
    denitrif += patch[0].soil_ns.nvolatilized_snk;
#endif

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic[] = "%4d %4d %4d %4d %3d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#else
    char out_basic[] = "%4d %4d %4d %4d %3d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#endif

#ifdef JMG_TRACKING
    char out_format[] = "%d %d %d ";
    strcat(out_format,out_basic);
#else
    char out_format[] = "";
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
        aleafc,
        aleafn,
        aleafc + afrootc + awoodc,
        aleafn + afrootn + awoodn,
    patch[0].litter_cs.litr1c + patch[0].litter_cs.litr2c + patch[0].litter_cs.litr3c + patch[0].litter_cs.litr4c,
        patch[0].soil_cs.soil1c + patch[0].soil_cs.soil2c + patch[0].soil_cs.soil3c + patch[0].soil_cs.soil4c,
    patch[0].litter_ns.litr1n + patch[0].litter_ns.litr2n + patch[0].litter_ns.litr3n + patch[0].litter_ns.litr4n,
        patch[0].soil_ns.soil1n + patch[0].soil_ns.soil2n + patch[0].soil_ns.soil3n + patch[0].soil_ns.soil4n,
    patch[0].soil_ns.nitrate,
    patch[0].soil_ns.sminn, patch[0].rootzone.depth*1000.0
#ifdef JMG_MORE_YEARLY_OUTPUT
        ,denitrif
#endif
             );
    return;
} /*end output_yearly_growth_patch*/
