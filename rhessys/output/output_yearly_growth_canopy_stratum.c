/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_growth_canopy_stratum						*/
/*																*/
/*	output_yearly_growth_canopy_stratum - creates output_growth files objects.		*/
/*																*/
/*	NAME														*/
/*	output_yearly_growth_canopy_stratum - output_growths  */
/*			current contents of a canopy_stratum.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_growth_canopy_stratum(int basinID,	*/
/*					int hillID, int zoneID,int patchID,         */
/*					struct	canopy_stratum_object	*canopy_stratum,				*/
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

void	output_yearly_growth_canopy_stratum( int basinID, int hillID, int zoneID,
			int patchID,
			struct	canopy_strata_object	*stratum,
			struct	date	current_date,
			struct	command_line_object *command_line,
            FILE *outfile
#ifdef JMG_TRACKING
            ,struct simtime *simtime
#endif
                                             )
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
#ifdef JMG_MORE_YEARLY_OUTPUT
    double AGBc = // above ground biomass carbon weight (kgC/m2)
            stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc + // leafc
                    stratum[0].cs.live_stemc + stratum[0].cs.dead_stemc + stratum[0].cs.livestemc_store + stratum[0].cs.deadstemc_store + stratum[0].cs.livestemc_transfer + stratum[0].cs.deadstemc_transfer + // stemc
                    stratum[0].cs.cwdc + stratum[0].cs.cpool; // remaining biomass
            /*AGBc += stratum[0].cs.live_stemc + stratum[0].cs.dead_stemc + stratum[0].cs.livestemc_store + stratum[0].cs.deadstemc_store + stratum[0].cs.livestemc_transfer + stratum[0].cs.deadstemc_transfer + // stemc
                    stratum[0].cs.cpool; // remaining biomass
                    */

    double LAI = // projected leaf area index (epc.proj_lai)
            (stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc) * stratum[0].defaults[0][0].epc.proj_sla;

    double rootc = // live croot + live froot + dead croot
            (stratum[0].cs.frootc + stratum[0].cs.live_crootc + stratum[0].cs.dead_crootc);

#endif

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic[] = "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#else
    char out_basic[] = "%d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
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
        	patchID,
        	stratum[0].ID,
#ifndef JMG_MORE_YEARLY_OUTPUT
            (stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc) * stratum[0].defaults[0][0].epc.proj_sla,
        	stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc,
        	stratum[0].ns.leafn + stratum[0].ns.leafn_store + stratum[0].ns.leafn_transfer + stratum[0].ns.dead_leafn,
        	stratum[0].cs.frootc + stratum[0].cs.frootc_store + stratum[0].cs.frootc_transfer,
        	stratum[0].ns.frootn + stratum[0].ns.frootn_store + stratum[0].ns.frootn_transfer,
        	(stratum[0].cs.live_crootc + stratum[0].cs.livecrootc_store + stratum[0].cs.livecrootc_transfer +
             stratum[0].cs.live_stemc + stratum[0].cs.livestemc_store + stratum[0].cs.livestemc_transfer +
             stratum[0].cs.dead_crootc + stratum[0].cs.deadcrootc_store + stratum[0].cs.deadcrootc_transfer +
             stratum[0].cs.dead_stemc + stratum[0].cs.deadstemc_store + stratum[0].cs.deadstemc_transfer) ,
        	(stratum[0].ns.live_crootn + stratum[0].ns.livecrootn_store + stratum[0].ns.livecrootn_transfer +
             stratum[0].ns.live_stemn + stratum[0].ns.livestemn_store + stratum[0].ns.livestemn_transfer +
             stratum[0].ns.dead_crootn + stratum[0].ns.deadcrootn_store + stratum[0].ns.deadcrootn_transfer +
             stratum[0].ns.dead_stemn + stratum[0].ns.deadstemn_store + stratum[0].ns.deadstemn_transfer) ,
        	stratum[0].cs.cwdc,
        	stratum[0].ns.cwdn,
            stratum[0].acc_year.psn,
            stratum[0].acc_year.mr,
            stratum[0].acc_year.gr,
            stratum[0].acc_year.minNSC,
            stratum[0].cs.mortality_fract,
            stratum[0].cs.snagc + stratum[0].cs.delay_snagc, // output the snag pool for beetle attack both carbon and nitrogen
            stratum[0].ns.snagn + stratum[0].ns.delay_snagn,
            stratum[0].cs.redneedlec + stratum[0].cs.delay_redneedlec,
            stratum[0].ns.redneedlen + stratum[0].ns.delay_redneedlen,
            stratum[0].cs.dead_rootc_beetle,
            stratum[0].ns.dead_rootn_beetle,
            stratum[0].epv.height, // the reason here height is different with fire.yearly, is fire.yearly is before burning but, stratum.yearly; if turn off the fire effect they should be the same
            stratum[0].rootzone.depth*1000.0
#else
            AGBc, // above ground biomass (kg-C/m2)
            stratum[0].epv.height, // canopy height (m)
            LAI, // Leaf Area Index (m2/m2)
            stratum[0].acc_year.psn, // yearly agg gross photosynthesis
            stratum[0].acc_year.mr, // yearly agg maintainance respiration
            stratum[0].acc_year.gr,  // yearly agg growth respiration
            stratum[0].cs.live_stemc + stratum[0].cs.dead_stemc + stratum[0].cs.livestemc_store + stratum[0].cs.deadstemc_store + stratum[0].cs.livestemc_transfer + stratum[0].cs.deadstemc_transfer, // stemc
            stratum[0].cs.leafc + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc, // leafc
            stratum[0].cs.frootc + stratum[0].cs.frootc_store + stratum[0].cs.frootc_transfer +
            stratum[0].cs.live_crootc + stratum[0].cs.livecrootc_store + stratum[0].cs.livecrootc_transfer +
            stratum[0].cs.dead_crootc + stratum[0].cs.deadcrootc_store + stratum[0].cs.deadcrootc_transfer, // rootc
            stratum[0].rootzone.depth*1000.0, // rootdepth (mm)
            rootc
#endif
            );

    stratum[0].acc_year.mr = 0.0;                                               //06012022LML
    stratum[0].acc_year.gr = 0.0;                                               //06012022LML
  if (command_line[0].f == NULL) { //If there is fire yearly growth output, set up set in the fire yearly growth output

    stratum[0].acc_year.psn = 0.0;
    stratum[0].acc_year.minNSC = -999;
#ifdef JMG_MORE_YEARLY_OUTPUT
    stratum[0].acc_year.LAIx = 0.0;
    stratum[0].acc_year.height = 0.0;
    stratum[0].acc_year.AGBc = 0.0;
    stratum[0].acc_year.BGBc = 0.0;
    stratum[0].acc_year.gpp = 0.0;
    stratum[0].acc_year.resp = 0.0;
    stratum[0].acc_year.npp = 0.0;
    stratum[0].acc_year.rootdepth = 0.0;
    stratum[0].acc_year.nuptake = 0.0;
    stratum[0].acc_year.length = 0.0;
#endif
  }

	return;
} /*end output_yearly_growth_canopy_stratum*/
