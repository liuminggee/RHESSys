/*--------------------------------------------------------------*/
/*                                                              */ 
/*		compute_annual_turnover				*/
/*                                                              */
/*  NAME                                                        */
/*		compute_annual_turnover				*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/* void compute_annual_turnover(				*/ 
/*			struct epconst_struct  epc,		*/
/*			struct epvar_struct	*epv,		*/
/*			struct phenology_struct *phen,		*/
/*			struct cstate_struct	*cs,		*/
/*			struct nstate_struct	*ns,		*/
/*			int grow_flag				*/
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*	computes livewood turnover (stem and coarse root)	*/
/*                                                              */
/*                                                              */
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*	P.Thornton (1998) version of BIOME_bgc			*/
/*                                                              */
/*                                                              */
/*--------------------------------------------------------------*/

#include <stdio.h>
#include "rhessys.h"

int compute_annual_turnover( 
							struct epconst_struct	epc,
							struct epvar_struct	*epv,
							struct cstate_struct *cs)
{
	/*------------------------------------------------------*/
	/*	Local function declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
	
	int ok=1;
	if ( (epc.veg_type == TREE) ){
		epv->day_livestem_turnover = (cs->live_stemc
			+ cs->livestemc_transfer + cs->livestemc_store)
			* epc.livewood_turnover / 365;
		epv->day_livecroot_turnover = (cs->live_crootc
			+ cs->livecrootc_transfer + cs->livecrootc_store)
			* epc.livewood_turnover	/ 365;
	}
	if (epc.veg_type == GRASS)
		epv->day_deadleaf_turnover = epc.deadleaf_turnover
		* (cs->dead_leafc ) / 365;

	return (!ok);
} /* end compute_annual_turnover */

//10012025LML add functions for check mass balance
double get_total_plant_C(struct cstate_struct *cs, bool isTREE) {
    //return kgC/m2 of strata; ignor beetle related pool NOW
    double total = 0.0;
    total += cs->cpool +
             cs->leafc +
             cs->dead_leafc +
             cs->frootc +
             cs->leafc_transfer +
             cs->frootc_transfer +
             cs->gresp_transfer +
             cs->leafc_store +
             cs->frootc_store +
             cs->gresp_store;
    if (isTREE) {
        total += cs->live_stemc +
                 cs->dead_stemc +
                 cs->live_crootc +
                 cs->dead_crootc +
                 cs->livestemc_transfer +
                 cs->deadstemc_transfer +
                 cs->livecrootc_transfer +
                 cs->deadcrootc_transfer +
                 cs->livestemc_store +
                 cs->deadstemc_store +
                 cs->livecrootc_store +
                 cs->deadcrootc_store +
                 cs->cwdc;
    }
    return total;
}

double get_total_plant_N(struct nstate_struct *ns, bool isTREE) {
    //return kgC/m2 of strata; ignor beetle related pool NOW
    double total = 0.0;
    total += ns->npool +
             ns->leafn +
             ns->dead_leafn +
             ns->frootn +
             ns->leafn_transfer +
             ns->frootn_transfer +
             ns->leafn_store +
             ns->frootn_store +
             ns->retransn;
    if (isTREE) {
        total += ns->live_stemn +
                 ns->dead_stemn +
                 ns->live_crootn +
                 ns->dead_crootn +
                 ns->livestemn_transfer +
                 ns->deadstemn_transfer +
                 ns->livecrootn_transfer +
                 ns->deadcrootn_transfer +
                 ns->livestemn_store +
                 ns->deadstemn_store +
                 ns->livecrootn_store +
                 ns->deadcrootn_store +
                 ns->cwdn;
    }
    return total;
}

double get_total_litter_C(struct litter_c_object *clitr) {
    //return kgC/m2 of patch
    double total = 0.0;
    total += clitr->litr1c +
             clitr->litr2c +
             clitr->litr3c +
             clitr->litr4c;
    return total;
}

double get_total_litter_N(struct litter_n_object *nlitr) {
    //return kgN/m2 of patch;
    double total = 0.0;
    total += nlitr->litr1n +
             nlitr->litr2n +
             nlitr->litr3n +
             nlitr->litr4n;
    return total;
}

double get_total_SOM_C(struct soil_c_object *somc) {
    //return kgC/m2 of patch
    //Note: including DOC
    double total = 0.0;
    total += somc->DOC +
             somc->soil1c +
             somc->soil2c +
             somc->soil3c +
             somc->soil4c;
    return total;
}

double get_total_SOM_N(struct soil_n_object *somn) {
    //return kgN/m2 of patch
    //Note: including DON
    double total = 0.0;
    total += somn->DON +
             somn->soil1n +
             somn->soil2n +
             somn->soil3n +
             somn->soil4n;
    return total;
}
