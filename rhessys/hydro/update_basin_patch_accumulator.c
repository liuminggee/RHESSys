/*--------------------------------------------------------------------------------------*/
/* 											*/
/*			update_basin_patch_accumulator					*/
/*											*/
/*	NAME										*/
/*	update_basin_patch_accumulator.c - update accumulator variables at the end of day	*/
/*					this process is taken from compute_subsurface_routing.c	*/
/*	SYNOPSIS									*/
/*	void update_basin_patch_accumulator( 						*/
/*					struct command_line_object *command_line,	*/
/*					struct basin_object *basin			*/
/*					struct date current_date)			*/
/*											*/
/* 											*/
/*											*/
/*	OPTIONS										*/
/*											*/
/*											*/
/*	DESCRIPTION									*/
/*	this function is called in basin_daily_F at the end of each day, it was in  	*/
/*	the compute_subsurface_routing, 										*/
/*											*/
/*											*/
/*	PROGRAMMER NOTES								*/
/*											*/
/*											*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"
#ifdef LIU_TRACKING_BASIN_LITTERC
int acumulate_carbon_flux(struct cdayflux_patch_struct *target, struct cdayflux_patch_struct *source, double scale);
#endif
void update_basin_patch_accumulator(
            struct command_line_object 	*command_line,
            struct basin_object 		*basin,
            struct date		 	current_date)
{
    /*----------------------------------------------------------------------*/
    /* Local variables definition                                           */
    /*-----------------------------------------------------------------------*/
    double scale;
    double tmp;
    //struct patch_object *patch;
    int b,h,p,z,c,s;
#ifdef JMG_MORE_YEARLY_OUTPUT
    int layer;
    //10072022LML moved into inside the parallel code block
    //struct  canopy_strata_object    *stratum;
    //double leafc, leafn, stemc, stemn, rootc, rootn, AGBc, AGBn, plantc, plantn, BGBc, BGBn;
#endif
    /*----------------------------------------------------------------------*/
    /* initializations		                                           */
    /*----------------------------------------------------------------------*/

    /*---------------------------------------------------------------------*/
    /*update accumulator variables                                            */
    /*-----------------------------------------------------------------------*/

    #pragma omp parallel for private(h,z,p)
    for (h=0; h < basin->num_hillslopes; ++h) {
      for(z=0; z < basin->hillslopes[h][0].num_zones; ++z) {
        for (p=0; p < basin->hillslopes[h][0].zones[z][0].num_patches; p++) {
            struct patch_object *patch=basin->hillslopes[h]->zones[z]->patches[p];

             patch[0].acc_year_trans += (patch[0].transpiration_unsat_zone
                    + patch[0].transpiration_sat_zone);
            if ((command_line[0].output_flags.monthly == 1)
                    && (command_line[0].p != NULL )) {
                patch[0].acc_month.theta += patch[0].rootzone.S;
                patch[0].acc_month.sm_deficit +=
                        max(0.0,
                                (patch[0].sat_deficit-patch[0].rz_storage-patch[0].unsat_storage));
                patch[0].acc_month.et += (patch[0].transpiration_unsat_zone
                        + patch[0].evaporation_surf
                        + patch[0].exfiltration_unsat_zone
                        + patch[0].exfiltration_sat_zone
                        + +patch[0].transpiration_sat_zone
                        + patch[0].evaporation);
                patch[0].acc_month.denitrif += patch[0].ndf.denitrif;
                patch[0].acc_month.nitrif += patch[0].ndf.sminn_to_nitrate;
                patch[0].acc_month.mineralized +=
                        patch[0].ndf.net_mineralized;
                patch[0].acc_month.uptake += patch[0].ndf.sminn_to_npool;
                patch[0].acc_month.DON_loss +=
                        (patch[0].soil_ns.DON_Qout_total
                                - patch[0].soil_ns.DON_Qout_total);
                patch[0].acc_month.DOC_loss +=
                        (patch[0].soil_cs.DOC_Qout_total
                                - patch[0].soil_cs.DOC_Qout_total);
                patch[0].acc_month.psn += patch[0].net_plant_psn;
                patch[0].acc_month.snowpack =
                        max(patch[0].snowpack.water_equivalent_depth, patch[0].acc_month.snowpack);
                patch[0].acc_month.lai =
                        max(patch[0].acc_month.lai, patch[0].lai);
                patch[0].acc_month.leach += (patch[0].soil_ns.leach
                        + patch[0].surface_ns_leach);
                patch[0].acc_month.length += 1;

            }
            if ((command_line[0].output_flags.yearly == 1)
                    && (command_line[0].p != NULL ) || command_line[0].grow_flag > 0) {// Here change condition for calculating ET for landclim model
                patch[0].acc_year.length += 1;
                if ((patch[0].sat_deficit - patch[0].unsat_storage)
                        > command_line[0].thresholds[SATDEF])
                    patch[0].acc_year.num_threshold += 1;
                patch[0].acc_year.theta += patch[0].rootzone.S;
                patch[0].acc_year.denitrif += patch[0].ndf.denitrif;
                patch[0].acc_year.nitrif += patch[0].ndf.sminn_to_nitrate;
                patch[0].acc_year.mineralized +=
                        patch[0].ndf.net_mineralized;
                patch[0].acc_year.uptake += patch[0].ndf.sminn_to_npool;
                patch[0].acc_year.leach += (patch[0].soil_ns.leach
                        + patch[0].surface_ns_leach);
                patch[0].acc_year.DON_loss +=
                        (patch[0].soil_ns.DON_Qout_total
                                - patch[0].soil_ns.DON_Qout_total);
                patch[0].acc_year.DOC_loss +=
                        (patch[0].soil_cs.DOC_Qout_total
                                - patch[0].soil_cs.DOC_Qout_total);
                patch[0].acc_year.streamflow += patch[0].streamflow;
                patch[0].acc_year.Qout_total += patch[0].Qout_total;
                patch[0].acc_year.Qin_total += patch[0].Qin_total;
                patch[0].acc_year.psn += patch[0].net_plant_psn;
                patch[0].acc_year.TPET += (patch[0].PE + patch[0].PET);
                patch[0].acc_year.PET += patch[0].PET;
                patch[0].acc_year.PE += patch[0].PE;
                patch[0].acc_year.burn += patch[0].burn;
                patch[0].acc_year.potential_recharge +=
                        patch[0].rain_throughfall;
                patch[0].acc_year.potential_recharge_wyd +=
                        patch[0].rain_throughfall
                                * round(patch[0].acc_year.length);
                patch[0].acc_year.recharge += patch[0].recharge;
                patch[0].acc_year.recharge_wyd += patch[0].recharge
                        * round(patch[0].acc_year.length);

                if (close_enough(patch[0].snowpack.water_equivalent_depth, 0)
                        && (patch[0].acc_year.snowpack > 0)) {
                    if (patch[0].acc_year.meltday
                            < patch[0].acc_year.peaksweday)
                        patch[0].acc_year.meltday = round(
                                patch[0].acc_year.length);
                }

                if (patch[0].snowpack.water_equivalent_depth
                        > patch[0].acc_year.snowpack) {
                    patch[0].acc_year.peaksweday = round(
                            patch[0].acc_year.length);
                }

                patch[0].acc_year.snowpack =
                        max(patch[0].snowpack.water_equivalent_depth,
                                patch[0].acc_year.snowpack);

                /* transpiration water stress computations */
                double tmp = (patch[0].transpiration_unsat_zone
                        + patch[0].exfiltration_unsat_zone
                        + patch[0].exfiltration_sat_zone
                        + patch[0].evaporation_surf
                        + patch[0].transpiration_sat_zone
                        + patch[0].evaporation);
                patch[0].acc_year.et += tmp;
                // for calculating decomposition

                patch[0].acc_year.num_days += 1;
            if (patch[0].acc_year.num_days <=365) {
                patch[0].acc_year.et_decom += tmp;
                patch[0].acc_year.et_decom_mean = patch[0].acc_year.et_decom/patch[0].acc_year.num_days;
                //printf("[num_days %d], [ET %lf] \n", patch[0].acc_year.num_days, patch[0].acc_year.et_decom_mean  );
                }
            else {
                patch[0].acc_year.et_decom_mean = (364 * patch[0].acc_year.et_decom_mean + tmp) / 365; //give more important ratio to ET
               // printf("[num_days %d], [ET %lf] \n", patch[0].acc_year.num_days, patch[0].acc_year.et_decom_mean  );
                                                                }

                tmp = (patch[0].transpiration_unsat_zone
                        + patch[0].transpiration_sat_zone);
                patch[0].acc_year.trans += tmp;

                //if(current_date.month == 9 & current_date.day == 1)
                //printf("\n The [acc_year.et %lf], [acc_year.trans %lf], [PET %lf], \n", patch[0].acc_year.et*1000, patch[0].acc_year.trans*1000, patch[0].acc_year.PET*1000);

                patch[0].acc_year.day7trans = (tmp / 14
                        + 13 / 14 * patch[0].acc_year.day7trans);
                patch[0].acc_year.day7pet = (patch[0].PET + patch[0].PE)
                        / 14 + 13 / 14 * patch[0].acc_year.day7pet;
                if (patch[0].acc_year.day7pet > patch[0].acc_year.maxpet) {
                    patch[0].acc_year.maxpet = patch[0].acc_year.day7pet;
                    patch[0].acc_year.rec_pet_wyd = 0;
                    patch[0].acc_year.max_pet_wyd = patch[0].acc_year.wyd;
                }

                if ((patch[0].acc_year.day7trans
                        > patch[0].acc_year.maxtrans)) {
                    patch[0].acc_year.maxtrans =
                            patch[0].acc_year.day7trans;
                    patch[0].acc_year.rec_wyd = 0;
                }

                if ((patch[0].acc_year.rec_wyd == 0)
                        && (patch[0].acc_year.day7trans
                                < patch[0].acc_year.maxtrans * 0.5)) {
                    patch[0].acc_year.rec_wyd = patch[0].acc_year.wyd;
                }

                if ((patch[0].acc_year.rec_pet_wyd == 0)
                        && (patch[0].acc_year.day7pet
                                < patch[0].acc_year.maxpet * 0.5)) {
                    patch[0].acc_year.rec_pet_wyd = patch[0].acc_year.wyd;
                }

                tmp = (patch[0].transpiration_unsat_zone
                        + patch[0].exfiltration_unsat_zone
                        + patch[0].exfiltration_sat_zone
                        + patch[0].evaporation_surf
                        + +patch[0].transpiration_sat_zone
                        + patch[0].evaporation);

                if ((patch[0].PET + patch[0].PE - tmp)
                        > patch[0].acc_year.sm_deficit)
                    patch[0].acc_year.sm_deficit = (patch[0].PET
                            + patch[0].PE - tmp);
                patch[0].acc_year.lai =
                        max(patch[0].acc_year.lai, patch[0].lai);

                tmp = patch[0].sat_deficit - patch[0].unsat_storage
                        - patch[0].rz_storage;
                if (tmp <= 0)
                    patch[0].acc_year.ndays_sat += 1;

                if (patch[0].rootzone.S > 0.7)
                    patch[0].acc_year.ndays_sat70 += 1;

                tmp =
                        max(0.0, (patch[0].rootzone.field_capacity/patch[0].rootzone.potential_sat -
                                        patch[0].wilting_point*patch[0].soil_defaults[0][0].porosity_0))
                                / 2.0
                                + patch[0].wilting_point
                                        * patch[0].soil_defaults[0][0].porosity_0;

                if ((patch[0].rootzone.S < tmp) && (current_date.month < 10)
                        && (patch[0].acc_year.midsm_wyd == 0)
                        && (patch[0].snowpack.water_equivalent_depth <= 0.0))
                    patch[0].acc_year.midsm_wyd = patch[0].acc_year.wyd;

                patch[0].acc_year.wyd = patch[0].acc_year.wyd + 1;

#ifdef JMG_MORE_YEARLY_OUTPUT
                patch[0].acc_year.n_deposition += basin->hillslopes[h][0].zones[z][0].ndep_NO3 + basin->hillslopes[h][0].zones[z][0].ndep_NH4; // JMG09272022
                patch[0].acc_year.soilc += patch[0].soil_cs.soil1c + patch[0].soil_cs.soil2c + patch[0].soil_cs.soil3c + patch[0].soil_cs.soil4c; // JMG09272022
                patch[0].acc_year.soiln += patch[0].soil_ns.soil1n + patch[0].soil_ns.soil2n + patch[0].soil_ns.soil3n + patch[0].soil_ns.soil4n; // JMG09272022
                patch[0].acc_year.litrc += patch[0].litter_cs.litr1c + patch[0].litter_cs.litr2c + patch[0].litter_cs.litr3c + patch[0].litter_cs.litr4c; // JMG09272022
                patch[0].acc_year.litrn += patch[0].litter_ns.litr1n + patch[0].litter_ns.litr2n + patch[0].litter_ns.litr3n + patch[0].litter_ns.litr4n; // JMG09272022

                for ( layer=0 ; layer<patch[0].num_layers; layer++ ){
                    for ( c=0 ; c<patch[0].layers[layer].count; c++ ){
                        struct  canopy_strata_object *stratum = patch[0].canopy_strata[(patch[0].layers[layer].strata[c])];

                        double leafc = stratum[0].cover_fraction	* (stratum[0].cs.leafc
                            + stratum[0].cs.leafc_store + stratum[0].cs.leafc_transfer + stratum[0].cs.dead_leafc);
                            //* patch[0].area;

                        double stemc = stratum[0].cover_fraction * (stratum[0].cs.live_stemc
                                + stratum[0].cs.dead_stemc
                                + stratum[0].cs.livestemc_store
                                + stratum[0].cs.deadstemc_store
                                + stratum[0].cs.livestemc_transfer
                                + stratum[0].cs.deadstemc_transfer
                                + stratum[0].cs.cwdc + stratum[0].cs.cpool);

                        double rootc = stratum[0].cover_fraction * (stratum[0].cs.live_crootc
                            + stratum[0].cs.dead_crootc
                            + stratum[0].cs.livecrootc_store
                            + stratum[0].cs.deadcrootc_store
                            + stratum[0].cs.livecrootc_transfer
                            + stratum[0].cs.deadcrootc_transfer
                            + stratum[0].cs.frootc
                            + stratum[0].cs.frootc_store
                            + stratum[0].cs.frootc_transfer);

                        double AGBc = (leafc + stemc);
                        double plantc = (AGBc + rootc);
                        double BGBc = (rootc);

                        double leafn = stratum[0].cover_fraction	* (stratum[0].ns.leafn
                            + stratum[0].ns.leafn_store + stratum[0].ns.leafn_transfer + stratum[0].ns.dead_leafn);
                            //* patch[0].area;

                        double stemn = stratum[0].cover_fraction * (stratum[0].ns.live_stemn
                                + stratum[0].ns.dead_stemn
                                + stratum[0].ns.livestemn_store
                                + stratum[0].ns.deadstemn_store
                                + stratum[0].ns.livestemn_transfer
                                + stratum[0].ns.deadstemn_transfer
                                + stratum[0].ns.cwdn + stratum[0].ns.npool);

                        double rootn = stratum[0].cover_fraction * (stratum[0].ns.live_crootn
                            + stratum[0].ns.dead_crootn
                            + stratum[0].ns.livecrootn_store
                            + stratum[0].ns.deadcrootn_store
                            + stratum[0].ns.livecrootn_transfer
                            + stratum[0].ns.deadcrootn_transfer
                            + stratum[0].ns.frootn
                            + stratum[0].ns.frootn_store
                            + stratum[0].ns.frootn_transfer);

                        double AGBn = (leafn + stemn);
                        // plantn = 0.0;
                        double plantn = (AGBn + rootn);
                        // plantn = (leafn + stemn + rootn);
                        double BGBn = (rootn);

                        patch[0].acc_year.plantc += plantc; // JMG09272022
                        patch[0].acc_year.plantn += plantn; // JMG09272022
                        patch[0].acc_year.AGBc += AGBc; // JMG09272022
                        patch[0].acc_year.BGBc += BGBc; // JMG09272022

                    }
                }
#endif
             } /* end if */
        } /* end of p*/
      } /* end of z*/
    } /* end of h*/
    if (command_line[0].b != NULL) {
      for (h=0; h < basin->num_hillslopes; ++h) {
        for(z=0; z < basin->hillslopes[h][0].num_zones; ++z) {
            for (p=0; p < basin->hillslopes[h][0].zones[z][0].num_patches; p++) {
                struct patch_object *patch=basin->hillslopes[h]->zones[z]->patches[p];
                if (command_line[0].output_flags.monthly == 1) {
                    scale = patch[0].area / basin[0].area;
                    basin[0].acc_month.streamflow += (patch[0].streamflow)
                            * scale;
                    basin[0].acc_month.et += (patch[0].transpiration_unsat_zone
                            + patch[0].evaporation_surf
                            + patch[0].exfiltration_unsat_zone
                            + patch[0].exfiltration_sat_zone
                            + patch[0].transpiration_sat_zone
                            + patch[0].evaporation) * scale;
                    basin[0].acc_month.denitrif += patch[0].ndf.denitrif
                            * scale;
                    basin[0].acc_month.nitrif += patch[0].ndf.sminn_to_nitrate
                            * scale;
                    basin[0].acc_month.mineralized +=
                            patch[0].ndf.net_mineralized * scale;
                    basin[0].acc_month.uptake += patch[0].ndf.sminn_to_npool
                            * scale;
                    basin[0].acc_month.DON_loss +=
                            (patch[0].soil_ns.DON_Qout_total
                                    - patch[0].soil_ns.DON_Qin_total) * scale;
                    basin[0].acc_month.DOC_loss +=
                            (patch[0].soil_cs.DOC_Qout_total
                                    - patch[0].soil_cs.DOC_Qin_total) * scale;
                    basin[0].acc_month.length += 1;
                    basin[0].acc_month.stream_NO3 += patch[0].streamflow_NO3
                            * scale;
                    basin[0].acc_month.stream_NH4 += patch[0].streamflow_NH4
                            * scale;
                    basin[0].acc_month.stream_DON += patch[0].streamflow_DON
                            * scale;
                    basin[0].acc_month.stream_DOC += patch[0].streamflow_DOC
                            * scale;
                    basin[0].acc_month.psn += patch[0].net_plant_psn * scale;
                    basin[0].acc_month.lai += patch[0].lai * scale;
                    basin[0].acc_month.leach += (patch[0].soil_ns.leach
                            + patch[0].surface_ns_leach) * scale;

#ifdef LIU_TRACKING_BASIN_LITTERC
                    acumulate_carbon_flux(&basin[0].acc_month.cdf,&patch[0].cdf,scale);
#endif

                }

                if (command_line[0].output_flags.yearly == 1) {
                    scale = patch[0].area / basin[0].area;
                    basin[0].acc_year.length += 1;
                    basin[0].acc_year.leach += (patch[0].soil_ns.leach
                            + patch[0].surface_ns_leach) * scale;
                    basin[0].acc_year.stream_NH4 += patch[0].streamflow_NH4
                            * scale;
                    basin[0].acc_year.stream_NO3 += patch[0].streamflow_NO3
                            * scale;
                    basin[0].acc_year.denitrif += patch[0].ndf.denitrif * scale;
                    basin[0].acc_year.nitrif += patch[0].ndf.sminn_to_nitrate
                            * scale;
                    basin[0].acc_year.mineralized +=
                            patch[0].ndf.net_mineralized * scale;
                    basin[0].acc_year.uptake += patch[0].ndf.sminn_to_npool
                            * scale;
                    basin[0].acc_year.DON_loss +=
                            (patch[0].soil_ns.DON_Qout_total
                                    - patch[0].soil_ns.DON_Qin_total) * scale;
                    basin[0].acc_year.DOC_loss +=
                            (patch[0].soil_cs.DOC_Qout_total
                                    - patch[0].soil_cs.DOC_Qin_total) * scale;
                    basin[0].acc_year.stream_DON += patch[0].streamflow_DON
                            * scale;
                    basin[0].acc_year.stream_DOC += patch[0].streamflow_DOC
                            * scale;
                    basin[0].acc_year.psn += patch[0].net_plant_psn * scale;

                    basin[0].acc_year.TPET += (patch[0].PE + patch[0].PET)
                            * scale;
                    basin[0].acc_year.PET += patch[0].PET * scale;
                    basin[0].acc_year.PE += patch[0].PE	* scale;

                    basin[0].acc_year.et += (patch[0].evaporation
                            + patch[0].evaporation_surf
                            + patch[0].exfiltration_unsat_zone
                            + patch[0].exfiltration_sat_zone
                            + patch[0].transpiration_unsat_zone
                            + patch[0].transpiration_sat_zone) * scale;
                    basin[0].acc_year.streamflow += (patch[0].streamflow)
                            * scale;
                    basin[0].acc_year.lai += patch[0].lai * scale;
                }
            } /* end of p*/
        } /* end of z*/
      } /* end of h*/
    } /* end if */


    return;
} /* end of update_basin_patch_accumulator.c */

#ifdef LIU_TRACKING_BASIN_LITTERC
int acumulate_carbon_flux(struct cdayflux_patch_struct *target, struct cdayflux_patch_struct *source, double scale)
{
    target->do_litr1c_loss += scale * source->do_litr1c_loss;
    target->do_litr2c_loss += scale * source->do_litr2c_loss;
    target->do_litr3c_loss += scale * source->do_litr3c_loss;
    target->do_litr4c_loss += scale * source->do_litr4c_loss;

    target->do_soil1c_loss += scale * source->do_soil1c_loss;
    target->do_soil2c_loss += scale * source->do_soil2c_loss;
    target->do_soil3c_loss += scale * source->do_soil3c_loss;
    target->do_soil4c_loss += scale * source->do_soil4c_loss;
    target->total_DOC_loss += scale * source->total_DOC_loss;
    target->DOC_to_gw += scale * source->DOC_to_gw;


    target->plitr1c_loss += scale * source->plitr1c_loss;
    target->plitr2c_loss += scale * source->plitr2c_loss;
    target->plitr3c_loss += scale * source->plitr3c_loss;
    target->plitr4c_loss += scale * source->plitr4c_loss;
    target->psoil1c_loss += scale * source->psoil1c_loss;
    target->psoil2c_loss += scale * source->psoil2c_loss;
    target->psoil3c_loss += scale * source->psoil3c_loss;
    target->psoil4c_loss += scale * source->psoil4c_loss;
    target->kl4 += scale * source->kl4;

    target->leafc_to_litr1c += scale * source->leafc_to_litr1c;
    target->leafc_to_litr2c += scale * source->leafc_to_litr2c;
    target->leafc_to_litr3c += scale * source->leafc_to_litr3c;
    target->leafc_to_litr4c += scale * source->leafc_to_litr4c;
    target->frootc_to_litr1c += scale * source->frootc_to_litr1c;
    target->frootc_to_litr2c += scale * source->frootc_to_litr2c;
    target->frootc_to_litr3c += scale * source->frootc_to_litr3c;
    target->frootc_to_litr4c += scale * source->frootc_to_litr4c;
    target->litr1c_to_soil1c += scale * source->litr1c_to_soil1c;
    target->litr2c_to_soil2c += scale * source->litr2c_to_soil2c;
    target->litr3c_to_litr2c += scale * source->litr3c_to_litr2c;
    target->litr4c_to_soil3c += scale * source->litr4c_to_soil3c;
    target->soil1c_to_soil2c += scale * source->soil1c_to_soil2c;
    target->soil2c_to_soil3c += scale * source->soil2c_to_soil3c;
    target->soil3c_to_soil4c += scale * source->soil3c_to_soil4c;
    target->cwdc_to_litr2c += scale * source->cwdc_to_litr2c;
    target->cwdc_to_litr3c += scale * source->cwdc_to_litr3c;
    target->cwdc_to_litr4c += scale * source->cwdc_to_litr4c;
    target->stemc_to_litr1c += scale * source->stemc_to_litr1c;
    target->stemc_to_cwdc += scale * source->stemc_to_cwdc;
    target->rootc_to_cwdc += scale * source->rootc_to_cwdc;

    target->snagc_to_cwdc += scale * source->snagc_to_cwdc;
    target->litr1c_hr += scale * source->litr1c_hr;
    target->litr2c_hr += scale * source->litr2c_hr;
    target->litr3c_hr += scale * source->litr3c_hr;
    target->litr4c_hr += scale * source->litr4c_hr;
    target->soil1c_hr += scale * source->soil1c_hr;
    target->soil2c_hr += scale * source->soil2c_hr;
    target->soil3c_hr += scale * source->soil3c_hr;
    target->soil4c_hr += scale * source->soil4c_hr;

    target->m_leafc_to_litr1c += scale * source->m_leafc_to_litr1c;
    target->m_leafc_to_litr2c += scale * source->m_leafc_to_litr2c;
    target->m_leafc_to_litr3c += scale * source->m_leafc_to_litr3c;
    target->m_leafc_to_litr4c += scale * source->m_leafc_to_litr4c;
    target->m_frootc_to_litr1c += scale * source->m_frootc_to_litr1c;
    target->m_frootc_to_litr2c += scale * source->m_frootc_to_litr2c;
    target->m_frootc_to_litr3c += scale * source->m_frootc_to_litr3c;
    target->m_frootc_to_litr4c += scale * source->m_frootc_to_litr4c;

    target->m_leafc_store_to_litr1c += scale * source->m_leafc_store_to_litr1c;
    target->m_frootc_store_to_litr1c += scale * source->m_frootc_store_to_litr1c;
    target->m_livestemc_store_to_litr1c += scale * source->m_livestemc_store_to_litr1c;
    target->m_deadstemc_store_to_litr1c += scale * source->m_deadstemc_store_to_litr1c;
    target->m_livecrootc_store_to_litr1c += scale * source->m_livecrootc_store_to_litr1c;
    target->m_deadcrootc_store_to_litr1c += scale * source->m_deadcrootc_store_to_litr1c;

    target->m_leafc_transfer_to_litr1c += scale * source->m_leafc_transfer_to_litr1c;
    target->m_frootc_transfer_to_litr1c += scale * source->m_frootc_transfer_to_litr1c;
    target->m_livestemc_transfer_to_litr1c += scale * source->m_livestemc_transfer_to_litr1c;
    target->m_deadstemc_transfer_to_litr1c += scale * source->m_deadstemc_transfer_to_litr1c;
    target->m_livecrootc_transfer_to_litr1c += scale * source->m_livecrootc_transfer_to_litr1c;
    target->m_deadcrootc_transfer_to_litr1c += scale * source->m_deadcrootc_transfer_to_litr1c;

    target->m_gresp_store_to_litr1c += scale * source->m_gresp_store_to_litr1c;
    target->m_gresp_transfer_to_litr1c += scale * source->m_gresp_transfer_to_litr1c;

    target->m_litr1c_to_atmos += scale * source->m_litr1c_to_atmos;
    target->m_litr2c_to_atmos += scale * source->m_litr2c_to_atmos;
    target->m_litr3c_to_atmos += scale * source->m_litr3c_to_atmos;
    target->m_litr4c_to_atmos += scale * source->m_litr4c_to_atmos;
    target->m_soil1c_to_atmos += scale * source->m_soil1c_to_atmos;
    target->m_soil2c_to_atmos += scale * source->m_soil2c_to_atmos;
    target->m_soil3c_to_atmos += scale * source->m_soil3c_to_atmos;
    target->m_soil4c_to_atmos += scale * source->m_soil4c_to_atmos;

    target->litterc_to_atmos += scale * source->litterc_to_atmos;
    target->litterc_to_soilc += scale * source->litterc_to_soilc;
    return 0;

}
#endif
