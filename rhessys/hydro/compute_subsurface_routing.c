/*--------------------------------------------------------------*/
/* 											*/
/*					compute_subsurface_routing			*/
/*											*/
/*	compute_subsurface_routing.c - creates a patch object				*/
/*											*/
/*	NAME										*/
/*	compute_subsurface_routing.c - creates a patch object				*/
/*											*/
/*	SYNOPSIS									*/
/*	struct routing_list_object compute_subsurface_routing( 				*/
/*							struct command_line_object command */
/*							struct hillslope_object *hillslopen)	*/
/*				 			int,			 	*/
/*							struct date *current_date)	*/
/*											*/
/* 											*/
/*											*/
/*	OPTIONS										*/
/*											*/
/*											*/
/*	DESCRIPTION									*/
/*											*/
/*											*/
/*											*/
/*											*/
/*	PROGRAMMER NOTES								*/
/*											*/
/*	June 16, 98 C.Tague								*/
/*	limit drainage to maximum saturation deficit defined by soil depth		*/
/*											*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "params.h"
#include "rhessys.h"
#include "functions.h"
#include <omp.h>

void compute_subsurface_routing(struct command_line_object *command_line,
		struct hillslope_object *hillslope, int n_timesteps, struct date current_date) {
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/

	void update_drainage_stream(struct patch_object *,
			struct command_line_object *, double, int);

	void update_drainage_road(struct patch_object *,
			struct command_line_object *, double, int);

	void update_drainage_land(struct patch_object *,
			struct command_line_object *, double, int);

	double compute_infiltration(int, double, double, double, double, double,
			double, double, double, double, double);

	double compute_z_final(int, double, double, double, double, double);

//	double compute_N_leached(int, double, double, double, double, double,
//			double, double, double, double, double, double, double,double *);

	double compute_layer_field_capacity(int, int, double, double, double,
			double, double, double, double, double, double);

	double compute_unsat_zone_drainage(int, int, double, double, double, double,
			double, double);

	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	int i, d;
	int j, k;
	int grow_flag, verbose_flag;
	double time_int, tmp;
    double theta, m, Ksat;
    //double NO3_out, NH4_out, DON_out, DOC_out;
    double return_flow;
	double water_balance, infiltration;
	double innundation_depth;
	double hillslope_outflow;
	double hillslope_rz_storage;
	double hillslope_unsat_storage;
	double hillslope_sat_deficit;
	double hillslope_return_flow;
	double hillslope_detention_store;
	double hillslope_area;
	double preday_hillslope_unsat_storage;
	double preday_hillslope_rz_storage;
	double preday_hillslope_sat_deficit;
	double preday_sat_deficit;
	double preday_hillslope_return_flow;
	double preday_hillslope_detention_store;
	double add_field_capacity, rz_drainage, unsat_drainage;
    double streamflow, Qin_total, Qstr_total;
    //struct patch_object *patch;
    //struct patch_object *neigh;
	/*--------------------------------------------------------------*/
	/*	initializations						*/
	/*--------------------------------------------------------------*/
	grow_flag = command_line[0].grow_flag;
	verbose_flag = command_line[0].verbose_flag;

	time_int = 1.0 / n_timesteps;
	hillslope_outflow = 0.0;
	hillslope_area = 0.0;
	hillslope_unsat_storage = 0.0;
	hillslope_rz_storage = 0.0;
	hillslope_sat_deficit = 0.0;
	hillslope_return_flow = 0.0;
	hillslope_detention_store = 0.0;
	preday_hillslope_rz_storage = 0.0;
	preday_hillslope_unsat_storage = 0.0;
	preday_hillslope_sat_deficit = 0.0;
	preday_hillslope_return_flow = 0.0;
	preday_hillslope_detention_store = 0.0;
	streamflow = 0.0;
	Qin_total = 0.0;
	Qstr_total = 0.0;
	d = 0;
	hillslope[0].hillslope_outflow = 0.0;
	hillslope[0].hillslope_area = 0.0;
	hillslope[0].hillslope_unsat_storage = 0.0;
	hillslope[0].hillslope_rz_storage = 0.0;
	hillslope[0].hillslope_sat_deficit = 0.0;
	hillslope[0].hillslope_return_flow = 0.0;
	hillslope[0].hillslope_detention_store = 0.0;
	hillslope[0].preday_hillslope_rz_storage = 0.0;
	hillslope[0].preday_hillslope_unsat_storage = 0.0;
	hillslope[0].preday_hillslope_sat_deficit = 0.0;
	hillslope[0].preday_hillslope_return_flow = 0.0;
	hillslope[0].preday_hillslope_detention_store = 0.0;

	// Note: this assumes that the set of patches in the surface routing table is identical to
	//       the set of patches in the subsurface flow table

    //#pragma omp parallel for reduction(+ : preday_hillslope_rz_storage,preday_hillslope_unsat_storage,preday_hillslope_sat_deficit,preday_hillslope_return_flow,preday_hillslope_detention_store,hillslope_area)
    for (int i = 0; i < hillslope->route_list->num_patches; i++) {
        struct patch_object *patch = hillslope->route_list->list[i];
		patch[0].streamflow = 0.0;
        //07202023LML patch[0].return_flow = 0.0;
		patch[0].base_flow = 0.0;
		patch[0].infiltration_excess = 0.0;
		preday_hillslope_rz_storage += patch[0].rz_storage * patch[0].area;
		preday_hillslope_unsat_storage += patch[0].unsat_storage * patch[0].area;
		preday_hillslope_sat_deficit += patch[0].sat_deficit * patch[0].area;
		preday_hillslope_return_flow += patch[0].return_flow * patch[0].area;
		preday_hillslope_detention_store += patch[0].detention_store * patch[0].area;
		hillslope_area += patch[0].area;
		patch[0].Qin_total = 0.0;
		patch[0].Qout_total = 0.0;
        //07202023LML patch[0].Qin = 0.0; //07202023LML the Qin might not been processed yet from neighboring patches
		patch[0].Qout = 0.0;
		patch[0].surface_Qin = 0.0;
		patch[0].surface_Qout = 0.0;

		patch[0].overland_flow = 0.0;

		patch[0].preday_sat_deficit = patch[0].sat_deficit;

        patch[0].preday_sat_deficit_z = compute_z_final_from_surface(
                patch[0].soil_defaults[0], -1.0 * patch[0].sat_deficit);

		patch[0].interim_sat = patch[0].sat_deficit - patch[0].unsat_storage;
		if ((patch[0].sat_deficit - patch[0].unsat_storage) < ZERO)
			patch[0].S = 1.0;
		else
			patch[0].S = patch[0].unsat_storage / patch[0].sat_deficit;

		if (grow_flag > 0) {
            //07202023LML patch[0].soil_ns.NO3_Qin = 0.0;
			patch[0].soil_ns.NO3_Qout = 0.0;
            //07202023LML patch[0].soil_ns.NH4_Qin = 0.0;
			patch[0].soil_ns.NH4_Qout = 0.0;
			patch[0].soil_ns.NO3_Qin_total = 0.0;
			patch[0].soil_ns.NO3_Qout_total = 0.0;
			patch[0].soil_ns.NH4_Qin_total = 0.0;
			patch[0].soil_ns.NH4_Qout_total = 0.0;
			patch[0].streamflow_DON = 0.0;
			patch[0].streamflow_DOC = 0.0;
			patch[0].streamflow_NO3 = 0.0;
			patch[0].streamflow_NH4 = 0.0;
			patch[0].soil_ns.DON_Qin_total = 0.0;
			patch[0].soil_ns.DON_Qout_total = 0.0;
			patch[0].soil_cs.DOC_Qin_total = 0.0;
			patch[0].soil_cs.DOC_Qout_total = 0.0;
			patch[0].surface_DON_Qin_total = 0.0;
			patch[0].surface_DON_Qout_total = 0.0;
			patch[0].surface_DOC_Qin_total = 0.0;
			patch[0].surface_DOC_Qout_total = 0.0;
			patch[0].soil_ns.leach = 0.0;
			patch[0].surface_ns_leach = 0.0;
			patch[0].soil_ns.DON_Qout = 0.0;
            //07202023LML patch[0].soil_ns.DON_Qin = 0.0;
			patch[0].soil_cs.DOC_Qout = 0.0;
            //07202023LML patch[0].soil_cs.DOC_Qin = 0.0;
			patch[0].surface_DON_Qout = 0.0;
            //07202023LML patch[0].surface_DON_Qin = 0.0;
			patch[0].surface_DOC_Qout = 0.0;
            //07202023LML patch[0].surface_DOC_Qin = 0.0;

			patch[0].streamNO3_from_surface	= 0.0;
			patch[0].streamNO3_from_sub = 0.0;


		}
	}
    // For openmp assign local variables back to main
    hillslope[0].preday_hillslope_rz_storage = preday_hillslope_rz_storage;
	hillslope[0].preday_hillslope_unsat_storage = preday_hillslope_unsat_storage;
	hillslope[0].preday_hillslope_sat_deficit = preday_hillslope_sat_deficit ;
	hillslope[0].preday_hillslope_return_flow = preday_hillslope_return_flow ;
	hillslope[0].preday_hillslope_detention_store = preday_hillslope_detention_store;
    hillslope[0].hillslope_area = hillslope_area;

	/*--------------------------------------------------------------*/
	/*	calculate Qout for each patch and add appropriate	*/
	/*	proportion of subsurface outflow to each neighbour	*/
	/*--------------------------------------------------------------*/
	for (k = 0; k < n_timesteps; k++) {

		/*patch[0].preday_sat_deficit_z = compute_z_final(verbose_flag,
				patch[0].soil_defaults[0][0].porosity_0,
				patch[0].soil_defaults[0][0].porosity_decay,
				patch[0].soil_defaults[0][0].soil_depth, 0.0,
				-1.0 * patch[0].sat_deficit);
		patch[0].preday_sat_deficit = patch[0].sat_deficit;*/

        //10172022LML seems no benifit #pragma omp parallel for
		for (i = 0; i < hillslope->route_list->num_patches; i++) {
            struct patch_object *patch = hillslope->route_list->list[i];
//#ifdef LIU_OMP_PATCH_LOCK
//            omp_set_lock(&locks_patch[0][patch[0].Unique_ID_index]);
//            printf("compute_subsurface_routing:locked_patch:%d,thread:%d num_streads:%d\n"
//                   ,patch[0].ID,omp_get_thread_num(),omp_get_num_threads());
//#endif
            patch[0].hourly_subsur2stream_flow = 0;
			patch[0].hourly_sur2stream_flow = 0;
			patch[0].hourly_stream_flow = 0;
			patch[0].hourly[0].streamflow_NO3 = 0;
			patch[0].hourly[0].streamflow_NO3_from_sub = 0;
			patch[0].hourly[0].streamflow_NO3_from_surface = 0;
			/*--------------------------------------------------------------*/
			/*	for roads, saturated throughflow beneath road cut	*/
			/*	is routed to downslope patches; saturated throughflow	*/
			/*	above the cut and overland flow is routed to the stream	*/
			/*								*/
			/*	for streams, no routing - all exported from hillslope	*/
			/*								*/
			/*	regular land patches - route to downslope neighbours    */
			/*--------------------------------------------------------------*/
			if ((patch[0].drainage_type == ROAD)
					&& (command_line[0].road_flag == 1)) {
				update_drainage_road(patch, command_line, time_int,
						verbose_flag);
			} else if (patch[0].drainage_type == STREAM) {
				update_drainage_stream(patch, command_line, time_int,
						verbose_flag);
			} else {
				update_drainage_land(patch, command_line, time_int,
						verbose_flag);
			}
//#ifdef LIU_OMP_PATCH_LOCK
//            omp_unset_lock(&locks_patch[0][patch[0].Unique_ID_index]);
//#endif
		} /* end i */

		/*--------------------------------------------------------------*/
		/*	update soil moisture and nitrogen stores		*/
		/*	check water balance					*/
		/*--------------------------------------------------------------*/
        //#pragma omp parallel for
        for (int i = 0; i < hillslope->route_list->num_patches; i++) {
            struct patch_object *patch = hillslope->route_list->list[i];

			/*--------------------------------------------------------------*/
			/*	update subsurface 				*/
			/*-------------------------------------------------------------------------*/

			/*-------------------------------------------------------------------------*/
			/*	Recompute current actual depth to water table				*/
			/*-------------------------------------------------------------------------*/
			patch[0].sat_deficit += (patch[0].Qout - patch[0].Qin);

            //if (patch[0].Qout * 1000 >= 2) {
            //    printf("patch[0].Qout(mm/day):%lf patch[0].Qin(mm/day):%lf sat_deficit_z:%lf gamma:%lf\n",
            //           patch[0].Qout*1000*24, patch[0].Qin*1000*24,patch[0].sat_deficit_z,patch[0].innundation_list->gamma);
            //}


            patch[0].sat_deficit_z = compute_z_final_from_surface(patch[0].soil_defaults[0],
					-1.0 * patch[0].sat_deficit);

			if (grow_flag > 0) {
				patch[0].soil_ns.nitrate += (patch[0].soil_ns.NO3_Qin
						- patch[0].soil_ns.NO3_Qout);

                //printf("0 nitrate:%lf\n",patch[0].soil_ns.nitrate);

				patch[0].soil_ns.sminn += (patch[0].soil_ns.NH4_Qin
						- patch[0].soil_ns.NH4_Qout);
				patch[0].soil_cs.DOC += (patch[0].soil_cs.DOC_Qin
						- patch[0].soil_cs.DOC_Qout);
				patch[0].soil_ns.DON += (patch[0].soil_ns.DON_Qin
						- patch[0].soil_ns.DON_Qout);
			}

			/*--------------------------------------------------------------*/
			/*      Recompute 	soil moisture storage                   */
			/*--------------------------------------------------------------*/

			if (patch[0].sat_deficit > patch[0].rootzone.potential_sat) {
				patch[0].rootzone.S =
						min(patch[0].rz_storage / patch[0].rootzone.potential_sat, 1.0);
				patch[0].S = patch[0].unsat_storage
						/ (patch[0].sat_deficit
								- patch[0].rootzone.potential_sat);
			} else {
				patch[0].rootzone.S =
						min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)
								/ patch[0].rootzone.potential_sat, 1.0);
				patch[0].S =
						min(patch[0].rz_storage / patch[0].sat_deficit, 1.0);
			}

			/*--------------------------------------------------------------*/
			/*	reset iterative  patch fluxes to zero			*/
			/*--------------------------------------------------------------*/
			patch[0].soil_ns.leach += (patch[0].soil_ns.DON_Qout
					+ patch[0].soil_ns.NH4_Qout + patch[0].soil_ns.NO3_Qout
					- patch[0].soil_ns.NH4_Qin - patch[0].soil_ns.NO3_Qin
					- patch[0].soil_ns.DON_Qin);
			patch[0].surface_ns_leach += ((patch[0].surface_NO3_Qout
					- patch[0].surface_NO3_Qin)
					+ (patch[0].surface_NH4_Qout - patch[0].surface_NH4_Qin)
					+ (patch[0].surface_DON_Qout - patch[0].surface_DON_Qin));
			patch[0].Qin_total += patch[0].Qin + patch[0].surface_Qin;
			patch[0].Qout_total += patch[0].Qout + patch[0].surface_Qout;

			patch[0].surface_Qin = 0.0;
			patch[0].surface_Qout = 0.0;
			patch[0].Qin = 0.0;
			patch[0].Qout = 0.0;
			if (grow_flag > 0) {
				patch[0].soil_cs.DOC_Qin_total += patch[0].soil_cs.DOC_Qin;
				patch[0].soil_cs.DOC_Qout_total += patch[0].soil_cs.DOC_Qout;
				patch[0].soil_ns.NH4_Qin_total += patch[0].soil_ns.NH4_Qin;
				patch[0].soil_ns.NH4_Qout_total += patch[0].soil_ns.NH4_Qout;
				patch[0].soil_ns.NO3_Qin_total += patch[0].soil_ns.NO3_Qin;
				patch[0].soil_ns.NO3_Qout_total += patch[0].soil_ns.NO3_Qout;
				patch[0].soil_ns.DON_Qin_total += patch[0].soil_ns.DON_Qin;
				patch[0].soil_ns.DON_Qout_total += patch[0].soil_ns.DON_Qout;
				patch[0].surface_DON_Qin_total += patch[0].surface_DON_Qin;
				patch[0].surface_DON_Qout_total += patch[0].surface_DON_Qout;
				patch[0].surface_DOC_Qin_total += patch[0].surface_DOC_Qin;
				patch[0].surface_DOC_Qout_total += patch[0].surface_DOC_Qout;

				patch[0].soil_ns.NH4_Qin = 0.0;
				patch[0].soil_ns.NH4_Qout = 0.0;
				patch[0].soil_ns.NO3_Qin = 0.0;
				patch[0].soil_ns.NO3_Qout = 0.0;
				patch[0].soil_ns.DON_Qout = 0.0;
				patch[0].soil_ns.DON_Qin = 0.0;
				patch[0].soil_cs.DOC_Qout = 0.0;
				patch[0].soil_cs.DOC_Qin = 0.0;
				patch[0].surface_NH4_Qout = 0.0;
				patch[0].surface_NH4_Qin = 0.0;
				patch[0].surface_NO3_Qout = 0.0;
				patch[0].surface_NO3_Qin = 0.0;
				patch[0].surface_DON_Qout = 0.0;
				patch[0].surface_DON_Qin = 0.0;
				patch[0].surface_DOC_Qout = 0.0;
				patch[0].surface_DOC_Qin = 0.0;

			}
			/*--------------------------------------------------------------*/
			/*	finalize streamflow and saturation deficits		*/
			/*								*/
			/*	move any saturation excess into detention store		*/
			/*	(i.e this needs to be routed on the next time step)	*/
			/* 	some streamflow may have already been accumulated from 	*/
			/* 	redirected streamflow					*/
			/*	water balance calculations				*/
			/* only on last iteration					*/
			/* **** note that streamflow is updated sequentially		*/
			/*	i.e not at the end; it is similar to Qout, in		*/
			/*	that it accumulates flux in from patches		*/
			/*	(roads) that direct water to the stream			*/
			/*--------------------------------------------------------------*/
			//if (k >=0){// (n_timesteps - 1)) //incorporate Tungs bug fix
                         patch[0].hourly_stream_flow += patch[0].hourly_subsur2stream_flow
								+ patch[0].hourly_sur2stream_flow;

				/*---------------------------------------------------------------------*/
				/*update daily output: the patch[0].base_flow and patch[0].return_flow
				 * is the summation of 24 hours return_flow and base_flow from previous
				 * calculation*/
				/*-----------------------------------------------------------------------*/
				/*if(k==n_timesteps-1){
				    if (patch[0].drainage_type == STREAM) {
					    patch[0].streamflow += patch[0].return_flow
							    + patch[0].base_flow;
				    }
				}*/

		} /* end i */


        if (k == (n_timesteps -1)) //end day
        {
            double leached[LEACH_ELEMENT_counts];
            //10172022LML seems no benifit #pragma omp parallel for
            for (int i = 0; i < hillslope->route_list->num_patches; i++) {
              struct patch_object *patch = hillslope->route_list->list[i];
              double excess = 0;
#ifdef LIU_CHECK_WATER_BALANCE
              //07192023LML check waterbalance
              double pre_patch_surface_water = patch[0].detention_store + patch[0].return_flow;
              double pre_patch_total_deficit = patch[0].sat_deficit - (patch[0].rz_storage + patch[0].unsat_storage);
              double pre_patch_total = pre_patch_surface_water - pre_patch_total_deficit;
              double total_Q_leave_patch = 0;
#endif





//#ifdef LIU_OMP_PATCH_LOCK
//              omp_set_lock(&locks_patch[0][patch[0].Unique_ID_index]);
              ////            printf("compute_subsurface_routing:locked_patch:%d,thread:%d num_streads:%d\n"
              ////                   ,patch[0].ID,omp_get_thread_num(),omp_get_num_threads());
//#endif

              if ((patch[0].sat_deficit
                    - (patch[0].unsat_storage + patch[0].rz_storage))
                    < -1.0 * ZERO) {
                excess = -1.0
                        * (patch[0].sat_deficit - patch[0].unsat_storage
                                - patch[0].rz_storage);
                patch[0].detention_store += excess;

                //if (patch[0].ID == 81497)
                //printf("detention_store:%lf excess:%lf \n",
                //       patch[0].detention_store*1000,
                //       excess*1000);


                patch[0].sat_deficit = 0.0;
                patch[0].unsat_storage = 0.0;
                patch[0].rz_storage = 0.0;
                //double Nout = 0;

                if (grow_flag > 0) {
                    double lpools[] = {patch[0].soil_ns.nitrate,
                                      patch[0].soil_ns.sminn,
                                      patch[0].soil_ns.DON,
                                      patch[0].soil_cs.DOC
                                      };
                    //double *Nout =
                    double t = compute_N_leached_from_soildef(verbose_flag,
                                    lpools, excess, 0.0, 0.0,
                                    //patch[0].m,
                                    //patch[0].innundation_list[d].gamma
                                    //        / patch[0].area * time_int,
                                    patch[0].soil_defaults[0],
                                    leached,
                                    LEACH_ELEMENT_counts
                                    //patch[0].transmissivity_profile
                                    );
                    patch[0].surface_DOC += leached[LDOC];
                    patch[0].soil_cs.DOC -= leached[LDOC];

                    patch[0].surface_DON += leached[LDON];
                    patch[0].soil_ns.DON -= leached[LDON];

                    patch[0].surface_NO3 += leached[LNO3];
                    patch[0].soil_ns.nitrate -= leached[LNO3];

                    patch[0].surface_NH4 += leached[LNH4];
                    patch[0].soil_ns.sminn -= leached[LNH4];
                    //free(Nout);
                }
            }

            /*--------------------------------------------------------------*/
            /*	final overland flow routing				*/
            /*--------------------------------------------------------------*/
            //06012023LML added max
            patch[0].overland_flow += max(patch[0].detention_store
                                - patch[0].soil_defaults[0][0].detention_store_size,0);

            if (((excess = patch[0].detention_store
                    - patch[0].soil_defaults[0][0].detention_store_size)
                    > ZERO) && (patch[0].detention_store > ZERO)) {

                if (patch[0].drainage_type == STREAM) {
                    if (grow_flag > 0) {
                        patch[0].streamflow_DON += (excess
                                / patch[0].detention_store)
                                * patch[0].surface_DON;
                        patch[0].streamflow_DOC += (excess
                                / patch[0].detention_store)
                                * patch[0].surface_DOC;

                        patch[0].streamflow_NO3 += (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NO3;
                        patch[0].streamNO3_from_surface +=(excess
                                / patch[0].detention_store)
                                * patch[0].surface_NO3;
                        patch[0].hourly[0].streamflow_NO3 += (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NO3;
                        patch[0].hourly[0].streamflow_NO3_from_surface +=(excess
                                / patch[0].detention_store)
                                * patch[0].surface_NO3;



                        patch[0].streamflow_NH4 += (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NH4;
                        patch[0].surface_DON -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_DON;
                        patch[0].surface_DOC -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_DOC;
                        patch[0].surface_NO3 -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NO3;
                        patch[0].surface_NH4 -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NH4;
                    }
                    patch[0].return_flow += excess;
                    patch[0].detention_store -= excess;
                    patch[0].Qout_total += excess;
                    patch[0].hourly_sur2stream_flow += excess;

                } else {
                    /*--------------------------------------------------------------*/
                    /* determine which innundation depth to consider		*/
                    /*--------------------------------------------------------------*/
                    if (patch[0].num_innundation_depths > 0) {
                        innundation_depth = patch[0].detention_store;
                        d = 0;
                        while ((innundation_depth
                                > patch[0].innundation_list[d].critical_depth)
                                && (d < patch[0].num_innundation_depths - 1)) {
                            d++;
                        }
                    } else {
                        d = 0;
                    }

                    for (int j = 0; j < patch->surface_innundation_list[d].num_neighbours; j++) {
                        struct patch_object *neigh = patch->surface_innundation_list[d].neighbours[j].patch;
                        double Qout = excess * patch->surface_innundation_list[d].neighbours[j].gamma; //m
                        double NO3_out,NH4_out,DON_out,DOC_out,Nout;
                        if (grow_flag > 0) {
                            NO3_out = Qout / patch[0].detention_store
                                    * patch[0].surface_NO3;
                            NH4_out = Qout / patch[0].detention_store
                                    * patch[0].surface_NH4;
                            DON_out = Qout / patch[0].detention_store
                                    * patch[0].surface_DON;
                            DOC_out = Qout / patch[0].detention_store
                                    * patch[0].surface_DOC;
                            Nout = NO3_out + NH4_out + DON_out;
                        }
//#ifdef LIU_OMP_PATCH_LOCK
//                        omp_set_lock(&locks_patch[0][neigh[0].Unique_ID_index]);
//#endif
                        if (neigh[0].drainage_type == STREAM) {
                            neigh[0].Qin_total += Qout * patch[0].area
                                    / neigh[0].area;
                            //07202023LML since it may not be processes at end of day on return_flow (will be cleared)
                            neigh[0].return_flow += Qout * patch[0].area
                                    / neigh[0].area;

                            if (grow_flag > 0) {
                                neigh[0].streamflow_DOC += (DOC_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].streamflow_DON += (DON_out
                                        * patch[0].area / neigh[0].area);

                                neigh[0].streamflow_NO3 += (NO3_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].streamNO3_from_surface +=(NO3_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].hourly[0].streamflow_NO3 += (NO3_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].hourly[0].streamflow_NO3_from_sub +=(NO3_out
                                        * patch[0].area / neigh[0].area);



                                neigh[0].streamflow_NH4 += (NH4_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].surface_ns_leach += (Nout
                                        * patch[0].area / neigh[0].area);
                            }
                        } else {
                            neigh[0].Qin_total += Qout * patch[0].area
                                    / neigh[0].area;
                            neigh[0].detention_store += Qout * patch[0].area
                                    / neigh[0].area;

                            //if (neigh[0].ID == 81497)
                            //printf("detention_store:%lf Qout:%lf \n",
                            //       neigh[0].detention_store*1000,
                            //       Qout*1000);


                            if (grow_flag > 0) {
                                neigh[0].surface_DOC += (DOC_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].surface_DON += (DON_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].surface_NO3 += (NO3_out
                                        * patch[0].area / neigh[0].area);
                                neigh[0].surface_ns_leach -= (Nout
                                        * patch[0].area / neigh[0].area);
                                neigh[0].surface_NH4 += (NH4_out
                                        * patch[0].area / neigh[0].area);

                            }

                        }
//#ifdef LIU_OMP_PATCH_LOCK
//                        omp_unset_lock(&locks_patch[0][neigh[0].Unique_ID_index]);
//#endif
                    }
                    if (grow_flag > 0) {
                        patch[0].surface_DOC -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_DOC;
                        patch[0].surface_DON -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_DON;
                        patch[0].surface_NO3 -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NO3;


                        patch[0].surface_NH4 -= (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NH4;
                        patch[0].surface_ns_leach += (excess
                                / patch[0].detention_store)
                                * patch[0].surface_NO3;
                    }
                    patch[0].detention_store -= excess;
                    patch[0].Qout_total += excess;
#ifdef LIU_CHECK_WATER_BALANCE
                    total_Q_leave_patch += excess; //07192023LML
#endif
                }
            }

            /*-------------------------------------------------------------------------*/
            /*Recompute current actual depth to water table				*/
            /*-------------------------------------------------------------------------*/
            patch[0].sat_deficit_z = compute_z_final_from_surface(patch[0].soil_defaults[0],
                    -1.0 * patch[0].sat_deficit);

            /*--------------------------------------------------------------*/
            /* 	leave behind field capacity			*/
            /*	if sat deficit has been lowered			*/
            /*	this should be an interactive process, we will use 	*/
            /*	0th order approximation					*/
            /* 	we do not do this once sat def is below 0.9 soil depth	*/
            /*     we use 0.9 to prevent numerical instability		*/
            /*--------------------------------------------------------------*/
            if ((patch[0].sat_deficit_z > patch[0].preday_sat_deficit_z)
                    && (patch[0].sat_deficit_z
                            < patch[0].soil_defaults[0][0].soil_depth * 0.9)) {
                add_field_capacity = compute_layer_field_capacity_from_soildef(
                        command_line[0].verbose_flag,
                        patch[0].soil_defaults[0],
                        patch[0].sat_deficit_z, patch[0].sat_deficit_z,
                        patch[0].preday_sat_deficit_z);

                add_field_capacity = max(add_field_capacity, 0.0);
                patch[0].sat_deficit += add_field_capacity;

                if ((patch[0].sat_deficit_z > patch[0].rootzone.depth)
                        && (patch[0].preday_sat_deficit_z
                                > patch[0].rootzone.depth))
                    patch[0].unsat_storage += add_field_capacity;
                else
                    patch[0].rz_storage += add_field_capacity;
            }

            if (patch[0].rootzone.depth > ZERO) {
                if ((patch[0].sat_deficit > ZERO)
                        && close_enough(patch[0].rz_storage, 0.0)) {
                    add_field_capacity = compute_layer_field_capacity_from_soildef(
                            command_line[0].verbose_flag,
                            patch[0].soil_defaults[0],
                            patch[0].sat_deficit_z, patch[0].sat_deficit_z,
                            0.0);
                    add_field_capacity = max(add_field_capacity, 0.0);
                    patch[0].sat_deficit += add_field_capacity;
                    patch[0].rz_storage += add_field_capacity;
                }
            } else {
                if ((patch[0].sat_deficit > ZERO)
                        && close_enough(patch[0].unsat_storage, 0.0)) {
                    add_field_capacity = compute_layer_field_capacity_from_soildef(
                            command_line[0].verbose_flag,
                            patch[0].soil_defaults[0],
                            patch[0].sat_deficit_z, patch[0].sat_deficit_z,
                            0.0);
                    add_field_capacity = max(add_field_capacity, 0.0);
                    patch[0].sat_deficit += add_field_capacity;
                    patch[0].unsat_storage += add_field_capacity;
                }
            }

            /*--------------------------------------------------------------*/
            /* try to infiltrate this water					*/
            /* use time_int as duration */
            /*--------------------------------------------------------------*/

            if (patch[0].detention_store > ZERO)
                if (patch[0].rootzone.depth > ZERO) {
                    infiltration = compute_infiltration_patch(verbose_flag,
                            patch, patch[0].rootzone.S,
                            (patch[0].detention_store), time_int);
                } else {
                    infiltration = compute_infiltration_patch(verbose_flag,
                            patch, patch[0].S,
                            (patch[0].detention_store), time_int);
                }
            else
                infiltration = 0.0;
            /*--------------------------------------------------------------*/
            /* added an surface N flux to surface N pool	and		*/
            /* allow infiltration of surface N				*/
            /*--------------------------------------------------------------*/
            if ((grow_flag > 0) && (infiltration > ZERO)) {
                patch[0].soil_ns.DON += ((infiltration
                        / patch[0].detention_store) * patch[0].surface_DON);
                patch[0].soil_cs.DOC += ((infiltration
                        / patch[0].detention_store) * patch[0].surface_DOC);
                patch[0].soil_ns.nitrate += ((infiltration
                        / patch[0].detention_store) * patch[0].surface_NO3);

                //printf("4 nitrate:%lf\n",patch[0].soil_ns.nitrate);

                patch[0].surface_NO3 -= ((infiltration
                        / patch[0].detention_store) * patch[0].surface_NO3);
                patch[0].soil_ns.sminn += ((infiltration
                        / patch[0].detention_store) * patch[0].surface_NH4);
                patch[0].surface_NH4 -= ((infiltration
                        / patch[0].detention_store) * patch[0].surface_NH4);
                patch[0].surface_DOC -= ((infiltration
                        / patch[0].detention_store) * patch[0].surface_DOC);
                patch[0].surface_DON -= ((infiltration
                        / patch[0].detention_store) * patch[0].surface_DON);
            }

            /*--------------------------------------------------------------*/
            /*	Determine if the infifltration will fill up the unsat	*/
            /*	zone or not.						*/
            /*	We use the strict assumption that sat deficit is the	*/
            /*	amount of water needed to saturate the soil.		*/
            /*--------------------------------------------------------------*/

            if (infiltration
                    > patch[0].sat_deficit - patch[0].unsat_storage
                            - patch[0].rz_storage) {
                /*--------------------------------------------------------------*/
                /*		Yes the unsat zone will be filled so we may	*/
                /*		as well treat the unsat_storage and infiltration*/
                /*		as water added to the water table.		*/
                /*--------------------------------------------------------------*/
                patch[0].sat_deficit -= (infiltration
                        + patch[0].unsat_storage + patch[0].rz_storage);
                /*--------------------------------------------------------------*/
                /*		There is no unsat_storage left.			*/
                /*--------------------------------------------------------------*/
                patch[0].unsat_storage = 0.0;
                patch[0].rz_storage = 0.0;
                patch[0].field_capacity = 0.0;
                patch[0].rootzone.field_capacity = 0.0;
            } else if ((patch[0].sat_deficit
                    > patch[0].rootzone.potential_sat)
                    && (infiltration
                            > patch[0].rootzone.potential_sat
                                    - patch[0].rz_storage)) {
                /*------------------------------------------------------------------------------*/
                /*		Just add the infiltration to the rz_storage and unsat_storage	*/
                /*------------------------------------------------------------------------------*/
                patch[0].unsat_storage += infiltration
                        - (patch[0].rootzone.potential_sat
                                - patch[0].rz_storage);
                patch[0].rz_storage = patch[0].rootzone.potential_sat;
            }
            /* Only rootzone layer saturated - perched water table case */
            else if ((patch[0].sat_deficit > patch[0].rootzone.potential_sat)
                    && (infiltration
                            <= patch[0].rootzone.potential_sat
                                    - patch[0].rz_storage)) {
                /*--------------------------------------------------------------*/
                /*		Just add the infiltration to the rz_storage	*/
                /*--------------------------------------------------------------*/
                patch[0].rz_storage += infiltration;
            }

            else if ((patch[0].sat_deficit
                    <= patch[0].rootzone.potential_sat)
                    && (infiltration
                            <= patch[0].sat_deficit - patch[0].rz_storage
                                    - patch[0].unsat_storage)) {
                patch[0].rz_storage += patch[0].unsat_storage;
                /* transfer left water in unsat storage to rootzone layer */
                patch[0].unsat_storage = 0;
                patch[0].rz_storage += infiltration;
                patch[0].field_capacity = 0;
            }

            if (patch[0].sat_deficit < 0.0) {
                patch[0].detention_store -= (patch[0].sat_deficit
                        - patch[0].unsat_storage);
                patch[0].sat_deficit = 0.0;
                patch[0].unsat_storage = 0.0;
            }
            patch[0].detention_store -= infiltration;
            /*--------------------------------------------------------------*/
            /* recompute saturation deficit					*/
            /*--------------------------------------------------------------*/
            patch[0].sat_deficit_z = compute_z_final_from_surface(patch[0].soil_defaults[0],
                    -1.0 * patch[0].sat_deficit);

            /*--------------------------------------------------------------*/
            /*	compute new field capacity				*/
            /*--------------------------------------------------------------*/
            if (patch[0].sat_deficit_z < patch[0].rootzone.depth) {
                patch[0].rootzone.field_capacity =
                        compute_layer_field_capacity_from_soildef(
                                command_line[0].verbose_flag,
                                patch[0].soil_defaults[0],
                                patch[0].sat_deficit_z,
                                patch[0].rootzone.depth, 0.0);

                patch[0].field_capacity = 0.0;
            } else {

                patch[0].rootzone.field_capacity =
                        compute_layer_field_capacity_from_soildef(
                                command_line[0].verbose_flag,
                                patch[0].soil_defaults[0],
                                patch[0].sat_deficit_z,
                                patch[0].rootzone.depth, 0.0);

                patch[0].field_capacity = compute_layer_field_capacity_from_soildef(
                        command_line[0].verbose_flag,
                        patch[0].soil_defaults[0],
                        patch[0].sat_deficit_z, patch[0].sat_deficit_z, 0)
                        - patch[0].rootzone.field_capacity;
            }

            /*--------------------------------------------------------------*/
            /*      Recompute patch soil moisture storage                   */
            /*--------------------------------------------------------------*/
            if (patch[0].sat_deficit < ZERO) {
                patch[0].S = 1.0;
                patch[0].rootzone.S = 1.0;
                rz_drainage = 0.0;
                unsat_drainage = 0.0;
            } else if (patch[0].sat_deficit_z > patch[0].rootzone.depth) { /* Constant vertical profile of soil porosity */

                /*-------------------------------------------------------*/
                /*	soil drainage and storage update	     	 */
                /*-------------------------------------------------------*/
                patch[0].rootzone.S =
                        min(patch[0].rz_storage / patch[0].rootzone.potential_sat, 1.0);
                rz_drainage = compute_unsat_zone_drainage_patch(
                        command_line[0].verbose_flag,
                        patch,
                        n_timesteps,
                        1,
                        0,
                        0);

                patch[0].rz_storage -= rz_drainage;
                patch[0].unsat_storage += rz_drainage;

                patch[0].S =
                        min(patch[0].unsat_storage / (patch[0].sat_deficit - patch[0].rootzone.potential_sat), 1.0);
                unsat_drainage = compute_unsat_zone_drainage_patch(
                        command_line[0].verbose_flag,
                        patch,
                        n_timesteps,
                        0,
                        0,
                        0);

                patch[0].unsat_storage -= unsat_drainage;
                patch[0].sat_deficit -= unsat_drainage;
            } else {
                patch[0].sat_deficit -= patch[0].unsat_storage; /* transfer left water in unsat storage to rootzone layer */
                patch[0].unsat_storage = 0.0;

                patch[0].S =
                        min(patch[0].rz_storage / patch[0].sat_deficit, 1.0);
                rz_drainage = compute_unsat_zone_drainage_patch(
                        command_line[0].verbose_flag,
                        patch,
                        n_timesteps,
                        0,
                        0,
                        1);

                unsat_drainage = 0.0;

                patch[0].rz_storage -= rz_drainage;
                patch[0].sat_deficit -= rz_drainage;
            }

            patch[0].unsat_drainage += unsat_drainage;
            patch[0].rz_drainage += rz_drainage;

            if (patch[0].sat_deficit > patch[0].rootzone.potential_sat)
                patch[0].rootzone.S =
                        min(patch[0].rz_storage / patch[0].rootzone.potential_sat, 1.0);
            else
                patch[0].rootzone.S =
                        min((patch[0].rz_storage + patch[0].rootzone.potential_sat - patch[0].sat_deficit)
                                / patch[0].rootzone.potential_sat, 1.0);

            /*-------------------c------------------------------------------------------*/
            /*	Recompute current actual depth to water table				*/
            /*-------------------------------------------------------------------------*/
            patch[0].sat_deficit_z = compute_z_final_from_surface(patch[0].soil_defaults[0],
                    -1.0 * patch[0].sat_deficit);



            /* ******************************** this is done by each hour*/
            /* patch[0].hourly_stream_flow += patch[0].hourly_subsur2stream_flow
                          + patch[0].hourly_sur2stream_flow;*/

            //hillslope[0].hillslope_return_flow += (patch[0].return_flow) * patch[0].area;
            if (patch[0].drainage_type == STREAM
                || patch[0].innundation_list[0].num_neighbours == 0             //09132022LML add this option
               ) {
                patch[0].streamflow += patch[0].return_flow + patch[0].base_flow; //LML note: patch base_flow will be diducted from later routine
#ifdef LIU_CHECK_WATER_BALANCE
                total_Q_leave_patch += patch[0].return_flow;
#endif
                patch[0].return_flow_printout = patch[0].return_flow;
                patch[0].return_flow = 0;                                       //07182023LML (returnflow will be flowout at end of the day, treated as storage)

                //if (patch[0].ID == 64301) {
                //    printf("patch detention_store(mm):%lf streamflow(mm):%lf return_flow(mm):%lf base_flow(mm):%lf\n"
                //           ,patch[0].detention_store*1000
                //           ,patch[0].streamflow*1000
                //           ,patch[0].return_flow*1000
                //           ,patch[0].base_flow*1000);
                //}    
            }
//#ifdef LIU_OMP_PATCH_LOCK
//            omp_unset_lock(&locks_patch[0][patch[0].Unique_ID_index]);
//#endif


            /*--------------------------------------------------------------*/
            /* final stream flow calculations				*/
            /*--------------------------------------------------------------*/

            #pragma omp critical (test1)
            {
            if (patch[0].drainage_type == STREAM || patch[0].innundation_list[0].num_neighbours == 0)
              hillslope[0].hillslope_return_flow += (patch[0].return_flow) * patch[0].area;
            hillslope[0].hillslope_outflow += (patch[0].streamflow) * patch[0].area;
            hillslope[0].hillslope_unsat_storage += patch[0].unsat_storage * patch[0].area;
            hillslope[0].hillslope_sat_deficit += patch[0].sat_deficit * patch[0].area;
            hillslope[0].hillslope_rz_storage += patch[0].rz_storage * patch[0].area;
            hillslope[0].hillslope_detention_store += patch[0].detention_store
                    * patch[0].area;
            }

#ifdef LIU_CHECK_WATER_BALANCE
            //07192023LML check water balance
            double post_patch_surface_water = patch[0].detention_store + patch[0].return_flow;
            double post_patch_total_deficit = patch[0].sat_deficit - (patch[0].rz_storage + patch[0].unsat_storage);
            double post_patch_total = post_patch_surface_water - post_patch_total_deficit;

            double total_fluxout = total_Q_leave_patch; //patch[0].Qout will be counted to watertable in later routine; this is surface flow out
            double water_balance = pre_patch_total - post_patch_total - total_fluxout;

            if (fabs(water_balance) >= 0.0001) {
            //if (fabs(total_fluxout) >= 0.0001) {
                printf("Error: patch outflow waterblance from compute_subsurface_routing:\n");
                printf(" waterbalance(mm):%.1f\n",water_balance*1000);
                printf("   surface_water change(%.1f): pre:%.1f post:%.1f\n"
                       ,(post_patch_surface_water-pre_patch_surface_water)*1000
                       ,pre_patch_surface_water*1000
                       ,post_patch_surface_water*1000);
                printf("   total_deficit change(%.1f): pre:%.1f post:%.1f\n"
                       ,(post_patch_total_deficit-pre_patch_total_deficit)*1000
                       ,pre_patch_total_deficit*1000
                       ,post_patch_total_deficit*1000);
                printf("   patch_total change(%.1f): pre:%.1f post:%.1f\n"
                       ,(post_patch_total-pre_patch_total)*1000
                       ,pre_patch_total*1000
                       ,post_patch_total*1000);
                printf("   patch fluxout: %.1f\n",total_fluxout*1000);
            }
#endif

            /*---------------------------------------------------------------------*/
            /*update accumulator variables                                            */
            /*-----------------------------------------------------------------------*/
            /* the accumulator is updated in update_hillslope_patch_accumulator.c in hillslope_daily_F.c*/
            } //end i
        } //last timestep


	} /* end k  */

	hillslope[0].hillslope_outflow /= hillslope_area;
	hillslope[0].preday_hillslope_rz_storage /= hillslope_area;
	hillslope[0].preday_hillslope_unsat_storage /= hillslope_area;
	hillslope[0].preday_hillslope_detention_store /= hillslope_area;
	hillslope[0].preday_hillslope_sat_deficit /= hillslope_area;
    hillslope[0].preday_hillslope_return_flow /= hillslope_area;                //07182023LML
	hillslope[0].hillslope_rz_storage /= hillslope_area;
	hillslope[0].hillslope_unsat_storage /= hillslope_area;
	hillslope[0].hillslope_detention_store /= hillslope_area;
	hillslope[0].hillslope_sat_deficit /= hillslope_area;
    hillslope[0].hillslope_return_flow /= hillslope_area;                       //07182023LML
	water_balance = hillslope[0].preday_hillslope_rz_storage + hillslope[0].preday_hillslope_unsat_storage
			+ hillslope[0].preday_hillslope_detention_store - hillslope[0].preday_hillslope_sat_deficit
            + hillslope[0].preday_hillslope_return_flow //07182023LML
			- (hillslope[0].hillslope_rz_storage + hillslope[0].hillslope_unsat_storage + hillslope[0].hillslope_detention_store
                    + hillslope[0].hillslope_return_flow //07182023LML
                    - hillslope[0].hillslope_sat_deficit)
            - hillslope[0].hillslope_outflow;

    //07202023LML this balance is not correct since some pataches may not have
    //been processed because of inflow from upstream patches. The overall method
    //may need to be rewrote
    //if (fabs(water_balance) >= 0.0001){
    if (FALSE) {
            printf("\n Water Balance(mm) is %12.8f on %ld %ld %ld for hillslope %d",
                water_balance*1000,
                current_date.day,
                current_date.month,
                current_date.year,
                hillslope[0].ID);
            printf("\n\tdelta_rz_storage:%lf (current:%lf preday:%lf)",(hillslope[0].hillslope_rz_storage - hillslope[0].preday_hillslope_rz_storage)*1000,hillslope[0].hillslope_rz_storage*1000,hillslope[0].preday_hillslope_rz_storage*1000);
            printf("\n\tdelta_unsat_storage:%lf (current:%lf preday:%lf)",(hillslope[0].hillslope_unsat_storage - hillslope[0].preday_hillslope_unsat_storage)*1000,hillslope[0].hillslope_unsat_storage*1000,hillslope[0].preday_hillslope_unsat_storage*1000);
            printf("\n\tdelta_detention_store:%lf (current:%lf preday:%lf)",(hillslope[0].hillslope_detention_store - hillslope[0].preday_hillslope_detention_store)*1000,hillslope[0].hillslope_detention_store*1000,hillslope[0].preday_hillslope_detention_store*1000);
            printf("\n\tdelta_sat_deficit:%lf (current:%lf preday:%lf)",(hillslope[0].hillslope_sat_deficit - hillslope[0].preday_hillslope_sat_deficit)*1000,hillslope[0].hillslope_sat_deficit*1000,hillslope[0].preday_hillslope_sat_deficit*1000);
            printf("\n\tdelta_returnflow:%lf (current:%lf preday:%lf)",(hillslope[0].hillslope_return_flow - hillslope[0].preday_hillslope_return_flow)*1000,hillslope[0].hillslope_return_flow*1000,hillslope[0].preday_hillslope_return_flow*1000);
            printf("\n\thillslope_outflow:%lf",hillslope[0].hillslope_outflow*1000);
    }


	if((command_line[0].output_flags.yearly == 1)
			&& (command_line[0].b != NULL )) {
		if (hillslope_outflow <= command_line[0].thresholds[STREAMFLOW])
			hillslope[0].acc_year.num_threshold += 1;
	}

	return;

} /*end compute_subsurface_routing.c*/

