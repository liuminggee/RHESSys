/*--------------------------------------------------------------*/
/* 											*/
/*					compute_subsurface_routing_hourly		*/
/*											*/
/*	compute_subsurface_routing_hourly.c - do subsurface computation during the end of each hour	*/
/*											*/
/*	NAME										*/
/*	compute_subsurface_routing_hourly.c - do subsurface computation during the end of each hour	*/
/*											*/
/*	SYNOPSIS									*/
/*	void compute_subsurface_routing_hourly( 								*/
/*						struct command_line_object command_line, */
/*							struct hillslope_object *hillslope)	*/
/*				 			int,			 	*/
/*							struct date *current_date)	*/
/*											*/
/* 											*/
/*											*/
/*	OPTIONS										*/
/*											*/
/*											*/
/*	DESCRIPTION									*/
/*	this function is called in hillslope_hourly_test at the end of each hour during	*/
/* 	hourly calculation, it is doing the compute_subsurface_routing for each hour	*/
/*											*/
/*											*/
/*											*/
/*	PROGRAMMER NOTES								*/
/*											*/
/*											*/
/*											*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "params.h"
#include "rhessys.h"
#include "functions.h"

void compute_subsurface_routing_hourly(
		struct command_line_object *command_line,
		struct hillslope_object *hillslope, 
		int n_timesteps, 
		struct date current_date) {
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

    //double compute_layer_field_capacity(int, int, double, double, double,
    //		double, double, double, double, double, double);

	double compute_unsat_zone_drainage(int, int, double, double, double, double,
			double, double);

	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
	int i, d;
	int j, k;
	int grow_flag, verbose_flag;
	double time_int, tmp;
	double theta, m, Ksat, Nout;
	double NO3_out, NH4_out, DON_out, DOC_out;
	double return_flow, excess;
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
	double streamflow, Qout, Qin_total, Qstr_total;
	struct patch_object *patch;
	struct patch_object *neigh;
    //struct litter_object *litter;
	d=0;
	/*--------------------------------------------------------------*/
	/*	initializations						*/
	/*--------------------------------------------------------------*/
	
		grow_flag = command_line[0].grow_flag;
		verbose_flag = command_line[0].verbose_flag;

		time_int = 1.0 / n_timesteps;

	if (current_date.hour==1)
	{
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
        hillslope_area = 0.0;

		hillslope[0].hillslope_outflow = 0.0;
		hillslope[0].hillslope_area = 0.0;
		hillslope[0].hillslope_unsat_storage = 0.0;
		hillslope[0].hillslope_rz_storage = 0.0;
		hillslope[0].hillslope_sat_deficit = 0.0;
		hillslope[0].hillslope_return_flow = 0.0;
		hillslope[0].hillslope_detention_store = 0.0;
        //hillslope[0].preday_hillslope_rz_storage = 0.0;
        //hillslope[0].preday_hillslope_unsat_storage = 0.0;
        //hillslope[0].preday_hillslope_sat_deficit = 0.0;
        //hillslope[0].preday_hillslope_return_flow = 0.0;
        //hillslope[0].preday_hillslope_detention_store = 0.0;
		streamflow = 0.0;
		Qin_total = 0.0;
		Qstr_total = 0.0;
		d = 0;
		// Note: this assumes that the set of patches in the surface routing table is identical to
		//       the set of patches in the subsurface flow table
        //#pragma omp parallel for reduction(+ : preday_hillslope_rz_storage,preday_hillslope_unsat_storage,preday_hillslope_sat_deficit,preday_hillslope_return_flow,preday_hillslope_detention_store,hillslope_area)
        for (int i = 0; i < hillslope->route_list->num_patches; i++) {
            struct patch_object *patch = hillslope->route_list->list[i];

			patch[0].streamflow = 0.0;
			patch[0].return_flow = 0.0;
			patch[0].base_flow = 0.0;
			patch[0].infiltration_excess = 0.0;
            //hillslope[0].preday_hillslope_rz_storage += patch[0].rz_storage * patch[0].area;
            //hillslope[0].preday_hillslope_unsat_storage += patch[0].unsat_storage * patch[0].area;
            //hillslope[0].preday_hillslope_sat_deficit += patch[0].sat_deficit * patch[0].area;
            //hillslope[0].preday_hillslope_return_flow += patch[0].return_flow * patch[0].area;
            //hillslope[0].preday_hillslope_detention_store += patch[0].detention_store
            //		* patch[0].area;
            //hillslope[0].hillslope_area += patch[0].area;
            preday_hillslope_rz_storage += patch[0].rz_storage * patch[0].area;
            preday_hillslope_unsat_storage += patch[0].unsat_storage * patch[0].area;
            preday_hillslope_sat_deficit += patch[0].sat_deficit * patch[0].area;
            preday_hillslope_return_flow += patch[0].return_flow * patch[0].area;
            preday_hillslope_detention_store += patch[0].detention_store
                    * patch[0].area;
            hillslope_area += patch[0].area;
			patch[0].Qin_total = 0.0;
			patch[0].Qout_total = 0.0;
			patch[0].Qin = 0.0;
			patch[0].Qout = 0.0;
			patch[0].surface_Qin = 0.0;
			patch[0].surface_Qout = 0.0;

			patch[0].overland_flow = 0.0;


			patch[0].interim_sat = patch[0].sat_deficit - patch[0].unsat_storage;
			if ((patch[0].sat_deficit - patch[0].unsat_storage) < ZERO)
				patch[0].S = 1.0;
			else
				patch[0].S = patch[0].unsat_storage / patch[0].sat_deficit;

			if (grow_flag > 0) {
				patch[0].soil_ns.NO3_Qin = 0.0;
				patch[0].soil_ns.NO3_Qout = 0.0;
				patch[0].soil_ns.NH4_Qin = 0.0;
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
				patch[0].soil_ns.DON_Qin = 0.0;
				patch[0].soil_cs.DOC_Qout = 0.0;
				patch[0].soil_cs.DOC_Qin = 0.0;
				patch[0].surface_DON_Qout = 0.0;
				patch[0].surface_DON_Qin = 0.0;
				patch[0].surface_DOC_Qout = 0.0;
				patch[0].surface_DOC_Qin = 0.0;
			
				patch[0].streamNO3_from_surface	= 0.0;
				patch[0].streamNO3_from_sub = 0.0;
			}
		}
        hillslope[0].preday_hillslope_rz_storage = preday_hillslope_rz_storage;
        hillslope[0].preday_hillslope_unsat_storage = preday_hillslope_unsat_storage;
        hillslope[0].preday_hillslope_sat_deficit = preday_hillslope_sat_deficit;
        hillslope[0].preday_hillslope_return_flow = preday_hillslope_return_flow;
        hillslope[0].preday_hillslope_detention_store = preday_hillslope_detention_store;
        hillslope[0].hillslope_area = hillslope_area;
	}
	/*--------------------------------------------------------------*/
	/*	calculate Qout for each patch and add appropriate	*/
	/*	proportion of subsurface outflow to each neighbour	*/
	/*--------------------------------------------------------------*/
        //#pragma omp parallel for
        for (int i = 0; i < hillslope->route_list->num_patches; i++) {
            struct patch_object *patch = hillslope->route_list->list[i];
						
			patch[0].preday_sat_deficit = patch[0].sat_deficit;
            patch[0].preday_sat_deficit_z = compute_z_final_from_surface(patch[0].soil_defaults[0],
					-1.0 * patch[0].sat_deficit);
			
            patch[0].hourly_subsur2stream_flow = 0;
			patch[0].hourly_sur2stream_flow = 0;
			patch[0].hourly_stream_flow = 0;
			patch[0].hourly[0].streamflow_NO3 = 0;
			patch[0].hourly[0].streamflow_NO3_from_sub = 0;
			patch[0].hourly[0].streamflow_NO3_from_surface = 0;
		}
        //10172022LML seems no benifit #pragma omp parallel for
        for (int i = 0; i < hillslope->route_list->num_patches; i++) {
            struct patch_object *patch = hillslope->route_list->list[i];
//#ifdef LIU_OMP_PATCH_LOCK
//            omp_set_lock(&locks_patch[patch[0].Unique_ID_index]);
//#endif
            struct litter_object *litter=&(patch[0].litter);
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
//            omp_unset_lock(&locks_patch[patch[0].Unique_ID_index]);
//#endif
		} /* end i */

		/*--------------------------------------------------------------*/
		/*	update soil moisture and nitrogen stores		*/
		/*	check water balance					*/
		/*--------------------------------------------------------------*/
        //10172022LML seems no benifit #pragma omp parallel for private(excess,rz_drainage,add_field_capacity,infiltration,unsat_drainage)
        for (int i = 0; i < hillslope->route_list->num_patches; i++) {
            struct patch_object *patch = hillslope->route_list->list[i];
//#ifdef LIU_OMP_PATCH_LOCK
//            omp_set_lock(&locks_patch[1][patch[1].Unique_ID_index]);
//#endif

			/*--------------------------------------------------------------*/
			/*	update subsurface 				*/
			/*-------------------------------------------------------------------------*/

			/*-------------------------------------------------------------------------*/
			/*	Recompute current actual depth to water table				*/
			/*-------------------------------------------------------------------------*/
			patch[0].sat_deficit += (patch[0].Qout - patch[0].Qin); // this part need to put into some where else

            patch[0].sat_deficit_z = compute_z_final_from_surface(patch[0].soil_defaults[0],
					-1.0 * patch[0].sat_deficit);

			if (grow_flag > 0) {
				patch[0].soil_ns.nitrate += (patch[0].soil_ns.NO3_Qin
						- patch[0].soil_ns.NO3_Qout);

                printf("3 nitrate:%lf\n",patch[0].soil_ns.nitrate);

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

				
					if (grow_flag > 0) {
                        double lpools[] = {patch[0].soil_ns.nitrate,
                                          patch[0].soil_ns.sminn,
                                          patch[0].soil_ns.DON,
                                          patch[0].soil_cs.DOC
                                          };
                        double leached[LEACH_ELEMENT_counts];
                        double t =
                                compute_N_leached_from_soildef(verbose_flag,
                                        lpools, excess, 0.0, 0.0,
                                        //patch[0].m,
                                        //patch[0].innundation_list[d].gamma
                                        //		/ patch[0].area * time_int,
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
                            double Qout = excess * patch->surface_innundation_list[d].neighbours[j].gamma;
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
//                            omp_set_lock(&locks_patch[1][neigh[1].Unique_ID_index]);
//#endif
							if (neigh[0].drainage_type == STREAM) {
								neigh[0].Qin_total += Qout * patch[0].area
										/ neigh[0].area;
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
                                //printf("detention_store:%lf Qout_2:%lf \n",
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
//                            omp_unset_lock(&locks_patch[1][neigh[1].Unique_ID_index]);
//#endif
                        } //j
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
                        } //if
						patch[0].detention_store -= excess;
						patch[0].Qout_total += excess;
                    } //else
                } //if
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
					add_field_capacity = compute_layer_field_capacity(
							command_line[0].verbose_flag,
							patch[0].soil_defaults[0][0].theta_psi_curve,
							patch[0].soil_defaults[0][0].psi_air_entry,
							patch[0].soil_defaults[0][0].pore_size_index,
							patch[0].soil_defaults[0][0].p3,
							patch[0].soil_defaults[0][0].p4,
							patch[0].soil_defaults[0][0].porosity_0,
							patch[0].soil_defaults[0][0].porosity_decay,
                            patch[0].soil_defaults[0][0].Dingman_coef,
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
						add_field_capacity = compute_layer_field_capacity(
								command_line[0].verbose_flag,
								patch[0].soil_defaults[0][0].theta_psi_curve,
								patch[0].soil_defaults[0][0].psi_air_entry,
								patch[0].soil_defaults[0][0].pore_size_index,
								patch[0].soil_defaults[0][0].p3,
								patch[0].soil_defaults[0][0].p4,
								patch[0].soil_defaults[0][0].porosity_0,
								patch[0].soil_defaults[0][0].porosity_decay,
                                patch[0].soil_defaults[0][0].Dingman_coef,
								patch[0].sat_deficit_z, patch[0].sat_deficit_z,
								0.0);
						add_field_capacity = max(add_field_capacity, 0.0);
						patch[0].sat_deficit += add_field_capacity;
						patch[0].rz_storage += add_field_capacity;
					}
				} else {
					if ((patch[0].sat_deficit > ZERO)
                            && close_enough(patch[0].unsat_storage, 0.0)) {
						add_field_capacity = compute_layer_field_capacity(
								command_line[0].verbose_flag,
								patch[0].soil_defaults[0][0].theta_psi_curve,
								patch[0].soil_defaults[0][0].psi_air_entry,
								patch[0].soil_defaults[0][0].pore_size_index,
								patch[0].soil_defaults[0][0].p3,
								patch[0].soil_defaults[0][0].p4,
								patch[0].soil_defaults[0][0].porosity_0,
								patch[0].soil_defaults[0][0].porosity_decay,
                                patch[0].soil_defaults[0][0].Dingman_coef,
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

                    printf("2 nitrate:%lf\n",patch[0].soil_ns.nitrate);

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
                patch[0].sat_deficit_z = compute_z_final_from_surface(
                        patch[0].soil_defaults[0],
						-1.0 * patch[0].sat_deficit);


			
				/*--------------------------------------------------------------*/
				/*	compute new field capacity				*/
				/*--------------------------------------------------------------*/
				if (patch[0].sat_deficit_z < patch[0].rootzone.depth) {
					patch[0].rootzone.field_capacity =
							compute_layer_field_capacity(
									command_line[0].verbose_flag,
									patch[0].soil_defaults[0][0].theta_psi_curve,
									patch[0].soil_defaults[0][0].psi_air_entry,
									patch[0].soil_defaults[0][0].pore_size_index,
									patch[0].soil_defaults[0][0].p3,
									patch[0].soil_defaults[0][0].p4,
									patch[0].soil_defaults[0][0].porosity_0,
									patch[0].soil_defaults[0][0].porosity_decay,
                                    patch[0].soil_defaults[0][0].Dingman_coef,
									patch[0].sat_deficit_z,
									patch[0].rootzone.depth, 0.0);

					patch[0].field_capacity = 0.0;
				} else {

					patch[0].rootzone.field_capacity =
							compute_layer_field_capacity(
									command_line[0].verbose_flag,
									patch[0].soil_defaults[0][0].theta_psi_curve,
									patch[0].soil_defaults[0][0].psi_air_entry,
									patch[0].soil_defaults[0][0].pore_size_index,
									patch[0].soil_defaults[0][0].p3,
									patch[0].soil_defaults[0][0].p4,
									patch[0].soil_defaults[0][0].porosity_0,
									patch[0].soil_defaults[0][0].porosity_decay,
                                    patch[0].soil_defaults[0][0].Dingman_coef,
									patch[0].sat_deficit_z,
									patch[0].rootzone.depth, 0.0);

					patch[0].field_capacity = compute_layer_field_capacity(
							command_line[0].verbose_flag,
							patch[0].soil_defaults[0][0].theta_psi_curve,
							patch[0].soil_defaults[0][0].psi_air_entry,
							patch[0].soil_defaults[0][0].pore_size_index,
							patch[0].soil_defaults[0][0].p3,
							patch[0].soil_defaults[0][0].p4,
							patch[0].soil_defaults[0][0].porosity_0,
							patch[0].soil_defaults[0][0].porosity_decay,
                            patch[0].soil_defaults[0][0].Dingman_coef,
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
                patch[0].sat_deficit_z = compute_z_final_from_surface(
                        patch[0].soil_defaults[0],
						-1.0 * patch[0].sat_deficit);



			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/
			/*--------------------------------------------------------------*/			
			patch[0].hourly_stream_flow += patch[0].hourly_subsur2stream_flow
		      				+ patch[0].hourly_sur2stream_flow;

            #pragma omp critical (test2)
            {
			hillslope[0].hillslope_return_flow += (patch[0].return_flow) * patch[0].area;
            }



			/* ******************************** */
			/* accumulate the daily returnflow and baseflow calculated from update_drainage*/
			/* The N calculation has been completed in update_drainage_***.c routing*/
			if (current_date.hour == n_timesteps){
				    if (patch[0].drainage_type == STREAM) {
					patch[0].streamflow += patch[0].return_flow
							+ patch[0].base_flow;
				    }
			}
		    

//#ifdef LIU_OMP_PATCH_LOCK
//            omp_unset_lock(&locks_patch[1][patch[1].Unique_ID_index]);
//#endif
		} /* end i */


	hillslope[0].hillslope_outflow /= hillslope[0].hillslope_area;
	hillslope[0].preday_hillslope_rz_storage /= hillslope[0].hillslope_area;
	hillslope[0].preday_hillslope_unsat_storage /= hillslope[0].hillslope_area;
	hillslope[0].preday_hillslope_detention_store /= hillslope[0].hillslope_area;
	hillslope[0].preday_hillslope_sat_deficit /= hillslope[0].hillslope_area;
	hillslope[0].hillslope_rz_storage /= hillslope[0].hillslope_area;
	hillslope[0].hillslope_unsat_storage /= hillslope[0].hillslope_area;
	hillslope[0].hillslope_detention_store /= hillslope[0].hillslope_area;
	hillslope[0].hillslope_sat_deficit /= hillslope[0].hillslope_area;
	water_balance = hillslope[0].preday_hillslope_rz_storage + hillslope[0].preday_hillslope_unsat_storage
			+ hillslope[0].preday_hillslope_detention_store - hillslope[0].preday_hillslope_sat_deficit
			- (hillslope[0].hillslope_rz_storage + hillslope[0].hillslope_unsat_storage + hillslope[0].hillslope_detention_store
					- hillslope[0].hillslope_sat_deficit) - hillslope[0].hillslope_outflow;

	if ((command_line[0].output_flags.yearly == 1)
			&& (command_line[0].b != NULL )) {
		if (hillslope[0].hillslope_outflow <= command_line[0].thresholds[STREAMFLOW])
			hillslope[0].acc_year.num_threshold += 1;
	}

	return;

} /*end compute_subsurface_routing_hourly.c*/


