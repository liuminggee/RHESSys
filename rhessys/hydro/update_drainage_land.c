/*--------------------------------------------------------------*/
/* 											*/
/*					update_drainage_land			*/
/*											*/
/*	update_drainage_land.c - creates a patch object				*/
/*											*/
/*	NAME										*/
/*	update_drainage_land.c - creates a patch object				*/
/*											*/
/*	SYNOPSIS									*/
/*	void update_drainage_land( 							*/
/*					struct patch_object *patch			*/
/*				 			double,			 	*/
/*				 			double,			 	*/
/*				 			double,			 	*/
/*							int,				*/
/*							int)				*/
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
/*											*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"
#include "functions.h"
#ifdef LIU_OMP_PATCH_LOCK
#include "params.h"
#endif

void  update_drainage_land(
					struct patch_object *patch,
					 struct command_line_object *command_line,
					 double time_int,
					 int verbose_flag)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	double  compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);
	
	double compute_varbased_returnflow(
		double,
		double,
		double,
		struct litter_object *);


	double compute_varbased_flow(
		int,
		double,
		double,
		double,
		double,
		double *,
		struct patch_object *);


//	double compute_N_leached(int,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double,
//		double *);
	
	double recompute_gamma(	
		struct patch_object *,
		double);


	double compute_infiltration( int,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double,
		double);
	
	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
    int d; //, idx;
    //double Ksat;
    double std_scale;
    double NH4_leached_to_patch;//, NH4_leached_to_stream;
    double NO3_leached_to_patch;//, NO3_leached_to_stream;
    double DON_leached_to_patch;//, DON_leached_to_stream;
    double DOC_leached_to_patch;//, DOC_leached_to_stream;
	double NO3_leached_to_surface; /* kg/m2 */
	double NH4_leached_to_surface; /* kg/m2 */
	double DON_leached_to_surface; /* kg/m2 */
	double DOC_leached_to_surface; /* kg/m2 */
    //double N_leached_total; /* kg/m2 */
    //double DON_leached_total; /* kg/m2 */
    //double DOC_leached_total; /* kg/m2 */
	double route_to_surface;  /* m3 */
	double return_flow,route_to_patch ;  /* m3 */
	double available_sat_water; /* m3 */
    double Qout = 0;  /* m */
    double innundation_depth; /* m */
	double total_gamma;
    double Nout; /* kg/m2 */
    //double t1,t2,t3;

    //struct patch_object *neigh;
	route_to_patch = 0.0;
	route_to_surface = 0.0;
	return_flow=0.0;

	DON_leached_to_patch = 0.0;
    //DON_leached_to_stream = 0.0;
	DOC_leached_to_patch = 0.0;
    //DOC_leached_to_stream = 0.0;
	NH4_leached_to_patch = 0.0;
    //NH4_leached_to_stream = 0.0;
	NO3_leached_to_patch = 0.0;
    //NO3_leached_to_stream = 0.0;
	NO3_leached_to_surface = 0.0;
	NH4_leached_to_surface = 0.0;
	DOC_leached_to_surface = 0.0;
	DON_leached_to_surface = 0.0;
    struct  soil_default *psoil_def = patch[0].soil_defaults[0];
#ifdef LIU_CHECK_WATER_BALANCE
    //07192023LML check waterbalance
    double pre_patch_surface_water = patch[0].detention_store + patch[0].return_flow;
    double pre_patch_total_deficit = patch[0].sat_deficit - (patch[0].rz_storage + patch[0].unsat_storage);
    double pre_patch_total = pre_patch_surface_water - pre_patch_total_deficit;
#endif
	/*--------------------------------------------------------------*/
	/*	m and K are multiplied by sensitivity analysis variables */
	/*--------------------------------------------------------------*/

    //m = patch[0].m ;
    //Ksat = patch[0].soil_defaults[0][0].Ksat_0 ;
	d=0;

	/*--------------------------------------------------------------*/
	/*	recalculate gamma based on current saturation deficits  */
	/*      to account the effect of changes in water table slope 	*/
	/*--------------------------------------------------------------*/

    total_gamma = recompute_gamma(patch,
                                  patch[0].innundation_list[d].gamma
                                  );

    //05312023LML note: it is total soil water content in saturated zone (m3)
    available_sat_water = max(((psoil_def[0].soil_water_cap
			- max(patch[0].sat_deficit,0.0))
			* patch[0].area),0.0);

	/*------------------------------------------------------------*/
	/*	calculate amuount of water output to patches			*/
	/*	this only computes subsurface flow, not overland flow	*/
	/*-----------------------------------------------------------*/

	std_scale = command_line[0].std_scale;

    //05312023LML note: it is subflow from this patch
	route_to_patch =  time_int * compute_varbased_flow(
		patch[0].num_soil_intervals,
		patch[0].std * std_scale, 
		patch[0].sat_deficit,
		total_gamma, 
        psoil_def[0].interval_size,
		patch[0].transmissivity_profile,
		patch);


    //printf("route_to_patch(mm/day):%lf total_gamma:%lf\n",route_to_patch*1000/time_int/patch[0].area,total_gamma);
    //double gamma = total_gamma / patch[0].area * time_int;


	if (route_to_patch < 0.0) route_to_patch = 0.0;
	if ( route_to_patch > available_sat_water) 
        route_to_patch *= available_sat_water/route_to_patch;
    double pQout = route_to_patch / patch[0].area;

	/*--------------------------------------------------------------*/
	/* compute Nitrogen leaching amount				*/
	/*--------------------------------------------------------------*/
    double lpools[] = {patch[0].soil_ns.nitrate,
                      patch[0].soil_ns.sminn,
                      patch[0].soil_ns.DON,
                      patch[0].soil_cs.DOC
                      };
    double leached[LEACH_ELEMENT_counts];

//#ifdef LIU_OMP_PATCH_LOCK
//        omp_set_lock(&locks_patch[0][patch[0].Unique_ID_index]);
//#endif

	if (command_line[0].grow_flag > 0) {
        //double *Nout =
        double t = compute_N_leached_from_soildef(
			verbose_flag,
            lpools,
            pQout,
			patch[0].sat_deficit,
            psoil_def[0].soil_water_cap,
            //m,
            //gamma,
            psoil_def,
            leached,
            LEACH_ELEMENT_counts
            //patch[0].transmissivity_profile
            );

        NO3_leached_to_patch = leached[LNO3] * patch[0].area;
        patch[0].soil_ns.NO3_Qout += leached[LNO3];

        NH4_leached_to_patch = leached[LNH4] * patch[0].area;
        patch[0].soil_ns.NH4_Qout += leached[LNH4];

        DON_leached_to_patch = leached[LDON] * patch[0].area;
        patch[0].soil_ns.DON_Qout += leached[LDON];

        DOC_leached_to_patch = leached[LDOC] * patch[0].area;
        patch[0].soil_cs.DOC_Qout += leached[LDOC];
        //free(Nout);

	}

	
    patch[0].Qout += pQout;                                                     //Will be counted to update water table in later routine
    //09132022LML
    if (patch[0].innundation_list[0].num_neighbours == 0) {
        patch[0].base_flow += pQout;
    }


	/*--------------------------------------------------------------*/
	/*	calculate any return flow associated with this patch	*/
	/*	and route any infiltration excess			*/
	/*	return flow is flow leaving patch (i.e surface_Qout)  	*/
	/*	note that return flow that becomes detention storage   */
	/*	is added to surface_Qin					*/
	/*	similarly with associated nitrogen			*/
	/* 	note we move unsat_storage into saturated storage in this case */
	/*	saturated zone will be updated in compute_subsurface_routing	*/
	/*	i.e becomes part of Qout				*/
	/*--------------------------------------------------------------*/
	if ((patch[0].sat_deficit-patch[0].rz_storage-patch[0].unsat_storage) < -1.0*ZERO) {
		return_flow = compute_varbased_returnflow(patch[0].std * std_scale, 
			patch[0].rz_storage+patch[0].unsat_storage,
			patch[0].sat_deficit, &(patch[0].litter));


        //fprintf(stderr,"return_flow(mm):%lf\n",return_flow*1000);


		patch[0].detention_store += return_flow;  


        //if (patch[0].ID == 81497)
        //printf("detention_store:%lf return_flow:%lf \n",
        //       patch[0].detention_store*1000,
        //       patch[0].return_flow*1000);

		patch[0].sat_deficit += (return_flow - (patch[0].unsat_storage+patch[0].rz_storage));
		patch[0].unsat_storage = 0.0;
		patch[0].rz_storage = 0.0;
	}
	/*--------------------------------------------------------------*/
	/*	calculated any N-transport associated with return flow  */
	/*	-note available N reduced by what has already been 	*/
	/*	we assume that only nitrate follows return flow		*/
	/*	lost in subsurface flow routing				*/
	/*--------------------------------------------------------------*/
    if (command_line[0].grow_flag > 0) {
            for (int i = 0; i < LEACH_ELEMENT_counts; i++) {
                switch(i) {
                case LNO3:
                    lpools[i] = patch[0].soil_ns.nitrate - NO3_leached_to_patch/patch[0].area;
                    break;
                case LNH4:
                    lpools[i] = patch[0].soil_ns.sminn - NH4_leached_to_patch/patch[0].area;
                    break;
                case LDON:
                    lpools[i] = patch[0].soil_ns.DON - DON_leached_to_patch/patch[0].area;
                    break;
                case LDOC:
                    lpools[i] = patch[0].soil_cs.DOC - DOC_leached_to_patch/patch[0].area;
                }
            }
            //double *pNout =
            double t = compute_N_leached_from_soildef(
				verbose_flag,
                lpools,
				return_flow,
				0.0,
				0.0,
                //m,
                //gamma,
                psoil_def,
                leached,
                LEACH_ELEMENT_counts
                //patch[0].transmissivity_profile
                    );

            patch[0].surface_NO3 += leached[LNO3];
            patch[0].soil_ns.NO3_Qout += leached[LNO3];

            patch[0].surface_NH4 += leached[LNH4];
            patch[0].soil_ns.NH4_Qout += leached[LNH4];

            patch[0].surface_DON += leached[LDON];
            patch[0].soil_ns.DON_Qout += leached[LDON];

            patch[0].surface_DOC += leached[LDOC];
            patch[0].soil_cs.DOC_Qout += leached[LDOC];
            //free(pNout);
    }
	
	/*--------------------------------------------------------------*/
	/*	route water and nitrogen lossed due to infiltration excess */
	/*--------------------------------------------------------------*/
    if ( (patch[0].detention_store > psoil_def[0].detention_store_size) &&
		(patch[0].detention_store > ZERO) ){

        Qout = patch[0].detention_store - psoil_def[0].detention_store_size;
        patch[0].overland_flow += Qout;

        //printf("Qout(mm):%.1f detention_store:%.1f\n",Qout*1000,patch[0].detention_store*1000);

        double ffrac = min(1.0, (Qout/ patch[0].detention_store));
		if (command_line[0].grow_flag > 0) {
            Nout = ffrac * patch[0].surface_DOC;
			DOC_leached_to_surface = Nout * patch[0].area;
			patch[0].surface_DOC -= Nout;
            Nout = ffrac * patch[0].surface_DON;
			DON_leached_to_surface = Nout * patch[0].area;
			patch[0].surface_DON -= Nout;
            Nout = ffrac * patch[0].surface_NO3;
			NO3_leached_to_surface = Nout * patch[0].area;
			patch[0].surface_NO3 -= Nout;
            Nout = ffrac * patch[0].surface_NH4;
			NH4_leached_to_surface = Nout * patch[0].area;
			patch[0].surface_NH4 -= Nout;
			}
        route_to_surface = Qout *  patch[0].area;
		patch[0].detention_store -= Qout;
		patch[0].surface_Qout += Qout;

        //09132022LML
        if (patch[0].innundation_list[0].num_neighbours == 0) {
            patch[0].return_flow += Qout;
        }

    }//if

#ifdef LIU_CHECK_WATER_BALANCE
    //07192023LML check water balance
    double post_patch_surface_water = patch[0].detention_store + patch[0].return_flow;
    double post_patch_total_deficit = patch[0].sat_deficit - (patch[0].rz_storage + patch[0].unsat_storage);
    double post_patch_total = post_patch_surface_water - post_patch_total_deficit;

    double total_fluxout = Qout; //patch[0].Qout will be counted to watertable in later routine; this is surface flow out
    double water_balance = pre_patch_total - post_patch_total - total_fluxout;

    if (fabs(water_balance) >= 0.0001) {
        printf("Error: patch outflow waterblance from update_drainage_land:\n");
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

//#ifdef LIU_OMP_PATCH_LOCK
//        omp_unset_lock(&locks_patch[0][patch[0].Unique_ID_index]);
//#endif


	if (NO3_leached_to_surface < 0.0)
		printf("WARNING %d %lf",patch[0].ID, NO3_leached_to_surface);

	/*--------------------------------------------------------------*/
	/*	route flow to neighbours				*/
	/*	route n_leaching if grow flag specfied			*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/* regular downslope routing */
	/*--------------------------------------------------------------*/
	if (command_line[0].noredist_flag == 0) {
      int numnb = patch[0].innundation_list[d].num_neighbours;
    //#pragma omp parallel for
      for (int j = 0; j < numnb; j++) {

        //07202023LML this routine seems problematic since the flow to downstream neighbors may not be able to be processed after the time loop but the flux has been cleared out at the beging of day

        struct  neighbour_object *nbo = &patch[0].surface_innundation_list[d].neighbours[j];
        struct  patch_object *neigh = nbo->patch;
//#ifdef LIU_OMP_PATCH_LOCK
//        omp_set_lock(&locks_patch[0][neigh[0].Unique_ID_index]);
//#endif
        double fgamma = nbo->gamma / neigh[0].area;
		/*--------------------------------------------------------------*/
		/* first transfer subsurface water and nitrogen */
		/*--------------------------------------------------------------*/
        double Qin =	fgamma * route_to_patch;
		if (Qin < 0) printf("\n warning negative routing from patch %d with gamma %lf", patch[0].ID, total_gamma);
        //06222023LML added surface waterbody option
        if (neigh[0].IsWaterBody) {
          if (command_line[0].grow_flag > 0) {
            neigh[0].surface_DON += fgamma * DON_leached_to_patch;
            neigh[0].surface_DOC += fgamma * DOC_leached_to_patch;
            neigh[0].surface_NO3 += fgamma * NO3_leached_to_patch;
            neigh[0].surface_NH4 += fgamma * NH4_leached_to_patch;
          }
          neigh[0].detention_store += Qin;
        } else {
          if (command_line[0].grow_flag > 0) {
            neigh[0].soil_ns.DON_Qin += fgamma * DON_leached_to_patch;
            neigh[0].soil_cs.DOC_Qin += fgamma * DOC_leached_to_patch;
            neigh[0].soil_ns.NO3_Qin += fgamma * NO3_leached_to_patch;
            neigh[0].soil_ns.NH4_Qin += fgamma * NH4_leached_to_patch;
          }
          neigh[0].Qin += Qin;
        }
//#ifdef LIU_OMP_PATCH_LOCK
//        omp_unset_lock(&locks_patch[0][neigh[0].Unique_ID_index]);
//#endif
      } //j

	/*--------------------------------------------------------------*/
	/* surface downslope routing */
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/* determine which innundation depth to consider		*/
	/*--------------------------------------------------------------*/
      d = 0;
      if (patch[0].num_innundation_depths > 0) {
		  innundation_depth = patch[0].detention_store + route_to_surface/patch[0].area; 
		  while ((innundation_depth > patch[0].innundation_list[d].critical_depth) 
			  && (d < patch[0].num_innundation_depths-1)) {
			  d++;}
      }
    //#pragma omp parallel for
      numnb = patch[0].innundation_list[d].num_neighbours;
      //test balance
      double sum_Qin = 0;
      for (int j = 0; j < numnb; j++) {
        struct  neighbour_object *nbo = &patch[0].surface_innundation_list[d].neighbours[j];
        struct  patch_object *neigh = nbo->patch;
        struct  soil_default *neigh_psoil = neigh[0].soil_defaults[0];
//#ifdef LIU_OMP_PATCH_LOCK
//        omp_set_lock(&locks_patch[0][neigh[0].Unique_ID_index]);
//#endif
        double fgamma = nbo->gamma / neigh[0].area;
		/*--------------------------------------------------------------*/
		/* now transfer surface water and nitrogen */
		/*	- first nitrogen					*/
		/*--------------------------------------------------------------*/
		if (command_line[0].grow_flag > 0) {
            neigh[0].surface_NO3 += fgamma * NO3_leached_to_surface;
            neigh[0].surface_NH4 += fgamma * NH4_leached_to_surface;
            neigh[0].surface_DON += fgamma * DON_leached_to_surface;
            neigh[0].surface_DOC += fgamma * DOC_leached_to_surface;
        }
		
		/*--------------------------------------------------------------*/
		/*	- now surface water 					*/
		/*	surface stores should be updated to facilitate transfer */
		/* added net surface water transfer to detention store		*/
		/*--------------------------------------------------------------*/

        double Qin = fgamma * route_to_surface;
        sum_Qin += Qin * neigh[0].area;

#ifdef LIU_CHECK_WATER_BALANCE
        //07192023LML check neigh waterbalance
        double pre_neigh_surface_water = neigh[0].detention_store + neigh[0].return_flow;
        double pre_neigh_total_deficit = neigh[0].sat_deficit - (neigh[0].rz_storage + neigh[0].unsat_storage);
        double pre_neigh_total = pre_neigh_surface_water - pre_neigh_total_deficit;
#endif


		neigh[0].detention_store += Qin;// need fix this
        //if (neigh[0].ID == 81497)
        //printf("detention_store:%lf Qin:%lf \n",
        //       neigh[0].detention_store*1000,
        //       Qin*1000);
		neigh[0].surface_Qin += Qin;
        double infiltration = 0;
		/*--------------------------------------------------------------*/
		/* try to infiltrate this water					*/ 
		/* use time_int as duration */
		/*--------------------------------------------------------------*/
		if (neigh[0].detention_store > ZERO) {
			if (neigh[0].rootzone.depth > ZERO) {
            infiltration = compute_infiltration_patch(
				verbose_flag,
                neigh,
				neigh[0].rootzone.S,
				(neigh[0].detention_store),	
                time_int);
			}
			else {
            infiltration = compute_infiltration_patch(
				verbose_flag,
                neigh,
				neigh[0].S,
                (neigh[0].detention_store),
                time_int);
			}
		}
		else infiltration = 0.0;
        double in_frac = infiltration / neigh[0].detention_store;

        //printf("infiltration:%f\n",infiltration);

		/*--------------------------------------------------------------*/
		/* added an surface N flux to surface N pool	and		*/
		/* allow infiltration of surface N				*/
		/*--------------------------------------------------------------*/
		if ((command_line[0].grow_flag > 0 ) && (infiltration > ZERO)) {
            neigh[0].soil_cs.DOC_Qin += in_frac * neigh[0].surface_DOC;
            neigh[0].surface_DOC -= in_frac * neigh[0].surface_DOC;
            neigh[0].soil_ns.DON_Qin += in_frac * neigh[0].surface_DON;
            neigh[0].surface_DON -= in_frac * neigh[0].surface_DON;
            neigh[0].soil_ns.NO3_Qin += in_frac * neigh[0].surface_NO3;
            neigh[0].surface_NO3 -= in_frac * neigh[0].surface_NO3;
            neigh[0].soil_ns.NH4_Qin += in_frac * neigh[0].surface_NH4;
            neigh[0].surface_NH4 -= in_frac * neigh[0].surface_NH4;
		}

		if (infiltration > neigh[0].sat_deficit - neigh[0].unsat_storage - neigh[0].rz_storage) {
			neigh[0].sat_deficit -= (infiltration + neigh[0].unsat_storage + neigh[0].rz_storage);
			neigh[0].unsat_storage = 0.0; 
			neigh[0].rz_storage = 0.0; 
			neigh[0].field_capacity = 0.0; 
			neigh[0].rootzone.field_capacity = 0.0; 
		}

		else if ((neigh[0].sat_deficit > neigh[0].rootzone.potential_sat) &&
			(infiltration > neigh[0].rootzone.potential_sat - neigh[0].rz_storage)) {
		/*------------------------------------------------------------------------------*/
		/*		Just add the infiltration to the rz_storage and unsat_storage	*/
		/*------------------------------------------------------------------------------*/
			neigh[0].unsat_storage += infiltration - (neigh[0].rootzone.potential_sat - neigh[0].rz_storage);
			neigh[0].rz_storage = neigh[0].rootzone.potential_sat;
		}								
		/* Only rootzone layer saturated - perched water table case */
		else if ((neigh[0].sat_deficit > neigh[0].rootzone.potential_sat) &&
			(infiltration <= neigh[0].rootzone.potential_sat - neigh[0].rz_storage)) {
			/*--------------------------------------------------------------*/
			/*		Just add the infiltration to the rz_storage	*/
			/*--------------------------------------------------------------*/
			neigh[0].rz_storage += infiltration;
		}
		else if ((neigh[0].sat_deficit <= neigh[0].rootzone.potential_sat) &&
			(infiltration <= neigh[0].sat_deficit - neigh[0].rz_storage - neigh[0].unsat_storage)) {
			neigh[0].rz_storage += neigh[0].unsat_storage;		
			/* transfer left water in unsat storage to rootzone layer */
			neigh[0].unsat_storage = 0;
			neigh[0].rz_storage += infiltration;
			neigh[0].field_capacity = 0;
		}

		neigh[0].detention_store -= infiltration;
//#ifdef LIU_OMP_PATCH_LOCK
//        omp_unset_lock(&locks_patch[0][neigh[0].Unique_ID_index]);
//#endif
#ifdef LIU_CHECK_WATER_BALANCE
        //07192023LML check water balance
        double post_neigh_surface_water = neigh[0].detention_store + neigh[0].return_flow;
        double post_neigh_total_deficit = neigh[0].sat_deficit - (neigh[0].rz_storage + neigh[0].unsat_storage);
        double post_neigh_total = post_neigh_surface_water - post_neigh_total_deficit;

        double neigh_total_fluxout = -Qin;
        double neigh_water_balance = pre_neigh_total - post_neigh_total - neigh_total_fluxout;

        if (fabs(neigh_water_balance) >= 0.0001) {
            printf("Error: neigh outflow waterblance from update_drainage_land:\n");
            printf(" waterbalance(mm):%.1f\n",neigh_water_balance*1000);
            printf("   surface_water change(%.1f): pre:%.1f post:%.1f\n"
                   ,(post_neigh_surface_water-pre_neigh_surface_water)*1000
                   ,pre_neigh_surface_water*1000
                   ,post_neigh_surface_water*1000);
            printf("   total_deficit change(%.1f): pre:%.1f post:%.1f\n"
                   ,(post_neigh_total_deficit-pre_neigh_total_deficit)*1000
                   ,pre_neigh_total_deficit*1000
                   ,post_neigh_total_deficit*1000);
            printf("   neigh_total change(%.1f): pre:%.1f post:%.1f\n"
                   ,(post_neigh_total-pre_neigh_total)*1000
                   ,pre_neigh_total*1000
                   ,post_neigh_total*1000);
            printf("   neigh fluxout: %.1f\n",total_fluxout*1000);
        }
#endif
      }//j
      if (fabs(sum_Qin - route_to_surface) >= 1) {
          printf("Balance error in update_drainage_land:\n");
          printf(" sum_Qin:%.1f route_to_surface:%.1f\n",sum_Qin,route_to_surface);
      }

    } /* end if redistribution flag */

    return;

} /*end update_drainage_land.c*/

