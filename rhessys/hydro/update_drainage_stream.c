/*--------------------------------------------------------------*/
/* 											*/
/*					update_drainage_stream			*/
/*											*/
/*	update_drainage_stream.c - creates a patch object				*/
/*											*/
/*	NAME										*/
/*	update_drainage_stream.c - creates a patch object				*/
/*											*/
/*	SYNOPSIS									*/
/*	void update_drainage_stream( 							*/
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


void  update_drainage_stream(
								 struct patch_object *patch,
								 struct command_line_object *command_line,
								 double time_int,
								 int verbose_flag)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	double  compute_return_flow(
		int,
		double  ,
		double  );
	
	double  compute_delta_water(
		int,
		double,
		double,
		double,
		double,
		double);
	
//	double compute_N_leached(
//		int,
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
//		double, double *);
	
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
		struct patch_object *patch);

	double recompute_gamma(	
		struct patch_object *,
		double);
	/*--------------------------------------------------------------*/ 
	/*	Local variable definition.				*/ 
	/*--------------------------------------------------------------*/ 
	int i, j,k, d; 
	double m, Ksat; 
	double return_flow;  /* m */ 
	double NO3_leached_total, NO3_leached_to_stream; /* kg/m2 */ 
	double NH4_leached_total, NH4_leached_to_stream; /* kg/m2 */ 
	double DON_leached_total, DON_leached_to_stream; /* kg/m2 */ 
	double DOC_leached_total, DOC_leached_to_stream; /* kg/m2 */ 
	double patch_int_depth;  /* m of H2O */
	double  route_to_stream; /* m3 */
	double route_to_surface;
	double  Qin, Qout,Qstr_total;  /* m */
	double gamma, total_gamma, percent_tobe_routed;
	double Nin, Nout;  /* kg/m2 */
	double t1,t2,t3;
	
	d=0;
	route_to_stream = 0.0;
	return_flow=0.0;
	NO3_leached_to_stream = 0.0;
	NH4_leached_to_stream = 0.0;
	DON_leached_to_stream = 0.0;
	DOC_leached_to_stream = 0.0;
	/*--------------------------------------------------------------*/
	/*	m and K are multiplied by sensitivity analysis variables */
	/*--------------------------------------------------------------*/
	m = patch[0].m ;
	Ksat = patch[0].soil_defaults[0][0].Ksat_0 ;

	/*--------------------------------------------------------------*/
	/*	for now there should be no recomputing of gamma for 	*/
	/*	streams because they do not route water to downslope	*/
	/*	neighbours						*/
	/*--------------------------------------------------------------*/
    total_gamma = patch[0].innundation_list[d].gamma;

	/*------------------------------------------------------------*/
	/*	calculate amuount of water output to stream as baseflow */
	/*-----------------------------------------------------------*/
	if (total_gamma < ZERO ) {
		gamma = Ksat * m * 2.0 * sqrt(patch[0].area)
			* time_int;
	}
	else {
		gamma = total_gamma * time_int;

	}

	route_to_stream = compute_varbased_flow(
		patch[0].num_soil_intervals,
		patch[0].std * command_line[0].std_scale,
		patch[0].sat_deficit,
		gamma,
		patch[0].soil_defaults[0][0].interval_size,
		patch[0].transmissivity_profile,
		patch);

	if (route_to_stream < 0.0) route_to_stream = 0.0;

		
	/*--------------------------------------------------------------*/
	/* compute Nitrogen leaching amount with baseflow		*/
	/*--------------------------------------------------------------*/
    double lpools[] = {patch[0].soil_ns.nitrate,
                      patch[0].soil_ns.sminn,
                      patch[0].soil_ns.DON,
                      patch[0].soil_cs.DOC
                      };
    double leached[LEACH_ELEMENT_counts];
	if (command_line[0].grow_flag > 0) {

            double t = compute_N_leached(
			verbose_flag,
            lpools,
			route_to_stream / patch[0].area,
			patch[0].sat_deficit,
			patch[0].soil_defaults[0][0].soil_water_cap,
            //m,
            //gamma / patch[0].area,
			patch[0].soil_defaults[0][0].porosity_0,
			patch[0].soil_defaults[0][0].porosity_decay,
            patch[0].soil_defaults[0][0].decay_rate,
			patch[0].soil_defaults[0][0].active_zone_z,
			patch[0].soil_defaults[0][0].soil_depth,
            patch[0].soil_defaults[0][0].adsorption_rate,
            leached
            //patch[0].transmissivity_profile
                );
        NO3_leached_to_stream = leached[LNO3];
        NH4_leached_to_stream = leached[LNH4];
        DON_leached_to_stream = leached[LDON];
        DOC_leached_to_stream = leached[LDOC];
        //free(leached_to_stream);

        patch[0].soil_ns.NO3_Qout += NO3_leached_to_stream;
        patch[0].soil_ns.NH4_Qout += NH4_leached_to_stream;
        patch[0].soil_ns.DON_Qout += DON_leached_to_stream;
        patch[0].soil_cs.DOC_Qout += DOC_leached_to_stream;
        patch[0].streamflow_NO3 += NO3_leached_to_stream;
        patch[0].streamNO3_from_sub += NO3_leached_to_stream;
        patch[0].hourly[0].streamflow_NO3 += NO3_leached_to_stream;
        patch[0].hourly[0].streamflow_NO3_from_sub += NO3_leached_to_stream;

        patch[0].streamflow_NH4 += NH4_leached_to_stream;
        patch[0].streamflow_DON += DON_leached_to_stream;
        patch[0].streamflow_DOC += DOC_leached_to_stream;
	}

	patch[0].Qout += (route_to_stream / patch[0].area);
	patch[0].base_flow += (route_to_stream / patch[0].area);
	patch[0].hourly_subsur2stream_flow += route_to_stream / patch[0].area;



	/*--------------------------------------------------------------*/
	/*	calculate any return flow to the stream in this patch   */
	/*	and route any infiltration excess			*/
	/*--------------------------------------------------------------*/
	if ((patch[0].sat_deficit-patch[0].rz_storage-patch[0].unsat_storage) < -1.0*ZERO) {
		return_flow = compute_varbased_returnflow(patch[0].std * command_line[0].std_scale, 
			patch[0].rz_storage+patch[0].unsat_storage,
			patch[0].sat_deficit, &(patch[0].litter));
		patch[0].detention_store += return_flow;  
		patch[0].sat_deficit += (return_flow - (patch[0].unsat_storage+patch[0].rz_storage));;
		patch[0].unsat_storage = 0.0;
		patch[0].rz_storage = 0.0;
	}

	/*--------------------------------------------------------------*/
	/*	calculated any N-transport associated with return flow  */
	/*	- note available N reduced by what has already been lost  */
	/*	due to subsurface routing above				*/
	/* 	note only nitrate is assumed to follow return flow		*/
	/*--------------------------------------------------------------*/
	if (return_flow > ZERO) {
        for (int i = 0; i < LEACH_ELEMENT_counts; i++) {
            switch(i) {
            case LNO3:
                lpools[i] = patch[0].soil_ns.nitrate - NO3_leached_to_stream;
                break;
            case LNH4:
                lpools[i] = patch[0].soil_ns.sminn - NH4_leached_to_stream;
                break;
            case LDON:
                lpools[i] = patch[0].soil_ns.DON - DON_leached_to_stream;
                break;
            case LDOC:
                lpools[i] = patch[0].soil_cs.DOC - DOC_leached_to_stream;
            }
        }
        //double *pNout =
        double t = compute_N_leached(
			verbose_flag,
            lpools,
			return_flow,
			0.0,
			0.0,
            //m,
            //gamma / patch[0].area,
			patch[0].soil_defaults[0][0].porosity_0,
			patch[0].soil_defaults[0][0].porosity_decay,
            patch[0].soil_defaults[0][0].decay_rate,
			patch[0].soil_defaults[0][0].active_zone_z,
			patch[0].soil_defaults[0][0].soil_depth,
            patch[0].soil_defaults[0][0].adsorption_rate,
            leached
            //patch[0].transmissivity_profile
                );
        patch[0].surface_NO3 += leached[LNO3];
        patch[0].soil_ns.NO3_Qout += leached[LNO3];

        patch[0].surface_NH4 += leached[LNH4];
        patch[0].soil_ns.NH4_Qout += leached[LNH4];

        patch[0].surface_DON += leached[LDON];
        patch[0].soil_ns.DON_Qout += leached[LDON];

        patch[0].surface_DOC += leached[LDON];
        patch[0].soil_cs.DOC_Qout += leached[LDON];
        //free(pNout);

	}

	/*--------------------------------------------------------------*/
	/*	route water and nitrogen lossed due to infiltration excess */
	/*	note we assume that this happens before return_flow losses */
	/*--------------------------------------------------------------*/

	if ((patch[0].detention_store > patch[0].soil_defaults[0][0].detention_store_size) &&
		(patch[0].detention_store > ZERO)) {
		Qout = (patch[0].detention_store - patch[0].soil_defaults[0][0].detention_store_size);
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_NO3;
		patch[0].surface_NO3  -= Nout;
		patch[0].streamflow_NO3 += Nout;
		patch[0].hourly[0].streamflow_NO3 += Nout;
		patch[0].streamNO3_from_surface +=Nout;
		patch[0].hourly[0].streamflow_NO3_from_surface +=Nout;

		patch[0].surface_ns_leach += Nout;
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_DOC;
		patch[0].surface_DOC  -= Nout;
		patch[0].streamflow_DOC += Nout;
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_DON;
		patch[0].surface_DON  -= Nout;
		patch[0].streamflow_DON += Nout;
		Nout = (min(1.0, Qout / patch[0].detention_store)) * patch[0].surface_NH4;
		patch[0].surface_NH4  -= Nout;
		patch[0].streamflow_NH4 += Nout;
		patch[0].detention_store -= Qout;
		patch[0].return_flow += Qout; 
		patch[0].hourly_sur2stream_flow += Qout;
		}

} /*end update_drainage_stream.c*/

