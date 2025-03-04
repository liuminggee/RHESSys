/*--------------------------------------------------------------*/
/* 																*/
/*					construct_routing_topology					*/
/*																*/
/*	construct_routing_topology.c - creates a patch object		*/
/*																*/
/*	NAME														*/
/*	construct_routing_topology.c - creates a patch object		*/
/*																*/
/*	SYNOPSIS													*/
/*	struct routing_list_object construct_routing_topology( 		*/
/*							struct basin_object *basin)			*/
/*																*/
/* 																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*  reads routing topology from input file						*/
/*	creates neighbourhood structure for each patch in the basin */
/*	returns a list giving order for patch-level routing			*/
/*																*/
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rhessys.h"

struct routing_list_object *construct_routing_topology(
      FILE *routing_file,
      struct hillslope_object *hillslope,
		  struct command_line_object *command_line,
		  bool surface)

{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	//struct patch_object *find_patch(int, int, int, struct basin_object *);
	struct patch_object *find_patch_in_hillslope(int, int, struct hillslope_object *);
    struct hillslope_object *find_hillslope_in_basin(int, struct basin_object *);

	int assign_neighbours_in_hillslope (struct neighbour_object *,
		int,
    struct hillslope_object *,
		FILE *);

	void *alloc(size_t, char *, char *);

	double * compute_transmissivity_curve( double, struct patch_object *, struct command_line_object *);

	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int		i,d,j;
	int		num_patches, num_neighbours;
	int		patch_ID, zone_ID, hill_ID;
	int		drainage_type;
	double	x,y,z, area, gamma, width;
	struct routing_list_object	*rlist;
	struct	patch_object	*patch;
	struct	patch_object	*stream;
	struct	innundation_object *innundation_list;


	rlist = (struct routing_list_object	*)alloc( sizeof(struct routing_list_object), "rlist", "construct_routing_topology");

	/*--------------------------------------------------------------*/
	/*  Try to open the routing file in read mode.                    */
	/*--------------------------------------------------------------*/
	//if ( (routing_file = fopen(routing_filename,"r")) == NULL ){
	//	fprintf(stderr,"FATAL ERROR:  Cannot open routing file %s\n",
	//		routing_filename);
	//	exit(EXIT_FAILURE);
	//} /*end if*/

	fscanf(routing_file,"%d",&num_patches);
	rlist->num_patches = num_patches;
	rlist->list = (struct patch_object **)alloc(
		num_patches * sizeof(struct patch_object *), "patch list",
		"construct_routing_topography");

	/*--------------------------------------------------------------*/
	/*	Read in  each patch record and find it		.				*/
	/*	if it is a stream add it to the basin level routing list	*/
	/*	otherwise add it to the hillslope level routing list		*/
	/*--------------------------------------------------------------*/
	for (i=0; i< num_patches; ++i) {
        patch = 0;
		fscanf(routing_file,"%d %d %d %lf %lf %lf %lf %lf %d %lf %d",
			&patch_ID,
			&zone_ID,
			&hill_ID,
			&x,&y,&z,
			&area,
			&area,
			&drainage_type,
            &gamma,                                                             //09092022LML: (m2) patch area X mean_slope (in terms of del_z/dist; or gamma_tot/perim_sum; gamma_tot = sum(tp_perimeter X slope_ij)))
                                                                                //09092022LML: SHOULD BE DOWNSTREAM BOUNDARY LENGTH x AVG_SLOPE (m)
			&num_neighbours);

        if  ( (hillslope != 0) && (patch_ID != 0) && (zone_ID != 0) && (hill_ID != 0) ) {
			patch = find_patch_in_hillslope(patch_ID, zone_ID, hillslope );
        } else {
            fprintf(stderr,"\nFATAL ERROR: in construct_routing_topology, could not findpatch with ID:%d in hillslope ID:%d.\n", patch_ID, hill_ID );
            exit(EXIT_FAILURE);
        }

		rlist->list[i] = patch;
        if (patch != 0) {
          if ((patch[0].soil_defaults[0][0].Ksat_0 < ZERO))
			printf("\n WARNING lateral Ksat (%lf) are close to zero for patch %d",
				patch[0].soil_defaults[0][0].Ksat_0, patch[0].ID);

          if (patch[0].soil_defaults[0][0].m < ZERO) {
		 	gamma = gamma * patch[0].soil_defaults[0][0].Ksat_0;
          } else {
#ifndef LIU_GAMMA_TRANSMISSIVITY_NEW
            gamma = gamma * patch[0].soil_defaults[0][0].m * patch[0].soil_defaults[0][0].Ksat_0;   //08312022LML: seems wrong if suppose to use average Ksat
#else

            //09012022LML
            //double avg_ksat_coefz = patch[0].soil_defaults[0][0].m_z
            //        * (1.0 - exp(-patch[0].soil_defaults[0][0].soil_depth)/patch[0].soil_defaults[0][0].m_z)
            //        / (patch[0].soil_defaults[0][0].soil_depth);
            double avg_ksat_coef = patch[0].soil_defaults[0][0].m
                    * (1.0 - exp(-patch[0].soil_defaults[0][0].soil_water_cap/patch[0].soil_defaults[0][0].m))
                    / (patch[0].soil_defaults[0][0].soil_water_cap);

            gamma = gamma * patch[0].soil_defaults[0][0].Ksat_0 * avg_ksat_coef; //(m3/day)
                                                                                 //SHOULD BE (m2/day) if original gamma input is in length X slope
#endif
          }

		/*--------------------------------------------------------------*/
		/*  Allocate innundation list array				*/
		/*	note for this routing there is only one innundation depth 	*/
		/*	however it is need to be compatablability 		*/
		/*--------------------------------------------------------------*/
          d=0;
          if ( surface ) {
			patch->surface_innundation_list = (struct innundation_object *)alloc( 1 *
								sizeof(struct innundation_object), "surface_innundation_list", "construct_routing_topology");
			innundation_list = patch->surface_innundation_list;
          } else {
			patch->innundation_list = (struct innundation_object *)alloc( 1 *
					sizeof(struct innundation_object), "innundation_list", "construct_routing_topology");
			innundation_list = patch->innundation_list;
          }

          if ( surface ) {
			patch[0].num_innundation_depths = 1;
          }

          innundation_list->num_neighbours = num_neighbours;

          innundation_list->gamma = gamma; //09092022LML (m)

		// TODO: what should critical depth be for a surface flow table?
          innundation_list->critical_depth = NULLVAL;

          if ( !surface ) {
			patch[0].stream_gamma = 0.0;
			patch[0].drainage_type = drainage_type;
			if ((patch[0].drainage_type != STREAM)
					&& (patch[0].innundation_list[d].gamma < ZERO)) {
				/*printf(
						"\n non-stream patches with zero gamma %d switched to stream for now",
						patch[0].ID);*/
				patch[0].drainage_type = STREAM;
			}
          }

		/*--------------------------------------------------------------*/
		/*  Allocate neighbour array									*/
		/*--------------------------------------------------------------*/
          innundation_list->neighbours = (struct neighbour_object *)alloc(num_neighbours *
				sizeof(struct neighbour_object), "neighbours", "construct_routing_topology");
          num_neighbours = assign_neighbours_in_hillslope(innundation_list->neighbours, num_neighbours, hillslope, routing_file);
          if ((num_neighbours == -9999) && (patch[0].drainage_type != STREAM)) {
			printf("\n WARNING sum of patch %d neigh gamma is not equal to 1.0", patch[0].ID);
          } else {
			innundation_list->num_neighbours = num_neighbours;
          }

          if ( drainage_type == ROAD ) {
			fscanf(routing_file,"%d %d %d %lf",
				&patch_ID,
				&zone_ID,
				&hill_ID,
				&width);
			// TODO: Decide if we need separate stream_gamma, road_cut_depth, and next_stream values for surface flow table
			if ( !surface ) {
				patch[0].stream_gamma = gamma;
				patch[0].road_cut_depth = width * tan(patch[0].slope);
                //09222023LML TODO: if the downstream of the road is another hillslope
                //, need find the patch from that hillslope.
                if (hill_ID != hillslope->ID) {
                    struct hillslope_object *hill = find_hillslope_in_basin(hill_ID,hillslope->basin);
                    stream = find_patch_in_hillslope(patch_ID, zone_ID, hill);
                } else {
                  stream = find_patch_in_hillslope(patch_ID, zone_ID, hillslope);
                }
				patch[0].next_stream = stream;
			}
          }

          if ( !surface ) {
			/*--------------------------------------------------------------*/
			/*	create a vector of transmssivities 			*/
			/*--------------------------------------------------------------*/
			patch[0].num_soil_intervals = (int) lround(patch[0].soil_defaults[0][0].soil_water_cap / patch[0].soil_defaults[0][0].interval_size);
			if (patch[0].num_soil_intervals > MAX_NUM_INTERVAL) {
				patch[0].num_soil_intervals = MAX_NUM_INTERVAL;
				patch[0].soil_defaults[0][0].interval_size = patch[0].soil_defaults[0][0].soil_water_cap / MAX_NUM_INTERVAL;
				}
			patch[0].transmissivity_profile = compute_transmissivity_curve(gamma, patch, command_line);
          }
        } else {//if patch!=0

        }
	}

	//fclose(routing_file);

	return(rlist);
} /*end construct_routing_topology.c*/
