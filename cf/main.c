/*--------------------------------------------------------------*/
/*                                                              */
/*		create_flowpaths										*/
/*                                                              */
/*  NAME                                                        */
/*		 create_flowpaths										*/
/*                                                              */
/*                                                              */
/*  SYNOPSIS                                                    */
/* 		 create_flowpaths( 								        */
/*                                                              */
/*  OPTIONS                                                     */
/*                                                              */
/*		-v 	Verbose Option										*/
/*		-l    roads to lowest flna interval 					*/
/*		-h  	roads to highest flna interval					*/
/*		-s print drainage statistics							*/
/*		-r 	road flag for drainage statistics					*/
/*		-stream	stream connectivity is assumed						*/
/*			1	random slope value								*/
/*			2	internal slope value							*/
/*			0	no connectivity (default)						*/
/*		-scaledem scale dem values by this amount                     */
/*		-scaletrans	scale streamside transmissivity						*/
/*		-bi	basin ID				*/
/*		-sw	use a sewer image to route water from roads     	*/
/*		-pst	print stream table	*/
/*		-pre	input image file name prefix					*/
/*		-w 	road width (default is 5m)							*/
/*		-a arcview ascii data files (default is GRASS ascii)	*/
/*		-slp slope flag - change use of slope		            */
/*			in computation of gamma			                    */
/*			0  standard				                            */
/*			1  internal slpe of patch		                    */
/*			2 max slope of patch		                    	*/
/*		-o output file name (default -pre opt + _flow_table.dat	*/
/*                                                              */
/*  DESCRIPTION                                                 */
/*                                                              */
/*  revision:  6.0  29 April, 2005                              */
/*  PROGRAMMER NOTES                                            */
/*                                                              */
/*--------------------------------------------------------------*/

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <grass/gis.h>

#include "grassio.h"
#include "blender.h"
#include "glb.h"
#include "sub.h"




main(int argc, char *argv[])  
{
/* local variable declarations */
    int		i/*j*/ /*r*/ /*c*/ /*area*/, num_stream, num_patches;   
	FILE	/**in1*/ *out1, *out2;    
	int		basinid, tmp, maxr, maxc;
	double	cell, width;
	int		pst_flag, f_flag, fl_flag, fh_flag, vflag;
	int		s_flag, r_flag, slp_flag, sc_flag, arc_flag;
	int		st_flag, sewer_flag;
	int		prefix_flag;
	int		suffix_flag;
	int		d_flag;
	double	scale_trans, scale_dem;
	char    input_prefix[MAXS];
	char    output_suffix[MAXS];
	float	scale_factor;



/* filenames for each image and file */

char    fndem[MAXS],  fnpartition[MAXS],    fnK[MAXS],     fnmpar[MAXS],
        fntable[MAXS], fnroot[MAXS], fnpatch[MAXS], fnhill[MAXS], fnzone[MAXS], fnslope[MAXS];
char	fnstream[MAXS];
char	fnroads[MAXS], fnsewers[MAXS];
char	fnflna[MAXS];
char	fnehr[MAXS], fnwhr[MAXS];
char	name[MAXS], name2[MAXS];

char*	fntemplate;
char	rnbasin[MAXS];
char	rnhillslope[MAXS];
char	rnzone[MAXS];
char	rnpatch[MAXS];
char*	rndem;
char*	rnK;
char*	rnmpar;
char*	rnroads;
char*	rnstream;

/* set pointers for images */

    double	    *dem;
    int          *patch;	
    float        *slope;	
    int          *hill;	
    int          *zone;	
    int			 *stream;
    int			 *roads;
    int			 *sewers;
    double		 *K;
	double		 *m_par;
	float		 *flna;
	float		 *ehr;
	float		 *whr;

	struct		flow_struct	*flow_table;

	d_flag 	 = 0;		/* debuf flag					 */
    vflag    = 0;		/* verbose flag					 */
    fl_flag  = 0;		/* roads to lowest flna			 */
    fh_flag  = 0;		/* roads to highest flna		 */
	s_flag	 = 0;		/* printing stats flag			 */
	r_flag 	 = 0;		/* road stats flag				 */
	sc_flag  = 1;		/* stream connectivity flag		 */
	slp_flag = 0;		/* slope use flag		         */
	st_flag  = 0;		/* scaling stream side patches 	 */
	sewer_flag = 0;		/* route through a sewer network */
    arc_flag = 0;		/* arcview input data flag		 */
	scale_trans = 1.0;
	scale_dem = 1.0;	/* scaling for dem values        */
	prefix_flag = 0;	/* input prefix flag             */
	suffix_flag = 0;	/* output suffix flag            */
	pst_flag = 0;		/* print stream table flag            */
	cell = 30.0;		/* default resolution            */
	width = 5;			/* default road width            */
	basinid = 1;

	// GRASS init
	G_gisinit(argv[0]);

	// GRASS module header
	struct GModule* module;
	module = G_define_module();
	module->keywords = "RHESSys";
	module->description = "Creates a flowpaths file for input into RHESSys";

	// GRASS arguments
	struct Flag* debug_flag = G_define_flag();
	debug_flag->key = 'g';
	debug_flag->description = "Enable printouts during compuation of flowpaths";

	struct Flag* lowest_flna_flag = G_define_flag();
	lowest_flna_flag->key = 'l';
	lowest_flna_flag->description = "Roads to lowest flna interval";

	struct Flag* highest_flna_flag = G_define_flag();
	highest_flna_flag->key = 'h';
	highest_flna_flag->description = "Roads to highest flna interval";

	struct Flag* drainage_stats_flag = G_define_flag();
	drainage_stats_flag->key = 'd';
	drainage_stats_flag->description = "Print drainage statistics";

	struct Flag* road_drainage_stats_flag = G_define_flag();
	road_drainage_stats_flag->key = 'r';
	road_drainage_stats_flag->description = "Road flag for drainage statistics";

	struct Option* stream_connectivity_opt = G_define_option();
	stream_connectivity_opt->key = "streamcon";
	stream_connectivity_opt->type = TYPE_STRING;
	stream_connectivity_opt->required = NO;
	stream_connectivity_opt->description = "Stream connectivity type: [random(default), internal, none]";
	stream_connectivity_opt->multiple=NO;

	struct Option* scale_dem_opt = G_define_option();
	scale_dem_opt->key = "scaledem";
	scale_dem_opt->type = TYPE_DOUBLE;
	scale_dem_opt->required = NO;
	scale_dem_opt->description = "DEM scaling factor";

	struct Option* scale_stream_trans = G_define_option();
	scale_stream_trans->key = "scaletrans";
	scale_stream_trans->type = TYPE_DOUBLE;
	scale_stream_trans->required = NO;
	scale_stream_trans->description = "Scaleing factor for streamside transmissivity";

	struct Flag* use_sewer_flag = G_define_flag();
	use_sewer_flag->key = 's';
	use_sewer_flag->description = "Use a sewer image to route water from roads";

	struct Flag* print_stream_table_flag = G_define_flag();
	print_stream_table_flag->key = 'p';
	print_stream_table_flag->description = "Print stream table";

	struct Option* prefix_opt = G_define_option();
	prefix_opt->key = "prefix";
	prefix_opt->type = TYPE_STRING;
	prefix_opt->required = YES;
	prefix_opt->description = "Prefix for reading in ASCII files, this will go away.";

	struct Option* road_width_opt = G_define_option();
	road_width_opt->key = "roadwidth";
	road_width_opt->type = TYPE_DOUBLE;
	road_width_opt->required = NO;
	road_width_opt->description = "Road width(m). [Default 5m]";
	
	struct Option* slope_use_opt = G_define_option();
	slope_use_opt->key = "slopeuse";
	slope_use_opt->type = TYPE_STRING;
	slope_use_opt->required = NO;
	slope_use_opt->description = "Change the use of slope in the compuation of gamma [standard(default), internal, max]";

	struct Option* basin_id_opt = G_define_option();
	basin_id_opt->key = "basinid";
	basin_id_opt->type = TYPE_INTEGER;
	basin_id_opt->required = NO;
	basin_id_opt->description = "Basin ID";

	struct Option* output_suffix_opt = G_define_option();
	output_suffix_opt->key = "output";
	output_suffix_opt->type = TYPE_STRING;
	output_suffix_opt->required = NO;
	output_suffix_opt->description = "Output suffix";

	// Arguments that specify the names of required raster maps
	struct Option* m_raster_opt = G_define_option();
	m_raster_opt->key = "m";
	m_raster_opt->type = TYPE_STRING;
	m_raster_opt->required = YES;
	m_raster_opt->description = "m";

	struct Option* K_raster_opt = G_define_option();
	K_raster_opt->key = "K";
	K_raster_opt->type = TYPE_STRING;
	K_raster_opt->required = YES;
	K_raster_opt->description = "K";

	struct Option* stream_raster_opt = G_define_option();
	stream_raster_opt->key = "stream";
	stream_raster_opt->type = TYPE_STRING;
	stream_raster_opt->required = YES;
	stream_raster_opt->description = "stream";

	struct Option* road_raster_opt = G_define_option();
	road_raster_opt->key = "road";
	road_raster_opt->type = TYPE_STRING;
	road_raster_opt->required = YES;
	road_raster_opt->description = "road";

	struct Option* dem_raster_opt = G_define_option();
	dem_raster_opt->key = "dem";
	dem_raster_opt->type = TYPE_STRING;
	dem_raster_opt->required = YES;
	dem_raster_opt->description = "dem";

	struct Option* slope_raster_opt = G_define_option();
	slope_raster_opt->key = "slope";
	slope_raster_opt->type = TYPE_STRING;
	slope_raster_opt->required = YES;
	slope_raster_opt->description = "slope";

	struct Option* template_opt = G_define_option();
	template_opt->key = "template";
	template_opt->type = TYPE_STRING;
	template_opt->required = YES;
	template_opt->description = "REHSSys template file from which to extract the basin, hill, zone, and patch GRASS rasters";
	
	// Parse GRASS arguments
	if (G_parser(argc, argv))
		exit(1);

	// Get values from GRASS arguments
	fl_flag = lowest_flna_flag->answer;
	fh_flag = highest_flna_flag->answer;
	if (fl_flag || fh_flag) {
		f_flag = 1;
	} else {
		f_flag = 0;
	}
	s_flag = drainage_stats_flag->answer;
	r_flag = road_drainage_stats_flag->answer;

	if (stream_connectivity_opt->answer != NULL) {  
		// Default is 1, random connectivity. See sc_flag declaration for setting default.
		if (strcmp("random",stream_connectivity_opt->answer)) {
			sc_flag = 1;
		} else if (strcmp("internal", stream_connectivity_opt->answer)) {
			sc_flag = 2;
		} else if (strcmp("none", stream_connectivity_opt->answer)) {
			sc_flag = 0;
		} else {
			G_fatal_error("\"%s\" is not a valid argument to stream", stream_connectivity_opt->answer);
		}
	}

	if (scale_dem_opt->answer != NULL) {	
		// Default is set at declaration, only modify if set
		if (sscanf(scale_dem_opt->answer, "%lf", &scale_dem) != 1) {
			G_fatal_error("Error setting the scale dem value");
		}
	}
	
	if (scale_stream_trans->answer != NULL) {
		// Default is set at declaration, only modify if set
		if (sscanf(scale_stream_trans->answer, "%lf", &scale_trans) != 1) {
			G_fatal_error("Error setting the scale trans value");
		}
	}

	sewer_flag = use_sewer_flag->answer;
	pst_flag = print_stream_table_flag->answer;
	strcpy(input_prefix, prefix_opt->answer);

	if (road_width_opt->answer != NULL) {
		// Default is set at declaration
		if (sscanf(road_width_opt->answer, "%lf", &width) != 1) {
			G_fatal_error("Error setting the road width value");
		}
	}

	if (slope_use_opt->answer != NULL) {
		if (strcmp("standard", slope_use_opt->answer)) {
			slp_flag = 0;
		} else if (strcmp("internal", slope_use_opt->answer)) {
			slp_flag = 1;
		} else if (strcmp("max", slope_use_opt->answer)) {
			slp_flag = 2;
		} else {
			G_fatal_error("\"%s\" is not a valid argument to slopeuse", slope_use_opt->answer);
		}
	}

	if (basin_id_opt->answer != NULL) {
		// Default set at declaration
		if (sscanf(basin_id_opt->answer, "%d", &basinid) != 1) {
			G_fatal_error("Error setting the basin ID value");
		}
	}

	if (output_suffix_opt->answer != NULL) {
		strcpy(output_suffix, output_suffix_opt->answer);	
		suffix_flag = 1;	
	} else {
		strcpy(output_suffix, "_flow_table.dat");
	}
	
		// Need to implement verbose	
	rndem = dem_raster_opt->answer;
	fntemplate = template_opt->answer;
	rnK = K_raster_opt->answer;
	rnmpar = m_raster_opt->answer;
	rnroads = road_raster_opt->answer;
	rnstream = stream_raster_opt->answer;



    printf("Create_flowpaths.C\n\n");

	
    input_prompt(&maxr, &maxc, input_prefix, fndem,fnslope,fnK,fnflna,
					fnpatch, fnzone,fnhill, fnstream, fnroads, fnsewers, fnmpar,fnpartition,
		 fntable,fnroot, fnehr, fnwhr, f_flag, sewer_flag, arc_flag);


	// Read in the names of the basin, hill, zone, and patch maps from the
	// template file.
	FILE* template_fp = fopen(fntemplate, "r");
	if (template_fp == NULL) {
		G_fatal_error("Can not open template file <%s>", fntemplate);
	}

	char	template_buffer[MAXS];
	char	first[MAXS];
	char	second[MAXS];
	while (fgets(template_buffer, sizeof(template_buffer), template_fp) != NULL) {
		sscanf(template_buffer, "%s %s", first, second);

		// Check if the token is anything we are looking for
		if (strcmp("_basin", first) == 0 ) {
			strcpy(rnbasin, second);
		} else if (strcmp("_hillslope", first) == 0 ) {
			strcpy(rnhillslope, second);
		} else if (strcmp("_zone", first) == 0 ) {
			strcpy(rnzone, second);
		} else if (strcmp("_patch", first) == 0 ) {
			strcpy(rnpatch, second);
		}
	}
	fclose (template_fp);
	printf("Basin: %s\n", rnbasin);
	printf("Hillslope: %s\n", rnhillslope);
	printf("Zone: %s\n", rnzone);
	printf("Patch: %s\n", rnpatch);



	/* open some diagnostic output files */

  strcpy(name, input_prefix);
  strcat(name, ".build"); 
  if ( (out1 = fopen(name, "w")) == NULL)
        {
        printf("cannot open file \n");
        exit(1);
       	} 

  strcpy(name2, input_prefix);
  strcat(name2, ".pit"); 
  if ( (out2 = fopen(name2, "w")) == NULL)
        {
        printf("cannot open file \n");
        exit(1);
       	} 


	/* allocate and input map images */
	// figure out what's happening with maxr, maxc, possible raster
	// size mismatch
	struct Cell_head dem_header;
	dem = raster2array(rndem, &dem_header, &maxr, &maxc);

	struct Cell_head patch_header;
	patch   = (int *) calloc(maxr*maxc, sizeof(int));      
    input_ascii_int(patch, fnpatch, maxr, maxc, arc_flag);	

	struct Cell_head zone_header;
	zone   = (int *) malloc(maxr*maxc*sizeof(int));
    input_ascii_int(zone, fnzone, maxr, maxc, arc_flag);	

	struct Cell_head hill_header;
    hill   = (int *) malloc(maxr*maxc*sizeof(int));
    input_ascii_int(hill, fnhill, maxr, maxc, arc_flag);	

	struct Cell_head stream_header;
    stream   = (int *) malloc(maxr*maxc*sizeof(int));
    input_ascii_int(stream, fnstream, maxr, maxc, arc_flag);	

   if ((sc_flag == 1) || (slp_flag > 0) ){ 

    	slope   = (float *) malloc(maxr*maxc*sizeof(int));
		scale_factor = 1;
    	input_ascii_float(slope, fnslope, maxr, maxc, arc_flag, scale_factor);	
	}

	struct Cell_head roads_header;
    roads = (int *) malloc(maxr*maxc*sizeof(int));
    input_ascii_int(roads, fnroads, maxr, maxc, arc_flag);	

    if (sewer_flag == 1) {
    	sewers = (int *) malloc(maxr*maxc*sizeof(int));
    	input_ascii_int(sewers, fnsewers, maxr, maxc, arc_flag);	
	}

	struct Cell_head K_header;
	K = raster2array(rnK, &K_header, NULL, NULL);

	struct Cell_head mpar_header;
	m_par = raster2array(rnmpar, &mpar_header, NULL, NULL);
	
	if (f_flag) {
		scale_factor = 1.0;
		flna = (float *)  malloc(maxr*maxc*sizeof(float));
		input_ascii_float(flna, fnflna, maxr, maxc, arc_flag, scale_factor);
		}
	else flna = NULL;


	/* allocate flow table */
	flow_table = (struct flow_struct *)calloc((maxr*maxc),sizeof(struct flow_struct));


	printf("\n Building flow table");
	num_patches = build_flow_table(flow_table, dem, slope, hill, zone, patch, 
					stream, roads, sewers, K, m_par, flna, out1, maxr, 
					maxc,f_flag, sc_flag, sewer_flag, slp_flag, cell, 
                    scale_dem);

	fclose(out1);

	/* processes patches - computing means and neighbour slopes and gammas */
	printf("\n Computing gamma");
	num_stream = compute_gamma(flow_table, num_patches, out2,scale_trans,cell,sc_flag,
                        slp_flag, d_flag);


	
	 /* remove pits and re-order patches appropriately */
	 
	printf("\n Removing pits");
	   remove_pits(flow_table, num_patches, sc_flag, slp_flag, cell, out2); 
	 
	/* add roads */
	printf("\n Adding roads");
	 add_roads(flow_table, num_patches, out2, cell);

	/* find_receiving patch for flna options */
	if (f_flag)
		route_roads_to_patches(flow_table, num_patches, fl_flag);
	

	printf("\n Computing upslope area");
 		tmp = compute_upslope_area(flow_table, num_patches, out2, r_flag,cell);

	if (s_flag == 1) {
		printf("\n Printing drainage stats");
		print_drain_stats(num_patches, flow_table);
		tmp = compute_dist_from_road(flow_table, num_patches, out2, cell);	
 		tmp = compute_drainage_density(flow_table, num_patches, cell);
		}
	

	printf("\n Printing flowtable");
	print_flow_table(num_patches,flow_table,sc_flag,slp_flag,cell,scale_trans,
                     input_prefix, output_suffix, width);

	if (pst_flag == 1) {
		printf("\n Printing  stream table");
		print_stream_table(num_patches,num_stream, flow_table,sc_flag,slp_flag,cell,scale_trans,
                     input_prefix, output_suffix, width, basinid);
		}

	
	fclose(out2);
	printf("\n Finished Createflowpaths \n\n");
    exit(0);


} /* end create_flowpaths.c */


