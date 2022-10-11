/*--------------------------------------------------------------*/
/* 								*/
/*	construct_soil_defaults					*/
/*								*/
/*	construct_soil_defaults.c - makes patch default		*/
/*			objects.				*/
/*								*/
/*	NAME							*/
/*	construct_soil_defaults.c - makes patch default		*/
/*			objects.				*/
/*								*/
/*	SYNOPSIS						*/
/*	struct soil_default *construct_soil_defaults(           */
/*			num_default_files,			*/
/*			default_files,				*/
/*			grow_flag,				*/
/*			default_object_list )			*/
/*								*/
/*	OPTIONS							*/
/*								*/
/*	DESCRIPTION						*/
/*								*/
/*	PROGRAMMER NOTES					*/
/*								*/
/*	Original code, January 15, 1996.			*/
/*	July 28, 1997	C.Tague					*/
/*	removed capillary rise soil variables i.e rooting depth */
/*	pore size index and suction				*/
/*								*/
/*	Sep 15 97 RAF						*/
/*	added cap rise variables back as we now do a cap rise	*/
/*	routine but removed wilting point which was never used.	*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "params.h"

struct soil_default *construct_soil_defaults(
			int	num_default_files,
			char	**default_files,
			struct command_line_object *command_line)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.				*/
	/*--------------------------------------------------------------*/
	void	*alloc(	size_t,
		char	*,
		char	*);

	double compute_delta_water(int, double, double,	double, double, double);
	int	parse_albedo_flag( char *);

	/*--------------------------------------------------------------*/
	/*	Local variable definition.				*/
	/*--------------------------------------------------------------*/
        int strbufLen = 256;
        int filenameLen = 1024;
	int	i;
	double 	ftmp,soil;
	FILE	*default_file;
        char	strbuf[strbufLen];
        char	outFilename[filenameLen];
	char	*newrecord;
	char	record[MAXSTR];
	struct 	soil_default *default_object_list;
	void	*alloc(	size_t, char *, char *);
        param *paramPtr = NULL;
        int paramCnt = 0;

	/*--------------------------------------------------------------*/
	/*	Allocate an array of default objects.			*/
	/*--------------------------------------------------------------*/
	default_object_list = (struct soil_default *)
		alloc(num_default_files * sizeof(struct soil_default),"default_object_list",
		"construct_soil_defaults");

	/*--------------------------------------------------------------*/
	/*	Loop through the default files list.			*/
	/*--------------------------------------------------------------*/
	for (i=0 ; i<num_default_files; i++){
		/*--------------------------------------------------------------*/
		/*		Try to open the ith default file.		*/
		/*--------------------------------------------------------------*/
		//if ( (default_file = fopen( default_files[i], "r")) == NULL ){
		//	fprintf(stderr,"FATAL ERROR:in construct_soil_defaults",
		//		"unable to open defaults file %d.\n",i);
		//	exit(EXIT_FAILURE);
		//} /*end if*/

		printf("Reading %s\n", default_files[i]);
                paramCnt = 0;
                if (paramPtr != NULL)
                    free(paramPtr);

                paramPtr = readParamFile(&paramCnt, default_files[i]);

		/*--------------------------------------------------------------*/
		/*		read the ith default file into the ith object.	*/
		/*--------------------------------------------------------------*/
        struct 	soil_default *defobj = &default_object_list[i];
        defobj[0].ID = 			getIntParam(&paramCnt, &paramPtr, "patch_default_ID", "%d", 3, 1);
        defobj[0].theta_psi_curve = 	getIntParam(&paramCnt, &paramPtr, "theta_psi_curve", "%d", 1.0, 1);
        defobj[0].Ksat_0 = 		getDoubleParam(&paramCnt, &paramPtr, "Ksat_0", "%lf", 3.0, 1);
        defobj[0].m = 			getDoubleParam(&paramCnt, &paramPtr, "m", "%lf", 0.12, 1);
        defobj[0].porosity_0 = 		getDoubleParam(&paramCnt, &paramPtr, "porosity_0", "%lf", 0.435, 1);
        defobj[0].porosity_decay = 	getDoubleParam(&paramCnt, &paramPtr, "porosity_decay", "%lf", 4000.0, 1);
        defobj[0].p3 = 			getDoubleParam(&paramCnt, &paramPtr, "P3", "%lf", 0.0, 1); // param name upper case in param file
        defobj[0].pore_size_index = 	getDoubleParam(&paramCnt, &paramPtr, "pore_size_index", "%lf", 0.204, 1);
        defobj[0].psi_air_entry = 		getDoubleParam(&paramCnt, &paramPtr, "psi_air_entry", "%lf", 0.218, 1);
        defobj[0].psi_max = 		getDoubleParam(&paramCnt, &paramPtr, "psi_max", "%lf", 0.01, 1);
        defobj[0].soil_depth = 		getDoubleParam(&paramCnt, &paramPtr, "soil_depth", "%lf", 200.0, 1);
         printf("soil_depth defs is: %lf\n",defobj[0].soil_depth);
        defobj[0].m_z = 			getDoubleParam(&paramCnt, &paramPtr, "m_z", "%lf", 0.4, 1);
        defobj[0].detention_store_size = 	getDoubleParam(&paramCnt, &paramPtr, "detention_store_size", "%lf", 0.0, 1);
        defobj[0].deltaz = 		getDoubleParam(&paramCnt, &paramPtr, "deltaZ", "%lf", 1.0, 1); // param name contains uppercase "Z" in param file
        defobj[0].active_zone_z = 		getDoubleParam(&paramCnt, &paramPtr, "active_zone_z", "%lf", 5.0, 1);
        if(defobj[0].active_zone_z > defobj[0].soil_depth){
            defobj[0].active_zone_z = defobj[0].soil_depth;
		}

        if (fabs(defobj[0].active_zone_z - defobj[0].soil_depth) > 0.5) {
			printf("\nNote that soil depth used for biogeochem cycling (active zone z)");
 			printf("\nis more than 0.5 meter different from hydrologic soil depth");
 			printf("\nfor soil default file: %s\n", default_files[i] );
			}

        defobj[0].maximum_snow_energy_deficit = 	getDoubleParam(&paramCnt, &paramPtr, "maximum_snow_energy_deficit", "%lf", -10.0, 1);
        defobj[0].snow_water_capacity = 		getDoubleParam(&paramCnt, &paramPtr, "snow_water_capacity", "%lf", 0.0, 1);
        defobj[0].snow_light_ext_coef = 		getDoubleParam(&paramCnt, &paramPtr, "snow_light_ext_coef", "%lf", 10000.0, 1);
        defobj[0].snow_melt_Tcoef = 		getDoubleParam(&paramCnt, &paramPtr, "snow_melt_Tcoef", "%lf", 0.05, 1);
         printf("snow melting coef: %lf\n",defobj[0].snow_melt_Tcoef);
        defobj[0].snow_albedo_flag = 	parse_albedo_flag(getStrParam(&paramCnt, &paramPtr, "snow_albedo_flag", "%s", "age", 1));
        defobj[0].bats_b = 		getDoubleParam(&paramCnt, &paramPtr, "bats_b", "%lf", 2.0, 1);
        defobj[0].bats_r3 = 		getDoubleParam(&paramCnt, &paramPtr, "bats_r3", "%lf", 0.3, 1);
/*
             defobj[0].snow_melt_Tcoef *= command_line[0].tmp_value;
*/
        defobj[0].max_heat_capacity = 	getDoubleParam(&paramCnt, &paramPtr, "max_heat_capacity", "%lf", 0.0, 1);
        defobj[0].min_heat_capacity = 	getDoubleParam(&paramCnt, &paramPtr, "min_heat_capacity", "%lf", 0.0, 1);
        defobj[0].albedo = 		getDoubleParam(&paramCnt, &paramPtr, "albedo", "%lf", 0.28, 1);
        defobj[0].NO3_adsorption_rate =	getDoubleParam(&paramCnt, &paramPtr, "NO3_adsorption_rate", "%lf", 0.0, 1);
        defobj[0].N_decay_rate = 	getDoubleParam(&paramCnt, &paramPtr, "N_decay", "%lf", 0.12, 1);
		 // this is for controling co2 effect on stomatal conductance
        defobj[0].m_CO2_effect = getIntParam(&paramCnt, &paramPtr, "m_CO2_effect", "%d", 0, 1); //default is off
        printf("\n The CO2 effect on plants conductance is : %d", defobj[0].m_CO2_effect);
		/*
		if (command_line[0].tmp_value > ZERO)
            defobj[0].N_decay_rate *= command_line[0].tmp_value;
		*/
        defobj[0].soil_type.sand =		getDoubleParam(&paramCnt, &paramPtr, "sand", "%lf", 0.7, 1);
        defobj[0].soil_type.silt =		getDoubleParam(&paramCnt, &paramPtr, "silt", "%lf", 0.2, 1);
        defobj[0].soil_type.clay = 	getDoubleParam(&paramCnt, &paramPtr, "clay", "%lf", 0.1, 1);
        soil =  defobj[0].soil_type.sand
            + defobj[0].soil_type.silt
            + defobj[0].soil_type.clay;
        if  (fabs(soil - 1.0) > ZERO) {
			fprintf(stderr,
				"FATAL ERROR:in construct_soil_defaults\n proportion sand, silt, clay = %f\n\n", soil);
			printf("\n %d -  %f %f %f %f \n",
                defobj[0].ID,
                defobj[0].N_decay_rate,
                defobj[0].soil_type.sand,
                defobj[0].soil_type.silt,
                defobj[0].soil_type.clay);
		} /*end if*/
		if (command_line[0].gw_flag > 0) {
            defobj[0].sat_to_gw_coeff = getDoubleParam(&paramCnt, &paramPtr, "sat_to_gw_coeff", "%lf", 1.0, 1);
            defobj[0].sat_to_gw_coeff *= command_line[0].sat_to_gw_coeff_mult;
			}

		/*-----------------------------------------------------------------------------
		 *  Fill and Spill parameters
		 *-----------------------------------------------------------------------------*/
        defobj[0].fs_spill =	getDoubleParam(&paramCnt, &paramPtr, "fs_spill", "%lf", 1, 1);
        defobj[0].fs_percolation =	getDoubleParam(&paramCnt, &paramPtr, "fs_percolation", "%lf", 1, 1);
        defobj[0].fs_threshold = 	getDoubleParam(&paramCnt, &paramPtr, "fs_threshold", "%lf", 0.2, 1);

		/*--------------------------------------------------------------*/
		/*	vertical soil m and K are initized using soil default	*/
		/*	but sensitivity analysis -s is not applied to them	*/
		/*      use -sv to change these parameters			*/
		/*--------------------------------------------------------------*/
            defobj[0].m_v = defobj[0].m;
            defobj[0].mz_v = defobj[0].m_z;
            defobj[0].Ksat_0_v = defobj[0].Ksat_0;


		/*--------------------------------------------------------------*/
		/* sensitivity adjustment of vertical drainage  soil paramters	*/
		/*--------------------------------------------------------------*/
		if (command_line[0].vsen_flag > 0) {
                defobj[0].m_v *= command_line[0].vsen[M];
                defobj[0].mz_v *= command_line[0].vsen[M];
                defobj[0].Ksat_0_v *= command_line[0].vsen[K];
		}

		/*--------------------------------------------------------------*/
		/* sensitivity adjustment of soil drainage paramters		*/
		/*--------------------------------------------------------------*/
		if (command_line[0].sen_flag > 0) {
                defobj[0].m *= command_line[0].sen[M];
                defobj[0].m_z *= command_line[0].sen[M];
                defobj[0].Ksat_0 *= command_line[0].sen[K];
                defobj[0].soil_depth *= command_line[0].sen[SOIL_DEPTH];
		}

		/*--------------------------------------------------------------*/
		/*      calculate water_equivalent depth of soil                */
		/*--------------------------------------------------------------*/
        defobj[0].soil_water_cap = compute_delta_water(
            0, defobj[0].porosity_0,
            defobj[0].porosity_decay,
            defobj[0].soil_depth,
            defobj[0].soil_depth,
			0.0);

		/*--------------------------------------------------------------*/
		/* initialization of optional default file parms		*/
		/*--------------------------------------------------------------*/
        defobj[0].theta_mean_std_p1 = 	getDoubleParam(&paramCnt, &paramPtr, "theta_mean_std_p1", "%lf", 0.0, 1);
        defobj[0].theta_mean_std_p2 = 	getDoubleParam(&paramCnt, &paramPtr, "theta_mean_std_p2", "%lf", 0.0, 1);
        defobj[0].gl_c = 			getDoubleParam(&paramCnt, &paramPtr, "gl_c", "%lf", 0.0062, 1);
        defobj[0].gsurf_slope = 		getDoubleParam(&paramCnt, &paramPtr, "gsurf_slope ", "%lf", 0.01, 1);
        defobj[0].gsurf_intercept = 	getDoubleParam(&paramCnt, &paramPtr, "gsurf_intercept", "%lf", 0.001, 1);
        defobj[0].p4 = 			getDoubleParam(&paramCnt, &paramPtr, "p4", "%lf", -1.5, 1);
        defobj[0].DOM_decay_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DOM_decay_rate", "%lf", 0.05, 1);
        defobj[0].NH4_adsorption_rate =    getDoubleParam(&paramCnt, &paramPtr, "NH4_adsorption_rate", "%lf", 0.000005, 1);
        defobj[0].DON_production_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DON_production_rate", "%lf", 0.03, 1);
        defobj[0].DOC_adsorption_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DOC_adsorption_rate", "%lf", 0.000023, 1);
        defobj[0].DON_adsorption_rate = 	getDoubleParam(&paramCnt, &paramPtr, "DON_adsorption_rate", "%lf", 0.000001, 1);
        defobj[0].interval_size = 		getDoubleParam(&paramCnt, &paramPtr, "interval_size", "%lf", INTERVAL_SIZE, 1);
        defobj[0].decom_model = getIntParam(&paramCnt, &paramPtr, "decom_model", "%d", 1, 1); // 1 is default, 2 is FireBGC 3 is LandClim

		/*--------------------------------------------------------------*/
		/* sensitivity adjustment of vertical drainage  soil paramters	*/
		/* an  scale Pore size index and psi air entry or other parameters	*/
		/* that control moisture retention (if curve 3 is used)		*/
		/*--------------------------------------------------------------*/
		if (command_line[0].vsen_alt_flag > 0) {
            if (defobj[0].theta_psi_curve != 3)  {
                defobj[0].psi_air_entry *= command_line[0].vsen_alt[PA];
                defobj[0].pore_size_index *= command_line[0].vsen_alt[PO];
                if (defobj[0].pore_size_index >= 1.0) {
					printf("\n Sensitivity analysis giving Pore Size Index > 1.0, not allowed, setting to 1.0\n");
                    defobj[0].pore_size_index = 0.999;
					}
			}
			else {
                defobj[0].p3 *= command_line[0].vsen_alt[PA];
                defobj[0].p4 *= command_line[0].vsen_alt[PO];
			}
		}

		/*--------------------------------------------------------------*/
		/*      Fire effect parameters                          	*/
		/*--------------------------------------------------------------*/
        defobj[0].overstory_height_thresh = getDoubleParam(&paramCnt, &paramPtr, "overstory_height_thresh", "%lf", 6, 1);
        printf("\n the fire effect model overstory height thresh: %lf\n",defobj[0].overstory_height_thresh);
        defobj[0].understory_height_thresh = getDoubleParam(&paramCnt, &paramPtr, "understory_height_thresh", "%lf", 4, 1);
        printf("\n the fire effect model understory height thresh: %lf\n",defobj[0].understory_height_thresh);


        //10112022LML
        defobj[0].Dingman_coef = pow((defobj[0].psi_air_entry / 3.4),defobj[0].pore_size_index);
        defobj[0].soil_deficit_threshold = defobj[0].porosity_0 * defobj[0].porosity_decay
                                                        * (1. - exp(-defobj[0].soil_depth/defobj[0].porosity_decay))
                                                        *(1. - defobj[0].fs_threshold);

		/*--------------------------------------------------------------*/
		/*		Close the ith default file.								*/
		/*--------------------------------------------------------------*/

                memset(strbuf, '\0', strbufLen);
                strcpy(strbuf, default_files[i]);
                char *s = strbuf;
                char *y = NULL;
                char *token = NULL;
                char filename[256];

                // Store filename portion of path in 't'
                while ((token = strtok(s, "/")) != NULL) {
                    // Save the latest component of the filename
                    strcpy(filename, token);
                    s = NULL;
                }

                // Remove the file extension, if one exists
                memset(strbuf, '\0', strbufLen);
                strcpy(strbuf, filename);
                free(s);
                s = strbuf;
                token = strtok(s, ".");
                if (token != NULL) {
                    strcpy(filename, token);
                }

                memset(outFilename, '\0', filenameLen);



            // Concatenate the output prefix with the filename of the input .def file
            // and "_soil.params"
            if (command_line[0].output_prefix != NULL) {
                strcat(outFilename, command_line[0].output_prefix);
                if (filename != NULL) {
                    strcat(outFilename, "_");
                    strcat(outFilename, filename);
                }
                strcat(outFilename, "_soil.params");
            }
            else {
                if (filename != NULL) {
                    strcat(outFilename, "_");
                    strcat(outFilename, filename);
                }
                strcat(outFilename, "soil.params");
            }

        printParams(paramCnt, paramPtr, outFilename);
	} /*end for*/


	if (paramPtr != NULL)
            free(paramPtr);
  return(default_object_list);
} /*end construct_soil_defaults*/
