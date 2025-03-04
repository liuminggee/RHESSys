/*--------------------------------------------------------------*/
/*                                                              */
/*      construct_beetle_defaults                               */
/*                                                              */
/*      construct_beetle_defaults.c - makes beetle default      */
/*                              model objects.                  */
/*                                                              */
/*      NAME                                                    */
/*      construct_beetle_defaults.c - makes beetle default      */
/*                                    objects.                  */
/*                                                              */
/*      SYNOPSIS                                                */
/*      struct beetle_default *construct_beetle_defaults(       */
/*                        num_default_files,                    */
/*                        default_files,                        */
/*                        beetlespread_flag,                    */
/*                                                              */
/*      OPTIONS                                                 */
/*                                                              */
/*      DESCRIPTION                                             */
/*      NR 2018-4-19                                            */
/*                                                              */
/*      PROGRAMMER NOTES                                        */
/*                                                              */
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "params.h"

struct spinup_default *construct_beetle_defaults(
        int     num_default_files,
        char    **default_files,
        struct command_line_object *command_line)

{
        /*--------------------------------------------------------------*/
        /*      Local function definition.                              */
        /*--------------------------------------------------------------*/
        void    *alloc( size_t, char *, char *);

        /*--------------------------------------------------------------*/
        /*      Local variable definition.                              */
        /*--------------------------------------------------------------*/
        int     i;
        int strbufLen = 256;
        int filenameLen = 1024;
        int paramCnt = 0;
        char    strbuf[strbufLen];
        char    outFilename[filenameLen];
        double  ftmp, soil;
        FILE    *default_file;
        char    *newrecord;
        char    record[MAXSTR];
        struct  beetle_default    *default_object_list;
        param *paramPtr = NULL;

        /*--------------------------------------------------------------*/
        /*      Allocate an array of default objects.                   */
        /*--------------------------------------------------------------*/
        default_object_list   = (struct beetle_default *)
                alloc(num_default_files *
                sizeof(struct beetle_default),"default_object_list",
                "construct_beetle_defaults");


        /*--------------------------------------------------------------*/
        /*      Loop through the default files list.                    */
        /*--------------------------------------------------------------*/
        for (i=0 ; i<num_default_files; i++){
                /*--------------------------------------------------------------*/
                /*              Try to open the ith default file.               */
                /*--------------------------------------------------------------*/
                printf("\nReading %s\n", default_files[i]);
                paramCnt = 0;
                if (paramPtr != NULL)
                    free(paramPtr);

                paramPtr = readParamFile(&paramCnt, default_files[i]);

                /*--------------------------------------------------------------*/
                /*              read the ith default file into the ith object   */
                /*--------------------------------------------------------------*/
                default_object_list[i].ID = getIntParam(&paramCnt, &paramPtr, "beetle_default_ID", "%d", 1, 1); // new param name
                /*--------------------------------------------------------------*/
                /*              assign parameters in  default and read the      */
                /*   optional parameter specification                           */
                /*--------------------------------------------------------------*/
                default_object_list[i].attack_mortality = getDoubleParam(&paramCnt, &paramPtr, "attack_mortality", "%lf", 0.95, 1); // the default mortality is 95%
		        printf("mortality: %lf\n",default_object_list[i].attack_mortality);
                default_object_list[i].year_delay = getDoubleParam(&paramCnt, &paramPtr, "year_delay", "%lf", 5, 1); // the default delay year is 5
		        printf("snag_year_delay: %d\n",default_object_list[i].year_delay);
                default_object_list[i].half_life = getDoubleParam(&paramCnt, &paramPtr, "half_life", "%lf", 10, 1); // the default delay year is 10
		        printf("snag_half_life: %d\n",default_object_list[i].half_life);
                default_object_list[i].n_rows=getIntParam(&paramCnt, &paramPtr, "n_rows", "%d", -1, 1);
                printf("n_rows: %d\n",default_object_list[i].n_rows);
                default_object_list[i].n_cols=getIntParam(&paramCnt, &paramPtr, "n_cols", "%d", -1, 1);
                printf("n_cols: %d\n",default_object_list[i].n_cols);
                default_object_list[i].year_attack=getIntParam(&paramCnt, &paramPtr, "year_attack", "%d", 1991, 1);
                printf("year_attack: %d\n",default_object_list[i].year_attack);
                default_object_list[i].beetle_in_buffer =getIntParam(&paramCnt, &paramPtr, "beetle_in_buffer", "%d", 0, 1);
                printf("beetle_in_buffer: %d\n",default_object_list[i].beetle_in_buffer);
                default_object_list[i].leaf_year_delay = getDoubleParam(&paramCnt, &paramPtr, "leaf_year_delay", "%lf", 3, 1); // the default delay year is 3
		        printf("leaf_year_delay: %d\n",default_object_list[i].leaf_year_delay);
                default_object_list[i].leaf_half_life = getDoubleParam(&paramCnt, &paramPtr, "leaf_half_life", "%lf", 2, 1); // the default delay year is 2
		        printf("leaf_half_life: %d\n",default_object_list[i].leaf_half_life);
		        default_object_list[i].calc_single_attack= getDoubleParam(&paramCnt, &paramPtr, "calc_single_attack", "%lf", 0, 1); // the default is not calculate single attack
		        printf("calculate single attack: %d\n",default_object_list[i].calc_single_attack);
		        default_object_list[i].min_abc= getDoubleParam(&paramCnt, &paramPtr, "min_abc", "%lf", 0.01, 1); // the min carbon g/m2 each patch have that attack happens before is 1000 change to 100 on 20190910
		        printf("the minimum carbon (g/m2..) in the patch that an attack will happen : %lf \n",default_object_list[i].min_abc);
		        default_object_list[i].mortality_type = getIntParam(&paramCnt, &paramPtr, "mortality_type", "%d", 1, 1); //type 1 is beetle type 2 is prescribed fire NR 2019/04/30
		        printf("the mortality type is: %d\n", default_object_list[i].mortality_type);
		        // the root go to litter pool or not NREN 20190908
		        default_object_list[i].root_alive= getIntParam(&paramCnt, &paramPtr, "root_alive", "%d", 1, 1); //1 is alive no touching of root after beetle attack 0 is root is dead, 2 is dead root go to decay pool combine with the half life below
		        printf("the root_alive after attack is: %d\n", default_object_list[i].root_alive);
		        default_object_list[i].harvest_dead_root= getIntParam(&paramCnt, &paramPtr, "harvest_dead_root", "%d", 1, 1); //1 is the dead root is harvested not go to litter pool,0 is dead root go to litter pool
		        printf("if harvest_dead_root is: %d\n", default_object_list[i].harvest_dead_root);
                default_object_list[i].deadroot_half_life = getDoubleParam(&paramCnt, &paramPtr, "deadroot_half_life", "%lf", 3, 1); // the default decay half life of deadroot is 3 years
		        printf("deadroot_half_life: %d\n",default_object_list[i].deadroot_half_life); //NREN 20190911
		        default_object_list[i].num_snag_sequence = getIntParam(&paramCnt, &paramPtr, "num_snag_sequence", "%d", 300, 1); //NREN this for dynamic memory allocation 20190914
		        printf("num_snag_sequence for memory allocation: %d\n",default_object_list[i].num_snag_sequence); //NREN 20190911
		        default_object_list[i].transfer_root_water = getIntParam(&paramCnt, &paramPtr, "transfer_root_water", "%d", 0, 1); //NREN this is for transfer root water is the rootdepth is changed
		        printf("need to transfer root water if the root depth is changed: %d\n",default_object_list[i].transfer_root_water); //NREN 20190914
		        // improve the beetle effect model to incoorperate with the harvest the dead stem and leaf
		        // if harvest leafc and harvest stem is true, here only control the lai and pai; you also need to change the delay time of stem and leaf to 100000
		        // at the same time, the harvest dead root should be true and root_alive is 0 too.
		        // if not harvest the leafc and stem, it will use the leaf_when_red and stem_when_red to calculate the LAI and PAI, and rootalive =0, harvest dead root is 0
		        default_object_list[i].lai_include_redneedle= getIntParam(&paramCnt, &paramPtr, "lai_include_redneedle", "%d", 0, 1); //1 is alive no touching of root after beetle attack 0 is root is dead
		        printf("calculating lai considering the redneedle: %d\n", default_object_list[i].lai_include_redneedle);
		        default_object_list[i].pai_include_snag = getIntParam(&paramCnt, &paramPtr, "pai_include_snag", "%d", 0, 1); //1 is alive no touching of root after beetle attack 0 is root is dead
		        printf("calculating pai considering the snag: %d\n", default_object_list[i].pai_include_snag);
		        default_object_list[i].height_include_snag = getIntParam(&paramCnt, &paramPtr, "height_include_snag", "%d", 1, 1); //1 is alive no touching of root after beetle attack 0 is root is dead
		        printf("calculating height considering the snag: %d\n", default_object_list[i].height_include_snag);




                /*--------------------------------------------------------------*/
                /*              Close the ith default file.                     */
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
                // and "_stratum.params"
                if (command_line[0].output_prefix != NULL) {
                    strcat(outFilename, command_line[0].output_prefix);
                    if (filename != NULL) {
                        strcat(outFilename, "_");
                        strcat(outFilename, filename);
                    }
                    strcat(outFilename, "_beetle.params");
                }
                else {
                    if (filename != NULL) {
                        strcat(outFilename, "_");
                        strcat(outFilename, filename);
                    }
                    strcat(outFilename, "beetle.params");
                }

            printParams(paramCnt, paramPtr, outFilename);
        } /*end for*/
        return(default_object_list);
} /*end construct_beetle_defaults*/



