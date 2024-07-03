/*--------------------------------------------------------------*/
/* 																*/
/*						handle_event							*/
/*																*/
/*	handle_event - handles a tec event.							*/
/*																*/
/*	NAME														*/
/*	handle_event - handles a tec event.							*/
/*																*/
/*	SYNOPSIS													*/
/*	(void) handle_event(										*/
/*					struct	tec_entry	*event,					*/
/*					struct command_line_object *command_line,	*/
/*					struct	date	current_date,				*/
/*					struct world_object *world)					*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	This routine adjusts state variables based on the current	*/
/*	tec event.													*/
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*	This rotine is called by execute_tec.c just before the 		*/
/*	next tec entry is read.  This routine is called continually	*/
/*	with successive tec entries until the next tec entry has a	*/
/*	date later than the current tec entry.  In this manner 		*/
/*	multiple events can be executed before another simulation	*/
/*	time step.													*/
/*																*/
/*	March 6, 97 - C.Tague										*/
/*	-added a state output event 								*/
/* 	 current date added to handle_event parameter list			*/
/*																*/
/*	Sept, 98 - C.Tague	*/
/*	added comma delimited output event */
/*																*/
/*	June, 2014 - X Chen	*/
/*	added hourly growth output event*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rhessys.h"
bool isNumeric(const char *array);
bool isArrayBlank(const char *array);
bool read_MTBS_burnt_severity(struct world_object *world,struct command_line_object *command_line);
int set_dominant_severity(int *data);
void	handle_event(
					 struct	tec_entry	*event,
					 struct command_line_object *command_line,
					 struct	date	current_date,
					 struct world_object *world)
{
	/*--------------------------------------------------------------*/
	/*	Local Function Declarations.								*/
	/*--------------------------------------------------------------*/
	void	execute_redefine_strata_event(
		struct world_object *,
		struct command_line_object *,
		struct	date);
	void	execute_redefine_world_event(
		struct world_object *,
		struct command_line_object *,
		struct	date);
	void	execute_redefine_world_mult_event(
		struct world_object *,
		struct command_line_object *,
		struct	date);
	void	execute_redefine_world_thin_event(
		struct world_object *,
		struct command_line_object *,
		struct	date,
		int);
	void	execute_road_construction_event(
		struct world_object *,
		struct command_line_object *,
		struct	date);
	void	execute_state_output_event(
		struct world_object *,
		struct date,
		struct date,
		struct command_line_object *);
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Below is a list of events which are caught if present.		*/
	/*	OTHERWISE a fatal error ensues.								*/
	/*--------------------------------------------------------------*/

	if ( !strcmp(event[0].command,"none") ){
		/* nothing here */
	}
	else if ( !strcmp(event[0].command,"print_yearly_on") ){
		command_line[0].output_flags.yearly= 1;
		command_line[0].output_yearly_date.month = current_date.month;
		command_line[0].output_yearly_date.day = current_date.day;
	}
	else if ( !strcmp(event[0].command,"print_yearly_off") ){
		command_line[0].output_flags.yearly= 0;
	}
	else if ( !strcmp(event[0].command,"print_yearly_growth_on")){
		command_line[0].output_flags.yearly_growth = 1;
		command_line[0].output_yearly_date.month = current_date.month;
		command_line[0].output_yearly_date.day = current_date.day;
	}
	else if (  !strcmp(event[0].command,"print_yearly_growth_off")){
		command_line[0].output_flags.yearly_growth = 0;
	}
	else if ( !strcmp(event[0].command,"print_monthly_on") ){
		command_line[0].output_flags.monthly= 1;
	}
	else if ( !strcmp(event[0].command,"print_monthly_off") ){
		command_line[0].output_flags.monthly= 0;
	}
	else if ( !strcmp(event[0].command,"print_daily_on") ){
		command_line[0].output_flags.daily= 1;
	}
	else if ( !strcmp(event[0].command,"print_daily_off") ){
		command_line[0].output_flags.daily= 0;
	}
	else if ( !strcmp(event[0].command,"print_daily_growth_on")){
		command_line[0].output_flags.daily_growth = 1;
	}
	else if ( !strcmp(event[0].command,"print_daily_csv_growth_on")){
		command_line[0].output_flags.daily_growth = 1;
	}
	else if ( !strcmp(event[0].command,"print_daily_csv_on") ){
		command_line[0].output_flags.daily = 1;
	}
	else if ( !strcmp(event[0].command,"print_yearly_csv_on") ){
		command_line[0].output_flags.yearly = 1;
		command_line[0].output_yearly_date.month = current_date.month;
		command_line[0].output_yearly_date.day = current_date.day;
	}
	else if ( !strcmp(event[0].command,"print_monthly_csv_on") ){
		command_line[0].output_flags.monthly = 1;
	}
	else if ( !strcmp(event[0].command,"print_daily_growth_off")){
		command_line[0].output_flags.daily_growth = 0;
	}
	else if ( !strcmp(event[0].command,"print_daily_csv_growth_off")){
		command_line[0].output_flags.daily_growth = 0;
	}
	else if ( !strcmp(event[0].command,"print_daily_csv_off") ){
		command_line[0].output_flags.daily = 0;
	}
	else if ( !strcmp(event[0].command,"print_yearly_csv_off") ){
		command_line[0].output_flags.yearly = 0;
	}
	else if ( !strcmp(event[0].command,"print_monthly_csv_off") ){
		command_line[0].output_flags.monthly = 0;
	}
	else if ( !strcmp(event[0].command,"print_hourly_on") ){
		command_line[0].output_flags.hourly= 1;
	}
	else if ( !strcmp(event[0].command,"print_hourly_off") ){
		command_line[0].output_flags.hourly= 0;
	}
	else if ( !strcmp(event[0].command,"print_hourly_growth_on")){
		command_line[0].output_flags.hourly_growth = 1;
	}
	else if ( !strcmp(event[0].command,"print_hourly_growth_off")){
		command_line[0].output_flags.hourly_growth = 0;
	}
//#ifdef LIU_BURN_ALL_AT_ONCE
    else if ( !strcmp(event[0].command,"burn_on")){
        command_line[0].burn_on_flag = 1;
        if(isNumeric(event[0].option)) {
           printf("option(%s)) is a number\n",event[0].option);
           int severity = atoi(event[0].option);
           if (severity != BURNT_SEVERITY_USE_COMMAND_ARGUMENTS &&
                           ((severity < MTBS_BURNT_SEVERITY_BACKGROUND)
                             || (severity >= MTBS_BURNT_SEVERITY_NUMCLASSES)))
           {
               printf("Warning: severity (%s) is out of range and use argument for mortality.\n",event[0].option);
               command_line[0].burn_on_severity = BURNT_SEVERITY_USE_COMMAND_ARGUMENTS;                        //Use specific command arguments to define severity
           } else {
               command_line[0].burn_on_severity = severity;
           }
        } else if (isArrayBlank(event[0].option)){
           printf("option(%s)) is blank\n",event[0].option);
           command_line[0].burn_on_severity = BURNT_SEVERITY_USE_COMMAND_ARGUMENTS;
        } else {
           command_line[0].burn_on_severity = BURNT_SEVERITY_USE_BURNT_SEVERITY_MAP;
           strcpy(command_line[0].burn_on_severity_filename,event[0].option);
           read_MTBS_burnt_severity(world,command_line);
        }
    }
    else if ( !strcmp(event[0].command,"burn_off")){
        command_line[0].burn_on_flag = 0;
    }
//#endif
	else if ( !strcmp(event[0].command,"output_current_state") ){
		execute_state_output_event(world, current_date,
			world[0].end_date,command_line);
	}
	else if ( !strcmp(event[0].command,"redefine_strata") ){
		execute_redefine_strata_event(world, command_line, current_date);
	}
	else if ( !strcmp(event[0].command,"redefine_world") ){
		execute_redefine_world_event(world, command_line, current_date);
	}
	else if ( !strcmp(event[0].command,"redefine_world_multiplier") ){
		execute_redefine_world_mult_event(world, command_line, current_date);
	}	
	else if ( !strcmp(event[0].command,"redefine_world_thin_remain") ){
		execute_redefine_world_thin_event(world, command_line, current_date, 1);
	}		
	else if ( !strcmp(event[0].command,"redefine_world_thin_harvest") ){
		execute_redefine_world_thin_event(world, command_line, current_date, 2);
	}
	else if ( !strcmp(event[0].command,"redefine_world_thin_snags") ){
		execute_redefine_world_thin_event(world, command_line, current_date, 3);
	}			
	else if ( !strcmp(event[0].command,"roads_on") ){
		command_line[0].road_flag = 1;
		execute_road_construction_event(world, command_line, current_date);
	}
	else if ( !strcmp(event[0].command,"roads_off") ){
		command_line[0].road_flag = 0;
	}
	else{
		fprintf(stderr,"FATAL ERROR: in handle event - event %s not recognized.\n",
			event[0].command);
		exit(EXIT_FAILURE);
	}
} /*end handle_event*/

/*01092024LML from ChatGPT*/
bool isNumeric(const char *array) {
    int i = 0;
    // Check for optional sign character (+/-) at the beginning
    if (array[0] == '+' || array[0] == '-') {
        i = 1; // Skip the sign character
    }
    bool digitEncountered = false; // Flag to track at least one digit in the array
    // Check each character to ensure it is a digit
    while (array[i] != '\0') {
        if (isdigit(array[i])) {
            digitEncountered = true;
            i++;
        } else {
            return false; // If any non-digit character is found, return false
        }
    }
    return digitEncountered; // Return true only if at least one digit was encountered
}
/*01092024LML from ChatGPT*/
bool isArrayBlank(const char *array) {
    int i = 0;
    // Check each character in the array
    while (array[i] != '\0') {
        if (!isspace(array[i])) {
            return false; // If any non-whitespace character is found, return false
        }
        i++;
    }
    return true; // Return true if all characters are whitespace
}
/*01092024LML read MTBS burnseverity*/
bool read_MTBS_burnt_severity(struct world_object *world,struct command_line_object *command_line) {
    FILE *burnt_severity = fopen(command_line[0].burn_on_severity_filename,"r");
    if (burnt_severity != NULL) {
      //loops all patches
      for (int basin = 0; basin < world[0].num_basin_files; basin++) {
        #pragma omp parallel for
        for (int h = 0 ; h < world[0].basins[basin]->num_hillslopes; h++) {
          for (int z = 0 ; z < world[0].basins[basin]->hillslopes[h]->num_zones; z++) {
            for (int p = 0; p < world[0].basins[basin]->hillslopes[h]->zones[z]->num_patches; p++) {
              struct patch_object *pt = world[0].basins[basin]->hillslopes[h]->zones[z]->patches[p];
              for (int bs = 0; bs < MTBS_BURNT_SEVERITY_NUMCLASSES; bs++) {
                 pt->burnt_severity_acount[bs] = 0;
              }
              pt->dominant_burnt_severity_type = MTBS_BURNT_SEVERITY_BACKGROUND;
            }
          }
        }
      }

      for(int i=0; i<world[0].defaults[0].fire->n_rows; i++){
        for(int j=0; j<world[0].defaults[0].fire->n_cols; j++){
          int tmp_burnt_severity = MTBS_BURNT_SEVERITY_BACKGROUND;
          fscanf(burnt_severity,"%d\t",&tmp_burnt_severity);

          struct patch_fire_object *patch_fire = &world->patch_fire_grid[i][j];
          if (patch_fire->num_patches > 0) {
            struct patch_object *patch = patch_fire->patches[0];
            if (tmp_burnt_severity >= 0) patch_fire->MTBS_burnt_severity = tmp_burnt_severity;
            else patch_fire->MTBS_burnt_severity = MTBS_BURNT_SEVERITY_BACKGROUND;
              patch->burnt_severity_acount[patch_fire->MTBS_burnt_severity]++;
          }
        }
      }
      //loops all patches
      for (int basin = 0; basin < world[0].num_basin_files; basin++) {
        #pragma omp parallel for
        for (int h = 0 ; h < world[0].basins[basin]->num_hillslopes; h++) {
          for (int z = 0 ; z < world[0].basins[basin]->hillslopes[h]->num_zones; z++) {
              for (int p = 0; p < world[0].basins[basin]->hillslopes[h]->zones[z]->num_patches; p++) {
                  struct patch_object *pt = world[0].basins[basin]->hillslopes[h]->zones[z]->patches[p];
                  pt->dominant_burnt_severity_type = get_dominant_severity(pt->burnt_severity_acount);
                  if (world[0].defaults[0].fire[0].fire_verbose == 3)
                    printf("b:%d\th:%d\tz:%d\tp:%d\tseverity:%d\n",world[0].basins[basin]->ID
                          ,world[0].basins[basin]->hillslopes[h]->ID
                          ,world[0].basins[basin]->hillslopes[h]->zones[z]->ID
                          ,pt->ID,pt->dominant_burnt_severity_type);
              }
          }
        }
      }
      //Get the dominant severity type
      return 0;
    }
    return -1;
}
//get dominant type
int get_dominant_severity(int *data){
  int maximum_acount = 0;
  int dominant_type = MTBS_BURNT_SEVERITY_BACKGROUND;
  for (int i = 0; i < MTBS_BURNT_SEVERITY_NUMCLASSES; i++) {
      if (data[i] >= maximum_acount && i != MTBS_BURNT_SEVERITY_NONMAPPING) {
          maximum_acount = data[i];
          dominant_type = i;
      }
  }
  if (maximum_acount == 0) dominant_type = MTBS_BURNT_SEVERITY_BACKGROUND;
  return dominant_type;
}
