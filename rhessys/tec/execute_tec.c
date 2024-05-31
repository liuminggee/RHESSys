/*--------------------------------------------------------------*/
/* 																*/
/*						execute_tec_file						*/
/*																*/
/*	execute_tec_file.c - executes the tec file with the world	*/
/*																*/
/*	NAME														*/
/*	execute_tec_file.c - executes the tec file with the world	*/
/*																*/
/*	SYNOPSIS													*/
/*	void *execute_tec(											*/
/*			struct	tec_object,									*/
/*			struct command_line_object,					 		*/
/*			struct	world_output_file_object,					*/
/*			struct	world_output_file_object,					*/
/*			struct world_object );								*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	To encourage future parellelism, operations on groups of	*/
/*	objects are grouped together (e.g. all canopy processes		*/
/*	for a given patch are put together) in an event loop.		*/
/*																*/
/*	The event loop cycles through an a priori temporal and		*/
/*	spatial hierarchy.  At the moment the hierarchies are:		*/
/*																*/
/*	SPATIAL						TEMPORAL						*/
/*	world						yearly							*/
/*	basin						monthly							*/
/*	hillslope					daily							*/
/*	zone						hourly							*/
/*	patch														*/
/*	stratum 													*/
/*																*/
/*	These levels are somewhat arbitrary but changes should be	*/
/*	documented.  Refer to the rhessys.h file for documentation 	*/
/*	of the spatial object hierarchy.							*/
/*																*/
/*	The temporal levels are handeled as implicit events.		*/
/*	For simplicity no sub-daily events are handled at the   	*/
/*	moment.														*/
/*																*/
/*	Note: It is important to start at the beginning of a day	*/
/*			and end before the beginning of a day since there	*/
/*			are day start and day end simulation calls.  So, 	*/
/*			make sure that the start and end dates have hour =1	*/
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	There is only one tec event that is kept track of.  		*/
/*	This is the next event.  									*/
/*																*/
/*	This implementation of the event loop allows a comprimise 	*/
/*	between a clean spatial and temporal layout of processes	*/
/*	and the desire to reduce RAM allocation.					*/
/*																*/
/*	The event loop starts at the world start date and continues	*/
/*	up to , but not including , the world end date.  For the	*/
/*	purposes of this code, the first hour in a day is hour 1 and*/
/*	the last hour is hour 24.  									*/
/*																*/
/*	The event loop has an hourly time step in which is does the	*/
/*	following:													*/
/*																*/
/*	1.  Initialze the event queue with a null event and set the	*/
/*			date to the start date.								*/
/*	2.  check if the current date is less than the end date		*/
/*		(if not the the event loop ends.						*/
/*	3.  	Perform the next event in the event queue.			*/
/*	4.  	Read the next event from tec file into event queue.	*/
/*			(if there are no more events set the next event to	*/
/*			 a null event).										*/
/*	5.  	check if  current date is less than the date of the	*/
/*			next event (if not go to 13).						*/
/*	6.			if hour = 1 call world_daily_I 			        */
/*	7.			call world_hourly								*/
/*	8.  		if hour = 24 	call world_daily_F				*/
/*								increment the day by one		*/
/*	9.			if the month of the new date does not equal that*/
/*					of the last day's date increment the month 	*/
/*	10.			if the year of the new date does not equal that*/
/*					of the last day's date increment the year 	*/
/*	11. 		increment the hour								*/
/*	12. 		goto 5											*/
/*	13.		At the moment there is no event handling here.		*/
/*	14.		goto 2												*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
#include "functions.h"

#ifndef MAX_LINE_LENGTH
#define MAX_LINE_LENGTH 1024
#endif
int numer_sections_char_array(char* line);
void	execute_tec(
					struct	tec_object *tecfile ,
					struct command_line_object *command_line,
					struct	world_output_file_object *outfile,
					struct	world_output_file_object *growth_outfile,
					struct world_object *world)
{
	/*--------------------------------------------------------------*/
	/*	Local Function Declarations.								*/
	/*--------------------------------------------------------------*/
	int		cal_date_lt(struct date, struct date );

	long	julday( struct date );

	struct	date	caldat( long );

	struct	tec_entry	*construct_tec_entry( struct date, char * );

	void	world_daily_I(
		long,
		struct world_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);

	void	world_hourly(
		struct world_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);

	void	world_daily_F(
		long,
		struct world_object *,
		struct command_line_object *,
		struct tec_entry *,
		struct date);

	void	handle_event(
		struct	tec_entry	*,
		struct	command_line_object	*,
		struct	date,
		struct	world_object	*);

	void	execute_yearly_growth_output_event(
		struct	world_object	*,
		struct	command_line_object	*,
		struct	date,
		struct	world_output_file_object *);

	void	execute_yearly_output_event(
		int,
		struct	world_object	*,
		struct	command_line_object	*,
		struct	date,
		struct	world_output_file_object *);

	void	execute_daily_output_event(
		struct	world_object	*,
		struct	command_line_object	*,
		struct	date,
		struct	world_output_file_object *);

	void	execute_daily_growth_output_event(
		struct	world_object	*,
		struct	command_line_object	*,
		struct	date,
		struct	world_output_file_object *);

	void	execute_monthly_output_event(
		struct	world_object	*,
		struct	command_line_object	*,
		struct	date,
		struct	world_output_file_object *);

	void	execute_hourly_output_event(
		struct	world_object	*,
		struct	command_line_object	*,
		struct	date,
		struct	world_output_file_object *);

	void	execute_firespread_event(
		struct	world_object	*,
		struct	command_line_object	*,
		struct	date);

	void    execute_beetlespread_event(
						world,
						command_line,
						current_date


	);
	void	execute_state_output_event(
		struct world_object *,
		struct date,
		struct date,
		struct command_line_object *);

	/*--------------------------------------------------------------*/
	/*	Local Variable Definition. 									*/
	/*--------------------------------------------------------------*/
	int check;
        int reset_flag;
	long	day;
	long	hour;
	long	month;
	long	year;
	struct	date	current_date;
	struct	date	next_date;
	struct	tec_entry	*event;

#ifdef JMG_TRACKING
    world[0].track_simtime.syr = 1;
    world[0].track_simtime.smth = 1;
    world[0].track_simtime.sday = 1;
#endif

	/*--------------------------------------------------------------*/
	/*	Initialize the indices into the base station clime sequences*/
	/*--------------------------------------------------------------*/
	reset_flag = 0.0;
	year = 0;
	month = 0;
	day = 0;
	hour = 0;

	/*--------------------------------------------------------------*/
	/*	Initialize the tec event									*/
	/*--------------------------------------------------------------*/

	event =  construct_tec_entry(world[0].end_date,"none");

	/*--------------------------------------------------------------*/
	/*	Loop from the start of the world to the end of the world.	*/
	/*--------------------------------------------------------------*/
	current_date = world[0].start_date;
	next_date = current_date;
//#ifdef LIU_BURN_ALL_AT_ONCE
    int spins_index = 0;                                                    //which spin
    int spin_year_index = 0;                                                //which year index within this spin
//#endif


	while ( cal_date_lt(current_date,world[0].end_date)){
		/*--------------------------------------------------------------*/
		/*		Perform the tec event.									*/
		/*--------------------------------------------------------------*/
        handle_event(event,command_line,current_date,world);
//#ifdef LIU_BURN_ALL_AT_ONCE
        printf("Event: year:%d mon:%d day:%d command:%s option:%s\n",event[0].cal_date.year,
                event[0].cal_date.month,event[0].cal_date.day,event[0].command,event[0].option);
//#endif
		/*--------------------------------------------------------------*/
		/*		read the next tec file entry.							*/
		/*		if we are not at the end of the tec file.				*/
		/*--------------------------------------------------------------*/
		if ( !(feof(tecfile[0].tfile))){
			/*--------------------------------------------------------------*/
			/*			read in the next tec line.							*/
			/*--------------------------------------------------------------*/
            char line[MAX_LINE_LENGTH];
            char *token;
            // Read a line from the file
            if (fgets(line, MAX_LINE_LENGTH, tecfile[0].tfile) == NULL) {
              printf("End of file reached or an error occurred while reading.\n");
            }
            //printf("Tokens: %d\n", numer_sections_char_array(line));
            if (numer_sections_char_array(line) == 5) {
              //printf("line:%s\n",line);
              memset(event[0].option, '\0', sizeof(event[0].option));
              //check = fscanf(tecfile[0].tfile,"%d %d %d %d %s\n",
              check = sscanf(line,"%d %d %d %d %s\n",
				&(event[0].cal_date.year),
				&(event[0].cal_date.month),
				&(event[0].cal_date.day),
				&(event[0].cal_date.hour),
				event[0].command);
            } else if (numer_sections_char_array(line) == 6) {
                check = sscanf(line,"%d %d %d %d %s %s\n",
                  &(event[0].cal_date.year),
                  &(event[0].cal_date.month),
                  &(event[0].cal_date.day),
                  &(event[0].cal_date.hour),
                  event[0].command,
                  event[0].option);
            }
			/*--------------------------------------------------------------*/
			/*			report an error if for some reason we cant read it	*/
			/*--------------------------------------------------------------*/
			if ( !check ){
				fprintf(stderr,"\nERROR:  the tec file is corrupted.");
				fclose(tecfile[0].tfile);
				exit(EXIT_FAILURE);
			} /*end if*/
		} /*end if*/
		/*--------------------------------------------------------------*/
		/* 	if end of tec file next event is the end of the world		*/
		/*--------------------------------------------------------------*/
		else{
			event =  construct_tec_entry(world[0].end_date, "none");
		} /*end if-else*/
		/*--------------------------------------------------------------*/
		/*		If the next event's date exceeds the end_date then		*/
		/*		set the next event to nothing and at a time at the 		*/
		/*		end of the simulation.									*/
		/*--------------------------------------------------------------*/
		if ( cal_date_lt(event[0].cal_date,	world[0].end_date) == 0  ){
			event =  construct_tec_entry(world[0].end_date,	"none");
		} /*end if*/
		/*--------------------------------------------------------------*/
		/*		Do stuff until the next tec event.						*/
		/*--------------------------------------------------------------*/
		while ( cal_date_lt(current_date, event[0].cal_date)){
			/*--------------------------------------------------------------*/
			/*			Simulate the world for the start of this day e		*/
			/*--------------------------------------------------------------*/
#ifdef LIU_DISPLY_RUN_INFO
            if ( current_date.hour == 1 ) printf("Current_date runyear = %d clim_year = %d mon = %d day = %d\n",
                   year,current_date.year,current_date.month,current_date.day); /*ning*/
#endif
            //fflush(stdout);
			if ( current_date.hour == 1 ){
                world_daily_I(
					day,
					world,
					command_line,
					event,
                    current_date);
			} /*end if*/
			/*--------------------------------------------------------------*/
			/*          Do hourly stuff for the day.                        */
			/*--------------------------------------------------------------*/
            world_hourly( world,
				command_line,
				event,
                current_date);

			/*--------------------------------------------------------------*/
			/*			Perform any requested hourly output					*/
			/*--------------------------------------------------------------*/
			if (command_line[0].output_flags.hourly == 1){
				execute_hourly_output_event(
							  world,
							  command_line,
							  current_date,
							  outfile);
			}

			if(command_line[0].output_flags.hourly_growth ==1 &&
					(command_line[0].grow_flag > 0) ){
				  execute_hourly_growth_output_event(
							      world,
							      command_line,
							      current_date,
							      growth_outfile);

				};
			/*--------------------------------------------------------------*/
			/*			Increment to the next hour.							*/
			/*--------------------------------------------------------------*/
			current_date.hour++;
			/*--------------------------------------------------------------*/
			/*			Check if this is a day end.							*/
			/*--------------------------------------------------------------*/
			if ( current_date.hour == 25 ){                
				/*--------------------------------------------------------------*/
				/*			Simulate the world for the end of this day e		*/
				/*--------------------------------------------------------------*/
//#ifdef JMG_TRACKING
//        /* Update simulation day count */
//                world[0].track_simtime.sday += 1;
//#endif
                world_daily_F(
					day,
					world,
					command_line,
					event,
                    current_date);
				/*--------------------------------------------------------------*/
				/*			Perform any requested daily output					*/
				/*--------------------------------------------------------------*/
                if ((command_line[0].output_flags.daily_growth == 1) &&
							(command_line[0].grow_flag > 0) ) {
						execute_daily_growth_output_event(
						world,
						command_line,
						current_date,
						growth_outfile);
				}
				if (command_line[0].output_flags.daily == 1) {
						execute_daily_output_event(
						world,
						command_line,
						current_date,
						outfile);
                               }
#ifdef LIU_WMFIRE_OUTPUT
                //05242022LML always print daily basin output
                if ((command_line[0].output_flags.daily != 1) || (command_line[0].b == NULL)) {
                    for (int b=0; b < world[0].num_basin_files; ++ b ) {
                            int basinID = command_line[0].b->basinID;
                            if (( world[0].basins[b][0].ID == basinID) || (basinID == -999)) {
                                output_basin(
                                    command_line[0].routing_flag,
                                    world[0].basins[b],
                                    current_date,
                                    outfile->basin->daily);
                                if (command_line[0].grow_flag != 0)
                                output_growth_basin(
                                    world[0].basins[b],
                                    current_date,
                                    growth_outfile->basin->daily
            #ifdef JMG_TRACKING
                                    ,&world[0].track_simtime
            #endif
                                    );
                            }
                    }
                }
#endif
				/*--------------------------------------------------------------*/
        /*  Output world state in spinup mode if targets met            */
				/*--------------------------------------------------------------*/

				if((command_line[0].vegspinup_flag > 0) && (world[0].target_status > 0)) {
		      execute_state_output_event(world, current_date, world[0].end_date,command_line);
          printf("\nSpinup completed YEAR %d MONTH %d DAY %d \n", current_date.year,current_date.month,current_date.day);
          exit(0);
        }

				/*--------------------------------------------------------------*/
				/*			Perform any requested yearly output					*/
				/*--------------------------------------------------------------*/
				 if (command_line[0].output_flags.yearly_growth == 1) {reset_flag=0;}

                //06192023LML note: this yearly output is conducted end of day and might cause one day's shift
				if ((command_line[0].output_flags.yearly == 1) &&
					(command_line[0].output_yearly_date.month==current_date.month)&&
                    (command_line[0].output_yearly_date.day == current_date.day)) {

                    //printf("current_date year%d month:%d day:%d\n",
                    //        current_date.year,
                    //        command_line[0].output_yearly_date.month,
                    //        command_line[0].output_yearly_date.day);

							execute_yearly_output_event(
							reset_flag,
							world,
							command_line,
							current_date,
							outfile);
                }

				if ((command_line[0].output_flags.yearly_growth == 1) &&
					(command_line[0].output_yearly_date.month==current_date.month)&&
					(command_line[0].output_yearly_date.day == current_date.day) &&
					(command_line[0].grow_flag > 0) )
					execute_yearly_growth_output_event(
					world,
					command_line,
					current_date,
					growth_outfile);
				/*--------------------------------------------------------------*/
				/*				Determine the new calendar date if we add 1 day.*/
				/*				Do this by first conversting the current cal	*/
				/* 				endar date into a julian day.  Then adding one	*/
				/*				to the julian day and the converting back to	*/
				/*				get tomorrows calendar date.					*/
				/*				We assume that it starts at hour 1.				*/
				/*--------------------------------------------------------------*/
				day = day + 1;
				next_date = caldat(julday(current_date)+1);
				current_date.day = next_date.day;
				current_date.hour = next_date.hour;
				if (command_line[0].verbose_flag > 0)
					fprintf(stderr,"\n\nYEAR %d MONTH %d DAY %d\n\n",
					current_date.year,current_date.month,current_date.day);

#ifdef JMG_TRACKING
                /* Update simulation day count */
                world[0].track_simtime.sday += 1;
#endif
			} /*end if*/
			/*--------------------------------------------------------------*/
			/*			Check if this is a month end.						*/
			/*--------------------------------------------------------------*/
			if ( next_date.month !=	current_date.month ){                
				/*--------------------------------------------------------------*/
				/*				Do monthly stuff.								*/
				/*--------------------------------------------------------------*/
//#ifdef JMG_TRACKING
//        /* Update simulation month count */
//                world[0].track_simtime.smth += 1;
//#endif
				/*--------------------------------------------------------------*/
				/* if fire spread is called - initiate fire spread routine 	*/
				/*--------------------------------------------------------------*/
                if (command_line[0].firespread_flag == 1)
                {
//#ifdef LIU_BURN_ALL_AT_ONCE
                    if ((!command_line[0].user_defined_fire_event_flag) ||
                         (command_line[0].user_defined_fire_event_flag &&
                          command_line[0].burn_on_flag == 1))
//#endif
					execute_firespread_event(
						world,
						command_line,
						current_date);
				}

				/*--------------------------------------------------------------*/
				/*			Perform any requested monthly output				*/
				/*--------------------------------------------------------------*/
                if (command_line[0].output_flags.monthly == 1){
						execute_monthly_output_event(
						world,
						command_line,
						current_date,
						outfile);
                }
				/*--------------------------------------------------------------*/
				/*				increment month 								*/
				/*--------------------------------------------------------------*/
				month = month + 1;
				current_date.month = next_date.month;

#ifdef JMG_TRACKING
                /* Update simulation month count */
                world[0].track_simtime.smth += 1;
#endif
			} /* end if */
			/*--------------------------------------------------------------*/
			/*			Check if this is a year end.						*/
			/*--------------------------------------------------------------*/
			if ( next_date.year != current_date.year ){
				/*--------------------------------------------------------------*/
				/*				Do yearly stuff.								*/
				/*--------------------------------------------------------------*/
//#ifdef JMG_TRACKING
//        /* Update simulation */
//                world[0].track_simtime.syr += 1;
//#endif
				/*--------------------------------------------------------------*/
				/* if beetle spread is called - initiate beetle spread routine 	*/
				/*--------------------------------------------------------------*/
	     		if (command_line[0].beetlespread_flag == 1 && current_date.year ==world[0].defaults[0].beetle[0].year_attack && world[0].defaults[0].beetle[0].calc_single_attack>0) {
					execute_beetlespread_event(
						world,
						command_line,
						current_date);
				}

				/*--------------------------------------------------------------*/
				/*				increment year  								*/
				/*-------------------------------------------------------------*/
//#ifdef LIU_DISPLY_RUN_INFO
                if (command_line[0].verbose_flag == -7) printf("\nYear %d\n", current_date.year);
//#endif
				year = year + 1;
                current_date.year= next_date.year;

#ifdef JMG_TRACKING
                /*Update simulation */
                world[0].track_simtime.syr += 1;
#endif

//#ifdef LIU_BURN_ALL_AT_ONCE
                if (command_line[0].fire_spin_flag != 0) {
                    if(spins_index < command_line[0].fire_spins) {
//#ifdef LIU_DISPLY_RUN_INFO
                        if (command_line[0].verbose_flag == -7) printf("Spins:%d Spin_year:%d\n",spins_index,spin_year_index);
//#endif
                        if (spin_year_index < (command_line[0].fire_spin_period - 1)) {
                            spin_year_index++;

                            //07272022LML select random climate data year from spin-up period
                            current_date.year = world[0].start_date.year
                                    + rand() % command_line[0].fire_spin_period;
                            next_date.year = current_date.year;
                            //int t1 = get_indays(current_date.year,current_date.month,current_date.day,1901,1);
                            //int t2 = get_indays(world[0].start_date.year,world[0].start_date.month,world[0].start_date.day,1901,1);
                            //day = t1 - t2;
                            //if (day < 0) {
                            //    printf("Date error for spin-up!\n");
                            //    exit(-1);
                            //}


                        } else {
                            spins_index++;
                            spin_year_index = 0;
                            current_date = world[0].start_date;
                            next_date.year = current_date.year;
                            day = 0;
                        }
                    }
                }
//#endif

/*#ifdef JMG_TRACKING*/
                /* Update simulation */
/*                world[0].track_simtime.syr += 1;
#endif */

            }  /*end if*/
        } /*end while*/
		} /*end while*/
		return;
} /*end execute_tec.c*/

//get total number of sections from a char array
int numer_sections_char_array(char* line) {
    char *token;
    // Read a line from the file
    char *cline = malloc(strlen(line) + 1);
    int tnum = 0;
    if (cline != NULL) {
        strcpy(cline, line);
        cline[strcspn(cline, "\n")] = '\0';
        token = strtok(cline, " \t");

        while (token != NULL) {
          //printf("Token: %s\n", token);
          token = strtok(NULL, " \t");
          tnum++;
        }
        free(cline);
    }
    return tnum;
}
