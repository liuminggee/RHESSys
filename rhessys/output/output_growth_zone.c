/*--------------------------------------------------------------*/
/* 																*/
/*					output_growth_zone						*/
/*																*/
/*	output_growth_zone - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_growth_zone - outputs current contents of a zone.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_growth_zone(										*/
/*					struct	zone_object	*zone,				*/
/*					struct	date	date,  						*/
/*					FILE 	*outfile)							*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	outputs spatial structure according to commandline			*/
/*	specifications to specific files							*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	We only permit one fileset per spatial modelling level.     */
/*	Each fileset has one file for each timestep.  				*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include "rhessys.h"

void	output_growth_zone(	int basinID, int hillID,
						   struct	zone_object	*zone,
						   struct	date	current_date,
                           FILE *outfile
#ifdef JMG_TRACKING
                           ,struct simtime *simtime
#endif
                            )
{
	/*------------------------------------------------------*/
	/*	Local Function Declarations.						*/
	/*------------------------------------------------------*/
	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/

    char out_basic[] = "%4d %4d %4d %3d %3d %3d %8.5f %8.5f %8.3f %8.3f %8.5f %f %f %f %f \n ";

#ifdef JMG_TRACKING
    char out_format[] = "%d %d %d ";
    strcat(out_format,out_basic);
#else
    char out_format[] = "";
    strcat(out_format, out_basic);
#endif
	
    fprintf(outfile, out_format,

#ifdef JMG_TRACKING
            simtime->sday,
            simtime->smth,
            simtime->syr,
#endif

        current_date.day,
		current_date.month,
		current_date.year,
		basinID,
		hillID,
		zone[0].ID,
		zone[0].rain * 1000.0,
		zone[0].snow * 1000.0,
		zone[0].metv.tday,
		zone[0].metv.tavg,
		zone[0].metv.vpd,
		zone[0].Kdown_direct,
		zone[0].Kdown_diffuse,
		zone[0].PAR_direct,
		zone[0].PAR_diffuse);
	return;
} /*end output_growth_zone*/
