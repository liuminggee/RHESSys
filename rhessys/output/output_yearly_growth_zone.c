/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_growth_zone						*/
/*																*/
/*	output_yearly_growth_zone - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_yearly_growth_zone - outputs current contents of a zone.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_growth_zone(int basinID, int hillID,	*/
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

void	output_yearly_growth_zone(	int basinID, int hillID,
				struct	zone_object	*zone,
				struct	date	current_date,
                FILE *outfile
#ifdef JMG_TRACKING
                ,struct simtime *simtime
#endif
                                    )
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	/*Not Implemented yet*/
} /*end output_yearly_growth_zone*/
