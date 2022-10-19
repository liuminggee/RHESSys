/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_canopy_stratum				*/
/*																*/
/*	output_yearly_canopy_stratum - creates output files objects.*/
/*																*/
/*	NAME														*/
/*	output_yearly_canopy_stratum - outputs current contents of a canopy_stratum*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_canopy_stratum(int bainsID, int hillID,*/
/*					int zoneID, int patchID,
/*					struct	canopy_stratum_object	*canopy_stratum,				*/
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

void	output_yearly_canopy_stratum( int basinID, int hillID,
									 int zoneID,
									 int patchID,
									 int reset_flag,
									 struct	canopy_strata_object	*stratum,
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
	
	if (stratum[0].acc_year.length == 0)
		stratum[0].acc_year.length = 1;
	/*--------------------------------------------------------------*/
	/*	output variables					*/
	/*--------------------------------------------------------------*/

    char out_basic[] = "%4d %d %d %d %d %d %lf %lf %lf\n";

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

		current_date.year,
		basinID,
		hillID,
		zoneID,
		patchID,
		stratum[0].ID,
		(stratum[0].acc_year.psn) / stratum[0].acc_year.length,
		(stratum[0].acc_year.lwp) / stratum[0].acc_year.length,
		stratum[0].rootzone.depth*1000.0);
	/*--------------------------------------------------------------*/
	/*	reset accumulator variables				*/
	/*--------------------------------------------------------------*/
	if (reset_flag == 1) {
	stratum[0].acc_year.psn = 0.0;
	stratum[0].acc_year.lwp = 0.0;
	stratum[0].acc_year.length = 0;
	}
	return;
	
} /*end output_yearly_canopy_stratum*/
