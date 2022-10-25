/*--------------------------------------------------------------*/
/* 																*/
/*					output_yearly_hillslope						*/
/*																*/
/*	output_yearly_hillslope - creates output files objects.		*/
/*																*/
/*	NAME														*/
/*	output_yearly_hillslope - outputs current contents of a hillslope.			*/
/*																*/
/*	SYNOPSIS													*/
/*	void	output_yearly_hillslope(										*/
/*					struct	hillslope_object	*hillslope,				*/
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

void	output_yearly_hillslope(	int basinID,
							 struct	hillslope_object	*hillslope,
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
	int check;


	if (hillslope[0].acc_year.length == 0) hillslope[0].acc_year.length = 1;

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic[] = "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n";
#else
    char out_basic[] = "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n";
#endif

#ifdef JMG_TRACKING
    char out_format[1000] = "%d %d %d ";
    strcat(out_format,out_basic);
#else
    char out_format[1000] = "";
    strcat(out_format, out_basic);
#endif

    check = fprintf(outfile, out_format,

#ifdef JMG_TRACKING
        simtime->sday,
        simtime->smth,
        simtime->syr,
#endif

		current_date.year-1,
		basinID,
		hillslope[0].ID,
#ifndef JMG_MORE_YEARLY_OUTPUT
        hillslope[0].acc_year.streamflow * 1000.0,
        hillslope[0].acc_year.stream_NO3 * 1000.0,
        hillslope[0].acc_year.denitrif * 1000.0,
        hillslope[0].acc_year.DOC_loss * 1000.0,
        hillslope[0].acc_year.DON_loss * 1000.0,
        hillslope[0].acc_year.et * 1000.0,
        hillslope[0].acc_year.psn * 1000.0,
        hillslope[0].acc_year.lai/ hillslope[0].acc_year.length,
        hillslope[0].acc_year.num_threshold,
        hillslope[0].acc_year.lai/ hillslope[0].acc_year.length,
        hillslope[0].acc_year.nitrif * 1000.0,
        hillslope[0].acc_year.mineralized * 1000.0,
        hillslope[0].acc_year.uptake * 1000.0,
#else
        hillslope[0].acc_year.pch_pcp * 1000.0,
        hillslope[0].acc_year.pch_et * 1000.0,
        hillslope[0].acc_year.pch_streamflow * 1000.0,
        hillslope[0].acc_year.pch_return_flow * 1000.0,
        hillslope[0].acc_year.pch_base_flow * 1000.0,
        hillslope[0].acc_year.hill_base_flow * 1000.0,
        hillslope[0].acc_year.pch_gw_drainage * 1000.0,
        (hillslope[0].acc_year.pch_rz_storage / hillslope[0].acc_year.length) * 1000.0, // average daily rz_storage
        (hillslope[0].acc_year.pch_unsat_storage / hillslope[0].acc_year.length) * 1000.0, // average daily unsat_storage
        (hillslope[0].acc_year.hill_gw_storage / hillslope[0].acc_year.length) * 1000.0, // average daily gw_storage
#endif
		hillslope[0].area
		);
	if (check <= 0) {
		fprintf(stdout,
			"\nWARNING: output error has occured in output_yearly_hillslope");
	}
	/*--------------------------------------------------------------*/
	/*	reset accumulator variables				*/
	/*--------------------------------------------------------------*/
	hillslope[0].acc_year.streamflow = 0.0;
	hillslope[0].acc_year.stream_NO3 = 0.0;
	hillslope[0].acc_year.et = 0.0;
	hillslope[0].acc_year.psn = 0.0;
	hillslope[0].acc_year.lai = 0.0;
	hillslope[0].acc_year.length = 0;
	hillslope[0].acc_year.DOC_loss = 0.0;
	hillslope[0].acc_year.DON_loss = 0.0;
	hillslope[0].acc_year.denitrif= 0.0;
	hillslope[0].acc_year.nitrif= 0.0;
	hillslope[0].acc_year.num_threshold = 0;
	hillslope[0].acc_year.lai = 0.0;
	hillslope[0].acc_year.mineralized = 0.0;
	hillslope[0].acc_year.uptake = 0.0;

#ifdef JMG_MORE_YEARLY_OUTPUT
    hillslope[0].acc_year.pch_pcp = 0.0; // JMG09122022
    hillslope[0].acc_year.pch_et = 0.0; // JG09122022
    hillslope[0].acc_year.pch_streamflow = 0.0; // JMG09122022
    hillslope[0].acc_year.pch_base_flow = 0.0; // JMG09122022
    hillslope[0].acc_year.pch_return_flow = 0.0; // JMG09122022
    hillslope[0].acc_year.pch_gw_drainage = 0.0; // JMG09122022
    hillslope[0].acc_year.pch_rz_storage = 0.0; // JMG09122022
    hillslope[0].acc_year.pch_unsat_storage = 0.0; // JMG09122022
    hillslope[0].acc_year.hill_gw_storage = 0.0; // JMG09122022
    hillslope[0].acc_year.hill_base_flow = 0.0; // JMG09122022
#endif

	return;
} /*end output_yearly_hillslope*/
