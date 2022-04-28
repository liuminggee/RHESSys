/*--------------------------------------------------------------*/
/* 																*/
/*					assign_base_station.c 						*/
/*																*/
/*																*/
/*	NAME														*/
/*																*/
/*	SYNOPSIS													*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*	Original code, January 15, 1996.							*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include "rhessys.h"
bool is_close_to_station(const double x, const double y, const base_station_object *station,
                         const base_station_ncheader_object *ncheader);
#ifdef LIU_FIND_GEO_CLOSEST_STATION
int get_closest_station_index(const double x, const double y, const base_station_object **stations,
                         const base_station_ncheader_object *ncheader, int num_base_stations);
int get_station_index_from_ID(const double x, const double y, const base_station_object **stations,
                              int num_base_stations, int search_ID);
#endif
struct base_station_object
		*assign_base_station_xy(
					 float		x,
					 float		y,
					 int		num_base_stations,
					 int		*notfound,
                     struct	base_station_object	**base_stations,
                     const struct base_station_ncheader_object *ncheader       
                     //#ifdef LIU_NETCDF_READER
                     //double dist_tol
                     //#endif
                     #ifdef FIND_STATION_BASED_ON_ID
                     ,const int basestation_id
                     #endif
        )
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/
	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	int	i;
	struct	base_station_object *base_station;
	/*--------------------------------------------------------------*/
	/*	Loop through all of the basestations available.			*/
	/*	and find the record which holds the matching base station	*/
	/*--------------------------------------------------------------*/
	i = 0;
	if (num_base_stations < 1) {
		*notfound = 1;
		return 0;
	}
	else {
#ifndef LIU_FIND_GEO_CLOSEST_STATION
        #ifdef FIND_STATION_BASED_ON_ID
        while ( (*(base_stations[i])).ID != basestation_id ) {
                
        #else
        while (!is_close_to_station(x,y,base_stations[i],ncheader)) {
        /*!is_approximately(x,(*(base_stations[i])).x,dist_tol) || !is_approximately(y,(*(base_stations[i])).y,dist_tol) */
        #endif
            //printf("\n      Assign: Starting while loop: i:%d\tzone_baseid:%d\tbstation_id:%d\tbstation_x:%lf\tbstation_y:%lf\n",i,basestation_id,(*(base_stations[i])).ID,(*(base_stations[i])).x,(*(base_stations[i])).y);
			i++;
			/*--------------------------------------------------------------*/
			/*	Report an error if no match was found.  Otherwise assign	*/
			/*	the base_station_pointer to point to this base_station.		*/
			/*--------------------------------------------------------------*/

            //printf("i:%d x_diff:%lf y_diff:%lf\n",i,x-base_stations[i]->lon,y-base_stations[i]->lat);

			if ( i >= num_base_stations ){
				//fprintf(stderr,"\n      Assign: NOT FOUND. Adding new base station for %lf %lf",x,y);
				*notfound = 1;
				return 0;
			}
		}  /* end-while */
        base_station = base_stations[i];
        return(base_station);
#else
        int find_index =
#ifdef FIND_STATION_BASED_ON_ID
            get_station_index_from_ID(x, y, base_stations, num_base_stations, basestation_id);
#else
            get_closest_station_index(x, y, base_stations, ncheader, num_base_stations);
#endif
        if (find_index < 0) {
            fprintf(stderr,"\n      Assign: NOT FOUND. Adding new base station for %lf %lf",x,y);
            return 0;
        } else {
            base_station = base_stations[find_index];
            return base_station;
        }
#endif
	}
} /*end assign_base_station*/
//_____________________________________________________________________
bool is_close_to_station(const double x, const double y, const base_station_object *station,
                         const base_station_ncheader_object *ncheader)
{ //Check if the site is close to station or not
    bool is_close = 0;
    if ((is_approximately(x,station->proj_x,ncheader->resolution_meter / 2.0) && is_approximately(y,station->proj_y,ncheader->resolution_meter / 2.0)) ||
        (is_approximately(x,station->lon,   ncheader->resolution_dd    / 2.0) && is_approximately(y,station->lat,   ncheader->resolution_dd    / 2.0)))
        is_close = 1;
    return is_close;
}
#ifdef LIU_FIND_GEO_CLOSEST_STATION
//_____________________________________________________________________
int get_closest_station_index(const double x, const double y, const base_station_object **stations,
                         const base_station_ncheader_object *ncheader, int num_base_stations)
{ //Check if the site is close to station or not

    //100821LML WARNING: assume zone & station are using geographic projection!!!

    int closest_index = -1;
    if (num_base_stations > 0) {
        double mindist = powf(x - stations[0]->lon,2) + powf(y - stations[0]->lat,2);
        for (int i = 0; i < num_base_stations; i++) {
            double dist = powf(x - stations[i]->lon,2) + powf(y - stations[i]->lat,2);
            if (dist <= mindist) {
                mindist = dist;
                closest_index = i;
            }
        }
        if (sqrt(mindist) >= 2.0 * ncheader->resolution_dd) {
            closest_index = -1;
        }
    }
    return closest_index;
}
//_____________________________________________________________________
int get_station_index_from_ID(const double x, const double y, const base_station_object **stations,
                              int num_base_stations, int search_ID)
{ //Check if the site is close to station or not

    //100821LML WARNING: assume zone & station are using geographic projection!!!

    int index = -1;
    if (num_base_stations > 0) {
        for (int i = 0; i < num_base_stations; i++) {
            if (search_ID == stations[i]->ID) {
                index = i;
            }
        }
    }
    return index;
}
#endif
