/*--------------------------------------------------------------*/
/* 																*/
/*					execute_firespread_event								*/
/*																*/
/*	execute_firespread_event.c - creates a patch object					*/
/*																*/
/*	NAME														*/
/*	execute_firespread_event.c - creates a patch object					*/
/*																*/
/*	SYNOPSIS													*/
/*																*/
/* 																*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*																*/
/*																*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*--------------------------------------------------------------*/
#include <string.h>
#include <stdio.h>
#include "rhessys.h"
#include "util/rhessys_fire.h"

// test comment
void execute_firespread_event(
									 struct	world_object *world,
									 struct	command_line_object	*command_line,
									 struct date	current_date)
{
	/*--------------------------------------------------------------*/
	/*	Local function definition.									*/
	/*--------------------------------------------------------------*/


	void *alloc(size_t, char *, char *);


	void	compute_fire_effects(
		struct patch_object *,
		double,
		struct	command_line_object	*);


	/*--------------------------------------------------------------*/
	/*	Local variable definition.									*/
	/*--------------------------------------------------------------*/
	struct fire_object **fire_grid;
	struct patch_fire_object **patch_fire_grid;
	struct patch_object *patch;
//	struct node_fire_wui_dist *tmp_node;
	int i,j,p,c,layer;
	double pspread;
	double mean_fuel_veg=0,mean_fuel_litter=0,mean_soil_moist=0,mean_fuel_moist=0,mean_relative_humidity=0,
		mean_wind_direction=0,mean_wind=0,mean_z=0,mean_temp=0,mean_et=0,mean_pet=0,mean_understory_et=0,mean_understory_pet=0, mean_trans=0;
	double denom_for_mean=0;

	patch_fire_grid=world[0].patch_fire_grid;
	fire_grid = world[0].fire_grid;

// add code here to define understory et and pet to calculate understory deficit. The definition of understory here
// will be the same as that for fire effects below

// default value to determine if you are in the understory, compare layer height to stratum-level height threshold (canopy_strata_upper[0].defaults[0][0].understory_height_thresh,
// compared to patch[0].canopy_strata[(patch[0].layers[layer+1].strata[c])].epv.height

	/*--------------------------------------------------------------*/
	/* update fire grid variables			*/
	/* first reset the values				*/
	/*--------------------------------------------------------------*/
    if(world[0].defaults[0].fire[0].fire_verbose == 3) { //NREN 20190912
	printf("In WMFire\n");}
    #pragma omp parallel for reduction(+ : denom_for_mean,mean_fuel_veg,mean_fuel_litter \
        ,mean_fuel_moist,mean_soil_moist,mean_relative_humidity,mean_wind_direction \
        ,mean_wind,mean_temp,mean_et,mean_pet,mean_understory_et,mean_understory_pet,mean_trans)
    for  (int i=0; i< world[0].num_fire_grid_row; i++) {
      for (int j=0; j < world[0].num_fire_grid_col; j++) {
        struct  fire_object *pfire = &fire_grid[i][j];
        pfire->fire_size=0; // reset grid to no fire
        if(patch_fire_grid[i][j].occupied_area==0)
		{
			  if(world[0].defaults[0].fire[0].fire_in_buffer==0)
			  {
                    pfire->fuel_veg = 0.0; // this should work to initialize the grid, so if none of the patches overlap a grid point the fuel is zero and fire doesn't spread
                    pfire->fuel_litter = 0.0;
                    pfire->fuel_moist = 100.0;
                    pfire->soil_moist = 100.0;
                    pfire->relative_humidity = 100.0;
                    pfire->wind_direction = 0;
                    pfire->wind = 0.0;
                    pfire->z=0.0; // mk: add so we can calculate the current elevation as the weighted mean elevation of the patches
                    pfire->temp=0.0;
                    pfire->et=1.0; // mk: should be 1 so that the deficit in the buffer is zero
                    pfire->pet=1.0;
                    pfire->understory_et=0.0;
                    pfire->understory_pet=1.0;
                    pfire->trans=0.0; //NR add transpiration
                    pfire->ign_available=0;
			  }
			  else // if denom_for_mean==0, then this initializes the buffer, otherwise the mean is filled in below
			  {
                    pfire->fuel_veg = 10.0; // this should work to initialize the grid, so if none of the patches overlap a grid point the fuel is zero and fire doesn't spread
                    pfire->fuel_litter = 10.0;
                    pfire->fuel_moist = 0;
                    pfire->soil_moist = 0;
                    pfire->relative_humidity = 0;
                    pfire->wind_direction = world[0].basins[0][0].hillslopes[0][0].zones[0][0].wind_direction;//patches[0].zone[0].wind_direction;// or pull the wind direction from the default wind
                    pfire->wind = 50.0;
                    pfire->z=patch_fire_grid[i][j].elev; // mk: add so we can calculate the current elevation as the weighted mean elevation of the patches
                    pfire->temp=0.0;
                    pfire->et=0.0;
                    pfire->pet=0.0;
                    pfire->understory_et=0.0;
                    pfire->understory_pet=0.0;
                    pfire->trans=0.0; //NR add transpiration
                    pfire->ign_available=0;
			  }
		}
		else
		{
            pfire->fuel_veg = 0.0; // this should work to initialize the grid, so if none of the patches overlap a grid point the fuel is zero and fire doesn't spread
            pfire->fuel_litter = 0.0;
            pfire->fuel_moist = 0.0;
            pfire->soil_moist = 0.0;
            pfire->relative_humidity = 0.0;
            pfire->wind_direction = 0.0;
            pfire->wind = 0.0;
            pfire->z=0.0; // mk: add so we can calculate the current elevation as the weighted mean elevation of the patches
            pfire->temp=0.0;
            pfire->et=0.0;
            pfire->pet=0.0;
            pfire->understory_et=0.0;
            pfire->understory_pet=0.0;
            pfire->trans=0.0; //NR add transpiration
            pfire->ign_available=1;	/* then make this available for ignition */
            pfire->burned = 0;
            //printf("col:%d\trow:%d\tign_available=1\n",j,i);
		}
        //if (patch_fire_grid[i][j].num_patches > 0)
        //  printf("checking num patches. row %d col %d numPatches %d\n",i,j,patch_fire_grid[i][j].num_patches);
        int num_patches = patch_fire_grid[i][j].num_patches;
        for (int p=0; p < num_patches; ++p) { // should just be 1 now...
//printf("Patch p: %d\n",p);
            double f_patch = patch_fire_grid[i][j].prop_patch_in_grid[p];
            struct patch_object *patch = patch_fire_grid[i][j].patches[p]; //So this is patch family now? points to patch family
//printf("Patch p1 %lf\n", patch[0].litter_cs.litr1c);

//printf("Patch p2: %d\n",p);
            if (world[0].defaults[0].fire[0].calc_above_ground_litter == 1) {
                pfire->fuel_litter += (patch[0].litter_cs.litr1c +	patch[0].litter_cs.litr2c +	// This sums the litter pools
                    patch[0].litter_cs.litr3c +	patch[0].litter_cs.litr4c) * f_patch * patch[0].prop_litrc_above_ground;
                }
            else
                pfire->fuel_litter += (patch[0].litter_cs.litr1c +	patch[0].litter_cs.litr2c +	// This sums the litter pools
                    patch[0].litter_cs.litr3c +	patch[0].litter_cs.litr4c) * f_patch;


			if( patch[0].litter.rain_capacity!=0)	// then update the fuel moisture, otherwise don't change it
                pfire->fuel_moist += (patch[0].litter.rain_stored / patch[0].litter.rain_capacity) * f_patch;
/*			fire_grid[i][j].fuel_moist += (patch[0].litter.rain_stored / patch[0].litter.rain_capacity) *
                        f_patch;
*/
//printf("Patch p: %d\n",p);

	// this is the canopy fuels

            for (int layer=0 ; layer<patch[0].num_layers; layer++ ){
                for (int c=0 ; c<patch[0].layers[layer].count; c++ ){
		//printf("Layers: %d\n",layer);

                    pfire->fuel_veg += (patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction
                        * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cs.leafc) * f_patch ;
                    pfire->fuel_litter +=(patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cover_fraction // add the dead leaf to litter pool
                        * patch[0].canopy_strata[(patch[0].layers[layer].strata[c])][0].cs.dead_leafc) * f_patch ;

                }
            }
//			printf("pixel veg and prop patch in grid: %lf\t%lf\n",pfire->fuel_veg,f_patch);


            pfire->soil_moist += patch[0].rootzone.S * f_patch;	//soil moisture, divided by proportion of the patch in that grid cell;

            pfire->wind += patch[0].zone[0].wind * f_patch;
            pfire->wind_direction += patch[0].zone[0].wind_direction * f_patch;
            pfire->relative_humidity += patch[0].zone[0].relative_humidity * f_patch;
            pfire->z += patch[0].z*f_patch; // elevation
            pfire->temp += patch[0].zone[0].metv.tavg*f_patch;// temperature? mk

            patch[0].fire.et = avgvalue_FIFO_Queue(&patch[0].fire.Q_et) * 1000;                         //m->mm
            patch[0].fire.pet = avgvalue_FIFO_Queue(&patch[0].fire.Q_pet) * 1000;                       //m->mm
            patch[0].fire.trans = avgvalue_FIFO_Queue(&patch[0].fire.Q_trans) * 1000;                   //m->mm
            patch[0].fire.understory_et = avgvalue_FIFO_Queue(&patch[0].fire.Q_understory_et) * 1000;   //m->mm
            patch[0].fire.understory_pet = avgvalue_FIFO_Queue(&patch[0].fire.Q_understory_pet) * 1000; //m->mm

            pfire->et += patch[0].fire.et * f_patch;
            pfire->pet += patch[0].fire.pet * f_patch;
            pfire->trans += patch[0].fire.trans * f_patch; //where is fire.trans
			if(patch[0].fire.understory_et==0 && patch[0].fire.understory_pet==0) // means no understory present, then use overall et and pet for deficit calculation for this patch
			{
                pfire->understory_et = pfire->et;
                pfire->understory_pet = pfire->pet;

			}
			else
			{
                pfire->understory_et += patch[0].fire.understory_et * f_patch;
                pfire->understory_pet += patch[0].fire.understory_pet * f_patch;

			}
	//printf("patch pet, patch et: %lf\t%lf\n",patch[0].fire.pet,patch[0].fire.et);

		}
        if(patch_fire_grid[i][j].occupied_area>0 && world[0].defaults[0].fire[0].fire_in_buffer==1) // if allowing fire into the buffer (on raster grid outside of watershed boundaries), then fill with mean field values within watershed boundary
		{ // this loop fills sums to calculate the mean value across watershed
			denom_for_mean+=1;
            mean_fuel_veg+=pfire->fuel_veg; // this should work to initialize the grid, so if none of the patches overlap a grid point the fuel is zero and fire doesn't spread
            mean_fuel_litter+=pfire->fuel_litter;
            mean_fuel_moist+=pfire->fuel_moist;
            mean_soil_moist+=pfire->soil_moist;
            mean_relative_humidity+=pfire->relative_humidity;
            mean_wind_direction+=pfire->wind_direction;
            mean_wind+=pfire->wind;
            mean_temp+=pfire->temp;
            mean_et+=pfire->et;
            mean_pet+=pfire->pet;
            mean_understory_et+=pfire->understory_et;
            mean_understory_pet+=pfire->understory_pet;
            mean_trans += pfire->trans; //NEW NR
        //	printf("et: %f  pet: %f  ",pfire->et,pfire->pet);
		}

        //07122023LML commented out
        //pfire->et *= 1000; // convert to mm
        //pfire->pet *= 1000; // convert to mm
        //pfire->understory_et *= 1000; // convert to mm
        //pfire->understory_pet *= 1000; // convert to mm
        //pfire->trans *= 1000;//NEW NR

        if (num_patches >= 1 && pfire->understory_et >= 1e-12 && pfire->understory_pet > 1e-12 && pfire->pet > 1e-12) {
          //printf("year:%d month:%d day:%d",current_date.year,current_date.month,current_date.day);
          ////printf(" num_patches:%d",num_patches);
          //printf(" col:%d,row:%d ",j,i);
          //printf(" understory et(%.1fmm)/pet(%.1fmm):%.2f",pfire->understory_et,pfire->understory_pet,(pfire->understory_et/(pfire->understory_pet)));
          //printf(" overcanopy et(%.1fmm)/pet(%.1fmm):%.2f\n",pfire->et,pfire->pet,(pfire->et/(pfire->pet)));

          //hard code for testing...
          //FILE *pFile2;
          //char *filename;
          //sprintf(filename,"/home/liuming/temp/check_pmoisture_fire_%d.txt",current_date.year);
          //pFile2=fopen(filename, "a");
          #ifdef LIU_CHECK_FIRE_OCCURENCE
          #pragma omp critical
          {
            fprintf(global_debug,"year:%d month:%d day:%d",current_date.year,current_date.month,current_date.day);
            //printf(" num_patches:%d",num_patches);
            fprintf(global_debug," col:%d,row:%d ",j,i);
            fprintf(global_debug," understory et(%.1fmm)/pet(%.1fmm):%.2f",pfire->understory_et,pfire->understory_pet,(pfire->understory_et/(pfire->understory_pet)));
            fprintf(global_debug," overcanopy et(%.1fmm)/pet(%.1fmm):%.2f\n",pfire->et,pfire->pet,(pfire->et/(pfire->pet)));
          }
          #endif
          //fclose(pFile2);
        }

      } //j
    } //i
//	printf("denom: %lf\t",denom_for_mean);
	if(denom_for_mean>0&&world[0].defaults[0].fire[0].fire_in_buffer==1) // so here we calculate the mean value
	{
//		printf("in denom if\n");
        mean_fuel_veg /= denom_for_mean;
        mean_fuel_litter /= denom_for_mean;
        mean_fuel_moist /= denom_for_mean;
        mean_soil_moist /= denom_for_mean;
        mean_relative_humidity /= denom_for_mean;
        mean_wind_direction /= denom_for_mean;
        mean_wind /= denom_for_mean;
        mean_temp /= denom_for_mean;
        mean_et /= denom_for_mean;
        mean_pet /= denom_for_mean;
        mean_understory_et /= denom_for_mean;
        mean_understory_pet /= denom_for_mean;
        mean_trans /= denom_for_mean;

	//	printf("mean et: %f  mean pet: %f  ",mean_et,mean_pet);

	//	printf("mean pet, mean et: %lf\t%lf\n",mean_pet,mean_et);
	//	printf("mean wind: %lf, mean direction %lf \n",mean_wind,mean_wind_direction);
        #pragma omp parallel for
        for  (int i=0; i< world[0].num_fire_grid_row; i++) {
          for (int j=0; j < world[0].num_fire_grid_col; j++) {
              if(patch_fire_grid[i][j].occupied_area==0) // and here we fill in the buffer
			  {
                struct  fire_object *pfire = &fire_grid[i][j];
                pfire->fuel_veg = mean_fuel_veg; // this should work to initialize the grid, so if none of the patches overlap a grid point the fuel is zero and fire doesn't spread
                pfire->fuel_litter = mean_fuel_litter;
                pfire->fuel_moist = mean_fuel_moist;
                pfire->soil_moist = mean_soil_moist;
                pfire->relative_humidity = mean_relative_humidity;
                pfire->wind_direction = mean_wind_direction;
                pfire->wind = mean_wind;
                pfire->temp=mean_temp;
                pfire->z=patch_fire_grid[i][j].elev;
                pfire->et=mean_et; //*1000; // convert to mm
                pfire->pet=mean_pet; //*1000; // convert to mm
                pfire->understory_et=mean_understory_et; //*1000; // convert to mm
                pfire->understory_pet=mean_understory_pet; //*1000; // convert to mm
                pfire->trans=mean_trans; //*1000;//NEW NR
	//	printf("in denom if take 2 update values\n");
			  }
             }//j
        }//i
	}

	/*--------------------------------------------------------------*/
	/* Call WMFire	 						*/
	/*--------------------------------------------------------------*/
	if(world[0].defaults[0].fire[0].fire_verbose == 3) { //NREN 20190912
	printf("calling WMFire: month %ld year %ld  cell res %lf  nrow %d ncol % d\n",current_date.month,current_date.year,command_line[0].fire_grid_res,world[0].num_fire_grid_row,world[0].num_fire_grid_col);}
// needs to return fire size, not just grid--create structure that includes fire size, or a 12-member array of fire sizes, and/or a tally of fires > 1000 acres

//#ifndef LIU_BURN_ALL_AT_ONCE
    if (!command_line->user_defined_fire_event_flag) {
      world[0].fire_grid=WMFire(command_line[0].output_prefix,command_line[0].fire_grid_res,world[0].num_fire_grid_row,world[0].num_fire_grid_col,current_date.year,current_date.month,world[0].fire_grid,*(world[0].defaults[0].fire));
    }
//#endif
    if(world[0].defaults[0].fire[0].fire_verbose == 3) {
 	printf("Finished calling WMFire\n"); }
	/*--------------------------------------------------------------*/
	/* update biomass after fire					*/
	/*--------------------------------------------------------------*/

	// if(world[0].fire_grid[0][0].fire_size>0) // only do this if there was a fire

//#ifdef LIU_BURN_ALL_AT_ONCE
//    double total_pspread = 0;                                                   //No real meaning, just for counting
//#endif
    #pragma omp parallel for
    for  (int i=0; i< world[0].num_fire_grid_row; i++) {
        for (int j=0; j < world[0].num_fire_grid_col; j++) {
            //if(patch_fire_grid[i][j].num_patches > 0)
            //    printf("i:%d j:%d num_patches:%d",i,j,patch_fire_grid[i][j].num_patches);
            struct  fire_object *pfire = &fire_grid[i][j];
            struct patch_fire_object *ppatch_fire = &patch_fire_grid[i][j];
            for (int p=0; p < patch_fire_grid[i][j].num_patches; ++p) {
                struct patch_object *patch = ppatch_fire->patches[p];
//#ifndef LIU_BURN_ALL_AT_ONCE
                if (!command_line->user_defined_fire_event_flag) {
                  patch[0].burn = pfire->burn * ppatch_fire->prop_grid_in_patch[p];
                  pspread = patch[0].burn;
                //if(patch_fire_grid[i][j].num_patches > 0)
                //    printf("\tp:%d prop_grid_in_patch:%f\n",p,patch_fire_grid[i][j].prop_grid_in_patch[p]);
                } else {
//#else
                //01092024 LML be careful that one same patch may have multiple
                //         fire pixels!
                //total_pspread += pspread;
                //05062022LML pspread = 1.0;                                                  //Handle the patch once for all pixels since they share same C&N pools
                  pspread = command_line[0].fire_pspread * ppatch_fire->prop_grid_in_patch[p]; //05182022LML
                  if (pspread < 0) {
                    if (command_line[0].fire_understory_mortality_rate > 0)     //use predefined understory mortality rate
                        pspread = command_line[0].fire_understory_mortality_rate
                                * ppatch_fire->prop_grid_in_patch[p];
                    else pspread = ppatch_fire->prop_grid_in_patch[p];          //assume firespread rate is 1
                  }
                }
//#endif
// so I think here we could flag whether to turn salient fire on in wui; convert fire size in pixels to ha, assuming the cell_res is in m
				/* (if pspread>0&world[0].fire_grid[0][0].fire_size*command_line[0].fire_grid_res*command_line[0].fire_grid_res*0.0001>=400) // also need a flag with the fire size to trigger event, because fire > 400 ha
				{
					// linked list loop
					for(w=0;w<3;w++) # for each level of salience, 1 = <= 3 km, 2 = <=5 km; 3=<=10 km
					{
                        tmp_node=fire_grid[i][j].wuiList[w] // where wuiList[w] is the linked list of patches within w index of this pixel
						while(tmp_node!=NULL)
						{
							if(tmp_node->patch.wuiFire==0||(i+1)<tmp_node.patch.wuiFire) // then this pixel is closer to the wui and should activate a more salient fire
								tmp_node->patch.wuiFire=i+1; // then flag this wuiPatch with a salient fire event
							tmp_node=tmp_node->next;
						}
					}
				}

				*/

                if(world[0].defaults[0].fire[0].calc_fire_effects==1)
				{
                    //printf("calling WMFire: pspread is %lf \n, burn is %lf \n", pspread, fire_grid[i][j].burn);
                    //printf("row:%d col:%d patch_ID:%d\n",i,j,patch->ID);
					compute_fire_effects(
						patch,
						pspread,
						command_line);
                }
            } //p
        } //j
    } //i
//#ifdef LIU_BURN_ALL_AT_ONCE
//    if (total_pspread > 0.1)
//        printf("Has burn event! total_pspread:%f\n",total_pspread);
//#endif
	return;
} /*end execute_firespread_event.c*/


// ----------printing code samples -----------
//	printf("in execute_firespread_event_1\n");
//	printf("*********************\n");
//	printf("Understory Height = %f\n", height_under);






/*

// For 	Local variable definition.
struct fire_veg_loss_struct	fire_veg_effects;

// RHESSys.h
/*----------------------------------------------------------*/
/*      Define post fire vegetation loss structure.    	    */
/*----------------------------------------------------------*/
/*struct fire_veg_loss_struct
{
	double pspread
	double layer_upper_height
	double layer_lower_height
	double layer_lower_c
	double understory_litter_c
};*/



