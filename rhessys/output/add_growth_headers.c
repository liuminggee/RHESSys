/*--------------------------------------------------------------*/
/* 																*/
/*					add_growth_headers					*/
/*																*/
/*	add_growth_headers - 												    	*/
/*	NAME														*/
/*	add_growth_headers    													*/
/*	SYNOPSIS													*/
/*	void add_growth_headers(struct world output_file_object *,				*/
/*			struct command_line_object *)					*/
/*																*/
/*	OPTIONS														*/
/*																*/
/*	DESCRIPTION													*/
/*																*/
/*	Adds headers for yearly, monthly, daily and	*/
/*	hourly basin output 					*/
/*																*/
/*	PROGRAMMER NOTES											*/
/*																*/
/*																*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "rhessys.h"


void add_growth_headers(struct world_output_file_object *world_output_files,
            struct command_line_object *command_line)
{
    /*--------------------------------------------------------------*/
    /*	Local function definition.									*/
    /*--------------------------------------------------------------*/
    /*--------------------------------------------------------------*/
    /*	Local variable definition.									*/
    /*--------------------------------------------------------------*/
    FILE *outfile;
    int check;
    /*--------------------------------------------------------------*/
    /*	Basin file headers					*/
    /*--------------------------------------------------------------*/

    if (command_line[0].b != NULL) {

    /*--------------------------------------------------------------*/
    /*	Hourly							*/
    /*--------------------------------------------------------------*/
    outfile = world_output_files[0].basin[0].hourly;
    fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
        "hour",
        "day",
        "month",
        "year",
        "basinID",
        //"lai",
        //"gpsn",
        //"plant_resp",
        //"soil_resp",
        //"nitrate",
        //"sminn",
        //"surfaceN",
        //"plantc",
        //"plantn",
        //"npool",
        //"litrc",
        //"litrn",
        //"soilc",
        //"soiln",
        "gwNO3",
        "gwNH4",
        "gwDON",
        "gwDOC",
        "streamflow_NO3",
        "streamflow_NH4",
        "streamflow_DON",
        "streamflow_DOC",
        "gwNO3out",
        "gwNH4out",
        "gwDONout",
        "gwDOCout",
        //"denitrif",
        //"nitrif",
        //"DOC",
        //"DON",
        //"root_depth",
        //"nfix",
        //"nuptake",
        //"grazingC",
        "StreamNO3_from_surface",
        "StreamNO3_from_sub");
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

    char out_basic_basin_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_basin_daily[1000] = "%s %s %s ";
    strcat(out_format_basin_daily,out_basic_basin_daily);
#else
    char out_format_basin_daily[1000] = "";
    strcat(out_format_basin_daily, out_basic_basin_daily);
#endif

    outfile = world_output_files[0].basin[0].daily;

    fprintf(outfile, out_format_basin_daily,
            //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "day",
        "month",
        "year",
        "basinID",
        "lai",
        "lai_b",
        "pai",
        "pai_b",
        "gpsn",
        "plant_resp",
        "leaf_resp",
        "soil_resp",
        "nitrate",
        "sminn",
        "surfaceN",
        "plantc",
        "plantn",
        "cpool",
        "npool",
        "litrc",
        "litrn",
        "soilc",
        "soiln",
        "soiln_noslow",
        "gwNO3",
        "gwNH4",
        "gwDON",
        "gwDOC",
        "streamflow_NO3",
        "streamflow_NH4",
        "streamflow_DON",
        "streamflow_DOC",
        "gwNO3out",
        "gwNH4out",
        "gwDONout",
        "gwDOCout",
        "denitrif",
        "nitrif",
        "DOC",
        "DON",
        "root_depth",
        "nfix",
        "nuptake",
        "grazingC",
        "StreamNO3_from_surface",
        "StreamNO3_from_sub", "N_dep","fertilizer_store",
        "understory_leafc",
        "understory_stemc",
        "understory_biomassc",
        "understory_height",
        "overstory_leafc",
        "overstory_stemc",
        "overstory_biomassc",
        "overstory_height",
        //fire
        "burn",
        "litterc_burned",
        "cwdc_to_atoms",
        "overstory_biomassc_consumed",
        "overstory_leafc_consumed",
        "overstory_stemc_consumed",
        "overstory_biomassc_mortality",
        "overstory_leafc_mortality",
        "overstory_stemc_mortality",
        "understory_biomassc_consumed",
        "understory_leafc_consumed",
        "understory_stemc_consumed",

        "total_snagc",
        "total_snagn",
        "total_redneedlec",
        "total_redneedlen",
        "deadrootc_beetle",
        "deadrootn_beetle",
        "understory_gpsn","understory_resp", "understory_rootdepth", "understory_npp", "ratio_abg_litter", "litr1_hr", "litr2_hr", "litr3_hr", "litr4_hr");
    /*--------------------------------------------------------------*/
    /*	Yearly 							*/
    /*--------------------------------------------------------------*/

    char out_basic_basin_yearly[] = "%s %s %s %s %s %s %s %s %s %s \n";

#ifdef JMG_TRACKING
    char out_format_basin_yearly[1000] = "%s %s %s ";
    strcat(out_format_basin_yearly,out_basic_basin_yearly);
#else
    char out_format_basin_yearly[1000] = "";
    strcat(out_format_basin_yearly, out_basic_basin_yearly);
#endif

    outfile = world_output_files[0].basin[0].yearly;
    fprintf(outfile, out_format_basin_yearly,
        //"%s %s %s %s %s %s %s %s %s %s \n",
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "year",
        "basinID",
        "gpsn",
        "plantresp",
        "newC",
        "soilhr",
        "strN",
        "denitrif","root_depth","mortf"
            );

    }

    /*--------------------------------------------------------------*/
    /*	Hillslope file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].h != NULL) {
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

    char out_basic_hillslope_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_hillslope_daily[1000] = "%s %s %s ";
    strcat(out_format_hillslope_daily,out_basic_hillslope_daily);
#else
    char out_format_hillslope_daily[1000] = "";
    strcat(out_format_hillslope_daily, out_basic_hillslope_daily);
#endif

    outfile = world_output_files[0].hillslope[0].daily;

    fprintf(outfile, out_format_hillslope_daily,
            //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "day",
        "month",
        "year",
        "basinID",
        "lai",
        "gpsn",
        "plant_resp",
        "soil_resp",
        "nitrate",
        "sminn",
        "surfaceN",
        "plantc",
        "plantn",
        "cpool",
        "npool",
        "litrc",
        "litrn",
        "soilc",
        "soiln",
        "gwNO3",
        "gwNH4",
        "gwDON",
        "gwDOC",
        "streamflow_NO3",
        "streamflow_NH4",
        "streamflow_DON",
        "streamflow_DOC",
        "gwNO3out",
        "gwNH4out",
        "gwDONout",
        "gwDOCout",
        "denitrif",
        "nitrif",
        "DOC",
        "DON",
        "root_depth",
        "nfix",
        "nuptake",
        "grazingC",
        "StreamNO3_from_surface",
        "StreamNO3_from_sub", "litr_decomp", "decom_landclim_daily");
    }

    /*--------------------------------------------------------------*/
    /*	Zone file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].z != NULL) {
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

    char out_basic_zone_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_zone_daily[1000] = "%s %s %s ";
    strcat(out_format_zone_daily,out_basic_zone_daily);
#else
    char out_format_zone_daily[1000] = "";
    strcat(out_format_zone_daily, out_basic_zone_daily);
#endif

    outfile = world_output_files[0].zone[0].daily;
    fprintf(outfile, out_format_zone_daily,
            //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n ",
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "day",
        "month",
        "year",
        "basinID",
        "hillID",
        "ID",
        "rain",
        "snow",
        "tday",
        "tavg",
        "vpd",
        "Kdown_direct",
        "Kdown_diffuse",
        "PAR_direct",
        "PAR_diffuse");
    }

    /*--------------------------------------------------------------*/
    /*	Patch file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].p != NULL) {
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

//#ifdef LIU_BURN_ALL_AT_ONCE
    char out_basic_patch_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
//#else
//    char out_basic_patch_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
//#endif

#ifdef JMG_TRACKING
    char out_format_patch_daily[1000] = "%s %s %s ";
    strcat(out_format_patch_daily,out_basic_patch_daily);
#else
    char out_format_patch_daily[1000] = "";
    strcat(out_format_patch_daily, out_basic_patch_daily);
#endif

    outfile = world_output_files[0].patch[0].daily;
    check = fprintf(outfile, out_format_patch_daily,

#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "day",
        "month",
        "year",
        "basinID",
        "hillID",
        "zoneID",
        "patchID",
        "lai",
        "plantc",
//#ifdef LIU_BURN_ALL_AT_ONCE
        "tot_leafc",
        "tot_frootc",
        "tot_woodc",
//#endif
        "plantn",
        "net_psn",
        "plant_resp",
        "soil_resp",
        "litr1c",
        "litr2c",
        "litr3c",
        "litr4c",
        "litr1n",
        "litr2n",
        "litr3n",
        "litr4n",
        "lit.rain_cap",
        "soil1c",
        "soil2c",
        "soil3c",
        "soil4c",
        "soil1n",
        "soil2n",
        "soil3n",
        "soil4n",
        "soilDON",
        "soilDOC",
        "denitrif",
        "netleach",
        "DON_loss",
        "DOC_loss",
        "soilNO3",
        "soilNH4",
        "streamNO3",
        "streamNH4",
        "streamDON",
        "streamDOC",
        "surfaceNO3",
        "surfaceNH4",
        "surfaceDOC",
        "surfaceDON",
        "height",
        "nuptake",
        "root_depth",
        "nfix",
        "grazingC",
        "area");
    /*--------------------------------------------------------------*/
    /*	Yearly 							*/
    /*--------------------------------------------------------------*/

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic_patch_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#else
    char out_basic_patch_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#endif

#ifdef JMG_TRACKING
    char out_format_patch_yearly[1000] = "%s %s %s ";
    strcat(out_format_patch_yearly,out_basic_patch_yearly);
#else
    char out_format_patch_yearly[1000] = "";
    strcat(out_format_patch_yearly, out_basic_patch_yearly);
#endif

    outfile = world_output_files[0].patch[0].yearly;
    fprintf(outfile, out_format_patch_yearly,

#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "year",
        "basinID",
        "hillID",
        "zoneID",
        "patchID",
        "leaf_c",
        "leaf_n",
        "plant_c",
        "plant_n",
        "litter_c",
        "soil_c",
        "litter_n",
        "soil_n",
        "nitrate",
        "sminn","root_depth"
#ifdef JMG_MORE_YEARLY_OUTPUT
        ,"denitrif"
#endif
            );
    }

    /*--------------------------------------------------------------*/
    /*	Stratum file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].c != NULL) {
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

    char out_basic_stratum_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_stratum_daily[1000] = "%s %s %s ";
    strcat(out_format_stratum_daily,out_basic_stratum_daily);
#else
    char out_format_stratum_daily[1000] = "";
    strcat(out_format_stratum_daily, out_basic_stratum_daily);
#endif

    outfile = world_output_files[0].canopy_stratum[0].daily;
    fprintf(outfile, out_format_stratum_daily,
        //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "day",
        "month",
        "year",
        "basinID",
        "hillID",
        "zoneID",
        "patchID",
        "stratumID",
        "proj_lai",
        "proj_lai_when_red", //Beetle effect
        "proj_pai",
        "proj_pai_when_red",
        "toal_snag_c", //decaying pool + delaying pool
        "toal_snag_n", //decaying pool +delaying pool
        "total_redneedle_c",
        "total_redneedle_n",  //end NREN
        "deadrootc_beetle",
        "deadrootn_beetle", //NREN 20190910
        "tleafc",
        "leafc",
        "leafc_store",
        "leafc_transfer",
        "tleafn",
        "leafn",
        "leafn_store",
        "leafn_transfer",
        "cpool",
        "npool",
        "dead_leafc",
        "frootc",
        "frootn",
        "live_stemc",
        "live_stemn",
        "leafc_store",
        "leafn_store",
        "dead_stemc",
        "dead_stemn",
        "live_crootc",
        "live_crootn",
        "dead_crootc",
        "dead_crootn",
        "cwdc",
        "mresp",
        "gresp",
        "psn_to_cpool",
        "age",
        "root_depth",
        "gwseasonday",
        "lfseasonday",
        "gsi",
        "nlimit",
        "fleaf",
        "froot",
        "fwood",
        "Nuptake",
        "smin2pl",
        "retrans2pl",
        "mort_fract",
        "assim_sunlit",
        "assim_shade",
        "trans_sunlit",
        "trans_shade",
        "proj_lai_sunlit",
        "proj_lai_shade");

  /*--------------------------------------------------------------*/
    /* Shadow	Daily 			                                   				*/
    /*--------------------------------------------------------------*/
    if (command_line[0].vegspinup_flag > 0){

    char out_basic_shadow_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_shadow_daily[1000] = "%s %s %s ";
    strcat(out_format_shadow_daily,out_basic_shadow_daily);
#else
    char out_format_shadow_daily[1000] = "";
    strcat(out_format_shadow_daily, out_basic_shadow_daily);
#endif

    outfile = world_output_files[0].shadow_strata[0].daily;
      fprintf(outfile, out_format_shadow_daily,
        //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "day",
        "month",
        "year",
        "basinID",
          "hillID  ",
          "zoneID  ",
          "patchID  ",
          "stratumID  ",
          "proj_lai  ",
          "leafc  ",
          "leafn  ",
          "cpool  ",
          "npool  ",
          "dead_leafc  ",
          "frootc  ",
          "frootn  ",
          "live_stemc  ",
          "live_stemn  ",
          "leafc_store  ",
          "leafn_store  ",
          "dead_stemc  ",
          "dead_stemn  ",
          "live_crootc  ",
          "live_crootn  ",
          "dead_crootc  ",
          "dead_crootn  ",
          "cwdc  ",
          "mresp  ",
          "gresp  ",
          "psn_to_cpool","age","root_depth","gwseasonday","lfseasonday","gsi", "nlimit",
          "fleaf","froot","fwood","Nuptake","smin2pl","retrans2pl","mort_fract");
  }

    /*--------------------------------------------------------------*/
    /*	Yearly 							*/
    /*--------------------------------------------------------------*/

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic_stratum_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#else
    char out_basic_stratum_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#endif

#ifdef JMG_TRACKING
    char out_format_stratum_yearly[1000] = "%s %s %s ";
    strcat(out_format_stratum_yearly,out_basic_stratum_yearly);
#else
    char out_format_stratum_yearly[1000] = "";
    strcat(out_format_stratum_yearly, out_basic_stratum_yearly);
#endif

    outfile = world_output_files[0].canopy_stratum[0].yearly;
    fprintf(outfile, out_format_stratum_yearly,

#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

         "year",
         "basinID",
         "hillID",
         "zoneID",
         "patchID",
         "stratumID",
#ifndef JMG_MORE_YEARLY_OUTPUT
         "proj_lai",
         "leafc",
         "leafn",
         "frootc",
         "frootn",
         "stemc",
         "stemn",
         "cwdc",
         "cwdn",
         "psn_net",
         "mr",
         "gr",
         "minNSC",
         "mortfract",
         "snagc",
         "snagn",
         "redneedlec",
         "redneedlen",
         "dead_rootc_beetle",
         "dead_rootn_beetle",
         "height",
         "rootdepth_mm"
#else
        "AGBc",
        "height",
        "LAI", // actually epv.proj_lai (projected LAI to surface plane)
        "psn_net",
        "mr",
        "gr",
        "stemc",
        "leafc",
        "rootc",
        "rootdepth_mm",
        "rootc_rd"
#endif
            );
    }


    if (command_line[0].f !=NULL) {
      /*--------------------------------------------------------------*/
      /*	Yearly 							*/
      /*--------------------------------------------------------------*/
      outfile = world_output_files[0].fire[0].yearly;
      fprintf(outfile, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
           "year",
           "basinID",
           "hillID",
           "zoneID",
           "patchID",
           "stratumID",
           "proj_lai",
           "leafc",
           "leafn",
           "frootc",
           "frootn",
           "stemc",
           "stemn",
           "cwdc",
           "cwdn",
           "psn","cpool", "mortfract","height","rootdepth");
    }




    return;
} /*end add_growth_headers*/
