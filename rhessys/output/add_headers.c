/*--------------------------------------------------------------*/
/* 																*/
/*					add_headers					*/
/*																*/
/*	add_headers -
                                                            */
/*	NAME														*/
/*	add_headers
                                                                */
/*	SYNOPSIS													*/
/*	void add_headers(struct world output_file_object *,				*/
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


void add_headers(struct world_output_file_object *world_output_files,
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
    outfile = world_output_files[0].basin[0].hourly;
    fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n",
    // the unit is based on mm and day
        "hour",
        "day",
        "month",
        "year",
        "basinID",
        "pot_surface_infil",
        //asnow_throughfall * 1000.0,
        "sat_def_z",
        "sat_def",
        "rz_stor",
        "unsat_stor",
        "rz_drainage",
        "unsat_drainage",
        //acap_rise * 1000.0,
        //aevaporation * 1000.0,
        //asnowpack * 1000.0,
        //atranspiration * 1000.0,
        "subsur2stream_flow",
        "sur2stream_flow",
        "streamflow",
        //apsn,
        //alai,
        "gw.Qout",
        "gw.storage",
        "detention_store",
        "%sat_area",
        "litter_store",
        "canopy_store",
        //aperc_snow *100,
        //asublimation * 1000.0,
        //var_trans,
        //aacctrans*1000,
        //var_acctrans,
        //aPET*1000,
        //adC13,
        "precip",
        //amortality_fract*100,
        //atmax,
        //atmin,
        //asnow*1000.0 ,
        "routedstreamflow");











    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

    char out_basic_basin_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_basin_daily[1000] = "%s %s %s ";
    strcat(out_format_basin_daily,out_basic_basin_daily);
#else
    char out_format_basin_daily[1000] = "";
    strcat(out_format_basin_daily, out_basic_basin_daily);
#endif

    outfile = world_output_files[0].basin[0].daily;
    fprintf(outfile, out_format_basin_daily,
            //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
#ifdef JMG_TRACKING
            "sday",
            "smth",
            "syr",
#endif
        "day",
        "month",
        "year",
        "basinID",
        "pot_surface_infil",
        "snow_thr",
        "sat_def_z",
        "sat_def",
        "rz_storage",
        "unsat_stor",
        "rz_drainage",
        "unsat_drain",
        "cap",
        "evap",
        "snowpack",
        "trans",
        "baseflow",
        "return",
        "streamflow",
        "psn",
        "nppcum",
        "lai",
        "gw.Qout",
        "gw.storage",
        "detention_store",
        "%sat_area",
        "litter_store",
        "canopy_store",
        "%snow_cover",
        "snow_subl",
        "trans_var",
        "acc_trans",
        "acctransv_var",
        "tpet",
        "pet",
        "pe",
        "dC13",
        "precip",
        "pcp_assim",
        "mortf",
        "tmax",
        "tmin",
        "tavg",
        "vpd",
        "snowfall",
        "recharge",
        "gpsn",
        "resp",
        "gs",
        "rootdepth",
        "plantc",
        "snowmelt",
        "canopysubl",
        "routedstreamflow",
        "canopy_snow",
        "height",
        "evap_can","evap_lit","evap_soil",
        "litrc",
        "Kdown","Ldown","Kup","Lup",
        "Kstar_can","Kstar_soil","Kstar_snow",
        "Lstar_can","Lstar_soil","Lstar_snow",
        "LE_canopy","LE_soil","LE_snow","Lstar_strat","canopydrip","ga", "litr_decomp", "decom_landclim_daily");

    /*--------------------------------------------------------------*/
    /*	Monthly							*/
    /*--------------------------------------------------------------*/

#ifdef LIU_TRACKING_BASIN_LITTERC
    char out_basic_basin_monthly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#else
    char out_basic_basin_monthly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#endif

#ifdef JMG_TRACKING
    char out_format_basin_monthly[1000] = "%s %s %s ";
    strcat(out_format_basin_monthly,out_basic_basin_monthly);
#else
    char out_format_basin_monthly[1000] = "";
    strcat(out_format_basin_monthly, out_basic_basin_monthly);
#endif

    outfile = world_output_files[0].basin[0].monthly;
    check = fprintf(outfile, out_format_basin_monthly,

#ifdef JMG_TRACKING
                    "sday",
                    "smth",
                    "syr",
#endif

        "month",
        "year",
        "basinID",
        "streamflow",
        "streamflow_NO3",
        "denitrif",
        "DOC",
        "DON",
        "et",
        "psn",
        "lai",
        "nitrif",
        "mineralized",
        "uptake"
#ifdef LIU_TRACKING_BASIN_LITTERC
        ,
        "leafc_to_litrc",
        "frootc_to_litrc",
        "cwdc_to_litrc",
        "stemc_to_litrc",
        "mort_to_litrc",
        "do_litrc_loss",
        "m_litrc_to_atmos",
        "litrc_to_atmos",
        "litrc_to_soilc"
#endif
        );
    /*--------------------------------------------------------------*/
    /*	Yearly 							*/
    /*--------------------------------------------------------------*/

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic_basin_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#else
    char out_basic_basin_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#endif

#ifdef JMG_TRACKING
    char out_format_basin_yearly[1000] = "%s %s %s ";
    strcat(out_format_basin_yearly,out_basic_basin_yearly);
#else
    char out_format_basin_yearly[1000] = "";
    strcat(out_format_basin_yearly, out_basic_basin_yearly);
#endif

    outfile = world_output_files[0].basin[0].yearly;
    check = fprintf(outfile, out_format_basin_yearly,

#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "year",
        "basinID",
#ifndef JMG_MORE_YEARLY_OUTPUT
        "streamflow",
        "streamflow_NO3",
        "denitrif",
        "DOC",
        "DON",
        "et",
        "psn","lai","nitrif",
        "mineralized", "uptake", "num_thresh","tpet","pet","pe"
#else
                    "pcp",
                    "et",
                    "TPET",
                    "PET",
                    "PE",
                    "streamflow",
                    "hill_base_flow",
                    "gw_drainage",
                    "rz_storage",
                    "unsat_storage",
                    "hill_gw_storage",
                    "n_deposition",
                    "denitrif",
                    "soilc",
                    "soiln",
                    "litrc",
                    "litrn",
                    "plantc",
                    "plantn",
                    "AGBc",
                    "BGBc",
                    "lai",
                    "sat_deficit",
                    "sat_deficit_z"
#endif
                    );
    }

    /*--------------------------------------------------------------*/
    /*	Hillslope file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].h != NULL) {
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

    char out_basic_hillslope_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_hillslope_daily[1000] = "%s %s %s ";
    strcat(out_format_hillslope_daily,out_basic_hillslope_daily);
#else
    char out_format_hillslope_daily[1000] = "";
    strcat(out_format_hillslope_daily, out_basic_hillslope_daily);
#endif

    outfile = world_output_files[0].hillslope[0].daily;
    fprintf(outfile, out_format_hillslope_daily,
            //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
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
        "rain_thr",
        "snow_thr",
        "sat_def_z",
        "sat_def",
        "unsat_stor",
        "unsat_drain",
        "cap",
        "evap",
        "snow",
        "trans",
        "baseflow",
        "return",
        "streamflow",
        "psn",
        "lai",
        "gw.Qout",
        "gw.storage",
        "precip", //NREN 2018/12/7
        "evap_surface",
        "soil_evap",
        "rz_storage",
        "detention_stor",
        "rain_stor",
        "litter_stor",
        "area",
        "pet", //REN 2019/11/29
        "snowpack_sublim",
        "canopy_subl",
        "height",
        "woodc",
        "lai_red",
        "gsurf",
        "potential_exfil",
        "soil_potential_evap",
        "rootzone.S",
        "gs",
        "ga",
        "exfiltration_sat",
        "exfiltration_unsat", "tmax", "tmin"
        );

    /*--------------------------------------------------------------*/
    /*	Monthly							*/
    /*--------------------------------------------------------------*/

    char out_basic_hillslope_monthly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_hillslope_monthly[1000] = "%s %s %s ";
    strcat(out_format_hillslope_monthly,out_basic_hillslope_monthly);
#else
    char out_format_hillslope_monthly[1000] = "";
    strcat(out_format_hillslope_monthly, out_basic_hillslope_monthly);
#endif

    outfile = world_output_files[0].hillslope[0].monthly;
    check = fprintf(outfile, out_format_hillslope_monthly,
        //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif
        "month",
        "year",
        "basinID",
        "hillslopeID",
        "streamflow",
        "streamflow_NO3",
        "snowpack",
        "denitrif",
        "DOC",
        "DON",
        "et",
        "psn",
        "lai",
        "nitrif",
        "mineralized",
        "uptake","area");
    /*--------------------------------------------------------------*/
    /*	Yearly 							*/
    /*--------------------------------------------------------------*/

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic_hillslope_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#else
    char out_basic_hillslope_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#endif

#ifdef JMG_TRACKING
    char out_format_hillslope_yearly[1000] = "%s %s %s ";
    strcat(out_format_hillslope_yearly,out_basic_hillslope_yearly);
#else
    char out_format_hillslope_yearly[1000] = "";
    strcat(out_format_hillslope_yearly, out_basic_hillslope_yearly);
#endif

    outfile = world_output_files[0].hillslope[0].yearly;
    check = fprintf(outfile, out_format_hillslope_yearly,

#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "year",
        "basinID",
        "hillslopeID",
#ifndef JMG_MORE_YEARLY_OUTPUT
        "streamflow",
        "streamflow_NO3",
        "denitrif",
        "DOC",
        "DON",
        "et",
        "psn","lai","nitrif",
                    "mineralized", "uptake",
#else
        "pch_pcp",
        "pch_et",
        "pch_streamflow",
        "pch_return_flow",
        "pch_base_flow",
        "hill_base_flow",
        "pch_gw_drainage",
        "pch_rz_storage",
        "pch_unsat_storage",
        "hill_gw_storage",
#endif
        "area");
    }
    /*--------------------------------------------------------------*/
    /*	Zone file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].z != NULL) {
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/

    char out_basic_zone_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_zone_daily[1000] = "%s %s %s ";
    strcat(out_format_zone_daily,out_basic_zone_daily);
#else
    char out_format_zone_daily[1000] = "";
    strcat(out_format_zone_daily, out_basic_zone_daily);
#endif

    outfile = world_output_files[0].zone[0].daily;
    fprintf(outfile, out_format_zone_daily,
            //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
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
        "tmax",
        "tmin",
        "vpd",
        "Kdown_direct",
        "Kdown_diffuse",
        "PAR_direct",
        "PAR_diffuse",
        "Ldown",
        "relH","aspect","z","slope","ehr","whr",
        "tdew","edew",
        "transmis",
        "wind",
        "deltaT","clearskytransmis","tcoeff1","cloudfrac","CO2_ppm");

    /*--------------------------------------------------------------*/
    /*	Monthly							*/
    /*--------------------------------------------------------------*/

    char out_basic_zone_monthly[] = "%s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_zone_monthly[1000] = "%s %s %s ";
    strcat(out_format_zone_monthly,out_basic_zone_monthly);
#else
    char out_format_zone_monthly[1000] = "";
    strcat(out_format_zone_monthly, out_basic_zone_monthly);
#endif

    outfile = world_output_files[0].zone[0].monthly;
    check = fprintf(outfile, out_format_zone_monthly,
        //"%s %s %s %s %s %s %s %s %s %s\n" ,
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "month",
        "year",
        "basinID",
        "hillID",
        "zoneID",
        "precip",
        "K_direct",
        "K_diffuse",
        "tmax", "tmin");

    /*--------------------------------------------------------------*/
    /*	Hourly 							*/
    /*--------------------------------------------------------------*/
    outfile = world_output_files[0].zone[0].hourly;
    fprintf(outfile,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n " ,
        "day",
        "month",
        "year",
        "hour",
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

    char out_basic_patch_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_patch_daily[1000] = "%s %s %s ";
    strcat(out_format_patch_daily,out_basic_patch_daily);
#else
    char out_format_patch_daily[1000] = "";
    strcat(out_format_patch_daily, out_basic_patch_daily);
#endif

    outfile = world_output_files[0].patch[0].daily;
        check = fprintf(outfile, out_format_patch_daily,
                        //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
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
                        "rain_thr",
                        "detention_store",
                        "sat_def_z",
                        "sat_def",
                        "rz_storage",
                        "potential_rz_store",
                        "rz_field_capacity",
                        "rz_wilting_point",
                        "unsat_stor",
                        "rz_drainage",
                        "unsat_drain",
                        "sublimation",
                        "return",
                        "evap",
                        "evap_surface",
                        "soil_evap",
                        "snow",
                        "snow_melt",
                        "trans_sat",
                        "trans_unsat",
                        "Qin",
                        "Qout",
                        "psn",
                        "root_zone.S",
                        "root.depth",
                        "litter.rain_stor",
                        "litter.S","area","tpet","pet","pe","lai","baseflow","streamflow","pcp","recharge",
                        "Kdowndirpch","Kdowndiffpch",
                        "Kupdirpch","Kupdifpch","Luppch",
                        "Kdowndirsubcan","Kdowndifsubcan","Ldownsubcan",
                        "Kstarcan","Kstardirsno","Kstardiffsno",
                        "Lstarcanopy","Lstarsnow","Lstarsoil",
                        "wind","windsnow","windzone","ga","gasnow","trans_reduc_perc","pch_field_cap",
                        "overland_flow","height","ustar","snow_albedo",
                        "Kstarsoil","Kdowndirsurf","Kdowndifsurf","exfil_unsat",
                        "snow_Rnet","snow_QLE","snow_QH","snow_Qrain","snow_Qmelt",
                        "LEcanopy",
                        "SED","snow_age", "psi", "t_scalar", "w_scalar",
                         "rate_scalar", "litr_decomp", "rate_landclim_year","rate_landclim_daily", "et_decom_mean", "Tsoil");

    /*--------------------------------------------------------------*/
    /*	Monthly							*/
    /*--------------------------------------------------------------*/

    char out_basic_patch_monthly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_patch_monthly[1000] = "%s %s %s ";
    strcat(out_format_patch_monthly,out_basic_patch_monthly);
#else
    char out_format_patch_monthly[1000] = "";
    strcat(out_format_patch_monthly, out_basic_patch_monthly);
#endif

    outfile = world_output_files[0].patch[0].monthly;
    check = fprintf(outfile, out_format_patch_monthly,
        //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",
#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "month",
        "year",
        "basinID",
        "hillID",
        "zoneID",
        "patchID",
        "leach",
        "denitrif",
        "soil_moist_deficit",
        "et",
        "psn",
        "DOC",
        "DON","lai","nitrif","mineralized","uptake","theta","snow","area","nitrate","sminn", "burn");
    /*--------------------------------------------------------------*/
    /*	Yearly							*/
    /*--------------------------------------------------------------*/

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic_patch_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#else
    char out_basic_patch_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
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
#ifndef JMG_MORE_YEARLY_OUTPUT
            "num_threshold_sat_def",
            "peaksweday",
            "meltday",
            "pklai",
            "peaklaiday",
            "leach",
            "denitrif",
            "DOC_loss",
            "DON_loss",
            "psn", "trans",
            "et","lai","nitrif","mineralized",
            "uptake","Theta","sd",
            "pkswe", "pktrans", "pkpet", "streamflow", "Qin","Qout","rec.wyd","rec.pet.wyd",
            "ndays_sat", "ndays_sat70", "midsm_wyd",
            "area","tpet","pet","pe","snowin_pcp","pot_recharge","recharge","recharge.wyd","pot_recharge.wyd"
#else
            "precip",
            "streamflow",
            "base_flow",
            "return_flow",
            "rz_storage",
            "unsat_storage",
            "gw_drainage",
            "overland_flow", // JMG09122022
            "et",
            "nday_sat",
            "nday_sat70",
            "plantc",
            "plantn",
            "litrc",
            "litrn",
            "soilc",
            "soiln",
            "lai",
            "psn",
            "uptake",
            "leach",
            "DON_loss",
            "denitrif",
            "nitrif",
            "mineralized",
            "sat_deficit",
            "sat_deficit_z"
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

    char out_basic_stratum_daily[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_stratum_daily[1000] = "%s %s %s ";
    strcat(out_format_stratum_daily,out_basic_stratum_daily);
#else
    char out_format_stratum_daily[1000] = "";
    strcat(out_format_stratum_daily, out_basic_stratum_daily);
#endif

    outfile = world_output_files[0].canopy_stratum[0].daily;
    fprintf(outfile, out_format_stratum_daily,
        //"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
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
        "lai",
        "proj_lai_when_red",
        "evap",
        "APAR_direct",
        "APAR_diffuse",
        "sublim",
        "trans",
        "ga",
        "gsurf",
        "gs",
        "psi",
        "leaf_day_mr",
        "psn_to_cpool",
        "rain_stored",
        "snow_stored",
        "rootzone.S",
        "m_APAR","m_tavg","m_LWP","m_CO2","m_tmin","m_vpd","dC13",
        "Kstar_dir","Kstar_dif",
        "Lstar","surf_heat",
        "height","covfrac","vegID", "wstress_days", "potential_psn_to_cpool");
    /*--------------------------------------------------------------*/
    /*	Monthly							*/
    /*--------------------------------------------------------------*/

    char out_basic_stratum_monthly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s\n";

#ifdef JMG_TRACKING
    char out_format_stratum_monthly[1000] = "%s %s %s ";
    strcat(out_format_stratum_monthly,out_basic_stratum_monthly);
#else
    char out_format_stratum_monthly[1000] = "";
    strcat(out_format_stratum_monthly, out_basic_stratum_monthly);
#endif

    outfile = world_output_files[0].canopy_stratum[0].monthly;
    fprintf(outfile, out_format_stratum_monthly,
            //"%s %s %s %s %s %s %s %s %s %s %s %s %s\n",

#ifdef JMG_TRACKING
        "sday",
        "smth",
        "syr",
#endif

        "month",
        "year",
        "basinID",
        "hillID",
        "zoneID",
        "patchID",
        "stratumID",
        "veg_parm_ID",
        "above_plantc",
        "cover_fraction",
        "height",
        "all_lai",
        "proj_lai");
    /*--------------------------------------------------------------*/
    /*	Yearly							*/
    /*--------------------------------------------------------------*/

#ifdef JMG_MORE_YEARLY_OUTPUT
    char out_basic_stratum_yearly[] = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n";
#else
    char out_basic_stratum_yearly[] = "%s %s %s %s %s %s %s %s %s\n";
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
            //"%s %s %s %s %s %s %s %s %s\n",
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
        "stratumID"

#ifndef JMG_MORE_YEARLY_OUTPUT
        ,"lai",
        "psn",
        "lwp"
#else
            ,"gpp",
            "resp",
            "npp",
            "AGBc",
            "BGBc",
            "LAI",
            "height",
            "rootdepth"
#endif
            );
}


    /*--------------------------------------------------------------*/
    /*	Fire file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].f != NULL) {
    /*--------------------------------------------------------------*/
    /*	Daily 							*/
    /*--------------------------------------------------------------*/
/*	outfile = world_output_files[0].fire[0].daily;
    fprintf(outfile,
        "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
        "day",
        "month",
        "year",
        "basinID",
        "hillID",
        "zoneID",
        "patchID",
        "stratumID",
        "vegID",
        "m_cwdc_to_atmos",
        "m_cwdn_to_atmos",
        "canopy_target_height",
        "canopy_target_height_u_prop",
        "canopy_target_prop_mort",
        "canopy_target_prop_mort_consumed",
        "canopy_target_prop_mort_u_component",
        "canopy_target_prop_mort_o_component",
        "canopy_target_prop_c_consumed",
        "canopy_target_prop_c_remain",
        "canopy_target_prop_c_remain_adjusted",
        "canopy_target_prop_c_remain_adjusted_leafc",
        "canopy_subtarget_height",
        "canopy_subtarget_height_u_prop",
        "canopy_subtarget_prop_mort",
        "canopy_subtarget_prop_mort_consumed",
        "canopy_subtarget_prop_c_consumed",
        "canopy_subtarget_c",
        "understory_c_consumed"); */


    /*--------------------------------------------------------------*/
    /*	yearly 							*/
    /*--------------------------------------------------------------*/
    outfile = world_output_files[0].fire[0].yearly;
    fprintf(outfile,
         "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" ,
         "year",
         "basinID",
         "hillID",
         "zoneID",
         "patchID",
         "stratumID",
         "vegID",
         "m_cwdc_to_atmos",
         "m_cwdn_to_atmos",
         "canopy_target_height",
         "canopy_target_height_u_prop",
         "canopy_target_prop_mort",
         "canopy_target_prop_mort_consumed",
         "canopy_target_prop_mort_u_component",
         "canopy_target_prop_mort_o_component",
         "canopy_target_prop_c_consumed",
         "canopy_target_prop_c_remain",
         "canopy_target_prop_c_remain_adjusted",
         "canopy_target_prop_c_remain_adjusted_leafc",
         "canopy_subtarget_height",
         "canopy_subtarget_height_u_prop",
         "canopy_subtarget_prop_mort",
         "canopy_subtarget_prop_mort_consumed",
         "canopy_subtarget_prop_c_consumed",
         "canopy_subtarget_biomassc", "litter_c_consumed",
         "understory_biomassc_consumed", "undersotry_leafc_consumed", "understory_stemc_consumed", "understory_rootc_consumed",
         "overstory_biomssc_consumed", "overstory_leafc_consumed", "overstory_stemc_consumed", "overstory_rootc_consumed",
         "overstory_biomassc_mortality", "overstory_leafc_mortality", "overstory_stemc_mortality", "overstory_rootc_mortality",
          "acc_length", "acc_length_overstory","acc_length_understory");

    }


    /*--------------------------------------------------------------*/
    /*	Stream routing file headers					*/
    /*--------------------------------------------------------------*/
    if (command_line[0].stro != NULL) {
        /*--------------------------------------------------------------*/
        /*	Daily 							*/
        /*--------------------------------------------------------------*/

        outfile = world_output_files[0].stream_routing[0].daily;
        fprintf(outfile, "%s %s %s %s %s %s %s %s %s\n",
                "day",
                "month",
                "year",
                "reachID",
                "Qout",
                "lateralinput",
                "Qin",
                "waterdepth",
                "reservoir.store");
    }
    return;
} /*end add_headers*/
