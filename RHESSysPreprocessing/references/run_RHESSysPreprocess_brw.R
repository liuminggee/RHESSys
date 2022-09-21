#### Bull Run Model - Calibration - generate new worldfile and flow table ####

# run_RHESSysPreprocess
# Will Burke 3/5/18

# Instructions
# ------------
# This is an example script showing how the RHESSysPreprocess.R function should be run.
# 1) Install or source the RHESSysPreprocessing package
# 2) Copy this script, and edit where indicated.
# 3) Run the RHESSysPreprocess.R function at the bottom.
# 4) The funciton will produce:
#       - worldfile
#       - flowtable
#       - metadata (if not supressed)

# Load Package
# ------------
# Assuming you installed the RHESSysPreprocessing package, you will need to load it
library(RHESSysPreprocessing)

# Filepaths
# ----------------
# This script uses relative filepaths. This means that it will look folders and files relative to your current working directory.
# If needed, set your current working directory to the folder of your project:
setwd( "C:/Users/PETBUser/OneDrive - Washington State University (email.wsu.edu)/Research/BRW_project/Modeling/Create/BRW")
#dir()
#dir("./input")

# This script also uses the "~", which is a shorthand method of navigating to your "home" user directory - typically the folder named for your username.

# Spatial Data
# ------------
# You will need to select your method of geospatial data input.  This is the means by which the spatial data referenced in your template

# Currently there are a two supported methods:
# 1) Raster - spatial data in any raster format supported by R GDAL will be read in from a folder.
# 2) GRASS GIS - GRASS 6 or 7, spatial data will be imported from the specified GRASS location and mapset.

# NOTES:
# - Due to a variety of factors, spatial data import via the raster method is both more robust and faster.
# - Regardless of import method, good practice for spatial data should be followed.
# - Input data should have the same projections, extents, and cell sizes. THIS MAY RESULT IN ERRORS IF NOT FOLLOWED.

# Raster
# ------
# To import spatial data from a folder of rasters:
# 1) Set type to "raster"
type = "raster"
# 2) Set typepars to the path of the folder containing your rasters
typepars = "./input/rasters"

# Template
# --------
# The worldfile template is the key document that outlines how your worldfile will be built.
# The template variable should point to the name and location of your template.
#template = "./templates/br_cal_template.txt"
template = "./input/templates/brw_statefile_template.txt"

# Name
# ----
# Set the name and path for all function outputs.
# Suffixes of .world, .flow, and .meta will be appended to the worldfile, flowtable, and metadata files respectively.
name = "./output/brw_test"

# Overwrite
# ---------
# TRUE/FALSE if an existing worldfile and flowtable be overwritten.
overwrite = TRUE

# Streams > also uses relative path as typepars
# -------
# Streams map to be used in creation of the flowtable - this is just the name of the map, to be found via the method indicated with "type"
streams = "streams.tif"

# Optional Flowtable Spatial Data > also use relative path as typepars
# -------------------------------
# These maps are optional inputs in flowtable creation
roads = "roads.tif"
impervious = "impervious_flip.tif"
# roofs = "roofs_map"

# Header
# ------
# TRUE/FALSE to produce a header file. Header file will be have same name(and location) set by "name", with the ".hdr" suffix.
header = TRUE

# Finally, run the function.  Depending on size, it may take a minute or two.
RHESSysPreprocess(
  template = template,
  name = name,
  type = type,
  typepars = typepars,
  streams = streams,
  overwrite = overwrite,
  header = header)
  #2


# CreateFlownet(name = name, readin = template, type = "raster",
#               typepars = typepars, streams = streams, roads = roads, road_width = 1,
#               impervious = impervious, overwrite = T)



