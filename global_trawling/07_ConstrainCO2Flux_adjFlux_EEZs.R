# 07_ConstrainCO2Flux_adjFlux_EEZs.R
# Created August 4, 2022
# Purpose: Seventh in series of scripts used to constrain the estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This seventh script performs the same calculations as in
# 05_ConstrainCO2Flux_adjFlux_extended.R (adjustment of the Sala et al
# benthic CO2 flux data, with time-integrated estimates of total emissions to
# the atmosphere for a variety of time horizons), but for specific
# exclusive economic zones (EEZs)

# *** Assumes user has already run the previous scripts in this series, and that
# the R objects generated using those scripts are in the environment already

# set the working directory

setwd("~/Code/global-trawling-CO2/") # for my laptop
# setwd("~/global-trawling-CO2/") # for AWS

# libraries (if not already loaded)

options("rgdal_show_exportToProj4_warnings"="none")
library(sp) # needs to be installed first
options("rgdal_show_exportToProj4_warnings"="none")
library(sf)
library(rgdal)
library(raster) # assumes you have some version of GDAL up and running, and 
# will force you to load several dependencies, including terra 
library(data.table)
library(parallel) # part of base; doesn't need to be installed
library(R.matlab) # to read .mat file

# if not already loaded, load in the file containing the coordinate matches,
# generated in previous script 02_ConstrainCO2Flux_coordMatch.R and saved to 
# global-trawling-CO2/data/global_trawling/derived/output/

load("data/global_trawling/derived/output/coord.matches.RData")

# need to load some shapefiles that define the EEZs

FMI_EEZ_polygons_v11_raw <- st_read("data/global_trawling/raw/FMI_world_EEZ_boundaries_v11/eez_v11.shp")

# subset to countries of interest
# good reference to make sure we capture the top countries, in additon to the others
# we are curious about for various reasons:
# Steadman et al., 2021, New perspectives on an old fishing practice: Scale, 
# context and impacts of bottom trawling; available at
# https://oursharedseas.com/wp-content/uploads/2021/12/HI-RES-REPORT-‘New-perspectives-on-an-old-fishing-practice.pdf

# define a list of countries we're interested in

trawlEEZs <- c("Argentina",
               "Canada",
               "China",
               "Denmark",
               "France",
               "Germany",
               "Greenland",
               "Iceland",
               "India",
               "Indonesia",
               "Ireland",
               "Italy",
               "Japan",
               "Malaysia",
               "Morocco",
               "Netherlands",
               "Norway",
               "Philippines",
               "Russia",
               "South Korea",
               "Sweden",
               "United Kingdom",
               "United States",
               "Vietnam")

# define function to subset the data by EEZ

EEZ_Argentina <- FMI_EEZ_polygons_v11_raw[FMI_EEZ_polygons_v11_raw$SOVEREIGN1=="Argentina",]

st_crop(harv_dtm, harv_boundary)

# run the calculations

# set up structure to hold results

predicted.PgCO2_per_year_to_atmos_EEZs <- as.data.frame(matrix(data = NA, 
                                                          nrow = length(seqFracYears.raw),
                                                          ncol = 1+length(trawlEEZs)))
colnames(predicted.PgCO2_per_year_to_atmos_EEZs) = c("Year",
                                                     trawlEEZs)
predicted.PgCO2_per_year_to_atmos_EEZs[,1] <- unlist(seqFracYears.raw)

# iterate

for (i in 1:nrow(predicted.PgCO2_per_year_to_atmos_EEZs)) {
  
  print(predicted.PgCO2_per_year_to_atmos_EEZs[i,1])
  
  time0 <- Sys.time()
  
  thisYear <- predicted.PgCO2_per_year_to_atmos_EEZs[i,1]
  EffluxFracs.thisyear <- genEffluxFracs(thisYear)
  adjCO2efflux.thisyear.global <- unlist(lapply(ind.nonZeroCO2, constrainFlux, EffluxFracs.thisyear))
  
  
  
  predicted.PgCO2_per_year_to_atmos[i,2] <- sum(adjCO2efflux.thisyear*SalaModel_cell_area, na.rm=T)*(1/10^9)
  
  time1 <- Sys.time()
  print(time1 - time0)
  
}
