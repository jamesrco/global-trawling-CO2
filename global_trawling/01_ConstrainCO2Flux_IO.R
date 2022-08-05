# 01_ConstrainCO2Flux_IO.R
# Created June 7, 2022
# Purpose: First in series of scripts used to constrain the crude estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This first script reads in the necessary data files

# set the working directory; create directory for output

setwd("~/Code/global-trawling-CO2/")

# libraries

options("rgdal_show_exportToProj4_warnings"="none")
library(sp) # needs to be installed first
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal) # needs to be installed first
library(raster) # assumes you have some version of GDAL up and running, and 
                # will force you to load several dependencies, including terra 
library(data.table) # needs to be installed first
library(parallel) # part of base; doesn't need to be installed

# library(ncdf4) # not needed anymore

# load datasets

# # Attempt to load Siegel et al. "sequestration fractions" from NetCDF
# 
# OCIM2_fseq_48L <- nc_open("data/global_trawling/raw/siegel_et_al_2021_v2/fseq_OCIM2_48L.nc", verbose=TRUE) # currently returns an error

# instead we'll use a modified version of the MATLAB script provided by Siegel et al. to get what we need: see "gen_fracs_to_constrain_trawlCO2.m" which should be in this same directory

# take the necessary detour into MATLAB at this point, if the output hasn't been generated; then return to R

# the output from gen_fracs_to_constrain_trawlCO2.m -- sequestration fractions for ocean bottom depths, plus the necessary metadata -- should now be in several .csv files found in data/global_trawling/derived/benthic_seqfractions

# load the benthic sequestration fractions for entire Siegel et al. model domain,
# from 1-200 y, then 200-1000 y in 100 y increments
# also load the years

fseq_bottom_multYears.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears.mat")
fseq_bottom.multyears <- fseq_bottom_multYears.raw$fseq.bottom.multyears # clean up a bit
fseq_bottom.multyears[fseq_bottom.multyears>=1] <- 1

seqFracYears.raw <- read.csv(file = "data/global_trawling/derived/benthic_seqfractions/benthic_years.csv",
                             header = FALSE)

# load metadata

fseq_bottom_depth_m.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/bottom_depth_m.csv",
                                 header = FALSE)
fseq_bottom_lat_degN.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/lat_degN.csv",
                                 header = FALSE)
fseq_bottom_long_degE.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/long_degE.csv",
                                 header = FALSE)

# load Sala et al. pCO2 data, as GeoTIFF, then convert to data frames
# with some help from https://datacarpentry.org/r-raster-vector-geospatial/01-raster-structure/ and
# https://www.neonscience.org/resources/learning-hub/tutorials/raster-data-r

# good instructions on enabling multithreading on Mac (for data.table) here:
# https://github.com/Rdatatable/data.table/wiki/Installation and (even more helpful)
# here: https://firas.io/post/data.table_openmp/
# *** removing & recompiling data.table from source is critical if you already have it 
# installed

Sala_bottomtrawl_Ia.raw <- 
  raster("data/global_trawling/raw/sala_et_al_2021/bottom_trawling_Ia.tif")
Sala_bottomtrawl_Ia.df <- as.data.frame(Sala_bottomtrawl_Ia.raw, xy = TRUE)

Sala_carbon_ranking.raw <- 
  raster("data/global_trawling/raw/sala_et_al_2021/carbon_ranking.tif")
Sala_carbon_ranking.df <- as.data.frame(Sala_carbon_ranking.raw, xy = TRUE)

Sala_CO2_efflux.raw <- 
  raster("data/global_trawling/raw/sala_et_al_2021/co2_efflux.tif")
Sala_CO2_efflux.df <- as.data.frame(Sala_CO2_efflux.raw, xy = FALSE) 
# # note that creating a data frame from this last object will require a lot of memory ... I had to up 
# # the R_MAX_VSIZE variable in .Renviron according to the directions here:
# # https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos

# load the year-by-year predicted global CO2 sediment remineralization rates from Sala et al.
# (this is a static .csv version of the object "results," generated from the simple
# model beginning on line 119 in
# https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/ancillary_analyses/timing_of_trawling_impacts.Rmd)

Sala_et_al_trawlTiming_results.raw <- read.csv(file = "data/global_trawling/derived/sala_et_al_2021_model/Sala_et_al_trawlTiming_results.csv",
                                               header = TRUE) 

colnames(Sala_et_al_trawlTiming_results.raw)[1] <- c("Year")

# since the underlying Sala et al. dataset is very sparse, can create a subset
# and run our functions on just that subset

ind.nonZeroCO2 <- which(!is.na(Sala_CO2_efflux.df$co2_efflux))
length(ind.nonZeroCO2)/length(Sala_CO2_efflux.df$co2_efflux) # values that aren't NA represent < 1% of the total dataset

# need to define cell area per email from jmayorga@bren.ucsb.edu:
# "the fluxes in the Geotiff are per km2, so please make sure to multiply times
# the pixel’s area before summing up."

SalaModel_cell_area <- 934.4789^2/1000000 # from https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/data_prep/update_bottom_trawling_impact.Rmd

# now, can find most appropriate (nearest) sequestration fraction for each point in the CO2 flux dataset

# create data frames of the coordinates of the points in the two different datasets

# CO2 flux data
Sala_CO2_efflux.coords <- xyFromCell(Sala_CO2_efflux.raw,c(1:length(Sala_CO2_efflux.raw))) # easy, can just use the xyFromCell function in raster package 

# sequestration fractions
Siegel_fseq.coords.df <- data.frame(as.numeric(rep(fseq_bottom_long_degE.raw[1,],91)),rep(fseq_bottom_lat_degN.raw[,1],180))
colnames(Siegel_fseq.coords.df) <- c("x","y")
# convert to eastings & westings rather than just eastings, for compatibility 
Siegel_fseq.coords.df$x[Siegel_fseq.coords.df$x>180] <- Siegel_fseq.coords.df$x[Siegel_fseq.coords.df$x>180]-360

# continues onto 02_ConstrainCO2Flux_coordMatch.R, which is set up to perform
# some intensive calculations using an AWS ECS instance