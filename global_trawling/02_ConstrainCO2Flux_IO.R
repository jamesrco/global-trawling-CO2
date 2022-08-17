# 02_ConstrainCO2Flux_IO.R
# Created June 7, 2022
# Purpose: Second in series of scripts used to constrain the crude estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This second script reads in the necessary starting data files, does some basic
# data characterization/exploration

# set the working directory; create directory for output

setwd("~/Code/global-trawling-CO2/")

# libraries

options("rgdal_show_exportToProj4_warnings"="none")
library(sp) # needs to be installed first
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal) # needs to be installed first
library(sf)
library(raster) # assumes you have some version of GDAL up and running, and 
                # will force you to load several dependencies, including terra 
library(data.table) # needs to be installed first
library(parallel) # part of base; doesn't need to be installed
library(R.matlab) # to read .mat file

# load datasets

# metadata from the Siegel et al. OCIM48 model; these files should have been
# generated using the MATLAB script "01_ConstrainCO2Flux_genSiegelMetaData.m,"
# which should be in this same directory

# load metadata and model domain characteristics

fseq_bottom_depth_m.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/bottom_depth_m.csv",
                                 header = FALSE)
fseq_lat_degN.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/lat_degN.csv",
                                 header = FALSE)
fseq_long_degE.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/long_degE.csv",
                                 header = FALSE)
fseq_surfaceLandMask.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/mask_surface.csv",
                                      header = FALSE)
fseq_OCIM_modelDepths.raw <- read.csv("data/global_trawling/derived/benthic_seqfractions/OCIM_modelDepths.csv",
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

# some exploration of the data

# load depths corresponding to the locations of all non-zero data points in the 
# Sala et al. efflux data (Sala_CO2_efflux.df)

Sala_CO2efflux_nonZeroValueDepths_raw <- st_read("data/global_trawling/derived/Sala_CO2efflux_nonZeroValueDepths/Sala_CO2efflux_nonZeroValueDepths.shp")

# adjust for values > 0 (something strange happened with the very near coastal
# points, perhaps as a result of the reprojection)
Sala_CO2efflux_nonZeroValueDepths_adj <- Sala_CO2efflux_nonZeroValueDepths_raw
Sala_CO2efflux_nonZeroValueDepths_adj$rvalue_1[Sala_CO2efflux_nonZeroValueDepths_adj$rvalue_1>0] <- 0

# take a quick look at a histogram; at which depths are most of the non-zero values?
effluxEstHist_numObs <- hist(Sala_CO2efflux_nonZeroValueDepths_adj$rvalue_1, breaks = seq(-3400,0,100),
     xlab = "Water column depth (m)",
     ylab = expression(paste("No. observations of est. CO"["2"]," efflux")),
     main = expression(paste("Histogram of CO"["2"]," efflux estimates from Sala et al. trawling disturbance model")))

# what about a similar plot that sums things by mass instead
bins <- cut(Sala_CO2efflux_nonZeroValueDepths_adj$rvalue_1, breaks = seq(-3400,0,100))
sums_PgCO2 <- tapply(Sala_CO2_efflux.df$co2_efflux[ind.nonZeroCO2], bins, sum)*
SalaModel_cell_area*(1/10^9)

effluxEstHist_mass <- effluxEstHist_numObs
effluxEstHist_mass$counts <- sums_PgCO2

plot(effluxEstHist_mass,
     xlim = c(-3500,0),
     ylim = c(0,1),
     xlab = "Water column depth at location of estimate (m)",
     ylab = expression(paste("Sala et al. est. CO"["2"]," efflux (Pg CO"["2"],")")),
     main = expression(paste("Total estimated CO"["2"]," efflux from benthic trawling disturbance, by water column depth")),
     yaxt = "n")

axis(2, at = seq(0,1,0.2))

sum(sums[34])/sum(sums, na.rm = T) # % total CO2 efflux (by mass) in 0-100 m depth bin
sum(sums[33:34])/sum(sums, na.rm = T) # % total CO2 efflux (by mass) in 0-200 m depth bin
sum(sums[1:32], na.rm = T)/sum(sums, na.rm = T) # % total CO2 efflux (by mass) originating in water depths > 200 m

# now, can find most appropriate (nearest) sequestration fraction for each point in the CO2 flux dataset

# create data frames of the coordinates of the points in the two different datasets

# CO2 flux data
Sala_CO2_efflux.coords <- xyFromCell(Sala_CO2_efflux.raw,c(1:length(Sala_CO2_efflux.raw))) # easy, can just use the xyFromCell function in raster package 
# make a subset since these data are sparse
Sala_CO2_efflux.coords.nonZero <- Sala_CO2_efflux.coords[ind.nonZeroCO2,]

# sequestration fractions
Siegel_fseq.coords.df <- data.frame(as.vector(apply(fseq_long_degE.raw,2,rep,91)),rep(fseq_lat_degN.raw[,1],180))
colnames(Siegel_fseq.coords.df) <- c("x","y")
# convert to eastings & westings rather than just eastings, for compatibility 
Siegel_fseq.coords.df$x[Siegel_fseq.coords.df$x>180] <- Siegel_fseq.coords.df$x[Siegel_fseq.coords.df$x>180]-360
# create another column that includes the index
Siegel_fseq.coords.df$ind <- seq(1:nrow(Siegel_fseq.coords.df))
# subset to just points which are in the OCIM ocean
ind_oceanPoints <- which(as.numeric(unlist(fseq_surfaceLandMask.raw))==1)
Siegel_fseq.coords.df.oceanOnly <- Siegel_fseq.coords.df[ind_oceanPoints,]

# continues onto 03_ConstrainCO2Flux_coordMatch.R