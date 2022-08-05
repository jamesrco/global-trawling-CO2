# 04_ConstrainCO2Flux_adjFlux_EEZs.R
# Created August 4, 2022
# Purpose: Fourth in series of scripts used to constrain the estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This fourth script performs the same calculations as in
# 03_ConstrainCO2Flux_adjFlux_extended.R (adjustment of the Sala et al
# benthic CO2 flux data, with time-integrated estimates of total emissions to
# the atmosphere for a variety of time horizons), but for specific
# exclusive economic zones (EEZs)

# *** Assumes user has already run 01_ConstrainCO2Flux_IO.R (in current session)
# and that the object "coord.matches.RData" generated using
# 02_ConstrainCO2Flux_coordMatch.R (likely via AWS) is present in
# global-trawling-CO2/data/global_trawling/derived/output/

# set the working directory

setwd("~/Code/global-trawling-CO2/") # for my laptop
# setwd("~/global-trawling-CO2/") # for AWS

# libraries (if not already loaded)

options("rgdal_show_exportToProj4_warnings"="none")
library(sp) # needs to be installed first
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal) # needs to be installed first
library(raster) # assumes you have some version of GDAL up and running, and 
# will force you to load several dependencies, including terra 
library(data.table) # needs to be installed first
library(parallel) # part of base; doesn't need to be installed
library(R.matlab) # to read .mat file

# if not already loaded, load in the file containing the coordinate matches,
# generated in previous script 02_ConstrainCO2Flux_coordMatch.R and saved to 
# global-trawling-CO2/data/global_trawling/derived/output/

load("data/global_trawling/derived/output/coord.matches.RData")

# need to load some shapefiles that define the EEZs



# generalized version of functions

genEffluxFracs <- function(year){
  yearInd <- which(seqFracYears.raw==year)
  fseq_bottom.thisyear <- fseq_bottom.multyears[,,yearInd]
  effluxFrac_bottom.thisyear <- 1-as.numeric(unlist(fseq_bottom.thisyear))
  return(effluxFrac_bottom.thisyear)
}

constrainFlux <- function(dataIndex, effluxFracs){
  adjFlux <- Sala_CO2_efflux.df$co2_efflux[dataIndex]*effluxFracs[coord.matches$Closest[dataIndex]]
  return(adjFlux)
}

# run the calculations

# *** these still blew up memory running with the entire dataset; even worse
# when using mclapply ... but much better with a reduced dataset

# set up structure to hold results

predicted.PgCO2_per_year_to_atmos_EEZs <- as.data.frame(matrix(data = NA, 
                                                          nrow = length(seqFracYears.raw),
                                                          ncol = 2))
colnames(predicted.PgCO2_per_year_to_atmos_EEZs) = c("Year",
                                                "PgCO2_per_year_to_atmos_China")
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
