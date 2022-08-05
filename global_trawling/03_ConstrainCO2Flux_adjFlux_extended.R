# 03_ConstrainCO2Flux_adjFlux_extended.R
# Created July 12, 2022
# Purpose: Third in series of scripts used to constrain the estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This third script performs the actual adjustment of the Sala et al
# benthic CO2 flux data and includes time-integrated estimates of total
# global emissions to the atmosphere for a variety of time horizons

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

# we are ready at this point to adjust the fluxes in the Sala et al dataset
# using the appropriate (nearest match) benthic sequestration fractions in the
# Siegel et al dataset

# if not already loaded, load in the file containing the coordinate matches,
# generated in previous script 02_ConstrainCO2Flux_coordMatch.R and saved to 
# global-trawling-CO2/data/global_trawling/derived/output/

load("data/global_trawling/derived/output/coord.matches.RData")

# define functions

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

predicted.PgCO2_per_year_to_atmos <- as.data.frame(matrix(data = NA, 
                                                          nrow = length(seqFracYears.raw),
                                                          ncol = 2))
colnames(predicted.PgCO2_per_year_to_atmos) = c("Year",
                                                "PgCO2_per_year_to_atmos_global")
predicted.PgCO2_per_year_to_atmos[,1] <- unlist(seqFracYears.raw)

# iterate

for (i in 1:nrow(predicted.PgCO2_per_year_to_atmos)) {

  print(predicted.PgCO2_per_year_to_atmos[i,1])
  
  time0 <- Sys.time()
  
  thisYear <- predicted.PgCO2_per_year_to_atmos[i,1]
  EffluxFracs.thisyear <- genEffluxFracs(thisYear)
  adjCO2efflux.thisyear <- unlist(lapply(ind.nonZeroCO2, constrainFlux, EffluxFracs.thisyear))
  predicted.PgCO2_per_year_to_atmos[i,2] <- sum(adjCO2efflux.thisyear*SalaModel_cell_area, na.rm=T)*(1/10^9)
  
  time1 <- Sys.time()
  print(time1 - time0)
  
}

# save

write.csv(predicted.PgCO2_per_year_to_atmos, file = "data/global_trawling/derived/output/adjCO2efflux_global_PgCO2_yr.csv",
          row.names = FALSE)

# now we can make time-integrated estimates of total emissions to the atmosphere

# assumes:
# 1. annual CO2 flux from sediments follows the trend posited by Sala et al. 2021,
# contained in the object "results" that is generated beginning on line 119 in
# https://github.com/emlab-ucsb/ocean-conservation-priorities/blob/master/ancillary_analyses/timing_of_trawling_impacts.Rmd
# *** this is where the assertion in the Sala et al. paper that sediment emissions
# in the first year are 1.47 Pg CO2, declining after 10 y to 0.58 Pg CO2, comes from 
# 2. the fraction of the sediment CO2 emissions in year n that will have reached the
# atmosphere 100-n years later is described by the benthic emissions fraction
# for the nth year based on Siegel et al. 2021 (1 - fseq_bottom)

# automate; plot

adjCO2efflux_PgCO2_cumulative <- as.data.frame(matrix(data = NA, 
                                                          nrow = 200,
                                                          ncol = 2))
colnames(adjCO2efflux_PgCO2_cumulative) = c("Year",
                                            "PgCO2_to_atmos_cumulative_global")
adjCO2efflux_PgCO2_cumulative[,1] <- unlist(seqFracYears.raw)[1:200]

for (i in 1:nrow(adjCO2efflux_PgCO2_cumulative)) {
  
  adjCO2efflux_PgCO2_cumulative[i,2] <- 
    sum(rev(predicted.PgCO2_per_year_to_atmos[1:i,2])*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  
}

write.csv(adjCO2efflux_PgCO2_cumulative, file = "data/global_trawling/derived/output/adjCO2efflux_global_PgCO2_cumulative.csv",
          row.names = FALSE)

plot(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative,
     type = "l", col = "black", lwd = 1, xlab="Year", 
     ylab = expression(paste("Pg CO"["2"]," emitted to atmosphere from global benthic trawling activity (cumulative)")))

# after 100 y of continuous trawling, a cumulative 1.27 Pg CO2 will have reached the atmosphere 
# after 200 y of continuous trawling, a cumulative 3.95 Pg CO2 will have reached the atmosphere 

# # send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.
# 
# system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-c/aws_provisioning/ssmtp/notification_email.txt"))