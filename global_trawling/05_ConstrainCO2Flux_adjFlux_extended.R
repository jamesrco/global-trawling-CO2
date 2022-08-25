# 05_ConstrainCO2Flux_adjFlux_extended.R
# Created July 12, 2022
# Purpose: Fifth in series of scripts used to constrain the estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This fifth script performs the actual adjustment of the Sala et al
# benthic CO2 flux data and includes time-integrated estimates of total
# global emissions to the atmosphere for a variety of time horizons

# *** Assumes user has already run the previous scripts in this series, and that
# the R objects generated using those scripts are in the environment already

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
# generated in previous script 03_ConstrainCO2Flux_coordMatch.R and saved to 
# global-trawling-CO2/data/global_trawling/derived/output/

load("data/global_trawling/derived/output/coord.matches.NonZero.RData")

# load the applicable sequestration fractions, generated in 04_ConstrainCO2Flux_genTrawlSeqFracs.m
# these are segmented so will need to be stitched together
# also load the years for which fractions were generated

fseq_bottom_multYears_1of7.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_1of7.mat")
fseq_bottom_multYears_2of7.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_2of7.mat")
fseq_bottom_multYears_3of7.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_3of7.mat")
fseq_bottom_multYears_4of7.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_4of7.mat")
fseq_bottom_multYears_5of7.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_5of7.mat")
fseq_bottom_multYears_6of7.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_6of7.mat")
fseq_bottom_multYears_7of7.raw <- readMat("data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_7of7.mat")

fseq_bottom_multYears.raw <- rbind(fseq_bottom_multYears_1of7.raw$fseq.bottom.multyears.1of7,
                                   fseq_bottom_multYears_2of7.raw$fseq.bottom.multyears.2of7,
                                   fseq_bottom_multYears_3of7.raw$fseq.bottom.multyears.3of7,
                                   fseq_bottom_multYears_4of7.raw$fseq.bottom.multyears.4of7,
                                   fseq_bottom_multYears_5of7.raw$fseq.bottom.multyears.5of7,
                                   fseq_bottom_multYears_6of7.raw$fseq.bottom.multyears.6of7,
                                   fseq_bottom_multYears_7of7.raw$fseq.bottom.multyears.7of7)

rm(fseq_bottom_multYears_1of7.raw,
   fseq_bottom_multYears_2of7.raw,
   fseq_bottom_multYears_3of7.raw,
   fseq_bottom_multYears_4of7.raw,
   fseq_bottom_multYears_5of7.raw,
   fseq_bottom_multYears_6of7.raw,
   fseq_bottom_multYears_7of7.raw)

fseq_bottom.multyears <- fseq_bottom_multYears.raw

fseq_bottom.multyears[fseq_bottom.multyears>=1] <- 1

# save a copy
save(fseq_bottom.multyears, file = "data/global_trawling/derived/benthic_seqfractions/fseq_bottom.multyears.RData")

seqFracYears.raw <- read.csv(file = "data/global_trawling/derived/benthic_seqfractions/trawlYears.csv",
                             header = FALSE)

# define functions

genEffluxFracs <- function(year){
  yearInd <- which(seqFracYears.raw==year)
  fseq_bottom.thisyear <- fseq_bottom.multyears[,yearInd]
  effluxFrac_bottom.thisyear <- 1-as.numeric(unlist(fseq_bottom.thisyear))
  return(effluxFrac_bottom.thisyear)
}

# run the calculations

# *** these still blew up memory running with the entire dataset; even worse
# when using mclapply ... but much better with a reduced dataset

# set up structure to hold results

predicted.PgCO2_per_year_to_atmos <- as.data.frame(matrix(data = NA, 
                                                          nrow = length(seqFracYears.raw),
                                                          ncol = 5))
colnames(predicted.PgCO2_per_year_to_atmos) = c("Year",
                                                "PgCO2_per_year_to_atmos_global_all_depths",
                                                "PgCO2_per_year_to_atmos_global_200m_less",
                                                "PgCO2_per_year_to_atmos_global_200m_greater",
                                                "PgCO2_per_year_to_atmos_global_400m_greater")
predicted.PgCO2_per_year_to_atmos[,1] <- unlist(seqFracYears.raw)

# iterate

for (i in 1:nrow(predicted.PgCO2_per_year_to_atmos)) {

#  print(predicted.PgCO2_per_year_to_atmos[i,1])
  
#  time0 <- Sys.time()
  
  thisYear <- predicted.PgCO2_per_year_to_atmos[i,1]
  EffluxFracs.thisyear <- genEffluxFracs(thisYear)
  adjCO2efflux.thisyear <- EffluxFracs.thisyear*Sala_CO2_efflux.df.nonZeroCO2$co2_efflux
  predicted.PgCO2_per_year_to_atmos[i,2] <- 
    sum(adjCO2efflux.thisyear*SalaModel_cell_area, na.rm=T)*(1/10^9)
  predicted.PgCO2_per_year_to_atmos[i,3] <- 
    sum(adjCO2efflux.thisyear[which(Sala_CO2_efflux.df.nonZeroCO2$bottom_depth > -200)]*SalaModel_cell_area, na.rm=T)*(1/10^9)
  predicted.PgCO2_per_year_to_atmos[i,4] <- 
    sum(adjCO2efflux.thisyear[which(Sala_CO2_efflux.df.nonZeroCO2$bottom_depth <= -200)]*SalaModel_cell_area, na.rm=T)*(1/10^9)
  predicted.PgCO2_per_year_to_atmos[i,5] <- 
    sum(adjCO2efflux.thisyear[which(Sala_CO2_efflux.df.nonZeroCO2$bottom_depth <= -400)]*SalaModel_cell_area, na.rm=T)*(1/10^9)
  
#  time1 <- Sys.time()
#  print(time1 - time0)
  
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
                                                          ncol = 8))
colnames(adjCO2efflux_PgCO2_cumulative) = c("Year",
                                            "PgCO2_to_atmos_cumulative_global_alldepths",
                                            "PgCO2_to_atmos_cumulative_global_200m_less",
                                            "PgCO2_to_atmos_cumulative_global_200m_greater",
                                            "PgCO2_to_atmos_cumulative_global_400m_greater",
                                            "PgCO2_to_atmos_cumulative_global_alldepths_unadjusted",
                                            "PgCO2_to_atmos_cumulative_global_200m_greater_unadjusted",
                                            "PgCO2_to_atmos_cumulative_global_400m_greater_unadjusted")
adjCO2efflux_PgCO2_cumulative[,1] <- unlist(seqFracYears.raw)[1:200]

for (i in 1:nrow(adjCO2efflux_PgCO2_cumulative)) {
  
  adjCO2efflux_PgCO2_cumulative[i,2] <-
    sum(rev(predicted.PgCO2_per_year_to_atmos[1:i,2])*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  adjCO2efflux_PgCO2_cumulative[i,3] <- 
    sum(rev(predicted.PgCO2_per_year_to_atmos[1:i,3])*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  adjCO2efflux_PgCO2_cumulative[i,4] <- 
    sum(rev(predicted.PgCO2_per_year_to_atmos[1:i,4])*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  adjCO2efflux_PgCO2_cumulative[i,5] <- 
    sum(rev(predicted.PgCO2_per_year_to_atmos[1:i,5])*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  adjCO2efflux_PgCO2_cumulative[i,6] <- 
    sum(Sala_et_al_trawlTiming_results.raw$C_remin[1]*(1/10^9)*(44/12)*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  adjCO2efflux_PgCO2_cumulative[i,7] <- 
    sum(sum(sums_PgCO2[1:32], na.rm = T)*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  adjCO2efflux_PgCO2_cumulative[i,8] <- 
    sum(sum(sums_PgCO2[1:30], na.rm = T)*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  }

write.csv(adjCO2efflux_PgCO2_cumulative, file = "data/global_trawling/derived/output/adjCO2efflux_global_PgCO2_cumulative.csv",
          row.names = FALSE)
# adjCO2efflux_PgCO2_cumulative <- read.csv("data/global_trawling/derived/output/adjCO2efflux_global_PgCO2_cumulative.csv")

# some exploratory plots

plot(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_alldepths,
     type = "l", col = "black", lwd = 1, xlab="Year", 
     ylab = expression(paste("Pg CO"["2"]," emitted to atmosphere from global benthic trawling activity (cumulative, all depths)")))

plot(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_less,
     type = "l", col = "black", lwd = 1, xlab="Year", 
     ylab = expression(paste("Pg CO"["2"]," emitted to atmosphere from global benthic trawling activity (cumulative, depths < 200 m)")))

plot(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater,
     type = "l", col = "black", lwd = 1, xlab="Year", xlim = c(0,30), ylim = c(0,4.4), 
     ylab = expression(paste("Pg CO"["2"]," emitted to atmosphere from global benthic trawling activity (cumulative, depths > 200 m)")))
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater_unadjusted,
      col = "red", lwd = 1)

plot(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater,
     type = "l", col = "black", lwd = 1, xlab="Year", xlim = c(0,200), ylim = c(0,27), 
     ylab = expression(paste("Pg CO"["2"]," emitted to atmosphere from global benthic trawling activity (cumulative, depths > 200 m)")))
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater_unadjusted,
      col = "red", lwd = 1)

plot(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_400m_greater,
     type = "l", col = "black", lwd = 1, xlab="Year", xlim = c(0,30), ylim = c(0,2.25), 
     ylab = expression(paste("Pg CO"["2"]," emitted to atmosphere from global benthic trawling activity (cumulative, depths > 400 m)")))
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_400m_greater_unadjusted,
      col = "red", lwd = 1)

plot(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_alldepths_unadjusted,
     type = "l", col = "red", lwd = 1, xlab="Year", xlim = c(0,30), ylim = c(0,20),
     ylab = expression(paste("Pg CO"["2"]," emitted to atmosphere from global benthic trawling activity (unadjusted)")))
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_alldepths,
      col = "black", lwd = 1)

# # send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.
# 
# system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-c/aws_provisioning/ssmtp/notification_email.txt"))