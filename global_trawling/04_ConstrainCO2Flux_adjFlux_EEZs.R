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

