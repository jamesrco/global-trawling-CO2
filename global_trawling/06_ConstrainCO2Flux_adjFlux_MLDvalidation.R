# 06_ConstrainCO2Flux_adjFlux_MLDvalidation.R
# Created August 19, 2022
# Purpose: Sixth in series of scripts used to constrain the estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This sixth script performs a simple validation of the primary approach to
# adjustment of the Sala fluxes using the Siegel et al. model output: Does 
# the adjusted estimate in year one agree with a simple estimate based on
# mixed layer depth

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

# load in MLD data
# using the Holte, Talley, et al., Argo-derived MLD climatology dataset, here:
# http://mixedlayer.ucsd.edu/data/Argo_mixedlayers_monthlyclim_04142022.mat

Argo_mixedlayers_monthlyclim_04142022.raw <- readMat("data/global_trawling/raw/UCSD_Argo_MLD_climatology/Argo_mixedlayers_monthlyclim_04142022.mat")

# extract annual max for each point; subset to only (ocean) locations that
# have max MLD estimates for them

annualMaxMLD <- apply(Argo_mixedlayers_monthlyclim_04142022.raw$mld.da.max, c(2,3), max, na.rm = T)
annualMaxMLD[which(is.infinite(annualMaxMLD))] <- NaN

Argo_MLD.max <- as.data.frame(cbind(as.vector(Argo_mixedlayers_monthlyclim_04142022.raw$latm),
                      as.vector(Argo_mixedlayers_monthlyclim_04142022.raw$lonm),
                      as.vector(annualMaxMLD)))
colnames(Argo_MLD.max) <- c("y","x","annualmaxMLD_m")

Argo_MLD.max.oceanOnly <- Argo_MLD.max[!is.na(Argo_MLD.max$annualmaxMLD_m),]
Argo_MLD.max.oceanOnly.dt <- data.table(Argo_MLD.max.oceanOnly)
                      
# specify a function to do the coordinate matching

distMLD <- function(a, b){
  dt <- data.table((Argo_MLD.max.oceanOnly.dt$x-a)^2+(Argo_MLD.max.oceanOnly.dt$y-b)^2)
  ind <- which.min(dt$V1)
  return(as.vector(Argo_MLD.max.oceanOnly.dt[ind,]))
}

# # do the matching, subset first
# 
# Sala_CO2_efflux.coords.nonZero.dt.sub <- Sala_CO2_efflux.coords.nonZero.dt[1:10000,]
# 
# time0 <- Sys.time()
# 
# MLDmatches.nonZero.sub <- Sala_CO2_efflux.coords.nonZero.dt.sub[, j = distMLD(x, y), by = 1:nrow(Sala_CO2_efflux.coords.nonZero.dt.sub)]
# 
# time1 <- Sys.time()
# print(time1 - time0)

# whole thing

time0 <- Sys.time()

MLDmatches.nonZero <- Sala_CO2_efflux.coords.nonZero.dt[, j = distMLD(x, y), by = 1:nrow(Sala_CO2_efflux.coords.nonZero.dt)]

time1 <- Sys.time()
print(time1 - time0)

# save this object

save(MLDmatches.nonZero, file = "data/global_trawling/derived/output/MLDmatches.nonZero.RData")

# now, can do some analysis of the data using bottom depth and MLD

# calculate how much mass is associated with emissions from waters in which

# (a) MLD >= bottom depth (no adjustment)

sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth))]*
      SalaModel_cell_area*(1/10^9))

# [1] 0.7981076 Pg CO2

# (b) MLD >= bottom depth ± 10 m

sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth-10))]*
      SalaModel_cell_area*(1/10^9))

# [1] 0.9225066 Pg CO2

# (c) MLD >= bottom depth ± 15 m

sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth-15))]*
      SalaModel_cell_area*(1/10^9))

# [1] 0.9774574 Pg CO2

# # send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.
# 
# system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-c/aws_provisioning/ssmtp/notification_email.txt"))