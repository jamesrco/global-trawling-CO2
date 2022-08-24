# 08_ConstrainCO2Flux_plots_and_tables.R
# Created August 24, 2022
# Purpose: Eighth in series of scripts used to constrain the estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This eighth script generates figures, statistics and tables for publication

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
library(dplyr)
library(rgdal)
library(raster) # assumes you have some version of GDAL up and running, and 
# will force you to load several dependencies, including terra 
library(data.table)
library(parallel) # part of base; doesn't need to be installed
library(R.matlab) # to read .mat file

# distribution by depth of mass CO2 remineralized from benthos under Sala et
# al model
# from 02_ConstrainCO2Flux_IO.R

# what about a similar plot that sums things by mass instead

pdf("img_output/Fig1.pdf",
    width = 3.54, height = 3.54, # 90 mm = 3.54 in
    bg = "white", pointsize = 9,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(lwd = lwd.thisplot)

plot(effluxEstHist_mass,
     xlim = c(-2500,0),
     ylim = c(0,1),
     main = NULL,
     xlab = "Depth bin (m)",
     ylab = expression(paste("Est. efflux from sediment (Pg CO"["2"],")")),
     yaxt = "n",
     lwd = lwd.thisplot)

axis(2, at = seq(0,1,0.2), lwd = lwd.thisplot)

dev.off() 

sum(sums_PgCO2[34])/sum(sums_PgCO2, na.rm = T) # fraction CO2 efflux (by mass) in 0-100 m depth bin
sum(sums_PgCO2[33:34])/sum(sums_PgCO2, na.rm = T) # fraction CO2 efflux (by mass) in 0-200 m depth bin
sum(sums_PgCO2[1:32], na.rm = T)/sum(sums_PgCO2, na.rm = T) # fraction CO2 efflux (by mass) originating in water depths > 200 m

