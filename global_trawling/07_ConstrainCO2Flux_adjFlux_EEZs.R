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
library(dplyr)
library(rgdal)
library(raster) # assumes you have some version of GDAL up and running, and 
# will force you to load several dependencies, including terra 
library(data.table)
library(parallel) # part of base; doesn't need to be installed
library(R.matlab) # to read .mat file
library(maptools)

# if not already loaded, load in the file containing the coordinate matches,
# generated in previous script 03_ConstrainCO2Flux_coordMatch.R and saved to 
# global-trawling-CO2/data/global_trawling/derived/output/

# also load some other necessary inputs, from earlier scripts

load("data/global_trawling/derived/output/coord.matches.NonZero.RData")
load("data/global_trawling/derived/benthic_seqfractions/fseq_bottom.multyears.RData")
seqFracYears.raw <- read.csv(file = "data/global_trawling/derived/benthic_seqfractions/trawlYears.csv",
                             header = FALSE)

# define functions

genEffluxFracs <- function(year){
  yearInd <- which(seqFracYears.raw==year)
  fseq_bottom.thisyear <- fseq_bottom.multyears[,yearInd]
  effluxFrac_bottom.thisyear <- 1-as.numeric(unlist(fseq_bottom.thisyear))
  return(effluxFrac_bottom.thisyear)
}

# load shapefile with EEZ polygons
FMI_EEZ_polygons_v11_raw <- st_read("data/global_trawling/raw/FMI_world_EEZ_boundaries_v11/eez_v11.shp")

# transform to match coordinate system we're using
FMI_EEZ_polygons_v11.54009 <- st_transform(FMI_EEZ_polygons_v11_raw, "ESRI:54009")
 
# define a list of countries we're interested in

# these are the countries with the 30 highest average landings from benthic trawling
# from 2009-2018 (per Sea Around Us data ... see calculations in 9_ConstrainCO2Flux_getSAUdata.R)
# and top 10 countries ranked by % of their total EEZ catch from bottom trawling gears 
# (per Steadman et al., 2021, New perspectives on an old fishing practice: Scale, 
# context and impacts of bottom trawling; available at
# https://oursharedseas.com/wp-content/uploads/2021/12/HI-RES-REPORT-‘New-perspectives-on-an-old-fishing-practice.pdf )

trawlEEZs <- c("China",
                   "Vietnam",
                   "Indonesia",
                   "India",
                   "Morocco",
                   "Japan",
                   "United States",
                   "South Korea",
                   "United Kingdom",
                   "Malaysia",
                   "Argentina",
                   "Myanmar",
                   "Mexico",
                   "Thailand",
                   "Russia",
                   "Norway",
                   "Guinea",
                   "Guinea-Bissau",
                   "Namibia",
                   "New Zealand",
                   "Denmark",
                   "Angola",
                   "Canada",
                   "Iceland",
                   "Pakistan",
                   "Turkey",
                   "Brazil",
                   "Bangladesh",
                   "South Africa",
                   "Ireland",
                   "Netherlands",
                   "Côte d'Ivoire",
                   "Germany",
                   "Republ. of Congo",
                   "Guyana")

# load in, format complete set of non-zero value points for compatibility
Sala_CO2_efflux.df.nonZeroCO2.raw <- read.csv("data/global_trawling/derived/output/Sala_CO2_efflux_nonZero.csv")
Sala_CO2_efflux.df.nonZeroCO2 <- data.frame(Sala_CO2_efflux.df.nonZeroCO2.raw)
rm(Sala_CO2_efflux.df.nonZeroCO2.raw)

points.sf <- st_as_sf(Sala_CO2_efflux.df.nonZeroCO2[,c("Sala_x","Sala_y")]*100000,
                      coords = c("Sala_x","Sala_y"),
                      crs = "ESRI:54009")

# make list object to hold EEZ points
EEZ.nonZeroCO2points <- vector(mode = "list", length = length(trawlEEZs))
names(EEZ.nonZeroCO2points) <- trawlEEZs

# iterate to subset points by EEZ

for (i in 1:length(EEZ.nonZeroCO2points)) {
  
  if (trawlEEZs[i]=="Côte d'Ivoire") {
    
    this.EEZ <- FMI_EEZ_polygons_v11.54009[FMI_EEZ_polygons_v11.54009$SOVEREIGN1=="Ivory Coast",]
    
    } else if (trawlEEZs[i]=="Republ. of Congo") {
      
      this.EEZ <- FMI_EEZ_polygons_v11.54009[FMI_EEZ_polygons_v11.54009$SOVEREIGN1=="Republic of the Congo",]
    
    } else {
      
      this.EEZ <- FMI_EEZ_polygons_v11.54009[FMI_EEZ_polygons_v11.54009$SOVEREIGN1==trawlEEZs[i],]
      
    }
  
      thisEEZ.CO2points <- as.data.frame(st_intersects(x = points.sf, y = this.EEZ))
      EEZ.nonZeroCO2points[[i]] <- thisEEZ.CO2points$row.id
  
}

# save this object
# worth noting it's pretty unlikely that there's *no* remineralization going on in
# Bangladeshi, Guyanese, Myanmar and Pakistani EEZs considering how much catch these
# countries land from benthic trawling (based on SAU data)
# maybe these countries are lacking in AIS ground stations?

save(EEZ.nonZeroCO2points, file = "data/global_trawling/derived/output/EEZ.nonZeroCO2points.RData")
# load("data/global_trawling/derived/output/EEZ.nonZeroCO2points.RData")

# now, can revisit the adjusted emissions calculations by EEZ

# first, year-by-year

# structure to hold results
predicted.PgCO2_per_year_to_atmos.byEEZ <- as.data.frame(matrix(data = NA, 
                                                          nrow = length(seqFracYears.raw),
                                                          ncol = 2 + length(trawlEEZs)))
colnames(predicted.PgCO2_per_year_to_atmos.byEEZ) = c("Year",
                                                "PgCO2_per_year_to_atmos_global_all_depths",
                                                trawlEEZs)
predicted.PgCO2_per_year_to_atmos.byEEZ[,1] <- unlist(seqFracYears.raw)

# iterate

for (i in 1:nrow(predicted.PgCO2_per_year_to_atmos.byEEZ)) {
  
  print(predicted.PgCO2_per_year_to_atmos.byEEZ[i,1])
  
  time0 <- Sys.time()
  
  thisYear <- predicted.PgCO2_per_year_to_atmos.byEEZ[i,1]
  
  EffluxFracs.thisyear <- genEffluxFracs(thisYear)
  adjCO2efflux.thisyear <- EffluxFracs.thisyear*Sala_CO2_efflux.df.nonZeroCO2$co2_efflux
  
  # store global calculation in column 2
  predicted.PgCO2_per_year_to_atmos.byEEZ[i,2] <- 
    sum(adjCO2efflux.thisyear, na.rm=T)*SalaModel_cell_area*(1/10^9)
  
  # now, do the calcs by EEZ
  
  for (j in 1:length(trawlEEZs)) {
    
    predicted.PgCO2_per_year_to_atmos.byEEZ[i,j+2] <- 
      sum(adjCO2efflux.thisyear[EEZ.nonZeroCO2points[[j]]], na.rm=T)*SalaModel_cell_area*(1/10^9)
    
  }
  
  time1 <- Sys.time()
  print(time1 - time0)
  
}

# save

write.csv(predicted.PgCO2_per_year_to_atmos.byEEZ, file = "data/global_trawling/derived/output/adjCO2efflux_global_PgCO2_yr_byEEZ.csv",
          row.names = FALSE)

# second, cumulative calculation

# structure to hold results
adjCO2efflux_PgCO2_cumulative.byEEZ <- vector(mode = "list", length = 2)
names(adjCO2efflux_PgCO2_cumulative.byEEZ) <- c("adjusted","unadjusted")

adjCO2efflux_PgCO2_cumulative.byEEZ[[1]] <- as.data.frame(matrix(data = NA, 
                                                      nrow = 200,
                                                      ncol = 2 + length(trawlEEZs)))
colnames(adjCO2efflux_PgCO2_cumulative.byEEZ[[1]]) = c("Year",
                                            "PgCO2_to_atmos_cumulative_global_alldepths",
                                            trawlEEZs)
adjCO2efflux_PgCO2_cumulative.byEEZ[[1]][,1] <- unlist(seqFracYears.raw)[1:200]

adjCO2efflux_PgCO2_cumulative.byEEZ[[2]] <- as.data.frame(matrix(data = NA, 
                                                                 nrow = 200,
                                                                 ncol = 2 + length(trawlEEZs)))
colnames(adjCO2efflux_PgCO2_cumulative.byEEZ[[2]]) = c("Year",
                                                       "PgCO2_to_atmos_cumulative_global_alldepths",
                                                       trawlEEZs)
adjCO2efflux_PgCO2_cumulative.byEEZ[[2]][,1] <- unlist(seqFracYears.raw)[1:200]

# iterate

for (i in 1:nrow(adjCO2efflux_PgCO2_cumulative.byEEZ[[1]])) {
  
  adjCO2efflux_PgCO2_cumulative.byEEZ[[1]][i,2] <- 
    sum(rev(predicted.PgCO2_per_year_to_atmos.byEEZ[1:i,2])*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  
  adjCO2efflux_PgCO2_cumulative.byEEZ[[2]][i,2] <- 
    sum(Sala_et_al_trawlTiming_results.raw$C_remin[1]*(1/10^9)*(44/12)*
          (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
             Sala_et_al_trawlTiming_results.raw$C_remin[1]))
  
  for (j in 1:length(trawlEEZs)) {
  
    adjCO2efflux_PgCO2_cumulative.byEEZ[[1]][i,j+2] <- 
      sum(rev(predicted.PgCO2_per_year_to_atmos.byEEZ[1:i,j+2])*
            (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
               Sala_et_al_trawlTiming_results.raw$C_remin[1]))
    
    adjCO2efflux_PgCO2_cumulative.byEEZ[[2]][i,j+2] <- 
      sum(sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[EEZ.nonZeroCO2points[[j]]], na.rm=T)*SalaModel_cell_area*(1/10^9)*
            (Sala_et_al_trawlTiming_results.raw$C_remin[1:i]/
               Sala_et_al_trawlTiming_results.raw$C_remin[1]))
    
  }
  
}

# save

write.csv(adjCO2efflux_PgCO2_cumulative.byEEZ, file = "data/global_trawling/derived/output/adjCO2efflux_global_PgCO2_cumulative_byEEZ.csv",
          row.names = FALSE)
# adjCO2efflux_PgCO2_cumulative.byEEZ <- read.csv(file = "data/global_trawling/derived/output/adjCO2efflux_global_PgCO2_cumulative_byEEZ.csv")

# quick look at % differences in the first year between adjusted and unadjusted
# estimates

adjCO2efflux_PgCO2_cumulative.byEEZ[[1]][1,]/
  adjCO2efflux_PgCO2_cumulative.byEEZ[[2]][1,]

# what explains this pattern? depth?

# as a first guess, let's take a look at weighted average depth of trawled area, by EEZ

wt_avg_EEZtrawldepths <- data.frame(vector(mode = "numeric", length = length(trawlEEZs)))
rownames(wt_avg_EEZtrawldepths) <- trawlEEZs
colnames(wt_avg_EEZtrawldepths) <- "weighted_avg_Depth_m"

for (i in 1:nrow(wt_avg_EEZtrawldepths)) {
  
  wt_avg_EEZtrawldepths[i,1] <- sum(Sala_CO2_efflux.df.nonZeroCO2$bottom_depth[EEZ.nonZeroCO2points[[i]]]*
                                   (Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[EEZ.nonZeroCO2points[[i]]]/
                                  sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[EEZ.nonZeroCO2points[[i]]])))
  
}

EEZtrawldepths.plotdata <- as.data.frame(cbind(as.numeric(unlist(wt_avg_EEZtrawldepths)),
                                as.numeric(adjCO2efflux_PgCO2_cumulative.byEEZ[1,3:37]/
                                adjCO2efflux_PgCO2_cumulative.byEEZ[1,40:74]),
                                as.numeric((as.numeric(adjCO2efflux_PgCO2_cumulative.byEEZ[30,3:37]-
                                   adjCO2efflux_PgCO2_cumulative.byEEZ[30,40:74])/
                                   adjCO2efflux_PgCO2_cumulative.byEEZ[30,40:74])*100)))
 # EEZtrawldepths.plotdata <- as.data.frame(cbind(as.numeric(unlist(wt_avg_EEZtrawldepths)),
 #                   as.numeric(unlist(adjCO2efflux_PgCO2_cumulative.byEEZ[[1]][1,]/
 #                     adjCO2efflux_PgCO2_cumulative.byEEZ[[2]][1,]))[3:37]))
colnames(EEZtrawldepths.plotdata) <- c("Weighted avg. depth by mass",
                        "Adjusted as % of unadjusted","Adjusted-Sala % deviation")
rownames(EEZtrawldepths.plotdata) <- trawlEEZs

EEZtrawldepths.plotdata <- EEZtrawldepths.plotdata[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar")),]

# make points of relative size
# one set, based on fraction of total est. global sediment remineralization
rel.sizeEEZpts.log <-
  -1/log10(adjCO2efflux_PgCO2_cumulative.byEEZ[1,3:37]/
             adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.PgCO2_to_atmos_cumulative_global_alldepths[1])
rel.sizeEEZpts.log <-
  rel.sizeEEZpts.log[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))]

rel.sizeEEZpts <-
  adjCO2efflux_PgCO2_cumulative.byEEZ[1,3:37]/
  adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.PgCO2_to_atmos_cumulative_global_alldepths[1]
rel.sizeEEZpts <-
  rel.sizeEEZpts[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))]

# another alternate set, based on avg landing tonnage from benthic trawling
# assumes user has run 09_ConstrainCO2Flux_getSAUdata.R

# make vector of landings data in right order
meanBenthicLandings_metrictons <- vector(mode = "numeric", length = length(trawlEEZs))
names(meanBenthicLandings_metrictons) <- trawlEEZs
for (i in 1:length(meanBenthicLandings_metrictons)) {
  
  meanBenthicLandings_metrictons[i] <- SAU_benthicCatchdata.sorted.aggregated$`Mean_2009-2018`[SAU_benthicCatchdata.sorted.aggregated$Country==names(meanBenthicLandings_metrictons)[i]]
  
}

rel.sizeEEZpts.landings <-
  meanBenthicLandings_metrictons/10^6
rel.sizeEEZpts.landings <-
  rel.sizeEEZpts.landings[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))]

# now, plot: with symbol size based on fraction of total est. global sediment remineralization 
plot(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
     pch = NA,
     xlim = c(0,500)) 

points(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
       EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
     pch = 21,
     col= "black",
     bg = "lightgrey",
     cex = as.numeric(rel.sizeEEZpts.log))

text(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`+20+rel.sizeEEZpts*40,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`+0.025,
     trawlEEZs[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))],
     cex = 0.7)

# set up some symbol sizes for the legend
leg.pchSizes <- c(.7,.5,.25,.1,.05)
leg.pchSizes.log <- -1/log10(leg.pchSizes)

legend(325, 0.95,
       legend = leg.pchSizes,
       title = paste("Fraction of total est.\nglobal sediment\nC remineralization"),
       pch = 21,
       pt.bg = "lightgrey",
       pt.cex = leg.pchSizes.log,
       box.lty=0,
       x.intersp = 3,
       y.intersp = c(3,2,1,1,1))

# pointLabel(-plotdata$`Weighted avg. depth by mass`,
#            plotdata$`Adjusted as % of unadjusted`,
#            labels = trawlEEZs[trawlEEZs!=c("Greenland","Philippines")],
#            method = "SANN",
#            offset = 0,
#            cex = 0.7,
#            allowSmallOverlap = FALSE,
#            trace = FALSE,
#            doPlot = TRUE)
           
# generally, depth appears to be a good explainer of the deviation 

# plot with symbol size based on avg landing tonnage from benthic trawling 
plot(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
     pch = NA,
     xlim = c(0,500)) 

points(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
       EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
       pch = 21,
       col= "black",
       bg = "lightgrey",
       cex = as.numeric(rel.sizeEEZpts.landings))

text(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`+20+rel.sizeEEZpts*40,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`+0.025,
     trawlEEZs[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))],
     cex = 0.7)

# set up some symbol sizes for the legend
leg.pchSizes <- c(5,2.5,1,0.5,0.25,.1)

legend(325, 0.95,
       legend = leg.pchSizes,
       title = paste("Avg. landings from\nbenthic trawling,\n2009-2018\n(Mt biomass)"),
       pch = 21,
       pt.bg = "lightgrey",
       pt.cex = leg.pchSizes,
       box.lty=0,
       x.intersp = 3,
       y.intersp = c(3,2,1,1,1))
