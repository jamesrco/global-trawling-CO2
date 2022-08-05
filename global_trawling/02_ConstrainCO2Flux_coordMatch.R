# 02_ConstrainCO2Flux_coordMatch.R
# Created June 7, 2022
# Purpose: Second in series of scripts used to constrain the crude estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This second script is set up to farm some intensive calculations out to AWS
# *** Assumes user has run the four shell provisioning scripts in
# zoonotic-C/aws-provisioning to get the remote computer provisioned properly

# set the working directory; create directory for output

setwd("~/global-trawling-CO2") # for AWS
dir.create("data/global_trawling/derived/output")

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

# library(ncdf4) # not needed anymore

# now that we have (more or less) apples to oranges, can find best match for
# each point in the CO2 flux dataset

# ******************************************************************************
# one approach: using just data.table
# ******************************************************************************

# *** seemed to take too long on my laptop and didn't lend itself to obvious parallelization
# but, given the error with approach #2, below, trying this with a 32 vCPU EC2
# instance on AWS (64 GB RAM plus a 35 GB disk swap)

# *** results from run on a c5a.8xlarge 32-vCPU AWS EC2 machine with 64 GB of
# RAM (and a 35 GB disk swap): Time difference of 2.307836 days; generated 
# a 1.6 GB data object

# first, set number of threads (need to change depending on cores or vCPUS)

setDTthreads(32) # shouldn't allow you to exceed actual # of cores or vCPUS,
# at least on linux

# convert to data tables; add decimal where necessary
Sala_CO2_efflux.coords.dt <- data.table(Sala_CO2_efflux.coords/10^5)
Siegel_fseq.coords.dt <- data.table(Siegel_fseq.coords.df)

# define function to find nearest match
# adapted from https://stackoverflow.com/questions/40211948/finding-closest-point-from-other-data-frame

dist1 <- function(a, b){
  dt <- data.table((Siegel_fseq.coords.dt$x-a)^2+(Siegel_fseq.coords.dt$y-b)^2)
  return(which.min(dt$V1))}

# find matches

# # test with a subset first; with benchmarking
# 
# Sala_CO2_efflux.coords.dt.sub <- Sala_CO2_efflux.coords.dt[1:10000,]
# 
# # find matches
# 
# time0 <- Sys.time()
# 
# coord.matches.test1 <- Sala_CO2_efflux.coords.dt.sub[, j = list(Closest =  dist1(x, y)), by = 1:nrow(Sala_CO2_efflux.coords.dt.sub)]
# 
# time1 <- Sys.time()
# print(time1 - time0)

# now with the whole enchilada
# *** with whole dataset, just ran and ran and ran, at least on my 2015 quad-core Intel i7 MBP

time0 <- Sys.time()

coord.matches <- Sala_CO2_efflux.coords.dt[, j = list(Closest =  dist1(x, y)), by = 1:nrow(Sala_CO2_efflux.coords.dt)]

time1 <- Sys.time()
print(time1 - time0)

# # ******************************************************************************
# # second approach: using data.table and apply to take advantage of parallelization
# # ******************************************************************************
# 
# # *** unfortunately, this approach didn't seem to work so well when I sent it to
# # a c5a.8xlarge 32-vCPU AWS EC2 machine with 64 GB of RAM (and a 35 GB disk swap):
# # ran for about half a day (benchmark time difference of 13.84486 hours),
# # which seemed right, but then returned this error:
# #
# # Warning message:                                                        
# # In parallel::mclapply(X = X, FUN = FUN, ...) :
# #  scheduled cores 1, 3, 5, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32 did not deliver results, all 
# #  values of the jobs will be affected
# #
# # the output object did save, with the correct number of rows, but I don't think
# # it worked correctly ... will try approach #1, above, on a similar EC2 instance
# 
# # convert Sala data to a data.frame, move decimal
# Sala_CO2_efflux.coords.df <- as.data.frame(Sala_CO2_efflux.coords/10^5)
# 
# # create a parallel version of sapply
# 
# mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
#   FUN <- match.fun(FUN)
#   answer <- parallel::mclapply(X = X, FUN = FUN, ...)
#   if (USE.NAMES && is.character(X) && is.null(names(answer))) 
#     names(answer) <- X
#   if (!isFALSE(simplify) && length(answer)) 
#     simplify2array(answer, higher = (simplify == "array"))
#   else answer
# }
# 
# # define function to find nearest match 
# # adapted from https://stackoverflow.com/questions/40211948/finding-closest-point-from-other-data-frame
# 
# dist2 <- function(df){
#   dt <- data.table((Siegel_fseq.coords.df$x-df$x)^2+(Siegel_fseq.coords.df$y-df$y)^2)
#   return(which.min(dt$V1))}
# 
# # test with a subset first; with benchmarking
# 
# Sala_CO2_efflux.coords.df.sub <- Sala_CO2_efflux.coords.df[1:10000,]
# 
# # find matches; may need to adjust no. of cores
# 
# time0 <- Sys.time()
# 
# coord.matches.test2 <- mcsapply(1:nrow(Sala_CO2_efflux.coords.df.sub), function(x) return(dist2(Sala_CO2_efflux.coords.df.sub[x,])), mc.cores = 4)
# 
# time1 <- Sys.time()
# print(time1 - time0)
# 
# # *** definitely faster than the approach #1 above, at least on my 2015 quad-core Intel i7 MBP
# 
# # save output
# 
# save(coord.matches.test2, file = "output/coord.matches.test2.RData")
# 
# # now with the whole enchilada
# 
# time0 <- Sys.time()
# 
# coord.matches <- mcsapply(1:nrow(Sala_CO2_efflux.coords.df), function(x) return(dist2(Sala_CO2_efflux.coords.df[x,])), mc.cores = 32)
# 
# time1 <- Sys.time()
# print(time1 - time0)

# save output

save(coord.matches, file = "data/global_trawling/derived/output/coord.matches.RData")

# send email when done ... assumes SSMTP has been installed and config file and text file for the email are in right place, etc.

system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-C/aws_provisioning/ssmtp/notification_email.txt"))

# continues onto 03b_ConstrainCO2Flux_adjFlux_extended.R, where we will perform
# the actual adjustment of the Sala et al. CO2 flux data