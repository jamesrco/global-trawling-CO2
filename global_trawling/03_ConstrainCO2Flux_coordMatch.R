# 03_ConstrainCO2Flux_coordMatch.R
# Created June 7, 2022
# Purpose: Third in series of scripts used to constrain the crude estimate of
# benthic CO2 flux in Sala  et al. 2021 using the sequestration fractions in
# Siegel et al. 2021
# Author: Jamie Collins, jcollins@edf.org

# This third script is set up to farm some calculations out to AWS, if needed/
# desired

# *** If using AWS, assumes user has run the four shell provisioning scripts in
# zoonotic-C/aws-provisioning to get the remote computer provisioned properly

# set the working directory; create directory for output

setwd("~/global-trawling-CO2") # for AWS
dir.create("data/global_trawling/derived/output_v2")

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
library(geodist) # calculate distance cheaply
library(pbapply) # parallelize distance calc
library(parallel) # parallelize distance calc

# now that we have (more or less) apples to oranges, can find best coordinate
# match for each point in the CO2 flux dataset

# # ******************************************************************************
# # one approach: using just data.table on the full (and very sparse) Sala data
# # ******************************************************************************
# 
# # *** seemed to take too long on my laptop and didn't lend itself to obvious parallelization
# # but, given the error with approach #2, below, trying this with a 32 vCPU EC2
# # instance on AWS (64 GB RAM plus a 35 GB disk swap)
# 
# # *** results from run on a c5a.8xlarge 32-vCPU AWS EC2 machine with 64 GB of
# # RAM (and a 35 GB disk swap): Time difference of 2.307836 days; generated 
# # a 1.6 GB data object
# 
# # first, set number of threads (need to change depending on cores or vCPUS)
# 
# setDTthreads(32) # shouldn't allow you to exceed actual # of cores or vCPUS,
# # at least on linux
# 
# # convert to data tables; add decimal where necessary
# Sala_CO2_efflux.coords.dt <- data.table(Sala_CO2_efflux.coords/10^5)
# 
# Siegel_fseq.coords.dt.oceanOnly <- data.table(Siegel_fseq.coords.df.oceanOnly)
# 
# # define function to find nearest lat/long match
# # adapted from https://stackoverflow.com/questions/40211948/finding-closest-point-from-other-data-frame
# 
# dist1 <- function(a, b){
#   dt <- data.table((Siegel_fseq.coords.dt.oceanOnly$x-a)^2+(Siegel_fseq.coords.dt.oceanOnly$y-b)^2)
#   return(which.min(dt$V1))}
# 
# # # find matches
# 
# # # test with a subset first; with benchmarking
# # 
# # Sala_CO2_efflux.coords.dt.sub <- Sala_CO2_efflux.coords.dt[1:10000,]
# # 
# # # find matches
# # 
# # time0 <- Sys.time()
# # 
# # coord.matches.test1 <- Sala_CO2_efflux.coords.dt.sub[, j = list(Closest =  dist1(x, y)), by = 1:nrow(Sala_CO2_efflux.coords.dt.sub)]
# 
# # time1 <- Sys.time()
# # print(time1 - time0)
# 
# # now with the whole enchilada
# # *** with whole dataset, just ran and ran and ran, at least on my 2015 quad-core Intel i7 MBP
# 
# time0 <- Sys.time()
# 
# coord.matches <- Sala_CO2_efflux.coords.dt[, j = list(Closest =  dist1(x, y)), by = 1:nrow(Sala_CO2_efflux.coords.dt)]
# 
# time1 <- Sys.time()
# print(time1 - time0)

# ******************************************************************************
# alternate first approach: using just data.table on a subset of the Sala data
# ******************************************************************************

# *** the coordinate matching took considerably less time (only about a half
# hour on my new MacBook Air M2) when I used a subset of the Sala et al
# dataset that contained only the coordinates of non-zero data points ... 
# which I should have done from the very start since the entire data set is 
# just so sparse to begin with

# # first, set number of threads (need to change depending on cores or vCPUS)
# 
# setDTthreads(32) # shouldn't allow you to exceed actual # of cores or vCPUS,
# # at least on linux

# convert to data tables; add decimal where necessary
Sala_CO2_efflux.coords.nonZero.dt <- data.table(Sala_CO2_efflux.coords.nonZero/10^5)

Siegel_fseq.coords.dt.oceanOnly <- data.table(Siegel_fseq.coords.df.oceanOnly)

# define function to find nearest lat/long match
# adapted from https://stackoverflow.com/questions/40211948/finding-closest-point-from-other-data-frame

dist1 <- function(a, b){
  dt <- data.table((Siegel_fseq.coords.dt.oceanOnly$x-a)^2+(Siegel_fseq.coords.dt.oceanOnly$y-b)^2)
  bestmatch.ind <- which.min(dt$V1)
  ind <- Siegel_fseq.coords.dt.oceanOnly$ind[bestmatch.ind]
  x <- Siegel_fseq.coords.dt.oceanOnly$x[bestmatch.ind]
  y <- Siegel_fseq.coords.dt.oceanOnly$y[bestmatch.ind]
  return(list(ind, x, y))
}
  
# find matches

# # test with a subset first; with benchmarking
# 
# Sala_CO2_efflux.coords.nonZero.dt.sub <- Sala_CO2_efflux.coords.nonZero.dt[1:10000,]
# 
# # find matches
# 
# time0 <- Sys.time()
# 
# coord.matchesNonZero.test1 <- Sala_CO2_efflux.coords.nonZero.dt.sub[, j = dist1(x, y), by = 1:nrow(Sala_CO2_efflux.coords.nonZero.dt.sub)]
# 
# time1 <- Sys.time()
# print(time1 - time0)

# now with the whole enchilada
# *** with whole dataset, just ran and ran and ran, at least on my 2015 quad-core Intel i7 MBP

time0 <- Sys.time()

coord.matchesNonZero <- Sala_CO2_efflux.coords.nonZero.dt[, j = dist1(x, y), by = 1:nrow(Sala_CO2_efflux.coords.nonZero.dt)]

time1 <- Sys.time()
print(time1 - time0)

# # ******************************************************************************
# # mrc approach: geodist and pbapply
# # ******************************************************************************

#relevant dataframes
# Siegel_fseq.coords.dt.oceanOnly # x, y, ind, 10442 x 3
# Sala_CO2_efflux.coords.nonZero.dt # x, y, 6864846 x 2

salan=data.matrix(Sala_CO2_efflux.coords.nonZero.dt)

#each row as list for parallel ops, remove cruise, date
salanl=lapply(seq_len(nrow(salan)), function(i) salan[i,])

# #get all combinations of cytograms
# nnn<-t(combn(1:nrow(dist1hbn),2))
# print('start wass')
# print(date())

#parlapply with progress bar
# https://cran.r-project.org/web/packages/pbapply/pbapply.pdf
#100 with progress bar, run 2
#   user  system elapsed
# 1.241   0.107  61.251
# nnnl=lapply(seq_len(nrow(nnn)), function(i) nnn[i,])

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(geodist))
parallel::clusterExport(cl, "Siegel_fseq.coords.dt.oceanOnly")

pbo = pboptions(type="txt")
siegel_ind <-pblapply(salanl, function(i) {
  which.min(geodist(i,Siegel_fseq.coords.dt.oceanOnly, measure = "haversine"))
}, cl = cl)

stopCluster(cl)

print('end distance')
print(date())

coord.matchesNonZero <- Siegel_fseq.coords.dt.oceanOnly[unlist(siegel_ind),]
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

# assign names, save output

colnames(coord.matchesNonZero) <- c("nrow","ind","x","y")
save(coord.matchesNonZero, file = "data/global_trawling/derived/output_v2/coord.matches.NonZero.RData")

# now, we can do some depth matching; ultimate goal: find nearest modeled depth in
# the OCIM model (from Siegel et al) to the ocean bottom depth at the locations of
# the data points in the Sala et al. dataset

# first we'll need to get the right bottom depth, in the right order, for each
# point in the Sala et al efflux data set

# create data table from coordinates of the depth data (these are not in the right order right now)
nonZeroValueDepths.coords.dt <- data.table(st_coordinates(Sala_CO2efflux_nonZeroValueDepths_adj)/10^5)

# create some character vectors (basically factors) that can be easily indexed and
# referenced for matching purposes
Sala_CO2_efflux.coords.nonZero.dt.rounded <- round(Sala_CO2_efflux.coords.nonZero.dt,6)
Sala_CO2_efflux.coords.nonZero.dt.rounded.char <- apply(Sala_CO2_efflux.coords.nonZero.dt.rounded, c(1,2), as.character)
Sala_CO2_efflux.coords.nonZero.dt.rounded.char <- apply(Sala_CO2_efflux.coords.nonZero.dt.rounded.char, 1, paste, collapse = "_")

nonZeroValueDepths.coords.dt.rounded <- round(nonZeroValueDepths.coords.dt,6)
nonZeroValueDepths.coords.dt.rounded.char <- apply(nonZeroValueDepths.coords.dt.rounded, c(1,2), as.character)
nonZeroValueDepths.coords.dt.rounded.char <- apply(nonZeroValueDepths.coords.dt.rounded.char, 1, paste, collapse = "_")

# reorder the bottom depth data to match the order of the efflux dataset
nonZeroValueDepths.ordered <- Sala_CO2efflux_nonZeroValueDepths_adj[order(match(nonZeroValueDepths.coords.dt.rounded.char, Sala_CO2_efflux.coords.nonZero.dt.rounded.char)),]   

# add the depths (and lats and longs, and the index to the OCIM48 lat/long 
# while we're at it too) back into the Sala_CO2_efflux.df data frame

Sala_CO2_efflux.df$bottom_depth <- rep(NA, nrow(Sala_CO2_efflux.df))
Sala_CO2_efflux.df$Sala_x <- rep(NA, nrow(Sala_CO2_efflux.df))
Sala_CO2_efflux.df$Sala_y <- rep(NA, nrow(Sala_CO2_efflux.df))
Sala_CO2_efflux.df$Siegel_ind <- rep(NA, nrow(Sala_CO2_efflux.df))
Sala_CO2_efflux.df$Siegel_x <- rep(NA, nrow(Sala_CO2_efflux.df))
Sala_CO2_efflux.df$Siegel_y <- rep(NA, nrow(Sala_CO2_efflux.df))

Sala_CO2_efflux.df$bottom_depth[ind.nonZeroCO2] <- nonZeroValueDepths.ordered$rvalue_1
Sala_CO2_efflux.df$Sala_x[ind.nonZeroCO2] <- Sala_CO2_efflux.coords.nonZero.dt$x
Sala_CO2_efflux.df$Sala_y[ind.nonZeroCO2] <- Sala_CO2_efflux.coords.nonZero.dt$y
Sala_CO2_efflux.df$Siegel_ind[ind.nonZeroCO2] <- coord.matchesNonZero$ind
Sala_CO2_efflux.df$Siegel_x[ind.nonZeroCO2] <- coord.matchesNonZero$x
Sala_CO2_efflux.df$Siegel_y[ind.nonZeroCO2] <- coord.matchesNonZero$y

Sala_CO2_efflux.df.nonZeroCO2 <- Sala_CO2_efflux.df[ind.nonZeroCO2,]

# save this object so we can reimport it in MATLAB; will just save a matrix
# containing values for the non-zero data points
write.csv(Sala_CO2_efflux.df.nonZeroCO2, "data/global_trawling/derived/output/Sala_CO2_efflux_nonZero.csv",
          row.names = FALSE)

system(paste0("ssmtp -v jcollins2139@gmail.com < ~/zoonotic-C/aws_provisioning/ssmtp/notification_email.txt"))

# continues onto 04_ConstrainCO2Flux_genTrawlSeqFracs.m, where we will perform
# depth matching and extract the best-fit sequestration fractions from the 
# Siegel mode output