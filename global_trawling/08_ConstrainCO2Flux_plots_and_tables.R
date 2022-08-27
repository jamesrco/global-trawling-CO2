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

# Fig. 1: distribution by depth of mass CO2 remineralized from benthos under Sala et al model
# from 02_ConstrainCO2Flux_IO.R

pdf("img_output/Fig1.pdf",
    width = 3.54, height = 2, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,4,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(effluxEstHist_mass,
     xlim = c(-2500,0),
     ylim = c(0,1),
     main = NULL,
     xlab = "Depth bin (m)",
     ylab = expression(paste("Est. efflux from sediment (Pg CO"["2"],")")),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     lwd = lwd.thisplot)

axis(2, at = seq(0,1,0.2), lwd = lwd.thisplot, tck = -0.02)
axis(1, lwd = lwd.thisplot, tck = -0.02)

dev.off() 

sum(sums_PgCO2[34])/sum(sums_PgCO2, na.rm = T) # fraction CO2 efflux (by mass) in 0-100 m depth bin
sum(sums_PgCO2[33:34])/sum(sums_PgCO2, na.rm = T) # fraction CO2 efflux (by mass) in 0-200 m depth bin
sum(sums_PgCO2[1:32], na.rm = T)/sum(sums_PgCO2, na.rm = T) # fraction CO2 efflux (by mass) originating in water depths > 200 m

# MLD validation data

# calculate how much mass is associated with emissions from waters in which

# (a) MLD >= bottom depth (no adjustment)

sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth))]*
      SalaModel_cell_area*(1/10^9))
sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth))])/
  sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux)

# 0.7981076 Pg CO2; fraction total mass: 0.5404034

# (b) MLD >= bottom depth ± 10 m

sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth-10))]*
      SalaModel_cell_area*(1/10^9))
sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth-10))])/
  sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux)

# 0.9225066 Pg CO2; fraction total mass: 0.6246347

# (c) MLD >= bottom depth ± 15 m

sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth-15))]*
      SalaModel_cell_area*(1/10^9))
sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth-15))])/
  sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux)

# 0.9774574 Pg CO2; fraction total mass: 0.6618422

# Fig. 2: plots of cumulative emissions over time, comparing adjusted and unadjusted data
# Time scale for these: 30 years
# 3 plots that will be panels a-c

pdf("img_output/Fig2a.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative$Year,
     adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_alldepths_unadjusted,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,30), ylim = c(0,20),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, all depths)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_alldepths,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

pdf("img_output/Fig2b.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative$Year,
     adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater_unadjusted,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,30), ylim = c(0,4.4),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, depths > 200 m)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

pdf("img_output/Fig2c.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative$Year,
     adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_400m_greater_unadjusted,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,30), ylim = c(0,2.25),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, depths > 400 m)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_400m_greater,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

pdf("img_output/Fig2d.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Year,
     adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.China,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,30), ylim = c(0,10.5),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, Chinese EEZ)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Year, adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.China,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

pdf("img_output/Fig2e.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Year,
     adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Morocco,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,30), ylim = c(0,0.11),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, Moroccan EEZ)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Year, adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Morocco,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

# Fig. S1: plots of cumulative emissions over time, comparing adjusted and unadjusted data
# Time scale for these: 200 years
# 3 plots that will be panels a-c

pdf("img_output/FigS1a.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative$Year,
     adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_alldepths_unadjusted,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,200), ylim = c(0,120),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, all depths)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_alldepths,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

pdf("img_output/FigS1b.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative$Year,
     adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater_unadjusted,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,200), ylim = c(0,27),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, depths > 200 m)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_200m_greater,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

pdf("img_output/FigS1c.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative$Year,
     adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_400m_greater_unadjusted,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,200), ylim = c(0,14),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, depths > 400 m)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative$Year, adjCO2efflux_PgCO2_cumulative$PgCO2_to_atmos_cumulative_global_400m_greater,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off()

pdf("img_output/FigS1d.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Year,
     adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.China,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,200), ylim = c(0,65),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, Chinese EEZ)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Year, adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.China,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

pdf("img_output/FigS1e.pdf",
    width = 3.54, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,1),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot)

plot(adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Year,
     adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Morocco,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,200), ylim = c(0,0.65),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, Moroccan EEZ)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Year, adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Morocco,
      col = "black", lwd = lwd.thisplot)

# add legend
legend("topleft",
       legend = c("Sala et al. estimate","Adjusted estimate"),
       col=c("darkgrey", "black"), lty=c(5,1), lwd = lwd.thisplot,
       box.lty=0,
       bg = NULL)

dev.off() 

# Fig. 3: relationship between adjusted-unadjusted deviation and some measure of depth of EEZ
# from 07_ConstrainCO2Flux_adjFlux_EEZs.R

# first version, with symbol size based on fraction of total est. global sediment remineralization 
pdf("img_output/Fig3.pdf",
    width = 4.75, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,10),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot,
    xpd=TRUE)

plot(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
     lwd = lwd.thisplot,
     pch = NA,
     xlim = c(0,550),
     ylim = c(0,1.1),
     yaxs = "i",
     xaxs = "i",
     yaxt = "n",
     xaxt = "n",
     xlab = "Weighted mean depth of EEZ trawled area (m)",
     ylab = paste("Adjusted emissions estimate as fraction\nof original Sala et al. estimate"))

# perform linear regression; plot
EEZtrawldepths.depthinv <- -EEZtrawldepths.plotdata[,1]
EEZ_fit <- lm(EEZtrawldepths.plotdata[,2] ~ EEZtrawldepths.depthinv)

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

points(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
       EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
       pch = 21,
       lwd = lwd.thisplot,
       col= "black",
       bg = "lightgrey",
       cex = as.numeric(rel.sizeEEZpts.log))

#add fitted regression line to scatterplot
EEZ_fit_predval <- predict(EEZ_fit, data.frame(EEZtrawldepths.depthinv))
lines(EEZ_fit_predval ~ EEZtrawldepths.depthinv, col = "light coral",
      lwd = 1,
      lty = 1)

text(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`+rel.sizeEEZpts*30-5,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`+0.025,
     trawlEEZs[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))],
     cex = 0.7,
     pos = 4)

r2 = summary(EEZ_fit)$adj.r.squared
my.p = summary(EEZ_fit)$coefficients[2,4]

rp = vector('expression',2)
rp[1] = substitute(expression(R^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=2)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, scientific = FALSE, digits = 1)))[2]

legend('topright', legend = rp, bty = 'n',
     cex = 0.8)
     
# set up some symbol sizes for the legend
# will need to be moved in Illustrator later
leg.pchSizes <- c(.7,.5,.25,.1,.05)
leg.pchSizes.log <- -1/log10(leg.pchSizes)

legend(600, 0.95,
       legend = leg.pchSizes,
       title = paste("Fraction of total est.\nglobal sediment\nC remineralization"),
       pch = 21,
       pt.bg = "lightgrey",
       pt.cex = leg.pchSizes.log,
       box.lty=0,
       bg = NULL,
       x.intersp = 3,
       y.intersp = c(3,2,1,1,1))

dev.off()

# alternate version, with symbol size based on avg landing tonnage from benthic trawling 
pdf("img_output/Fig3_alt.pdf",
    width = 4.75, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,10),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot,
    xpd=TRUE)

plot(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
     lwd = lwd.thisplot,
     pch = NA,
     xlim = c(0,550),
     ylim = c(0,1.1),
     yaxs = "i",
     xaxs = "i",
     yaxt = "n",
     xaxt = "n",
     xlab = "Weighted mean depth of EEZ trawled area (m)",
     ylab = paste("Adjusted emissions estimate as fraction\nof original Sala et al. estimate"))

# perform linear regression; plot
EEZtrawldepths.depthinv <- -EEZtrawldepths.plotdata[,1]
EEZ_fit <- lm(EEZtrawldepths.plotdata[,2] ~ EEZtrawldepths.depthinv)

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

points(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
       EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`,
       pch = 21,
       lwd = lwd.thisplot,
       col= "black",
       bg = "lightgrey",
       cex = as.numeric(rel.sizeEEZpts.landings))

#add fitted regression line to scatterplot
EEZ_fit_predval <- predict(EEZ_fit, data.frame(EEZtrawldepths.depthinv))
lines(EEZ_fit_predval ~ EEZtrawldepths.depthinv, col = "light coral",
      lwd = 1,
      lty = 1)

text(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`+rel.sizeEEZpts*30-5,
     EEZtrawldepths.plotdata$`Adjusted as % of unadjusted`+0.025,
     trawlEEZs[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))],
     cex = 0.7,
     pos = 4)

r2 = summary(EEZ_fit)$adj.r.squared
my.p = summary(EEZ_fit)$coefficients[2,4]

rp = vector('expression',2)
rp[1] = substitute(expression(R^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=2)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, scientific = FALSE, digits = 1)))[2]

legend('topright', legend = rp, bty = 'n',
       cex = 0.8)

# set up some symbol sizes for the legend
# will need to be moved in Illustrator later
leg.pchSizes <- c(5,2.5,1,0.5,0.25,.1)

legend(600, 0.95,
       legend = leg.pchSizes,
       title = paste("Avg. landings from\nbenthic trawling,\n2009-2018\n(Mt biomass)"),
       pch = 21,
       pt.bg = "lightgrey",
       pt.cex = leg.pchSizes,
       box.lty=0,
       bg = NULL,
       x.intersp = 3,
       y.intersp = c(3,2,1,1,1))

dev.off() 

# second alternate version, with symbol size based on avg landing tonnage from benthic trawling,
# and the y-axis using a different measure of deviation
pdf("img_output/Fig3_alt2.pdf",
    width = 4.75, height = 2.5, # 90 mm = 3.54 in
    bg = "white", pointsize = 7,
    colormodel = "cmyk",
    paper = "A4")

lwd.thisplot = 0.75

par(mar=c(3,5,1,10),
    mgp=c(2,0.5,0),
    lwd = lwd.thisplot,
    xpd=TRUE)

plot(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
     EEZtrawldepths.plotdata$`Adjusted-Sala % deviation`,
     lwd = lwd.thisplot,
     pch = NA,
     xlim = c(0,550),
     ylim = rev(c(5,-85)),
     yaxs = "i",
     xaxs = "i",
     yaxt = "n",
     xaxt = "n",
     xlab = "Weighted mean depth of EEZ trawled area (m)",
     ylab = paste("Adjusted estimate, % deviation from original\nSala et al. estimate\n(cumulative emissions after 30 years)"))

# perform linear regression; plot
EEZtrawldepths.depthinv <- -EEZtrawldepths.plotdata[,1]
EEZ_fit.alt <- lm(EEZtrawldepths.plotdata[,3] ~ EEZtrawldepths.depthinv)

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

points(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`,
       EEZtrawldepths.plotdata$`Adjusted-Sala % deviation`,
       pch = 21,
       lwd = lwd.thisplot,
       col= "black",
       bg = "lightgrey",
       cex = as.numeric(rel.sizeEEZpts.landings))

#add fitted regression line to scatterplot
EEZ_fit_predval.alt <- predict(EEZ_fit.alt, data.frame(EEZtrawldepths.depthinv))
lines(EEZ_fit_predval.alt ~ EEZtrawldepths.depthinv, col = "light coral",
      lwd = 1,
      lty = 1)

text(-EEZtrawldepths.plotdata$`Weighted avg. depth by mass`+rel.sizeEEZpts*30-5,
     EEZtrawldepths.plotdata$`Adjusted-Sala % deviation`+0.025,
     trawlEEZs[!(trawlEEZs %in% c("Guyana","Bangladesh","Pakistan","Myanmar"))],
     cex = 0.7,
     pos = 4)

r2.alt = summary(EEZ_fit.alt)$adj.r.squared
my.p.alt = summary(EEZ_fit.alt)$coefficients[2,4]

rp = vector('expression',2)
rp[1] = substitute(expression(R^2 == MYVALUE), 
                   list(MYVALUE = format(r2.alt,dig=2)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p.alt, scientific = FALSE, digits = 1)))[2]

legend('topright', legend = rp, bty = 'n',
       cex = 0.8)

# set up some symbol sizes for the legend
# will need to be moved in Illustrator later
leg.pchSizes <- c(5,2.5,1,0.5,0.25,.1)

legend(600, 0.95,
       legend = leg.pchSizes,
       title = paste("Avg. landings from\nbenthic trawling,\n2009-2018\n(Mt biomass)"),
       pch = 21,
       pt.bg = "lightgrey",
       pt.cex = leg.pchSizes,
       box.lty=0,
       bg = NULL,
       x.intersp = 3,
       y.intersp = c(3,2,1,1,1))

dev.off() 

# get some information concerning the regression

# summary(EEZ_fit)
# 
# > summary(EEZ_fit)
# 
# Call:
#   lm(formula = EEZtrawldepths.plotdata[, 2] ~ EEZtrawldepths.depthinv)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.49652 -0.08074  0.02328  0.13844  0.20553 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              0.8325444  0.0524483  15.874 7.73e-16 ***
#   EEZtrawldepths.depthinv -0.0017577  0.0002276  -7.722 1.63e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1805 on 29 degrees of freedom
# Multiple R-squared:  0.6728,	Adjusted R-squared:  0.6615 
# F-statistic: 59.62 on 1 and 29 DF,  p-value: 1.628e-08

summary(EEZ_fit.alt)

# > summary(EEZ_fit.alt)
# 
# Call:
#   lm(formula = EEZtrawldepths.plotdata[, 3] ~ EEZtrawldepths.depthinv)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -33.360  -3.532   0.888   5.550  19.954 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              1.46458    3.16352   0.463    0.647    
# EEZtrawldepths.depthinv -0.12887    0.01373  -9.386 2.72e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 10.89 on 29 degrees of freedom
# Multiple R-squared:  0.7523,	Adjusted R-squared:  0.7438 
# F-statistic: 88.09 on 1 and 29 DF,  p-value: 2.725e-10

# Finally, a Table S1 with some relevant data

# create table structure
TableS1.years <- c(1,5,10,25,30,50,75,100,200)
TableS1.numcols <- 3+3*2+length(TableS1.years)+3
TableS1.numrows <- 3+4+4+4+4+4+4+length(trawlEEZs)*4

TableS1 <- as.data.frame(matrix(data = NA,
                                nrow = TableS1.numrows,
                                ncol = TableS1.numcols))
                                     
TableS1[sort(c(seq(5,TableS1.numrows,4),seq(6,TableS1.numrows,4),
                         seq(7,TableS1.numrows,4))),2] <-
rep(c("Original estimate (Sala et al.)","Adjusted estimate (this analysis)",
                           "Percent deviation"),TableS1.numrows/4)

TableS1[seq(4,TableS1.numrows,4),1] <- c("Global trawled ocean area, all depths",
                                         "Global trawled ocean area, bottom depth║ <= 100 m",
                                         "Global trawled ocean area, bottom depth║ <= max. MLD¶",
                                         "Global trawled ocean area, bottom depth║ <= 200 m",
                                         "Global trawled ocean area, bottom depth║ <= 400 m",
                                         "Global trawled ocean area, bottom depth║ > 400 m",
                                         sort(trawlEEZs))

TableS1[1,c(3,5,7,9)] <- c("Avg. landings from benthic trawling, 2009-2018*",
                           "Weighted mean depth of area subjected to benthic trawling†",
                           "Fraction by mass of total est. global sediment C remineralization‡",
                           "Cumulative emissions after years of continuous trawling (Pg CO2)")
TableS1[2,c(3:8)] <- c("Mt biomass","Rank by EEZ§",
                       "m","Rank by EEZ§",
                       "—","Rank by EEZ§")
TableS1[2,c(9,11,12,13,14,16,17,18,19)] <- TableS1.years
TableS1[3,c(9,11,12,13,14,16,17,18,19)] <- "Pg CO2"
TableS1[3,c(10,15,20)] <- "Rank by EEZ§"

# populate with data

# landings data
TableS1[seq(28,TableS1.numrows,4),3] <- formatC(signif(meanBenthicLandings_metrictons/10^6,digits=2), digits=2,format="fg", flag="#")[sort(trawlEEZs, index.return = TRUE)$ix]
TableS1[seq(28,TableS1.numrows,4),4] <- order(meanBenthicLandings_metrictons/10^6, decreasing = T)[sort(trawlEEZs, index.return = TRUE)$ix]

# weighted mean trawled area depth
trawlDepthsforTable <- wt_avg_EEZtrawldepths$weighted_avg_Depth_m
trawlDepthsforTable[trawlDepthsforTable==0] <- NA
TableS1[seq(28,TableS1.numrows,4),5] <- round(-trawlDepthsforTable,digits=0)[sort(trawlEEZs, index.return = TRUE)$ix]
TableS1[seq(28,TableS1.numrows,4),6] <- rank(trawlDepthsforTable, na.last = T)[sort(trawlEEZs, index.return = TRUE)$ix]
TableS1[which(as.numeric(TableS1[,6])>length(trawlDepthsforTable)-sum(is.na(trawlDepthsforTable))),6] <- NA

# fraction total mass
TableS1[4,7] <- formatC(signif(sum(sums_PgCO2, na.rm = T)/sum(sums_PgCO2, na.rm = T),digits=2), digits=2,format="fg", flag="#")
TableS1[8,7] <- formatC(signif(sum(sums_PgCO2[34])/sum(sums_PgCO2, na.rm = T),digits=2), digits=2,format="fg", flag="#") # fraction CO2 efflux (by mass) in 0-100 m depth bin
TableS1[12,7] <- formatC(signif(sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux[which(MLDmatches.nonZero$annualmaxMLD_m>=(-Sala_CO2_efflux.df.nonZeroCO2$bottom_depth-10))])/
                                  sum(Sala_CO2_efflux.df.nonZeroCO2$co2_efflux),digits=2), digits=2,format="fg", flag="#") # fraction CO2 efflux (by mass) where bottom depth <= max. MLD ± 10 m
TableS1[16,7] <- formatC(signif(sum(sums_PgCO2[33:34])/sum(sums_PgCO2, na.rm = T),digits=2), digits=2,format="fg", flag="#") # fraction CO2 efflux (by mass) in 0-200 m depth bins
TableS1[20,7] <- formatC(signif(sum(sums_PgCO2[31:34])/sum(sums_PgCO2, na.rm = T),digits=2), digits=2,format="fg", flag="#") # fraction CO2 efflux (by mass) in 0-400 m depth bins
TableS1[24,7] <- formatC(signif(sum(sums_PgCO2[1:30], na.rm = T)/sum(sums_PgCO2, na.rm = T),digits=2), digits=2,format="fg", flag="#") # fraction CO2 efflux (by mass) originating in water depths > 400 m
fractMassbyEEZforTable <- 
  (adjCO2efflux_PgCO2_cumulative.byEEZ[1,40:74]/
  adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.PgCO2_to_atmos_cumulative_global_alldepths[1])[sort(trawlEEZs, index.return = TRUE)$ix]
fractMassbyEEZforTable[fractMassbyEEZforTable==0] <- NA
TableS1[seq(28,TableS1.numrows,4),7] <- formatC(signif(as.numeric(fractMassbyEEZforTable),digits=2), digits=2,format="fg", flag="#")
TableS1[which(as.numeric(TableS1[,7])>length(fractMassbyEEZforTable)-sum(is.na(fractMassbyEEZforTable))),7] <- NA
TableS1[seq(28,TableS1.numrows,4),8] <- as.integer(rank(-fractMassbyEEZforTable, na.last = T))
TableS1[which(as.numeric(TableS1[,8])>length(fractMassbyEEZforTable)-sum(is.na(fractMassbyEEZforTable))),8] <- NA

# cumulative emissions estimates

# year data, without ranks first

S1yearInsert.ind <- c(9,11,12,13,14,16,17,18,19)
  
for (i in 1:length(TableS1.years)) {
  
  thisYear <- TableS1.years[i]
  thisCol <- S1yearInsert.ind[i]
  
  TableS1[c(5,seq(29,TableS1.numrows,4)),thisCol] <-
    c(formatC(signif(as.numeric(adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,39]),digits=3), digits=3,format="fg", flag="#"),
      formatC(signif(as.numeric(adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,40:74]),digits=3), digits=3,format="fg", flag="#")[sort(trawlEEZs, index.return = TRUE)$ix])
  TableS1[which(as.numeric(TableS1[,thisCol])==0),thisCol] <- NA
  
  TableS1[c(6,seq(30,TableS1.numrows,4)),thisCol] <-
    c(formatC(signif(as.numeric(adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,2]),digits=3), digits=3,format="fg", flag="#"),
      formatC(signif(as.numeric(adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,3:37]),digits=3), digits=3,format="fg", flag="#")[sort(trawlEEZs, index.return = TRUE)$ix])
  TableS1[which(as.numeric(TableS1[,thisCol])==0),thisCol] <- NA
  
  TableS1[7,thisCol] <-
    formatC(signif(as.numeric(((adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,2]-adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,39])/
                                 adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,39])*100)
                   ,digits=2), digits=2,format="fg", flag="#")
  TableS1[seq(31,TableS1.numrows,4),thisCol] <-
    formatC(signif(as.numeric(((adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,3:37]-adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,40:74])/
                                 adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,40:74])*100)
                   ,digits=2), digits=2,format="fg", flag="#")[sort(trawlEEZs, index.return = TRUE)$ix]
  
}

# now do the rank data

S1rankInsert.ind <- c(10,15,20)
rankYears <- c(1,30,200)

for (i in 1:length(rankYears)) {
  
  thisYear <- rankYears[i]
  thisCol <- S1rankInsert.ind[i]
  
  TableS1[seq(29,TableS1.numrows,4),thisCol] <-
    as.integer(rank(-adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,40:74], na.last = T))[sort(trawlEEZs, index.return = TRUE)$ix]
  TableS1[which(as.numeric(TableS1[,thisCol])>length(fractMassbyEEZforTable)-sum(is.na(fractMassbyEEZforTable))),thisCol] <- NA
  
  TableS1[seq(30,TableS1.numrows,4),thisCol] <-
    as.integer(rank(-adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,3:37], na.last = T))[sort(trawlEEZs, index.return = TRUE)$ix]
  TableS1[which(as.numeric(TableS1[,thisCol])>length(fractMassbyEEZforTable)-sum(is.na(fractMassbyEEZforTable))),thisCol] <- NA
  
  TableS1[seq(31,TableS1.numrows,4),thisCol] <-
    rank(as.numeric(((adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,3:37]-adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,40:74])/
                       adjCO2efflux_PgCO2_cumulative.byEEZ[thisYear,40:74])*100), na.last = TRUE)[sort(trawlEEZs, index.return = TRUE)$ix]
  TableS1[which(as.numeric(TableS1[,thisCol])>length(fractMassbyEEZforTable)-sum(is.na(fractMassbyEEZforTable))),thisCol] <- NA
  
}

# now, save the table

write.csv(TableS1, file = "data/global_trawling/derived/output/TableS1_raw.csv",
          row.names = FALSE, col.names = FALSE, na = "")

