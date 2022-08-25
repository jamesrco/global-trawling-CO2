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
     adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Canada,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,30), ylim = c(0,0.15),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, Canadian EEZ)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Year, adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Canada,
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
     adjCO2efflux_PgCO2_cumulative.byEEZ$unadjusted.Canada,
     type = "l", col = "darkgrey", lty = 5, lwd = lwd.thisplot,
     xlab="Year", xlim = c(0,200), ylim = c(0,0.9),
     yaxt = "n",
     xaxt = "n",
     yaxs = "i",
     xaxs = "i",
     ylab = paste("Pg CO2 emitted to atmosphere\n(cumulative, Canadian EEZ)"))

axis(side = 1, lwd = lwd.thisplot, tck = -0.02)
axis(side = 2, lwd = lwd.thisplot, tck = -0.02)

# add adjusted estimate
lines(adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Year, adjCO2efflux_PgCO2_cumulative.byEEZ$adjusted.Canada,
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
     xlim = c(0,500),
     ylim = c(0,1.1),
     yaxs = "i",
     xaxs = "i",
     yaxt = "n",
     xaxt = "n",
     xlab = "Mass-weighted depth of EEZ trawled area (m)",
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
     trawlEEZs[trawlEEZs!=c("Greenland","Philippines")],
     cex = 0.7,
     pos = 4)

r2 = summary(EEZ_fit)$adj.r.squared
my.p = summary(EEZ_fit)$coefficients[2,4]

rp = vector('expression',2)
rp[1] = substitute(expression(R^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, scientific = FALSE, digits = 1)))[2]

legend('topright', legend = rp, bty = 'n',
     cex = 0.8)
     
# set up some symbol sizes for the legend
# will need to be moved in Illustrator later
leg.pchSizes <- c(.7,.5,.25,.1,.05)
leg.pchSizes.log <- -1/log10(leg.pchSizes)

legend(550, 0.95,
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

# get some information concerning the regression

summary(EEZ_fit)

# > summary(EEZ_fit)
# 
# Call:
#   lm(formula = EEZtrawldepths.plotdata[, 2] ~ EEZtrawldepths.depthinv)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.53470 -0.08247  0.05543  0.14212  0.23272 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              0.8833679  0.0639742  13.808 1.10e-11 ***
#   EEZtrawldepths.depthinv -0.0019380  0.0003466  -5.591 1.79e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1914 on 20 degrees of freedom
# Multiple R-squared:  0.6099,	Adjusted R-squared:  0.5904 
# F-statistic: 31.26 on 1 and 20 DF,  p-value: 1.793e-05

