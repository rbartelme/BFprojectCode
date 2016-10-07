### Ammonia Oxidation Calculation, Based on Mussman et al 2010 ##
#
# Version 3: 8-24-15, added wet weight biofilter sand conversion factor based on 
# 9-25-13 Wet:Dry Weight Experiment
#
omega_x <- 5 #Mean Cell Residence Time in days
#
omega <- 0.006611 #Hydraulic Retention time in Days
#
#Yao <- 0.34 #kg VSS/kg NH4+-N; assumes 20% loss of NH4 due to Heterotrophic assimilation; Mußman et al 2011
Yao <- 0.0525 #g organic carbon/g NH4+-N; Timmons and Ebeling, 2013
#
#
bao <-  0.15 #engodenous respiration rate constant per day****
# defined as oxygen consumption rate in the absence of substrate
#nearly indepdendent of substrate concentration
#governed by the growth and decay of heterotrophs and nitrifiers
#measured by BOD
DeltaAmmonia <- 5.7 # Delta Ammonia Across component in mg/L
cell.volume <- (4/3)*pi*(0.96/2)^3 
#0.96 um diameter from nichols & prosser paper for single C. Nitrosocosmicus (close relative based on amoA)

#cell volume in cubic micrometers
cell.mass <- 310*cell.volume*(10^-12)
#used Mussman et al conversion of 310 fg*C/um^3; 10^-12 term converts fg to mg
#
biomass <- (omega_x/omega)*(Yao/(1+bao*omega_x))*DeltaAmmonia #Mussman modified Rittman equation
biomass <- biomass/630 #630g sand in 1 liter of biofilter water
biomass.y <- biomass/cell.mass
#2 zeros removed from y lims; changed xlims to a day
plot(x=omega_x, y=biomass.y,xlim=c(0,30),ylim=c(1e6,1e10), log="y", col="red",xlab="MCRT in Days",  ylab="Biomass in cells/g sand")
title(main="Ammonia Oxidizer Biomass\nAs A Function of Mean Cell Residence Time")
for(i in 1:10000){
  x<-runif(1, 0.01, 30)
  DeltaAmmonia.V<-runif(1,2.67,8.73)
  y<-(x/omega)*(Yao/(1+bao*x))*DeltaAmmonia.V
  y<-y/630
  y<-y/cell.mass
  points(x=x,y=y,pch=42,col="Grey")
    }
for(i in 1:1000){
  x<-runif(1, 0.01, 30)
  y<-(x/omega)*(Yao/(1+bao*x))*DeltaAmmonia
  y<-y/630
  y<-y/cell.mass
  points(x=x,y=y,pch=42,col="Red")
    }
##############################################################
# Approximate Ammmonia Oxidizer Biomass Values based on qPCR #
##############################################################
#
#
abline(h=383162990, lty=2,col="Black")#Approximate Total Ammonia oxidizer biomass at surface based on qPCR estimates
text(10,1.1e8,font=2,"3.8e8 AOM cells via qPCR",cex=0.8, pos=3,col='black')
#
#Add legend for modeled biomass
legend(locator(1),legend=c("Mean Modeled AOM Biomass","Maximum/Minimum\nModeled AOM Biomass"),pch=c(15,15),col=c("Red","Grey"),cex=0.8,bty="n")
#
#Sum of surface archaea and AOB amoA genes 111849985
#Comammox amoA: 271313005
#Total approximate ammonia oxidizier biomass via qPCR = 383162990 cells/g



