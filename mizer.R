# I. Dynamic size spectrum model for the 2013 climate

rm(list=ls())
# Load packages and some additional functions

# # If run on Mac laptop
setwd("~/Google Drive/hybrid-SDM")

# If run on Windows Desktop
setwd("G:/My Drive/hybrid-SDM")

# Load necessary packages and functions
source("mizer_misc_functions.R")
library(mizerHowTo)
library(mizerExperimental)
library(tidyverse)

# Load speces parameters
nefsc_params<-read_csv("nefsc_species_params_scenario3.csv")

# Load interaction matrix (2013)
inter<-read_csv("interactions_2013_3.csv")

rownames(inter)<-colnames(inter)

# Load time-averaged catch data
avgCatch<-read_csv("catch2.csv")

# Set up parameters for mizer
params_uncalibrated<-newMultispeciesParams(nefsc_params,inter,kappa = 1e11)

species_params(params_uncalibrated)

# Specify fishing gear for each species.Target species information from NOAA website
# and mizer tutorial
gear_params(params_uncalibrated)$gear<-c("Pelagic","Pelagic","Industrial","Pelagic","Pelagic","Otter","Beam",
                                         "Beam","Otter","Otter","Otter","Otter","Beam","Beam","Otter","Otter","Otter","Otter","Industrial")

summary(params_uncalibrated)

# Run the model with defaultRmax values (Inf)
sim_uncalibrated<-project(object = params_uncalibrated,effort = 0.2,t_max = 100,t_save = 1)

# Normally, some species will go extinct so plotSummary function may not work here
plotlyBiomass(sim_uncalibrated)


# Roughly estimate R_max from carrying capacity and max weight of species
params_guessed<-params_uncalibrated
params_guessed@species_params$R_max <- params_guessed@resource_params$kappa*params_guessed@species_params$w_inf^-1

params_guessed<-setParams(params_guessed)

# Run the model with guessed Rmax
sim_guessed<-project(params_guessed,effort = 0.2,t_max = 100)
plotSummary(sim_guessed)
plotlyBiomass(sim_guessed)

RDI_RDD <- as.data.frame(getRDI(sim_guessed@params)/getRDD(sim_guessed@params))
RDI_RDD$w_inf<-params_uncalibrated@species_params$w_inf
RDI_RDD<-RDI_RDD[order(RDI_RDD$w_inf),]
colnames(RDI_RDD)[1] <- "ratio"
RDI_RDD

# Make a plot
ggplot(RDI_RDD) +
  geom_bar(aes(x = rownames(RDI_RDD), y = ratio, fill = rownames(RDI_RDD)),stat="identity") + scale_fill_manual(name = "Species", values = sim_guessed@params@linecolour) + scale_y_continuous(name = "RDI/RDD", trans = "log10") +
  scale_x_discrete(name = "Species") +
  theme(panel.background = element_rect(fill = "white", color = "black"), legend.position = "none")+
  geom_hline(aes(yintercept = 10))

# 1st optimization of Rmax to get the estimated biomass to be in the right ballpark
params_optim <- params_guessed
vary <-  log10(params_optim@species_params$R_max) # variable to explore
params_optim<-setParams(params_optim)
# set up workers
noCores <- detectCores() - 1 # keep some spare core
cl <- makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
  library(mizerExperimental)
  library(optimParallel)
})
optim_result <- optimParallel::optimParallel(par=vary,getError,params=params_optim, dat = avgCatch$catch, method   ="L-BFGS-B",lower=c(rep(3,19)),upper= c(rep(15,19)),
                                             parallel=list(loginfo=TRUE, forward=TRUE))
stopCluster(cl)
saveRDS(optim_result,"optimParallel_Rmax1.RDS")

# Update Rmax based on optimization results
params_optim@species_params$R_max <- 10^optim_result$par 
params_optim <-setParams(params_optim)
sim_optim <- project(params_optim, effort = 0.2, t_max = 100, dt=0.1,initial_n = sim_guessed@n[100,,],initial_n_pp = sim_guessed@n_pp[100,])
# saveRDS(sim_optim,"sim_optim1.RDS")
plotCalibration(sim_optim,catch_dat = avgCatch$catch)

plot_dat <- plotFmsy(params_optim, returnData = T)
# saveRDS(plot_dat, "Fmsy1.rds")

# Create another mizer parameter object for erepro-tweaking
params_optim2<-sim_optim@params

# Tweaking erepro based on Fmsy curves
## 1st round
params_optim2@species_params$erepro<-10^-1
params_optim2@species_params$erepro[-c(1,2)]<-10^-3
params_optim2@species_params$erepro[c(6,9,11:15,18)]<-10^-3.5
params_optim2@species_params$erepro[c(6,9,13,15)]<-10^-4
params_optim2@species_params$erepro[c(6,15)]<-10^-4.5

## 2nd round
params_optim2@species_params$erepro[c(1,2)]<-10^-0.7
params_optim2@species_params$erepro[2]<-10^-2
params_optim2@species_params$erepro[c(3,7,19)]<-10^-2.5
params_optim2@species_params$erepro[c(5,8)]<-10^-3.5
params_optim2@species_params$erepro[c(11,12,14)]<-10^-3.7
params_optim2@species_params$erepro[10]<-10^-4
params_optim2@species_params$erepro[13]<-10^-4.2
params_optim2@species_params$erepro[c(9,15)]<-10^-4.5

## 3rd round
params_optim2@species_params$erepro[1]<-10^-0.2
params_optim2@species_params$erepro[3]<-10^-2
params_optim2@species_params$erepro[4]<-10^-1.5
params_optim2@species_params$erepro[c(7,8)]<-10^-3
params_optim2@species_params$erepro[13]<-10^-4.5
params_optim2@species_params$erepro[19]<-10^-3

## 4th round
params_optim2@species_params$erepro[3]<-10^-2.3
params_optim2@species_params$erepro[4]<-10^-2
params_optim2@species_params$erepro[7]<-10^-2.5
params_optim2@species_params$erepro[8]<-10^-3.2
params_optim2@species_params$erepro[19]<-10^-2.7

## 5th round
params_optim2@species_params$erepro[7]<-10^-2.7
params_optim2@species_params$erepro[19]<-10^-3

plotFmsy(params_optim2)

# Run the model with tweaked erepro values
sim_optim2 <- project(params_optim2, effort = 0.2, t_max = 100, progress_bar = T)
plotSummary(sim_optim2)

# 2nd optimization: this time repeat 5 times to make sure we do get the optimal Rmax
sim_loop <- sim_optim2
params_loop <- sim_loop@params

for(i in 1:5)
{
  # saving the last time step in the param object
  params_loop@initial_n <- sim_loop@n[dim(sim_loop@n)[1],,]
  params_loop@initial_n_pp <- sim_loop@n_pp[dim(sim_loop@n_pp)[1],]
  params_calibration <- params_loop
  vary <-  log10(params_calibration@species_params$R_max)
  params_calibration<-setParams(params_calibration)
  noCores <- detectCores() - 1 # keep a spare core
  cl <- makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, as.list(ls()))
  clusterEvalQ(cl, {
    library(mizerExperimental)
    library(optimParallel)
  })
  optim_loop <-optimParallel::optimParallel(par=vary,getError,params=params_calibration, dat = avgCatch$catch, 
                                            method   ="L-BFGS-B",lower=c(rep(3,19)),upper= c(rep(15,19)),
                                            parallel=list(loginfo=TRUE, forward=TRUE))
  stopCluster(cl)
  
  # optim values:
  params_loop@species_params$R_max <- 10^optim_loop$par 
  # set the param object 
  params_loop <-setParams(params_loop)
  sim_loop <- project(params_loop, effort = 0.2, t_max = 100, dt=0.1, initial_n = params_loop@initial_n ,
                      initial_n_pp = params_loop@initial_n_pp, progress_bar = F)
}
saveRDS(optim_loop,"optimParallel_Rmax2.RDS")

sim_loop <- project(params_loop, effort = 0.2, t_max = 300, progress_bar = T)
plotSummary(sim_loop)

params_loop@initial_n <- sim_loop@n[dim(sim_loop@n)[1],,]
params_loop@initial_n_pp <- sim_loop@n_pp[dim(sim_loop@n_pp)[1],]
sim_loop <- project(params_loop, effort = 0.2, progress_bar = F)

plotSummary(sim_loop)
plotPredObsYield(sim_loop,dat = avgCatch$catch)
df<-plotPredObsYield(sim_loop,dat = avgCatch$catch,returnData = TRUE)
reg<-lm(df$obs~df$value)
summary(reg)

RDI_RDD <- as.data.frame(getRDI(sim_loop@params)/getRDD(sim_loop@params))
RDI_RDD$w_inf<-params_loop@species_params$w_inf
colnames(RDI_RDD)[1]<-"ratio"
plot(x = log10(RDI_RDD$w_inf),y = log10(RDI_RDD$ratio))

saveRDS(sim_loop,"sim_optim_final.RDS")

# Find the appropriate background resource level by examining growth curves
## Compare growth in the model vs growth according k_vb
params_growth<-params_loop
params_growth<-setParams(params_growth,kappa=1.5e11)

sim_growth<-project(params_growth, effort = 0.2, t_max = 100, 
                              initial_n = params_loop@initial_n ,
                              initial_n_pp = params_loop@initial_n_pp)

plotSummary(sim_growth)

checkGrowth(sim_growth)

## Check feeding level
plotFeedingLevel(sim_growth)
plotDiet2(sim_growth)

## Tweaking species gamma values to further calibrate model growth curves
shiny_gamma(sim_growth@params,dat = avgCatch$catch)

params_growth2<-setParams(params_growth)
params_growth2@species_params$gamma[15]<-5e-11

params_growth2<-setParams(params_growth2)


sim_growth2<-project(params_growth2, effort = 0.2, t_max = 100, 
                    initial_n = params_loop@initial_n ,
                    initial_n_pp = params_loop@initial_n_pp)

plotSummary(sim_growth2)
checkGrowth(sim_growth2)
plotFeedingLevel(sim_growth2)
plotDiet2(sim_growth2)
plotFmsy(sim_growth2@params)

plotPredObsYield(sim_growth2,dat = avgCatch$catch)
df<-plotPredObsYield(sim_loop,dat = avgCatch$catch,returnData = TRUE)
reg<-lm(df$obs~df$value)
summary(reg)

saveRDS(sim_growth2,"sim_2013.RDS")


# II. Dynamic size spectrum model for RCP45 and RCP85
rm(list=ls())
# Load packages and some additional functions

# If run on Mac laptop
setwd("~/Google Drive/hybrid-SDM")

# If run on Windows Desktop
setwd("G:/?ڪ????ݵw??/hybrid-SDM")

# Load necessary packages and functions
source("mizer_misc_functions.R")
library(mizerHowTo)
library(mizerExperimental)
library(tidyverse)

# Load the 2013 model
sim_2013<-readRDS("sim_2013.RDS")
avgCatch<-read_csv("catch2.csv")
inter_2013<-read_csv("interactions_2013_3.csv")
rownames(inter_2013)<-colnames(inter_2013)

# RCP45 model
# Read in species parameters under rcp45
nefsc_params_rcp45<-read_csv("nefsc_species_params_rcp45_scenario3.csv")

# Read in interaction matrix under rcp45
inter_rcp45<-read_csv("interactions_rcp45_3.csv")

rownames(inter_rcp45)<-colnames(inter_rcp45)

# Set up mizer parameters, using background resources in sim_Rmax_final
params_rcp45<-newMultispeciesParams(nefsc_params_rcp45,inter_rcp45,
                                    kappa = 1.12*sim_2013@params@resource_params$kappa)

params_rcp45@species_params$gamma<-sim_2013@params@species_params$gamma

params_rcp45<-setParams(params_rcp45)

# Specify fishing gear for each species.Target species information from NOAA website
# and mizer tutorial
gear_params(params_rcp45)$gear<-c("Pelagic","Pelagic","Industrial","Pelagic","Pelagic","Otter","Beam","Beam","Otter",
                                  "Otter","Otter","Otter","Beam","Beam","Otter","Otter","Otter","Otter","Industrial")

# Calculate R_max based on kappa and w_inf: RCP45
params_rcp45@species_params$R_max <- params_rcp45@resource_params$kappa*params_rcp45@species_params$w_inf^-1

params_rcp45@species_params$erepro<-sim_2013@params@species_params$erepro

# Run simulation, using the end species abundance in sim_Rmax_final as initial abundance 
sim_rcp45<- project(params_rcp45, effort = 0.2, t_max = 300,
                    initial_n = sim_2013@n[100,,])

plotSummary(sim_rcp45)

plotlyBiomass(sim_rcp45)

# Tweak erepro to get reasonable Fmsy curves
params_optim2_rcp45<-sim_rcp45@params

# 1st round 
params_optim2_rcp45@species_params$erepro[c(3,4)]<-10^-3
params_optim2_rcp45@species_params$erepro[c(9,15)]<-10^-4.7
params_optim2_rcp45@species_params$erepro[c(11,12,14,18)]<-10^-4
params_optim2_rcp45@species_params$erepro[c(17,19)]<-10^-3.5

# 2ng round
params_optim2_rcp45@species_params$erepro[14]<-10^-4.2
params_optim2_rcp45@species_params$erepro[15]<-10^-5
params_optim2_rcp45@species_params$erepro[18]<-10^-4.5

plotFmsy(params_optim2_rcp45)

sim_optim2_rcp45<-project(params_optim2_rcp45,effort = 0.2,t_max = 100)
plotSummary(sim_optim2_rcp45)

# Optimization of Rmax with estimated time-averaged catch data
sim_loop <- sim_optim2_rcp45
params_loop <- sim_loop@params

for(i in 1:5)
{
  # saving the last time step in the param object
  params_loop@initial_n <- sim_loop@n[dim(sim_loop@n)[1],,]
  params_loop@initial_n_pp <- sim_loop@n_pp[dim(sim_loop@n_pp)[1],]
  params_calibration <- params_loop
  vary <-  log10(params_calibration@species_params$R_max)
  params_calibration<-setParams(params_calibration)
  noCores <- detectCores() - 1 # keep a spare core
  cl <- makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, as.list(ls()))
  clusterEvalQ(cl, {
    library(mizerExperimental)
    library(optimParallel)
  })
  optim_loop <-optimParallel::optimParallel(par=vary,getError,params=params_calibration, dat = avgCatch$catch_inf_rcp45, 
                                            method   ="L-BFGS-B",lower=c(rep(3,19)),upper= c(rep(15,19)),
                                            parallel=list(loginfo=TRUE, forward=TRUE))
  stopCluster(cl)
  
  # optim values:
  params_loop@species_params$R_max <- 10^optim_loop$par 
  # set the param object 
  params_loop <-setParams(params_loop)
  sim_loop <- project(params_loop, effort = 0.2, t_max = 100, dt=0.1, initial_n = params_loop@initial_n ,
                      initial_n_pp = params_loop@initial_n_pp, progress_bar = F)
}
saveRDS(optim_loop,"optimParallel_Rmax_rcp45.RDS")

sim_loop_rcp45 <- project(params_loop, effort = 0.2, t_max = 700, progress_bar = T)
plotSummary(sim_loop_rcp45)
plotlyBiomass(sim_loop_rcp45)

params_loop@initial_n <- sim_loop_rcp45@n[dim(sim_loop@n)[1],,]
params_loop@initial_n_pp <- sim_loop_rcp45@n_pp[dim(sim_loop@n_pp)[1],]
sim_loop_rcp45 <- project(params_loop, effort = 0.2, progress_bar = T)

plotSummary(sim_loop_rcp45)
plotlyBiomass(sim_loop_rcp45)
plotFmsy(sim_loop_rcp45@params)
checkGrowth(sim_loop_rcp45)
plotDiet2(sim_loop_rcp45)

plotPredObsYield(sim_loop_rcp45,dat = avgCatch$catch_inf_rcp45)
df<-plotPredObsYield(sim_loop_rcp45,dat = avgCatch$catch_inf_rcp45,returnData = TRUE)
reg<-lm(df$obs~df$value)
summary(reg)

RDI_RDD <- as.data.frame(getRDI(sim_loop_rcp45@params)/getRDD(sim_loop_rcp45@params))
RDI_RDD$w_inf<-sim_loop_rcp45@params@species_params$w_inf
colnames(RDI_RDD)[1]<-"ratio"
plot(x = log10(RDI_RDD$w_inf),y = log10(RDI_RDD$ratio))

saveRDS(sim_loop_rcp45,"sim_rcp45.RDS")

# rcp85 model
# Read in species parameters under rcp85
sim_2013<-readRDS("sim_2013.RDS")
sim_rcp45<-readRDS("sim_rcp45.RDS")
nefsc_params_rcp85<-read_csv("nefsc_species_params_rcp85_scenario3.csv")

# Read in interaction matrix under rcp85
inter_rcp85<-read_csv("interactions_rcp85_3.csv")
rownames(inter_rcp85)<-colnames(inter_rcp85)

# Set up mizer parameters, using background resources in sim_Rmax_final
params_rcp85<-newMultispeciesParams(nefsc_params_rcp85,inter_rcp85,
                                    kappa = 1.12*sim_2013@params@resource_params$kappa)

params_rcp85@species_params$gamma<-sim_2013@params@species_params$gamma

params_rcp85<-setParams(params_rcp85)

# Specify fishing gear for each species.Target species information from NOAA website
# and mizer tutorial
gear_params(params_rcp85)$gear<-c("Pelagic","Pelagic","Industrial","Pelagic","Pelagic","Otter","Beam","Beam","Otter",
                                  "Otter","Otter","Otter","Beam","Beam","Otter","Otter","Otter","Otter","Industrial")

# Calculate R_max based on kappa and w_inf: rcp85
params_rcp85@species_params$R_max <- params_rcp85@resource_params$kappa*params_rcp85@species_params$w_inf^-1

# Use erepro from the RCP45 model as starting points
params_rcp85@species_params$erepro<-sim_rcp45@params@species_params$erepro

# Run simulation, using the end species abundance in sim_Rmax_final as initial abundance 
sim_rcp85<- project(params_rcp85, effort = 0.2, t_max = 300,
                    initial_n = sim_2013@n[100,,])

plotSummary(sim_rcp85)

plotlyBiomass(sim_rcp85)

# Tweak erepro to get reasonable Fmsy curves
params_optim_rcp85<-sim_rcp85@params

# 1st round 
params_optim_rcp85@species_params$erepro[6]<-10^-4.7
params_optim_rcp85@species_params$erepro[7]<-10^-3
params_optim_rcp85@species_params$erepro[8]<-10^-3.5
params_optim_rcp85@species_params$erepro[9]<-10^-5
params_optim_rcp85@species_params$erepro[12]<-10^-4.2

# 2nd round
params_optim_rcp85@species_params$erepro[7]<-10^-3.2
params_optim_rcp85@species_params$erepro[12]<-10^-4.5

plotFmsy(params_optim_rcp85)

sim_optim_rcp85<-project(params_optim_rcp85,effort = 0.2,t_max = 100)
plot(sim_optim_rcp85)

# Optimization of Rmax with estimated time-averaged catch data
sim_loop <- sim_optim_rcp85
params_loop <- sim_loop@params

for(i in 1:5)
{
  # saving the last time step in the param object
  params_loop@initial_n <- sim_loop@n[dim(sim_loop@n)[1],,]
  params_loop@initial_n_pp <- sim_loop@n_pp[dim(sim_loop@n_pp)[1],]
  params_calibration <- params_loop
  vary <-  log10(params_calibration@species_params$R_max)
  params_calibration<-setParams(params_calibration)
  noCores <- detectCores() - 1 # keep a spare core
  cl <- makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, as.list(ls()))
  clusterEvalQ(cl, {
    library(mizerExperimental)
    library(optimParallel)
  })
  optim_loop <-optimParallel::optimParallel(par=vary,getError,params=params_calibration, dat = avgCatch$catch_inf_rcp85, 
                                            method   ="L-BFGS-B",lower=c(rep(3,19)),upper= c(rep(15,19)),
                                            parallel=list(loginfo=TRUE, forward=TRUE))
  stopCluster(cl)
  
  # optim values:
  params_loop@species_params$R_max <- 10^optim_loop$par 
  # set the param object 
  params_loop <-setParams(params_loop)
  sim_loop <- project(params_loop, effort = 0.2, t_max = 100, dt=0.1, initial_n = params_loop@initial_n ,
                      initial_n_pp = params_loop@initial_n_pp, progress_bar = F)
}
saveRDS(optim_loop,"optimParallel_Rmax_rcp85.RDS")

params_loop@species_params$R_max <- 10^optim_loop$par 
params_loop <-setParams(params_loop)

sim_loop_rcp85 <- project(params_loop, effort = 0.2, t_max = 300, progress_bar = TRUE)
plotSummary(sim_loop_rcp85)
plotlyBiomass(sim_loop_rcp85)

params_loop@initial_n <- sim_loop_rcp85@n[dim(sim_loop@n)[1],,]
params_loop@initial_n_pp <- sim_loop_rcp85@n_pp[dim(sim_loop@n_pp)[1],]
sim_loop_rcp85 <- project(params_loop, effort = 0.2, progress_bar = TRUE)

plotSummary(sim_loop_rcp85)
plotlyBiomass(sim_loop_rcp85)
plotFmsy(sim_loop_rcp85@params)
checkGrowth(sim_loop_rcp85)
plotDiet2(sim_loop_rcp85)

plotPredObsYield(sim_loop_rcp85,dat = avgCatch$catch_inf_rcp85)
df<-plotPredObsYield(sim_loop_rcp85,dat = avgCatch$catch_inf_rcp85,returnData = TRUE)
reg<-lm(df$obs~df$value)
summary(reg)

RDI_RDD <- as.data.frame(getRDI(sim_loop_rcp85@params)/getRDD(sim_loop_rcp85@params))
RDI_RDD$w_inf<-sim_loop_rcp85@params@species_params$w_inf
colnames(RDI_RDD)[1]<-"ratio"
plot(x = log10(RDI_RDD$w_inf),y = log10(RDI_RDD$ratio))

saveRDS(sim_loop_rcp85,"sim_rcp85.RDS")

### Plotting results
# Plotting size spectra of each species in 2013, rcp45, and rcp85
rm(list=ls())

setwd("~/Google Drive/hybrid-SDM")
source("mizer_misc_functions.R")

species_param<-read_csv("nefsc_species_params_scenario3.csv")

species_param<-species_param[order(species_param$w_inf),]

avgCatch<-read_csv("catch2.csv")

sim_2013<-readRDS("sim_2013.RDS")
sim_rcp45<-readRDS("sim_rcp45.RDS")
sim_rcp85<-readRDS("sim_rcp85.RDS")

final_n_2013<-N(sim_2013)[idxFinalT(sim_2013),,,drop=FALSE]
str(final_n_2013)
nf_2013<-melt(final_n_2013)
nf_2013<-filter(nf_2013,value>0)
nf_2013$scenario<-"rcp2013"

final_n_rcp45<-N(sim_rcp45)[idxFinalT(sim_rcp45),,,drop=FALSE]
str(final_n_rcp45)
nf_rcp45<-melt(final_n_rcp45)
nf_rcp45<-filter(nf_rcp45,value>0)
nf_rcp45$scenario<-"rcp45"

final_n_rcp85<-N(sim_rcp85)[idxFinalT(sim_rcp85),,,drop=FALSE]
str(final_n_rcp85)
nf_rcp85<-melt(final_n_rcp85)
nf_rcp85<-filter(nf_rcp85,value>0)
nf_rcp85$scenario<-"rcp85"

nf<-rbind(nf_2013,nf_rcp45,nf_rcp85)

write_csv(nf,"size_spectra.csv")

nf<-read_csv(paste0("~/Google Drive/hybrid-SDM/size_spectra_fishing",fishing,".csv"))

# ## Abundance ~ size 
# for (i in 1:19){
# 
# focal_sp<-sim_2013@params@species_params$species[i]
# 
# df<-nf[nf$sp==focal_sp,]
# 
# Abnc_2013<-sum(df[df$scenario=="rcp2013",]$value)
# Abnc_rcp45<-sum(df[df$scenario=="rcp45",]$value)
# Abnc_rcp85<-sum(df[df$scenario=="rcp85",]$value)
# 
# avgMass_2013<-sum(df[df$scenario=="rcp2013",]$w*df[df$scenario=="rcp2013",]$value)/sum(df[df$scenario=="rcp2013",]$value)
# avgMass_rcp45<-sum(df[df$scenario=="rcp45",]$w*df[df$scenario=="rcp45",]$value)/sum(df[df$scenario=="rcp45",]$value)
# avgMass_rcp85<-sum(df[df$scenario=="rcp85",]$w*df[df$scenario=="rcp85",]$value)/sum(df[df$scenario=="rcp85",]$value)
# 
# pdf(file = paste0("~/Desktop/",focal_sp,".pdf"))
# 
# plot<-ggplot(df)+
#   theme_bw()+
#   theme(axis.text.x = element_text(face="bold",size=20),
#         axis.text.y = element_text(face="bold",size=20),
#         legend.position = "none")+
#   labs(x="",y="")+
#   geom_line(aes(x = w,y = value,color=scenario),size=1,linetype=1)+
#   geom_hline(aes(yintercept = Abnc_2013),color = "red",size=2,linetype=2,alpha = 0.4)+
#   geom_hline(aes(yintercept = Abnc_rcp45),color = "green",size=2,linetype=2,alpha = 0.4)+
#   geom_hline(aes(yintercept = Abnc_rcp85),color = "blue",size=2,linetype=2,alpha = 0.4)+
#   geom_vline(aes(xintercept = avgMass_2013),color = "red",size = 2, linetype = 4, alpha = 0.4)+
#   geom_vline(aes(xintercept = avgMass_rcp45),color = "green",size = 2, linetype = 4, alpha = 0.4)+
#   geom_vline(aes(xintercept = avgMass_rcp85),color = "blue",size = 2, linetype = 4, alpha = 0.4)+
#   scale_x_log10()+
#   scale_y_log10(limits = c(0.001,max(Abnc_2013,Abnc_rcp45,Abnc_rcp85)))
# 
# print(plot)
# 
# dev.off()
# 
# }
# 
# ## Relative abundance ~ size
# 
# for (i in 1:19){
#   
#   focal_sp<-sim_2013@params@species_params$species[i]
#   
#   df<-nf[nf$sp==focal_sp,]
#   
#   df$relative<--99
#   
#   nrow_2013<-nrow(df[df$scenario=="rcp2013",])
#   nrow_rcp45<-nrow(df[df$scenario=="rcp45",])
#   nrow_rcp85<-nrow(df[df$scenario=="rcp85",])
#   
#   df$relative[1:nrow_2013]<-df$value[1:nrow_2013]/df$value[1]
#   df$relative[(nrow_2013+1):(nrow_2013+nrow_rcp45)]<-df$value[(nrow_2013+1):(nrow_2013+nrow_rcp45)]/df$value[nrow_2013+1]
#   df$relative[(nrow(df)-nrow_rcp85+1):nrow(df)]<-df$value[(nrow(df)-nrow_rcp85+1):nrow(df)]/df$value[nrow(df)-nrow_rcp85+1]
#   
#   pdf(file = paste0("~/Desktop/",focal_sp,".pdf"))
# 
#   plot<-ggplot(df)+
#     theme_bw()+
#     theme(axis.text.x = element_text(face="bold",size=20),
#           axis.text.y = element_text(face="bold",size=20),
#           legend.position = "none")+
#     labs(x="",y="")+
#     geom_line(aes(x = w,y = relative,color=scenario),size=1,linetype=1)+
#     scale_x_log10()+
#     scale_y_log10(limits=c(1e-7,1))
# 
#   print(plot)
# 
#   dev.off()
#   
# }

# Biomass density ~ size 
rm(list=ls())

setwd("~/Google Drive/hybrid-SDM")
source("mizer_misc_functions.R")
library(gridExtra)

species_param<-read_csv("nefsc_species_params_scenario3.csv")

species_param<-species_param[order(species_param$w_inf),]

avgCatch<-read_csv("catch2.csv")

sim_2013<-readRDS("sim_2013.RDS")
sim_rcp45<-readRDS("sim_rcp45.RDS")
sim_rcp85<-readRDS("sim_rcp85.RDS")

df_2013<-plotSpectra(sim_2013)$data
df_rcp45<-plotSpectra(sim_rcp45)$data
df_rcp85<-plotSpectra(sim_rcp85)$data

nf<-rbind(df_2013,df_rcp45,df_rcp85)
nf$scenario<-rep(c("rcp2013","rcp45","rcp85"),times = c(1667,1635,1631))

write_csv(nf,"size_spectra_biomassdensity.csv")

nf<-read_csv("size_spectra_biomassdensity.csv")
nf2<-read_csv("size_spectra_abundance.csv")

plot_sizeSpectra<-function(focal_sp){
  
  df<-nf[nf$Species==focal_sp,]
  df2<-nf2[nf2$sp==focal_sp,]
  
  biomass_2013<-sum(df[df$scenario=="rcp2013",]$value)
  biomass_rcp45<-sum(df[df$scenario=="rcp45",]$value)
  biomass_rcp85<-sum(df[df$scenario=="rcp85",]$value)
  
  avgMass_2013<-sum(df2[df2$scenario=="rcp2013",]$w*df2[df2$scenario=="rcp2013",]$value)/sum(df2[df2$scenario=="rcp2013",]$value)
  avgMass_rcp45<-sum(df2[df2$scenario=="rcp45",]$w*df2[df2$scenario=="rcp45",]$value)/sum(df2[df2$scenario=="rcp45",]$value)
  avgMass_rcp85<-sum(df2[df2$scenario=="rcp85",]$w*df2[df2$scenario=="rcp85",]$value)/sum(df2[df2$scenario=="rcp85",]$value)
  
  ggplot(df)+
    theme_bw()+
    theme(axis.text.x = element_text(face="bold",size=10),
          axis.text.y = element_text(face="bold",size=10),
          legend.position = "none")+
    labs(x="",y="")+
    geom_line(aes(x = w,y = value,color=scenario),size=1,linetype=1)+
    geom_hline(aes(yintercept = biomass_2013),color = "red",size=1,linetype=2,alpha = 0.4)+
    geom_hline(aes(yintercept = biomass_rcp45),color = "green",size=1,linetype=2,alpha = 0.4)+
    geom_hline(aes(yintercept = biomass_rcp85),color = "blue",size=1,linetype=2,alpha = 0.4)+
    geom_vline(aes(xintercept = avgMass_2013),color = "red",size = 1, linetype = 4, alpha = 0.4)+
    geom_vline(aes(xintercept = avgMass_rcp45),color = "green",size = 1, linetype = 4, alpha = 0.4)+
    geom_vline(aes(xintercept = avgMass_rcp85),color = "blue",size = 1, linetype = 4, alpha = 0.4)+
    scale_x_log10()+
    scale_y_log10(limits = c(NA,max(biomass_2013,biomass_rcp45,biomass_rcp85)))
  
}

focal_sp<-species_param$species

p<-list()

for (i in 1:19){
  
  p[[i]]<-plot_sizeSpectra(focal_sp[i])
  
}

do.call(grid.arrange,p)

# Getting data of total biomass and average mass 
df_all<-data.frame(matrix(ncol=7,nrow=length(focal_sp)))

df_all[,1]<-focal_sp
colnames(df_all)[1]<-"species"

for (i in 1:length(focal_sp)){
  
  df<-nf[nf$Species==focal_sp[i],]
  df2<-nf2[nf2$sp==focal_sp[i],]
  
  df_all[i,2]<-sum(df[df$scenario=="rcp2013",]$value)
  df_all[i,3]<-sum(df[df$scenario=="rcp45",]$value)
  df_all[i,4]<-sum(df[df$scenario=="rcp85",]$value)
  
  df_all[i,5]<-sum(df2[df2$scenario=="rcp2013",]$w*df2[df2$scenario=="rcp2013",]$value)/sum(df2[df2$scenario=="rcp2013",]$value)
  df_all[i,6]<-sum(df2[df2$scenario=="rcp45",]$w*df2[df2$scenario=="rcp45",]$value)/sum(df2[df2$scenario=="rcp45",]$value)
  df_all[i,7]<-sum(df2[df2$scenario=="rcp85",]$w*df2[df2$scenario=="rcp85",]$value)/sum(df2[df2$scenario=="rcp85",]$value)  
  
}

colnames(df_all)[2:7]<-c("biomass_2013","biomass_rcp45","biomass_rcp85",
                         "avgMass_2013","avgMass_rcp45","avgMass_rcp85")

write_csv(df_all,"~/Desktop/df.csv")

# Plotting relative biomass density ~ size
nf<-read_csv("size_spectra_biomassdensity.csv")

plot_relativeSizeSpectra<-function(focal_sp){
  
  df<-nf[nf$Species==focal_sp,]
  
  df$relative<--99
  
  nrow_2013<-nrow(df[df$scenario=="rcp2013",])
  nrow_rcp45<-nrow(df[df$scenario=="rcp45",])
  nrow_rcp85<-nrow(df[df$scenario=="rcp85",])
  
  df$relative[1:nrow_2013]<-df$value[1:nrow_2013]/df$value[1]
  df$relative[(nrow_2013+1):(nrow_2013+nrow_rcp45)]<-df$value[(nrow_2013+1):(nrow_2013+nrow_rcp45)]/df$value[nrow_2013+1]
  df$relative[(nrow(df)-nrow_rcp85+1):nrow(df)]<-df$value[(nrow(df)-nrow_rcp85+1):nrow(df)]/df$value[nrow(df)-nrow_rcp85+1]
  
  ggplot(df)+
    theme_bw()+
    theme(axis.text.x = element_text(face="bold",size=10),
          axis.text.y = element_text(face="bold",size=10),
          legend.position = "none")+
    labs(x="",y="")+
    geom_line(aes(x = w,y = relative,color=scenario),size=1,linetype=1)+
    scale_x_log10()+
    scale_y_log10()
  
}

focal_sp<-species_param$species

p<-list()

for (i in 1:19){
  
  p[[i]]<-plot_relativeSizeSpectra(focal_sp[i])
  
}

do.call(grid.arrange,p)

# Plotting community size spectrum
## Biomass density ~ size
nf<-read_csv("size_spectra_biomassdensity.csv")
df<-nf[nf$Species!="Resource",]

now<-subset(df,subset = scenario == "rcp2013")
rcp45<-subset(df,subset = scenario == "rcp45")
rcp85<-subset(df,subset = scenario == "rcp85")

plotdata_now<-aggregate(now$value~now$w,FUN = sum)
plotdata_now$relative<-plotdata_now$`now$value`/plotdata_now$`now$value`[1]
colnames(plotdata_now)<-c("size","biomass_density","relative_biomass_density")

plotdata_rcp45<-aggregate(rcp45$value~rcp45$w,FUN = sum)
plotdata_rcp45$relative<-plotdata_rcp45$`rcp45$value`/plotdata_rcp45$`rcp45$value`[1]
colnames(plotdata_rcp45)<-c("size","biomass_density","relative_biomass_density")

plotdata_rcp85<-aggregate(rcp85$value~rcp85$w,FUN = sum)
plotdata_rcp85$relative<-plotdata_rcp85$`rcp85$value`/plotdata_rcp85$`rcp85$value`[1]
colnames(plotdata_rcp85)<-c("size","biomass_density","relative_biomass_density")

plotdata<-rbind(plotdata_now,plotdata_rcp45,plotdata_rcp85)
plotdata$scenario<-rep(c("rcp2013","rcp45","rcp85"),each = 100)

ggplot(data = plotdata)+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",size=25),
        axis.text.y = element_text(face="bold",size=25),
        legend.position = "none")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x = "",y = "")+
  geom_line(aes(x = size,y = biomass_density,color = scenario),size = 2)

## Relative biomass density ~ size
ggplot(data = plotdata)+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",size=25),
        axis.text.y = element_text(face="bold",size=25),
        legend.position = "none")+
  scale_x_log10()+
  scale_y_log10()+
  labs(x = "",y = "")+
  geom_line(aes(x = size,y = relative_biomass_density,color = scenario),size = 2)

# Plotting relative change in total abundance and mean body size
rm(list=ls())
setwd("~/Google Drive/hybrid-SDM/")

library(tidyverse)
library(ggplot2)

spectra<-read_csv("size_spectra.csv")
speciesdata<-read_csv("nefsc_species_params_scenario3.csv")

sp<-levels(as.factor(spectra$sp))

dAbnc_rcp45<-NULL
dAbnc_rcp85<-NULL
dMass_rcp45<-NULL
dMass_rcp85<-NULL

for (i in 1:19){
  
  df<-spectra[spectra$sp==sp[i],]
  df_2013<-df[df$scenario=="rcp2013",]
  df_rcp45<-df[df$scenario=="rcp45",]
  df_rcp85<-df[df$scenario=="rcp85",]
  
  dAbnc_rcp45[i]<-(sum(df_rcp45$value)-sum(df_2013$value))/sum(df_2013$value)
  dAbnc_rcp85[i]<-(sum(df_rcp85$value)-sum(df_2013$value))/sum(df_2013$value)
  
  avgMass_2013<-sum(df_2013$w*df_2013$value)/sum(df_2013$value)
  avgMass_rcp45<-sum(df_rcp45$w*df_rcp45$value)/sum(df_rcp45$value)
  avgMass_rcp85<-sum(df_rcp85$w*df_rcp85$value)/sum(df_rcp85$value)
  
  dMass_rcp45[i]<-(avgMass_rcp45-avgMass_2013)/avgMass_2013
  dMass_rcp85[i]<-(avgMass_rcp85-avgMass_2013)/avgMass_2013
  
}

data<-data.frame(cbind(dAbnc_rcp45,dAbnc_rcp85,dMass_rcp45,dMass_rcp85))
impact<-cbind(sp,speciesdata$w_inf,data)
colnames(impact)[2]<-"w_inf"

# Manually add habitat and range variables before saving it
write_csv(impact,"warming_impact.csv")

# Plotting relative change in total abundance and mean body size against max size
rm(list=ls())

library(tidyverse)
library(ggplot2)

setwd("~/Google Drive/hybrid-SDM/")

df<-read_csv("warming_impact.csv")

## Relative change in abundance
# RCP45
reg1<-lm(log10(dAbnc_rcp45+1)~log10(w_inf),data = df)

plot_abnc_rcp45<-ggplot(data = df)+
  theme_bw()+
  theme(text = element_text(size=20),legend.position = "none")+
  labs(x = "",y = "")+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(x = w_inf,y = dAbnc_rcp45+1,
                 color = as.factor(habitat_type),shape = as.factor(range)),size=3)+
  scale_shape_manual(values= c(1,19))+
  scale_color_manual(values = c("purple","black","orange"))+
  geom_hline(yintercept = 1,color="red",linetype=2)+
  geom_abline(slope = reg1$coefficients[2],intercept = reg1$coefficients[1]+1,color = "gray")

# plot_abnc_rcp45
summary(reg1)

# RCP85
reg1<-lm(log10(dAbnc_rcp85+1)~log10(w_inf),data = df)

plot_abnc_rcp85<-ggplot(data = df)+
  theme_bw()+
  theme(text = element_text(size=20),legend.position = "none")+
  labs(x = "",y = "")+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(x = w_inf,y = dAbnc_rcp85+1,
                 color = as.factor(habitat_type),shape = as.factor(range)),size=3)+
  scale_shape_manual(values= c(1,19))+
  scale_color_manual(values = c("purple","black","orange"))+
  geom_hline(yintercept = 1,color="red",linetype=2)+
  geom_abline(slope = reg1$coefficients[2],intercept = reg1$coefficients[1]+1,color = "gray")

# plot_abnc_rcp85
summary(reg1)

## Relative change in mean size
# RCP45
reg2<-lm(log10(dMass_rcp45+1)~log10(w_inf),data = df)

plot_mass_rcp45<-ggplot(data = df)+
  theme_bw()+
  theme(text = element_text(size=20),legend.position = "none")+
  labs(x = "",y = "")+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(x = w_inf,y = dMass_rcp45+1,
                 color = as.factor(habitat_type),shape = as.factor(range)),size=3)+
  scale_shape_manual(values= c(1,19))+
  scale_color_manual(values = c("purple","black","orange"))+
  geom_hline(yintercept = 1,color="red",linetype=2)+
  geom_abline(slope = reg2$coefficients[2],intercept = reg2$coefficients[1]+1,color = "gray")

# plot_mass_rcp45
summary(reg2)

# RCP85
reg2<-lm(log10(dMass_rcp85+1)~log10(w_inf),data = df)

plot_mass_rcp85<-plot_abnc<-ggplot(data = df)+
  theme_bw()+
  theme(text = element_text(size=20),legend.position = "none")+
  labs(x = "",y = "")+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(x = w_inf,y = dMass_rcp85+1,
                 color = as.factor(habitat_type),shape = as.factor(range)),size=3)+
  scale_shape_manual(values= c(1,19))+
  scale_color_manual(values = c("purple","black","orange"))+
  geom_hline(yintercept = 1,color="red",linetype=2)+
  geom_abline(slope = reg2$coefficients[2],intercept = reg2$coefficients[1]+1,color = "gray")

# plot_mass_rcp85
summary(reg2)

cowplot::plot_grid(plot_abnc_rcp45,plot_abnc_rcp85,plot_mass_rcp45,plot_mass_rcp85,ncol=2,align = "hv",rel_widths = c(0.5,0.5))


# Plot relative change in biomass ~ changes in w_inf, vb, and suitable grids
rm(list=ls())

setwd("~/Google Drive/hybrid-SDM")

library(mizerHowTo)
library(tidyverse)
library(gridExtra)
library(car)

biomass<-read_csv("biomass.csv")
speciesdata<-read_csv("speciesdata_SDM.csv")

sp<-levels(as.factor(speciesdata$species))

dBiomass_rcp45<-NULL
dBiomass_rcp85<-NULL

dW_inf_rcp45<-NULL
dW_inf_rcp85<-NULL

dVb_rcp45<-NULL
dVb_rcp85<-NULL

dGrid_rcp45<-NULL
dGrid_rcp85<-NULL

for (i in 1:19){
  
  biomass_2013<-biomass[biomass$species==sp[i]&biomass$scenario=="2013",]$value
  biomass_rcp45<-biomass[biomass$species==sp[i]&biomass$scenario=="rcp45",]$value
  biomass_rcp85<-biomass[biomass$species==sp[i]&biomass$scenario=="rcp85",]$value
  
  dBiomass_rcp45[i]<-(biomass_rcp45-biomass_2013)/biomass_2013
  dBiomass_rcp85[i]<-(biomass_rcp85-biomass_2013)/biomass_2013
  
  w_inf_2013<-speciesdata[speciesdata$species==sp[i],]$maxweight_2013
  w_inf_rcp45<-speciesdata[speciesdata$species==sp[i],]$maxweight_RCP45
  w_inf_rcp85<-speciesdata[speciesdata$species==sp[i],]$maxweight_RCP85
  
  dW_inf_rcp45[i]<-(w_inf_rcp45-w_inf_2013)/w_inf_2013
  dW_inf_rcp85[i]<-(w_inf_rcp85-w_inf_2013)/w_inf_2013
  
  vb_2013<-speciesdata[speciesdata$species==sp[i],]$vonBert_2013
  vb_rcp45<-speciesdata[speciesdata$species==sp[i],]$vonBert_RCP45
  vb_rcp85<-speciesdata[speciesdata$species==sp[i],]$vonBert_RCP85
  
  dVb_rcp45[i]<-(vb_rcp45-vb_2013)/vb_2013
  dVb_rcp85[i]<-(vb_rcp85-vb_2013)/vb_2013
  
  grid_2013<-speciesdata[speciesdata$species==sp[i],]$suitable_2013
  grid_rcp45<-speciesdata[speciesdata$species==sp[i],]$suitable_rcp45
  grid_rcp85<-speciesdata[speciesdata$species==sp[i],]$suitable_rcp85
  
  dGrid_rcp45[i]<-(grid_rcp45-grid_2013)/grid_2013
  dGrid_rcp85[i]<-(grid_rcp85-grid_2013)/grid_2013
  
}

df<-data.frame(cbind(sp,dBiomass_rcp45,dBiomass_rcp85,dW_inf_rcp45,dW_inf_rcp85,
                     dVb_rcp45,dVb_rcp85,dGrid_rcp45,dGrid_rcp85))

df[,2:9]<-as.numeric(df[,2:9])

reg1_rcp45<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dW_inf_rcp45)+1))
reg1_rcp85<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dW_inf_rcp85)+1))


p1<-ggplot(data = df)+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.title = element_text(size=15,face = "bold"))+
  labs(x = "Relative change in maximum size",y = "Relative change in biomass",size=10)+
  geom_point(aes(x = as.numeric(dW_inf_rcp45)+1,y = as.numeric(dBiomass_rcp45)+1),size=3)+
  geom_point(aes(x = as.numeric(dW_inf_rcp85)+1,y = as.numeric(dBiomass_rcp85)+1),pch = 1,size=3)+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(aes(yintercept = 1),color = "red",linetype=2)+
  geom_abline(intercept = coef(reg1_rcp45)[1]+1,slope = coef(reg1_rcp45)[2],color = "black",alpha=0.4,size=1)+
  geom_abline(intercept = coef(reg1_rcp85)[1]+1,slope = coef(reg1_rcp85)[2],linetype=2,color = "black",alpha=0.4,size=1)

reg2_rcp45<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dVb_rcp45)+1))
reg2_rcp85<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dVb_rcp85)+1))

p2<-ggplot(data = df)+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.title = element_text(size=15,face = "bold"))+
  labs(x = "Relative change in growth coefficient",y = "Relative change in biomass",size=10)+
  geom_point(aes(x = as.numeric(dVb_rcp45)+1,y = as.numeric(dBiomass_rcp45)+1),size=3)+
  geom_point(aes(x = as.numeric(dVb_rcp85)+1,y = as.numeric(dBiomass_rcp85)+1),pch = 1,size=3)+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(aes(yintercept = 1),color = "red",linetype=2)+
  geom_abline(intercept = coef(reg2_rcp45)[1]+1,slope = coef(reg2_rcp45)[2],color = "black",alpha=0.4,size=1)+
  geom_abline(intercept = coef(reg2_rcp85)[1]+1,slope = coef(reg2_rcp85)[2],linetype=2,color = "black",alpha=0.4,size=1)


reg3_rcp45<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dGrid_rcp45)+1))
reg3_rcp85<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dGrid_rcp85)+1))

p3<-ggplot(data = df)+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.title = element_text(size=15,face = "bold"))+
  labs(x = "Relative change in habitat availability",y = "Relative change in biomass",size=10)+
  geom_point(aes(x = as.numeric(dGrid_rcp45)+1,y = as.numeric(dBiomass_rcp45)+1),size=3)+
  geom_point(aes(x = as.numeric(dGrid_rcp85)+1,y = as.numeric(dBiomass_rcp85)+1),pch = 1,size=3)+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(aes(yintercept = 1),color = "red",linetype=2)+
  geom_abline(intercept = coef(reg3_rcp45)[1]+1,slope = coef(reg3_rcp45)[2],color = "black",alpha=0.4,size=1)+
  geom_abline(intercept = coef(reg3_rcp85)[1]+1,slope = coef(reg3_rcp85)[2],linetype=2,color = "black",alpha=0.4,size=1)

grid.arrange(p1,p2,p3,ncol=3)

# Some formal stats
# RCP45
reg_rcp45_all<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dW_inf_rcp45)+1)+
                log10(as.numeric(dVb_rcp45)+1)+log10(as.numeric(dGrid_rcp45)+1),data = df)

reg_rcp45_red1<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dW_inf_rcp45)+1),data = df)

reg_rcp45_red2<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dVb_rcp45)+1),data = df)

reg_rcp45_red3<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dGrid_rcp45)+1),data = df)

# Full - dW_inf
reg_rcp45_red4<-lm(log10(as.numeric(dBiomass_rcp45)+1)~
                     log10(as.numeric(dVb_rcp45)+1)+log10(as.numeric(dGrid_rcp45)+1),data = df)

# Full - dVb
reg_rcp45_red5<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dW_inf_rcp45)+1)
                     +log10(as.numeric(dGrid_rcp45)+1),data = df)

# Full - dGrid
reg_rcp45_red6<-lm(log10(as.numeric(dBiomass_rcp45)+1)~log10(as.numeric(dW_inf_rcp45)+1)+
                     log10(as.numeric(dVb_rcp45)+1),data = df)


AIC(reg_rcp45_all,reg_rcp45_red1,reg_rcp45_red2,reg_rcp45_red3,reg_rcp45_red4,reg_rcp45_red5,reg_rcp45_red6)

Anova(reg_rcp45_all,type = 3)


# RCP85
reg_rcp85_all<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dW_inf_rcp85)+1)+
                    log10(as.numeric(dVb_rcp85)+1)+log10(as.numeric(dGrid_rcp85)+1),data = df)

reg_rcp85_red1<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dW_inf_rcp85)+1),data = df)

reg_rcp85_red2<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dVb_rcp85)+1),data = df)

reg_rcp85_red3<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dGrid_rcp85)+1),data = df)

# Full - dW_inf
reg_rcp85_red4<-lm(log10(as.numeric(dBiomass_rcp85)+1)~
                     log10(as.numeric(dVb_rcp85)+1)+log10(as.numeric(dGrid_rcp85)+1),data = df)

# Full - dVb
reg_rcp85_red5<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dW_inf_rcp85)+1)
                   +log10(as.numeric(dGrid_rcp85)+1),data = df)

# Full - dGrid
reg_rcp85_red6<-lm(log10(as.numeric(dBiomass_rcp85)+1)~log10(as.numeric(dW_inf_rcp85)+1)+
                     log10(as.numeric(dVb_rcp85)+1),data = df)




AIC(reg_rcp85_all,reg_rcp85_red1,reg_rcp85_red2,reg_rcp85_red3,reg_rcp85_red4,reg_rcp85_red5,reg_rcp85_red6)

summary(reg_rcp85_all)
