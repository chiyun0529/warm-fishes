### Run "Control" models where life-history changes and interaction are not taken into account
rm(list=ls())
# Load packages and some additional functions

# If run on Mac laptop
setwd("~/Google Drive/hybrid-SDM")

# If run on Windows Desktop
setwd("G:/My Drive/hybrid-SDM")

# Load necessary packages and functions
source("mizer_misc_functions.R")
library(mizerHowTo)
library(mizerExperimental)
library(tidyverse)

sim_2013<-readRDS("sim_2013.RDS")
sim_rcp45<-readRDS("sim_rcp45.RDS")
sim_rcp85<-readRDS("sim_rcp85.RDS")

avgCatch<-read_csv("catch2.csv")

nefsc_params_2013<-read_csv("nefsc_species_params_scenario3.csv")
nefsc_params_rcp45<-read_csv("nefsc_species_params_rcp45_scenario3.csv")
nefsc_params_rcp85<-read_csv("nefsc_species_params_rcp85_scenario3.csv")

inter_2013<-read_csv("interactions_2013_3.csv")
rownames(inter_2013)<-colnames(inter_2013)

inter_rcp45<-read_csv("interactions_rcp45_3.csv")
rownames(inter_rcp45)<-colnames(inter_rcp45)

inter_rcp85<-read_csv("interactions_rcp85_3.csv")
rownames(inter_rcp85)<-colnames(inter_rcp45)

## Control model 1: w/ life-history changes, w/o in interaction changes
## RCP45
params_rcp45_nointer<-sim_rcp45_nointer@params
## sim_rcp85_nointer was a good parameter starting point
params_rcp45_nointer@species_params<-sim_rcp85_nointer@params@species_params

plotFmsy(params_rcp45_nointer)

sim_rcp45_nointer<-project(params_rcp45_nointer,effort = 0.2,t_max = 100)

plotSummary(sim_rcp45_nointer)
plotlyBiomass(sim_rcp45_nointer)
plotPredObsYield(sim_rcp45_nointer,avgCatch$catch_inf_rcp45,returnData = TRUE)

sim_loop <- sim_rcp45_nointer
params_loop <- sim_rcp45_nointer@params

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
saveRDS(optim_loop,"optimParallel_Rmax2_rcp45nointer.RDS")

sim_loop<-project(params_loop,effort = 0.2,t_max = 300)
plotSummary(sim_loop)

params_loop@initial_n <- sim_loop@n[dim(sim_loop@n)[1],,]
params_loop@initial_n_pp <- sim_loop@n_pp[dim(sim_loop@n_pp)[1],]
sim_loop <- project(params_loop, effort = 0.2)

plotSummary(sim_loop)
plotPredObsYield(sim_loop,dat = avgCatch$catch_inf_rcp45,returnData = TRUE)
plotFmsy(sim_loop@params)
checkGrowth(sim_loop)

saveRDS(sim_loop,"sim_rcp45_nointer.RDS")

## Control model 1: w/ life-history changes, w/o in interaction changes
## RCP85
## Everything looked normal so no optimization is necessary
params_rcp85_nointer<-setParams(sim_rcp85@params,interaction = inter_2013)

sim_rcp85_nointer<-project(params_rcp85_nointer,effort = 0.2,t_max = 100)

plotSummary(sim_rcp85_nointer)
plotPredObsYield(sim_rcp85_nointer,avgCatch$catch_inf_rcp85,returnData = TRUE)

sim_loop <- sim_rcp85_nointer
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
saveRDS(optim_loop,"optimParallel_Rmax_rcp85nointer.RDS")

params_loop@species_params$gamma<-sim_rcp85@params@species_params$gamma

sim_loop_nointer<-project(params_loop,effort = 0.2,t_max = 300,progress_bar = TRUE)

plotSummary(sim_loop_nointer)
plotlyBiomass(sim_loop_nointer)

checkGrowth(sim_loop_nointer)

saveRDS(sim_loop_nointer,"sim_rcp85_nointer.RDS")


## Control model 2: w/o life-history changes, w/ in interaction changes
## RCP45
params_rcp45_nolife<-newMultispeciesParams(nefsc_params_2013,interaction = inter_rcp45,kappa = sim_rcp45@params@resource_params$kappa)

params_rcp45_nolife@species_params$R_max<-params_rcp45_nolife@resource_params$kappa*params_rcp45_nolife@species_params$w_inf^-1

gear_params(params_rcp45_nolife)$gear<-c("Pelagic","Pelagic","Industrial","Pelagic","Pelagic","Otter","Beam",
                                    "Beam","Otter","Otter","Otter","Otter","Beam","Beam","Otter","Otter","Otter","Otter","Industrial")


sim_rcp45_nolife<-project(params_rcp45_nolife,effort = 0.2,t_max = 100)
plotlyBiomass(sim_rcp45_nolife)

# 1st optimization to get biomass to the right ballpark
params_optim<-params_rcp45_nolife
vary<-log10(params_optim@species_params$R_max)
params_optim<-setParams(params_optim)

noCores <- detectCores() - 1 # keep some spare core
cl <- makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
  library(mizerExperimental)
  library(optimParallel)
})
optim_result <- optimParallel::optimParallel(par=vary,getError,params=params_optim, dat = avgCatch$catch_inf_rcp45, method   ="L-BFGS-B",lower=c(rep(3,19)),upper= c(rep(15,19)),
                                             parallel=list(loginfo=TRUE, forward=TRUE))
stopCluster(cl)
saveRDS(optim_result,"optimParallel_Rmax1_rcp45nolife.RDS")

params_optim@species_params$R_max<-10^optim_result$par

params_rcp45_nolife<-setParams(params_optim)

sim_rcp45_nolife2<-project(params_rcp45_nolife,effort = 0.2,t_max = 100,dt = 0.1,
                          initial_n = sim_rcp45_nolife@n[100,,],initial_n_pp = sim_rcp45_nolife@n_pp[100,])

plotSummary(sim_rcp45_nolife2)
plotPredObsYield(sim_rcp45_nolife2,avgCatch$catch_inf_rcp45)

# Tweak erepro to get reasonable Fmsy curves
params_optim2<-sim_rcp45_nolife2@params

# 1st round
params_optim2@species_params$erepro<-sim_2013@params@species_params$erepro
params_optim2@species_params$erepro[2]<-10^-1.5
params_optim2@species_params$erepro[4]<-10^-1
params_optim2@species_params$erepro[5]<-10^-3
params_optim2@species_params$erepro[7]<-10^-2
params_optim2@species_params$erepro[8]<-10^-3.5
params_optim2@species_params$erepro[11]<-10^-4
params_optim2@species_params$erepro[12]<-10^-3
params_optim2@species_params$erepro[16]<-10^-3.5

# 2nd round
params_optim2@species_params$erepro[2]<-10^-1
params_optim2@species_params$erepro[4]<-10^-0.2
params_optim2@species_params$erepro[5]<-10^-3.5
params_optim2@species_params$erepro[7]<-10^-1
params_optim2@species_params$erepro[8]<-10^-4
params_optim2@species_params$erepro[12]<-10^-2

# 3rd round
params_optim2@species_params$erepro[2]<-10^-0.2
params_optim2@species_params$erepro[7]<-10^-0.2
params_optim2@species_params$erepro[12]<-10^-0.5

plotFmsy(params_optim2)

sim_rcp45_nolife3<-project(params_optim2,effort = 0.2,t_max = 100)

plotSummary(sim_rcp45_nolife3)

# 2nd optimization
sim_loop <- sim_rcp45_nolife3
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
saveRDS(optim_loop,"optimParallel_Rmax2_rcp45nolife.RDS")

sim_loop<-project(params_loop,effort = 0.2,t_max = 300)
plotSummary(sim_loop)

params_loop@initial_n <- sim_loop@n[dim(sim_loop@n)[1],,]
params_loop@initial_n_pp <- sim_loop@n_pp[dim(sim_loop@n_pp)[1],]
sim_loop <- project(params_loop, effort = 0.2)

plotSummary(sim_loop)
plotPredObsYield(sim_loop,dat = avgCatch$catch_inf_rcp45)

RDI_RDD <- as.data.frame(getRDI(sim_loop@params)/getRDD(sim_loop@params))
RDI_RDD$w_inf<-params_loop@species_params$w_inf
colnames(RDI_RDD)[1]<-"ratio"
plot(x = log10(RDI_RDD$w_inf),y = log10(RDI_RDD$ratio))

saveRDS(sim_loop,"sim_rcp45_nolife.RDS")


## RCP85
params_rcp85_nolife<-newMultispeciesParams(nefsc_params_2013,interaction = inter_rcp85,kappa = sim_rcp85@params@resource_params$kappa)

params_rcp85_nolife@species_params$R_max<-params_rcp85_nolife@resource_params$kappa*params_rcp85_nolife@species_params$w_inf^-1

gear_params(params_rcp85_nolife)$gear<-c("Pelagic","Pelagic","Industrial","Pelagic","Pelagic","Otter","Beam",
                                         "Beam","Otter","Otter","Otter","Otter","Beam","Beam","Otter","Otter","Otter","Otter","Industrial")


sim_rcp85_nolife<-project(params_rcp85_nolife,effort = 0.2,t_max = 100)
plotlyBiomass(sim_rcp85_nolife)

# 1st optimization to get biomass to the right ballpark
params_optim<-params_rcp85_nolife
vary<-log10(params_optim@species_params$R_max)
params_optim<-setParams(params_optim)

noCores <- detectCores() - 1 # keep some spare core
cl <- makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
  library(mizerExperimental)
  library(optimParallel)
})
optim_result <- optimParallel::optimParallel(par=vary,getError,params=params_optim, dat = avgCatch$catch_inf_rcp85, method   ="L-BFGS-B",lower=c(rep(3,19)),upper= c(rep(15,19)),
                                             parallel=list(loginfo=TRUE, forward=TRUE))
stopCluster(cl)
saveRDS(optim_result,"optimParallel_Rmax1_rcp85nolife.RDS")

params_optim@species_params$R_max<-10^optim_result$par

params_rcp85_nolife<-setParams(params_optim)

sim_rcp85_nolife2<-project(params_rcp85_nolife,effort = 0.2,t_max = 100,dt = 0.1,
                           initial_n = sim_rcp85_nolife@n[100,,],initial_n_pp = sim_rcp85_nolife@n_pp[100,])

plotSummary(sim_rcp85_nolife2)
plotPredObsYield(sim_rcp85_nolife2,avgCatch$catch_inf_rcp85)

# Tweak erepro to get reasonable Fmsy curves
params_optim2<-sim_rcp85_nolife2@params

# 1st round
params_optim2@species_params$erepro<-sim_2013@params@species_params$erepro
params_optim2@species_params$erepro[4]<-10^-1
params_optim2@species_params$erepro[5]<-10^-2.5
params_optim2@species_params$erepro[11]<-10^-4
params_optim2@species_params$erepro[12]<-10^-3
params_optim2@species_params$erepro[15]<-10^-4.7
params_optim2@species_params$erepro[16]<-10^-3.5
params_optim2@species_params$erepro[19]<-10^-2

# 2nd round
params_optim2@species_params$erepro[4]<-10^-0.2
params_optim2@species_params$erepro[5]<-10^-3
params_optim2@species_params$erepro[11]<-10^-4.5
params_optim2@species_params$erepro[12]<-10^-2
params_optim2@species_params$erepro[16]<-10^-4
params_optim2@species_params$erepro[19]<-10^-1

# 3rd round
params_optim2@species_params$erepro[5]<-10^-3.2
params_optim2@species_params$erepro[12]<-10^-0.2
params_optim2@species_params$erepro[19]<-10^-0.2
  

plotFmsy(params_optim2)

sim_rcp85_nolife3<-project(params_optim2,effort = 0.2,t_max = 100)

plotSummary(sim_rcp85_nolife3)

# 2nd optimization
sim_loop <- sim_rcp85_nolife3
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
saveRDS(optim_loop,"optimParallel_Rmax2_rcp85nolife.RDS")

sim_loop<-project(params_loop,effort = 0.2,t_max = 300)
plotSummary(sim_loop)

params_loop@initial_n <- sim_loop@n[dim(sim_loop@n)[1],,]
params_loop@initial_n_pp <- sim_loop@n_pp[dim(sim_loop@n_pp)[1],]
sim_loop <- project(params_loop, effort = 0.2)

plotSummary(sim_loop)
plotPredObsYield(sim_loop,dat = avgCatch$catch_inf_rcp85)

RDI_RDD <- as.data.frame(getRDI(sim_loop@params)/getRDD(sim_loop@params))
RDI_RDD$w_inf<-params_loop@species_params$w_inf
colnames(RDI_RDD)[1]<-"ratio"
plot(x = log10(RDI_RDD$w_inf),y = log10(RDI_RDD$ratio))

saveRDS(sim_loop,"sim_rcp85_nolife.RDS")

## Plotting the results
rm(list=ls())

# On mac
setwd("~/Google Drive/hybrid-SDM/")

# On Windows
setwd("~/My Drive/hybrid-SDM")

# Original model
sim_rcp45<-readRDS("sim_rcp45.RDS")
sim_rcp85<-readRDS("sim_rcp85.RDS")

# Life-history only
sim_rcp45_nointer<-readRDS("sim_rcp45_nointer.RDS")
sim_rcp85_nointer<-readRDS("sim_rcp85_nointer.RDS")

# Interaction changes only
sim_rcp45_nolife<-readRDS("sim_rcp45_nolife.RDS")
sim_rcp85_nolife<-readRDS("sim_rcp85_nolife.RDS")

# Plotting for RCP45
## Size spectra from original model 
plotdata<-plotSpectra(sim_rcp45,total = TRUE)
df1<-data.frame(cbind(plotdata$data$w,plotdata$data$value))
colnames(df1)<-c("w","value")
df2<-aggregate(df1$value~df1$w,FUN = sum)

# Size spectra from static interaction model
plotdata2<-plotSpectra(sim_rcp45_nointer,total=TRUE)
df1_2<-data.frame(cbind(plotdata2$data$w,plotdata2$data$value))
colnames(df1_2)<-c("w","value")
df2_2<-aggregate(df1_2$value~df1_2$w,FUN = sum)

# Size spectra from static life history model
plotdata3<-plotSpectra(sim_rcp45_nolife,total = TRUE)
df1_3<-data.frame(cbind(plotdata3$data$w,plotdata3$data$value))
colnames(df1_3)<-c("w","value")
df2_3<-aggregate(df1_3$value~df1_3$w,FUN = sum)

# Lifehis + interactions (black)
# Static interaction (red)
# Static life history (blue)
plot<-ggplot(data = df2)+
  theme_bw()+
  theme(text = element_text(size =20))+
  scale_x_log10(limits=c(1,NA))+
  scale_y_log10()+
  labs(x = "Size (g)",y = "Biomass density")+
  geom_line(aes(x = `df1$w`,y = `df1$value`),size = 1)+
  geom_line(data = df2_2,aes(x = `df1_2$w`,y = `df1_2$value`),color = "red",linetype =2, size = 1)+
  geom_line(data = df2_3,aes(x = `df1_3$w`,y = `df1_3$value`),color = "blue",linetype = 4, size = 1)

plot

# Plotting for RCP85
## Size spectra from original model 
plotdata<-plotSpectra(sim_rcp85,total = TRUE)
df1<-data.frame(cbind(plotdata$data$w,plotdata$data$value))
colnames(df1)<-c("w","value")
df2<-aggregate(df1$value~df1$w,FUN = sum)

# Size spectra from static interaction model
plotdata2<-plotSpectra(sim_rcp85_nointer,total=TRUE)
df1_2<-data.frame(cbind(plotdata2$data$w,plotdata2$data$value))
colnames(df1_2)<-c("w","value")
df2_2<-aggregate(df1_2$value~df1_2$w,FUN = sum)

# Size spectra from static life history model
plotdata3<-plotSpectra(sim_rcp85_nolife,total = TRUE)
df1_3<-data.frame(cbind(plotdata3$data$w,plotdata3$data$value))
colnames(df1_3)<-c("w","value")
df2_3<-aggregate(df1_3$value~df1_3$w,FUN = sum)

# Plot for RCP85
# Lifehis + interactions (black)
# Static interaction (red)
# Static life history (blue)
plot<-ggplot(data = df2)+
  theme_bw()+
  theme(text = element_text(size =20))+
  scale_x_log10(limits=c(1,NA))+
  scale_y_log10()+
  labs(x = "Size (g)",y = "Biomass density")+
  geom_line(aes(x = `df1$w`,y = `df1$value`),size = 1)+
  geom_line(data = df2_2,aes(x = `df1_2$w`,y = `df1_2$value`),color = "red",linetype =2, size = 1)+
  geom_line(data = df2_3,aes(x = `df1_3$w`,y = `df1_3$value`),color = "blue",linetype = 4, size = 1)

plot

# Plotting species size spectra for RCP45 and 85, comparing results when
# life-history changes and distribution shifts are/are not taken into account

rm(list=ls())

setwd("~/Google Drive/hybrid-SDM")
source("mizer_misc_functions.R")
library(gridExtra)

species_param<-read_csv("nefsc_species_params_scenario3.csv")
species_param<-species_param[order(species_param$w_inf),]

sim_rcp45<-readRDS("sim_rcp45.RDS")
sim_rcp45_nolife<-readRDS("sim_rcp45_nolife.RDS")
sim_rcp45_nointer<-readRDS("sim_rcp45_nointer.RDS")

sim_rcp85<-readRDS("sim_rcp85.RDS")
sim_rcp85_nolife<-readRDS("sim_rcp85_nolife.RDS")
sim_rcp85_nointer<-readRDS("sim_rcp85_nointer.RDS")

df_rcp45<-plotSpectra(sim_rcp45)$data
df_rcp45_nolife<-plotSpectra(sim_rcp45_nolife)$data
df_rcp45_nointer<-plotSpectra(sim_rcp45_nointer)$data

df_rcp85<-plotSpectra(sim_rcp85)$data
df_rcp85_nolife<-plotSpectra(sim_rcp85_nolife)$data
df_rcp85_nointer<-plotSpectra(sim_rcp85_nointer)$data

nf_rcp45<-rbind(df_rcp45,df_rcp45_nolife,df_rcp45_nointer)
nf_rcp45$scenario<-rep(c("normal","nolife","nointer"),times=c(1635,1667,1633))

nf_rcp85<-rbind(df_rcp85,df_rcp85_nolife,df_rcp85_nointer)
nf_rcp85$scenario<-rep(c("normal","nolife","nointer"),times=c(1631,1666,1631))


nf<-nf_rcp45
# or
nf<-nf_rcp85

plot_sizeSpectra<-function(focal_sp){
  
  df<-nf[nf$Species==focal_sp,]
  
  biomass_normal<-sum(df[df$scenario=="normal",]$value)
  biomass_nolife<-sum(df[df$scenario=="nolife",]$value)
  biomass_nointer<-sum(df[df$scenario=="nointer",]$value)
  
  ggplot(df)+
    theme_bw()+
    theme(axis.text.x = element_text(face="bold",size=10),
          axis.text.y = element_text(face="bold",size=10),
          legend.position = "none")+
    labs(x="",y="")+
    geom_line(aes(x = w,y = value,color=scenario),size=1,linetype=1)+
    scale_color_manual(values = c("black","dark gray","light gray"))+
    scale_x_log10()+
    scale_y_log10(limits = c(NA,max(biomass_normal,biomass_nolife,biomass_nointer)))
  
}

# plot_sizeSpectra("Clupea harengus")

focal_sp<-species_param$species

p<-list()

for (i in 1:19){
  
  p[[i]]<-plot_sizeSpectra(focal_sp[i])
  
}

do.call(grid.arrange,p)
