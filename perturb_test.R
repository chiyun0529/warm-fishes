rm(list=ls())

setwd("~/Google Drive/hybrid-SDM")
setwd("G:/My Drive/hybrid-SDM")
source("mizer_misc_functions.R")

sim_2013<-readRDS("sim_2013.RDS")
params_2013<-sim_2013@params

initial_n_new_2013<-cbind(sim_2013@n[101,,1:54],sim_2013@n[101,,55:100]*0.5)

sim_2013_perturb<-project(params_2013,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                           initial_n = initial_n_new_2013,initial_n_pp = sim_2013@n_pp[1,])

perturb_2013<-getBiomass(sim_2013_perturb)
steady_2013<-getBiomass(sim_2013)

total_biomass_perturb_2013<-data.frame(apply(perturb_2013,MARGIN = 1,FUN = sum))
total_biomass_steady_2013<-sum(steady_2013[101,])

relative_biomass_2013<-total_biomass_perturb_2013/total_biomass_steady_2013
colnames(relative_biomass_2013)<-"relative_biomass"

sim_rcp45<-readRDS("sim_rcp45.RDS")
sim_rcp85<-readRDS("sim_rcp85.RDS")

params_rcp45<-sim_rcp45@params
params_rcp85<-sim_rcp85@params

initial_n_new_rcp45<-cbind(sim_rcp45@n[101,,1:53],sim_rcp45@n[101,,54:100]*0.5)

sim_rcp45_perturb<-project(params_rcp45,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                          initial_n = initial_n_new_rcp45,initial_n_pp = sim_rcp45@n_pp[1,])

perturb_rcp45<-getBiomass(sim_rcp45_perturb)
steady_rcp45<-getBiomass(sim_rcp45)

total_biomass_perturb_rcp45<-data.frame(apply(perturb_rcp45,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp45<-sum(steady_rcp45[101,])

relative_biomass_rcp45<-total_biomass_perturb_rcp45/total_biomass_steady_rcp45
colnames(relative_biomass_rcp45)<-"relative_biomass"

initial_n_new_rcp85<-cbind(sim_rcp85@n[101,,1:53],sim_rcp85@n[101,,54:100]*0.5)

sim_rcp85_perturb<-project(params_rcp85,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                           initial_n = initial_n_new_rcp85,initial_n_pp = sim_rcp85@n_pp[1,])

perturb_rcp85<-getBiomass(sim_rcp85_perturb)
steady_rcp85<-getBiomass(sim_rcp85)

total_biomass_perturb_rcp85<-data.frame(apply(perturb_rcp85,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp85<-sum(steady_rcp85[101,])

relative_biomass_rcp85<-total_biomass_perturb_rcp85/total_biomass_steady_rcp85
colnames(relative_biomass_rcp85)<-"relative_biomass"

## Estimate resilience if either interaction or life-history changes was not considered
sim_rcp45_nolife<-readRDS("sim_rcp45_nolife.RDS")
sim_rcp45_nointer<-readRDS("sim_rcp45_nointer.RDS")
sim_rcp85_nolife<-readRDS("sim_rcp85_nolife.RDS")
sim_rcp85_nointer<-readRDS("sim_rcp85_nointer.RDS")

params_rcp45_nolife<-sim_rcp45_nolife@params
params_rcp45_nointer<-sim_rcp45_nointer@params

params_rcp85_nolife<-sim_rcp85_nolife@params
params_rcp85_nointer<-sim_rcp85_nointer@params

initial_n_new_rcp45_nolife<-cbind(sim_rcp45_nolife@n[101,,1:54],sim_rcp45_nolife@n[101,,55:100]*0.5)

sim_rcp45_nolife_perturb<-project(params_rcp45_nolife,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                           initial_n = initial_n_new_rcp45_nolife,initial_n_pp = sim_rcp45_nolife@n_pp[1,])

plotlyBiomass(sim_rcp45_nolife_perturb)


initial_n_new_rcp45_nointer<-cbind(sim_rcp45_nointer@n[101,,1:53],sim_rcp45_nointer@n[101,,54:100]*0.5)


sim_rcp45_nointer_perturb<-project(params_rcp45_nointer,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                                  initial_n = initial_n_new_rcp45_nointer,initial_n_pp = sim_rcp45_nointer@n_pp[1,])

plotlyBiomass(sim_rcp45_nointer_perturb)


initial_n_new_rcp85_nolife<-cbind(sim_rcp85_nolife@n[101,,1:54],sim_rcp85_nolife@n[101,,55:100]*0.5)

sim_rcp85_nolife_perturb<-project(params_rcp85_nolife,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                                  initial_n = initial_n_new_rcp85_nolife,initial_n_pp = sim_rcp85_nolife@n_pp[1,])

plotlyBiomass(sim_rcp85_nolife_perturb)


initial_n_new_rcp85_nointer<-cbind(sim_rcp85_nointer@n[101,,1:53],sim_rcp85_nointer@n[101,,54:100]*0.5)


sim_rcp85_nointer_perturb<-project(params_rcp85_nointer,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                                   initial_n = initial_n_new_rcp85_nointer,initial_n_pp = sim_rcp85_nointer@n_pp[1,])

plotlyBiomass(sim_rcp85_nointer_perturb)

# RCP45_nolife
perturb_rcp45_nolife<-getBiomass(sim_rcp45_nolife_perturb)
steady_rcp45_nolife<-getBiomass(sim_rcp45_nolife)

total_biomass_perturb_rcp45_nolife<-data.frame(apply(perturb_rcp45_nolife,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp45_nolife<-sum(steady_rcp45_nolife[101,])

relative_biomass_rcp45_nolife<-total_biomass_perturb_rcp45_nolife/total_biomass_steady_rcp45_nolife
colnames(relative_biomass_rcp45_nolife)<-"relative_biomass"

# RCP45_nointer
perturb_rcp45_nointer<-getBiomass(sim_rcp45_nointer_perturb)
steady_rcp45_nointer<-getBiomass(sim_rcp45_nointer)

total_biomass_perturb_rcp45_nointer<-data.frame(apply(perturb_rcp45_nointer,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp45_nointer<-sum(steady_rcp45_nointer[101,])

relative_biomass_rcp45_nointer<-total_biomass_perturb_rcp45_nointer/total_biomass_steady_rcp45_nointer
colnames(relative_biomass_rcp45_nointer)<-"relative_biomass"

# RCP85_nolife
perturb_rcp85_nolife<-getBiomass(sim_rcp85_nolife_perturb)
steady_rcp85_nolife<-getBiomass(sim_rcp85_nolife)

total_biomass_perturb_rcp85_nolife<-data.frame(apply(perturb_rcp85_nolife,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp85_nolife<-sum(steady_rcp85_nolife[101,])

relative_biomass_rcp85_nolife<-total_biomass_perturb_rcp85_nolife/total_biomass_steady_rcp85_nolife
colnames(relative_biomass_rcp85_nolife)<-"relative_biomass"

# RCP85_nointer
perturb_rcp85_nointer<-getBiomass(sim_rcp85_nointer_perturb)
steady_rcp85_nointer<-getBiomass(sim_rcp85_nointer)

total_biomass_perturb_rcp85_nointer<-data.frame(apply(perturb_rcp85_nointer,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp85_nointer<-sum(steady_rcp85_nointer[101,])

relative_biomass_rcp85_nointer<-total_biomass_perturb_rcp85_nointer/total_biomass_steady_rcp85_nointer
colnames(relative_biomass_rcp85_nointer)<-"relative_biomass"

df<-data.frame(c(0:1000),relative_biomass_2013,relative_biomass_rcp45,relative_biomass_rcp85,
               relative_biomass_rcp45_nointer,relative_biomass_rcp45_nolife,
               relative_biomass_rcp85_nointer,relative_biomass_rcp85_nolife)

colnames(df)<-c("time","rcp2013","rcp45","rcp85")

df$time<-c(seq(from = 0.1, to = 100.1, by = 0.1))

df$res_2013<-1/sum(df$rcp2013^2*0.1)
df$res_rcp45<-1/sum(df$rcp45^2*0.1)
df$res_rcp85<-1/sum(df$rcp85^2*0.1)
df$res_rcp45_nointer<-1/sum(df$rcp45_nointer^2*0.1)
df$res_rcp45_nolife<-1/sum(df$rcp45_nolife^2*0.1)
df$res_rcp85_nointer<-1/sum(df$rcp85_nointer^2*0.1)
df$res_rcp85_nolife<-1/sum(df$rcp85_nolife^2*0.1)

write_csv(df,"resilience_top-down.csv")

## Testing for stabilization
rm(list=ls())

resilience<-read_csv("resilience_top-down.csv")

library(tseries)

# Augmented Dickey-Fuller test to see if the time series has stabilized 
adf.test(resilience$rcp2013,alternative = "e")
adf.test(resilience$rcp45,alternative = "e")
adf.test(resilience$rcp85,alternative = "e")
adf.test(resilience$rcp45_nointer,alternative = "e")
adf.test(resilience$rcp45_nolife,alternative = "e")
adf.test(resilience$rcp85_nointer,alternative = "e")
adf.test(resilience$rcp85_nolife,alternative = "e")
## The results show that they all have by the end of the time series

## Bottom-up perturbations
rm(list=ls())

setwd("~/Google Drive/hybrid-SDM")

setwd("G:/My Drive/hybrid-SDM")

source("mizer_misc_functions.R")

# 2013
sim_2013<-readRDS("sim_2013.RDS")
params_2013<-sim_2013@params

initial_n_pp_new_2013<-c(sim_2013@n_pp[101,1:119]*1.5,sim_2013@n_pp[101,120:220])

sim_2013_perturb<-project(params_2013,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                          initial_n = sim_2013@n[101,,],initial_n_pp = initial_n_pp_new_2013)

plotlyBiomass(sim_2013_perturb)
plotlyBiomass(sim_2013)

perturb_2013<-getBiomass(sim_2013_perturb)
steady_2013<-getBiomass(sim_2013)

total_biomass_perturb_2013<-data.frame(apply(perturb_2013,MARGIN = 1,FUN = sum))
total_biomass_steady_2013<-sum(steady_2013[101,])

relative_biomass_2013<-total_biomass_perturb_2013/total_biomass_steady_2013
colnames(relative_biomass_2013)<-"relative_biomass"

# RCP45
sim_rcp45<-readRDS("sim_rcp45.RDS")
params_rcp45<-sim_rcp45@params

initial_n_pp_new_rcp45<-c(sim_rcp45@n_pp[101,1:119]*1.5,sim_rcp45@n_pp[101,120:219])

sim_rcp45_perturb<-project(params_rcp45,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                           initial_n = sim_rcp45@n[101,,],initial_n_pp = initial_n_pp_new_rcp45)

plotlyBiomass(sim_rcp45_perturb)
plotlyBiomass(sim_rcp45)

perturb_rcp45<-getBiomass(sim_rcp45_perturb)
steady_rcp45<-getBiomass(sim_rcp45)

total_biomass_perturb_rcp45<-data.frame(apply(perturb_rcp45,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp45<-sum(steady_rcp45[101,])

relative_biomass_rcp45<-total_biomass_perturb_rcp45/total_biomass_steady_rcp45
colnames(relative_biomass_rcp45)<-"relative_biomass"

# RCP85
sim_rcp85<-readRDS("sim_rcp85.RDS")
params_rcp85<-sim_rcp85@params

initial_n_pp_new_rcp85<-c(sim_rcp85@n_pp[101,1:118]*1.5,sim_rcp85@n_pp[101,119:218])

sim_rcp85_perturb<-project(params_rcp85,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                           initial_n = sim_rcp85@n[101,,],initial_n_pp = initial_n_pp_new_rcp85)

plotlyBiomass(sim_rcp85_perturb)
plotlyBiomass(sim_rcp85)

perturb_rcp85<-getBiomass(sim_rcp85_perturb)
steady_rcp85<-getBiomass(sim_rcp85)

total_biomass_perturb_rcp85<-data.frame(apply(perturb_rcp85,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp85<-sum(steady_rcp85[101,])

relative_biomass_rcp85<-total_biomass_perturb_rcp85/total_biomass_steady_rcp85
colnames(relative_biomass_rcp85)<-"relative_biomass"

## Estimate resilience if either interaction or life-history changes was not considered
sim_rcp45_nolife<-readRDS("sim_rcp45_nolife.RDS")
sim_rcp45_nointer<-readRDS("sim_rcp45_nointer.RDS")
sim_rcp85_nolife<-readRDS("sim_rcp85_nolife.RDS")
sim_rcp85_nointer<-readRDS("sim_rcp85_nointer.RDS")

# RCP45_nolife
params_rcp45_nolife<-sim_rcp45_nolife@params

initial_n_pp_new_rcp45nolife<-c(sim_rcp45_nolife@n_pp[101,1:120]*1.5,sim_rcp45_nolife@n_pp[101,121:220])

sim_rcp45_nolife_perturb<-project(params_rcp45_nolife,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                           initial_n = sim_rcp45_nolife@n[101,,],initial_n_pp = initial_n_pp_new_rcp45nolife)

plotlyBiomass(sim_rcp45_nolife_perturb)

perturb_rcp45_nolife<-getBiomass(sim_rcp45_nolife_perturb)
steady_rcp45_nolife<-getBiomass(sim_rcp45_nolife)

total_biomass_perturb_rcp45_nolife<-data.frame(apply(perturb_rcp45_nolife,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp45_nolife<-sum(steady_rcp45_nolife[101,])

relative_biomass_rcp45_nolife<-total_biomass_perturb_rcp45_nolife/total_biomass_steady_rcp45_nolife
colnames(relative_biomass_rcp45_nolife)<-"relative_biomass"

# RCP45_nointer
params_rcp45_nointer<-sim_rcp45_nointer@params

initial_n_pp_new_rcp45nointer<-c(sim_rcp45_nointer@n_pp[101,1:119]*1.5,sim_rcp45_nointer@n_pp[101,120:219])

sim_rcp45_nointer_perturb<-project(params_rcp45_nointer,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                                  initial_n = sim_rcp45_nointer@n[101,,],initial_n_pp = initial_n_pp_new_rcp45nointer)


plotlyBiomass(sim_rcp45_nointer_perturb)

perturb_rcp45_nointer<-getBiomass(sim_rcp45_nointer_perturb)
steady_rcp45_nointer<-getBiomass(sim_rcp45_nointer)

total_biomass_perturb_rcp45_nointer<-data.frame(apply(perturb_rcp45_nointer,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp45_nointer<-sum(steady_rcp45_nointer[101,])

relative_biomass_rcp45_nointer<-total_biomass_perturb_rcp45_nointer/total_biomass_steady_rcp45_nointer
colnames(relative_biomass_rcp45_nointer)<-"relative_biomass"

# RCP85_nolife
params_rcp85_nolife<-sim_rcp85_nolife@params

initial_n_pp_new_rcp85nolife<-c(sim_rcp85_nolife@n_pp[101,1:120]*1.5,sim_rcp85_nolife@n_pp[101,121:220])

sim_rcp85_nolife_perturb<-project(params_rcp85_nolife,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                                  initial_n = sim_rcp85_nolife@n[101,,],initial_n_pp = initial_n_pp_new_rcp85nolife)

plotlyBiomass(sim_rcp85_nolife_perturb)

perturb_rcp85_nolife<-getBiomass(sim_rcp85_nolife_perturb)
steady_rcp85_nolife<-getBiomass(sim_rcp85_nolife)

total_biomass_perturb_rcp85_nolife<-data.frame(apply(perturb_rcp85_nolife,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp85_nolife<-sum(steady_rcp85_nolife[101,])

relative_biomass_rcp85_nolife<-total_biomass_perturb_rcp85_nolife/total_biomass_steady_rcp85_nolife
colnames(relative_biomass_rcp85_nolife)<-"relative_biomass"

# RCP85_nointer
params_rcp85_nointer<-sim_rcp85_nointer@params

initial_n_pp_new_rcp85nointer<-c(sim_rcp85_nointer@n_pp[101,1:118]*1.5,sim_rcp85_nointer@n_pp[101,119:218])

sim_rcp85_nointer_perturb<-project(params_rcp85_nointer,effort = 0.2,dt = 0.1,t_save = 0.1,t_max = 100,
                                   initial_n = sim_rcp85_nointer@n[101,,],initial_n_pp = initial_n_pp_new_rcp85nointer)

plotlyBiomass(sim_rcp85_nointer_perturb)

perturb_rcp85_nointer<-getBiomass(sim_rcp85_nointer_perturb)
steady_rcp85_nointer<-getBiomass(sim_rcp85_nointer)

total_biomass_perturb_rcp85_nointer<-data.frame(apply(perturb_rcp85_nointer,MARGIN = 1,FUN = sum))
total_biomass_steady_rcp85_nointer<-sum(steady_rcp85_nointer[101,])

relative_biomass_rcp85_nointer<-total_biomass_perturb_rcp85_nointer/total_biomass_steady_rcp85_nointer
colnames(relative_biomass_rcp85_nointer)<-"relative_biomass"

df<-data.frame(c(0:1000),relative_biomass_2013,relative_biomass_rcp45,relative_biomass_rcp85,
               relative_biomass_rcp45_nointer,relative_biomass_rcp45_nolife,
               relative_biomass_rcp85_nointer,relative_biomass_rcp85_nolife)

colnames(df)<-c("time","rcp2013","rcp45","rcp85","rcp45_nointer","rcp45_nolife",
                "rcp85_nointer","rcp85_nolife")

df$time<-rownames(df)

df$res_2013<-1/sum(df$rcp2013^2*0.1)
df$res_rcp45<-1/sum(df$rcp45^2*0.1)
df$res_rcp85<-1/sum(df$rcp85^2*0.1)
df$res_rcp45_nointer<-1/sum(df$rcp45_nointer^2*0.1)
df$res_rcp45_nolife<-1/sum(df$rcp45_nolife^2*0.1)
df$res_rcp85_nointer<-1/sum(df$rcp85_nointer^2*0.1)
df$res_rcp85_nolife<-1/sum(df$rcp85_nolife^2*0.1)

write_csv(df,"resilience_bottom-up.csv")

rm(list=ls())

resilience<-read_csv("resilience_bottom-up.csv")

# Augmented Dickey-Fuller test to see if the time series has stabilized 
adf.test(resilience$rcp2013,alternative = "e")
adf.test(resilience$rcp45,alternative = "e")
adf.test(resilience$rcp85,alternative = "e")
adf.test(resilience$rcp45_nointer,alternative = "e")
adf.test(resilience$rcp45_nolife,alternative = "e")
adf.test(resilience$rcp85_nointer,alternative = "e")
adf.test(resilience$rcp85_nolife,alternative = "e")
## The results show that they all have by the end of the time series
