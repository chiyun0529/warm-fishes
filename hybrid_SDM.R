rm(list=ls())

# If run on Mac
setwd("~/Google Drive/hybrid-SDM")

# If run on Windows
setwd("G:/?ڪ????ݵw??/hybrid-SDM")

source("Bayesian_SDM.r")

nefsc<-read_csv("nefsc_SDM.csv")
speciesdata<-read_csv("speciesdata_SDM.csv")
world <- ne_countries(scale = "medium", returnclass = "sf")

# When run on mac
gm_auth_configure(path = "~/Google Drive/GmailCredentials.json")

# When run on windows
gm_auth_configure(path = "G:/?ڪ????ݵw??/GmailCredentials.json")

sender<-"kuochiyun@gmail.com"

email <- gm_mime() %>%
  gm_to(sender) %>%
  gm_from(sender) %>%
  gm_subject("Analysis complete") %>%
  gm_text_body("as title")


################################
# 1. Cross-validation

for (i in 1:nrow(speciesdata)){

species<-as.character(speciesdata[i,1])
dir.create(paste0("~/Desktop/hybrid- SDM/",species))

source("hybrid_SDM_globals.r")

SDM_prelim(species = species,presence_ratio = as.numeric(speciesdata[i,4]),ctmax = as.numeric(speciesdata[i,5]))

}

################################
# 2. SDM using all data after cross-validation
rm(list=ls())

nefsc<-read_csv("~/Desktop/hybrid-SDM/nefsc_iap_t.csv")
speciesdata<-read_csv("~/Desktop/hybrid-SDM/speciesdata.csv")
nefsc_SDM<-read_csv("~/Desktop/hybrid-SDM/nefsc_SDM.csv")
unique_locations<-data.frame(unique(nefsc_SDM$grid_id))

temps<-data.frame(cbind(unique_locations,matrix(0,nrow = nrow(unique_locations),ncol=9)))

colnames(temps)<-c("grid_id","2013_s","2013_b","2013_i","rcp45_s","rcp45_b","rcp45_i","rcp85_s","rcp85_b","rcp85_i")

for (i in 1:nrow(temps)){
  
  temps[i,2]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$sst_2013[1]
  temps[i,3]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$sbt_2013[1]
  temps[i,4]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$sit_2013[1]
  
  temps[i,5]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$RCP45SMean[1]
  temps[i,6]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$RCP45BMean[1]
  temps[i,7]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$RCP45IMean[1]
  
  temps[i,8]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$RCP85SMean[1]
  temps[i,9]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$RCP85BMean[1]
  temps[i,10]<-nefsc_SDM[nefsc_SDM$grid_id==temps[i,1],]$RCP85IMean[1]
}

write_csv(temps,"~/Desktop/hybrid-SDM/temps_SDM.csv")

# Create presence-absence data for SDM
for (i in 1:nrow(speciesdata)){
  
  species<-as.character(speciesdata[i,2])
  presence_ratio<-as.numeric(speciesdata[i,5])
  setwd(paste0("~/Google Drive/hybrid-SDM/",species))
  
  nefsc$presence<-ifelse(nefsc$scientificName == species,1,0)
  
  presence<-subset(nefsc,subset = presence ==1)
  absence<-subset(nefsc,subset = presence ==0)
  
  true_absence<-anti_join(x = absence,y = presence,by = c("grid_lon","grid_lat","depth0"))
  
  if (presence_ratio>=0.2){
    
    occurrence<-rbind(presence,true_absence)
  } else {
    true_absence_rows<-sample(c(1:nrow(true_absence)),size = nrow(presence)*5,replace = FALSE)
    true_absence_sub<-true_absence[true_absence_rows,]
    
    occurrence<-rbind(presence,true_absence_sub)
    
  }
  
  write_csv(occurrence,paste0("presence_absence_",species,".csv"))
  
}

temps<-read_csv("~/Desktop/hybrid-SDM/temps_SDM.csv")
source("~/Desktop/hybrid-SDM/hybrid_SDM_globals.r")
source("~/Desktop/hybrid-SDM/Bayesian_SDM.r")

settings <- list(
  chains=4, 
  tuning = 1000, 
  burnin = 4000, 
  samples=4000, 
  thin=1, 
  seed=NA,
  diagnostics=TRUE,
  quantiles=c(0.05, 0.95))


for (i in 1:nrow(speciesdata)){
  
  SDM(species = as.character(speciesdata[i,2]),habitat = as.character(speciesdata[i,3]),SDM_type = as.character(speciesdata[i,8]))
  
}

################################

# 3.1 Estimate the area of suitable habitats by counting "suitable" grids

rm(list = ls())

speciesdata<-read_csv("speciesdata_SDM.csv")

speciesdata$suitable_2013<--NA
speciesdata$suitable_rcp45<--NA
speciesdata$suitable_rcp85<--NA

for (i in 1:nrow(speciesdata)){
  
  species<-speciesdata[i,]$species
  
  if (is.na(speciesdata[i,]$SDM_type)==FALSE){
  
  k_max<-speciesdata[i,]$k_max
    
  load(paste0("~/SDM results/",species,"/",species,"_SDM.rdata"))
  
  if (speciesdata[i,]$SDM_type=="N"){  
    
  m1Predictions_2013$tempPredict$presence<-ifelse(m1Predictions_2013$tempPredict$mean>=k_max,1,0)
  m1Predictions_rcp45$tempPredict$presence<-ifelse(m1Predictions_rcp45$tempPredict$mean>=k_max,1,0)
  m1Predictions_rcp85$tempPredict$presence<-ifelse(m1Predictions_rcp85$tempPredict$mean>=k_max,1,0)
  
  suitable_2013<-sum(m1Predictions_2013$tempPredict$presence)
  suitable_rcp45<-sum(m1Predictions_rcp45$tempPredict$presence)
  suitable_rcp85<-sum(m1Predictions_rcp85$tempPredict$presence)
  
  } else if (speciesdata[i,]$SDM_type=="H"){
  
  mmPredictions_2013$tempPredict$presence<-ifelse(mmPredictions_2013$tempPredict$mean>=k_max,1,0)
  mmPredictions_rcp45$tempPredict$presence<-ifelse(mmPredictions_rcp45$tempPredict$mean>=k_max,1,0)
  mmPredictions_rcp85$tempPredict$presence<-ifelse(mmPredictions_rcp85$tempPredict$mean>=k_max,1,0)
    
  suitable_2013<-sum(mmPredictions_2013$tempPredict$presence)
  suitable_rcp45<-sum(mmPredictions_rcp45$tempPredict$presence)
  suitable_rcp85<-sum(mmPredictions_rcp85$tempPredict$presence)   
      
  }
    
  speciesdata[i,]$suitable_2013<-suitable_2013
  speciesdata[i,]$suitable_rcp45<-suitable_rcp45
  speciesdata[i,]$suitable_rcp85<-suitable_rcp85
  }
}

temps<-read_csv("grid_temps.csv")

temps$grid_lat<-0 
temps$grid_lon<-0

for (i in 1:nrow(temps)){
  
  temps[i,]$grid_lat<-nefsc[nefsc$grid_id==temps[i,]$grid_id,]$grid_lat[1]
  temps[i,]$grid_lon<-nefsc[nefsc$grid_id==temps[i,]$grid_id,]$grid_lon[1]
  
}

# 3.2 Constructing interaction matrix for mizer based on SDM results

rm(list=ls())

temps<-read_csv("grid_temps.csv")

speciesdata<-read_csv("speciesdata_SDM.csv")

species_list<-speciesdata[is.na(speciesdata$k_max)==FALSE,]$species

# Removing "Hippoglossoides platessoides", "Pollachius virens", and "Sebastes fasciatus"  
species_list<-species_list[-c(4,7,10,11,18,21)]

occurrence_2013<-data.frame(matrix(0,nrow = 677,ncol = 25))

rownames(occurrence_2013)<-temps$grid_id

occurrence_rcp45<-occurrence_rcp85<-occurrence_2013

species_list<-speciesdata[is.na(speciesdata$k_max)==FALSE,]$species

for (i in 1:length(species_list)){
  
  species<-species_list[i]  
  
  colnames(occurrence_2013)[i]<-species
  colnames(occurrence_rcp45)[i]<-species
  colnames(occurrence_rcp85)[i]<-species
  
  k_max<-speciesdata[speciesdata$species==species,]$k_max
  
  
  load(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/",species,"_SDM.rdata"))
  
  # On Windows
  # load("G:/?ڪ????ݵw??/hybrid-SDM/SDM results/Clupea harengus/Clupea harengus_SDM.rdata")
  
  if (speciesdata[speciesdata$species==species,]$SDM_type=="H"){
    
    mmPredictions_2013$tempPredict$presence<-ifelse(mmPredictions_2013$tempPredict$mean>=k_max,1,0)
    mmPredictions_rcp45$tempPredict$presence<-ifelse(mmPredictions_rcp45$tempPredict$mean>=k_max,1,0)
    mmPredictions_rcp85$tempPredict$presence<-ifelse(mmPredictions_rcp85$tempPredict$mean>=k_max,1,0)
    
    occurrence_2013[,i]<-mmPredictions_2013$tempPredict$presence
    occurrence_rcp45[,i]<-mmPredictions_rcp45$tempPredict$presence
    occurrence_rcp85[,i]<-mmPredictions_rcp85$tempPredict$presence
    
  } else if (speciesdata[speciesdata$species==species,]$SDM_type=="N"){
    
    m1Predictions_2013$tempPredict$presence<-ifelse(m1Predictions_2013$tempPredict$mean>=k_max,1,0)
    m1Predictions_rcp45$tempPredict$presence<-ifelse(m1Predictions_rcp45$tempPredict$mean>=k_max,1,0)
    m1Predictions_rcp85$tempPredict$presence<-ifelse(m1Predictions_rcp85$tempPredict$mean>=k_max,1,0)  
    
    occurrence_2013[,i]<-m1Predictions_2013$tempPredict$presence
    occurrence_rcp45[,i]<-m1Predictions_rcp45$tempPredict$presence
    occurrence_rcp85[,i]<-m1Predictions_rcp85$tempPredict$presence
    
  }
}

write_csv(occurrence_2013,"occurrence_2013.csv")
write_csv(occurrence_rcp45,"occurrence_rcp45.csv")
write_csv(occurrence_rcp85,"occurrence_rcp85.csv")


# For each species, generate a matrix of degrees of spatial overlap
# with the other species, taking into account predicted presence from SDM
# and habitat preference.

# Read in predicted presence of each species in each grid

climate<-c("2013","rcp45","rcp85")

scenario<-3

# Choose the scenario
if (scenario==1){

occurrence<-read_csv("occurrence_2013.csv")

} else if (scenario==2){
  
occurrence<-read_csv("occurrence_rcp45.csv")

} else if (scenario==3){

occurrence<-read_csv("occurrence_rcp85.csv")

}

occurrence<-occurrence[,-c(4,7,10,11,18,21)]

# Calculate the degree of spatial overlap for each species pair
iteration<-length(species_list)

for (i in 1:(iteration-1)){
  
  df<-matrix(-99,nrow = 677,ncol = iteration-i)
  
  j<-1
  
  while (i+j<=iteration){
    
    species1<-species_list[i]
    species2<-species_list[i+j]
    
    habitat1<-speciesdata[speciesdata$species==species1,]$habitat
    habitat2<-speciesdata[speciesdata$species==species2,]$habitat
    
    
    if (habitat1==habitat2){
      
      h<-1
      
    } else if (habitat1!=habitat2) {
      
      if (habitat1!="benthopelagic"& habitat2!="benthopelagic"){
        
        h<-0.5
        # h<-0
        
      } else {
        
        h<-1
        # h<-0.5
      }
    }
    
    df2<-cbind(occurrence[,i],occurrence[,i+j])
    
    for (k in 1:677){
      
      if ((df2[k,1]==1) & (df2[k,2]==1)){
        
        df2[k,3]<-h
        
      } else {df2[k,3]<-0}
      
    } 
    
    df[,j]<-df2[,3]
    
    j<-j+1
    
  }
  
  df<-data.frame(df)
  colnames(df)<-species_list[(i+1):iteration]
  write_csv(df,paste0(species1,"_co_occur.csv"))
}

# Calculate overall spatial overlap for each species pair
# Setting up different climate scenarios
climate<-c("2013","rcp45","rcp85")

scenario<-3

# On Mac
dir<-paste0("~/Google Drive/hybrid-SDM/spatial overlap_",climate[scenario],"_3/")

# On Windows
dir<-paste0("G:/我的雲端硬碟/hybrid-SDM/spatial overlap_",climate[scenario],"_3/")


files<-list.files(dir)

# Setting up the interactions matrix

interactions<-matrix(-99,nrow = iteration,ncol = iteration)

for (i in 1:length(files)){
  
  # sp1<-species_list[i]
  # df_sp1<-occurrence[,i]

  interactions[i,i]<-1
  
  df<-read_csv(paste0(dir,files[i]))
  
  no_pairs<-ncol(df)
  
  for (j in 1:no_pairs){
  
  # Fill in the corresponding cells in the interaction matrix with overall overlap,
  # which is calculated as 
    # df_spj<-occurrence[,j+1]
    # df_join<-cbind(df_sp1,df_spj)  
    
    # double_absence<-nrow(df_join[which(df_join[,1]==0 & df_join[,2]==0),])
    
  interactions[i,min(j+i,iteration)]<-interactions[min(j+i,iteration),i]<-sum(df[,j])/677
  # interactions[i,min(j+i,iteration)]<-interactions[min(j+i,iteration),i]<-sum(df[,j])/(677-double_absence)
  
  }
  
}

interactions[iteration,iteration]<-1

interactions<-data.frame(interactions)

dimnames(interactions)[[1]]<-dimnames(interactions)[[2]]<-species_list

write_csv(interactions,paste0("interactions_",climate[scenario],".csv"))

################################
# 4 calculating mean temperatures of inhabited grids in 2013, rcp45, and rcp85 based on SDM results
rm(list=ls())
temps<-read_csv("grid_temps.csv")
speciesdata<-read_csv("speciesdata_SDM.csv")

climate<-c("2013","rcp45","rcp85")

scenario<-3

species_list<-speciesdata$species

df_temp<-data.frame(matrix(-99,nrow = 25,ncol = 2))

colnames(df_temp)<-c("species",paste0("temp_",climate[scenario]))

for (i in 1:length(species_list)){

# Change species names as necessary
species<-species_list[i]

df_temp[i,1]<-species_list[i]

# On Mac
load(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/",species,"_SDM.rdata"))

# On Windows
# load(paste0("G:/我的雲端硬碟/Google Drive/hybrid-SDM/SDM results/",species,"/",species,"_SDM.rdata"))

k_max<-speciesdata[speciesdata$species==species,]$k_max

SDM_type<-speciesdata[speciesdata$species==species,]$SDM_type

grid_coords<-cbind(temps$grid_lat,temps$grid_lon)

if (SDM_type=="H"){

  if (scenario==1){
  
    mmPredictions_2013$tempPredict$presence<-ifelse(mmPredictions_2013$tempPredict$mean>=k_max,1,0)
    grids_all<-cbind(mmPredictions_2013$tempPredict,grid_coords)
    
    grid_temp<-cbind(mmPredictions_2013$tempDomain,grid_coords)
    
  } else if (scenario==2){
  
    mmPredictions_rcp45$tempPredict$presence<-ifelse(mmPredictions_rcp45$tempPredict$mean>=k_max,1,0)
    grids_all<-cbind(mmPredictions_rcp45$tempPredict,grid_coords)
    
    grid_temp<-cbind(mmPredictions_rcp45$tempDomain,grid_coords)
  
  } else if (scenario==3){
  
    mmPredictions_rcp85$tempPredict$presence<-ifelse(mmPredictions_rcp85$tempPredict$mean>=k_max,1,0)
    grids_all<-cbind(mmPredictions_rcp85$tempPredict,grid_coords)
    
    grid_temp<-cbind(mmPredictions_rcp85$tempDomain,grid_coords)
  }


} else if (SDM_type=="N"){
  
m1Predictions_2013$tempPredict$presence<-ifelse(m1Predictions_2013$tempPredict$mean>=k_max,1,0)
m1Predictions_rcp45$tempPredict$presence<-ifelse(m1Predictions_rcp45$tempPredict$mean>=k_max,1,0)
m1Predictions_rcp85$tempPredict$presence<-ifelse(m1Predictions_rcp85$tempPredict$mean>=k_max,1,0)
  
  if (scenario==1){
  
  m1Predictions_2013$tempPredict$presence<-ifelse(m1Predictions_2013$tempPredict$mean>=k_max,1,0)
  grids_all<-cbind(m1Predictions_2013$tempPredict,grid_coords)
  
  grid_temp<-cbind(m1Predictions_2013$tempDomain,grid_coords)
  
} else if (scenario==2){
  
  m1Predictions_rcp45$tempPredict$presence<-ifelse(m1Predictions_rcp45$tempPredict$mean>=k_max,1,0)
  grids_all<-cbind(m1Predictions_rcp45$tempPredict,grid_coords)
  
  grid_temp<-cbind(m1Predictions_rcp45$tempDomain,grid_coords)
  
} else if (scenario==3){
  
  m1Predictions_rcp85$tempPredict$presence<-ifelse(m1Predictions_rcp85$tempPredict$mean>=k_max,1,0)
  grids_all<-cbind(m1Predictions_rcp85$tempPredict,grid_coords)
  
  grid_temp<-cbind(m1Predictions_rcp85$tempDomain,grid_coords)
} 
  
}

colnames(grid_temp)[c(2,3)]<-c("lat","lon")

presence<-subset(grids_all,subset=presence==1)
absence<-subset(grids_all,subset = presence==0)

colnames(presence)[6:7]<-c("lat","lon")
colnames(absence)[6:7]<-c("lat","lon")

# Calculate mean temperature of presence grid in 2013, RCP45, and RCP85

grid_id_presence<-rownames(presence)

temp_presence<-NULL

for (j in 1:length(grid_id_presence)){
  
  temp_presence[j]<-grid_temp[rownames(grid_temp)==grid_id_presence[j],]$temp
  mean_temp<-mean(temp_presence)
  
}

write_csv(presence,paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/presence_",climate[scenario],".csv"))
write_csv(absence,paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/absence_",climate[scenario],".csv"))

df_temp[i,2]<-mean_temp
}

################################
# 5 Plotting presence and absence
rm(list=ls())
speciesdata<-read_csv("speciesdata_SDM.csv")
species_list<-speciesdata$species

world <- ne_countries(scale = "medium", returnclass = "sf")

for (i in 1:length(species_list)){
  
  presence<-absence<-NULL
  
  # Change species names as necessary
  species<-species_list[i]
      
  presence_2013<-read_csv(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/presence_2013.csv"))
  presence_rcp45<-read_csv(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/presence_rcp45.csv"))
  
  # absence<-read_csv(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/absence_",climate[scenario],".csv"))
    
  png(paste0("~/Desktop/",species,".png"))
  
  plot<-ggplot(data = world)+
    theme_bw()+
    labs(x = "",y = "",title = species)+
    geom_sf()+
    coord_sf(xlim = c(-85,-60),ylim = c(25,50))+
    geom_point(data = presence,aes(x=lon,y = lat),size=0.2)+
    geom_point(data = absence,aes(x=lon,y = lat),size=0.2,color="red")
  
  print(plot)
  dev.off()
}

## Plotting presences in 2013, rcp45, rcp85
rm(list=ls())

setwd("~/Google Drive/hybrid-SDM/")

library(tidyverse)
library(ggConvexHull)
source("Bayesian_SDM.r")
source("geom_bag.r")

speciesdata<-read_csv("nefsc_species_params_scenario3.csv")

sp<-speciesdata$species

world <- ne_countries(scale = "medium", returnclass = "sf")

for (i in 1:19){
  
  species<-sp[i]
  
  presence_2013<-read_csv(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/presence_2013.csv"))
  presence_2013$scenario<-"rcp2013"
  
  presence_rcp45<-read_csv(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/presence_rcp45.csv"))
  presence_rcp45$scenario<-"rcp45"
  
  presence_rcp85<-read_csv(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/presence_rcp85.csv"))
  presence_rcp85$scenario<-"rcp85"
  
  df_all<-rbind(presence_2013,presence_rcp45,presence_rcp85)
  df<-df_all[df_all$lon!=-74.625 & df_all$lat!=28.875,]
  
  hull<-df %>% 
    group_by(scenario) %>%
    slice(chull(lon,lat))
  
  # absence<-read_csv(paste0("~/Google Drive/hybrid-SDM/SDM results/",species,"/absence_",climate[scenario],".csv"))
  
  pdf(paste0("~/Desktop/",species,".pdf"))
  
  plot<-ggplot(data = world)+
    theme_bw()+
    theme(legend.position = "none")+
    labs(x = "",y = "")+
    geom_sf()+
    coord_sf(xlim = c(-85,-60),ylim = c(25,50))+
    # geom_bag(data = hull,aes(x=lon,y = lat,fill = scenario),prop = 0.70,alpha = 0.3)+
    # scale_fill_manual(values = c("red","dark green","blue"))
    geom_point(data = df,aes(x=lon,y = lat,color = scenario,shape=scenario),size=2)+
    scale_shape_manual(values = c(1,4,15))+
    scale_color_manual(values = c("#FF000080","#00640080","#0000FF33"))
  
  
  print(plot)
  dev.off()
}



speciesdata$vonBert_2013<-speciesdata$vonBert_intercept+speciesdata$vonBert_slope*(speciesdata$temp_2013-speciesdata$mean_temp)
speciesdata$vonBert_RCP45<-speciesdata$vonBert_intercept+speciesdata$vonBert_slope*(speciesdata$temp_rcp45-speciesdata$mean_temp)
speciesdata$vonBert_RCP85<-speciesdata$vonBert_intercept+speciesdata$vonBert_slope*(speciesdata$temp_rcp85-speciesdata$mean_temp)

speciesdata$maxweight_2013<-speciesdata$maxweight_intercept+speciesdata$maxweight_slope*(speciesdata$temp_2013-speciesdata$mean_temp)
speciesdata$maxweight_RCP45<-speciesdata$maxweight_intercept+speciesdata$maxweight_slope*(speciesdata$temp_rcp45-speciesdata$mean_temp)
speciesdata$maxweight_RCP85<-speciesdata$maxweight_intercept+speciesdata$maxweight_slope*(speciesdata$temp_rcp85-speciesdata$mean_temp)

speciesdata$matureweight_2013<-speciesdata$matureweight_intercept+speciesdata$matureweight_slope*(speciesdata$temp_2013-speciesdata$mean_temp)
speciesdata$maxtureweight_RCP45<-speciesdata$matureweight_intercept+speciesdata$matureweight_slope*(speciesdata$temp_rcp45-speciesdata$mean_temp)
speciesdata$maxtureweight_RCP85<-speciesdata$matureweight_intercept+speciesdata$matureweight_slope*(speciesdata$temp_rcp85-speciesdata$mean_temp)

write_csv(speciesdata,"speciesdata_SDM.csv")
