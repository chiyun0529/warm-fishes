library(tidyverse)
library(ggplot2)
library(raster)
library(rasterVis)
library(lme4)
library(RColorBrewer)

# Population-level life history data

rm(list=ls())

setwd("~/Google Drive/hybrid-SDM")

lifehis<-read_csv("life_history_nefsc.csv")

genus_all<-levels(as.factor(lifehis$genus))

genus_data<-matrix(0,ncol=3,nrow=length(genus_all))

for (i in 1:length(genus_all)){
  
  data<-subset(lifehis,genus==genus_all[i])
  genus_data[i,3]<-mean(data$temp)
  genus_data[i,1]<-genus_all[i]
  genus_data[i,2]<-nrow(data)
  
}

genus_data<-data.frame(genus_data)
colnames(genus_data)<-c("genus","num_pops","mean_temp")

lifehis$temp_anomaly<--99
lifehis$max_temp<--99
lifehis$min_temp<--99

for (i in 1:length(genus_all)){
  
  mean_temp<-as.numeric(genus_data$mean_temp[i])
  lifehis[lifehis$genus==genus_all[i],]$temp_anomaly<-lifehis[lifehis$genus==genus_all[i],]$temp-mean_temp
  lifehis[lifehis$genus==genus_all[i],]$max_temp<-max(lifehis[lifehis$genus==genus_all[i],]$temp)
  lifehis[lifehis$genus==genus_all[i],]$min_temp<-min(lifehis[lifehis$genus==genus_all[i],]$temp)
  
}

write_csv(lifehis,"life_history_nefsc.csv")

# Exclude species with less than 3 data points
lifehis<-read_csv("life_history_nefsc.csv")
lifehis_part<-subset(lifehis,subset=genus!="Peprilus"&genus!="Stenotomus"&genus!="Amblyraja"&genus!="Enchelyopus")

# Model selection indicates that glmm with only genus as a random factor is the best model for max_size
glm.maxsize_full<-lmer(max_weight~temp_anomaly+(temp_anomaly|genus),data = lifehis_part)
glm.maxsize_reduced1<-lmer(max_weight~temp_anomaly+(1|genus),data = lifehis_part)
glm.maxsize_reduced2<-lm(max_weight~temp_anomaly,data=lifehis_part)
glm.maxsize_reduced3<-lmer(max_weight~(1|genus),data=lifehis_part)

AIC(glm.maxsize_full,glm.maxsize_reduced1,glm.maxsize_reduced2,glm.maxsize_reduced3)

# Get genus-specific intercepts
coef(glm.maxsize_reduced3)


# Model selection indicates that the full model is the best model for max_weight
glm.maxweight_full<-lmer(max_weight~temp_anomaly+(temp_anomaly|genus),data = lifehis_part)
glm.maxweight_reduced1<-lmer(max_weight~temp_anomaly+(1|genus),data = lifehis_part)
glm.maxweight_reduced2<-lm(max_weight~temp_anomaly,data=lifehis_part)
glm.maxweight_reduced3<-lmer(max_weight~(1|genus),data=lifehis_part)

AIC(glm.maxweight_full,glm.maxweight_reduced1,glm.maxweight_reduced2,glm.maxweight_reduced3)

maxweight_coef<-coef(glm.maxweight_full)$genus

# Model selection indicates that the full model is the best model for weight_mature
glm.matureweight_full<-lmer(weight_maturity~temp_anomaly+(temp_anomaly|genus),data = lifehis_part)
glm.matureweight_reduced1<-lmer(weight_maturity~temp_anomaly+(1|genus),data = lifehis_part)
glm.matureweight_reduced2<-lm(weight_maturity~temp_anomaly,data=lifehis_part)
glm.matureweight_reduced3<-lmer(weight_maturity~(1|genus),data=lifehis_part)

AIC(glm.matureweight_full,glm.matureweight_reduced1,glm.matureweight_reduced2,glm.matureweight_reduced3)

matureweight_coef<-coef(glm.matureweight_full)$genus

# Model selection indicates that the full model is the best model for vonBert
glm.vonBert_full<-lmer(vonBert~temp_anomaly+(temp_anomaly|genus),data = lifehis_part)
glm.vonBert_reduced1<-lmer(vonBert~temp_anomaly+(1|genus),data = lifehis_part)
glm.vonBert_reduced2<-lm(vonBert~temp_anomaly,data = lifehis_part)
glm.vonBert_reduced3<-lmer(vonBert~(1|genus),data = lifehis_part)

AIC(glm.vonBert_full,glm.vonBert_reduced1,glm.vonBert_reduced2, glm.vonBert_reduced3)

# Get genus-specific slopes and intercepts for mature_weight

speciesdata<-read_csv("~/Desktop/hybrid-SDM/speciesdata_SDM.csv")

speciesdata$matureweight_intercept<--99
speciesdata$matureweight_slope<--99

genus<-rownames(matureweight_coef)

for (i in 1:length(genus)){
  
  speciesdata[speciesdata$genus==genus[i],]$matureweight_intercept<-matureweight_coef[rownames(matureweight_coef)==genus[i],1]
  speciesdata[speciesdata$genus==genus[i],]$matureweight_slope<-matureweight_coef[rownames(matureweight_coef)==genus[i],2]
  
}

# Get genus-specific slopes and intercepts for vonBert
vonBert_coef<-coef(glm.vonBert_full)$genus

peciesdata$vonBert_intercept<--99
speciesdata$vonBert_slope<--99

for (i in 1:length(genus)){
  
  speciesdata[speciesdata$genus==genus[i],]$vonBert_intercept<-vonBert_coef[rownames(vonBert_coef)==genus[i],1]
  speciesdata[speciesdata$genus==genus[i],]$vonBert_slope<-vonBert_coef[rownames(vonBert_coef)==genus[i],2]
  
}

speciesdata[speciesdata$vonBert_intercept==-99,]$vonBert_intercept<-NA
speciesdata[speciesdata$vonBert_slope==-99,]$vonBert_slope<-NA

write_csv(speciesdata,"~/Desktop/hybrid-SDM/speciesdata.csv")

# Get genus specific slopes and intercepts for max weight
speciesdata<-read_csv("~/Desktop/hybrid-SDM/speciesdata.csv")

genus<-rownames(maxweight_coef)

for (i in 1:length(genus)){
  
  speciesdata[speciesdata$genus==genus[i],]$maxweight_intercept<-maxweight_coef[rownames(maxweight_coef)==genus[i],1]
  speciesdata[speciesdata$genus==genus[i],]$maxweight_slope<-maxweight_coef[rownames(maxweight_coef)==genus[i],2]
  
}

speciesdata$mean_temp<--99

genus<-levels(as.factor(lifehis_part$genus))

for (i in 1:length(genus)){
  
  speciesdata2[speciesdata2$genus==genus[i],]$mean_temp<-mean(lifehis_part[lifehis_part$genus==genus[i],]$temp)
  
}

speciesdata[speciesdata$mean_temp==-99,]$mean_temp<-NA

speciesdata$temp_2013<--99

species<-levels(as.factor(nefsc_SDM$species))

temp_2013<-data.frame(matrix(-99,ncol=2,nrow=length(species)))

nefsc_SDM<-read_csv("~/Desktop/hybrid-SDM/nefsc_SDM.csv")

for (i in 1:length(species)){
  
  df<-nefsc_SDM[nefsc_SDM$species==species[i],]
  temp_2013[i,1]<-species[i]
  
  if (df$habitat[1]=="pelagic"){
  temp_2013[i,2]<-mean(df$sst_2013)
  } else if (df$habitat[1]=="demersal"){
  temp_2013[i,2]<-mean(df$sbt_2013)  
  } else {
  temp_2013[i,2]<-mean(df$sit_2013)  
  }
}

species<-levels(as.factor(speciesdata$species))

for (i in 1:length(species)){
  
  speciesdata[speciesdata$species==species[i],]$temp_2013<-temp_2013[temp_2013$X1==species[i],]$X2
  
}

speciesdata2<-speciesdata[is.na(speciesdata$vonBert_intercept)==FALSE,]

write_csv(speciesdata,"~/Desktop/hybrid-SDM/speciesdata.csv")

speciesdata$vonBert_2013<-speciesdata$vonBert_slope*(speciesdata$temp_2013-speciesdata$mean_temp)+speciesdata$vonBert_intercept

speciesdata$vonBert_RCP45<-speciesdata$vonBert_slope*2+speciesdata$vonBert_intercept

speciesdata$vonBert_RCP85<-speciesdata$vonBert_slope*4+speciesdata$vonBert_intercept

speciesdata$maxweight_2013<-speciesdata$maxweight_slope*(speciesdata$temp_2013-speciesdata$mean_temp)+speciesdata$maxweight_intercept

speciesdata$maxweight_RCP45<-speciesdata$maxweight_slope*2+speciesdata$maxweight_intercept

speciesdata$maxweight_RCP85<-speciesdata$maxweight_slope*4+speciesdata$maxweight_intercept

speciesdata$matureweight_2013<-speciesdata$matureweight_slope*(speciesdata$temp_2013-speciesdata$mean_temp)+speciesdata$matureweight_intercept

speciesdata$maxtureweight_RCP45<-speciesdata$matureweight_slope*2+speciesdata$matureweight_intercept

speciesdata$maxtureweight_RCP85<-speciesdata$matureweight_slope*4+speciesdata$matureweight_intercept

# Plotting life history traits ~ temp separately for each genus
rm(list=ls())

lifehis<-read_csv("life_history_nefsc.csv")
genus_list<-c("Alosa","Ammodytes","Anchoa","Clupea","Gadus","Glyptocephalus",
         "Limanda","Lophius","Melanogrammus","Merluccius","Myoxocephalus",
         "Paralichthys","Pseudopleuronectes","Raja","Squalus","Urophycis",
         "Zoarces")

lifehis_part<-subset(lifehis,subset=genus %in% genus_list)


# Max_weight ~ temp_anomaly with species specific regression lines
reg_all<-lm(log10(lifehis_part$max_weight)~lifehis_part$temp_anomaly)

plot(x = lifehis_part$temp_anomaly,y = log10(lifehis_part$max_weight),
     type = "n",cex.axis =2,cex.lab = 2,
     xlab = "Temperature anomaly (degree C)",
     ylab = "Log10 Maximal Size")

slope<-NULL
intercept<-NULL

for (i in 1:length(genus_list)){
  
  habitat<-lifehis_part[lifehis_part$genus==genus_list[i],]$habitat[1]
  
  reg<-lm(log10(lifehis_part[lifehis_part$genus==genus_list[i],]$max_weight)~lifehis_part[lifehis_part$genus==genus_list[i],]$temp_anomaly)
  intercept[i]<-reg$coefficients[1]
  slope[i]<-reg$coefficients[2]
  points(x = lifehis_part[lifehis_part$genus==genus_list[i],]$temp_anomaly,y = log10(lifehis_part[lifehis_part$genus==genus_list[i],]$max_weight),xlab="",ylab="",pch=16,cex.lab=1.5,cex.axis=1.5)
  
  if (habitat == "pelagic"){
  abline(reg,col = alpha(colour = "orange",alpha = 0.3),lwd=2)
    } else if (habitat == "benthopelagic"){
  abline(reg,col = alpha(colour = "purple",alpha = 0.3),lwd=2)
    } else {
  abline(reg,col = alpha(colour = "dark gray",alpha = 0.3),lwd=2)
  }
}

abline(reg_all,lty=2,lwd=2)

coeffs<-data.frame(cbind(genus_list,slope,intercept))
write_csv(coeffs,"~/Desktop/coeffs.csv")

# vonBert ~ temp_anomaly with species specific regression lines
reg_all<-lm(log10(lifehis_part$vonBert)~lifehis_part$temp_anomaly)
plot(x = lifehis_part$temp_anomaly,y = log10(lifehis_part$vonBert),
     type = "n",cex.axis =2,cex.lab = 2,
     xlab = "Temperature anomaly (degree C)",
     ylab = "Log10 von Bertalanffy coefficient")

slope<-intercept<-NULL

for (i in 1:length(genus_list)){
  
  habitat<-lifehis_part[lifehis_part$genus==genus_list[i],]$habitat[1]
  
  reg<-lm(log10(lifehis_part[lifehis_part$genus==genus_list[i],]$vonBert)~lifehis_part[lifehis_part$genus==genus_list[i],]$temp_anomaly)
  intercept[i]<-reg$coefficients[1]
  slope[i]<-reg$coefficients[2]
  points(x = lifehis_part[lifehis_part$genus==genus_list[i],]$temp_anomaly,y = log10(lifehis_part[lifehis_part$genus==genus_list[i],]$vonBert),xlab="",ylab="",pch=16,cex.lab=1.5,cex.axis=1.5)
  
  if (habitat == "pelagic"){
    abline(reg,col = alpha(colour = "orange",alpha = 0.3),lwd=2)
  } else if (habitat == "benthopelagic"){
    abline(reg,col = alpha(colour = "purple",alpha = 0.3),lwd=2)
  } else {
    abline(reg,col = alpha(colour = "dark gray",alpha = 0.3),lwd=2)
  }
  
}

abline(reg_all,lty=2,lwd=2)

coeffs<-data.frame(cbind(genus_list,slope,intercept))
write_csv(coeffs,"~/Desktop/coeffs.csv")


# weight_maturity ~ temp_anomaly with species specific regression lines
reg_all<-lm(log10(lifehis_part$weight_maturity)~lifehis_part$temp_anomaly)
plot(x = lifehis_part$temp_anomaly,y = log10(lifehis_part$weight_maturity),
     type = "n",cex.axis =2,cex.lab = 2,
     xlab = "Temperature anomaly (degree C)",
     ylab = "Log10 Weight at maturity")

slope<-intercept<-NULL

for (i in 1:length(genus_list)){
  
  habitat<-lifehis_part[lifehis_part$genus==genus_list[i],]$habitat[1]
  
  reg<-lm(log10(lifehis_part[lifehis_part$genus==genus_list[i],]$weight_maturity)~lifehis_part[lifehis_part$genus==genus_list[i],]$temp_anomaly)
  intercept[i]<-reg$coefficients[1]
  slope[i]<-reg$coefficients[2]
  points(x = lifehis_part[lifehis_part$genus==genus_list[i],]$temp_anomaly,y = log10(lifehis_part[lifehis_part$genus==genus_list[i],]$weight_maturity),xlab="",ylab="",pch=16,cex.lab=1.5,cex.axis=1.5)
  
  if (habitat == "pelagic"){
    abline(reg,col = alpha(colour = "orange",alpha = 0.3),lwd=2)
  } else if (habitat == "benthopelagic"){
    abline(reg,col = alpha(colour = "purple",alpha = 0.3),lwd=2)
  } else {
    abline(reg,col = alpha(colour = "dark gray",alpha = 0.3),lwd=2)
  }

}

abline(reg_all,lty=2,lwd=2)

coeffs<-data.frame(cbind(genus_list,slope,intercept))
write_csv(coeffs,"~/Desktop/coeffs.csv")




