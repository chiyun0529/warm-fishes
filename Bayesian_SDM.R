library(readr)
library(stringr)
library(tidyverse)
library(gmailr)
library(ncdf4)
library(sp)
library(rjags)
library(coda)
library(dplyr)
library(ecospat)
library(parallel)
library(runjags)
library(ggplot2)              
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)


SDM<-function(species,presence_ratio,ctmax){
  
  setwd(paste0("C:/users/USER/Desktop/hybrid-SDM/",species))
  
  nefsc$presence<-ifelse(nefsc$scientificName == species,1,0)
  
  presence<-subset(nefsc,subset = presence ==1)
  absence<-subset(nefsc,subset = presence ==0)
  
  # Create "true" absences 
  true_absence<-anti_join(x = absence,y = presence,by = c("grid_lon","grid_lat","depth0"))
  
  # Full dataset with presences and true absences
  occurrence<-rbind(presence,true_absence)
  
  
  ################################
  # Create training and testing datasets for cross-validation
  for (i in 1:20) {
    
    # Randomly select 70% of presence as training data
    training_rows <- sample(x = c(1:nrow(presence)), size = round(0.7 * nrow(presence)), replace = FALSE)
   
    training <- presence[training_rows, ]
    testing <- presence[-training_rows, ]
    
    # If presence:absence >=1:5, randomly select 70% of absence as training data
    if (presence_ratio>=0.2){
    training_rows_others <- sample(x = c(1:nrow(true_absence)), size = round(0.7 * nrow(true_absence)), replace = FALSE)
    others_training <- true_absence[training_rows_others, ]
    others_testing <- true_absence[-training_rows_others, ]
    
    } else {
      
    # If presence:absence < 1:5, randomly select 5 x nrow(training) from absence as training data
    training_rows_others <- sample(x = c(1:nrow(true_absence)), size = 5 * nrow(training), replace = FALSE)
    others_training <- true_absence[training_rows_others, ]
    
    others_leftover<-true_absence[-training_rows_others,]
    
    # Randomly select 5 x nrow(testing) from the leftover absence as testing data
    others_testing_rows<-sample(x = c(1:nrow(others_leftover)),size = 5*nrow(testing),replace = FALSE)
    others_testing <- others_leftover[others_testing_rows,]
    }
    
    # Make training and testing data
    training_data <- rbind(training, others_training)
    testing_data <- rbind(testing, others_testing)
    
    write_csv(training_data, paste0("C:/users/USER/Desktop/cross-validation/training_data_", 
                                    i, ".csv"))
    write_csv(testing_data, paste0("C:/users/USER/Desktop/cross-validation/testing_data_", 
                                   i, ".csv"))
    }
  
  ################################
  # If ctmax is available, create a physiology model
  if (is.na(ctmax)==FALSE){
    # The physiological model
    T <- seq(from = 0, to = 40, length.out = 100)

    Tmin = min(presence$iap_t)
    Tmax = max(presence$iap_t)

    CTmax<-ctmax

    b <- (Tmin + Tmax)/2

    c_ctmax <- sqrt(((CTmax - b)^2)/(2 * abs(log(0.05))))
    
    gaussian <- function(x, c) {
      exp(-(x - b)^2/(2 * c^2))
    }

    psi <- gaussian(x = T, c = c_ctmax)

    m2Data <- list(psi = psi, temp = T)
    m2Pars <- c("a0", "a1", "a2", "phi")

    model2 <- jags.model(paste0(species,"_model2.jags"), data = m2Data,
                         n.chains = settings$chains, n.adapt = settings$tuning)

    update(model2, settings$burnin)

    m2Results <- coda.samples(model2, m2Pars, settings$samples, thin = settings$thin)

    if (settings$diagnostics) {
      print(gelman.diag(m2Results))
    }

    save(m2Results, file = paste0("C:/users/USER/Desktop/hybrid-SDM/",species,"/",species,"_m2.rdata"))
  }

  
  ################################
  boyce.results<-data.frame(matrix(0,nrow=20,ncol=4))
  colnames(boyce.results)<-c("Boyce_index_SDM","Kmax_SDM","Boyce_index_hybrid-SDM","Kmax_hybrid-SDM")
  
  for (i in 1:20){
    
    training_data<-read_csv(paste0("C:/users/USER/Desktop/cross-validation/training_data_",i,".csv"))  
    testing_data<-read_csv(paste0("C:/users/USER/Desktop/cross-validation/testing_data_",i,".csv"))
    
    # Normal SDM
    bPrior <- matrix(c(
      0.0, 1.0E-4,    # mean and tau for b0 (intercept)
      0.0, 1.0E-4,    # b1 (first temp parameter)
      0.0, 1.0E-4),    # b2 (temp^2 parameter)
      byrow=TRUE, ncol=2)
    
    m1Data_tr <- list(temp = training_data$iap_t, presence = training_data$presence, bPrior=bPrior)
    m1Pars <- c('b0', 'b1', 'b2')
    
    model1<-jags.model(paste0(species,"_metamodel.jags"), data = m1Data_tr, n.chains = settings$chains, n.adapt = settings$tuning)
    
    update(model1, settings$burnin)
    
    m1Results_tr <- coda.samples(model1, m1Pars, settings$samples, thin=settings$thin)
    
    save(m1Results_tr,file = paste0("C:/users/USER/Desktop/hybrid-SDM/",species,"/",species,"_m1.rdata"))
    
    # Use temperature from the testing dataset for making predictions
    tempPrediction <- unique(testing_data$iap_t)
    
    m1Predictions_tr <- list(
      tempPredict = predict_psi(m1Results_tr, tempPrediction, quantiles=settings$quantiles),
      tempDomain = data.frame(temp = tempPrediction))
    
    testing_data$suitability_m1<-0
    
    for (j in 1:length(tempPrediction)){
    
    testing_data[testing_data$iap_t==tempPrediction[j],]$suitability_m1<-m1Predictions_tr$tempPredict$mean[j]
    }
    
    boyce_obs<-subset(testing_data,subset = presence == 1)
    
    boyce.results[i,1]<-ecospat.boyce(fit = testing_data$suitability_m1,obs = boyce_obs$suitability_m1)$Spearman.cor
    boyce.results[i,2]<-ecospat.max.kappa(Pred = testing_data$suitability_m1,Sp.occ = testing_data$presence)$max.Kappa
  
    ################################
    # Build the metamodel if ctmax is available
    if (is.na(ctmax)==FALSE){
    load(paste0(species,"_m2.rdata"))
    
    m2Stats <- as.data.frame(summary(m2Results)$statistics)
    m2Stats$tau <- 1/(m2Stats$SD^2)
    
    precisionReduction <- 1
    
    bPrior <- matrix(c(
      m2Stats$Mean[1],m2Stats$tau[1]/precisionReduction, # mean and tau for b0 (intercept)
      m2Stats$Mean[2],m2Stats$tau[2]/precisionReduction,  # b1 (first temp parameter)
      m2Stats$Mean[3],m2Stats$tau[3]/precisionReduction), # b2 (temp^2 parameter) - this and the one above are uninformative because we only have info on precipitation
      byrow=TRUE, ncol=2)
    
    # Build the metamodel
    mmData_tr <- list(
      presence = training_data$presence,
      temp = training_data$iap_t,
      bPrior = bPrior)
    mmPars <- c('b0', 'b1', 'b2')
    
    metamodel <- jags.model(paste0(species,"_metamodel.jags"), data = mmData_tr,  n.chains=settings$chains, n.adapt=settings$tuning)
    update(metamodel, settings$burnin)
    
    mmResults_tr <- coda.samples(metamodel, mmPars, settings$samples, thin=settings$thin)
    
    save(mmResults_tr,file = paste0("C:/users/USER/Desktop/hybrid-SDM/",species,"/",species,"_mm.rdata"))
    
    tempPrediction <- unique(testing_data$iap_t)
    
    mmPredictions_tr <- list(
      # This "tempPredict" predicts the probabilities of occurrence for every known herring location
      tempPredict = predict_psi(mmResults_tr, tempPrediction, quantiles=settings$quantiles),
      tempDomain = data.frame(temp = tempPrediction))
    
    testing_data$suitability_mm<-0
    
    for (j in 1:length(tempPrediction)){
      
      testing_data[testing_data$iap_t==tempPrediction[j],]$suitability_mm<-mmPredictions_tr$tempPredict$mean[j]
    }
    
    
    boyce_obs<-subset(testing_data,subset = presence ==1)
    
    boyce.results[i,3]<-ecospat.boyce(fit = testing_data$suitability_mm,obs = boyce_obs$suitability_mm)$Spearman.cor
    boyce.results[i,4]<-ecospat.max.kappa(Pred = testing_data$suitability_mm,Sp.occ = testing_data$presence)$max.Kappa
    }
  }
  write_csv(boyce.results,"boyce_results.csv")
}
