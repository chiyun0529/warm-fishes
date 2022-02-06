# warm-fishes
Investigating the effect of warming on large marine fish communities by combining species distribution modeling, life-history modeling, and size spectrum modeling (mizer).

Below is a list of files in this repository along with brief descriptions:

1. **Bayesian_SDM.R:** R code for running Bayesian SDMs. From Talluto et al. (2016) 
2. **hybrid_SDM.R:** R code for running the Bayesian SDMs and performing cross-validations
3. **interactions_2013.csv:** the interacton matrix for the 2013 mizer model
4. **interactions_rcp45.csv:** the interaction matricx for the 2100-RCP4.5 mizer model
5. **interactions_rcp85.csv:** the interaction matricx for the 2100-RCP8.5 mizer model
6. **life_history.R:** R code for performing statistical life-history analyses
7. **life_history_nefsc.csv:** population-level life history data for the statistical life-history analysis
8. **lifehis_interaction.R:** R code for building and calibrating mizer model under future climate without accounting for changes in life history and species interactions
9. **mizer.R:** R code for building mizer models for the 2013 and the two 2100 climates
10. **mizer_misc_functions.R:** miscellaneous functions written by the mizer development team to aid visual inspection of model results
11. **nefsc_species_params_2013.csv:** multispecies parameters for the 2013 mizer model
12. **nefsc_species_params_rcp45.csv:** multispecies parameters for the 2100-RCP4.5 mizer model
13. **nefsc_species_params_rcp85.csv:** multispecies parameters for the 2100-RCP8.5 mizer model
14. **perturb_test.R:** R code for evaluating community resilience
15. **resilience_bottom-up.csv:** community resilience against bottom-up perturbations
16. **resilience_top-down.csv:** community resilience against top-down perturbations
17. **sim_2013.RDS:** calibrated 2013 mizer model
18. **sim_rcp45.RDS:** calbrated 2100-RCP4.5 mizer model
19. **sim_rcp45_nointer.RDS:** calibrated 2100-RCP4.5 mizer model without accounting for changes in species interactions
20. **sim_rcp45_nolife.RDS:** calibrated 2100-RCP4.5 mizer model without accounting for changes in species life histories
21. **sim_rcp85.RDS:** calibrated 2100-RCP8.5 mizer model
22. **sim_rcp85_nointer.RDS:** calibrated 2100-RCP8.5 mizer model without accounting for changes in species interactions
23. **sim_rcp85_nolife.RDS:** calibrated 2100-RCP8.5 mizer model without accounting for changes in species life histories
24. **warming_impact.csv:** changes in abundance and total biomass for each species under the two warming scenarios
