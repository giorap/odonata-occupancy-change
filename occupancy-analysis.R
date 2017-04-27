####################################################################################################
################################# odonata-occupancy-change ######################################### 
####################################################################################################
################################### OCCUPANCY ANALYSIS #############################################
#### R code originally developed by Daniel Turek
### Load packages
require("nimble")
require("coda")
# require("methods") # needed for running on SB cluster
# require("foreach") # needed for running for loop in parallel
### Source nimble definitions
source("src/nimble-functions.R")
### Load data set for model fitting. Load dataset_full for Full dataset and dataset_opportunistic for Opportunistic dataset
inputData <- readRDS("data/dataset_opportunistic.rds")

# #### Increasing memory limit for R
# memory.limit()
# memory.limit(size = 7000) # Size in Mb
# memory.limit()

### Define model in BUGS/NIMBLE language
code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:7) {
        beta[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[1]*year[i] + beta[2]*min_temp[i] + beta[3]*total_precip[i]
        logit(p_obs[i]) <- beta[4] + beta[5]*list_length[i] + beta[6]*jdn[i] + beta[7]*jdn2[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

### Initiate for loop to run through all odonata species
## Generate character vector of species names
species_names <- names(inputData)[-c(1:9)]
## Create output directory
dir.create("occupancy_models", showWarnings = FALSE)
## Run occupancy model across all species
for (species_name in species_names) {
        source("run-occupancy-model.R")
        saveRDS(model_output, paste("occupancy_models/", species_name, "_model.rds", sep = ""))
}
