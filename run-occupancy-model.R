####################################################################################################
################################# odonata-occupancy-change ######################################### 
####################################################################################################
############################## RUN SINGLE OCCUPANCY MODEL ##########################################
#### R code originally developed by Daniel Turek

print(paste("Initializing ", species_name, " model", sep = ""))

#### Set number of iterations (niter) and number of iterations to discard as burn-in (burnin)
niter <- 500000
burnin <- 100000

#### Set up model constants
constants <- with(inputData,
                  list(N=N, nsite=nsite, 
                       year=year, 
                       min_temp=min_temp,
                       total_precip=total_precip,
                       list_length=list_length, 
                       jdn=jdn,
                       jdn2=jdn2, 
                       siteID=siteID)
                  )
data <- list(y = inputData[[species_name]])
inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0, 7))
modelInfo <- list(code=code, constants=constants, data=data, inits=inits, name=species_name)
        
#### Set up model and samplers
Rmodel <- nimbleModel(modelInfo$code,
                      modelInfo$constants,
                      modelInfo$data,
                      modelInfo$inits)
Cmodel <- compileNimble(Rmodel)
spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:7]')
spec$addSampler('beta[1:4]', 'RW_block') # detection sub-model sampler
spec$addSampler('beta[5:7]', 'RW_block') # occupancy sub-model sampler
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha')) # random effect sampler
spec$getSamplers() # Check samplers
# I have commented these last two lines out to reduce memory requirements for saved output. 
# However, they might have to be added back in for the final run.
# spec$addMonitors(c('p_occ')) # add a monitor to get p_occ in output
# spec$addMonitors(c('p_obs')) # add a monitor to get p_obs in output

        
#‘posterior predictive’,
#‘RW’, ‘RW block’, ‘binary’, ‘slice’, ‘crossLevel’, and ‘RW llFunction’

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
        
        
#### Run MCMC with 500,000 iterations and 100,000 burn-in
system.time(samplesList <- lapply(1:3, mcmcClusterFunction))
        
print("MODELS SUCCESSFULLY RUN!")
		
# Make mcmc.list with only covariates and mu.alpha
mcmcs <- mcmc.list(lapply(1:length(samplesList), function(x) as.mcmc(samplesList[[x]][,1:7])))

## Rhat
coda::gelman.diag(mcmcs, autoburnin = FALSE)
gelman.plot(mcmcs)
        
## Effective sample size
effectiveSize(mcmcs)
        
## Posterior Density Plots
# pdf(paste("output/trace_and_posterior_density_plots_", species_name, ".pdf", sep = ""))
plot(mcmcs[[1]], ask = FALSE)
# dev.off()
        
#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[1]], 2, mean),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.975))))
names(results) <- c("mean", "cil", "ciu")
results$params <- row.names(results)

#### Create and return model output
model_output <- list(samplesList = samplesList, mcmcs = mcmcs, results = results)