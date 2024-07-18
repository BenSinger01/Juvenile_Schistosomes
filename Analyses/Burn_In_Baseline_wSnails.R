## Baseline States
## Benjamin John Singer, July 2024.
## Script to burn in baseline endemic states for MDA simulations.

source('Analyses/Schisto_Simulation_Functions.R')

# load parameters from shell call
args <- commandArgs(trailingOnly = TRUE)
# random seed
seed <- args[1]
set.seed(as.numeric(seed))

setting <- "High"
epgpwp <- 4
juvenile_duration <- 6
community_size <- 500
lambda_param <- 1

# choose calibrated parameters for baseline state
preschool_ratio <- 0.61 # susceptibility ratio preschool:SAC
adult_ratio <- 0.82 # suceptibility ration adult:SAC
nbinom_param <- 0.0175 # overdispersion parameter for infections
# parameters for shape of lognormal susceptibility distribution:
dispersion_param <- 2.5
mu_param <- -3.9

# snail model prameters
K_param <- 5
q_param <- 23
eta_param <- 1.55e-7
start_prev <- 0.63
start_worms <- 13

# generate baseline state
baseline_state <- burn_in_sim(lambda_param,dispersion_param,mu_param,epgpwp,2000,
    N=community_size,infect_dispersion=nbinom_param,juvenile_duration=juvenile_duration,
    snails=TRUE,K=K_param,eta=eta_param,q=q_param,start_prev=start_prev,start_worms=start_worms)

# save in baselines folder and appropriate subfolder depending on parameters
dir.create("Data/baselines/", showWarnings = FALSE)
dir.create("Data/baselines/Snails/", showWarnings = FALSE)
saveRDS(baseline_state,file=paste("Data/baselines/Snails/",as.character(seed),".RData",sep=''))