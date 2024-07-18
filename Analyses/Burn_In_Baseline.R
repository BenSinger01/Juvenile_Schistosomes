## Baseline States
## Benjamin John Singer, July 2024.
## Script to burn in baseline endemic states for MDA simulations.

source('Analyses/Schisto_Simulation_Functions.R')

# load parameters from shell call
args <- commandArgs(trailingOnly = TRUE)
# choose transmission setting (Low, Moderate, High)
setting <- args[1]
# Average baseline fecundity of worm pairs (1, 4, 15)
epgpwp <- as.numeric(args[2])
# length of juvenile life stage (4 to 10)
juvenile_duration <- as.numeric(args[3])
# degree of seasonality (0 to 1)
seasonality <- as.numeric(args[4])/100
# start time for seasonality (28 for low transmission season, 54 for high transmission season)
start_t <- as.numeric(args[5])
# random seed
seed <- args[6]
set.seed(as.numeric(seed))

community_size <- 500
lambda_param <- 1

# choose calibrated parameters for baseline state, dependent on epgpwp
if (epgpwp==4){
  preschool_ratio <- 0.61 # susceptibility ratio preschool:SAC
  adult_ratio <- 0.82 # suceptibility ration adult:SAC
  nbinom_param <- 0.0175 # overdispersion parameter for infections
  # parameters for shape of lognormal susceptibility distribution:
  if (setting=="Low"){
    dispersion_param <- 2.5
    mu_param <- -7.0
  } else if (setting=="Moderate"){
    dispersion_param <- 2.4
    mu_param <- -5.4
  } else if (setting=="High"){
    dispersion_param <- 2.5
    mu_param <- -3.9
  }
} else if (epgpwp==15){
  preschool_ratio <- 0.65
  adult_ratio <- 0.85
  nbinom_param <- 1.5^(-13)
  if (setting=="Low"){
    dispersion_param <- 1.49
    mu_param <- -6.41
  } else if (setting=="Moderate"){
    dispersion_param <- 1.41
    mu_param <- -5.51
  } else if (setting=="High"){
    dispersion_param <- 1.5
    mu_param <- -4.3
  }
} else if (epgpwp==1){
  preschool_ratio <- 0.55
  adult_ratio <- 0.77
  nbinom_param <- 1.5^(-7)
  if (setting=="Low"){
    dispersion_param <- 3.125
    mu_param <- -6.66
  } else if (setting=="Moderate"){
    dispersion_param <- 2.93
    mu_param <- -4.77
  } else if (setting=="High"){
    dispersion_param <- -2.84
    mu_param <- -3.02
  }
}

# generate baseline state
baseline_state <- burn_in_sim(lambda_param,dispersion_param,mu_param,epgpwp,
  infect_dispersion=nbinom_param,juvenile_duration=juvenile_duration,seasonality=seasonality,start_t=start_t)

# save in baselines folder and appropriate subfolder depending on parameters
dir.create(paste("Data/baselines/"), showWarnings = FALSE)
if (seasonality==0){
  if ((epgpwp==4)&(juvenile_duration==6)){
    dir.create(paste("Data/baselines/",setting,"/",sep=''), showWarnings = FALSE)
    saveRDS(baseline_state,file=paste("Data/baselines/",setting,"/",as.character(seed),".RData",sep=''))
  } else if ((epgpwp!=4)&(juvenile_duration==6)){
    dir.create(paste("Data/baselines/",setting,"_",as.character(epgpwp),"epgpwp/"), showWarnings = FALSE)
    saveRDS(baseline_state,file=paste("Data/baselines/",setting,"_",as.character(epgpwp),"epgpwp/",as.character(seed),".RData",sep=''))
  } else if ((epgpwp==4)&(juvenile_duration!=6)){
    dir.create(paste("Data/baselines/",setting,"_",as.character(juvenile_duration),"juve/"), showWarnings = FALSE)
    saveRDS(baseline_state,file=paste("Data/baselines/",setting,"_",as.character(juvenile_duration),"juve/",as.character(seed),".RData",sep=''))
  }
} else {
  dir.create(paste("Data/baselines/",setting,"_",as.character(seasonality*100),"seasonal",as.character(start_t),"/"), showWarnings = FALSE)
  saveRDS(baseline_state,file=paste("Data/baselines/",setting,"_",as.character(seasonality*100),"seasonal",as.character(start_t),"/",as.character(seed),".RData",sep=''))
}
