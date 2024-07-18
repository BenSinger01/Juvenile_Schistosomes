## MDA simulation for grid figures
## Benjamin John Singer, July 2024.
## Code to simulate MDA over a five year period, with variable efficacy against adults and juveniles.

source('Analyses/Schisto_Simulation_Functions.R')
require(parallel)

# Load parameters from shell call
args <- commandArgs(trailingOnly = TRUE)
# Average baseline fecundity of worm pairs (1, 4, 15)
epgpwp <- as.numeric(args[1]) # epgpwp = eggs per gram per worm pair
# Length of juvenile life stage (4 to 10)
juvenile_duration <- as.numeric(args[2])
# Endemic setting (Low, Moderate, High)
setting <- args[3]
# How much does infection risk change seasonally? (0 to 1)
seasonality <- as.numeric(args[4])/100
# At what week does the simulation start, i.e. what point on the cosine curve
start_t <- as.numeric(args[5])

# Proportion of force of infection that does not change with treatment (0 to 1)
static_force <- as.numeric(args[6])/100
# How does force of infection depend on current infections? (saturation, no_saturation, prev)
dyn_type <- args[7]
# Proportion of the population treated (0 to 1)
coverage <- as.numeric(args[8])/100
# Proportion of the population systematically nonadherent (0 to 1-coverage)
nonadherance <- as.numeric(args[9])/100

# Average proportion of adults killed by treatment
eff_a <- as.numeric(args[10])/100
# Average proportion of juveniles killed by treatment
eff_j <- as.numeric(args[11])/100

# How many slots available to run code in parallel?
nslots <- as.numeric(args[12])
# How many simulations?
N <- as.numeric(args[13])
seeds <- 1:N
nseeds <- length(seeds)

community_size <- 500
lambda_param <- 1

# period of simulation. Here, five years + time for follow-up.
T=5*52+11

# choose calibrated parameters for baseline, dependent on epgpwp
if (epgpwp==4){
  preschool_ratio <- 0.61 # susceptibility ratio preschool:SAC
  adult_ratio <- 0.82 # susceptibility ratio adult:SAC
  nbinom_param <- 0.0175 # overdispersion parameter for infections
  beta_param <- 0.49 # overdispersion parameter for treatment
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
  beta_param <- 0.51
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
  beta_param <- 0.555
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

treat_weeks <- seq(from=1,to=T,by=52)
followup <- 1

# define parallel function
parallel_plots <- function(seed){
  set.seed(seed)
  # load pre-run baseline states
  if (seasonality==0){
    if ((epgpwp==4)&(juvenile_duration==6)){
      baseline_state <- readRDS(paste("Data/baselines/",setting,"/",as.character(seed),".RData",sep=''))
    } else if ((epgpwp!=4)&(juvenile_duration==6)){
      baseline_state <- readRDS(paste("Data/baselines/",setting,"_",as.character(epgpwp),"epgpwp/",as.character(seed),".RData",sep=''))
    } else if ((epgpwp==4)&(juvenile_duration!=6)){
      baseline_state <- readRDS(paste("Data/baselines/",setting,"_",as.character(juvenile_duration),"juve/",as.character(seed),".RData",sep=''))
    }
  } else {
    baseline_state <- readRDS(paste("Data/baselines/",setting,"_",as.character(seasonality*100),"seasonal",as.character(t_start),"/",as.character(seed),".RData",sep=''))
  }
  sim <- strategy_sim(baseline_state,lambda_param,treat_weeks,eff_a,eff_j,T,epgpwp,juvenile_duration,
    nbinom_param,beta_param,community_size,coverage,nonadherance,static_force,dyn_type,seasonality,start_t)
  return(sim$plots)
}

# run simulations in parallel
analysis <- mclapply(seeds,parallel_plots,mc.cores=nslots)

drug_strategy <- 'single'
# save to dataframe. Make sure to edit filename if performing alternate analysis
df <- data.frame(matrix(ncol=22,nrow=T*nseeds))
colnames(df) <- c('Setting','Strategy','Coverage','Static Force','Dynamic Type','epgpwp','t','Prevalence True','Light True','Moderate True','High True','AMI True','GMI True','Prevalence 3dx2s','Light 3dx2s','Moderate 3dx2s','High 3dx2s','AMI 3dx2s','GMI 3dx2s','Juvenile Prevalence','Juveniles','Adults')
for (i in 1:nseeds){
  for (t in 1:T){
    df[(i-1)*T+t,1:6] <- c(setting,drug_strategy,coverage*100,static_force*100,dyn_type,epgpwp)
    df[(i-1)*T+t,7] <- t
    df[(i-1)*T+t,8:22] <- analysis[[i]][,t]
  }
}

dir.create("Data/MDA_sims/", showWarnings = FALSE)
dir.create("Data/MDA_sims/Grid_sims", showWarnings = FALSE)
if (seasonality==0){
  if ((epgpwp==4)&(juvenile_duration==6)){
    dir.create(paste("Data/MDA_sims/Grid_sims/",setting,"/",sep=''), showWarnings = FALSE)
    filename <- paste("Data/MDA_sims/Grid_sims/",setting,"/",as.character(eff_a*100),"a",as.character(eff_j*100),"j","_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),".csv",sep='')
  } else if ((epgpwp!=4)&(juvenile_duration==6)){
    dir.create(paste("Data/MDA_sims/Grid_sims/",setting,"_",as.character(epgpwp),"epgpwp/",sep=''), showWarnings = FALSE)
    filename <- paste("Data/MDA_sims/Grid_sims/",setting,"_",as.character(epgpwp),"epgpwp/",as.character(eff_a*100),"a",as.character(eff_j*100),"j","_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),".csv",sep='')
  } else if ((epgpwp==4)&(juvenile_duration!=6)){
    dir.create(paste("Data/MDA_sims/Grid_sims/",setting,"_",as.character(juvenile_duration),"juve/",sep=''), showWarnings = FALSE)
    filename <- paste("Data/MDA_sims/Grid_sims/",setting,"_",as.character(juvenile_duration),"juve/",as.character(eff_a*100),"a",as.character(eff_j*100),"j","_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),".csv",sep='')
  }
} else {
  dir.create(paste("Data/MDA_sims/Grid_sims/",setting,"_",as.character(seasonality*100),"seasonal",as.character(start_t),"/",sep=''), showWarnings = FALSE)
  filename <- paste("Data/MDA_sims/Grid_sims/",setting,"_",as.character(seasonality*100),"seasonal",as.character(start_t),"/",as.character(eff_a*100),"a",as.character(eff_j*100),"j","_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),".csv",sep='')
}
write.table(df, file = filename, append = FALSE, quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
