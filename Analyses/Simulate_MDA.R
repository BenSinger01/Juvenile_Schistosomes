## MDA simulation
## Benjamin John Singer, July 2024.
## Code to simulate MDA over a five year period.

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
t_start <- as.numeric(args[5])

# Proportion of force of infection that does not change with treatment (0 to 1)
static_force <- as.numeric(args[6])/100
# How does force of infection depend on current infections? (saturation, no_saturation, prev)
dyn_type <- args[7]
# Treatment efficacy estimate (Mean, Upper, Lower)
case <- args[8]
# Treatment strategy (Single, Double, NovelA, NovelB, NovelC, Arm1, Arm2, Arm3)
drug_strategy <- args[9]
# Proportion of the population treated (0 to 1)
coverage <- as.numeric(args[10])/100
# Proportion of the population systematically nonadherent (0 to 1-coverage)
nonadherance <- as.numeric(args[11])/100
# Scale effectiveness of treatment (0+)
eff_scale <- as.numeric(args[12])/100
# Degree of transmission reduction attributable to WASH intervention (0 to 1)
wash <- as.numeric(args[13])/100

# How many slots available to run code in parallel?
nslots <- as.numeric(args[14])
# How many simulations?
N <- as.numeric(args[15])
seeds <- 1:N
nseeds <- length(seeds)

community_size <- 500
lambda_param <- 1*(1-wash)

# period of simulation. Here, five years + time for follow-up.
T=5*52+11

# choose calibrated parameters for baseline and pzq treatment, dependent on epgpwp
if (epgpwp==4){
  preschool_ratio <- 0.61 # susceptibility ratio preschool:SAC
  adult_ratio <- 0.82 # susceptibility ratio adult:SAC
  nbinom_param <- 0.0175 # overdispersion parameter for infections
  eff_novel <- 0.999 # high-efficacy value for novel drugs
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
  # parameters for treatment
  if (case=="Mean"){ # best guess at praziquantel behaviour
    eff_base <- 0.95 # average proportion of adult schistosomes killed by treatment
    beta_param <- 0.49 # overdispersion of schistosome death in treatment
  } else if (case=="Upper"){ # higher efficacy praziquantel scenario
    eff_base <- 0.991
    beta_param <- 0.54
  } else if (case=="Lower"){ # lower efficacy praziquantel scenario
    eff_base <- 0.929
    beta_param <- 0.5
  }
} else if (epgpwp==15){
  preschool_ratio <- 0.65
  adult_ratio <- 0.85
  nbinom_param <- 1.5^(-13)
  eff_novel <- 0.999
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
  if (case=="Mean"){
    eff_base <- 0.909
    beta_param <- 0.51
  } else if (case=="Upper"){
    eff_base <- 0.976
    beta_param <- 0.37
  } else if (case=="Lower"){
    eff_base <- 0.887
    beta_param <- 0.355
  }
} else if (epgpwp==1){
  preschool_ratio <- 0.55
  adult_ratio <- 0.77
  nbinom_param <- 1.5^(-7)
  eff_novel <- 0.988
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
  if (case=="Mean"){
    eff_base <- 0.991
    beta_param <- 0.555
  } else if (case=="Upper"){
    eff_base <- 0.9965
    beta_param <- 0.495
  } else if (case=="Lower"){
    eff_base <- 0.975
    beta_param <- 0.44
  }
}

# Choose drug efficacies based on drug/strategy choice.
# Arm 1, 2, 3 are arms simulating SCORE study conditions.
if (drug_strategy=='Single'){ # treat annually with single dose of praziquantel
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_base
  eff_j <- 0
  followup <- 1
} else if (drug_strategy=='NovelA'){ # treat annually with drug effective against juveniles
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_base
  eff_j <- 1
  followup <- 1
} else if (drug_strategy=='NovelC'){ # treat annually with drug with extra efficacy against adults and juveniles
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_novel
  eff_j <- 1
  followup <- 1
} else if (drug_strategy=='NovelB'){ # treat annually with drug with extra efficacy against adults
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_novel
  eff_j <- 0
  followup <- 1
} else if (drug_strategy=='Double'){ # treat annually with spaced doses of praziquantel
  treat_weeks <- c(seq(from=1,to=T,by=52),seq(from=7,to=T,by=52))
  eff_a <- eff_base
  eff_j <- 0
  followup <- 2
} else if (drug_strategy=='Arm1'){ # treat annually with praziquantel
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_base
  eff_j <- 0
  followup <- 1
  if (setting=='Low'){
    coverage <- 0.62
  } else if (setting=='Moderate'){
    coverage <- 0.89
  }
} else if (drug_strategy=='Arm2'){ # treat first two years with praziquantel
  treat_weeks <- c(1,53)
  eff_a <- eff_base
  eff_j <- 0
  followup <- 1
  if (setting=='Low'){
    coverage <- 0.88
  } else if (setting=='Moderate'){
    coverage <- 0.95
  }
} else if (drug_strategy=='Arm3'){ # treat first and third year with praziquantel
  treat_weeks <- c(1,105)
  eff_a <- eff_base
  eff_j <- 0
  followup <- 1
  if (setting=='Moderate'){
    coverage <- 0.76
  } else if (setting=='High'){
    coverage <- 0.58
  }
}

# scale effectiveness if eff_scale != 1
if (eff_scale>=1){
  eff_a_scaled <- 1-(1-eff_a)*(1/eff_scale)
  if (eff_j>0){
    eff_j_scaled <- 1-(1-eff_j)*(1/eff_scale)
  } else {
    eff_j_scaled <- 0
  }
} else {
  eff_a_scaled <- eff_a*eff_scale
  if (eff_j>0){
    eff_j_scaled <- eff_j*eff_scale
  } else {
    eff_j_scaled <- 0
  }
}

# define parallelized version of the strategy simulation function
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

  sim <- strategy_sim(baseline_state,lambda_param,treat_weeks,eff_a_scaled,eff_j_scaled,T,epgpwp,juvenile_duration,
    nbinom_param,beta_param,community_size,coverage,nonadherance,static_force,dyn_type,seasonality,t_start)
  return(sim$plots)
}

# run simulations in parallel
analysis <- mclapply(seeds,parallel_plots,mc.cores=nslots)

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

# Write to csv
dir.create("Data/MDA_sims/", showWarnings = FALSE)
if (seasonality==0){
  if ((epgpwp==4)&(juvenile_duration==6)){
    dir.create(paste("Data/MDA_sims/",setting,"/",sep=''), showWarnings = FALSE)
    filename <- paste("Data/MDA_sims/",setting,"/",drug_strategy,"_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),"_",as.character(eff_scale*100),"scale_",as.character(wash*100),"wash_",case,".csv",sep='')
  } else if ((epgpwp!=4)&(juvenile_duration==6)){
    dir.create(paste("Data/MDA_sims/",setting,"_",as.character(epgpwp),"epgpwp/",sep=''), showWarnings = FALSE)
    filename <- paste("Data/MDA_sims/",setting,"_",as.character(epgpwp),"epgpwp/",drug_strategy,"_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),"_",as.character(eff_scale*100),"scale_",as.character(wash*100),"wash_",case,".csv",sep='')
  } else if ((epgpwp==4)&(juvenile_duration!=6)){
    dir.create(paste("Data/MDA_sims/",setting,"_",as.character(juvenile_duration),"juve/",sep=''), showWarnings = FALSE)
    filename <- paste("Data/MDA_sims/",setting,"_",as.character(juvenile_duration),"juve/",drug_strategy,"_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),"_",as.character(eff_scale*100),"scale_",as.character(wash*100),"wash_",case,".csv",sep='')
  }
} else {
  dir.create(paste("Data/MDA_sims/",setting,"_",as.character(seasonality*100),"seasonal",as.character(start_t),"/",sep=''), showWarnings = FALSE)
  filename <- paste("Data/MDA_sims/",setting,"_",as.character(seasonality*100),"seasonal",as.character(start_t),"/",drug_strategy,"_",as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),"_",as.character(eff_scale*100),"scale_",as.character(wash*100),"wash_",case,".csv",sep='')
}
write.table(df, file = filename, append = FALSE, quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
