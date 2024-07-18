## Run Snail Simulation Table
## B J Singer
## Code to run long-term MDA simulations in parallel with snails

source('Analyses/Schisto_Simulation_Functions.R')
require(parallel)

# load parameters from shell call
args <- commandArgs(trailingOnly = TRUE)
# Treatment efficacy estimate (Mean, Upper, Lower)
case <- args[1]
# Treatment strategy (Single, Double, NovelA, NovelB, NovelC, Arm1, Arm2, Arm3)
drug_strategy <- args[2]
# Proportion of the population treated
coverage <- as.numeric(args[3])/100
# Proportion of the population systematically nonadherent
nonadherance <- as.numeric(args[4])/100
# Scale effectiveness of treatment (0+)
eff_scale <- as.numeric(args[5])/100

# How many slots available to run code in parallel?
nslots <- as.numeric(args[6])
# How many simulations?
N <- as.numeric(args[7])
seeds <- 1:N
nseeds <- length(seeds)

setting <- "High"
epgpwp <- 4
juvenile_duration <- 6
community_size <- 500
lambda_param <- 1

# period of simulation. Here, five years + time for follow-up.
T=5*52+11

# calibrated setting paramters
preschool_ratio <- 0.61
adult_ratio <- 0.82
nbinom_param <- 0.0175
dispersion_param <- 2.5
mu_param <- -3.9

# snail model prameters
K_param <- 5
q_param <- 23
eta_param <- 1.55e-7
start_prev <- 0.63
start_worms <- 13

# calibrated treatment parameters
eff_novel <- 0.999
if (case=="Mean"){
  eff_base <- 0.95
  beta_param <- 0.49
} else if (case=="Upper"){
  eff_base <- 0.991
  beta_param <- 0.54
} else if (case=="Lower"){
  eff_base <- 0.929
  beta_param <- 0.5
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
  baseline_state <- readRDS(paste("Data/baselines/Snails/",as.character(seed),".RData",sep=''))
  sim <- strategy_sim(baseline_state$humans,lambda_param,treat_weeks,eff_a_scaled,eff_j_scaled,T,
    epgpwp,juvenile_duration,nbinom_param,beta_param,community_size,coverage,nonadherance,
    sn=baseline_state$snails)
  return(sim$plots)
}

# run simulations in parallel
analysis <- mclapply(seeds,parallel_plots,mc.cores=nslots)

# save to dataframe. Make sure to edit filename if performing alternate analysis
df <- data.frame(matrix(ncol=20,nrow=T*nseeds))
colnames(df) <- c('Setting','Strategy','Coverage','epgpwp','t','Prevalence True','Light True','Moderate True','High True','AMI True','GMI True','Prevalence 3dx2s','Light 3dx2s','Moderate 3dx2s','High 3dx2s','AMI 3dx2s','GMI 3dx2s','Juvenile Prevalence','Juveniles','Adults')
for (i in 1:nseeds){
  for (t in 1:T){
    df[(i-1)*T+t,1:4] <- c(setting,drug_strategy,coverage*100,epgpwp)
    df[(i-1)*T+t,5] <- t
    df[(i-1)*T+t,6:20] <- analysis[[i]][,t]
  }
}

dir.create("Data/MDA_sims/", showWarnings = FALSE)
dir.create("Data/MDA_sims/Snails/", showWarnings = FALSE)
filename=paste("Data/MDA_sims/Snails/",drug_strategy,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),"_",as.character(eff_scale*100),"scale_",case,".csv",sep='')
write.table(df, file = filename, append = FALSE, quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
