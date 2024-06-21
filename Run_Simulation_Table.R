## Run Simulation Table
## B J Singer
## Code to run long-term MDA simulations in parallel

# choose calibrated parameters for baseline and pzq treatment, dependent on epgpwp
if (epgpwp==4){
  preschool_ratio <- 0.61
  adult_ratio <- 0.82
  nbinom_param <- 0.0175
  eff_novel <- 0.999
  if (setting=="Elimination"){
    dispersion_param <- 1.7
    mu_param <- -7.2
  } else if (setting=="Low"){
    dispersion_param <- 2.5
    mu_param <- -7.0
  } else if (setting=="Moderate"){
    dispersion_param <- 2.4
    mu_param <- -5.4
  } else if (setting=="High"){
    dispersion_param <- 2.5
    mu_param <- -3.9
  }
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
} else if (epgpwp==15){
  preschool_ratio <- 0.65
  adult_ratio <- 0.85
  nbinom_param <- 1.5^(-13)
  eff_novel <- 0.999
  if (setting=="Elimination"){
    dispersion_param <- 0.3
    mu_param <- -6.3
  } else if (setting=="Low"){
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
  if (setting=="Elimination"){
    dispersion_param <- 2.4
    mu_param <- -7.5
  } else if (setting=="Low"){
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
# Labelling for novel drugs changed during write-up:
#   novel1 here corresponds to Novel Drug A in manuscript
#   novelC corresponds to Novel Drug B
#   novelB correspons to Novel Drug C
# Arm 1, 2, 3 are arms simulating SCORE study conditions.
if (drug_strategy=='single'){
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_base
  eff_j <- 0
  followup <- 1
} else if (drug_strategy=='novel1'){
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_base
  eff_j <- 1
  followup <- 1
} else if (drug_strategy=='novel2'){
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- 1-(1-eff_base)/2
  eff_j <- 1
  followup <- 1
} else if (drug_strategy=='novelB'){
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_novel
  eff_j <- 1
  followup <- 1
} else if (drug_strategy=='novelC'){
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_novel
  eff_j <- 0
  followup <- 1
} else if (drug_strategy=='double'){
  treat_weeks <- c(seq(from=1,to=T,by=52),seq(from=7,to=T,by=52))
  eff_a <- eff_base
  eff_j <- 0
  followup <- 2
} else if (drug_strategy=='arm1'){
  treat_weeks <- seq(from=1,to=T,by=52)
  eff_a <- eff_base
  eff_j <- 0
  followup <- 1
  if (setting=='Low'){
    coverage <- 0.62
  } else if (setting=='Moderate'){
    coverage <- 0.89
  }
} else if (drug_strategy=='arm2'){
  treat_weeks <- c(1,53)
  eff_a <- eff_base
  eff_j <- 0
  followup <- 1
  if (setting=='Low'){
    coverage <- 0.88
  } else if (setting=='Moderate'){
    coverage <- 0.95
  }
} else if (drug_strategy=='arm3'){
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
  if (epgpwp==4){
    baseline_state <- readRDS(paste("baselines/",setting,"/",as.character(seed),".RData",sep=''))
  } else {
    baseline_state <- readRDS(paste("baselines/",setting,"_",as.character(epgpwp),"epgpwp/",as.character(seed),".RData",sep=''))
  }
  sim <- strategy_sim(baseline_state,lambda_param
    ,treat_weeks,eff_a_scaled,eff_j_scaled,T,delay_efficacy=eff_a_scaled
    ,N=community_size,coverage=coverage,nonadherance=nonadherance
    ,static_force=static_force,dynamic_egg=dyn_type
    ,binomial_treatment=TRUE,treatment_dispersion=beta_param,treatment_redraw=redraw
    ,no_cap=nocap,wormpair_EPG_conversion=epgpwp
    ,infect_dispersion=nbinom_param,nbinom_trunc=100)
  return(sim$plots)
}

# run simulations in parallel
analysis <- mclapply(seeds,parallel_plots,mc.cores=nslots)

# save to dataframe. Make sure to edit filename if performing alternate analysis
df <- data.frame(matrix(ncol=20,nrow=T*nseeds))
colnames(df) <- c('Setting','Strategy','Coverage','Static Force','Dynamic Type','epgpwp','t','Prevalence True','Light True','Moderate True','High True','AMI True','GMI True','Prevalence 3dx2s','Light 3dx2s','Moderate 3dx2s','High 3dx2s','AMI 3dx2s','GMI 3dx2s','Juvenile Prevalence')
for (i in 1:nseeds){
  for (t in 1:T){
    df[(i-1)*T+t,1:6] <- c(setting,drug_strategy,coverage*100,static_force*100,dyn_type,epgpwp)
    df[(i-1)*T+t,7] <- t
    df[(i-1)*T+t,8:20] <- analysis[[i]][,t]
  }
}

# output to folder Tables_smallsample if N<500
if (N<500){
  tag = "_smallsample"
} else {
  tag = ""
}

# add epgpwp to end of folder name if not 4
if (epgpwp==4){
  tag2 = ""
} else {
  tag2 = as.character(epgpwp)
}

# output simulation as csv
filename=paste('Tables',tag,tag2,'/',drug_strategy,'_',setting,'_',as.character(static_force*100),dyn_type,"_",as.character(coverage*100),"cov",as.character(nonadherance*100),"_",epgpwp,"epgpwp_",case,".csv",sep='')
write.table(df, file = filename, append = TRUE, quote = FALSE, sep = ',', row.names = FALSE, col.names = !file.exists(filename))
