
# choose calibrated parameters for baseline state, dependent on epgpwp
if (epgpwp==4){
  preschool_ratio <- 0.61
  adult_ratio <- 0.82
  nbinom_param <- 0.0175
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
} else if (epgpwp==15){
  preschool_ratio <- 0.65
  adult_ratio <- 0.85
  nbinom_param <- 1.5^(-13)
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
} else if (epgpwp==1){
  preschool_ratio <- 0.55
  adult_ratio <- 0.77
  nbinom_param <- 1.5^(-7)
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
}

# run and save basline state
baseline_state <- burn_in_sim(lambda_param,dispersion_param,mu_param=mu_param,N=community_size,wormpair_EPG_conversion=epgpwp,no_cap=nocap,infect_dispersion=nbinom_param,nbinom_trunc=100,juvenile_duration=juvenile_duration,
  preschool_ratio=preschool_ratio,adult_ratio=adult_ratio)
if ((epgpwp==4)&(juvenile_duration==6)){
  saveRDS(baseline_state,file=paste("baselines/",setting,"/",as.character(seed),".RData",sep=''))
} else if ((epgpwp!=4)&(juvenile_duration==6)){
  saveRDS(baseline_state,file=paste("baselines/",setting,"_",as.character(epgpwp),"epgpwp/",as.character(seed),".RData",sep=''))
} else if ((epgpwp==4)&(juvenile_duration!=6)){
  saveRDS(baseline_state,file=paste("baselines/",setting,"_",as.character(juvenile_duration),"juve/",as.character(seed),".RData",sep=''))
}
