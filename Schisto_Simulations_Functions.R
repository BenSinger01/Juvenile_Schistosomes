source("KK_sensitivity.R")

create_population <- function(lambda_all,dispersion_param, mu_param,N=500) {
  
  # Create columns
  Number <- seq(1:N)
  Age <- 18 # [] Consider distribution later 
  Juvenile_worm1 <- rep(0,N)
  Juvenile_worm2 <- rep(0,N)
  Juvenile_worm3 <- rep(0,N)
  Juvenile_worm4 <- rep(0,N)
  Juvenile_worm5 <- rep(0,N)
  Juvenile_worm6 <- rep(0,N)
  Juvenile_worm7 <- rep(0,N)
  Juvenile_worm8 <- rep(0,N)
  Adult_worm <- rep(0,N)
  Infection_binary <- rep(0,N)
  # Test_1dx2s <- rep(0,N)
  Test_3dx2s <- rep(0,N)
  # Test_1dx1s <- rep(0,N)
  EPG <- rep(0,N)
  # Cure_binary <- rep(0,N)
  # EPG_reduction <- rep(0,N)
  
  Person_risk_beta <- rlnorm(n=N,sdlog=dispersion_param,meanlog=mu_param)

  # Create dataframe
  Population <- data.frame(Number, Person_risk_beta, Age, Juvenile_worm1, Juvenile_worm2, Juvenile_worm3, Juvenile_worm4, Juvenile_worm5, Juvenile_worm6, Juvenile_worm7, Juvenile_worm8, Adult_worm, Infection_binary, Test_3dx2s, EPG)
  
  Population
}


schisto_drug_sim <- function(df, lambda_all, worm_death,wormpair_EPG_conversion,no_cap=FALSE,infect_dispersion=NULL,nbinom_trunc=NULL,juvenile_duration=6) { 
  N <- dim(df)[1]

  for (i in 1:N) {
    
    # Step 1 (Worm aging step)
    worm_stage <- df[paste("Juvenile_worm",as.character(juvenile_duration),sep='')]
    Adult_worm_aging <- worm_stage[[i,1]]
    df$Juvenile_worm8[i]<- df$Juvenile_worm7[i]
    df$Juvenile_worm7[i]<- df$Juvenile_worm6[i]
    df$Juvenile_worm6[i]<- df$Juvenile_worm5[i]
    df$Juvenile_worm5[i]<- df$Juvenile_worm4[i]
    df$Juvenile_worm4[i]<- df$Juvenile_worm3[i]
    df$Juvenile_worm3[i]<- df$Juvenile_worm2[i]
    df$Juvenile_worm2[i]<- df$Juvenile_worm1[i]
    df$Juvenile_worm1[i]<-0
    
    # Step 2 (Worm death step)
    Adult_worm_death <- rbinom(size=df$Adult_worm[i],n=1, prob=worm_death) # 208 (week) vs 1460 (day)
    
    # Step 3 (Infection step)
    lambda_i <- lambda_all * df$Person_risk_beta[i]
    
    if (is.null(infect_dispersion)){
      Juvenile_worm_infect <- rpois(1,lambda_i)
    } else {
      Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
      if (!(is.null(nbinom_trunc))){
        while (max(Juvenile_worm_infect)>nbinom_trunc){
          Juvenile_worm_infect<- rnbinom(1,mu=lambda_i,size=infect_dispersion)
        }
      }
    }
    
    # Update dataframe
    df$Juvenile_worm1[i] <- Juvenile_worm_infect 
    df$Adult_worm[i] <- df$Adult_worm[i] + Adult_worm_aging - Adult_worm_death

    # Step 4 (Compute EPG)
    Pairs_i <- df$Adult_worm[i]%/%2
    max_epg <- 5000
    if (no_cap){
      df$EPG[i] <- wormpair_EPG_conversion*Pairs_i
    } else {
      df$EPG[i] <- round(max_epg*(1-exp(-wormpair_EPG_conversion*Pairs_i/max_epg)))
    }

    # # Step 5 (Define infection status)
    # if (df$EPG[i]==0 & df$Infection_binary[i]==1) {
    #   df$Cure_binary[i] <- 1
    # }
    df$Infection_binary[i] <- ifelse(df$EPG[i]>0, 1, 0)

    # Step 6 (Test positivity)
    # Simulates results of two different testing regimes. Each uses two
    # slides per stool. One just uses one stool, the other uses three stools
    # sampled on separate days. Both approaches are used in the SCORE data.
    # test_1dx2s <- KK_result(df$EPG[i],2,1)
    test_3dx2s <- KK_result(df$EPG[i],2,3)
    # test_1dx1s <- KK_result(df$EPG[i],1,1)
    # df$Test_1dx2s[i] <- test_1dx2s
    df$Test_3dx2s[i] <- test_3dx2s
    # df$Test_1dx1s[i] <- test_1dx1s
  }
  df
}


burn_in_sim <- function(lambda_all,dispersion_param,mu_param,worm_death=1/208,wormpair_EPG_conversion=5,n_burn_sim=700,N=500,no_cap=FALSE,infect_dispersion=NULL,nbinom_trunc=NULL,juvenile_duration=6){
  # Create initial population
  Population_iter <- create_population(lambda_all,dispersion_param, mu_param,N=N)
  
  for (burn in 1:n_burn_sim) {
    Population_iter <- schisto_drug_sim(Population_iter, lambda_all, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration)
  }
  
  return(Population_iter)
}

strategy_sim <- function(Initial_state,lambda_all,Treatment_weeks,Adult_efficacy,Juvenile_efficacy,Time_outcome,delay_efficacy=0.5,worm_death=1/208,wormpair_EPG_conversion=5,juvenile_duration=6,N=500, coverage=1, nonadherance=0, static_force=1,binomial_treatment=FALSE,treatment_dispersion=NULL,treatment_redraw=FALSE,dynamic_egg='no_saturation',no_cap=FALSE,infect_dispersion=NULL,nbinom_trunc=NULL,followup1_week=5,followup2_week=11){
  if ((1-nonadherance)<coverage){
    print("Coverage not possible with given level of nonadherance")
  }
  # Initialize
  n_treated <- floor(nrow(Initial_state)*coverage)
  Strategy <- Initial_state
  if (binomial_treatment&!(is.null(treatment_dispersion))){
    alpha <- Adult_efficacy/treatment_dispersion
    beta <- (1-Adult_efficacy)/treatment_dispersion
    ratio_juvenile <- Juvenile_efficacy/Adult_efficacy
    if (!treatment_redraw){
      Strategy["personal_efficacy"] <- rbeta(nrow(Initial_state),alpha,beta)
    }
  }
  Strategy_plot <- matrix(0, nrow=13, ncol=(Time_outcome+1))
  Strategy_plot[1,1] <- sum(Initial_state$EPG>0)/length(Initial_state$EPG)*100 # Any
  Strategy_plot[2,1] <- sum(Initial_state$EPG>0 & Initial_state$EPG<100)/sum(Initial_state$EPG>0)*100 # Light
  Strategy_plot[3,1] <- sum(Initial_state$EPG>99 & Initial_state$EPG<400)/sum(Initial_state$EPG>0)*100 #Moderate
  Strategy_plot[4,1] <- sum(Initial_state$EPG>399)/sum(Initial_state$EPG>0)*100 # Heavy
  Strategy_plot[5,1] <- mean(Initial_state[Initial_state$EPG>0,]$EPG) # mean EPG
  Strategy_plot[6,1] <- exp(mean(log(Initial_state[Initial_state$EPG>0,]$EPG+1))) # GMI
  Strategy_plot[7,1] <- sum(Initial_state$Test_3dx2s>0)/length(Initial_state$Test_3dx2s)*100 # Any 3dx2s
  Strategy_plot[8,1] <- sum(Initial_state$Test_3dx2s>0 & Initial_state$EPG<100)/sum(Initial_state$Test_3dx2s>0)*100 # Light 3dx2s
  Strategy_plot[9,1] <- sum(Initial_state$Test_3dx2s>0 & (Initial_state$EPG>99 & Initial_state$EPG<400))/sum(Initial_state$Test_3dx2s>0)*100 #Moderate 3dx2s
  Strategy_plot[10,1] <- sum(Initial_state$Test_3dx2s>0 & Initial_state$EPG>399)/sum(Initial_state$Test_3dx2s>0)*100 # Heavy 3dx2s
  Strategy_plot[11,1] <- mean(Initial_state[Initial_state$Test_3dx2s>0,]$EPG) # mean EPG
  Strategy_plot[12,1] <- exp(mean(log(Initial_state[Initial_state$Test_3dx2s>0,]$EPG+1))) # GMI 3dx2s
  Strategy_plot[13,1] <- sum((Initial_state$Juvenile_worm1+Initial_state$Juvenile_worm2+Initial_state$Juvenile_worm3+Initial_state$Juvenile_worm4+Initial_state$Juvenile_worm5+Initial_state$Juvenile_worm6)>0)/length(Initial_state$EPG)*100 #Juvenile prevalence
  Strategy_plot[is.na(Strategy_plot[,1]),1] <- 0

  # Forecast
  for (sim_time in 1:Time_outcome) {
    if (sim_time %in% Treatment_weeks & coverage>0){
      # if the treatment time is not an exact number of years from the first
      # treatment, then it's the second round perscription for that year.
      efficacy <- Adult_efficacy
      if ((sim_time-min(Treatment_weeks))%%26!=0) {efficacy<-delay_efficacy}
      # Simulate treatment
      if ((sim_time-min(Treatment_weeks))%%26==0) {treat_rows <- sample(floor(nrow(Strategy)*(1-nonadherance)),n_treated)}
      Treated <- Strategy[treat_rows,]
      if (binomial_treatment){
        if (is.null(treatment_dispersion)){
          for (i in 1:n_treated){
            Treated$Adult_worm[i] <- rbinom(1,Treated$Adult_worm[i],(1-efficacy))
            Treated$Juvenile_worm1[i] <- rbinom(1,Treated$Juvenile_worm1[i],(1-Juvenile_efficacy))
            Treated$Juvenile_worm2[i] <- rbinom(1,Treated$Juvenile_worm2[i],(1-Juvenile_efficacy))
            Treated$Juvenile_worm3[i] <- rbinom(1,Treated$Juvenile_worm3[i],(1-Juvenile_efficacy))
            Treated$Juvenile_worm4[i] <- rbinom(1,Treated$Juvenile_worm4[i],(1-Juvenile_efficacy))
            Treated$Juvenile_worm5[i] <- rbinom(1,Treated$Juvenile_worm5[i],(1-Juvenile_efficacy))
            Treated$Juvenile_worm6[i] <- rbinom(1,Treated$Juvenile_worm6[i],(1-Juvenile_efficacy))
            if (juvenile_duration>6){
              Treated$Juvenile_worm7[i] <- rbinom(1,Treated$Juvenile_worm7[i],(1-Juvenile_efficacy))
              Treated$Juvenile_worm8[i] <- rbinom(1,Treated$Juvenile_worm8[i],(1-Juvenile_efficacy))
            }
          }
        } else {
          if (treatment_redraw){
            Strategy["personal_efficacy"] <- rbeta(nrow(Initial_state),alpha,beta)
          }
          for (i in 1:n_treated){
            pefficacy <- min(1,Strategy$personal_efficacy[i]*(efficacy/Adult_efficacy))
            personal_Juvenile_efficacy <- min(1,ratio_juvenile*pefficacy)
            Treated$Adult_worm[i] <- rbinom(1,Treated$Adult_worm[i],(1-pefficacy))
            Treated$Juvenile_worm1[i] <- rbinom(1,Treated$Juvenile_worm1[i],(1-personal_Juvenile_efficacy))
            Treated$Juvenile_worm2[i] <- rbinom(1,Treated$Juvenile_worm2[i],(1-personal_Juvenile_efficacy))
            Treated$Juvenile_worm3[i] <- rbinom(1,Treated$Juvenile_worm3[i],(1-personal_Juvenile_efficacy))
            Treated$Juvenile_worm4[i] <- rbinom(1,Treated$Juvenile_worm4[i],(1-personal_Juvenile_efficacy))
            Treated$Juvenile_worm5[i] <- rbinom(1,Treated$Juvenile_worm5[i],(1-personal_Juvenile_efficacy))
            Treated$Juvenile_worm6[i] <- rbinom(1,Treated$Juvenile_worm6[i],(1-personal_Juvenile_efficacy))
            if (juvenile_duration>6){
              Treated$Juvenile_worm7[i] <- rbinom(1,Treated$Juvenile_worm7[i],(1-personal_Juvenile_efficacy))
              Treated$Juvenile_worm8[i] <- rbinom(1,Treated$Juvenile_worm8[i],(1-personal_Juvenile_efficacy))
            }
          }
        }
      } else {
        Treated$Adult_worm <- floor(Treated$Adult_worm*(1-efficacy))
        Treated$Juvenile_worm1 <- floor(Treated$Juvenile_worm1*(1-Juvenile_efficacy))
        Treated$Juvenile_worm2 <- floor(Treated$Juvenile_worm2*(1-Juvenile_efficacy))
        Treated$Juvenile_worm3 <- floor(Treated$Juvenile_worm3*(1-Juvenile_efficacy))
        Treated$Juvenile_worm4 <- floor(Treated$Juvenile_worm4*(1-Juvenile_efficacy))
        Treated$Juvenile_worm5 <- floor(Treated$Juvenile_worm5*(1-Juvenile_efficacy))
        Treated$Juvenile_worm6 <- floor(Treated$Juvenile_worm6*(1-Juvenile_efficacy))
        if (juvenile_duration>6){
          Treated$Juvenile_worm7 <- floor(Treated$Juvenile_worm7*(1-Juvenile_efficacy))
          Treated$Juvenile_worm8 <- floor(Treated$Juvenile_worm8*(1-Juvenile_efficacy))
        }
      }
      for (r in 1:n_treated) {Strategy[treat_rows[r],] <- Treated[r,]}
      state_treat <- data.frame(Strategy)
      for (i in 1:N) { # update epg, infection binary, and test results
        # print(state_treat[i,])
        Pairs_i <- state_treat$Adult_worm[i]%/%2
        # print(Pairs_i)
        max_epg <- 5000
        if (no_cap){
          state_treat$EPG[i] <- wormpair_EPG_conversion*Pairs_i
        } else {
          state_treat$EPG[i] <- round(max_epg*(1-exp(-wormpair_EPG_conversion*Pairs_i/max_epg)))
        }
        # if (state_treat$EPG[i]==0 & state_treat$Infection_binary[i]==1) {
        #   state_treat$Cure_binary[i] <- 1
        # }
        state_treat$Infection_binary[i] <- ifelse(state_treat$EPG[i]>0, 1, 0)
        # test_1dx2s <- KK_result(state_treat$EPG[i],2,1)
        test_3dx2s <- KK_result(state_treat$EPG[i],2,3)
        # test_1dx1s <- KK_result(state_treat$EPG[i],1,1)
        # state_treat$Test_1dx2s[i] <- test_1dx2s
        state_treat$Test_3dx2s[i] <- test_3dx2s
        # state_treat$Test_1dx1s[i] <- test_1dx1s
      }
    }
    
    if (dynamic_egg=='no_saturation'){
      delay_prev <- mean(c(Strategy_plot[5,max(c(1,sim_time-4))]*Strategy_plot[1,max(c(1,sim_time-4))],Strategy_plot[5,max(c(1,sim_time-5))]*Strategy_plot[1,max(c(1,sim_time-5))],Strategy_plot[5,max(c(1,sim_time-6))]*Strategy_plot[1,max(c(1,sim_time-6))],Strategy_plot[5,max(c(1,sim_time-7))]*Strategy_plot[1,max(c(1,sim_time-7))],Strategy_plot[5,max(c(1,sim_time-8))]*Strategy_plot[1,max(c(1,sim_time-8))]))
      base_prev <- Strategy_plot[5,1]*Strategy_plot[1,1]
    } else if (dynamic_egg=='saturation'){
      base_prev <- Strategy_plot[5,1]*Strategy_plot[1,1]
      unsat_delay_prev <- mean(c(Strategy_plot[5,max(c(1,sim_time-4))]*Strategy_plot[1,max(c(1,sim_time-4))],Strategy_plot[5,max(c(1,sim_time-5))]*Strategy_plot[1,max(c(1,sim_time-5))],Strategy_plot[5,max(c(1,sim_time-6))]*Strategy_plot[1,max(c(1,sim_time-6))],Strategy_plot[5,max(c(1,sim_time-7))]*Strategy_plot[1,max(c(1,sim_time-7))],Strategy_plot[5,max(c(1,sim_time-8))]*Strategy_plot[1,max(c(1,sim_time-8))]))
      mag <- base_prev/(1-exp(-base_prev))
      delay_prev <- mag*(1 - exp(-unsat_delay_prev))
    } else {
      delay_prev <- mean(c(Strategy_plot[1,max(c(1,sim_time-4))],Strategy_plot[1,max(c(1,sim_time-5))],Strategy_plot[1,max(c(1,sim_time-6))],Strategy_plot[1,max(c(1,sim_time-7))],Strategy_plot[1,max(c(1,sim_time-8))]))
      base_prev <- Strategy_plot[1,1]
    }
    lambda_dyn <- lambda_all*(static_force + ((1-static_force) * (delay_prev/base_prev)))
    if (is.na(lambda_dyn)){
      lambda_dyn <- 0
    }
    Strategy <- schisto_drug_sim(Strategy, lambda_dyn, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration)
    if (sim_time==followup1_week){state_followup1 <- data.frame(Strategy)}
    if (sim_time==followup2_week){state_followup2 <- data.frame(Strategy)}
    
    Strategy_plot[1,sim_time+1] <- sum(Strategy$EPG>0)/length(Strategy$EPG)*100 # Any
    Strategy_plot[2,sim_time+1] <- sum(Strategy$EPG>0 & Strategy$EPG<100)/sum(Strategy$EPG>0)*100 # Light
    Strategy_plot[3,sim_time+1] <- sum(Strategy$EPG>99 & Strategy$EPG<400)/sum(Strategy$EPG>0)*100 #Moderate
    Strategy_plot[4,sim_time+1] <- sum(Strategy$EPG>399)/sum(Strategy$EPG>0)*100 # Heavy
    Strategy_plot[5,sim_time+1] <- mean(Strategy[Strategy$EPG>0,]$EPG) # mean EPG
    Strategy_plot[6,sim_time+1] <- exp(mean(log(Strategy[Strategy$EPG>0,]$EPG+1))) # GMI
    Strategy_plot[7,sim_time+1] <- sum(Strategy$Test_3dx2s>0)/length(Strategy$Test_3dx2s)*100 # Any 3dx2s
    Strategy_plot[8,sim_time+1] <- sum(Strategy$Test_3dx2s>0 & Strategy$EPG<100)/sum(Strategy$Test_3dx2s>0)*100 # Light 3dx2s
    Strategy_plot[9,sim_time+1] <- sum(Strategy$Test_3dx2s>0 & (Strategy$EPG>99 & Strategy$EPG<400))/sum(Strategy$Test_3dx2s>0)*100 #Moderate 3dx2s
    Strategy_plot[10,sim_time+1] <- sum(Strategy$Test_3dx2s>0 & Strategy$EPG>399)/sum(Strategy$Test_3dx2s>0)*100 # Heavy 3dx2s
    Strategy_plot[11,sim_time+1] <- mean(Strategy[Strategy$Test_3dx2s>0,]$EPG) # mean EPG
    Strategy_plot[12,sim_time+1] <- exp(mean(log(Strategy[Strategy$Test_3dx2s>0,]$EPG+1))) # GMI 3dx2s
    Strategy_plot[13,sim_time+1] <- sum((Strategy$Juvenile_worm1+Strategy$Juvenile_worm2+Strategy$Juvenile_worm3+Strategy$Juvenile_worm4+Strategy$Juvenile_worm5+Strategy$Juvenile_worm6)>0)/length(Strategy$EPG)*100 #Juvenile prevalence
    Strategy_plot[is.na(Strategy_plot[,sim_time+1]),sim_time+1] <- 0  
  }
  
  return(list('final_state'=Strategy,'treated'=treat_rows,'state_treat'=state_treat,'state_followup1'=state_followup1,'state_followup2'=state_followup2,'plots'=Strategy_plot))
}

outcomes_row <- function(analysis,drug, coverage, initial, followup){
  if (followup==1){
    state <- analysis$state_followup1
  } else if (followup==2){
    state <- analysis$state_followup2
  }
  n_adults <- sum(state$Adult_worm)
  n_juveniles <- sum(state$Juvenile_worm1+state$Juvenile_worm2+state$Juvenile_worm3+state$Juvenile_worm4+state$Juvenile_worm5+state$Juvenile_worm6)
  mean_epg <- mean(state[state$Infection_binary>0,]$EPG)
  mean_epg_3dx2s <- mean(state[state$Test_3dx2s>0,]$EPG)
  logepg <- log(state[state$Infection_binary>0,]$EPG)
  logepg_3dx2s <- log(state[state$Test_3dx2s>0,]$EPG)
  gmean_epg <- exp(mean(logepg[is.finite(logepg)]))
  gmean_epg_3dx2s <- exp(mean(logepg_3dx2s[is.finite(logepg_3dx2s)]))
  base_prevalence <- round(sum(initial$EPG>0)/length(initial$EPG)*100,digits = 1)
  base_prevalence_3dx2s <- round(sum(initial$Test_3dx2s>0)/length(initial$EPG)*100,digits = 1)
  prevalence <- round(sum(state$EPG>0)/length(state$EPG)*100,digits = 1)
  prevalence_3dx2s <- round(sum(state$Test_3dx2s>0)/length(state$Test_3dx2s)*100,digits = 1)
  light <- round(sum(state$EPG>0 & state$EPG<100)/sum(state$EPG>0)*100)
  light_3dx2s <- round(sum(state$Test_3dx2s & state$EPG<100)/sum(state$Test_3dx2s)*100)
  moderate <- round(sum(state$EPG>=100 & state$EPG<400)/sum(state$EPG>0)*100)
  heavy <- round(sum(state$EPG>=400)/sum(state$EPG>0)*100)
  
  initial_treated <- initial[analysis$treated,]
  is_infected <- initial_treated$Infection_binary==1
  is_light <- is_infected & initial_treated$EPG<100
  is_moderate <- is_infected & initial_treated$EPG>=100 & is_infected & initial_treated$EPG<400
  is_heavy <- is_infected & initial_treated$EPG>400
  is_infected_3dx2s <- initial_treated$Test_3dx2s==1
  is_light_3dx2s <- is_infected_3dx2s & initial_treated$EPG<100
  is_moderate_3dx2s <- is_infected_3dx2s & initial_treated$EPG>=100 & initial_treated$EPG<400
  is_heavy_3dx2s <- is_infected_3dx2s & initial_treated$EPG>400
  
  initial_cured <- initial_treated[is_infected,]
  initial_cured_light <- initial_treated[is_light,]
  initial_cured_moderate <- initial_treated[is_moderate,]
  initial_cured_heavy <- initial_treated[is_heavy,]
  initial_cured_3dx2s <- initial_treated[is_infected_3dx2s,]
  initial_cured_light_3dx2s <- initial_treated[is_light_3dx2s,]
  initial_cured_moderate_3dx2s <- initial_treated[is_moderate_3dx2s,]
  initial_cured_heavy_3dx2s <- initial_treated[is_heavy_3dx2s,]

  final_treated <- state[analysis$treated,]

  final_cured <- final_treated[is_infected,]
  final_cured_light <- final_treated[is_light,]
  final_cured_moderate <- final_treated[is_moderate,]
  final_cured_heavy <- final_treated[is_heavy,]
  final_cured_3dx2s <- final_treated[is_infected_3dx2s,]
  final_cured_light_3dx2s <- final_treated[is_light_3dx2s,]
  final_cured_moderate_3dx2s <- final_treated[is_moderate_3dx2s,]
  final_cured_heavy_3dx2s <- final_treated[is_heavy_3dx2s,]
  
  cure_rate <- (1 - (sum(final_cured$Infection_binary)/sum(initial_cured$Infection_binary)))*100
  cure_rate_light <- (1 - (sum(final_cured_light$Infection_binary)/sum(initial_cured_light$Infection_binary)))*100
  cure_rate_moderate <- (1 - (sum(final_cured_moderate$Infection_binary)/sum(initial_cured_moderate$Infection_binary)))*100
  cure_rate_heavy <- (1 - (sum(final_cured_heavy$Infection_binary)/sum(initial_cured_heavy$Infection_binary)))*100
  cure_rate_3dx2s <- (1 - (sum(final_cured_3dx2s$Test_3dx2s)/sum(initial_cured_3dx2s$Test_3dx2s)))*100
  cure_rate_light_3dx2s <- (1 - (sum(final_cured_light_3dx2s$Test_3dx2s)/sum(initial_cured_light_3dx2s$Test_3dx2s)))*100
  cure_rate_moderate_3dx2s <- (1 - (sum(final_cured_moderate_3dx2s$Test_3dx2s)/sum(initial_cured_moderate_3dx2s$Test_3dx2s)))*100
  cure_rate_heavy_3dx2s <- (1 - (sum(final_cured_heavy_3dx2s$Test_3dx2s)/sum(initial_cured_heavy_3dx2s$Test_3dx2s)))*100
  epg_red <- (initial_cured$EPG - final_cured$EPG)/initial_cured$EPG
  epg_red_3dx2s <- (initial_cured_3dx2s$EPG - final_cured_3dx2s$EPG)/initial_cured_3dx2s$EPG
  gmean_epg_treat <- exp(mean(log(final_cured$EPG)))
  gmean_epg_treat_3dx2s <- exp(mean(log(final_cured_3dx2s$EPG)))
  mean_epg_reduction <- mean(epg_red)*100
  mean_epg_reduction_3dx2s <- mean(epg_red_3dx2s)*100
  logepg_red <- log(epg_red)
  logepg_red_3dx2s <- log(epg_red_3dx2s)
  gmean_epg_reduction <- exp(mean(logepg_red[is.finite(logepg_red)]))*100
  gmean_epg_reduction_3dx2s <- exp(mean(logepg_red_3dx2s[is.finite(logepg_red_3dx2s)]))*100
  
  is_cured <- final_cured$Infection_binary==0
  is_cured_3dx2s <- final_cured$Test_3dx2s==0
  treatment_cured <- analysis$state_treat[is_cured,]
  treatment_cured_3dx2s <- analysis$state_treat[is_cured_3dx2s,]
  survived_juveniles <- treatment_cured[(treatment_cured$Juvenile_worm1 + treatment_cured$Juvenile_worm2 +treatment_cured$Juvenile_worm3 +treatment_cured$Juvenile_worm4 +treatment_cured$Juvenile_worm5 +treatment_cured$Juvenile_worm6) > 0,]
  survived_juveniles_3dx2s <- treatment_cured_3dx2s[(treatment_cured_3dx2s$Juvenile_worm1 + treatment_cured_3dx2s$Juvenile_worm2 +treatment_cured_3dx2s$Juvenile_worm3 +treatment_cured_3dx2s$Juvenile_worm4 +treatment_cured_3dx2s$Juvenile_worm5 +treatment_cured_3dx2s$Juvenile_worm6) > 0,]
  juvenile_cured <- 100*nrow(survived_juveniles)/nrow(treatment_cured)
  juvenile_cured_3dx2s <- 100*nrow(survived_juveniles_3dx2s)/nrow(treatment_cured_3dx2s)

  baseline_juveniles <- initial$Juvenile_worm1 + initial$Juvenile_worm2 + initial$Juvenile_worm3 + initial$Juvenile_worm4 + initial$Juvenile_worm5 + initial$Juvenile_worm6
  juvenile_infections <- initial[baseline_juveniles > 0,]
  juvenile_ratio_baseline <- 100*sum(juvenile_infections$Infection_binary == 1) / sum(initial$Infection_binary == 1)
  
  final_juveniles <- state$Juvenile_worm1 + state$Juvenile_worm2 + state$Juvenile_worm3 + state$Juvenile_worm4 + state$Juvenile_worm5 + state$Juvenile_worm6
  juvenile_prevalence <- 100*sum(final_juveniles>0)/length(final_juveniles)

  return(list('base_prevalence'=base_prevalence,'base_prevalence_3dx2s'=base_prevalence_3dx2s,'drug'=drug,'coverage'=coverage,'n_adults'=n_adults,'n_juveniles'=n_juveniles,'mean_epg_3dx2s'=mean_epg_3dx2s,'gmean_epg_3dx2s'=gmean_epg_3dx2s,'prevalence_3dx2s'=prevalence_3dx2s,'light_3dx2s'=light_3dx2s,'cure_rate_3dx2s'=cure_rate_3dx2s,'cure_rate_light_3dx2s'=cure_rate_light_3dx2s,'cure_rate_moderate_3dx2s'=cure_rate_moderate_3dx2s,'cure_rate_heavy_3dx2s'=cure_rate_heavy_3dx2s,'mean_epg_reduction_3dx2s'=mean_epg_reduction_3dx2s,'gmean_epg_reduction_3dx2s'=gmean_epg_reduction_3dx2s,'juvenile_cured_3dx2s'=juvenile_cured_3dx2s,'mean_epg'=mean_epg,'gmean_epg'=gmean_epg,'prevalence'=prevalence,'intensity_prev'=paste(light,moderate,heavy,sep='/'), 'cure_rate'=cure_rate,'cure_rate_light'=cure_rate_light,'cure_rate_moderate'=cure_rate_moderate,'cure_rate_heavy'=cure_rate_heavy, 'mean_epg_reduction'=mean_epg_reduction,'gmean_epg_reduction'=gmean_epg_reduction,'base_cases_juveniles'=juvenile_ratio_baseline,'juvenile_cured'=juvenile_cured,'juvenile_prevalence'=juvenile_prevalence))
}
