## Schistosomiasis Simulation Functions
## B J Singer
## Functions to simulate schistosomiasis infection, tracking juvenile schistosomes as they age

# Function to create population of uninfected individuals as base for simulations
create_population <- function(lambda_all,dispersion_param,mu_param,preschool_ratio=1,adult_ratio=1,N=500) {
  
  # Create columns
  Number <- seq(1:N)
  # Juvenile worm compartments, only first six used in main analysis
  Juvenile_worm1 <- rep(0,N)
  Juvenile_worm2 <- rep(0,N)
  Juvenile_worm3 <- rep(0,N)
  Juvenile_worm4 <- rep(0,N)
  Juvenile_worm5 <- rep(0,N)
  Juvenile_worm6 <- rep(0,N)
  Juvenile_worm7 <- rep(0,N)
  Juvenile_worm8 <- rep(0,N)
  Juvenile_worm9 <- rep(0,N)
  Juvenile_worm10 <- rep(0,N)
  # Adult worm compartment
  Adult_worm <- rep(0,N)
  # Keeping track of outcomes
  Infection_binary <- rep(0,N)
  Test_2dx2s <- rep(0,N)
  Test_3dx2s <- rep(0,N)
  # Test_1dx1s <- rep(0,N)
  EPG <- rep(0,N)
  # Cure_binary <- rep(0,N)
  # EPG_reduction <- rep(0,N)

  # Age groups: preschool-age children, school-age children, adolescents/adults
  preschool_proportion <- 0.1
  school_proportion <- 0.25
  adult_proportion <- 0.65
  Age <- sample(c('Preschool','School','Adult'),N,replace=TRUE,prob=c(preschool_proportion,school_proportion,adult_proportion))

  # Balance relative susceptability so that average susceptability is constant
  school_relative_beta <- 1/(school_proportion+preschool_proportion*preschool_ratio+adult_proportion*adult_ratio)
  preschool_relative_beta <- preschool_ratio*school_relative_beta
  adult_relative_beta <- adult_ratio*school_relative_beta
          
  # Calculate mu and dispersion so that mean of lognormal susceptiblity distribuions are proportioned according to relative risk
  mu_preschool <- log (preschool_relative_beta^2*exp(2*mu_param+dispersion_param^2)/sqrt((preschool_relative_beta^2+exp(dispersion_param^2)-1)*exp(2*mu_param+dispersion_param^2)))
  dispersion_preschool <- sqrt(2*log(sqrt((preschool_relative_beta^2+exp(dispersion_param^2)-1)*exp(2*mu_param+dispersion_param^2))/(preschool_relative_beta*exp(mu_param+.5*dispersion_param^2))))

  mu_school <- log(school_relative_beta^2*exp(2*mu_param+dispersion_param^2)/sqrt((school_relative_beta^2+exp(dispersion_param^2)-1)*exp(2*mu_param+dispersion_param^2)))
  dispersion_school <- sqrt(2*log(sqrt((school_relative_beta^2+exp(dispersion_param^2)-1)*exp(2*mu_param+dispersion_param^2))/(school_relative_beta*exp(mu_param+.5*dispersion_param^2))))

  mu_adult <- log(adult_relative_beta^2*exp(2*mu_param+dispersion_param^2)/sqrt((adult_relative_beta^2+exp(dispersion_param^2)-1)*exp(2*mu_param+dispersion_param^2)))
  dispersion_adult <- sqrt(2*log(sqrt((adult_relative_beta^2+exp(dispersion_param^2)-1)*exp(2*mu_param+dispersion_param^2))/(adult_relative_beta*exp(mu_param+.5*dispersion_param^2))))


  Person_risk_beta <- rep(0,N)
  # Assign random susceptibility from lognormal distribution to each individual, with dependence on age group
  for (i in 1:N){
    if (Age[i]=='Preschool'){
      mu <- mu_preschool
      dispersion <- dispersion_preschool
    } else if (Age[i]=="School"){
      mu <- mu_school
      dispersion <- dispersion_school
    } else if (Age[i]=="Adult"){
      mu <- mu_adult
      dispersion <- dispersion_adult
    }
    Person_risk_beta[i] <- rlnorm(n=1,sdlog=dispersion,meanlog=mu)
  }

  # Create dataframe
  Population <- data.frame(Number, Person_risk_beta, Age, Juvenile_worm1, Juvenile_worm2, Juvenile_worm3, Juvenile_worm4, Juvenile_worm5, Juvenile_worm6, Juvenile_worm7, Juvenile_worm8, Juvenile_worm9, Juvenile_worm10, Adult_worm, Infection_binary, Test_2dx2s, Test_3dx2s, EPG)
  
  Population
}

# Function to create snail compartments and keep track of snail model parameters
create_snail_population <- function(B_0=7/10,K_0=5,mu=7/100,mu_I=7/25,incubation=1/4,recovery=0,eta=1.8e-7,q=15,start=0){
  Snail_Population <- list(
    # Snail birth rate
    B_0=B_0,
    # Carrying capacity density
    K_0=K_0,
    # Snail death rate
    mu=mu,
    # Additional mortality for infected snails
    mu_I=mu_I,
    # Rate of transfer from E to I compartment
    incubation=incubation,
    # Rate of transfer from E to S compartment
    recovery=recovery,
    # Force of infection for snails for each unit egg shedding
    eta=eta,
    # Force of infection for humans for each infected snail
    q=q,
    # Time step for seasonal model
    t=start,
    # Initialize compartments to reasonable guesses at equilibrium values
    S=K_0*(1-mu/B_0)*0.9,
    E=K_0*(1-mu/B_0)*0.05,
    I=K_0*(1-mu/B_0)*0.05
    )
}

# Function to simulate on week's worth of infection and aging of schistosomes
schisto_drug_sim <- function(df, lambda_all, worm_death,wormpair_EPG_conversion,no_cap=FALSE,max_epg=5000,infect_dispersion=NULL,nbinom_trunc=NULL,juvenile_duration=6,seasonality=0,t=0,maturation_0=NULL,maturation_alpha=NULL,skip_infect=NULL,sn=NULL,sn_seasonality=FALSE) { 
  N <- dim(df)[1]

  # if a snail state has been passed, update according to snail model
  if (!is.null(sn)){
    # if carrying capacity is seasonal, update that
    if (sn_seasonality){
      K <- sn$K_0*(1+0.5*sin(2*3.14*sn$t/52))
    } else {
      K <- sn$K_0
    }

    # update snail compartments according to equations in SI
    N_sn <- sn$S + sn$E + sn$I
    egg_shedding <- sum(df$EPG)
    new_S <- sn$S + sn$B_0*(1-N_sn/K)*(sn$S+sn$E) + sn$recovery*sn$E -
     sn$eta*egg_shedding*sn$S - sn$mu*sn$S
    new_E <- sn$E + sn$eta*egg_shedding*sn$S -
     sn$incubation*sn$E - sn$mu*sn$E - sn$recovery*sn$E
    new_I <- sn$I + sn$incubation*sn$E - (sn$mu+sn$mu_I)*sn$I

    sn$S <- new_S
    sn$E <- new_E
    sn$I <- new_I
    # iterate snail time
    sn$t <- sn$t + 1

    # update FOI according to snail shedding
    lambda_all <- sn$q*sn$I*lambda_all
  }

  for (i in 1:N) {
    # Step 1 (Worm aging step)
    # From oldest to youngest compartment, age juveniles into next category
    Adult_worm_aging <- df[i,paste("Juvenile_worm",as.character(juvenile_duration),sep='')]
    for (jc in juvenile_duration:3){
      df[i,paste("Juvenile_worm",as.character(jc),sep='')] <- df[i,paste("Juvenile_worm",as.character(jc-1),sep='')]
    }
    # Check for density-dependent worm maturation (not included in final analyses)
    if (!(is.null(maturation_0))){
      # density-dependent maturation from initial compartment
      aging_worms <- rbinom(size=df$Juvenile_worm1[i],n=1,prob=maturation_0*exp(-maturation_alpha*df$Adult_worm[i]))
      df$Juvenile_worm1[i] <- df$Juvenile_worm1[i] - aging_worms
      df$Juvenile_worm2[i] <- aging_worms
    } else {
      # fixed maturation from initial compartment
      df$Juvenile_worm2[i] <- df$Juvenile_worm1[i]
      df$Juvenile_worm1[i] <- 0
    }
    
    # Step 2 (Worm death step)
    Adult_worm_death <- rbinom(size=df$Adult_worm[i],n=1, prob=worm_death) # 208 (week) vs 1460 (day)
    
    # Step 3 (Infection step)
    # define force of infection, with optional seasonality
    lambda_i <- lambda_all * df$Person_risk_beta[i] * (1 + seasonality*cos(2*pi*t/52))

    # Null or False skip_infect do same thing
    if (is.null(skip_infect)){
      if (is.null(infect_dispersion)){
        # poisson distributed infection with cercariae from the environment
        Juvenile_worm_infect <- rpois(1,lambda_i)
      } else {
        # overdispersed infection, with truncation to prevent unrealistic extreme values
        Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
        # while loop implements truncation
        if (!(is.null(nbinom_trunc))){
          while (max(Juvenile_worm_infect)>nbinom_trunc){
            Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
          }
        }
      }
    } else if (!skip_infect[i]) {
      if (is.null(infect_dispersion)){
        Juvenile_worm_infect <- rpois(1,lambda_i)
      } else {
        Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
        if (!(is.null(nbinom_trunc))){
          while (max(Juvenile_worm_infect)>nbinom_trunc){
            Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
          }
        }
      }
    } else {
      # if infection step is skipped in this time step (ie DHP-like treatment), don't infect
      Juvenile_worm_infect <- 0
    }
    
    # Update dataframe
    df$Juvenile_worm1[i] <- df$Juvenile_worm1[i] + Juvenile_worm_infect 
    df$Adult_worm[i] <- df$Adult_worm[i] + Adult_worm_aging - Adult_worm_death

    # Step 4 (Compute EPG)
    Pairs_i <- df$Adult_worm[i]%/%2
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
    # slides per stool. One uses two stools, the other uses three stools
    # sampled on separate days. The former is used for Mnkugwe et al., the
    # latter for SCORE data.
    # test_1dx2s <- KK_result(df$EPG[i],2,1)
    test_3dx2s <- KK_result(df$EPG[i],2,3)
    test_2dx2s <- KK_result(df$EPG[i],2,2)
    # test_1dx1s <- KK_result(df$EPG[i],1,1)
    # df$Test_1dx2s[i] <- test_1dx2s
    df$Test_3dx2s[i] <- test_3dx2s
    df$Test_2dx2s[i] <- test_2dx2s
    # df$Test_1dx1s[i] <- test_1dx1s
  }

  # iterate plain datafra if no snails, otherwise list
  if (is.null(sn)){
    df
  } else {
    outlist <- list("humans" = df, "snails" = sn)
    outlist
  }

}

# Function to burn in baseline state
burn_in_sim <- function(lambda_all,dispersion_param,mu_param,worm_death=1/208,wormpair_EPG_conversion=5,n_burn_sim=700,preschool_ratio=1,adult_ratio=1,N=500,no_cap=FALSE,max_epg=5000,infect_dispersion=NULL,nbinom_trunc=NULL,juvenile_duration=6,seasonality=0,start_t=28,maturation_0=NULL,maturation_alpha=NULL,snails=FALSE,K=5,eta=1.8e-7,q=15,start_prev=0.5,start_worms=20,sn_seasonality=FALSE){
  # Create initial population
  Population_iter <- create_population(lambda_all,dispersion_param,mu_param,preschool_ratio=preschool_ratio,adult_ratio=adult_ratio,N=N)
  
  # If snails are to be included, create a snail population
  if (snails){
    sn <- create_snail_population(K=K,eta=eta,q=q)

    # initialize human population with some infections so that snail infections don't go to zero
    susceptibility_quantile <- quantile(Population_iter$Person_risk_beta,1-start_prev)
    Population_iter[Population_iter$Person_risk_beta>susceptibility_quantile,]$Adult_worm <- start_worms
    
    iter_list <- list(humans=Population_iter,snails=sn)
  }

  # Simulate n_burn_sim weeks of infection dynamics
  for (burn in 1:n_burn_sim) {
    if (snails){
      iter_list <- schisto_drug_sim(iter_list$humans,lambda_all, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,max_epg=max_epg,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration,maturation_0=maturation_0,maturation_alpha=maturation_alpha,sn=iter_list$snails,sn_seasonality=sn_seasonality)
    } else {
      Population_iter <- schisto_drug_sim(Population_iter, lambda_all, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,max_epg=max_epg,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration,seasonality=seasonality,t=start_t+burn,maturation_0=maturation_0,maturation_alpha=maturation_alpha)
    }
  }
  
  # Return baseline state at equilibrium
  if (snails){
    return(iter_list)
  } else {
    return(Population_iter)
  }
}

# Function to simulate mass drug administration and following weeks of reinfections and evaluations
strategy_sim <- function(Initial_state,lambda_all,Treatment_weeks,
                         Adult_efficacy,Juvenile_efficacy,
                         Time_outcome,delay_efficacy=0.5,
                         worm_death=1/208,
                         wormpair_EPG_conversion=5,
                         juvenile_duration=6,
                         maturation_0=NULL,maturation_alpha=NULL,N=500,
                         coverage=1, nonadherance=0, 
                         static_force=1,
                         binomial_treatment=FALSE,treatment_dispersion=NULL,
                         treatment_redraw=FALSE,
                         dynamic_egg='no_saturation',seasonality=0,
                         t_start=28,no_cap=FALSE,max_epg=5000,
                         infect_dispersion=NULL,nbinom_trunc=NULL,
                         followup1_week=5,followup2_week=11,dhp_like=FALSE,sn=NULL,sn_seasonality=FALSE){
  # Check that parameters are valid
  if ((1-nonadherance)<coverage){
    print("Coverage not possible with given level of nonadherance")
  } 
  if ((Adult_efficacy>1)|(Adult_efficacy<0)){
    print("Invalid adult efficacy")
  } 
  if ((Juvenile_efficacy>1)|(Juvenile_efficacy<0)){
    print("Invalid juvenile efficacy")
  }

  # Initialize
  n_treated <- floor(nrow(Initial_state)*coverage)
  Strategy <- Initial_state

  # Overdispersed treatment requires setting treatment efficacy in eacy individual according to beta distribution
  if (binomial_treatment&!(is.null(treatment_dispersion))){
    alpha <- Adult_efficacy/treatment_dispersion
    beta <- (1-Adult_efficacy)/treatment_dispersion
    ratio_juvenile <- Juvenile_efficacy/Adult_efficacy
    if (!treatment_redraw){
      Strategy["personal_efficacy"] <- rbeta(nrow(Initial_state),alpha,beta)
    }
  }

  # Initialize Strategy_plot, an object to keep track of certain variables over time 
  Strategy_plot <- matrix(0, nrow=13, ncol=(Time_outcome+1))
  Strategy_plot[1,1] <- sum(Initial_state$EPG>0)/length(Initial_state$EPG)*100 # Infection prevalence
  Strategy_plot[2,1] <- sum(Initial_state$EPG>0 & Initial_state$EPG<100)/sum(Initial_state$EPG>0)*100 # Proportion light infections
  Strategy_plot[3,1] <- sum(Initial_state$EPG>99 & Initial_state$EPG<400)/sum(Initial_state$EPG>0)*100 # Moderate
  Strategy_plot[4,1] <- sum(Initial_state$EPG>399)/sum(Initial_state$EPG>0)*100 # Heavy
  Strategy_plot[5,1] <- mean(Initial_state[Initial_state$EPG>0,]$EPG) # mean EPG
  Strategy_plot[6,1] <- exp(mean(log(Initial_state[Initial_state$EPG>0,]$EPG+1))) # GMI
  Strategy_plot[7,1] <- sum(Initial_state$Test_3dx2s>0)/length(Initial_state$Test_3dx2s)*100 # 3dx2s observed infection prevalence
  Strategy_plot[8,1] <- sum(Initial_state$Test_3dx2s>0 & Initial_state$EPG<100)/sum(Initial_state$Test_3dx2s>0)*100 # Light 3dx2s
  Strategy_plot[9,1] <- sum(Initial_state$Test_3dx2s>0 & (Initial_state$EPG>99 & Initial_state$EPG<400))/sum(Initial_state$Test_3dx2s>0)*100 #Moderate 3dx2s
  Strategy_plot[10,1] <- sum(Initial_state$Test_3dx2s>0 & Initial_state$EPG>399)/sum(Initial_state$Test_3dx2s>0)*100 # Heavy 3dx2s
  Strategy_plot[11,1] <- mean(Initial_state[Initial_state$Test_3dx2s>0,]$EPG) # mean EPG
  Strategy_plot[12,1] <- exp(mean(log(Initial_state[Initial_state$Test_3dx2s>0,]$EPG+1))) # GMI 3dx2s
  
  # Juvenile prevalence, depending on life stage duration:
  if (juvenile_duration==6){
      Strategy_plot[13,1] <- sum((Initial_state$Juvenile_worm1+Initial_state$Juvenile_worm2+Initial_state$Juvenile_worm3+Initial_state$Juvenile_worm4+Initial_state$Juvenile_worm5+Initial_state$Juvenile_worm6)>0)/length(Initial_state$EPG)*100 #Juvenile prevalence
  } else if (juvenile_duration==4){
      Strategy_plot[13,1] <- sum((Initial_state$Juvenile_worm1+Initial_state$Juvenile_worm2+Initial_state$Juvenile_worm3+Initial_state$Juvenile_worm4)>0)/length(Initial_state$EPG)*100 #Juvenile prevalence
  } else if (juvenile_duration==8){
      Strategy_plot[13,1] <- sum((Initial_state$Juvenile_worm1+Initial_state$Juvenile_worm2+Initial_state$Juvenile_worm3+Initial_state$Juvenile_worm4+Initial_state$Juvenile_worm5+Initial_state$Juvenile_worm6+Initial_state$Juvenile_worm7+Initial_state$Juvenile_worm8)>0)/length(Initial_state$EPG)*100 #Juvenile prevalence
  } else if ((juvenile_duration==10)|juvenile_duration=="variable"){
      Strategy_plot[13,1] <- sum((Initial_state$Juvenile_worm1+Initial_state$Juvenile_worm2+Initial_state$Juvenile_worm3+Initial_state$Juvenile_worm4+Initial_state$Juvenile_worm5+Initial_state$Juvenile_worm6+Initial_state$Juvenile_worm7+Initial_state$Juvenile_worm8+Initial_state$Juvenile_worm9+Initial_state$Juvenile_worm10)>0)/length(Initial_state$EPG)*100 #Juvenile prevalence
  }
  Strategy_plot[is.na(Strategy_plot[,1]),1] <- 0

  # Simulate enough time steps to cover treatment and followups
  for (sim_time in 1:Time_outcome) {
    # If in a follow-up week, assign states to these observations
    if (sim_time==followup1_week){state_followup1 <- data.frame(Strategy)}
    if (sim_time==followup2_week){state_followup2 <- data.frame(Strategy)}

    # Update dynamic force of infection according to the specified model
    if (dynamic_egg=='no_saturation'){
      # FOI depends on mean egg shedding in last 4–8 weeks
      delay_prev <- mean(c(Strategy_plot[5,max(c(1,sim_time-4))]*Strategy_plot[1,max(c(1,sim_time-4))],Strategy_plot[5,max(c(1,sim_time-5))]*Strategy_plot[1,max(c(1,sim_time-5))],Strategy_plot[5,max(c(1,sim_time-6))]*Strategy_plot[1,max(c(1,sim_time-6))],Strategy_plot[5,max(c(1,sim_time-7))]*Strategy_plot[1,max(c(1,sim_time-7))],Strategy_plot[5,max(c(1,sim_time-8))]*Strategy_plot[1,max(c(1,sim_time-8))]))
      base_prev <- Strategy_plot[5,1]*Strategy_plot[1,1]
    } else if (dynamic_egg=='saturation'){
      # FOI depends on non-linear saturation function of mean egg shedding in last 4–8 weeks
      base_prev <- Strategy_plot[5,1]*Strategy_plot[1,1]/100
      unsat_delay_prev <- mean(c(Strategy_plot[5,max(c(1,sim_time-4))]*Strategy_plot[1,max(c(1,sim_time-4))],Strategy_plot[5,max(c(1,sim_time-5))]*Strategy_plot[1,max(c(1,sim_time-5))],Strategy_plot[5,max(c(1,sim_time-6))]*Strategy_plot[1,max(c(1,sim_time-6))],Strategy_plot[5,max(c(1,sim_time-7))]*Strategy_plot[1,max(c(1,sim_time-7))],Strategy_plot[5,max(c(1,sim_time-8))]*Strategy_plot[1,max(c(1,sim_time-8))]))/100
      mag <- base_prev/(1-exp(-base_prev/10))
      delay_prev <- mag*(1 - exp(-unsat_delay_prev/10))
    } else {
      # FOI depends on mean prevalence in last 4–8 weeks
      delay_prev <- mean(c(Strategy_plot[1,max(c(1,sim_time-4))],Strategy_plot[1,max(c(1,sim_time-5))],Strategy_plot[1,max(c(1,sim_time-6))],Strategy_plot[1,max(c(1,sim_time-7))],Strategy_plot[1,max(c(1,sim_time-8))]))
      base_prev <- Strategy_plot[1,1]
    }
    # update FOI
    lambda_dyn <- lambda_all*(static_force + ((1-static_force) * (delay_prev/base_prev)))
    if (is.na(lambda_dyn)){
      lambda_dyn <- 0
    }

    # Simulate treatment, then infection
    if (sim_time %in% Treatment_weeks & coverage>0){
      # If the treatment time is not an exact number of years from the first
      # treatment, then it's the second round of treatment for that year, and
      # efficacy should be delay_efficacy. Otherwise it's just Adult_efficacy.
      efficacy <- Adult_efficacy
      if ((sim_time-min(Treatment_weeks))%%26!=0) {efficacy<-delay_efficacy}
      # Simulate treatment
      # Select random individuals to be treated according to coverage, and assign last rows as nonadherent individuals
      if ((sim_time-min(Treatment_weeks))%%26==0) {treat_rows <- sample(floor(nrow(Strategy)*(1-nonadherance)),n_treated)}
      Treated <- Strategy[treat_rows,]
      if (binomial_treatment){
        # Assuming some dispersion in treatment
        if (is.null(treatment_dispersion)){
          # Assuming treatment is not overdispersed, just use binomial
          for (i in 1:n_treated){
            # Reduce worm compartments according to efficacy
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
            if ((juvenile_duration>8)){
              Treated$Juvenile_worm9[i] <- rbinom(1,Treated$Juvenile_worm9[i],(1-Juvenile_efficacy))
              Treated$Juvenile_worm10[i] <- rbinom(1,Treated$Juvenile_worm10[i],(1-Juvenile_efficacy))}
          }
        } else {
          if (treatment_redraw){
            # Redraw means that efficacy in individuals is not stable over time, i.e. overdispersion is not
            # due to differences in physiology or behaviour between individuals.
            Strategy["personal_efficacy"] <- rbeta(nrow(Initial_state),alpha,beta)
          }
          for (i in 1:n_treated){
            # Overdispersed treatment is beta-binomial distributed
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
            if (juvenile_duration>8){
              Treated$Juvenile_worm9[i] <- rbinom(1,Treated$Juvenile_worm9[i],(1-personal_Juvenile_efficacy))
              Treated$Juvenile_worm10[i] <- rbinom(1,Treated$Juvenile_worm10[i],(1-personal_Juvenile_efficacy))
            }
          }
        }
      } else {
        # This assumes that there is no dispersion in treatment except from rounding
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
        if (juvenile_duration>8){
          Treated$Juvenile_worm9 <- floor(Treated$Juvenile_worm9*(1-Juvenile_efficacy))
          Treated$Juvenile_worm10 <- floor(Treated$Juvenile_worm10*(1-Juvenile_efficacy))
        }
      }
      for (r in 1:n_treated) {Strategy[treat_rows[r],] <- Treated[r,]}

      ##store information about what occurs at treatment
      state_treat <- data.frame(Strategy)
      for (i in 1:N) { # update epg, infection binary, and test results
        Pairs_i <- state_treat$Adult_worm[i]%/%2
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
        test_2dx2s <- KK_result(state_treat$EPG[i],2,2)
        # test_1dx1s <- KK_result(state_treat$EPG[i],1,1)
        # state_treat$Test_1dx2s[i] <- test_1dx2s
        state_treat$Test_3dx2s[i] <- test_3dx2s
        state_treat$Test_2dx2s[i] <- test_2dx2s
        # state_treat$Test_1dx1s[i] <- test_1dx1s
      }

      ## Infection step, including maturation, updating EPGs etc.
      # Optionally with treated individuals protected from infection
      if (dhp_like){
        skip_infect <- rep(c(FALSE),N)
        skip_infect[treat_rows] <- TRUE
      } else {
        skip_infect <- NULL
      }
      # Infection with snails, or no snails
      if (!is.null(sn)){
        slist <- schisto_drug_sim(Strategy, lambda_dyn, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,
        max_epg=max_epg,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration,
        maturation_0=maturation_0,maturation_alpha=maturation_alpha,
        skip_infect=skip_infect,
        sn=sn,sn_seasonality=sn_seasonality)
        Strategy <- slist$humans
        sn <- slist$snails
      } else {
        Strategy <- schisto_drug_sim(Strategy, lambda_dyn, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,
          max_epg=max_epg,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration,
          maturation_0=maturation_0,maturation_alpha=maturation_alpha,
          seasonality=seasonality,t=700+t_start+sim_time,
          skip_infect=skip_infect)
      }
    } else {
      # If no treatment just run the model normally for one time step
      if (!is.null(sn)){
        slist <- schisto_drug_sim(Strategy, lambda_dyn, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,
          max_epg=max_epg,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration,
          maturation_0=maturation_0,maturation_alpha=maturation_alpha,
          sn=sn,sn_seasonality=sn_seasonality)
        Strategy <- slist$humans
        sn <- slist$snails
      } else {
        Strategy <- schisto_drug_sim(Strategy, lambda_dyn, worm_death,wormpair_EPG_conversion,no_cap,infect_dispersion,
          max_epg=max_epg,nbinom_trunc=nbinom_trunc,juvenile_duration=juvenile_duration,
          seasonality=seasonality,t=700+t_start+sim_time,
          maturation_0=maturation_0,maturation_alpha=maturation_alpha)
      }
    }

    # Update Strategy_plot variables
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
    if (juvenile_duration==6){
      Strategy_plot[13,sim_time+1] <- sum((Strategy$Juvenile_worm1+Strategy$Juvenile_worm2+Strategy$Juvenile_worm3+Strategy$Juvenile_worm4+Strategy$Juvenile_worm5+Strategy$Juvenile_worm6)>0)/length(Strategy$EPG)*100 #Juvenile prevalence
    } else if (juvenile_duration==4){
      Strategy_plot[13,sim_time+1] <- sum((Strategy$Juvenile_worm1+Strategy$Juvenile_worm2+Strategy$Juvenile_worm3+Strategy$Juvenile_worm4)>0)/length(Strategy$EPG)*100 #Juvenile prevalence
    } else if (juvenile_duration==8){
      Strategy_plot[13,sim_time+1] <- sum((Strategy$Juvenile_worm1+Strategy$Juvenile_worm2+Strategy$Juvenile_worm3+Strategy$Juvenile_worm4+Strategy$Juvenile_worm5+Strategy$Juvenile_worm6+Strategy$Juvenile_worm7+Strategy$Juvenile_worm8)>0)/length(Strategy$EPG)*100 #Juvenile prevalence
    } else if (juvenile_duration==10|(juvenile_duration=="variable")){
      Strategy_plot[13,sim_time+1] <- sum((Strategy$Juvenile_worm1+Strategy$Juvenile_worm2+Strategy$Juvenile_worm3+Strategy$Juvenile_worm4+Strategy$Juvenile_worm5+Strategy$Juvenile_worm6+Strategy$Juvenile_worm7+Strategy$Juvenile_worm8+Strategy$Juvenile_worm9+Strategy$Juvenile_worm10)>0)/length(Strategy$EPG)*100 #Juvenile prevalence
    }
    Strategy_plot[is.na(Strategy_plot[,sim_time+1]),sim_time+1] <- 0  
  }
  
  # return important objects and variables
  return(list('final_state'=Strategy,'treated'=treat_rows,'state_treat'=state_treat,'state_followup1'=state_followup1,'state_followup2'=state_followup2,'plots'=Strategy_plot))
}

# Calculate a variety of relevant outcomes of treatment
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
  mean_epg_preschool_3dx2s <- mean(state[state$Age=="Preschool" & state$Test_3dx2s>0,]$EPG)
  mean_epg_school_3dx2s <- mean(state[state$Age=="School" & state$Test_3dx2s>0,]$EPG)
  mean_epg_adult_3dx2s <- mean(state[state$Age=="Adult" & state$Test_3dx2s>0,]$EPG)
  logepg <- log(state[state$Infection_binary>0,]$EPG)
  logepg_3dx2s <- log(state[state$Test_3dx2s>0,]$EPG)
  logepg_preschool_3dx2s <- log(state[state$Age=="Preschool" & state$Test_3dx2s>0,]$EPG)
  logepg_school_3dx2s <- log(state[state$Age=="School" & state$Test_3dx2s>0,]$EPG)
  logepg_adult_3dx2s <- log(state[state$Age=="Adult" & state$Test_3dx2s>0,]$EPG)
  gmean_epg <- exp(mean(logepg[is.finite(logepg)]))
  gmean_epg_3dx2s <- exp(mean(logepg_3dx2s[is.finite(logepg_3dx2s)]))
  gmean_epg_preschool_3dx2s <- exp(mean(logepg_preschool_3dx2s[is.finite(logepg_preschool_3dx2s)]))
  gmean_epg_school_3dx2s <- exp(mean(logepg_school_3dx2s[is.finite(logepg_school_3dx2s)]))
  gmean_epg_adult_3dx2s <- exp(mean(logepg_adult_3dx2s[is.finite(logepg_adult_3dx2s)]))
  base_prevalence <- round(sum(initial$EPG>0)/length(initial$EPG)*100,digits = 1)
  base_prevalence_3dx2s <- round(sum(initial$Test_3dx2s>0)/length(initial$EPG)*100,digits = 1)
  base_prevalence_preschool_3dx2s <- round(sum(initial[initial$Age=="Preschool",]$Test_3dx2s>0)/
    length(initial[initial$Age=="Preschool",]$EPG)*100,digits = 1)
  base_prevalence_school_3dx2s <- round(sum(initial[initial$Age=="School",]$Test_3dx2s>0)/
    length(initial[initial$Age=="School",]$EPG)*100,digits = 1)
  base_prevalence_adult_3dx2s <- round(sum(initial[initial$Age=="Adult",]$Test_3dx2s>0)/
    length(initial[initial$Age=="Adult",]$EPG)*100,digits = 1)
  prevalence <- round(sum(state$EPG>0)/length(state$EPG)*100,digits = 1)
  prevalence_3dx2s <- round(sum(state$Test_3dx2s>0)/length(state$Test_3dx2s)*100,digits = 1)
  prevalence_preschool_3dx2s <- round(sum(state[state$Age=="Preschool",]$Test_3dx2s>0)/
    length(state[state$Age=="Preschool",]$EPG)*100,digits = 1)
  prevalence_school_3dx2s <- round(sum(state[state$Age=="School",]$Test_3dx2s>0)/
    length(state[state$Age=="School",]$EPG)*100,digits = 1)
  prevalence_adult_3dx2s <- round(sum(state[state$Age=="Adult",]$Test_3dx2s>0)/
    length(initial[initial$Age=="Adult",]$EPG)*100,digits = 1)
  light <- round(sum(state$EPG>0 & state$EPG<100)/sum(state$EPG>0)*100)
  light_3dx2s <- round(sum(state$Test_3dx2s & state$EPG<100)/sum(state$Test_3dx2s)*100)
  light_preschool_3dx2s <- round(sum(state$Age=="Preschool" & (state$Test_3dx2s & state$EPG<100))/sum(state[state$Age=="Preschool",]$Test_3dx2s)*100)
  light_school_3dx2s <- round(sum(state$Age=="School" & (state$Test_3dx2s & state$EPG<100))/sum(state[state$Age=="School",]$Test_3dx2s)*100)
  light_adult_3dx2s <- round(sum(state$Age=="Adult" & (state$Test_3dx2s & state$EPG<100))/sum(state[state$Age=="Adult",]$Test_3dx2s)*100)
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
  initial_cured_preschool <- initial_treated[is_infected & initial_treated$Age=="Preschool",]
  initial_cured_school <- initial_treated[is_infected & initial_treated$Age=="School",]
  initial_cured_adult <- initial_treated[is_infected & initial_treated$Age=="Adult",]
  initial_cured_light <- initial_treated[is_light,]
  initial_cured_moderate <- initial_treated[is_moderate,]
  initial_cured_heavy <- initial_treated[is_heavy,]
  initial_cured_3dx2s <- initial_treated[is_infected_3dx2s,]
  initial_cured_preschool_3dx2s <- initial_treated[is_infected_3dx2s & initial_treated$Age=="Preschool",]
  initial_cured_school_3dx2s <- initial_treated[is_infected_3dx2s & initial_treated$Age=="School",]
  initial_cured_adult_3dx2s <- initial_treated[is_infected_3dx2s & initial_treated$Age=="Adult",]
  initial_cured_light_3dx2s <- initial_treated[is_light_3dx2s,]
  initial_cured_moderate_3dx2s <- initial_treated[is_moderate_3dx2s,]
  initial_cured_heavy_3dx2s <- initial_treated[is_heavy_3dx2s,]

  final_treated <- state[analysis$treated,]

  final_cured <- final_treated[is_infected,]
  final_cured_preschool <- final_treated[is_infected & final_treated$Age=="Preschool",]
  final_cured_school <- final_treated[is_infected & final_treated$Age=="School",]
  final_cured_adult <- final_treated[is_infected & final_treated$Age=="Adult",]
  final_cured_light <- final_treated[is_light,]
  final_cured_moderate <- final_treated[is_moderate,]
  final_cured_heavy <- final_treated[is_heavy,]
  final_cured_3dx2s <- final_treated[is_infected_3dx2s,]
  final_cured_preschool_3dx2s <- final_treated[is_infected_3dx2s & final_treated$Age=="Preschool",]
  final_cured_school_3dx2s <- final_treated[is_infected_3dx2s & final_treated$Age=="School",]
  final_cured_adult_3dx2s <- final_treated[is_infected_3dx2s & final_treated$Age=="Adult",]
  final_cured_light_3dx2s <- final_treated[is_light_3dx2s,]
  final_cured_moderate_3dx2s <- final_treated[is_moderate_3dx2s,]
  final_cured_heavy_3dx2s <- final_treated[is_heavy_3dx2s,]
  
  cure_rate <- (1 - (sum(final_cured$Infection_binary)/sum(initial_cured$Infection_binary)))*100
  cure_rate_preschool <- (1 - (sum(final_cured_preschool$Infection_binary)/sum(initial_cured_preschool$Infection_binary)))*100
  cure_rate_school <- (1 - (sum(final_cured_school$Infection_binary)/sum(initial_cured_school$Infection_binary)))*100
  cure_rate_adult <- (1 - (sum(final_cured_adult$Infection_binary)/sum(initial_cured_adult$Infection_binary)))*100
  cure_rate_light <- (1 - (sum(final_cured_light$Infection_binary)/sum(initial_cured_light$Infection_binary)))*100
  cure_rate_moderate <- (1 - (sum(final_cured_moderate$Infection_binary)/sum(initial_cured_moderate$Infection_binary)))*100
  cure_rate_heavy <- (1 - (sum(final_cured_heavy$Infection_binary)/sum(initial_cured_heavy$Infection_binary)))*100
  cure_rate_3dx2s <- (1 - (sum(final_cured_3dx2s$Test_3dx2s)/sum(initial_cured_3dx2s$Test_3dx2s)))*100
  cure_rate_preschool_3dx2s <- (1 - (sum(final_cured_preschool_3dx2s$Infection_binary)/sum(initial_cured_preschool_3dx2s$Infection_binary)))*100
  cure_rate_school_3dx2s <- (1 - (sum(final_cured_school_3dx2s$Infection_binary)/sum(initial_cured_school_3dx2s$Infection_binary)))*100
  cure_rate_adult_3dx2s <- (1 - (sum(final_cured_adult_3dx2s$Infection_binary)/sum(initial_cured_adult_3dx2s$Infection_binary)))*100
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

  return(list('base_prevalence'=base_prevalence,'base_prevalence_3dx2s'=base_prevalence_3dx2s,
    'base_prevalence_preschool_3dx2s'=base_prevalence_preschool_3dx2s,
    'base_prevalence_school_3dx2s'=base_prevalence_school_3dx2s,'base_prevalence_adult_3dx2s'=base_prevalence_adult_3dx2s,
    'drug'=drug,'coverage'=coverage,'n_adults'=n_adults,'n_juveniles'=n_juveniles,
    'mean_epg_3dx2s'=mean_epg_3dx2s,'gmean_epg_3dx2s'=gmean_epg_3dx2s,
    'prevalence_3dx2s'=prevalence_3dx2s,'light_3dx2s'=light_3dx2s,
    'prevalence_preschool_3dx2s'=prevalence_preschool_3dx2s,
    'prevalence_school_3dx2s'=prevalence_school_3dx2s,'prevalence_adult_3dx2s'=prevalence_adult_3dx2s,
    'cure_rate_3dx2s'=cure_rate_3dx2s,'cure_rate_light_3dx2s'=cure_rate_light_3dx2s,
    'cure_rate_moderate_3dx2s'=cure_rate_moderate_3dx2s,'cure_rate_heavy_3dx2s'=cure_rate_heavy_3dx2s,
    'cure_rate_preschool_3dx2s'=cure_rate_preschool_3dx2s,
    'cure_rate_school_3dx2s'=cure_rate_school_3dx2s,'cure_rate_adult_3dx2s'=cure_rate_adult_3dx2s,
    'mean_epg_reduction_3dx2s'=mean_epg_reduction_3dx2s,'gmean_epg_reduction_3dx2s'=gmean_epg_reduction_3dx2s,
    'juvenile_cured_3dx2s'=juvenile_cured_3dx2s,
    'mean_epg'=mean_epg,'gmean_epg'=gmean_epg,'prevalence'=prevalence,
    'intensity_prev'=paste(light,moderate,heavy,sep='/'),
    'cure_rate'=cure_rate,'cure_rate_light'=cure_rate_light,
    'cure_rate_moderate'=cure_rate_moderate,'cure_rate_heavy'=cure_rate_heavy,
    'mean_epg_reduction'=mean_epg_reduction,'gmean_epg_reduction'=gmean_epg_reduction,
    'base_cases_juveniles'=juvenile_ratio_baseline,'juvenile_cured'=juvenile_cured,'juvenile_prevalence'=juvenile_prevalence))
}
