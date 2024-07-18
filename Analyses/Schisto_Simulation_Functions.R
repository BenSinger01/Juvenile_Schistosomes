## Schistosomiasis Simulation Functions
## B J Singer
## Functions to simulate schistosomiasis infection, tracking juvenile schistosomes as they age

source("Analyses/KK_sensitivity.R")

# Fixed parameters in all versions of the model
worm_death <- 1/208 # Worms live four years on average, 208 weeks
max_epg <- 5000 # Saturation value of worm fecundity
nbinom_trunc <- 100 # Maximum number of new juveniles in one individual in one week

# Function to create population of uninfected individuals as base for simulations
create_population <- function(lambda_all, # Overall force of infection, 1 in main analysis
                              dispersion_param, mu_param, # Susceptibility distribution parameters, vary by setting
                              preschool_ratio = 0.61, adult_ratio = 0.82, # Relative susceptability of pre-school and adults to SAC
                              N=500 # Size of the human population
                              )
{
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
  Infection_binary <- rep(0,N) # Is this individual infected?
  Test <- rep(0,N) # Do they test positive?
  EPG <- rep(0,N) # How many eggs per gram of feces?

  # Age groups: preschool-age children, school-age children, adolescents/adults
  # Source: Lo NC, Lai YS, Karagiannis-Voules DA, Bogoch II, Coulibaly JT, Bendavid E, et al. Assessment of global guidelines for preventive chemotherapy against schistosomiasis and soil-transmitted helminthiasis: a cost-effectiveness modelling study. Lancet Infect Dis. 2016 Sep 1;16(9):1065–75.
  preschool_proportion <- 0.1
  school_proportion <- 0.25
  adult_proportion <- 0.65
  # Assign age group randomly, with (approximately) correct proportions
  Age <- sample(c('Preschool','School','Adult'),N,replace=TRUE,prob=c(preschool_proportion,school_proportion,adult_proportion))

  # Balance relative susceptability so that average susceptability is constant
  school_relative_beta <- 1/(school_proportion+preschool_proportion*preschool_ratio+adult_proportion*adult_ratio)
  preschool_relative_beta <- preschool_ratio*school_relative_beta
  adult_relative_beta <- adult_ratio*school_relative_beta

  # Calculate mu and dispersion so that mean of lognormal susceptiblity distributions are proportioned according to relative risk
  # These expressions are derived by rearranging the expression for the mean of the lognormal distribution
  mu_preschool <- log(preschool_relative_beta^2*exp(2*mu_param+dispersion_param^2)/sqrt((preschool_relative_beta^2+exp(dispersion_param^2)-1)*exp(2*mu_param+dispersion_param^2)))
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
  Population <- data.frame(Number, Person_risk_beta, Age, Juvenile_worm1, Juvenile_worm2, Juvenile_worm3, Juvenile_worm4, Juvenile_worm5, Juvenile_worm6, Juvenile_worm7, Juvenile_worm8, Juvenile_worm9, Juvenile_worm10, Adult_worm, Infection_binary, Test, EPG)
  
  Population
}

# Function to create snail compartments and keep track of snail model parameters
create_snail_population <- function(B_0 = 7/10, # Snail birth rate
                                    K_0 = 5, # Carrying capacity density
                                    mu = 7/100, # Snail death rate
                                    mu_I = 7/25, # Additional mortality for infected snails
                                    incubation = 1/4, # Rate of transfer from E to I compartment
                                    eta = 1.55e-7, # Force of infection for snails for each unit egg shedding
                                    q = 23 # Force of infection for humans for each infected snail
                                    )
{
  Snail_Population <- list(B_0=B_0,K_0=K_0,mu=mu,mu_I=mu_I,incubation=incubation,eta=eta,q=q,
                            # Initialize compartments to reasonable guesses at equilibrium values
                            S=K_0*(1-mu/B_0)*0.9,
                            E=K_0*(1-mu/B_0)*0.05,
                            I=K_0*(1-mu/B_0)*0.05
                            )
  Snail_Population
}

# Function to simulate on week's worth of infection and aging of schistosomes
increment_sim <- function(df, # Current state
                           lambda_all = 1, # Overall force of infection
                           wormpair_EPG_conversion = 4, # Average baseline fecundity of worm pairs
                           infect_dispersion = 0.0175, # Degree of overdispersion in new infection intensity
                           juvenile_duration = 6, # How long do juveniles take to mature?
                           seasonality = 0, # How much does infection risk change seasonally? 0 to 1
                           t = 0, # What week of the simulation is it?
                           skip_infect = NULL, # Skip the infection step? e.g. DHP treatment
                           sn = NULL, # Snail state
                           test_samples = 3, # Number of Kato-Katz samples in test
                           test_slides = 2 # Number of Kato-Katz slides per sample
                           )
{ 
  N <- dim(df)[1]

  # if a snail state has been passed, update according to snail model
  if (!is.null(sn)){
    # update snail compartments according to equations (see SI pp7–8 for explanation of the model)
    N_sn <- sn$S + sn$E + sn$I
    egg_shedding <- sum(df$EPG)
    new_S <- sn$S + sn$B_0*(1-N_sn/sn$K_0)*(sn$S+sn$E) -
     sn$eta*egg_shedding*sn$S - sn$mu*sn$S
    new_E <- sn$E + sn$eta*egg_shedding*sn$S -
     sn$incubation*sn$E - sn$mu*sn$E
    new_I <- sn$I + sn$incubation*sn$E - (sn$mu+sn$mu_I)*sn$I

    sn$S <- new_S
    sn$E <- new_E
    sn$I <- new_I

    # update FOI proportionally to infected snails (proxy of free cercariae)
    lambda_all <- sn$q*sn$I*lambda_all
  }

  # Iterate over individual humans
  for (i in 1:N) {
    # Step 1 (Worm aging step)
    # From oldest to youngest compartment, age juveniles into next category
    Adult_worm_aging <- df[i,paste("Juvenile_worm",as.character(juvenile_duration),sep='')]
    for (jc in (juvenile_duration-1):1){
      df[i,paste("Juvenile_worm",as.character(jc+1),sep='')] <- df[i,paste("Juvenile_worm",as.character(jc),sep='')]
    }
    
    # Step 2 (Worm death step)
    Adult_worm_death <- rbinom(size=df$Adult_worm[i],n=1, prob=worm_death)
    
    # Step 3 (Infection step)
    # define force of infection, with optional seasonality
    lambda_i <- lambda_all * df$Person_risk_beta[i] * (1 + seasonality*cos(2*pi*t/52))

    # Infection step, possible to skip for certain treatment types
    # Null or False skip_infect do same thing
    if (is.null(skip_infect)){ # Null skip_infect
      # overdispersed infection, with truncation to prevent unrealistic extreme values
      Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
      # while loop implements truncation
      if (!(is.null(nbinom_trunc))){
        while (max(Juvenile_worm_infect)>nbinom_trunc){
          Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
        }
      }
    } else if (!skip_infect[i]) { # False skip_infect for this individual
      # Same infection step as above
      Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
      if (!(is.null(nbinom_trunc))){
        while (max(Juvenile_worm_infect)>nbinom_trunc){
          Juvenile_worm_infect <- rnbinom(1,mu=lambda_i,size=infect_dispersion)
        }
      }
    } else {
      # if infection step is skipped in this time step (ie DHP-like treatment), don't infect
      Juvenile_worm_infect <- 0
    }
    
    # Update adults and new juveniles
    df$Juvenile_worm1[i] <- Juvenile_worm_infect 
    df$Adult_worm[i] <- df$Adult_worm[i] + Adult_worm_aging - Adult_worm_death

    # Step 4 (Compute EPG)
    # Calculation of EPG based on density-dependent fecundity,
    # each worm pair generates fewer eggs if there are many worm pairs,
    # such that the maximum EPG is max_epg.
    Pairs_i <- df$Adult_worm[i]%/%2
    df$EPG[i] <- round(max_epg*(1-exp(-wormpair_EPG_conversion*Pairs_i/max_epg)))

    # # Step 5 (Define infection status)
    df$Infection_binary[i] <- ifelse(df$EPG[i]>0, 1, 0)

    # Step 6 (Test positivity)
    # Simulates results of two different testing regimes. Each uses two
    # slides per stool. One uses two stools, the other uses three stools
    # sampled on separate days. The former is used for Mnkugwe et al., the
    # latter for SCORE data. The KK_result function is defined in KK_sensitivity.R.
    df$Test[i] <- KK_result(df$EPG[i],test_slides,test_samples)
  }

  # Return single dataframe if no snails, or list if snails included
  if (is.null(sn)){
    df
  } else {
    outlist <- list("humans" = df, "snails" = sn)
    outlist
  }

}

# Function to burn in baseline state
burn_in_sim <- function(lambda_all, # Overall force of infection
                        dispersion_param, mu_param, # Susceptibility distribution parameters, vary by setting
                        wormpair_EPG_conversion = 4, # Average baseline fecundity of worm pairs
                        n_burn_sim = 700, # Number of weeks to run simulation
                        preschool_ratio = 0.61, adult_ratio = 0.82, # Relative susceptability of pre-school and adults to SAC
                        N = 500, # Size of population
                        infect_dispersion = 0.0175, # Degree of overdispersion in new infection intensity
                        juvenile_duration = 6, # How long do juveniles take to mature?
                        seasonality = 0, # How much does infection risk change seasonally? 0 to 1
                        start_t = 28, # At what week does the simulation start, i.e. what point on the cosine curve
                        snails = FALSE, # Include the snail model?
                        K=5,eta=1.55e-7,q=23, # Parameters for the snail model
                        start_prev=0.63,start_worms=13 # Initial values for the snail model
                        )
{
  # Create initial population
  Population_iter <- create_population(lambda_all,dispersion_param,mu_param,preschool_ratio=preschool_ratio,adult_ratio=adult_ratio,N=N)
  
  # If snails are to be included, create a snail population
  if (snails){
    sn <- create_snail_population(K=K,eta=eta,q=q)

    # Initialize human population with some infections so that snail infections don't go to zero
    # Each infected human has startm_worms adult worms
    susceptibility_quantile <- quantile(Population_iter$Person_risk_beta,1-start_prev)
    Population_iter[Population_iter$Person_risk_beta>susceptibility_quantile,]$Adult_worm <- start_worms
    
    iter_list <- list(humans=Population_iter,snails=sn)
  }

  # Simulate n_burn_sim weeks of infection dynamics
  for (burn in 1:n_burn_sim) {
    if (snails){
      iter_list <- increment_sim(iter_list$humans,lambda_all,wormpair_EPG_conversion,infect_dispersion,juvenile_duration,sn=iter_list$snails)
    } else {
      Population_iter <- increment_sim(Population_iter,lambda_all,wormpair_EPG_conversion,infect_dispersion,juvenile_duration,seasonality=seasonality,t=start_t+burn)
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
strategy_sim <- function(Initial_state, # Baseline state
                         lambda_all, # Overall force of infection
                         Treatment_weeks, # Weeks from simulation start in which treatment occurs
                         Adult_efficacy, # Average proportion of adult schistosomes killed by treatment
                         Juvenile_efficacy, # Average proportion of juvenile schistosomes killed by treatment
                         Time_outcome, # Time of final follow-up
                         wormpair_EPG_conversion = 4, # Average baseline fecundity of worm pairs
                         juvenile_duration = 6, # How long do juveniles take to mature?
                         infect_dispersion = 0.0175, # Degree of overdispersion in new infection intensity
                         treatment_dispersion = 0.49, # Overdispersion of schistosome death in treatment
                         N = 500, # Size of the human population
                         coverage = 1, # Proportion of the population treated
                         nonadherance = 0, # Proportion of the population systematically nonadherent
                         static_force = 1, # Proportion of force of infection that does not change with treatment
                         dynamic_egg = 'no_saturation', # How does force of infection depend on current infections?
                         seasonality = 0, # How much does infection risk change seasonally? 0 to 1
                         t_start = 28, # At what week does the simulation start, i.e. what point on the cosine curve
                         followup1_week = 5, followup2_week = 11, # When are initial outcomes observed?
                         dhp_like = FALSE, # Is the infection step skipped after treatment?
                         sn=NULL, # Snail state
                         test_samples = 3, # Number of Kato-Katz samples in test
                         test_slides = 2 # Number of Kato-Katz slides per sample
                         )
{
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

  # Treatment is oversdispersed, so each individual has an individual treatment efficacy drawn from a beta distribution
  # Beta distribution parameters:
  alpha <- Adult_efficacy/treatment_dispersion
  beta <- (1-Adult_efficacy)/treatment_dispersion
  alpha_juve <- Juvenile_efficacy/treatment_dispersion
  beta_juve <- (1-Juvenile_efficacy)/treatment_dispersion
  # generate uniform random numbers
  percentiles <- runif(nrow(Initial_state),0,1)
  # draw adult efficacy from beta distribution
  adult_ind_efficacy <- qbeta(percentiles,alpha,beta)
  Strategy["personal_efficacy"] <- adult_ind_efficacy
  # set juvenile efficacy to equivalent percentile in distribution with different mean
  # thus, adult and juvenile efficacy scale together but each fits a beta distribution
  juvenile_ind_efficacy <- qbeta(percentiles,alpha_juve,beta_juve)
  Strategy["juvenile_personal_efficacy"] <- juvenile_ind_efficacy

  # Initialize Strategy_plot, an object to keep track of certain variables over time 
  Strategy_plot <- matrix(0, nrow=15, ncol=(Time_outcome+1))
  Strategy_plot[1,1] <- sum(Initial_state$EPG>0)/length(Initial_state$EPG)*100 # Infection prevalence
  Strategy_plot[2,1] <- sum(Initial_state$EPG>0 & Initial_state$EPG<100)/sum(Initial_state$EPG>0)*100 # Proportion light infections
  Strategy_plot[3,1] <- sum(Initial_state$EPG>99 & Initial_state$EPG<400)/sum(Initial_state$EPG>0)*100 # Moderate
  Strategy_plot[4,1] <- sum(Initial_state$EPG>399)/sum(Initial_state$EPG>0)*100 # Heavy
  Strategy_plot[5,1] <- mean(Initial_state[Initial_state$EPG>0,]$EPG) # mean EPG
  Strategy_plot[6,1] <- exp(mean(log(Initial_state[Initial_state$EPG>0,]$EPG+1))) # GMI
  Strategy_plot[7,1] <- sum(Initial_state$Test>0)/length(Initial_state$Test)*100 # observed infection prevalence
  Strategy_plot[8,1] <- sum(Initial_state$Test>0 & Initial_state$EPG<100)/sum(Initial_state$Test>0)*100 # Light observed
  Strategy_plot[9,1] <- sum(Initial_state$Test>0 & (Initial_state$EPG>99 & Initial_state$EPG<400))/sum(Initial_state$Test>0)*100 #Moderate observed
  Strategy_plot[10,1] <- sum(Initial_state$Test>0 & Initial_state$EPG>399)/sum(Initial_state$Test>0)*100 # Heavy observed
  Strategy_plot[11,1] <- mean(Initial_state[Initial_state$Test>0,]$EPG) # mean EPG observed
  Strategy_plot[12,1] <- exp(mean(log(Initial_state[Initial_state$Test>0,]$EPG+1))) # GMI observed
  # Number of juveniles:
  for (jc in 1:juvenile_duration){ # iterate over juvenile compartments
    Strategy_plot[14,1] <- Strategy_plot[14,1] + sum(Initial_state[,paste("Juvenile_worm",as.character(jc),sep='')]) # add to total
  }
  Strategy_plot[13,1] <- Strategy_plot[14,1]/nrow(Initial_state)*100 # Juvenile prevalence
  Strategy_plot[15,1] <- sum(Initial_state$Adult_worm) #No. adults
  Strategy_plot[is.na(Strategy_plot[,1]),1] <- 0 # Replace NA prevalence with zero

  # Simulate enough time steps to cover treatment and followups
  for (sim_time in 1:Time_outcome) {
    # If in a follow-up week, assign states to these observations
    if (sim_time==followup1_week){state_followup1 <- data.frame(Strategy)}
    if (sim_time==followup2_week){state_followup2 <- data.frame(Strategy)}

    # Update dynamic force of infection according to the specified model
    if (dynamic_egg=='no_saturation'){
      # FOI depends on mean egg shedding in last 4–8 weeks
      # Multiply mean epg of infections by number of infections to get total egg shedding
      delay_prev <- mean(c(Strategy_plot[5,max(c(1,sim_time-4))]*Strategy_plot[1,max(c(1,sim_time-4))],Strategy_plot[5,max(c(1,sim_time-5))]*Strategy_plot[1,max(c(1,sim_time-5))],Strategy_plot[5,max(c(1,sim_time-6))]*Strategy_plot[1,max(c(1,sim_time-6))],Strategy_plot[5,max(c(1,sim_time-7))]*Strategy_plot[1,max(c(1,sim_time-7))],Strategy_plot[5,max(c(1,sim_time-8))]*Strategy_plot[1,max(c(1,sim_time-8))]))
      # Initial value scales dynamic model
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
      # Simulate treatment
      # Select random individuals to be treated according to coverage, and assign last rows as nonadherent individuals
      # If multiple treatments in a round, only select for first treatment then keep treating same people until next round
      # Rounds are assumed to start on the first week of the year (52 weeks), or half way through the year
      if ((sim_time-min(Treatment_weeks))%%26==0) {treat_rows <- sample(floor(nrow(Strategy)*(1-nonadherance)),n_treated)}
      Treated <- Strategy[treat_rows,]
      for (i in 1:n_treated){
        # Overdispersed treatment is beta-binomial distributed, i.e. mean worm reduction parameter in the binomial comes from the beta distribution defined above
        # Reduce size of each worm compartment randomly according to binomial distribution with beta-distributed parameter
        Treated$Adult_worm[i] <- rbinom(1,Treated$Adult_worm[i],(1-Strategy$personal_efficacy[i]))
        for (jc in 1:juvenile_duration){ # iterate over juvenile compartments
          Treated[i,paste("Juvenile_worm",as.character(jc),sep='')] <- rbinom(1,Treated[i,paste("Juvenile_worm",as.character(jc),sep='')],(1-Strategy$juvenile_personal_efficacy[i]))
        }
      }
      for (r in 1:n_treated) {Strategy[treat_rows[r],] <- Treated[r,]} # Update state

      ##store information about what occurs at treatment
      state_treat <- data.frame(Strategy)
      for (i in 1:N) { # update epg, infection binary, and test results
        Pairs_i <- state_treat$Adult_worm[i]%/%2
        state_treat$EPG[i] <- round(max_epg*(1-exp(-wormpair_EPG_conversion*Pairs_i/max_epg)))
        state_treat$Infection_binary[i] <- ifelse(state_treat$EPG[i]>0, 1, 0)
        state_treat$Test[i] <- KK_result(state_treat$EPG[i],test_slides,test_samples)
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
        slist <- increment_sim(Strategy,lambda_dyn,wormpair_EPG_conversion,infect_dispersion,juvenile_duration,
          skip_infect=skip_infect,sn=sn,test_samples=test_samples,test_slides=test_slides)
        Strategy <- slist$humans
        sn <- slist$snails
      } else {
        Strategy <- increment_sim(Strategy,lambda_dyn,wormpair_EPG_conversion,infect_dispersion,juvenile_duration,seasonality,700+t_start+sim_time,skip_infect,
          test_samples=test_samples,test_slides=test_slides)
      }
    } else {
      # If no treatment just run the model normally for one time step
      if (!is.null(sn)){
        slist <- increment_sim(Strategy,lambda_dyn,wormpair_EPG_conversion,infect_dispersion,juvenile_duration,
          sn=sn,test_samples=test_samples,test_slides=test_slides)
        Strategy <- slist$humans
        sn <- slist$snails
      } else {
        Strategy <- increment_sim(Strategy,lambda_dyn,wormpair_EPG_conversion,infect_dispersion,juvenile_duration,seasonality,700+t_start+sim_time,
          test_samples=test_samples,test_slides=test_slides)
      }
    }

    # Update Strategy_plot variables
    Strategy_plot[1,sim_time+1] <- sum(Strategy$EPG>0)/length(Strategy$EPG)*100 # Any
    Strategy_plot[2,sim_time+1] <- sum(Strategy$EPG>0 & Strategy$EPG<100)/sum(Strategy$EPG>0)*100 # Light
    Strategy_plot[3,sim_time+1] <- sum(Strategy$EPG>99 & Strategy$EPG<400)/sum(Strategy$EPG>0)*100 #Moderate
    Strategy_plot[4,sim_time+1] <- sum(Strategy$EPG>399)/sum(Strategy$EPG>0)*100 # Heavy
    Strategy_plot[5,sim_time+1] <- mean(Strategy[Strategy$EPG>0,]$EPG) # mean EPG
    Strategy_plot[6,sim_time+1] <- exp(mean(log(Strategy[Strategy$EPG>0,]$EPG+1))) # GMI
    Strategy_plot[7,sim_time+1] <- sum(Strategy$Test>0)/length(Strategy$Test)*100 # Any observed
    Strategy_plot[8,sim_time+1] <- sum(Strategy$Test>0 & Strategy$EPG<100)/sum(Strategy$Test>0)*100 # Light observed
    Strategy_plot[9,sim_time+1] <- sum(Strategy$Test>0 & (Strategy$EPG>99 & Strategy$EPG<400))/sum(Strategy$Test>0)*100 #Moderate observed
    Strategy_plot[10,sim_time+1] <- sum(Strategy$Test>0 & Strategy$EPG>399)/sum(Strategy$Test>0)*100 # Heavy observed
    Strategy_plot[11,sim_time+1] <- mean(Strategy[Strategy$Test>0,]$EPG) # mean EPG observed
    Strategy_plot[12,sim_time+1] <- exp(mean(log(Strategy[Strategy$Test>0,]$EPG+1))) # GMI observed
    # Juvenile prevalence:
    juve_presence_by_jc <- matrix(NA,nrow=N,ncol=juvenile_duration)
    juve_presence <- 0
    # Number of juveniles:
    for (jc in 1:juvenile_duration){ # iterate over juvenile compartments
      juve_presence_by_jc[,jc] <- Strategy[,paste("Juvenile_worm",as.character(jc),sep='')]>0
      Strategy_plot[14,sim_time+1] <- Strategy_plot[14,sim_time+1] + sum(Strategy[,paste("Juvenile_worm",as.character(jc),sep='')]) # add to total
    }
    for (individual in 1:N){
      juve_presence <- juve_presence + any(juve_presence_by_jc[individual,])
    }
    Strategy_plot[13,sim_time+1] <- juve_presence/nrow(Strategy)*100 # Juvenile prevalence
    Strategy_plot[15,sim_time+1] <- sum(Strategy$Adult_worm) #No. adults
    Strategy_plot[is.na(Strategy_plot[,sim_time+1]),sim_time+1] <- 0 # Replace NA prevalence with zero
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
  mean_epg_test <- mean(state[state$Test>0,]$EPG)
  mean_epg_preschool_test <- mean(state[state$Age=="Preschool" & state$Test>0,]$EPG)
  mean_epg_school_test <- mean(state[state$Age=="School" & state$Test>0,]$EPG)
  mean_epg_adult_test <- mean(state[state$Age=="Adult" & state$Test>0,]$EPG)
  logepg <- log(state[state$Infection_binary>0,]$EPG+1)
  logepg_test <- log(state[state$Test>0,]$EPG+1)
  logepg_preschool_test <- log(state[state$Age=="Preschool" & state$Test>0,]$EPG+1)
  logepg_school_test <- log(state[state$Age=="School" & state$Test>0,]$EPG+1)
  logepg_adult_test <- log(state[state$Age=="Adult" & state$Test>0,]$EPG+1)
  gmean_epg <- exp(mean(logepg))
  gmean_epg_test <- exp(mean(logepg_test))
  gmean_epg_preschool_test <- exp(mean(logepg_preschool_test))
  gmean_epg_school_test <- exp(mean(logepg_school_test))
  gmean_epg_adult_test <- exp(mean(logepg_adult_test))
  base_prevalence <- round(sum(initial$EPG>0)/length(initial$EPG)*100,digits = 1)
  base_prevalence_test <- round(sum(initial$Test>0)/length(initial$EPG)*100,digits = 1)
  base_prevalence_preschool_test <- round(sum(initial[initial$Age=="Preschool",]$Test>0)/
    length(initial[initial$Age=="Preschool",]$EPG)*100,digits = 1)
  base_prevalence_school_test <- round(sum(initial[initial$Age=="School",]$Test>0)/
    length(initial[initial$Age=="School",]$EPG)*100,digits = 1)
  base_prevalence_adult_test <- round(sum(initial[initial$Age=="Adult",]$Test>0)/
    length(initial[initial$Age=="Adult",]$EPG)*100,digits = 1)
  prevalence <- round(sum(state$EPG>0)/length(state$EPG)*100,digits = 1)
  prevalence_test <- round(sum(state$Test>0)/length(state$Test)*100,digits = 1)
  prevalence_preschool_test <- round(sum(state[state$Age=="Preschool",]$Test>0)/
    length(state[state$Age=="Preschool",]$EPG)*100,digits = 1)
  prevalence_school_test <- round(sum(state[state$Age=="School",]$Test>0)/
    length(state[state$Age=="School",]$EPG)*100,digits = 1)
  prevalence_adult_test <- round(sum(state[state$Age=="Adult",]$Test>0)/
    length(initial[initial$Age=="Adult",]$EPG)*100,digits = 1)
  light <- round(sum(state$EPG>0 & state$EPG<100)/sum(state$EPG>0)*100)
  light_test <- round(sum(state$Test & state$EPG<100)/sum(state$Test)*100)
  light_preschool_test <- round(sum(state$Age=="Preschool" & (state$Test & state$EPG<100))/sum(state[state$Age=="Preschool",]$Test)*100)
  light_school_test <- round(sum(state$Age=="School" & (state$Test & state$EPG<100))/sum(state[state$Age=="School",]$Test)*100)
  light_adult_test <- round(sum(state$Age=="Adult" & (state$Test & state$EPG<100))/sum(state[state$Age=="Adult",]$Test)*100)
  moderate <- round(sum(state$EPG>=100 & state$EPG<400)/sum(state$EPG>0)*100)
  heavy <- round(sum(state$EPG>=400)/sum(state$EPG>0)*100)
  
  initial_treated <- initial[analysis$treated,]
  is_infected <- initial_treated$Infection_binary==1
  is_light <- is_infected & initial_treated$EPG<100
  is_moderate <- is_infected & initial_treated$EPG>=100 & is_infected & initial_treated$EPG<400
  is_heavy <- is_infected & initial_treated$EPG>400
  is_infected_test <- initial_treated$Test==1
  is_light_test <- is_infected_test & initial_treated$EPG<100
  is_moderate_test <- is_infected_test & initial_treated$EPG>=100 & initial_treated$EPG<400
  is_heavy_test <- is_infected_test & initial_treated$EPG>400
  
  initial_cured <- initial_treated[is_infected,]
  initial_cured_preschool <- initial_treated[is_infected & initial_treated$Age=="Preschool",]
  initial_cured_school <- initial_treated[is_infected & initial_treated$Age=="School",]
  initial_cured_adult <- initial_treated[is_infected & initial_treated$Age=="Adult",]
  initial_cured_light <- initial_treated[is_light,]
  initial_cured_moderate <- initial_treated[is_moderate,]
  initial_cured_heavy <- initial_treated[is_heavy,]
  initial_cured_test <- initial_treated[is_infected_test,]
  initial_cured_preschool_test <- initial_treated[is_infected_test & initial_treated$Age=="Preschool",]
  initial_cured_school_test <- initial_treated[is_infected_test & initial_treated$Age=="School",]
  initial_cured_adult_test <- initial_treated[is_infected_test & initial_treated$Age=="Adult",]
  initial_cured_light_test <- initial_treated[is_light_test,]
  initial_cured_moderate_test <- initial_treated[is_moderate_test,]
  initial_cured_heavy_test <- initial_treated[is_heavy_test,]

  final_treated <- state[analysis$treated,]

  final_cured <- final_treated[is_infected,]
  final_cured_preschool <- final_treated[is_infected & final_treated$Age=="Preschool",]
  final_cured_school <- final_treated[is_infected & final_treated$Age=="School",]
  final_cured_adult <- final_treated[is_infected & final_treated$Age=="Adult",]
  final_cured_light <- final_treated[is_light,]
  final_cured_moderate <- final_treated[is_moderate,]
  final_cured_heavy <- final_treated[is_heavy,]
  final_cured_test <- final_treated[is_infected_test,]
  final_cured_preschool_test <- final_treated[is_infected_test & final_treated$Age=="Preschool",]
  final_cured_school_test <- final_treated[is_infected_test & final_treated$Age=="School",]
  final_cured_adult_test <- final_treated[is_infected_test & final_treated$Age=="Adult",]
  final_cured_light_test <- final_treated[is_light_test,]
  final_cured_moderate_test <- final_treated[is_moderate_test,]
  final_cured_heavy_test <- final_treated[is_heavy_test,]
  
  cure_rate <- (1 - (sum(final_cured$Infection_binary)/sum(initial_cured$Infection_binary)))*100
  cure_rate_preschool <- (1 - (sum(final_cured_preschool$Infection_binary)/sum(initial_cured_preschool$Infection_binary)))*100
  cure_rate_school <- (1 - (sum(final_cured_school$Infection_binary)/sum(initial_cured_school$Infection_binary)))*100
  cure_rate_adult <- (1 - (sum(final_cured_adult$Infection_binary)/sum(initial_cured_adult$Infection_binary)))*100
  cure_rate_light <- (1 - (sum(final_cured_light$Infection_binary)/sum(initial_cured_light$Infection_binary)))*100
  cure_rate_moderate <- (1 - (sum(final_cured_moderate$Infection_binary)/sum(initial_cured_moderate$Infection_binary)))*100
  cure_rate_heavy <- (1 - (sum(final_cured_heavy$Infection_binary)/sum(initial_cured_heavy$Infection_binary)))*100
  cure_rate_test <- (1 - (sum(final_cured_test$Test)/sum(initial_cured_test$Test)))*100
  cure_rate_preschool_test <- (1 - (sum(final_cured_preschool_test$Infection_binary)/sum(initial_cured_preschool_test$Infection_binary)))*100
  cure_rate_school_test <- (1 - (sum(final_cured_school_test$Infection_binary)/sum(initial_cured_school_test$Infection_binary)))*100
  cure_rate_adult_test <- (1 - (sum(final_cured_adult_test$Infection_binary)/sum(initial_cured_adult_test$Infection_binary)))*100
  cure_rate_light_test <- (1 - (sum(final_cured_light_test$Test)/sum(initial_cured_light_test$Test)))*100
  cure_rate_moderate_test <- (1 - (sum(final_cured_moderate_test$Test)/sum(initial_cured_moderate_test$Test)))*100
  cure_rate_heavy_test <- (1 - (sum(final_cured_heavy_test$Test)/sum(initial_cured_heavy_test$Test)))*100
  epg_red <- (initial_cured$EPG - final_cured$EPG)/initial_cured$EPG
  epg_red_test <- (initial_cured_test$EPG - final_cured_test$EPG)/initial_cured_test$EPG
  mean_epg_reduction <- mean(epg_red)*100
  mean_epg_reduction_test <- mean(epg_red_test)*100
  logepg_red <- log(epg_red)
  logepg_red_test <- log(epg_red_test)
  gmean_epg_reduction <- exp(mean(logepg_red[is.finite(logepg_red)]))*100
  gmean_epg_reduction_test <- exp(mean(logepg_red_test[is.finite(logepg_red_test)]))*100
  
  is_cured <- final_cured$Infection_binary==0
  is_cured_test <- final_cured$Test==0
  treatment_cured <- analysis$state_treat[is_cured,]
  treatment_cured_test <- analysis$state_treat[is_cured_test,]
  survived_juveniles <- treatment_cured[(treatment_cured$Juvenile_worm1 + treatment_cured$Juvenile_worm2 +treatment_cured$Juvenile_worm3 +treatment_cured$Juvenile_worm4 +treatment_cured$Juvenile_worm5 +treatment_cured$Juvenile_worm6) > 0,]
  survived_juveniles_test <- treatment_cured_test[(treatment_cured_test$Juvenile_worm1 + treatment_cured_test$Juvenile_worm2 +treatment_cured_test$Juvenile_worm3 +treatment_cured_test$Juvenile_worm4 +treatment_cured_test$Juvenile_worm5 +treatment_cured_test$Juvenile_worm6) > 0,]
  juvenile_cured <- 100*nrow(survived_juveniles)/nrow(treatment_cured)
  juvenile_cured_test <- 100*nrow(survived_juveniles_test)/nrow(treatment_cured_test)

  baseline_juveniles <- initial$Juvenile_worm1 + initial$Juvenile_worm2 + initial$Juvenile_worm3 + initial$Juvenile_worm4 + initial$Juvenile_worm5 + initial$Juvenile_worm6
  juvenile_infections <- initial[baseline_juveniles > 0,]
  juvenile_ratio_baseline <- 100*sum(juvenile_infections$Infection_binary == 1) / sum(initial$Infection_binary == 1)
  
  final_juveniles <- state$Juvenile_worm1 + state$Juvenile_worm2 + state$Juvenile_worm3 + state$Juvenile_worm4 + state$Juvenile_worm5 + state$Juvenile_worm6
  juvenile_prevalence <- 100*sum(final_juveniles>0)/length(final_juveniles)

  return(list('base_prevalence'=base_prevalence,'base_prevalence_test'=base_prevalence_test,
    'base_prevalence_preschool_test'=base_prevalence_preschool_test,
    'base_prevalence_school_test'=base_prevalence_school_test,'base_prevalence_adult_test'=base_prevalence_adult_test,
    'drug'=drug,'coverage'=coverage,'n_adults'=n_adults,'n_juveniles'=n_juveniles,
    'mean_epg_test'=mean_epg_test,'gmean_epg_test'=gmean_epg_test,
    'mean_epg_preschool_test'=mean_epg_preschool_test,
    'mean_epg_school_test'=mean_epg_school_test,'mean_epg_adult_test'=mean_epg_adult_test,
    'gmean_epg_preschool_test'=gmean_epg_preschool_test,
    'gmean_epg_school_test'=gmean_epg_school_test,'gmean_epg_adult_test'=gmean_epg_adult_test,
    'prevalence_test'=prevalence_test,'light_test'=light_test,
    'light_preschool_test'=light_preschool_test,
    'light_school_test'=light_school_test,'light_adult_test'=light_adult_test,
    'prevalence_preschool_test'=prevalence_preschool_test,
    'prevalence_school_test'=prevalence_school_test,'prevalence_adult_test'=prevalence_adult_test,
    'cure_rate_test'=cure_rate_test,'cure_rate_light_test'=cure_rate_light_test,
    'cure_rate_moderate_test'=cure_rate_moderate_test,'cure_rate_heavy_test'=cure_rate_heavy_test,
    'cure_rate_preschool_test'=cure_rate_preschool_test,
    'cure_rate_school_test'=cure_rate_school_test,'cure_rate_adult_test'=cure_rate_adult_test,
    'mean_epg_reduction_test'=mean_epg_reduction_test,'gmean_epg_reduction_test'=gmean_epg_reduction_test,
    'juvenile_cured_test'=juvenile_cured_test,
    'mean_epg'=mean_epg,'gmean_epg'=gmean_epg,'prevalence'=prevalence,
    'intensity_prev'=paste(light,moderate,heavy,sep='/'),
    'cure_rate'=cure_rate,'cure_rate_light'=cure_rate_light,
    'cure_rate_moderate'=cure_rate_moderate,'cure_rate_heavy'=cure_rate_heavy,
    'cure_rate_preschool'=cure_rate_preschool,
    'cure_rate_school'=cure_rate_school,'cure_rate_adult'=cure_rate_adult,
    'mean_epg_reduction'=mean_epg_reduction,'gmean_epg_reduction'=gmean_epg_reduction,
    'base_cases_juveniles'=juvenile_ratio_baseline,'juvenile_cured'=juvenile_cured,'juvenile_prevalence'=juvenile_prevalence))
}