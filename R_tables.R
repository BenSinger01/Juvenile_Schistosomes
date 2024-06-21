rm(list=ls())
gc()

library(here)
setwd(here())

source("KK_sensitivity.R")
source('Schisto_Simulations_Functions.R')

require(parallel)

coverage <- 75/100
nonadherance <- 10/100
dyn_type <- "no_saturation"
epgpwp <- 4
eff_scale <- 100/100
nslots <- 4 # number of cores to use
N <- 50 # number of simulations

nocap <- FALSE
nodelay <- FALSE
redraw <- FALSE
community_size <- 500

seeds <- 1:N
nseeds <- length(seeds)

# period of simulation. Here, five years + time for follow-up.
T=5*52+11
lambda_param <- 1

for (case in c("Mean", "Lower", "Upper")) {
  for (sett in c("Low", "Moderate", "High")) {
    for (static in c(0, 50, 100)) {
      for (strat in c("single", "double", "novel1", "novelC", "novelB")) {
        setting <- sett
        drug_strategy <- strat
        static_force <- static/100
        source("Run_Simulation_Table.R")
      }
    }
  }
}
