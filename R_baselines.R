rm(list=ls())
gc()

library(here)
setwd(here())

source("KK_sensitivity.R")
source('Schisto_Simulations_Functions.R')

num_baseline <- 50
epgpwp <- 4
juvenile_duration <- 6

nocap <- FALSE
community_size <- 500
lambda_param <- 1

for (seed in 1:num_baseline) {
  for (sett in c("Low", "Moderate", "High")) {
    setting <- sett
    set.seed(as.numeric(seed))
    source("baseline_states_nbinom.R")
  }
}
