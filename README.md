# Population level impact of mass drug administration against schistosomiasis with novel anthelmintic drugs targeting juvenile schistosomes: a modelling study

## Summary
This repository contains analytic code to simulate a human population exposed to _Schistosoma mansoni_ parasites, including explicit representation of the juvenile parasites. It includes an individual-based mechanistic model of _S. mansoni_ infection and simulated mass drug administration control programs in three different epidemiologic environments. It can simulate control programs with praziquantel (single-dose and two-dose regimen) and hypothetical novel drugs with varying assumptions on efficacy against adult and juvenile parasites. Outputs include figures illustrating infection prevalence over time.

Fixed model parameters were calibrated through a grid search (see Methods section of study). The data used for this analysis is publicly available  here: [SCORE data](https://clinepidb.org/ce/app/workspace/analyses/DS_d6a1141fbf/new/details#Contacts).

The study has been accepted for publication in the Lancet Microbe, with the final article in preparation.

This code was written using R version 4.4.1 and Python version 3.9.7.

## Structure
```
Juvenile_Schistosomes/  
└──Analysis/ - All code for simulations and plotting.  
└──Data/  
│  └─baselines/ - Simulated baseline states stored here.  
│  └─MDA_sims/ - Simulations of mass drug administration stored here.  
└──Figures/ - Figures output here.
```  

## Running the Main Analysis
To run the main analysis and reproduce figures 1–3 in the manuscript, simply call bash Run_Analysis.sh from the command line in the Juvenile_Schistosomes folder. To run sensitivity analyses, edit parameters in Run_Analysis.sh.

## Contact
Please direct any questions to the study authors, either by submitting an issue or emailing:

Benjamin J Singer, Stanford University, benjaminjohnsinger@gmail.com

Nathan Lo, Stanford University, nathan.lo@stanford.edu
