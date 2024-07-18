## Run Analyses
## Benjamin J Singer, July 2024.
## Shell script to run main analyses and generate figures 1–3 from Singer et al. 2024

# Number of cores available on machine for parallelization, set as high as possible so code runs fast
nslots=4

# Number of repeats to run of each simulation, 250 in original paper
N=4

# Vary these baseline parameters for sensitivity analyses
epgpwp=4
juvenile_duration=6
seasonality=0
t_start=0

# Burn in baselines
# for seed in $(seq 1 $N);
# 	do for setting in Low Moderate High;
# 		do Rscript Analyses/Burn_In_Baseline.R $setting $epgpwp $juvenile_durattion $seasonality $t_start $seed;
# 	done;
# done

# Vary these treatment parameters for sensitivity analyses
dyn_type=no_saturation
coverage=75
nonadherance=10
eff_scale=100
wash=0

# Simulate MDA
# for case in Mean Upper Lower
# 	do for setting in Low High Moderate
# 		do for static_force in 100 0 50
# 			do for strat in Single Double
# 				do Rscript Analyses/Simulate_MDA.R $epgpwp $juvenile_duration $setting $seasonality $t_start $static_force $dyn_type $case $strat $coverage $nonadherance $eff_scale $wash $nslots $N
# 			done
# 		done
# 	done
# done

# for setting in Low High Moderate
# 	do for static_force in 100 0 50
# 		do for strat in NovelA NovelB NovelC
# 			do Rscript Analyses/Simulate_MDA.R $epgpwp $juvenile_duration $setting $seasonality $t_start $static_force $dyn_type Mean $strat $coverage $nonadherance $eff_scale $wash $nslots $N
# 		done
# 	done
# done

# Recreate Figure 1
# python3 Analyses/Plot_MDA.py $epgpwp $juvenile_duration $seasonality $t_start $dyn_type $coverage $nonadherance $eff_scale $wash 4

# Simulate MDA with varying efficacies
for static_force in 100 0
	do for effa in 90 95 100
		do for effj in 0 50 90 95 100
			do for setting in Low High
				do Rscript Analyses/Simulate_MDA_Grids.R $epgpwp $juvenile_duration $setting $seasonality $t_start $static_force $dyn_type $coverage $nonadherance $effa $effj $nslots $N
			done
		done
	done
done

# Vary these parameters to generate alternate figures
n_years=4
delay=5

# Recreate Figures 2 and 3
python3 Analyses/Plot_Grid.py $epgpwp $juvenile_duration $seasonality $t_start 100 $dyn_type $coverage $nonadherance $n_years $delay
python3 Analyses/Plot_Grid.py $epgpwp $juvenile_duration $seasonality $t_start 0 $dyn_type $coverage $nonadherance $n_years $delay