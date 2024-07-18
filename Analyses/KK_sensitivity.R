# This function returns the result of a Kato-Katz test for schistosomiasis
# stochastically, depending on the number of slides and the number of stools
# taken. This is prameterised using the medians of Barenbold 2018's posterior
# distributions.

EPG_dispersion_median <- 5.537135243621125
EPG_daily_sd_median <- 1.2738676151429251

KK_result <- function(EPG,nslides=1,ndays=1){
	if (EPG==0){
		return(0)
	} else {
		intensity <- (EPG/24)/exp((EPG_daily_sd_median**2)/2)
		result <- 0
		for (day in 1:ndays){
			day_intensity <- intensity*exp(rnorm(n=1,mean=0,sd=EPG_daily_sd_median))
			p_negative_slide <- (EPG_dispersion_median/(day_intensity+EPG_dispersion_median))**EPG_dispersion_median
			result <- result|any(rbinom(n=nslides,size=1,prob=1-p_negative_slide))
		}
		return(result)
	}
}