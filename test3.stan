# This is a test of STAN in emacs using stan-mode
#
# looks like the color coding in stan-mode is working
#
# Steven Skates 
# 23 Jan 2014
#
#

data {
	int<lower=1> nMarkers; //number of markers
	int<lower=1> nPairs; //number of pairs
	real	cases[nMarkers,nPairs]; //1 vector of length nPairs per marker
	real	paired_controls[nMarkers,nPairs];
	real	unpaired_controls[nMarkers,nPairs];
}
parameters {
	//parameters
		real mu[nMarkers,2]; // mu[marker#,group]
		real<lower=0> sigma[nMarkers,2]; // sigma[marker#,group]
		real<lower=0,upper=1> rho[nMarkers]; //correlation coefficient
		real<lower=0,upper=1> theta[nMarkers]; //vector proportions of cases like controls
	//hyperparameters
		real mumu[2]; //LOC location hyperparameters
		real<lower=0> musig[2]; //LOC scale hyperparameter
		real sigmu[2]; //log(SCALE) location hyperparameters
		real<lower=0> sigsig[2]; //log(SCALE) scale hyperparameter
}
model {
	real ps;
	for (z in 1:2) {
	//mu hyperpriors
		mumu[z] ~ student_t(5,4,4);
		log(musig[z]) ~ student_t(5,0.5,0.5);
		increment_log_prob(-log(musig[z]));
	//sigma hyperpriors
		sigmu[z] ~ student_t(5,-1.5,0.5);
		log(sigsig[z]) ~ student_t(5,-1,0.75);
		increment_log_prob(-log(sigsig[z]));
	}
	for (i in 1:nMarkers) {
		//priors
			theta[i] ~ uniform(0,1);
			rho[i] ~ uniform(0,1);
		//hyperparameters
			for (z in 1:2) {
				mu[i,z] ~ student_t(5,mumu[z],musig[z]);
				log(sigma[i,z]) ~ student_t(5,sigmu[z],sigsig[z]);
				increment_log_prob(-log(sigma[i,z]));
			}
		//parameters
			for (j in 1:nPairs) {
				unpaired_controls[i,j] ~ student_t(5,mu[i,1],sigma[i,1]);
				paired_controls[i,j] ~ student_t(5,mu[i,1],sigma[i,1]);			
				ps <- log_sum_exp(
					log(theta[i]) + student_t_log(cases[i,j]
										,6
										,mu[i,1] + rho[i]*(paired_controls[i,j] - mu[i,1])
										,sqrt((1-rho[i]*rho[i])*(5+square((paired_controls[i,j]- mu[i,1])/sigma[i,1]))/6)*sigma[i,1])
					,log(1-theta[i]) + student_t_log(cases[i,j]
								,6
								,mu[i,2] + (sigma[i,2]/sigma[i,1])*rho[i]*(paired_controls[i,j] - mu[i,1])
								,sqrt((1-rho[i]*rho[i])*(5+square((paired_controls[i,j]- mu[i,1])/sigma[i,1]))/6)*sigma[i,2])
					); // cases are paired and mixture
				increment_log_prob(ps);
			}
	}
}
generated quantities {
	real diff[nMarkers]; //mu2-mu1
	real cut_point[nMarkers];
	real wp_cut_point[nMarkers]; //within-person cutpoint
	real<lower=0,upper=1> sens98[nMarkers]; //@spec <- 0.98
	real<lower=0,upper=1> wp_sens98[nMarkers]; //@spec <- 0.98, within-person sens98

	for (k in 1:nMarkers) {
		diff[k] <- mu[k,2] - mu[k,1];
		cut_point[k] <- mu[k,1] + 2.756509 * sigma[k,1];
		sens98[k] <- 1.0 - (
				theta[k] * student_t_cdf(cut_point[k],5,mu[k,1],sigma[k,1])
				+(1-theta[k]) * student_t_cdf(cut_point[k],5,mu[k,2],sigma[k,2])
				); //model based proportion of cases above 98th quantile of cntl mean
		wp_cut_point[k] <- 2.756509 * (sigma[k,1]*sqrt(2*(1-rho[k]))); // (mu[k,1]-mu[k,1])+2.756...
		wp_sens98[k] <- 1.0 - (
				theta[k] * student_t_cdf(wp_cut_point[k],5,0,sigma[k,1]*sqrt(2*(1-rho[k])))
				+(1-theta[k]) * student_t_cdf(wp_cut_point[k],5,diff[k],sqrt(square(sigma[k,2])+square(sigma[k,1])+2*sigma[k,2]*sigma[k,1]*rho[k]))
				); //within-person sens98
	}
}
