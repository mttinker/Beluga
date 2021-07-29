# Beluga model
# SOund quality detection probability 
require(readxl)
require(gtools)
require(parallel)
require(ggplot2)
require(dplyr)
require(fitdistrplus)
require(stats)
require(bayesplot)
require(gridExtra)
require(htmlTable)
require(loo)
require(cmdstanr)
require(posterior)
rstan::rstan_options(javascript=FALSE)
#
# Load data -----------------------------------------------


# Process data -------------------------------------------





# Fit model ---------------------------------------------
nburnin = 500                    # number warm-up (burn-in) samples
nsamples = 5000                 # desired total number samples
fitmodel = c("Beluga_fit.stan")    # name of file with stan code
stan.data <- list(N1=Nobs1,N2=Nobs2,Strikes=Strikes,Trials=Trials,
                  Detect=Detect,Nsite1=Nsite1,Nsite2=Nsite2,K=K,
                  S1=S1,S2=S2,X1=X1,X2=X2,ix1=ix1,ix2=ix2,Mic=Mic) # 
parms <- c("ppp","Tstat","Tstat_new","tau","sig1","sig2","kappa1","kappa2","Upsilon",
           "alpha","ynew","log_lik") # 
cores = detectCores()
ncore = min(20,cores-4)
Niter = round(nsamples/ncore)
mod <- cmdstan_model(fitmodel)   # compiles model (if necessary)
suppressMessages(                # Suppress messages/warnings (if desired)
  suppressWarnings (
    fit <- mod$sample(
      data = stan.data,
      seed = 1234,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      # init = init.fun,
      iter_warmup = nburnin,
      iter_sampling = Niter,
      max_treedepth = 12,
      adapt_delta = 0.8
    )
  )
)
source("cmdstan_sumstats.r")

# Diagnostic Plots -------------------------------------
mcmc_trace(fit$draws("tau"))

rr = sample(Nsims,pmin(1000,Nsims))
y = c(Strikes[ix1], Detect[ix2]); prop_zero <- function(x) mean(x == 0)
y_rep = as.matrix(mcmc[,startsWith(vn, "ynew[")])
ppc_stat(y, y_rep[rr,], stat = "prop_zero") + 
  ggtitle ("Posterior predictive distribution, observed vs out-of-sample predictions") +
  labs(x = "Proportion of records with 0 strikes",y = "Frequency") + theme_classic()

Bayes_P = sumstats[which(vns=="ppp"),1]
if(Bayes_P>0.75){Bayes_P=1-Bayes_P}
xx = as.matrix(mcmc[,vns=="Tstat"][rr])
yy = as.matrix(mcmc[,vns=="Tstat_new"][rr])
df_ppc = data.frame(x = log(xx), y = log(yy))
ggplot(df_ppc,aes(x=x,y=y)) +
  geom_point(color="blue") + 
  labs(x="Discrepancy measure for actual data set",
       y = "Discrepancy measure for new data",
       title= "Posterior predictive check, sum of squared Pearson residuals", 
       subtitle = paste0("Bayesian-P = ",Bayes_P)) +
  geom_abline(slope=1,intercept=0) +
  theme_classic()

mcmc_areas(fit$draws(variables = c("tau","sig1","sig2")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle("Parameter posterior distributions, variance params") +
  labs(x="Parameter value",y="Posterior density of predictor variables") +
  theme_classic()
