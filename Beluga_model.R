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
initvec = read_excel("data/Beluga_data.xlsx",sheet = 1)
Harvest = read_excel("data/Beluga_data.xlsx",sheet = 2)
Hv_ages = read_excel("data/Beluga_data.xlsx",sheet = 3)
df_Surv = read_excel("data/Beluga_data.xlsx",sheet = 4)
df_Strd = read_excel("data/Beluga_data.xlsx",sheet = 5)

# Process data -------------------------------------------
Year1 = min(Harvest$Year) ; 
YearT = max(df_Strd$Year) ;
Years = seq(Year1,YearT)  ;
YearsN = Years - Year1 + 1
Nyrs = max(YearsN) ;
NyrsH = max(Harvest$Year) - Year1 + 1
NStg = 12; NAge = 9
YrSv = df_Surv$Year - Year1 + 1
Nsrv = length(YrSv)
ObsS = df_Surv$`Pop size estimate`
VarS = (df_Surv$SE)^2
invSc = ObsS / VarS
PJ = df_Surv$Prop_juv
YrSt = df_Strd$Year - Year1 + 1
Nstr = length(YrSt)
StrNB = df_Strd$NB_shift
StrOA = df_Strd$Older_shift
Harv = pmax(0.01,Harvest$Removals)
#
# Determine appropriate precision param for Beta dist, proportion Juveniles
# (assuming binomial distrib. of raw counts of Juvs vs all animals)
# N_ob = round(df_Surv$`Pop size estimate`/(2.021*2.09))
# N_jv = round(df_Surv$Prop_juv * N_ob)
# N = round(mean(N_ob)); J = round(mean(N_jv)); P = J/N
# P_rnd = (rbinom(100000,N,P))/N
# coefs = coef(fitdist(P_rnd,"beta")); upsilon = round(coefs[1]/P)
# rm(P_rnd,coefs,N,P,J,N_ob,N_jv)
upsilon = 245
#
# *** STILL NEED TO ADD AGE DIST OF STRANDED ANIMALS ***
#
# Initialize Beluga matrix with 12 stages, 9 age classes: 
#  stages 1:8 = ages 0-7 (both sexes), stage 9 = 8+ males, 
#  stages 10 = 8+ fems "available", 11 = 8+ fems preg, 12 = 8+ fems w. calves
#    (Priors for vital rates based on estimates from Mosnier 2014)
Sn = .72     # newborn calf survival (from 0.5 yr old to 1.5 yr old at next survey)
Snn = sqrt(Sn) # neonate survival (birth to time of survey, Sept 1)
Sy = .9   # yearling survival (1.5 to 2.5 year old)
Sa = .95   # adult survival (annual survival for 2yr old and older)
Pr = .82    # pregnancy rate (for "available" females)
source("Make_matrix.R")
M = Make_matrix(Snn,Sn,Sy,Sa,Pr)
lambda = eigen(M); lambda = lambda$values[1]
ssd=abs(eigen(M)$vectors[,1]); ssd = ssd/sum(ssd)  # stable stage distribution 
#
stan.data <- list(Nyrs=Nyrs,NyrsH=NyrsH,NStg=NStg,NAge=NAge,
                  Nsrv=Nsrv,Nstr=Nstr,YrSv=YrSv,YrSt=YrSt,
                  ObsS=ObsS,invSc=invSc,PJ=PJ,StrNB=StrNB,StrOA=StrOA,
                  ssd=ssd,Harv=Harv,upsilon=upsilon) # 
parms <- c("ppp","Tstat","Tstat_new","sig_N","sig_P","sig_H","alpha","phi",
           "Pr_mn","S_A","S_Y","S_N_mn","gamma_A","gamma_Y","gamma_N","gamma_H_mn",
           "PD_NB","PD_OA","N","Pr","S_N","ynew") # 
#
init_fun <- function() {list(sig_N=runif(1, 1.45, 1.55),
                             sig_P=runif(1, 1, 2),
                             sig_H=runif(1, 2.8, 3),
                             logNinit=runif(1, .55, .65),
                             gammA=runif(1, .275, .325),
                             gamma_Y=runif(1, .35, .45),
                             gamma_N=runif(1, .45, .55),
                             alpha=runif(1, -1, 1),
                             phi=runif(1, .65, .75),
                             PD_NB=runif(1, .05, .07),
                             PD_OA=runif(1, .2, .21),
                             gamma_H_mn=runif(1, .15, .25)
)}   
#
# Fit model ---------------------------------------------
nburnin = 500                    # number warm-up (burn-in) samples
nsamples = 5000                 # desired total number samples
fitmodel = c("Beluga_fit.stan")    # name of file with stan code
mod <- cmdstan_model(fitmodel)   # compiles model (if necessary)
cores = detectCores()
ncore = min(20,cores-4)
Niter = round(nsamples/ncore)
suppressMessages(                # Suppress messages/warnings (if desired)
  suppressWarnings (
    fit <- mod$sample(
      data = stan.data,
      seed = 1234,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      # init = init_fun,
      iter_warmup = nburnin,
      iter_sampling = Niter,
      max_treedepth = 12,
      adapt_delta = 0.8
    )
  )
)
# tmp = fit$output(); tmp[[1]][40:60]
source("cmdstan_sumstats.r")
#
# Diagnostic Plots -------------------------------------

mcmc_trace(fit$draws("sig_N"))
mcmc_trace(fit$draws("gamma_A"))

set.seed(123)
rr = sample(Nsims,1000)
y = ObsS ;
y_rep <- (as.matrix(mcmc[,which(startsWith(vn,"ynew["))]))
ppc_dens_overlay(y, y_rep[rr[1:100],]) + 
  ggtitle ("Posterior predictive distribution, observed vs out-of-sample predictions") +
  labs(x = "Survey Counts",y = "Frequency") + theme_classic()

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

mcmc_areas(fit$draws(variables = c("sig_P","sig_N","sig_H")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle("Parameter posterior distributions, variance params") +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()

mcmc_areas(fit$draws(variables = c("gamma_A","gamma_Y","gamma_N")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle("Parameter posterior distributions, log hazard params") +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()

mcmc_areas(fit$draws(variables = c("S_A","S_Y","S_N_mn","Pr_mn")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle("Parameter posterior distributions, vital rates") +
  scale_y_discrete(labels=c("Survival, adult","Survival, yearling",
                            "Survival, newborn","Pregnancy rate")) +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()

mcmc_areas(fit$draws(variables = c("PD_NB","PD_OA")),
           area_method="equal height",
           prob = 0.8) + 
  scale_y_discrete(labels=c("Prob detect Newborn","Prob detect Older")) +
  ggtitle("Parameter posterior distributions, carcass detection probabilities") +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()

Npred = sumstats[which(startsWith(vns,"N[")),1]
Npred_lo = sumstats[which(startsWith(vns,"N[")),4]
Npred_hi = sumstats[which(startsWith(vns,"N[")),8]
df_Nplt = data.frame(Year=Years,Npred=Npred,
                       Npred_lo=Npred_lo,Npred_hi=Npred_hi)
df_Nplt$Survey_est = rep(NA,Nyrs); df_Nplt$Survey_est_SE = rep(NA,Nyrs)
df_Nplt$Survey_est[YrSv] = ObsS; df_Nplt$Survey_est_SE[YrSv] = sqrt(VarS)

ggplot(df_Nplt[which(Years>1925),],aes(x=Year,y=Npred)) +
  geom_ribbon(aes(ymin=Npred_lo,ymax=Npred_hi),alpha=0.3) +
  geom_line() +
  geom_point(aes(y=Survey_est)) +
  geom_errorbar(aes(ymin=Survey_est-Survey_est_SE,ymax=Survey_est+Survey_est_SE)) +
  labs(x="Year",y="Estimated Abundance") +
  ggtitle("Beluga population abundance, model projections (1925-2012)") +
  theme_classic()

Spred = sumstats[which(startsWith(vns,"S_N[")),1]
Spred_lo = sumstats[which(startsWith(vns,"S_N[")),5]
Spred_hi = sumstats[which(startsWith(vns,"S_N[")),7]
df_Splt = data.frame(Year=Years,Spred=Spred,
                     Spred_lo=Spred_lo,Spred_hi=Spred_hi)

ggplot(df_Splt[which(Years>1989),],aes(x=Year,y=Spred)) +
  geom_ribbon(aes(ymin=Spred_lo,ymax=Spred_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Estimated surviva rate, newborns") +
  ggtitle("Beluga newborn survival rate, model projections (1990-2012)") +
  theme_classic()

Ppred = sumstats[which(startsWith(vns,"Pr[")),1]
Ppred_lo = sumstats[which(startsWith(vns,"Pr[")),5]
Ppred_hi = sumstats[which(startsWith(vns,"Pr[")),7]
df_Pplt = data.frame(Year=Years,Ppred=Ppred,
                     Ppred_lo=Ppred_lo,Ppred_hi=Ppred_hi)

ggplot(df_Pplt[which(Years>1989),],aes(x=Year,y=Ppred)) +
  geom_ribbon(aes(ymin=Ppred_lo,ymax=Ppred_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Estimated pregancy rate") +
  ggtitle("Beluga adult pregnancy rate, model projections (1990-2012)") +
  theme_classic()

fit$save_object(file = paste0("results/beluga_rslt_",Sys.Date(),"_fit.RDS"))
# fit = readRDS(paste0(filename))

rm(fit,mod)
save.image(file=paste0("results/beluga_rslt_",Sys.Date(),".rdata"))

fit = readRDS(paste0("results/beluga_rslt_",Sys.Date(),"_fit.RDS"))

