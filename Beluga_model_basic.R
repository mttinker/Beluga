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
# df_Strd = read_excel("data/Strandings1983-2020_20210802.xlsx")
# Process data -------------------------------------------
Year1 = min(Harvest$Year) ; 
#YearT = max(df_Strd$YYYY) ;
YearT = max(df_Strd$Year) ;
Years = seq(Year1,YearT)  ;
YearsN = Years - Year1 + 1
Nyrs = max(YearsN) ;
NyrsH = max(Harvest$Year) - Year1 + 1
NStg = 12; NAge = 9 # (newborns + 8 year-classes)
YrSv = df_Surv$Year - Year1 + 1
Nsrv = length(YrSv)
ObsS = df_Surv$`Pop size estimate`
VarS = (df_Surv$SE)^2
invSc = ObsS / VarS
PJ = df_Surv$Prop_juv
Harv = pmax(0.01,Harvest$Removals)
# Process stranding data: group counts by survey year (Oct - Sept)
#  and make matrix of age at death vectors for 1+ ages
# source("Stranding_Age_process.R")
#
StrNB = df_Strd$Newborns
StrOA = df_Strd$Older
YrSt = df_Strd$Year - Year1 + 1
Nstr = length(YrSt)
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
# Initialize Beluga matrix with 12 stages, 9 age classes: 
#  stages 1:8 = ages 0-7 (both sexes), stage 9 = 8+ males, 
#  stages 10 = 8+ fems "available", 11 = 8+ fems preg, 12 = 8+ fems w. calves
#    (Priors for vital rates based on estimates from Mosnier 2014)
Snn = .9     # newborn calf survival(3 months)
# Snn = SnnA^(1/4) # neonate survival (birth to time of survey, Sept 1)
Sy = .9   # yearling survival 
Sn = Snn^2 * sqrt(Sy)
Sa = .94   # adult survival (annual survival for 2yr old and older)
Pr = .8    # pregnancy rate (for "available" females)
source("Make_matrix.R")
M = Make_matrix(Snn,Sn,Sy,Sa,Pr)
lambda = eigen(M); lambda = Re(lambda$values[1])
ssd=abs(eigen(M)$vectors[,1]); ssd = ssd/sum(ssd)  # stable stage distribution 
nt = ssd*1000; Nt = sum(nt)
for(t in 1:100){
  nt = M %*% nt
  Nt = c(Nt,sum(nt))
}
plot(Nt)
#
stan.data <- list(Nyrs=Nyrs,NyrsH=NyrsH,NStg=NStg,NAge=NAge,
                  Nsrv=Nsrv,Nstr=Nstr,YrSv=YrSv,YrSt=YrSt,
                  ObsS=ObsS,invSc=invSc,PJuv=PJ,StrNB=StrNB,StrOA=StrOA,
                  ssd=ssd,Harv=Harv,upsilon=upsilon) # 
parms <- c("ppp","Tstat","Tstat_new","logN1","alpha",
           "gamma_N","phi","gamma_Y","gamma_A","gamma_H_mn",
           "sig_P","sig_N","sig_H","theta_NB","theta_OA",
           "Pr_mn","S_NN_mn","S_N_mn","S_Y","S_A","N","Pr",
           "S_NN","S_N","ppn_J","ppn_Av","ppn_Pr","ppn_Wc","ynew")  
#
# Fit model ---------------------------------------------
nburnin = 750                    # number warm-up (burn-in) samples
nsamples = 10000                 # desired total number samples
fitmodel = c("Beluga_fit_basic.stan")    # name of file with stan code
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
#
plt_trace1 = mcmc_trace(fit$draws("logN1"))
plt_trace2 = mcmc_trace(fit$draws("gamma_A"))
grid.arrange(plt_trace1,plt_trace2)
#
set.seed(123)
rr = sample(Nsims,1000)
y = ObsS ;
y_rep <- (as.matrix(mcmc[,which(startsWith(vn,"ynew["))]))
plt_ppc1 = ppc_dens_overlay(y, y_rep[rr[1:100],]) + 
  ggtitle ("Posterior predictive distribution, observed vs out-of-sample predictions") +
  labs(x = "Survey Counts",y = "Frequency") + theme_classic()
#
Bayes_P = sumstats[which(vns=="ppp"),1]
if(Bayes_P>0.75){Bayes_P=1-Bayes_P}
xx = as.matrix(mcmc[,vns=="Tstat"][rr])
yy = as.matrix(mcmc[,vns=="Tstat_new"][rr])
df_ppc = data.frame(x = log(xx), y = log(yy))
plt_ppc2 = ggplot(df_ppc,aes(x=x,y=y)) +
  geom_point(color="blue") + 
  labs(x="Discrepancy measure for actual data set",
       y = "Discrepancy measure for new data",
       title= "Posterior predictive check, sum of squared Pearson residuals", 
       subtitle = paste0("Bayesian-P = ",Bayes_P)) +
  geom_abline(slope=1,intercept=0) +
  theme_classic()
grid.arrange(plt_ppc1,plt_ppc2)
#
mcmc_areas(fit$draws(variables = c("sig_P","sig_N","sig_H")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle("Posterior distributions, variance parameters") +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()
#
mcmc_areas(fit$draws(variables = c("phi")),
           area_method="equal height",
           prob = 0.8) +
  ggtitle("Posterior distribution, density dependent log hazard ratio") +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()
#
mcmc_areas(fit$draws(variables = c("gamma_A","gamma_Y","gamma_N")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle("Posterior distributions, log hazard params") +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()
#
mcmc_areas(fit$draws(variables = c("S_A","S_Y","S_N_mn","S_NN_mn","Pr_mn")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle("Posterior distributions, vital rates") +
  scale_y_discrete(labels=c("Survival, adult","Survival, yearling",
                            "Survival, newborn", "Survival, neonatal (3 mo)",
                            "Pregnancy rate")) +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()
#
mcmc_areas(fit$draws(variables = c("theta_NB","theta_OA")),
           area_method="equal height",
           prob = 0.8) + 
  scale_y_discrete(labels=c("Prob detect Newborn","Prob detect Older")) +
  ggtitle("Posterior distributions, carcass detection probabilities") +
  labs(x="Parameter value",y="Posterior sample density") +
  theme_classic()
#
Npred = sumstats[which(startsWith(vns,"N[")),1]
Npred_lo = sumstats[which(startsWith(vns,"N[")),4]
Npred_hi = sumstats[which(startsWith(vns,"N[")),8]
df_Nplt = data.frame(Year=Years,Npred=Npred,
                       Npred_lo=Npred_lo,Npred_hi=Npred_hi)
df_Nplt$Survey_est = rep(NA,Nyrs); df_Nplt$Survey_est_SE = rep(NA,Nyrs)
df_Nplt$Survey_est[YrSv] = ObsS; df_Nplt$Survey_est_SE[YrSv] = sqrt(VarS)
#
ggplot(df_Nplt[which(Years>1925),],aes(x=Year,y=Npred)) +
  geom_ribbon(aes(ymin=Npred_lo,ymax=Npred_hi),alpha=0.3) +
  geom_line() +
  geom_point(aes(y=Survey_est)) +
  geom_errorbar(aes(ymin=Survey_est-Survey_est_SE,ymax=Survey_est+Survey_est_SE)) +
  labs(x="Year",y="Estimated Abundance") +
  ggtitle("Beluga population abundance, model projections (1925-2012)") +
  theme_classic()
#
Spred = sumstats[which(startsWith(vns,"S_NN[")),1]
Spred_lo = sumstats[which(startsWith(vns,"S_NN[")),5]
Spred_hi = sumstats[which(startsWith(vns,"S_NN[")),7]
df_Splt = data.frame(Year=Years,Spred=Spred,
                     Spred_lo=Spred_lo,Spred_hi=Spred_hi)
#
plt_SvNb = ggplot(df_Splt[which(Years>1982 & Years<2013),],aes(x=Year,y=Spred)) +
  geom_ribbon(aes(ymin=Spred_lo,ymax=Spred_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Estimated survival rate, neonates") +
  ggtitle("Beluga neonatal survival, model projections (1990-2012)") +
  theme_classic()
#
Ppred = sumstats[which(startsWith(vns,"Pr[")),1]
Ppred_lo = sumstats[which(startsWith(vns,"Pr[")),5]
Ppred_hi = sumstats[which(startsWith(vns,"Pr[")),7]
df_Pplt = data.frame(Year=Years,Ppred=Ppred,
                     Ppred_lo=Ppred_lo,Ppred_hi=Ppred_hi)
#
plt_PrgRt = ggplot(df_Pplt[which(Years>1982 & Years<2013),],aes(x=Year,y=Ppred)) +
  geom_ribbon(aes(ymin=Ppred_lo,ymax=Ppred_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Estimated pregancy rate") +
  ggtitle("Beluga adult pregnancy rate, model projections (1990-2012)") +
  theme_classic()
#
grid.arrange(plt_SvNb,plt_PrgRt)

P_Av_pred = sumstats[which(startsWith(vns,"ppn_Av[")),1]
P_Av_pred_lo = sumstats[which(startsWith(vns,"ppn_Av[")),5]
P_Av_pred_hi = sumstats[which(startsWith(vns,"ppn_Av[")),7]
P_Pr_pred = sumstats[which(startsWith(vns,"ppn_Pr[")),1]
P_Pr_pred_lo = sumstats[which(startsWith(vns,"ppn_Pr[")),5]
P_Pr_pred_hi = sumstats[which(startsWith(vns,"ppn_Pr[")),7]
P_Wc_pred = sumstats[which(startsWith(vns,"ppn_Wc[")),1]
P_Wc_pred_lo = sumstats[which(startsWith(vns,"ppn_Wc[")),5]
P_Wc_pred_hi = sumstats[which(startsWith(vns,"ppn_Wc[")),7]
df_Femstg_plt = data.frame(Year=Years,Status=rep("Available",Nyrs),Prop=P_Av_pred,
                                               Prop_lo=P_Av_pred_lo,Prop_hi=P_Av_pred_hi)
df_Femstg_plt = rbind(df_Femstg_plt,data.frame(Year=Years,Status=rep("Pregnant",Nyrs),
                                               Prop=P_Pr_pred,
                                               Prop_lo=P_Pr_pred_lo,Prop_hi=P_Pr_pred_hi) )
df_Femstg_plt = rbind(df_Femstg_plt,data.frame(Year=Years,Status=rep("With_calf",Nyrs),
                                               Prop=P_Wc_pred,
                                               Prop_lo=P_Wc_pred_lo,Prop_hi=P_Wc_pred_hi) )                      
#
plt_pAv = ggplot(df_Femstg_plt[which(df_Femstg_plt$Year>1982 & df_Femstg_plt$Status=="Available" ),],
                 aes(x=Year,y=Prop)) +  # group=Status,color=Status,fill=Status
  geom_ribbon(aes(ymin=Prop_lo,ymax=Prop_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Estimate") +
  ggtitle(paste0("Beluga female status (1983-",as.character(YearT),"): Proportion Available")) +
  theme_classic()
plt_pPr = ggplot(df_Femstg_plt[which(df_Femstg_plt$Year>1982 & df_Femstg_plt$Status=="Pregnant" ),],
                 aes(x=Year,y=Prop)) +  # group=Status,color=Status,fill=Status
  geom_ribbon(aes(ymin=Prop_lo,ymax=Prop_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Estimate") +
  ggtitle(paste0("Beluga female status (1983-",as.character(YearT),"): Proportion Pregnant")) +
  theme_classic()
plt_pWc = ggplot(df_Femstg_plt[which(df_Femstg_plt$Year>1982 & df_Femstg_plt$Status=="With_calf" ),],
                 aes(x=Year,y=Prop)) +  # group=Status,color=Status,fill=Status
  geom_ribbon(aes(ymin=Prop_lo,ymax=Prop_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Estimate") +
  ggtitle(paste0("Beluga female status (1983-",as.character(YearT),"): Proportion With Calf")) +
  theme_classic()
#
PJpred = sumstats[which(startsWith(vns,"ppn_J[")),1]
PJpred_lo = sumstats[which(startsWith(vns,"ppn_J[")),4]
PJpred_hi = sumstats[which(startsWith(vns,"ppn_J[")),8]
df_PJplt = data.frame(Year=Years,PJpred=PJpred,
                      PJpred_lo=PJpred_lo,PJpred_hi=PJpred_hi)
df_PJplt$PJ_obs = rep(NA,Nyrs); df_PJplt$PJ_obs_se = rep(NA,Nyrs)
df_PJplt$PJ_obs[YrSv] = PJ
N_ob = round(df_Surv$`Pop size estimate`/(2.021*2.09))
N_jv = PJ*N_ob ; se = sqrt(N_ob*PJ*(1-PJ)); CV = se/N_jv
df_PJplt$PJ_obs_se[YrSv] = CV*PJ
plt_pJv = ggplot(df_PJplt[which(Years>1982),],aes(x=Year,y=PJpred)) +
  geom_ribbon(aes(ymin=PJpred_lo,ymax=PJpred_hi),alpha=0.3) +
  geom_line() +
  geom_point(aes(y=PJ_obs)) +
  geom_errorbar(aes(ymin=PJ_obs-1.96*PJ_obs_se,ymax=PJ_obs+1.96*PJ_obs_se)) +
  labs(x="Year",y="Estimate") +
  ggtitle(paste0("Beluga survey age structure (1983-",
                 as.character(YearT),"): Proportion Juvenile")) +
  theme_classic()
#
grid.arrange(grobs=list(plt_pAv,plt_pPr,plt_pWc,plt_pJv),nrow=4)
#
source("Compare_priors.R")
# Save Results -----------------------------------------------------------
fit$save_object(file = paste0("results/beluga_base_rslt_",Sys.Date(),"_fit.RDS"))
# fit = readRDS(paste0(filename))

rm(fit,mod)
save.image(file=paste0("results/beluga_base_rslt_",Sys.Date(),".rdata"))

fit = readRDS(paste0("results/beluga_base_rslt_",Sys.Date(),"_fit.RDS"))

