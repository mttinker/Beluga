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
# NOTE: need to load results file first
#
post_N1 = exp(mcmc[,which(vn=="logN1")])
pri_N1 = exp(rnorm(Nsims,8.5,1))
df_post_N = data.frame(Estimate = c(rep("Prior",length(pri_N1)),
                                     rep("Posterior",length(post_N1)) ),
                        Value = c(pri_N1,post_N1)) 
df_post_N$Estimate = factor(df_post_N$Estimate)
ggplot(df_post_N,aes(x=Value,group=Estimate,fill=Estimate)) +
  geom_density(adjust=1.5,alpha=0.3) +
  scale_x_continuous(limits = c(0,50000)) +
  ggtitle("Prior vs Posterior, Initial population size") +
  theme_classic() 
#
post_gmmaA = mcmc[,which(vn=="gamma_A")]
post_gmmaY = mcmc[,which(vn=="gamma_Y")]
post_S_A = mcmc[,which(vn=="S_A")]
tmp = rnorm(Nsims*2.5,0,2); tmp = sample(tmp[tmp>-2 & tmp<3],Nsims)
pri_S_A = exp(-exp(tmp - 2.25))
df_post_SA = data.frame(Estimate = c(rep("Prior",length(pri_S_A)),
                                     rep("Posterior",length(post_S_A)) ),
                        Value = c(pri_S_A,post_S_A)) 
df_post_SA$Estimate = factor(df_post_SA$Estimate)
plt1 = ggplot(df_post_SA,aes(x=Value,group=Estimate,fill=Estimate)) +
  geom_density(adjust=1.5,alpha=0.3) +
  scale_x_continuous(limits = c(0,1)) +
  ggtitle("Prior vs Posterior, Adult Survival") +
  theme_classic() + theme(legend.position = "none")
#
post_S_Y = mcmc[,which(vn=="S_Y")]
tmp = rnorm(Nsims*2.5,0,2); tmp = sample(tmp[tmp>0 & tmp<3],Nsims)
pri_S_Y = exp(-exp(-2.25 + tmp + post_gmmaA ))
df_post_SY = data.frame(Estimate = c(rep("Prior",length(pri_S_Y)),
                                     rep("Posterior",length(post_S_Y)) ),
                        Value = c(pri_S_Y,post_S_Y)) 
df_post_SY$Estimate = factor(df_post_SY$Estimate)
plt2 = ggplot(df_post_SY,aes(x=Value,group=Estimate,fill=Estimate)) +
  geom_density(adjust=1.5,alpha=0.3) +
  scale_x_continuous(limits = c(0,1)) +
  ggtitle("Prior vs Posterior, Yearling Survival") +
  theme_classic() + theme(legend.position = "none")
#
post_S_NN = mcmc[,which(vn=="S_NN_mn")]
tmp = rnorm(Nsims*2.5,0,1); tmp = sample(tmp[tmp>-3 & tmp<3],Nsims)
tmp2 = rnorm(Nsims*2.5,0,1); tmp2 = sample(tmp2[tmp2>0],Nsims)
pri_S_NN = (exp(-exp(-.85 + tmp + 0.1*tmp2)))
df_post_SN = data.frame(Estimate = c(rep("Prior",length(pri_S_NN)),
                                     rep("Posterior",length(post_S_NN)) ),
                        Value = c(pri_S_NN,post_S_NN)) 
df_post_SN$Estimate = factor(df_post_SN$Estimate)
plt3 = ggplot(df_post_SN,aes(x=Value,group=Estimate,fill=Estimate)) +
  geom_density(adjust=1.5,alpha=0.3) +
  scale_x_continuous(limits = c(0,1)) +
  ggtitle("Prior vs Posterior, Neonatal Survival (mean)") +
  theme_classic()
#
post_Pr = mcmc[,which(vn=="Pr_mn")]
tmp = rnorm(Nsims,0,1.5)
pri_Pr = inv.logit(tmp)
df_post_Pr = data.frame(Estimate = c(rep("Prior",length(pri_Pr)),
                                     rep("Posterior",length(post_Pr)) ),
                        Value = c(pri_Pr,post_Pr)) 
df_post_Pr$Estimate = factor(df_post_Pr$Estimate)
plt4 = ggplot(df_post_Pr,aes(x=Value,group=Estimate,fill=Estimate)) +
  geom_density(adjust=1.5,alpha=0.3) +
  scale_x_continuous(limits = c(0,1)) +
  ggtitle("Prior vs Posterior, Pregancy Rate (mean)") +
  theme_classic()
#
rm(tmp,tmp2)
grid.arrange(grobs=list(plt1,plt3,plt2,plt4),nrow=2)



