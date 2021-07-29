# Work space
# 
# Beluga matrix with 3 adult female compartments (available, pregnant, with calf)
M0 = matrix(0, nrow = 12,ncol = 12)
# Note: stages 1:8 = ages 0-7 (both sexes), 9 = 8+ males, 
#   10 = 8+ females "available", 11 = 8+ females preg, 12 = 8+ females w. calves 
Sn = .83   # early newborn survival (birth to time of survey, Sept 1)
Sc = .8    # calf survival (from 0.5 yr old to 1.5 yr old at next survey)
Sy = .87   # yearling survival (1.5 to 2.5 year old)
Sa = .95   # adult survival (annual survival for 2yr old and older)
Pr = .9    # pregnancy rate (for "available" females)

Make_matrix <- function(M0,Sn,Sc,Sy,Sa,Pr) {
  M = M0
  # or STAN version: matrix[12,12] M; M = rep_matrix(0,12,12) ;
  M[1,11] = Sa * Sn
  M[2,1] = Sc 
  M[3,2] = Sy
  M[4:8,3:7] = diag(1,5,5)*rep(Sa,5) 
  # STAN version: M[4:8,3:7] = add_diag(M[4:8,3:7], rep_vector(Sa,5)) 
  #  or repeat 5 steps, e.g.: M[4,3] = Sa
  M[9,8] = Sa/2     
  M[9,9] = Sa
  M[10,8] = Sa/2
  M[10,10] = Sa * (1-Pr)
  M[10,11] = Sa * (1-Sn)
  M[10,12] = Sa 
  M[11,10] = Sa * Pr
  M[12,11] = Sa * Sn
  return(M)
}
M = Make_matrix(M0,Sn,Sc,Sy,Sa,Pr)
# get stable state distribution (ssd)
eigs = eigen(M)
lambda = eigs$values[1]
W=abs(eigen(M)$vectors[,1]) # W=matrix of right eigenvectors 
ssd = W/sum(W)              # stable stage distribution
n = 1000 * ssd
# annual demographic transitions accomplished with single matrix multiplication
n = M %*% n
# Use "Death Matrix" to create expected death assemblage for 9 age classes (0yr, 1:8+)
DM = diag(1,9,9)*c(1-Sc,1-Sy,rep((1-Sa),7))
# STAN version: DM and nt and nd declared as local variables, and 
#  DM = add_diag(rep_matrix(0,9,9),[1-Sc,1-Sy,rep_vector((1-Sa),7)]) ;
nt = c(n[1:8],sum(n[9:12]))  # Age-based vector (9 classes instead of 12)
# Matrix multiplication to get # dead by age class 
nd = DM %*% nt
# adjust age0 to add newborns who are born and die before next survey
nd[1] = nd[1] + n[11] * Sa * (1-Sn)

# ADDITIONAL NOTES:
# - use beta distribution for observed proportion immature in survey
# - use multinomial (maybe Dirichlet-multinomial?) for age comp of carcasses
#   (relative prob of observing newborns reduced, maybe start as fixed ratio)
# - log hazards (gamma) for adults/yearling/age0, use incremental approach
#   where gmm0 ~ normal(0,1), gamma_A = gmm0 - 3, and each additional gamma term
#   is a log hazard ratio required to be > 0 (half-normal prior):
#    - haz_A = exp(gamma_A + eps_SA)
#    - haz_Y = exp(gamma_A + gamma_Y + eps_S1)
#    - haz_C = exp(gamma_A + gamma_Y + gamma_C + eps_S0)
#     S_A = exp(-haz_A)   : survival for immatures/adults
#     S_Y = exp(-haz_Y)   : survival for yearlings
#     S_C = exp(-haz_C)   :survival for newborns, ~ 6mo to ~ 1.5 yrs (yearling) 
#     S_N = exp(-haz_C/2) ; #(fractional mortality for newborns March to Sept) 
# - age0 log hazards include stochastic hierarchical random effect (eps_S0),
#    where eps0 ~ normal(0,1) and eps_S0 = eps0 * sig_S0
# - stochasticity for yearlings may be lower and independent from age0, 
#        so eps1 ~ normal(0,1), eps_S1 = eps1 * sig_S1 
#       and eps_SA = eps1 * sig_SA (separate variance term)
# - preg rate Pr = inv.logit(Beta0 - Beta1*N_t/1000 + eps_P), 
#    where epsP ~ normal(0,1) and eps_P = epsP * sig_P
# - add density-depend. term to preg AND haz_Y? other covariates for pregnancy/survival? 
#  



