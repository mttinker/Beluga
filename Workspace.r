# Work space
# 
# Beluga matrix with 3 adult female compartments (available, pregnant, with calf)
M0 = matrix(0, nrow = 12,ncol = 12)
# Note: stages 1:8 = ages 0-7 (both sexes), 9 = 8+ males, 
#   10 = 8+ females "available", 11 = 8+ females preg, 12 = 8+ females w. calves 
Sn = .7   # newborn survival (birth to time of survey, Sept 1)
Sc = .8   # calf survival (from 0 to 1 yr old)
Sy = .85   # yearling survival (1 to 2 year old)
Sa = .95  # adult survival (annual survival for 2yr old and older)
Pr = .9   # pregnancy rate (for "available" females)

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
# annual demographic transitions accomplished with single matrix mult
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
# - use beta distribution for proportion immatures in survey
# - use multinomial (maybe Dirichlet-multinomial) for age comp of carcasses
#   (relative prob of observing newborns reduced, maybe start as fixed ratio)
# - log hazards (gamma) for adults/yearling/age0, incremental approach:
#    - haz_A = exp(-5 + gamma_A + epsA)
#    - haz_Y = exp(-5 + gamma_A + gamma_Y + epsA)
#    - haz_0 = exp(-5 + gamma_A + gamma_Y + gamma_0 + eps0)
# - age0 log hazards include stochastic hierarchical random effect (eps0))
# - for yearling/adults, also stochastic? use scaled value of eps0, or independent epsA?
# - preg rate stochastic also? may not be identifiable from newborn mortality...
# - spit age0 hazards between newborn and calf survival, i.e. Sn = exp(-1*haz_0*fraction_nb)  
# - add density dependence term to haz_0? or other covariates to pregnancy and survival? 
#  
# Example of user function in STAN with vector (should also work with matrix):
#  (the example function standardizes a vector variable to have a mean of 0 
#    and standard deviation of 1, or just center it if scale=0.)
# functions {
#   vector stdized(int N, vector x, int scale) {
#     vector[N] x_sc;
#     
#     x_sc = scale ? x-mean(x) : (x-mean(x))/sd(x);
#     
#     return x_sc;
#   }
# }


