// Beluga model
// This Stan code executes an age-structured population model for St Lawrence
//  Belugas, fit to surey data, mortaliy data and harvest data
//
functions {
  // user function to create projection matrix
  matrix makemat(int Ns, vector G, vector S) {
     matrix[Ns,Ns] M; M = rep_matrix(0,Ns,Ns) ;
     M[1:(Ns-1),1:(Ns-1)] = add_diag(M[1:(Ns-1),1:(Ns-1)], S[1:(Ns-1)] .* (1 -  G[1:(Ns-1)])) ;
     M[2:Ns,1:(Ns-1)] = add_diag(M[2:Ns,1:(Ns-1)], S[1:(Ns-1)] .* G[1:(Ns-1)]) ;
     return M;
  }
}
// Section 1. Data inputs for model
data {
  int<lower=1> Nyrs ;            // Total years of counts (if 2011-2020 then 10)
  int<lower=1> YrB ;             // Year before temporal break (if 2013 then 3)
  int<lower=1> Nsz ;             // Number size classes (assumed to be 10)
  vector[Nsz-1] Gr ;             // Growth transition probs for size classes
  simplex[Nsz] SZinit ;          // Initial size distribution
  vector[Nsz] lgtPD ;            // baseline logit detection fxn
  vector[Nyrs] Obs_Tot ;         // Observed total densities by year
  int Obs_Sz[Nyrs,Nsz] ;         // Observed "counts" by size class, by year
  real R0pri ;                   // prior for log mean Recruitment (R)
  real N0pri ;                   // prior for initial abundance
  vector[Nsz] Zeros ;            // vector of zeros (to make recruitment vector)
  real gamma0 ;                  // baseline log hazards
  real invscpri ;                // inverse scale param for gamma dist. densities
}
// Section 2. The parameters to be estimated 
parameters {
  real<lower=0> tau ;             // precision param for Dirichelet ditribution
  real<lower=0> invscale ;        // inverse scale param for Gamma ditribution
  real<lower=-10,upper=10> phiDmn;// mean change in detect prob post break, logit 
  real<lower=-5,upper=5> phiSmn ; // mean change in log hazards post break, ratio 
  real<lower=-2,upper=5> phiRmn ; // mean change in log recruit post break, ratio 
  real<lower=0> sigD1 ;           // variation in detect prob pre break, logit 
  real<lower=0> sigS1 ;           // variation in log hazards pre break, ratio 
  real<lower=0> sigR1 ;           // variation in log recruit pre break, ratio   
  real<lower=0> sigD   ;          // variation in detect prob post break, logit 
  real<lower=0> sigS   ;          // variation in log hazards post break, ratio 
  real<lower=0> sigR   ;          // variation in log recruit post break, ratio   
  vector[Nyrs] phiD ;             // random effect, detect prob by year
  vector[Nyrs] phiS ;             // random effect, baseline survival by year
  vector[Nyrs] phiR ;             // random effect, baseline recruitment by year
  simplex[Nsz] pie[Nyrs] ;        // predicted size dist. with overdispersion
} 
// Section 3. Additional transformed parameters, including key model dynamics
transformed parameters {
  vector[Nsz] R1 ;                // initial recruitment 
  vector[Nsz] S1;                 // initial survival rate
  vector[Nsz] P1 ;                // initial detection rate
  matrix[Nsz,Nsz] M1 ;            // initial matrix
  vector[Nsz] n0      ;           // initil pop vector
  simplex[Nsz] SzDst[Nyrs] ;      // predicted size dist, all years
  vector[Nyrs] N ;                // true density, all years
  vector[Nsz] n[Nyrs] ;           // true population vector, all years
  vector[Nyrs] D ;                // detectable density, all years
  vector[Nsz] d[Nyrs] ;           // detectable density by size, all years
  // Initial Prob detection vector
  P1 = inv_logit(lgtPD + phiD[1]) ;
  // Initial survival
  S1 = rep_vector(exp(-exp(gamma0 + phiS[1])),Nsz)  ;
  // Initial recruitment
  R1 = Zeros ;
  R1[1] = exp(R0pri + phiR[1]) ;
  // Projection matrix for first year
  M1 = makemat(Nsz,Gr,S1) ;
  // Initialize pop vector
  n0 = N0pri * SZinit ;
  // Demographic calcs for first year
  n[1] = M1 * (n0 + R1) ;
  N[1] = sum(n[1]) ;
  d[1] = P1 .* n[1] ;
  D[1] = sum(d[1]) ;
  SzDst[1] = (d[1] + .00000001) / (D[1] + Nsz * .00000001 ) ;
  // Loop through remaining years
  for (t in 2:Nyrs){
    vector[Nsz] R2 ; 
    vector[Nsz] S2 ;  
    vector[Nsz] P2 ;
    vector[Nsz] nt1 ;
    matrix[Nsz,Nsz] M2;
    R2 = Zeros ;
    P2 = inv_logit(lgtPD + phiD[t] ) ;
    S2 = rep_vector(exp(-exp(gamma0 + phiS[t] )),Nsz)  ;
    R2[1] = exp(R0pri + phiR[t]) ;
    M2 = makemat(Nsz,Gr,S2) ;
    n[t] = M2 * (n[t-1] + R2) ;
    N[t] = sum(n[t]) ;
    d[t] = P2 .* n[t] ;
    D[t] = sum(d[t]) ;
    SzDst[t] = (d[t] + .00000001) / (D[t] + Nsz * .00000001 ) ;
  }
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // -Uchin Counts, total by year 
  Obs_Tot ~ gamma(D * invscale, invscale) ; 
  // -Sizeclass distributions by year 
  //  NOTE: Use Dirichlet-multinomial distribution for overdispersed size freqs 
  for(t in 1:Nyrs){
    pie[t] ~ dirichlet(tau * SzDst[t]);
    Obs_Sz[t,] ~ multinomial(pie[t]) ;
  }
  // B) Prior distributions for model parameters:
  // Hierarchical random effects:  
  phiD[1:YrB] ~ normal(0, sigD1) ;
  phiS[1:YrB] ~ normal(0, sigS1) ;
  phiR[1:YrB] ~ normal(0, sigR1) ;
  phiD[(YrB+1):Nyrs] ~ normal(phiDmn, sigD) ;
  phiS[(YrB+1):Nyrs] ~ normal(phiSmn, sigS) ;
  phiR[(YrB+1):Nyrs] ~ normal(phiRmn, sigR) ;
  // Base parameter priors:
  phiDmn ~ normal(0,5) ;
  phiSmn ~ normal(0,2.5) ;
  phiRmn ~ normal(0,1) ;
  sigD ~ normal(0,2) ;
  sigS ~ normal(0,1) ;
  sigR ~ normal(0,1) ;
  sigD1 ~ normal(0,2) ;
  sigS1 ~ normal(0,1) ;
  sigR1 ~ normal(0,1) ;  
  tau ~ normal(35,2.5) ;
  invscale ~ normal(invscpri, .1*invscpri) ;
}
// Section 5. Derived parameters and statistics 
 generated quantities {
  real log_lik[3*Nyrs] ;        // Log liklihood of obs. data (for LooIC)
  real ynew[Nyrs]  ;            // New observations (out of sample) for ppc  
  vector[Nyrs] P_resid;         // Pearson residuals, observed data
  vector[Nyrs] P_resid_new;     // Pearson residuals, new data
  real Tstat ;                  // Test statistic, observed data (chi-2)
  real Tstat_new ;              // Test statistic, new data (chi-2)
  real ppp ;                    // posterior predictive P-value 
  int<lower=0> c;
  c = 0 ;
  for (t in 1:Nyrs) {
    // Summed Log-likelihood of observations
    c = c+1 ;
    log_lik[c] = gamma_lpdf(Obs_Tot[t] | D[t] * invscale, invscale) ;      
    c = c+1 ;
    log_lik[c] = dirichlet_lpdf(pie[t,] | tau * SzDst[t]) ;
    c = c+1 ;
    log_lik[c] = multinomial_lpmf(Obs_Sz[t,] | pie[t]) ;
    // Pearson residuals (for gamma distrib. densities only)
    ynew[t] = gamma_rng(D[t] * invscale, invscale) ;
    P_resid[t] = square(Obs_Tot[t] - D[t]) / (D[t] / invscale) ;  
    P_resid_new[t] = square(ynew[t] - D[t]) / (D[t] / invscale) ;  
  }
  Tstat = sum(P_resid) ;
  Tstat_new = sum(P_resid_new) ;
  ppp = Tstat > Tstat_new ? 1 : 0;
}
