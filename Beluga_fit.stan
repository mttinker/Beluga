// Beluga model
// This Stan code executes an age-structured population model for St Lawrence
//  Belugas, fit to surey data and stranding data 
//
functions {
  // user function to create demog. transition projection matrix
  matrix makemat(int NSt, real Pr, real Snn, real Sn, real Sy, real Sa) {
    // Snn = survival of neonatal calves from birth to Sept survey
    matrix[NSt,NSt] M = rep_matrix(0,NSt,NSt) ;
    M[1,11] = Sa * Snn ;
    M[2,1] = Sn ;
    M[3,2] = Sy ;
    M[4:8,3:7] = add_diag(M[4:8,3:7], rep_vector(Sa,5)) ;
    M[9,8] = Sa/2 ;    
    M[9,9] = Sa ;
    M[10,8] = Sa/2 * (1-Pr) ;
    M[11,8] = Sa/2 * Pr ;
    M[11,10] = Sa * Pr ;
    M[10,10] = Sa * (1-Pr) ;
    M[12,11] = Sa * Snn ;
    M[10,11] = Sa * (1-Snn) ;
    M[10,12] = Sa * (Sn + (1-Sn) * (1-Pr)) ;
    M[11,12] = Sa * (1-Sn) * Pr ;
    return M ;
  }
  // user function to create mortality projection matrix
  matrix makeDmat(int Nag, real Sn, real Sy, real Sa) {
    matrix[Nag,Nag] D = rep_matrix(0,Nag,Nag) ;
    D[1,1] = 1-Sn ;
    D[2,2] = 1-Sy ;
    D[3:Nag,3:Nag] = add_diag(D[3:Nag,3:Nag], rep_vector((1-Sa),(Nag-2) )) ;
    return D ;
  }
}
// Section 1. Data inputs for model
data {
  int<lower=1> Nyrs ;            // Years of dynamics 
  int<lower=1> NyrsH ;           // Years for which harvest occured
  int<lower=1> NStg ;            // Number stages
  int<lower=1> NAge ;            // Number age classes 
  int<lower=1> Nsrv ;            // Number surveys
  int<lower=1> Nstr ;            // Number stranding counts
  int<lower=1> YrSv[Nsrv] ;      // Year index for surveys
  int<lower=1> YrSt[Nstr] ;      // Year index for stranding counts
  vector<lower=0>[Nsrv] ObsS ;   // Observed survey estimates
  vector<lower=0>[Nsrv] invSc ;  // inverse scale (precision) for survey counts
  real<lower=0,upper=1> PJ[Nsrv];// Proportion juveniles in surveys
  int<lower=0> StrNB[Nstr] ;     // Stranding counts of newborns
  int<lower=0> StrOA[Nstr] ;     // Stranding counts of older animals
  //int<lower=0> Agects[Nstr,NAge];// matrix of adult counts by age class (cols) & year (rows)
  simplex[NStg] ssd ;            // Initial stage distribution
  vector<lower=0>[NyrsH] Harv ;  // Harvest totals by year
  real upsilon ;                 // precision param for ppn juveniles
}
// Section 2. The parameters to be estimated 
parameters {
  real<lower=0> logNinit ;      // initial log abundance
  real<lower=0,upper=1.5> gammA ;           // baseline log hazards + 3 (for adults)
  real<lower=0,upper=2> gamma_Y ;         // log hazard ratio, yearlings vs adults
  real<lower=0,upper=2> gamma_N ;         // log hazard ratio, newborns vs yearlings
  real alpha ;                    // logit param for pregancy rate
  real<lower=0,upper=3> phi ;             // density dependent param
  real<lower=0,upper=3.5> sig_P ;           // stochasticity in pregancy rates (logit)
  vector[Nyrs] epsP ;             // z-vec of devs for stoch. effect, Preg
  real<lower=0,upper=3.5> sig_N ;           // stochasticity in newborn survival (log haz)
  vector[Nyrs] epsN ;             // z-vec of devs for stoch. effect, NB
  real gamma_H_mn ;               // mean harvest log hazards  
  real<lower=0,upper=5> sig_H ;   // variance in harvest log hazards 
  vector[NyrsH] gamma_H ;         // annual harvest log hazards during early years  
  real<lower=0,upper=.5> PD_NB ;  // Probability detect newborn carcass
  real<lower=0,upper=.5> PD_OA ;  // Probability detect older animal carcass
  // real<lower=0.1,upper=15> upsilon ;// precision param for beta dist, proportion Juv
  //real<lower=0> tau ;             // precision param for dirichlet-multinomial age dist  
  //simplex[NAge] pie[Nstr];        // probs of stranding by age (for multinomial)
} 
// Section 3. Additional transformed parameters, including key model dynamics
transformed parameters {
  real gamma_A = gammA - 3;               // adult log hazards 
  real haz_A = exp(gamma_A) ;             // adult hazards (natural)
  real haz_Y = exp(gamma_A + gamma_Y) ;   // yearling hazards (natural)
  real S_A = exp(-haz_A) ;                // adult survival rate
  real S_Y = exp(-haz_Y) ;                // yearling survival rate
  // real Pr = inv_logit(3 + alpha) ;       // pregancy rate (if constant)
  vector[Nyrs] Pr ;                       // pregancy rate (if variable)
  vector[Nyrs] eps_P = epsP * sig_P;      // stoch. effects, Preg
  vector[Nyrs] eps_N = epsN * sig_N;      // stoch. effects, NB
  vector[NStg] n[Nyrs] ;
  vector[Nyrs] N ;
  vector[Nyrs] ppnJ ;
  vector[Nyrs] nd_NB ;
  vector[Nyrs] nd_OA ;
  vector[NyrsH] Harv_expect ;
  //simplex[NAge] Agvc[Nyrs] ;      // predicted carcass age dist, all years
  // Year 1
  N[1] = exp(8.5+logNinit) ;
  n[1] = N[1] * ssd ;
  Harv_expect[1] = Harv[1] ;
  Pr[1] = inv_logit(3 + alpha + eps_P[1]) ;
  ppnJ[1:NyrsH] = rep_vector(0,NyrsH) ;
  // Loop through early years with harvest
  for(t in 2:NyrsH){
    real haz_H ;
    real haz_N ;
    real S_N ;
    real S_Yh ;
    real S_Ah ;
    vector[NAge] nt ;
    vector[NAge] nd ;
    matrix[NStg,NStg] M;
    matrix[NAge,NAge] D;
    Pr[t] = inv_logit(3 + alpha + eps_P[t]) ;
    haz_H = exp(gamma_H[t]) ;
    S_Yh = exp(-1 * (haz_Y + haz_H) ) ;
    S_Ah = exp(-1 * (haz_A + haz_H) ) ;
    haz_N = exp(gamma_A + gamma_Y + gamma_N + phi * (N[t-1]/10000) + eps_N[t]) ;
    S_N = exp(-haz_N) ;
    M = makemat(NStg,Pr[t],sqrt(S_N),S_N,S_Yh,S_Ah) ;
    n[t] = M * n[t-1]  ;
    N[t] = sum(n[t]) ;
    nt[1:8] = n[t-1][1:8] ;  
    nt[9] = sum(n[t-1][9:12]) ;  
    D = makeDmat(NAge,S_N,S_Yh,S_Ah) ;
    nd = D * nt ;
    Harv_expect[t] = (nd[2] * (haz_H/(haz_Y + haz_H))) + 
                     (sum(nd[3:9]) * (haz_H/(haz_A + haz_H))) ; 
  }  
  // loop through remaining years (no harvest)
  for(t in (NyrsH+1):Nyrs){
    real haz_N ;
    real S_N ;
    vector[NAge] nt ;
    vector[NAge] nd ;
    matrix[NStg,NStg] M;
    matrix[NAge,NAge] D;
    Pr[t] = inv_logit(3 + alpha + eps_P[t]) ;
    haz_N = exp(gamma_A + gamma_Y + gamma_N + phi * (N[t-1]/10000) + eps_N[t]) ;
    S_N = exp(-haz_N) ;
    M = makemat(NStg,Pr[t],sqrt(S_N),S_N,S_Y,S_A) ;
    n[t] = M * n[t-1]  ;
    N[t] = sum(n[t]) ;
    ppnJ[t] = (0.0001 + sum(n[t][1:2])) / (0.0001 +  N[t]) ;
    nt[1:8] = n[t-1][1:8] ;  
    nt[9] = sum(n[t-1][9:12]) ;  
    D = makeDmat(NAge,S_N,S_Y,S_A) ;
    nd = D * nt ;
    nd[1] = nd[1] + n[t-1][11] * S_A * (1-sqrt(S_N)) ;
    nd_NB[t] = nd[1]+.001 ;
    nd_OA[t] = sum(nd[2:NAge])+.001 ;
    //Agvc[t] = (nd + .00000001) / (sum(nd) + NAge * .00000001 ) ; 
  }
}
// Section 4. Estimating model parameters (drawing from probability distributions)
model {
  // A) Observed nodes:
  // - Harvest numbers
  Harv ~ gamma(Harv_expect * 100, 100) ;
  // -Survey Estimates, by year. NOTE for gamma, invscale = M / V
  ObsS ~ gamma(N[YrSv] .* invSc, invSc) ; 
  // -Proportion juvs (<2yr-olds) 
  PJ ~ beta(upsilon*ppnJ[YrSv]+.0001,upsilon*(1-ppnJ[YrSv])+.0001 ) ;
  // -Stranding counts, newborns
  StrNB ~ poisson( PD_NB * nd_NB[YrSt]) ;     
  // -Stranding counts, older animals
  StrOA ~ poisson( PD_OA * nd_OA[YrSt]) ;     
  // -Age vectors of stranded animals
  //for(i in 1:Nstr){
    // NOTE: Use dirichlet-multinomial to handle error/variance in age counts 
  //  pie[i] ~ dirichlet(10 * tau * Agvc[YrSt[i]]);
  //  Agects[i,] ~ multinomial(pie[i]) ;
  //}
  //
  // B) Prior distributions for model parameters:
  // Hierarchical random effects:  
  epsP ~ normal(0, 1) ;
  epsN ~ normal(0, 1) ;
  gamma_H ~ normal(gamma_H_mn-5,sig_H) ;
  // Base parameter priors:
  gamma_H_mn ~ normal(0,1) ;
  sig_H ~ normal(0,2.5) ;
  sig_P ~ normal(0,2.5) ;
  sig_N ~ normal(0,2.5) ;
  logNinit ~ normal(0,1) ;
  gammA ~ normal(0,1) ;
  gamma_Y ~ normal(0,1) ;
  gamma_N ~ normal(0,1) ;
  alpha ~ normal(0,1) ;
  phi ~ normal(0,1) ;
  PD_NB ~ beta(1,5) ;
  PD_OA ~ beta(1,5) ;
  // upsilon ~ cauchy(0,2.5) ;
  // tau ~ cauchy(0,2.5) ;
}
// Section 5. Derived parameters and statistics 
 generated quantities {
  real ynew[Nsrv]  ;            // New observations (out of sample) for ppc  
  vector[Nsrv] P_resid;         // Pearson residuals, observed data
  vector[Nsrv] P_resid_new;     // Pearson residuals, new data
  real Tstat ;                  // Test statistic, observed data (chi-2)
  real Tstat_new ;              // Test statistic, new data (chi-2)
  real ppp ;                    // posterior predictive P-value 
  for (t in 1:Nsrv) {
    // Pearson residuals (for gamma distrib. densities only)
    ynew[t] = gamma_rng(N[YrSv[t]] * invSc[t], invSc[t]) ;
    P_resid[t] = square(ObsS[t] - N[YrSv[t]]) / (N[YrSv[t]]/ invSc[t]) ;  
    P_resid_new[t] = square(ynew[t] - N[YrSv[t]]) / (N[YrSv[t]]/ invSc[t]) ;  
  }
  Tstat = sum(P_resid) ;
  Tstat_new = sum(P_resid_new) ;
  ppp = Tstat > Tstat_new ? 1 : 0;
}
