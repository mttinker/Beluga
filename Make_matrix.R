Make_matrix <- function(Snn,Sn,Sy,Sa,Pr) {
  M = matrix(0, nrow = 12,ncol = 12)
  # or STAN version: matrix[12,12] M; M = rep_matrix(0,12,12) ;
  M[1,11] = Sa * Snn
  M[2,1] = Sn 
  M[3,2] = Sy
  M[4:8,3:7] = diag(1,5,5)*rep(Sa,5) 
  # STAN version: M[4:8,3:7] = add_diag(M[4:8,3:7], rep_vector(Sa,5)) 
  #  or repeat 5 steps, e.g.: M[4,3] = Sa
  M[9,8] = Sa/2     
  M[9,9] = Sa
  M[10,8] = Sa/2 * (1-Pr)
  M[11,8] = Sa/2 * Pr  
  M[11,10] = Sa * Pr
  M[10,10] = Sa * (1-Pr)
  M[12,11] = Sa * Snn
  M[10,11] = Sa * (1-Snn)  
  M[10,12] = Sa * (Sn + Sn * (1-Pr))
  M[11,12] = Sa * (1-Sn) * Pr
  return(M)
}