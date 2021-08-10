# Process Beluga stranding age data
Year1str = min(df_Strd$YYYY)
StrNB = numeric()
StrOA = numeric()
YrSt = numeric()
Agects = matrix(0,nrow = YearT - Year1str + 1, ncol = NAge - 1)
t = 0
for(y in Year1str:YearT){
  t = t+1
  yr = y - Year1 + 1
  YrSt[t] = yr
  # For newborns: this years births plus newborns who were alive at last survey
  ii = which( (df_Strd$YYYY==y & df_Strd$MM<9 & df_Strd$Best_GLG==0) | 
              (df_Strd$YYYY==y & df_Strd$MM>2 & df_Strd$MM<9 & df_Strd$Best_GLG==1) |
              (df_Strd$YYYY==(y-1) & df_Strd$MM>8 & df_Strd$Best_GLG==0))
  StrNB[t] = length(ii)
  # For ages 1 and up: includes age t in Sept to age t+1 following August
  ii = which( (df_Strd$YYYY==y & df_Strd$MM>2 & df_Strd$MM<9 & 
                   (is.na(df_Strd$Best_GLG) | df_Strd$Best_GLG > 1) ) |
                (df_Strd$YYYY==y & df_Strd$MM<3 & (is.na(df_Strd$Best_GLG) | df_Strd$Best_GLG > 0) ) |
                (df_Strd$YYYY==(y-1) & df_Strd$MM>8 & (is.na(df_Strd$Best_GLG) | df_Strd$Best_GLG > 0)))
  StrOA[t] = length(ii)
  # Loop through age classes to calculate # stranded by age
  for (a in 1:(NAge-1)){
    if(a<(NAge-1)){
      ii = which( (df_Strd$YYYY==y & df_Strd$MM>2 & df_Strd$MM<9 & df_Strd$Best_GLG == a+1) | 
                    (df_Strd$YYYY==y & df_Strd$MM<3 & df_Strd$Best_GLG == a) |
                    (df_Strd$YYYY==(y-1) & df_Strd$MM>8 & df_Strd$Best_GLG == a))     
    }else{
      ii = which( (df_Strd$YYYY==y & df_Strd$MM>2 & df_Strd$MM<9 & df_Strd$Best_GLG >= a+1) | 
                  (df_Strd$YYYY==y & df_Strd$MM<3 & df_Strd$Best_GLG >= a) |
                  (df_Strd$YYYY==(y-1) & df_Strd$MM>8 & df_Strd$Best_GLG >= a))     
    }
    Agects[t,a] = length(ii)
  }
}
Nstr = length(YrSt)
