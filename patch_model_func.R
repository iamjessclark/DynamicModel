# model function # patch model 


#library(deSolve)
#library(tidyverse)

# model function ####

RVFV_patch <- function(t, Y, pars){
  
  # states patch 1
  S_H1 <- Y[1]
  E_H1 <- Y[2]
  I_H1 <- Y[3]
  R_H1 <- Y[4]
  
  Eggs1 <- Y[5]
  Pupae1 <- Y[6]
  
  S_Va1 <- Y[7]
  E_Va1 <- Y[8]
  I_Va1 <- Y[9]
  
  # states patch 2
  S_H2 <- Y[10]
  E_H2 <- Y[11]
  I_H2 <- Y[12]
  R_H2 <- Y[13]
  
  Eggs2 <- Y[14]
  Pupae2 <- Y[15]
  
  S_Va2 <- Y[16]
  E_Va2 <- Y[17]
  I_Va2 <- Y[18]
  
  # patch 1 pars
  kv1 <- pars[1] 
  Reggs1 <- pars[2] 
  Heggs1 <- pars[3] # hatching rate eggs
  meggs1 <- pars[4] # death eggs 
  Peggs1 <- pars[5] # probability of survival eggs
  mpup1 <- pars[6] # deaths pupae
  pupmaturation1 <- pars[7] 
  pupsurvival1 <- pars[8]
  mv1 <- pars[9]
  lv1 <- pars[10]
  probVinfected1 <- pars[11]
  biterate1 <- pars[12]
  
  y1  <- pars[13]
  lh1 <- pars[14]
  mh1 <- pars[15]
  xh1 <- pars[16]
  kh1 <- pars[17]
  probHinfected1 <- pars[18]
  ah1 <- pars[19]
  
  # patch 2 pars
  kv2 <- pars[20]
  Reggs2 <- pars[21]
  Heggs2 <- pars[22]
  Peggs2 <- pars[23]
  meggs2 <- pars[24]
  mpup2 <- pars[25]
  pupmaturation2 <- pars[26]
  pupsurvival2 <- pars[27]
  mv2 <- pars[28]
  lv2 <- pars[29]
  probVinfected2 <- pars[30]
  biterate2 <- pars[31]
  
  y2  <- pars[32]
  lh2 <- pars[33]
  mh2 <- pars[34]
  xh2 <- pars[35]
  kh2 <- pars[36]
  probHinfected2 <- pars[37]
  ah2 <- pars[38]
  
  pH12 <- pars[39]
  pH21 <- pars[40]
  
  dYdt <- vector(length=18)
  
  # carrying capacity terms
  Kh1 <- (1-((S_H1 + E_H1 + I_H1 + R_H1)/kh1))
  Kh2 <- (1-((S_H2 + E_H2 + I_H2 + R_H2)/kh2))
  Kv1 <- (1-(Pupae1/kv1))
  Kv2 <- (1-(Pupae2/kv2))
  
  # beta terms
  BH1 <- biterate1*probHinfected1
  BH2 <- biterate2*probHinfected2
  BV1 <- biterate1*probVinfected1
  BV2 <- biterate2*probVinfected2
  
  Hostp1_births = (xh1*(S_H1 + E_H1 + I_H1 + R_H1)*Kh1) 
  hostp1_infection = (BH1*(S_H1/(S_H1 + E_H1 + I_H1 + R_H1))*I_Va1)
  hostp1_progression = lh1*E_H1
  I_Hp1_mort = ah1*I_H1
  I_Hp1_clear = y1*I_H1
  
  Hostp2_births = (xh2*(S_H2 + E_H2 + I_H2 + R_H2)*Kh2) 
  hostp2_infection = (BH2*(S_H2/(S_H2 + E_H2 + I_H2 + R_H2))*I_Va2)
  hostp2_progression = lh2*E_H2
  I_Hp2_mort = ah2*I_H2
  I_Hp2_clear = y2*I_H2
  
  Eggproductionp1 = Reggs1*(S_Va1 + E_Va1 + I_Va1)
  Hatchingp1 = Eggs1*Heggs1*Peggs1*Kv1
  Maturationp1 = Pupae1*pupmaturation1*pupsurvival1
  Vecp1_infection = (BV1*(I_H1/(S_H1 + E_H1 + I_H1 + R_H1))*S_Va1) 
  Vecp1_progression = lv1*E_Va1
  
  Eggproductionp2 = Reggs2*(S_Va2 + E_Va2 + I_Va2)
  Hatchingp2 = Eggs2*Heggs2*Peggs2*Kv2
  Maturationp2 = Pupae2*pupmaturation2*pupsurvival2
  Vecp2_infection = (BV2*(I_H2/(S_H2 + E_H2 + I_H2 + R_H2))*S_Va2) 
  Vecp2_progression = lv2*E_Va2
  
  #Hosts p1
  dYdt[1] <- Hostp1_births  - hostp1_infection - mh1*S_H1 + ((pH21*S_H2)-(pH12*S_H1)) # susceptible
  dYdt[2] <- hostp1_infection - hostp1_progression - mh1*E_H1 + ((pH21*E_H2)-(pH12*E_H1)) # exposed
  dYdt[3] <- hostp1_progression - I_Hp1_mort - I_Hp1_clear - mh1*I_H1 + ((pH21*I_H2)-(pH12*I_H1)) # infected 
  dYdt[4] <- I_Hp1_clear - mh1*R_H1 + ((pH21*R_H2)-(pH12*R_H1)) # recovered
  
  #Vectors p1
  dYdt[5] <- -Eggs1*meggs1
  dYdt[6] <- - Pupae1*mpup1
  dYdt[7] <- - Vecp1_infection - (mv1*S_Va1) # susceptible
  dYdt[8] <- Vecp1_infection - Vecp1_progression - (mv1*E_Va1) # exposed
  dYdt[9] <- Vecp1_progression - (mv1*I_Va1) # infected
  
  #H2
  dYdt[10] <- Hostp2_births  - hostp2_infection - mh2*S_H2 + ((pH12*S_H1)-(pH21*S_H2))
  dYdt[11] <- hostp2_infection - hostp2_progression - mh2*E_H2 + ((pH12*E_H1)-(pH21*E_H2))
  dYdt[12] <- hostp2_progression - I_Hp2_mort - I_Hp2_clear - mh2*I_H2 + ((pH12*I_H1)-(pH21*I_H2))
  dYdt[13] <- I_Hp2_clear - mh2*R_H2 + ((pH12*R_H1)-(pH21*R_H2))
  
  #V2
  dYdt[14] <- Eggproductionp2 - Hatchingp2 - Eggs2*meggs2
  dYdt[15] <- Hatchingp2  - Maturationp2 - Pupae2*mpup2
  dYdt[16] <- Maturationp2  - Vecp2_infection - (mv2*S_Va2)
  dYdt[17] <- Vecp2_infection - Vecp2_progression - (mv2*E_Va2)
  dYdt[18] <- Vecp2_progression - (mv2*I_Va2)
  
  
  return(list(dYdt))
  
}

# forcing functions ####

# hatching
eventhatch <- function(t, Y, parms) {
  with(as.list(c(Y, parms)),                 #local environment to evaluate derivatives
       {
         #seasonal egg hatching in rural environments (population 1)
         Y[5] <- Y[5] + ((Y[7]+Y[8]+Y[9])*Reggs1) - Y[5]*Heggs1 # rate eggs are laid - hatching
         Y[6] <- Y[6] + Y[5]*Heggs1 - Y[6]*pupmaturation1*pupsurvival1  # hatching - maturation 
         Y[7] <- Y[7] + Y[6]*pupmaturation1*pupsurvival1
         return(Y)
       }
  )
}
