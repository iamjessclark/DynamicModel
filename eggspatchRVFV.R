#source("parameters.R")
#source("rain_data.R")
#source("patch_model_func.R")

# Run Model ####
SolPatch <- as.data.frame(ode(init,t, RVFV_patch, pars, 
                              events = list(func = eventhatch, 
                                            time = hatchtimes)))


# Plot outputs ####

colnames(SolPatch) <- c("time", "Sh1", "Eh1", "Ih1", "Rh1", 
                        "Eggs1", "Pupae1", "Sv1", "Ev1", "Iv1",  
                        "Sh2", "Eh2", "Ih2", "Rh2",
                        "Eggs2", "Pupae2", "Sv2", "Ev2", "Iv2")

SolPatch <- SolPatch %>%
  mutate(totalh1 = (Sh1 + Eh1 + Ih1 + Rh1), 
         totalh2 = (Sh2 + Eh2 + Ih2 + Rh2), 
         totalv1 = (Sv1 + Ev1 + Iv1), 
         totalv2 = (Sv2 + Ev2 + Iv2)) %>% 
  pivot_longer(-time, names_to = "class", values_to = "value") %>%
  mutate(Agent = ifelse(class == "Sh1"|class == "Sh2"|
                                class == "Eh1"|class == "Eh2"|
                                class == "Ih1"|class == "Ih2"|
                                class == "Rh1"|class == "Rh2"|
                                class=="totalh1"| class=="totalh2", "host", "vector"), 
        Status = factor(case_when(class == "Sh1"|class == "Sh2"| class == "Sv1"| class == "Sv2" ~ "susceptible", 
                                  class == "Eh1"|class == "Eh2"| class == "Ev1"| class == "Ev2" ~ "exposed", 
                                  class == "Ih1"|class == "Ih2"| class == "Iv1"| class == "Iv2" ~ "infected", 
                                  class == "Rh1"|class == "Rh2" ~ "recovered", 
                                  class== "totalh1"| class=="totalh2"| class == "totalv1"| class == "totalv2" ~ "total", 
                                  class == "Eggs1"|class == "Eggs2" ~ "Eggs", 
                                  class == "Pupae1" | class == "Pupae2" ~ "Pupae"), 
                               levels = c("Eggs", "Pupae", "susceptible", "exposed", "infected", "recovered", "total")), 
        Patch = ifelse(class == "Sh1"| class == "Eh1"| class == "Ih1"| class ==  "Rh1"|
                                class == "Sv1"| class == "Ev1"| class == "Iv1"| class == "Eggs1" | class == "Pupae1" |
                                class == "totalh1" | class == "totalv1", "rural", "peri-urban"))
