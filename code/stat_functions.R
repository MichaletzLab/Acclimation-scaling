# stat_functions.R - functions fitting the best model for each of photo, resp, and growth

# Fits best NLME model to respiration as a function of temperature
fit.resp.nlme = function(resp_data) {
  
  resp_data$taxon = as.factor(resp_data$taxon)
  
  # I'm going with this model - the model with no species effect on E was preferred very slightly
  # by AIC, but not significantly different according to LRT, but I would like species level
  # estimates of E, so it seems reasonable.
  
  resp.nlme = nlme(Photo ~ arrhenius(Tleaf, J_ref, E, T_ref = 10),
                     data = resp_data,
                     fixed = J_ref + E ~ taxon - 1,
                     random = J_ref + E ~ 1,
                     groups = ~ curveID,
                     weights = varPower(),
                     method = "REML",
                     start = c(0.8, 0.6,0.8, 0.6,0.8, 0.6))
  
  return(resp.nlme)

}


fit.photo.nlme = function(photo_data) {
  
  photo_data$taxon = as.factor(photo_data$taxon)
  
  photo.nlme.1 = nlme(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                      data = photo_data,
                      fixed = r_tref + e + eh + topt ~ 1,
                      random = pdDiag(r_tref + e + eh + topt ~ 1),
                      groups = ~ curveID,
                      start = c(r_tref = 21, e = 0.71, eh = 1.46, topt = 22.7),
                      na.action = na.omit,
                      method = "ML",
                      verbose = F)
  
  photo.nlme.2 = update(photo.nlme.1, 
                        fixed = r_tref + e + eh + topt ~ taxon - 1,
                        start = c(14,14,14,0.5,0.5,0.5,1,1,1,25,25,25))
  
  # photo.nlme.3c = update(photo.nlme.1,
  #                        fixed = list(eh ~ 1, r_tref + e + topt ~ taxon - 1),
  #                        start = c(1.1, 14.8, 10.5, 22.9, 0.67,0.29,0.79, 23.4,24.2,21.4)) 
  # 
  photo.nlme.6 = nlme(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                      data = photo_data,
                      fixed = r_tref + e + eh + topt ~ taxon - 1,
                      #start = c(12.5, 10.5, 22.7, 0.51,0.24,0.78,1.1,1.5,1.1,23.5,24.9,21.3),
                      start = list(fixed = fixef(photo.nlme.2), random = ranef(photo.nlme.2)),
                      random = list(r_tref + e + eh + topt ~ 1),
                      groups = ~ curveID,
                      #start = c(r_tref = 21, e = 0.71, eh = 1.46, topt = 22.7),
                      na.action = na.omit,
                      method = "ML",
                      verbose = T,
                      control = nlmeControl(pnlsMaxIter = 50, msMaxIter = 100, maxIter = 100)
  )
  
#   # NOW WHY ISN"T THIS CONVERGING???
#   photo.nlme.7c = update(photo.nlme.6,
#                          fixed = list(eh ~ 1, r_tref + e + topt ~ taxon - 1),
# #                         start = list(fixed = fixef(photo.nlme.3c)))
#                          start = c(1.12, 11.24, 12.18, 22.65, 0.51, 0.28, 0.79, 22.65, 25.5, 21.2))
#   
  photo.nlme = update(photo.nlme.6, 
                      start = list(fixed = fixef(photo.nlme.6), random = ranef(photo.nlme.6)),
                      method = "REML")

  
  return(photo.nlme)
}


fit.photo.nlme.with.35 = function(photo_data, photo_data_unfiltered, photo.nlme) {

  photo_all = bind_rows(photo_data, subset(photo_data_unfiltered, treatment == 35))
  
  photo_all$taxon = as.factor(photo_all$taxon)
  # 
  # photo.nlme.1 = nlme(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
  #                     data = photo_all,
  #                     fixed = r_tref + e + eh + topt ~ 1,
  #                     random = pdDiag(r_tref + e + eh + topt ~ 1),
  #                     groups = ~ curveID,
  #                     start = c(r_tref = 21, e = 0.71, eh = 1.46, topt = 22.7),
  #                     na.action = na.omit,
  #                     method = "ML",
  #                     verbose = F)
  # 
  # photo.nlme.2 = update(photo.nlme.1, 
  #                       fixed = r_tref + e + eh + topt ~ taxon - 1,
  #                       start = c(15,10,23,0.5,0.24,0.79,1.1,1.5,1.1,20,25,20),
  #                       #control = nlmeControl(pnlsMaxIter = 50, msMaxIter = 100, maxIter = 100),
  #                       verbose = T)
  
  # photo.nlme.3c = update(photo.nlme.1,
  #                        fixed = list(eh ~ 1, r_tref + e + topt ~ taxon - 1),
  #                        start = c(1.1, 14.8, 10.5, 22.9, 0.67,0.29,0.79, 23.4,24.2,21.4)) 
  # 
  photo.nlme.6 = nlme(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                      data = photo_all,
                      fixed = r_tref + e + eh + topt ~ taxon - 1,
                      #start = c(12.5, 10.5, 22.7, 0.51,0.24,0.78,1.1,1.5,1.1,23.5,24.9,21.3),
                      start = list(fixed = fixef(photo.nlme)),
                      random = list(r_tref + e + eh + topt ~ 1),
                      groups = ~ curveID,
                      #start = c(r_tref = 21, e = 0.71, eh = 1.46, topt = 22.7),
                      na.action = na.omit,
                      method = "REML",
                      verbose = T,
                      control = nlmeControl(pnlsMaxIter = 50, msMaxIter = 100, maxIter = 100)
  )
  
  #   # NOW WHY ISN"T THIS CONVERGING???
  #   photo.nlme.7c = update(photo.nlme.6,
  #                          fixed = list(eh ~ 1, r_tref + e + topt ~ taxon - 1),
  # #                         start = list(fixed = fixef(photo.nlme.3c)))
  #                          start = c(1.12, 11.24, 12.18, 22.65, 0.51, 0.28, 0.79, 22.65, 25.5, 21.2))
  #   
  # photo.nlme.final = update(photo.nlme.6, 
  #                     start = list(fixed = fixef(photo.nlme.6), random = ranef(photo.nlme.6)),
  #                     method = "REML")
  # 
  # 
  # 
  
  return(photo.nlme.6)
}

# Fits best NLME model to biomass as a fn of time and temperature
fit.bm.nlme = function(bm_data, T_tissue_data) {
  biomass = bm_data %>% rename(treatment = temperature_C, taxon = species)
  biomass_thermal = merge(biomass, T_tissue_data, by=c("taxon", "treatment"))
  bm = biomass_thermal %>% select(taxon, treatment, temp = T_tissue, mass = biomass_total, time = elapsed_days)
  
  bm$taxon = as.factor(bm$taxon)
  bm$treatment = as.factor(bm$treatment)
  
  bm.nlme = nlme(mass ~ biomass_model(temp, time, J_ref, E, E_D, T_opt, M_0, alpha, tref = 10),
                   data = bm,
                   fixed = J_ref + E + E_D + T_opt + M_0 + alpha ~ taxon - 1,
                   random = M_0 ~ 1,
                   groups = ~ treatment,
                   weights = varPower(),
                   method = "REML",
                   start = c(0.06,0.06,0.06,0.3,0.3,0.3,3,3,3,21,21,21,1,1,1,0.9,0.9,0.9))
  
  return(bm.nlme)
}

# Model of biomass growth as a function of time and temperature
biomass_model = function(temp, time, r_tref, e, eh, topt, M_0, beta, tref=10) {
  
  r = pawar_2018(temp, r_tref, e, eh, topt, tref)
  biomass = ( M_0^(1-beta) + r*time*(1-beta) )^(1/(1-beta))
  return(biomass)
  
}


# Returns Arrhenius function for the given parameters
# Borrowed from curve_fitting_michaletz_2021
arrhenius <- function(temp, J_ref, E, T_ref = 25) { 
  # temp   : Temperature values to evaluate function at (C)
  # J_ref  : Rate at T_ref (units depend on rate)
  # E      : Activation energy (eV; E > 0)
  # T_ref  : Reference temperature for normalization (C)
  
  k <- 8.62e-5           # Boltzmann's constant (eV K-1)
  temp = temp + 273.15   # Convert to Kelvin
  T_ref = T_ref + 273.15 # Convert to Kelvin
  
  # Evaluate and return
  return( J_ref * exp( E * (1/(k*T_ref) - 1/(k*temp)) ) ) 
}

