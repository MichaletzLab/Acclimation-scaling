# Functions for making model predictions of acclimation effects on
# temperature response parameters

# Kattge and Knorr 2007 model
K2007 = function(Vcmax25, Tgrowth) {
  
  PPFDin = 600
  
  rJV = 2.59 - 0.035*Tgrowth # Fix
  Jmax25 = Vcmax25*rJV
  
  Ea_V = 71.513
  Ed_V = 200
  dS_V = 668.39 - 1.07*Tgrowth
  
  Ea_J = 49.884
  Ed_J = 200
  dS_J = 659.70 - 0.75*Tgrowth
  
  dat = Photosyn(VPD = 1,
                 Ca = 420,
                 PPFD = PPFDin,
                 Tleaf = Tgrowth,
                 Patm = 100,
                 Jmax = Jmax25,
                 Vcmax = Vcmax25,
                 EaV = Ea_V*1000,
                 EdVC = Ed_V*1000,
                 delsC = dS_V,
                 EaJ = Ea_J*1000,
                 EdVJ = Ed_J*1000,
                 delsJ = dS_J)
}

# Kumarathunge et al 2018 model
K2018 = function(Vcmax25, Thome, Tgrowth) {
  
  PPFDin = 600
  
  JVr = 2.56 - 0.0375*Thome - 0.0202*(Tgrowth - Thome)
  Jmax25 = JVr*Vcmax25
  
  Ea_V = 42.6 + 1.14*Tgrowth
  Ea_J = 40.71
  
  dS_V = 645.13-0.38*Tgrowth
  dS_J = 658.77 - 0.84*Thome - 0.52*(Tgrowth - Thome)
  
  Ed_V = 200 #kJ/mol
  Ed_J = 200 #kJ/mol
  
  dat = Photosyn(VPD = 1,
                 Ca = 420,
                 PPFD = PPFDin,
                 Tleaf = Tgrowth,
                 Patm = 100,
                 Jmax = Jmax25,
                 Vcmax = Vcmax25,
                 EaV = Ea_V*1000,
                 EdVC = Ed_V*1000,
                 delsC = dS_V,
                 EaJ = Ea_J*1000,
                 EdVJ = Ed_J*1000,
                 delsJ = dS_J)
  
  return(dat)
}


# Smith 2018 / 2019? model
S2018 = function(temps) {
  
  PPFDin = 600
  PPFDgrowth = 250
  
  res = calc_optimal_vcmax(pathway = "C3",     # Photosynthetic pathway
                           tg_c = temps,       # Growth temp
                           z = 85,             # Elevation (m)
                           vpdo = 1,#vpds,        # VPD at sea level
                           cao = 420,          # Atmospheric CO2
                           oao = 209460,       # Atmospheric O2
                           paro = PPFDgrowth,  # PAR (use growth PAR or sat PAR?) 
                           beta = 146,         # Cost parameter - use this standard value
                           theta = 0.85,       # Using assumed value
                           chi = NA,           # Not needed
                           q0 = 0.257,         # Using assumed value
                           q0_resp = "yes",    # Using assumed value
                           q0_int = -0.0805,   # Using assumed value
                           lma = NA,           # LMA (this doesn't affect A a all)
                           f = 1)              # Fraction of year in growing season (_)
  
  # At very low temperatures, S19 model predicts negative Vcmax and Jmax;
  # This is clearly not physically possible, so I am restricting the range
  # to > 0
  res = subset(res, vcmax > 0 & jmax > 0)
  
  dat = Photosyn(VPD = 1,
                 Ca = 420,
                 PPFD = PPFDin,
                 Tleaf = res$tg_c,#Tgrowth,
                 Patm = 100,
                 Jmax = res$jmax,#Jmax25,
                 Vcmax = res$vcmax,#Vcmax25,
                 GammaStar = res$gammastar,
                 Km = res$km,
                 Ci = 420*res$chi,
                 #EaV = 0, #Ea_V*1000,
                 #EdVC = 0, #Ed_V*1000,
                 #delsC = 0, #dS_V,
                 #EaJ = 0,#Ea_J*1000,
                 #EdVJ = 0,#Ed_J*1000,
                 #delsJ = 0,#dS_J,
                 Tcorrect = F)
  return(dat)
  
}

