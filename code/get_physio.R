# Get physiological traits
get_physio = function(photo_data, resp_data, T_tissue_data) {
  
  results = data.frame()
  
  for(i in unique(photo_data$curveID)) {
  
  # Grab a curve
    curAT = subset(photo_data, curveID == i)
  # Find growth T for curve and sat PPFD value
    species = unique(curAT$taxon)
    curtreatment = unique(curAT$treatment)
    tissueT = subset(T_tissue_data, taxon == species & treatment == curtreatment)$T_tissue
    PPFD = round(mean(curAT$Qin))
    
  # Start by fitting a model to each AT curve
  # Extract A25, Agrowth, Amax
    z.A = smooth.spline(x = curAT$Tleaf, y = curAT$Photo, nknots = 5)
    
    Sx = seq(min(curAT$Tleaf), max(curAT$Tleaf), 0.01)
    
    #plot(Sx, predict(z.A, Sx)$y)
    #points(curAT$Tleaf, curAT$Photo, col="red")
    
    Amax <- max(predict(z.A, Sx)$y)
    A25 = predict(z.A, 25)$y
    Agrowth = predict(z.A, tissueT)$y
    
  # Fit a model to gs-T
  # Extract gs25, gs_growth
    z.gs = smooth.spline(x = curAT$Tleaf, y = curAT$gsw, nknots = 5)
    
    #plot(Sx, predict(z.gs, Sx)$y)
    #points(curAT$Tleaf, curAT$gsw, col="red")
    
    gs25 = predict(z.gs, 25)$y
    gsgrowth = predict(z.gs, tissueT)$y

  # Fit a model to RT
  # Extract R25, Rgrowth
  # Note that we don't have a corresponding RT curve for each at curve so we either
  # have to use means here or skip this step
    curRT = subset(resp_data, taxon == species & treatment == curtreatment)
    z.R = smooth.spline(x = curRT$Tleaf, y = curRT$Photo, nknots = 5)
    
    Sx = seq(min(curRT$Tleaf), max(curRT$Tleaf), 0.01)
    #plot(Sx, predict(z.R, Sx)$y)
    #points(curRT$Tleaf, curRT$Photo, col="red")
    
    Rd25 = predict(z.R, 25)$y
    Rdgrowth = predict(z.R, tissueT)$y
  
  # Fit a model to Ci-T
  # Extract Ci25, Cigrowth
    z.Ci = smooth.spline(x = curAT$Tleaf, y=curAT$Ci, nknots = 5)
    #plot(Sx, predict(z.Ci, Sx)$y)
    #points(curAT$Tleaf, curAT$Ci, col="red")
    Ci25 = predict(z.Ci, 25)$y
    Cigrowth = predict(z.Ci, tissueT)$y
        
  # Use fitaci onepoint method
  # Get Vcmax25, Vcmaxgrowth
    Vcmax25 = fitaci(data.frame(Photo = A25, Tleaf = 25, Ci = Ci25, PARi = PPFD, Rd = Rd25),
                     fitmethod = "onepoint", Tcorrect = F)$Vcmax
    Vcmaxgrowth = fitaci(data.frame(Photo = Agrowth, Tleaf = tissueT, Ci = Cigrowth, PARi = PPFD, Rd = Rdgrowth),
                         fitmethod = "onepoint", Tcorrect = F)$Vcmax
    
  # Compute iWUE25, iWUEgrowth
    iWUE25 = A25/gs25
    iWUEgrowth = Agrowth/gsgrowth
    
  # Pull Topt_A, E_A, and E_R from models
  
  # Assemble dataframe
    results = bind_rows(results, 
                        data.frame(curveID = i, taxon = species, treatment = curtreatment, tissueT,
                                   Agrowth, A25, Amax,
                                   Rdgrowth, Rd25,
                                   gsgrowth, gs25,
                                   Cigrowth, Ci25,
                                   iWUEgrowth, iWUE25,
                                   Vcmaxgrowth, Vcmax25))
  }
  return(results)
}


