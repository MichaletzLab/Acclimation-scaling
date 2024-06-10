library(progress)
library(rTPC)

library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(scales)
library(gridExtra)
library(grDevices)

fit_curves_all <- function(Data, x = "Tleaf", y = "Photo", tref = 25, PLOT=TRUE) {

  set.seed(1)
  Data$tref = tref

  # Standardize names of x and y for use throughout
  if(x != "Tleaf") {
    Data$Tleaf = Data[[x]] # Rename X as Tleaf
    Data[[x]] = NULL      # Get rid of originals
  }
  if(y != "Photo") {
    Data$Photo = Data[[y]] # Rename Y as Photo
    Data[[y]] = NULL      # Get rid of originals
  }

  Data = subset(Data, Photo > 0)

  # Open output file if plotting desired
  if(PLOT==TRUE)
    grDevices::pdf("plots_all.pdf", width = 10, height = 5)

  #Data = assign_curve_type(Data)
  Data$curve_type = "Schoolfield"
  Data = Data %>% arrange(by_group = curveID)
  results = c()

  ######################
  # Curve fitting loop #
  ######################
  message("Fitting curves; this may take some time.")
  n_iter <- length(unique(Data$curveID))
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = n_iter,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar

  for (i in unique(Data$curveID)) {

    d_1 <- subset(Data, curveID == i)    # Grab only current curveID
    #    d_1 <-subset(d_1, Photo > 0)         # Drop all values with Y < 0

    if (is.na(d_1$curve_type)[1]) {     # If no curve should be fit, all NAs
      r_tref = r_tref_SE = e = e_SE = eh = eh_SE = topt = topt_SE = r_max = breadth = AIC = r_sq = NA

    } else if(d_1$curve_type[1] == "Schoolfield") { # If curve type is Schoolfield, fit one

      # Fit Schoolfield function using nls_multstart
      fit <- nls_multstart(Photo ~ pawar_2018(temp = Tleaf, r_tref, e, eh, topt, tref),
                           data = d_1,
                           iter = 1000,
                           start_lower = c(r_tref = 0, e = 0, eh = 0.2, topt = 0),
                           start_upper = c(r_tref = 20, e = 2, eh = 5, topt = 50),
                           supp_errors = 'Y',
                           na.action = na.omit,
                           lower = c(r_tref = -10, e = 0, eh = 0, topt = 0))

      if(typeof(fit) != "NULL") { # If the fit was successful, extract parameters estimates and SE
        r_tref <- coef(fit)["r_tref"]
        r_tref_SE <- summary(fit)$coefficients[,'Std. Error']['r_tref']
        e <- coef(fit)["e"]
        e_SE <- summary(fit)$coefficients[,'Std. Error']['e']
        eh <- coef(fit)["eh"]
        eh_SE <- summary(fit)$coefficients[,'Std. Error']['eh']
        topt <- coef(fit)["topt"]
        topt_SE <- summary(fit)$coefficients[,'Std. Error']['topt']
        AIC <- AIC(fit)
        r_sq = 1-sum(resid(fit)^2)/sum((d_1$Photo-mean(d_1$Photo))^2)

        r_max = get_rmax_2(fit)
        breadth = get_breadth_2(fit, level = 0.8)

        # if (T_opt < min(d_1$Tleaf)) {
        #   J_ref = J_ref_SE = E = E_SE = E_D = E_D_SE = T_opt = T_opt_SE = AIC = r_sq = NA
        #   d_1$failure_status = "Poor fit: T_opt out of bounds"
        # }

      } else { # Otherwise, likely a convergence failure. All set to NA
        #       lnB0 = lnB0_SE = E = E_SE = E_D = E_D_SE = T_h = T_h_SE = AIC = r_sq = NA
        r_tref = r_tref_SE = e = e_SE = eh = eh_SE = topt = topt_SE = r_max = breadth = AIC = r_sq = NA
        d_1$failure_status = "Convergence failure"
      }

    } else if(d_1$curve_type[1] == "Arrhenius") { # If curve type is Arrhenius, fit one

      # Fit Arrhenius function using nls_multstart
      fit <- nls_multstart(Photo ~ arrhenius(temp = Tleaf, r_tref, e, tref),
                           data = d_1,
                           iter = 1000,
                           start_lower = c(r_tref = 0, e = 0),
                           start_upper = c(r_tref = 20, e = 4),
                           supp_errors = 'Y',
                           na.action = na.omit,
                           lower = c(r_tref = -10, e = 0))


      if(typeof(fit) != "NULL") { # If the fit was successful, extract parameters estimates and SE
        r_tref <- coef(fit)["r_tref"]
        r_tref_SE <- summary(fit)$coefficients[,'Std. Error']['r_tref']
        e <- coef(fit)["e"]
        e_SE <- summary(fit)$coefficients[,'Std. Error']['e']
        eh = eh_SE = topt = topt_SE = r_max = breadth = NA
        AIC <- AIC(fit)
        r_sq = 1-sum(resid(fit)^2)/sum((d_1$Photo-mean(d_1$Photo))^2)

      } else { # Otherwise, likely a convergence failure. All set to NA
        #        lnB0 = lnB0_SE = E = E_SE = E_D = E_D_SE = T_h = T_h_SE = AIC = r_sq = NA
        r_tref = r_tref_SE = e = e_SE = eh = eh_SE = topt = topt_SE =  r_max = breadth = AIC = r_sq = NA
        d_1$failure_status = "Convergence failure"
      }
    }


    ##########################
    # Plot curves if desired #
    ##########################
    if (PLOT==TRUE) {

      # Pick a set of temperature values to compute model over
      tmp_temps <- seq(min(floor(d_1$Tleaf)), ceiling(max(d_1$Tleaf)), length = 200)

      # Compute model predictions depending on model type and select curve color
      if(is.na(d_1$curve_type[1])) {
        tmp_model = mean(d_1$Photo)
        color_model = "white"
      } else if(d_1$curve_type[1] == "Schoolfield") {
        tmp_model <- pawar_2018(tmp_temps, r_tref, e, eh, topt, tref)
        color_model = "red"
      } else if(d_1$curve_type[1] == "Arrhenius") {
        tmp_model <- arrhenius(tmp_temps, r_tref, e, tref)
        color_model = "blue"
      }

      # Assemble model predicitons and raw data to plot
      ModelToPlotS <- data.frame(Temperature = tmp_temps, TraitValue = tmp_model)
      DataToPlot <- na.omit(data.frame(Temperature = d_1$Tleaf, TraitValue = d_1$Photo))

      # Build plots - normal
      p <- ggplot() +
        geom_point(data = DataToPlot, aes(x = Temperature, y = TraitValue),
                   shape = 21, size = 2.75, col = "black", fill = "white") +
        geom_line(data = ModelToPlotS,aes(x = Temperature, y = TraitValue), colour = color_model) +
        ggtitle(paste("curveID=", i, ", ", d_1$Taxon[1], sep="")) +
        xlab(expression(paste("Leaf temperature (", degree, C, ")"))) +
        ylab(expression("Assimilation rate (" * mu ~ "mol" ~m^-2 ~s^-1 * ")")) +
        theme_bw(base_size=12) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) +
        annotate("text", size = 3,
                 label=paste("R^2 = ", sprintf("%.2f", r_sq),"\nE =", format(e, digits = 3),"\nAIC =",format(AIC,digits=3), "\nCurve type = ", d_1$curve_type[1]),
                 x = min(DataToPlot[, "Temperature"]),
                 y = mean(DataToPlot[, "TraitValue"]),
                 hjust=0,
                 fontface = 3)

      # Build plots - Arrhenius plot
      p3 <- ggplot() +
        #        geom_point(data = DataToPlot, aes(x = 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
        #                   shape = 21, size = 2.75, col = "black", fill = "white") +
        #        geom_line(data = ModelToPlotS, aes(x = 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
        #                  colour = color_model) +
        geom_point(data = DataToPlot, aes(x = 1/(0.00008617*(tref + 273.15)) - 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
                   shape = 21, size = 2.75, col = "black", fill = "white") +
        geom_line(data = ModelToPlotS, aes(x = 1/(0.00008617*(tref + 273.15)) - 1/(0.00008617*(Temperature + 273.15)), y = TraitValue),
                  colour = color_model) +
        ggtitle(paste("curveID=", i, ", ", d_1$Taxon[1], sep="")) +
        #        xlab(expression(paste('Leaf temperature ', '1/',italic('kT'), ' (',  eV^{-1}, ')'))) +
        xlab(expression(paste('Leaf temperature ', '[ 1/',italic('kT_ref'), ' - 1/',italic('kT'), '] (',  eV^{-1}, ')'))) +
        ylab(expression("Assimilation rate (" * mu ~ "mol" ~m^-2 ~s^-1 * ")")) +
        #        scale_x_continuous(sec.axis = sec_axis(trans = ~ (1/(.*0.00008617))-273.15 , name = expression(paste("Leaf temperature (", degree, C, ")")))) +
        scale_x_continuous(sec.axis = sec_axis(trans = ~ (tref+273.15)/(1 - 0.00008617*(tref+273.15)*.) - 273.15, name = expression(paste("Leaf temperature (", degree, C, ")")))) +
        scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3),
                           labels = trans_format("log", math_format(e^.x))) +
        theme_bw(base_size=12) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"))

      grid.arrange(p,p3,nrow=1) # Arrange plots next to each other
    } # End plotting

    # Increment progress bar
    pb$tick()

    # Build dataframe with results
    results = rbind(results, data.frame(d_1[1,] %>% select(-Tleaf, -Photo), r_tref, r_tref_SE, e, e_SE,
                                        eh, eh_SE, topt, topt_SE, r_max, breadth, AIC, r_sq))
  }

  # Close output file
  if(PLOT==TRUE)
    grDevices::dev.off()

  return(results)
}

# 
# 
# # Returns Arrhenius function for the given parameters
# arrhenius <- function(temp, r_tref, e, tref = 25) {
#   # temp   : Temperature values to evaluate function at (C)
#   # r_tref  : Rate at T_ref (units depend on rate)
#   # e      : Activation energy (eV; E > 0)
#   # tref  : Reference temperature for normalization (C)
# 
#   k <- 8.62e-5           # Boltzmann's constant (eV K-1)
#   temp = temp + 273.15   # Convert to Kelvin
#   tref = tref + 273.15 # Convert to Kelvin
# 
#   # Evaluate and return
#   return( r_tref * exp( e * (1/(k*tref) - 1/(k*temp)) ) )
# }
# 


# In the given dataset, assigns each curveID to be fit with Arrhenius, Schoolfield, or neither
assign_curve_type <- function(data) {

  filter_criteria = c() # Holder

  #new_dat = c() ###
  # Find metrics for filtering for each curve
  for (i in unique(data$curveID)) {

    cur_dat = subset(data, curveID == i) # Grab the current data set
    Tpk = cur_dat$Tleaf[which( cur_dat$Photo == max(cur_dat$Photo) )[1]] # Find peak
    #cur_dat = subset(cur_dat, Tleaf <= Tpk) ###

    nbefore = sum( cur_dat$Tleaf < Tpk ) # Find number of points before peak
    nafter = sum( cur_dat$Tleaf > Tpk )  # Find number of points before peak
    Tbefore = Tpk-min(cur_dat$Tleaf)     # Find temperature span before peak

    filter_criteria = rbind(filter_criteria, data.frame(curveID = i, nbefore, nafter, Tbefore))

    #new_dat = rbind(new_dat, cur_dat) ###
  }

  # Merge filtering metrics with dataset for ease of filtering
  data = merge(data, filter_criteria, by="curveID")

  data$curve_type = NA
  data$failure_status = NA

  # Reject curves if the are decreasing only, have too little T span, or too few points
  data$failure_status = ifelse(data$nbefore < 2, "Rejected: fewer than two points before peak", data$failure_status)
  data$failure_status = ifelse(data$Tbefore < 4, "Rejected: ascending region too small", data$failure_status)
  data$failure_status = ifelse(data$nbefore+data$nafter < 4, "Rejected: fewer than 5 data points", data$failure_status)

  # Assign curve type based on whether or not they have an identifiable peak, or NA if failed above
  data$curve_type = ifelse(data$nbefore > 0 & data$nafter > 0, "Schoolfield", data$curve_type)
  data$curve_type = ifelse(data$nafter < 1, "Arrhenius", data$curve_type)
  data$curve_type = ifelse(!is.na(data$failure_status), NA, data$curve_type)

  data = data %>% select(-nbefore, -nafter, -Tbefore)  # Remove junk
  return(data)
}

get_rmax_2 = function (model)
{
  x <- model$m$getEnv()
  formula <- stats::as.formula(model$m$formula())
  param_ind <- all.vars(formula[[3]])[!all.vars(formula[[3]]) %in%
                                        names(model$m$getPars())]
  vals <- x[[param_ind[1]]]
  vals2 <- x[[param_ind[2]]]
  newdata <- data.frame(x = seq(min(vals), max(vals), by = 0.001),
                        y = unique(vals2),
                        stringsAsFactors = FALSE)
  names(newdata) <- param_ind
  newdata$preds <- stats::predict(model, newdata)
  rmax = newdata[newdata$preds == max(newdata$preds), "preds"]
  return(mean(rmax))
}

get_breadth_2 = function (model, level = 0.8)
{
  x <- model$m$getEnv()
  formula <- stats::as.formula(model$m$formula())
  param_ind <- all.vars(formula[[3]])[!all.vars(formula[[3]]) %in%
                                        names(model$m$getPars())]
  vals <- 0:50 #x[[param_ind[1]]]
  vals2 <- x[[param_ind[2]]]
  newdata_extrap <- data.frame(x = seq(min(vals), max(vals), by = 0.001),
                               y = unique(vals2),
                               stringsAsFactors = FALSE)
  newdata <- data.frame(x = seq(min(vals), max(vals), by = 1),
                        y = unique(vals2),
                        stringsAsFactors = FALSE)
  names(newdata) <- param_ind
  names(newdata_extrap) <- param_ind
  newdata$preds <- stats::predict(model, newdata = newdata)
  newdata_extrap$preds <- stats::predict(model, newdata = newdata_extrap)
  newdata_extrap <- newdata_extrap[!is.nan(newdata_extrap$preds),
  ]
  topt = newdata[newdata$preds == max(newdata$preds, na.rm = TRUE),
                 param_ind[1]]
  rmax = newdata[newdata$preds == max(newdata$preds), "preds"]
  newdata_extrap_low <- newdata_extrap[newdata_extrap[, param_ind[1]] <=
                                         topt, ]
  newdata_extrap_high <- newdata_extrap[newdata_extrap[, param_ind[1]] >=
                                          topt, ]
  newdata_extrap_low <- newdata_extrap_low[newdata_extrap_low$preds >=
                                             (level * rmax), ]
  newdata_extrap_high <- newdata_extrap_high[newdata_extrap_high$preds >=
                                               (level * rmax), ]
  low_val <- suppressWarnings(min(newdata_extrap_low[, param_ind[1]]))
  high_val <- suppressWarnings(max(newdata_extrap_high[, param_ind[1]]))
  return(high_val - low_val)
}

