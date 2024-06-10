make_main_fig_3 = function() {#function(resp_data, resp.nlme, photo_data, photo.nlme, bm_data, lma.df) {
  
  
  # Get Resp parameters
  resp.coefs = coef(resp.nlme)
  resp.coefs$curveID = rownames(resp.coefs)
  resp.coefs = resp_data %>% select(curveID, taxon, treatment) %>% unique() %>% merge(resp.coefs, by="curveID")
  
  resp.coefs$J_ref = resp.coefs$E = 0
  
  whichBOOF = resp.coefs$taxon == "BOOF"
  whichRASA = resp.coefs$taxon == "RASA"
  whichHOVU = resp.coefs$taxon == "HOVU"
  
  resp.coefs$J_ref[whichBOOF] = resp.coefs$J_ref.taxonBOOF[whichBOOF]
  resp.coefs$J_ref[whichRASA] = resp.coefs$J_ref.taxonRASA[whichRASA]
  resp.coefs$J_ref[whichHOVU] = resp.coefs$J_ref.taxonHOVU[whichHOVU]
  
  resp.coefs$E[whichBOOF] = resp.coefs$E.taxonBOOF[whichBOOF]
  resp.coefs$E[whichRASA] = resp.coefs$E.taxonRASA[whichRASA]
  resp.coefs$E[whichHOVU] = resp.coefs$E.taxonHOVU[whichHOVU]
  
  resp.coefs$J_ref = resp.coefs$J_ref + resp.coefs$`J_ref.(Intercept)`
  resp.coefs$E = resp.coefs$E + resp.coefs$`E.(Intercept)`
  
  resp.coefs = resp.coefs %>% 
    select(taxon, treatment, R10 = J_ref, ER = E) %>% 
    gather( key = "variable", value = "value", c("R10", "ER"))
  
  # Get Photo Parameters
  photo.coefs = coef(photo.nlme)
  photo.coefs$curveID = rownames(photo.coefs)
  photo.coefs = photo_data %>% select(curveID, taxon, treatment) %>% unique() %>% merge(photo.coefs, by="curveID")
  
  photo.coefs$r_tref = photo.coefs$e = photo.coefs$eh = photo.coefs$topt = 0
  
  whichBOOF = photo.coefs$taxon == "BOOF"
  whichRASA = photo.coefs$taxon == "RASA"
  whichHOVU = photo.coefs$taxon == "HOVU"
  
  photo.coefs$r_tref[whichBOOF] = photo.coefs$r_tref.taxonBOOF[whichBOOF]
  photo.coefs$r_tref[whichRASA] = photo.coefs$r_tref.taxonRASA[whichRASA]
  photo.coefs$r_tref[whichHOVU] = photo.coefs$r_tref.taxonHOVU[whichHOVU]
  
  photo.coefs$e[whichBOOF] = photo.coefs$e.taxonBOOF[whichBOOF]
  photo.coefs$e[whichRASA] = photo.coefs$e.taxonRASA[whichRASA]
  photo.coefs$e[whichHOVU] = photo.coefs$e.taxonHOVU[whichHOVU]
  
  photo.coefs$eh[whichBOOF] = photo.coefs$eh.taxonBOOF[whichBOOF]
  photo.coefs$eh[whichRASA] = photo.coefs$eh.taxonRASA[whichRASA]
  photo.coefs$eh[whichHOVU] = photo.coefs$eh.taxonHOVU[whichHOVU]
  
  photo.coefs$topt[whichBOOF] = photo.coefs$topt.taxonBOOF[whichBOOF]
  photo.coefs$topt[whichRASA] = photo.coefs$topt.taxonRASA[whichRASA]
  photo.coefs$topt[whichHOVU] = photo.coefs$topt.taxonHOVU[whichHOVU]
  
  
  photo.coefs$r_tref = photo.coefs$r_tref + photo.coefs$`r_tref.(Intercept)`
  photo.coefs$e = photo.coefs$e + photo.coefs$`e.(Intercept)`
  photo.coefs$eh = photo.coefs$eh + photo.coefs$`eh.(Intercept)`
  photo.coefs$topt = photo.coefs$topt + photo.coefs$`topt.(Intercept)`
  
  photo.coefs = photo.coefs %>% 
    select(taxon, treatment, A10 = r_tref, EA = e, ED = eh, Topt = topt) %>% 
    gather( key = "variable", value = "value", c("A10", "EA", "ED", "Topt"))
  
  # Get BM parameters
  bm.coefs = bm_data %>% 
    mutate(variable = "mu") %>% 
    select(taxon = species, treatment = temperature_C, variable, value = prop_leafy)
  
  # Get LMA parameters
  lma.coefs = lma.df %>% 
    mutate(taxon = toupper(species), variable = "lma") %>%
    select(taxon, treatment, variable, value = lma_g_m2)
  
  # We could add a shitload of other stuff here - all the physio data,
  # which is like 13 more variables. Not sure all of this is interesting, though
  # Maybe I should choose the ones I think will be the most interesting
  # Like gs and iWUE might be nice for a "water use" category
  # Vcmax25 and Amax might be nice
  # EA and Topt are obvious
  # LMA and Mu are nice
  # Is A10 interesting, or R10? I guess kind of;
  # Agrowth and Rgrowth might also be interesting
  # That makes 12 parameters, which is a nice rectangle
  
  physio.coefs = physio_data %>% 
    select(taxon, treatment, iWUEgrowth, gsgrowth, Vcmax25, Agrowth, Rdgrowth) %>% 
    gather( key = "variable", value = "value", c("iWUEgrowth", "gsgrowth", "Vcmax25", "Agrowth", "Rdgrowth"))
  
  
  
  
  # Assemble dataframe
  plot.data = bind_rows(photo.coefs, resp.coefs, bm.coefs, lma.coefs, physio.coefs)
  #plot.data = subset(plot.data, variable != "ED")
  plot.data = subset(plot.data, variable != "Agrowth")
  
  #plot.data$variable = factor(plot.data$variable, levels = c("A10", "Agrowth", "EA","Topt","Vcmax25", "R10",
  #                                                           "iWUEgrowth", "Rdgrowth", "ER", "gsgrowth", "lma","mu"))
  plot.data$variable = factor(plot.data$variable, levels = c("A10", "EA","Topt","ED", "Vcmax25", "R10",
                                                             "iWUEgrowth", "Rdgrowth", "ER", "gsgrowth", "lma","mu"))
  
  
  # Now for each taxon and variable, we need to decide linear or quadratic
  plot.data = plot.data %>% group_by(taxon, variable) %>% mutate(curveID = group_indices())
  plot.data$type = NA
  for (i in unique(plot.data$curveID)) {
    subdat = subset(plot.data, curveID == i)
    subdat$treatment2 = subdat$treatment^2
    
    z1 = lm(value ~ treatment, subdat)
    z2 = lm(value ~ treatment + treatment2, subdat)
    
    if(AIC(z2) < AIC(z1)-2){
      plot.data$type[plot.data$curveID == i] = "quadratic"
    } else {
      plot.data$type[plot.data$curveID == i] = "linear"
    }
    
  }
  
  plot.data = plot.data %>% 
    subset(variable != "iWUEgrowth") %>% 
    subset(variable != "gsgrowth") %>% 
    subset(variable != "Rdgrowth")
  
 plot.data$variable =  factor(plot.data$variable, levels = c("A10", "EA", "Topt", "ED",  "R10", "ER", "Vcmax25", "lma" ,"mu"))
  
  anno = data.frame(label = c("a.","b.","d.","c.","e.","f.","i.","h.","g."),
                    x = rep(5,9),
                    y = rep(0,9),
                    variable = plot.data$variable %>% unique())
  

  
  
  # Make plots
  p = ggplot(data = plot.data, aes(x = treatment, y = value, color = taxon)) +
    my_theme +
    theme_classic() +
    scale_color_brewer(palette = "Dark2",labels = c(expression(italic("B. officinalis")), expression(italic("H. vulgare")), expression(italic("R. sativus")))) +
    #geom_smooth() +
    stat_smooth(data = subset(plot.data, type == "linear"), method = "lm", formula = y ~ x, size = 0.5) +
    stat_smooth(data = subset(plot.data, type == "quadratic"), method = "lm", formula = y ~ x + I(x^2), size = 0.5) +
    geom_jitter(size = 0.5,width = 1,height = 0) +
    theme(legend.position = "bottom") +
    theme_transparent +
    #theme(legend.text = element_markdown()) +
    theme(legend.title = element_blank()) +
    xlab("Treatment temperature (ºC)") +
    ylab(NULL) +
    #scale_color_discrete(labels = c("B.", "H.", "R.")) +
    theme(strip.background = element_blank(),
          strip.placement = "outside") +
    geom_text(data = anno, 
              #aes(xlabel=label),
              mapping = aes(x = -Inf, y = -Inf, label = label),
              hjust   = -0.2,
              #vjust   = -8.5,
              vjust = -9.3,
              color = "black") +
    #theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    facet_wrap(~variable, 
               scales="free", 
               ncol = 3,
               strip.position = "left",
               labeller = as_labeller(c(A10 = "A[10]~(μmol~m^-2~s^-1)",
                                        EA = "E[A]^A~(eV)",
                                        Topt = "T[opt]^A~(phantom()^o*C)",
                                        ED = "E[D]^A~(eV)",
                                        #Agrowth = "A[growth]~(μmol/m^2*s)", 
                                        R10 = "R[10]~(μmol~m^-2~s^-1)",
                                        ER = "E[A]^R~(eV)",
                                        lma = "ρ[A]~(g~m^-2)",
                                        mu = "μ~(dimensionless)", 
                                        Rdgrowth = "R[growth]~(μmol~m^-2~s^-1)", 
                                        Vcmax25 = "V[cmax,25]~(μmol~m^-2~s^-1)", 
                                        iWUEgrowth = "iWUE~(μmol~mol^-1)", 
                                        gsgrowth="g[sw]~(mol~m^-2~s^-1)"), 
                                      default = label_parsed)) #+

  p
  
  # Save plots to file
  svg("figures/main_fig_3.svg", width = 5.5, height = 5)
  grid.arrange(p)
  dev.off()
  
}
