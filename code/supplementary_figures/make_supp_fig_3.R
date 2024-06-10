make_supp_fig_3 = function() {
  
  
  #library(ggpattern)
  
  
  
  # Next, we need to evaluate the model at the tissue T for each curve
  # The easiest way to do this is probably to use predict()
  photo_data %>% select(curveID, taxon, treatment) %>% unique() %>% merge(T_tissue_data) %>% 
    rename(Tleaf = T_tissue) %>% mutate(Photo = predict(photo.nlme,.)) -> photo_composite
  
  # Now, to fit an Arrhenius function to each composite curve
  photo.nls.boof = nls(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                       data = subset(photo_composite, taxon == "BOOF"),
                       start = c(r_tref = 5, e = 0.5, eh = 1, topt = 21))
  
  photo.nls.hovu = nls(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                       data = subset(photo_composite, taxon == "HOVU"),
                       start = c(r_tref = 5, e = 0.5, eh = 1, topt = 21))
  
  # nls.rasa = nls(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt = 20, tref = 10),
  #                data = subset(photo_composite, taxon == "RASA"),
  #                start = c(r_tref = 10, e = 1, eh = 2.29))
  
  photo.nls.rasa = nls_multstart(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                                 data = subset(photo_composite, taxon == "RASA"),
                                 iter = 1000,
                                 start_lower = c(r_tref = 5, e = 0.1, eh = 1, topt = 15),
                                 start_upper = c(r_tref = 50, e = 5, eh = 10, topt = 50),
                                 lower = c(r_tref = 0, e = 0, eh = 0, topt = 0),
                                 supp_errors = "Y")
  
  strip_data = bind_rows(
    
    data.frame(low = intervals(photo.nlme, which = "fixed")$fixed[4:6,1], 
               E = intervals(photo.nlme, which = "fixed")$fixed[4:6,2], 
               high = intervals(photo.nlme, which = "fixed")$fixed[4:6,3],
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "An"), # Boof E
    
    data.frame(low = c(confint2(photo.nls.boof)[2,1],confint2(photo.nls.hovu)[2,1],confint2(photo.nls.rasa)[2,1]),
               E = c(coef(photo.nls.boof)["e"],coef(photo.nls.hovu)["e"],coef(photo.nls.rasa)["e"]),
               high = c(confint2(photo.nls.boof)[2,2],confint2(photo.nls.hovu)[2,2],confint2(photo.nls.rasa)[2,2]),
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "An'"), # Boof E
    
    data.frame(low = intervals(bm.nlme)$fixed[4:6,1], 
               E = intervals(bm.nlme)$fixed[4:6,2], 
               high = intervals(bm.nlme)$fixed[4:6,3],
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "G")
  )
  
  
  
  #strip_data$rate = as.factor(strip_data$rate)
  #strip_data$rate = factor(strip_data$rate, c("Growth", "Photo_net", "Photo_net_accl"))
  
  
  strip_data_E_D = bind_rows(
    
    
    data.frame(low = intervals(photo.nlme, which = "fixed")$fixed[7:9,1], 
               E_D = intervals(photo.nlme, which = "fixed")$fixed[7:9,2], 
               high = intervals(photo.nlme, which = "fixed")$fixed[7:9,3],
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "An"), # Boof E
    
    data.frame(low = c(confint2(photo.nls.boof)[3,1],confint2(photo.nls.hovu)[3,1],confint2(photo.nls.rasa)[3,1]),
               E_D = c(coef(photo.nls.boof)["eh"],coef(photo.nls.hovu)["eh"],coef(photo.nls.rasa)["eh"]),
               high = c(confint2(photo.nls.boof)[3,2],confint2(photo.nls.hovu)[3,2],confint2(photo.nls.rasa)[3,2]),
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "An'"), # Boof E
    
    data.frame(low = intervals(bm.nlme)$fixed[7:9,1], 
               E_D = intervals(bm.nlme)$fixed[7:9,2], 
               high = intervals(bm.nlme)$fixed[7:9,3],
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "G")
  )
  
  
  strip_data_T_opt = bind_rows(
    
    
    data.frame(low = intervals(photo.nlme, which = "fixed")$fixed[10:12,1], 
               T_opt = intervals(photo.nlme, which = "fixed")$fixed[10:12,2], 
               high = intervals(photo.nlme, which = "fixed")$fixed[10:12,3],
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "An"), # Boof E
    
    data.frame(low = c(confint2(photo.nls.boof)[4,1],confint2(photo.nls.hovu)[4,1],confint2(photo.nls.rasa)[4,1]),
               T_opt = c(coef(photo.nls.boof)["topt"],coef(photo.nls.hovu)["topt"],coef(photo.nls.rasa)["topt"]),
               high = c(confint2(photo.nls.boof)[4,2],confint2(photo.nls.hovu)[4,2],confint2(photo.nls.rasa)[4,2]),
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "An'"), # Boof E
    
    data.frame(low = intervals(bm.nlme)$fixed[10:12,1], 
               T_opt = intervals(bm.nlme)$fixed[10:12,2], 
               high = intervals(bm.nlme)$fixed[10:12,3],
               taxon = c("BOOF", "HOVU", "RASA"),
               model = "G")
  )
  
  
  
  
  
  
  
  
  
  
  # For each of three acclimation models, predict Vcmax, Jmax, etc.
  ## Kattge and Knorr 2007
  ## Kumarathunge 2018
  ## Smith ??
  # Generate An'(T) predictions with FvCB model
  # Fit SS to An' and extract T-dependence parameters
  # Compare to G and An' empirical
  # Make figure and possibly table
  
  # Predict photosynthetic capacity with Kumarathunge 2018
  ## For each species, we need Thome, Tgrowth (Vector), and Vcmax25
  ## We will use Topt of growth as Thome, and Vcmax25 from Thome I guess
  
  
  
  Thome_BOOF = coef(bm.nlme)$T_opt.taxonBOOF[1]
  Vc25_BOOF = subset(physio_data, taxon == "BOOF")$Vcmax25 %>% mean()
  Thome_HOVU = coef(bm.nlme)$T_opt.taxonHOVU[1]
  Vc25_HOVU = subset(physio_data, taxon == "HOVU")$Vcmax25 %>% mean()
  Thome_RASA = coef(bm.nlme)$T_opt.taxonRASA[1]
  Vc25_RASA = subset(physio_data, taxon == "RASA")$Vcmax25 %>% mean()
  
  Tgrowth = seq(2.5,35, 0.5)
  
  dat_K2018_BOOF = K2018(Vc25_BOOF, Thome_BOOF, Tgrowth) %>% mutate(taxon = "BOOF", model = "K18")
  dat_K2018_HOVU = K2018(Vc25_HOVU, Thome_HOVU, Tgrowth) %>% mutate(taxon = "HOVU", model = "K18")
  dat_K2018_RASA = K2018(Vc25_RASA, Thome_RASA, Tgrowth) %>% mutate(taxon = "RASA", model = "K18")
  
  # Let's try the Kattge and Knorr model
  dat_K2007_BOOF = K2007(Vc25_BOOF, Tgrowth) %>% mutate(taxon = "BOOF", model = "K07")
  dat_K2007_HOVU = K2007(Vc25_HOVU, Tgrowth) %>% mutate(taxon = "HOVU", model = "K07")
  dat_K2007_RASA = K2007(Vc25_RASA, Tgrowth) %>% mutate(taxon = "RASA", model = "K07")
  
  # Finally, the least-cost / coordination hypothesis
  dat_S2018_BOOF = S2018(Tgrowth) %>% mutate(taxon = "BOOF", model = "S19")
  dat_S2018_HOVU = S2018(Tgrowth) %>% mutate(taxon = "HOVU", model = "S19")
  dat_S2018_RASA = S2018(Tgrowth) %>% mutate(taxon = "RASA", model = "S19")

  
  # Once we have the data, we need to collate it
  data_all = bind_rows(dat_K2018_BOOF, dat_K2018_HOVU, dat_K2018_RASA,
                       dat_K2007_BOOF, dat_K2007_HOVU, dat_K2007_RASA,
                       dat_S2018_BOOF, dat_S2018_HOVU, dat_S2018_RASA)
  
  data_all = data_all %>% group_by(taxon, model) %>% mutate(curveID = group_indices())
  
  # Fit SS
  fits = fit_curves_all(data_all, x = "Tleaf", y = "ALEAF", tref = 10, PLOT = T)
  fits_small = fits %>% select(taxon, model, r_tref, e, eh, topt)
  
  # Now we gotta assemble the barchart comparison
  #fits_small
  fits_long = gather(fits_small, key = "parameter", value = "value", c("e","eh","topt"))
  
  # Get strip data from figure 3 code
  e_dat = strip_data %>% select(taxon, model, E) %>% rename(value = E) %>% mutate(parameter = "e")
  eh_dat = strip_data_E_D %>% select(taxon, model, E_D) %>% rename(value = E_D) %>% mutate(parameter = "eh")
  topt_dat = strip_data_T_opt %>% select(taxon, model, T_opt) %>% rename(value = T_opt) %>% mutate(parameter = "topt")
  
  dat_all = bind_rows(e_dat, eh_dat, topt_dat) 
  dat_all$taxon = gsub("B. officinalis", "BOOF", dat_all$taxon)
  dat_all$taxon = gsub("H. vulgare", "HOVU", dat_all$taxon)
  dat_all$taxon = gsub("R. sativus", "RASA", dat_all$taxon)
  # dat_all$model = gsub("Growth", "G", dat_all$model)
  # dat_all$model = gsub("An", "bquote(A[n])", dat_all$model)
  # dat_all$model = gsub("An_accl", "An'", dat_all$model)
  
  
  
  fits_long = fits_long %>% bind_rows(dat_all) %>% mutate(newvar = paste(taxon, parameter))
  
  # I think we are back to a 3x3 panel grid; each panel has to be a species x parameter comparison
  labels = c(`BOOF e` = "EA (eV)", `BOOF eh` = "ED (eV)", `BOOF topt` = "Topt (C)",
             `HOVU e` = "EA (eV)", `HOVU eh` = "ED (eV)", `HOVU topt` = "Topt (C)",
             `RASA e` = "EA (eV)", `RASA eh` = "ED (eV)", `RASA topt` = "Topt (C)")
  
  color_vals = c("#00BA38", "#01611e", "#619CFF","purple1", "purple3", "purple4")
  #fits_long$model = factor("Growth", "Photo_net", )
  # 
  # 
  # ggplot(data = fits_long, aes(x = model, y = value, fill = model)) +
  #   geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = "white") +
  #   #my_theme +
  #   facet_wrap(~ newvar, scales="free",
  #              strip.position = "left",  #,
  #              labeller = as_labeller(labels)) +
  #              #labeller = as_labeller(c(e = "EA (eV)", eh = "ED (eV)", topt = "Topt (C)", BOOF = "element_blank()", HOVU = "", RASA = "" ))) +
  #   my_theme +
  #   theme_classic() +
  #   scale_fill_manual(values = color_vals) +
  #   theme(strip.background = element_blank(),
  #         strip.placement = "outside",
  #         axis.title.y = element_blank(),
  #         axis.title.x = element_blank(),
  #         legend.position = "none")

  # stuff = geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = "white") +
  #   #my_theme +
  #   # facet_wrap(~ newvar, scales="free",
  #   #            strip.position = "left",  #,
  #   #            labeller = as_labeller(labels)) +
  #   #labeller = as_labeller(c(e = "EA (eV)", eh = "ED (eV)", topt = "Topt (C)", BOOF = "element_blank()", HOVU = "", RASA = "" ))) +
  #   my_theme +
  #   theme_classic() +
  #   scale_fill_manual(values = color_vals) +
  #   theme(strip.background = element_blank(),
  #         strip.placement = "outside",
  #         axis.title.y = element_blank(),
  #         axis.title.x = element_blank(),
  #         legend.position = "none")
  # 
  
  a = fits_long %>% subset(taxon == "BOOF" & parameter == "e")
  a = strip_data %>% subset(taxon == "BOOF") %>% merge(a,., by=c("model", "taxon"), all=T)
  p1 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    geom_abline(slope = 0, intercept = 0.32, lty=2) +
    geom_abline(slope = 0, intercept = 0.65, lty=3) +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = expression(paste("a. ", italic("Borago officinalis")))) +
    ylab(bquote("E"[A]~"(eV)")) +
    ylim(0,2.0) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  a = fits_long %>% subset(taxon == "HOVU" & parameter == "e")
  a = strip_data %>% subset(taxon == "HOVU") %>% merge(a,., by=c("model", "taxon"), all=T)
  p2 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    geom_abline(slope = 0, intercept = 0.32, lty=2) +
    geom_abline(slope = 0, intercept = 0.65, lty=3) +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = expression(paste("b. ", italic("Hordeum vulgare")))) +
    ylab(bquote("E"[A]~"(eV)")) +
    ylim(0,2.0) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  a = fits_long %>% subset(taxon == "RASA" & parameter == "e")
  a = strip_data %>% subset(taxon == "RASA") %>% merge(a,., by=c("model", "taxon"), all=T)
  p3 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    geom_abline(slope = 0, intercept = 0.32, lty=2) +
    geom_abline(slope = 0, intercept = 0.65, lty=3) +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = expression(paste("c. ", italic("Raphanus sativus")))) +
    ylab(bquote("E"[A]~"(eV)")) +
    ylim(0,2.0) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  a = fits_long %>% subset(taxon == "BOOF" & parameter == "eh")
  a = strip_data_E_D %>% subset(taxon == "BOOF") %>% merge(a,., by=c("model", "taxon"), all=T)
  p4 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = "d.") +
    ylab(bquote("E"[D]~"(eV)")) +
    ylim(0,5.4) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  
  a = fits_long %>% subset(taxon == "HOVU" & parameter == "eh")
  a = strip_data_E_D %>% subset(taxon == "HOVU") %>% merge(a,., by=c("model", "taxon"), all=T)
  p5 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = "e.") +
    ylab(bquote("E"[D]~"(eV)")) +
    ylim(0,5.4) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  
  a = fits_long %>% subset(taxon == "RASA" & parameter == "eh")
  a = strip_data_E_D %>% subset(taxon == "RASA") %>% merge(a,., by=c("model", "taxon"), all=T)
  p6 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = "f.") +
    ylab(bquote("E"[D]~"(eV)")) +
    ylim(0,5.4) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  
  a = fits_long %>% subset(taxon == "BOOF" & parameter == "topt")
  a = strip_data_T_opt %>% subset(taxon == "BOOF") %>% merge(a,., by=c("model", "taxon"), all=T)
  p7 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = "g.") +
    ylab(bquote("T"[opt]~"(°C)")) +
    ylim(0,26.5) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  a = fits_long %>% subset(taxon == "HOVU" & parameter == "topt")
  a = strip_data_T_opt %>% subset(taxon == "HOVU") %>% merge(a,., by=c("model", "taxon"), all=T)
  p8 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = "h.") +
    ylab(bquote("T"[opt]~"(°C)")) +
    ylim(0,26.5) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  
  
  a = fits_long %>% subset(taxon == "RASA" & parameter == "topt")
  a = strip_data_T_opt %>% subset(taxon == "RASA") %>% merge(a,., by=c("model", "taxon"), all=T)
  p9 = ggplot(data = a, aes(x = model, y = value, fill = model)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = NA) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    my_theme +
    theme_classic() +
    scale_fill_manual(values = color_vals) +
    scale_x_discrete(labels = c(bquote(A[n]), bquote(A*"′"[n]), bquote(G[M]), "K07", "K19", "S19")) +
    labs(title = "i.") +
    ylab(bquote("T"[opt]~"(°C)")) +
    ylim(0,26.5) +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank(),
          legend.position = "none")
  # 
  # svg("figures/main_fig_3.svg", 7, 7)
  # plot_grid(p1,p4,p7,p2,p5,p8,p3,p6,p9, ncol=3)
  # dev.off()
  # 
  
  # Add code here to plot supplementary figure 5
  
  
  # For loop
  # Step through the fits_small data frame
  # For each row, generate SS model predictions on [2.5, 35]
  # Save to a big dataframe with model and taxon
  
  data_preds = data.frame()
  for(i in 1:dim(fits_small)[1]) {
    
    curdat = fits_small[i,]
    temps = seq(2.5,35,0.5)
    preddat = pawar_2018(temps, 
                         r_tref = curdat$r_tref, 
                         e = curdat$e, 
                         eh = curdat$eh,
                         topt = curdat$topt,
                         tref = 10)
    
    newdf = data.frame(taxon = curdat$taxon,
                       model = curdat$model,
                       Tleaf = temps,
                       ALEAF = preddat)
    
    data_preds = bind_rows(data_preds, newdf)
  }
  
  p_sup = ggplot(data_all, aes(x = Tleaf, y = ALEAF, color = taxon)) + 
    geom_point() +
    geom_line(data = data_preds, color = "black") +
    scale_color_brewer(palette = "Dark2") +
    my_theme +
    xlab("Growth temperaure (ºC)") +
    ylab(bquote("Assimilation rate (µmol "~m^"-2"~s^"-1"*")")) +
    #theme_classic() +
    facet_grid(rows = vars(taxon), cols = vars(model), scales="free",
               labeller = as_labeller(c(K07 = "K07",
                                        K18 = "K19",
                                        S19 = "S19",
                                        BOOF = "B. officinalis",
                                        HOVU = "H. vulgare",
                                        RASA = "R. sativus")))
  
  svg("figures/supp_fig_3.svg", width = 6, height = 5)
  grid.arrange(p_sup)
  dev.off()
  
}
