make_main_fig_2 = function() {#function(photo.nlme, photo_data, bm.nlme, bm_data, T_tissue_data) {
  
  cols = c("#00D0FF", "#282AB7", "#A200ED", "#FF0000", "#660000")
  
  photo_temp=photo_data %>% 
    group_by(curveID) %>% 
    mutate(Tmin = min(Tleaf), Tmax = max(Tleaf)) %>% 
    select(curveID, taxon, treatment, Tmin, Tmax) %>% 
    unique()
  
  photo_pred = c()
  for (i in unique(photo_temp$curveID)) {
    curdat = subset(photo_temp, curveID == i)
    df = data.frame(curveID = curdat$curveID, 
                    treatment = curdat$treatment, 
                    taxon = curdat$taxon, 
                    Tleaf = seq(curdat$Tmin, curdat$Tmax, 0.2))
    photo_pred = bind_rows(photo_pred, df)
  }
  
  # Now let's make predictions from the model
  photo_pred$Photo = predict(photo.nlme, photo_pred)
  
  # Next, we need to evaluate the model at the tissue T for each curve
  # The easiest way to do this is probably to use predict()
  photo_data %>% select(curveID, taxon, treatment) %>% unique() %>% merge(T_tissue_data) %>% 
    rename(Tleaf = T_tissue) %>% mutate(Photo = predict(photo.nlme,.)) -> photo_composite
  
  # Now, to fit an Arrhenius function to each composite curve
  nls.boof = nls(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                 data = subset(photo_composite, taxon == "BOOF"),
                 start = c(r_tref = 5, e = 0.5, eh = 1, topt = 21))
  
  nls.hovu = nls(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                 data = subset(photo_composite, taxon == "HOVU"),
                 start = c(r_tref = 5, e = 0.5, eh = 1, topt = 21))
  
  # nls.rasa = nls(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt = 20, tref = 10),
  #                data = subset(photo_composite, taxon == "RASA"),
  #                start = c(r_tref = 10, e = 1, eh = 2.29))
  
  nls.rasa = nls_multstart(Photo ~ pawar_2018(Tleaf, r_tref, e, eh, topt, tref = 10),
                           data = subset(photo_composite, taxon == "RASA"),
                           iter = 1000,
                           start_lower = c(r_tref = 5, e = 0.1, eh = 1, topt = 15),
                           start_upper = c(r_tref = 50, e = 5, eh = 10, topt = 50),
                           lower = c(r_tref = 0, e = 0, eh = 0, topt = 0),
                           supp_errors = "Y")
  
    trange = seq(3.5, 33.5, 0.5)
  
  comp.boof = data.frame(Tleaf = trange) %>% mutate(Photo = predict(nls.boof, .))
  comp.hovu = data.frame(Tleaf = trange) %>% mutate(Photo = predict(nls.hovu, .))
  comp.rasa = data.frame(Tleaf = trange) %>% mutate(Photo = predict(nls.rasa, .))
  
  #pdf("figures/fig10_boof.pdf", width = 3, height = 2.65)
  spec = "BOOF"
  a = subset(photo_data, taxon == spec)
  b = subset(photo_pred, taxon == spec)
  c = subset(photo_composite, taxon == spec)
  #d = subset(composite_fit_pred, taxon == spec)
  p1a = ggplot(data = a, aes(x = Tleaf, y = Photo, color = treatment)) +
    geom_point(size = 0.5) + 
    #scale_colour_gradient(low = "#0000ff80", high = "#ff000080", name = "Treatment (ºC)") +
    #scale_colour_gradient(limits = c(2.5,35),low = "blue", high = "red", name = "Treatment (ºC)") +
    scale_colour_gradientn(colours = cols, limits = c(2.5, 35), name = "Treatment (ºC)") +
    #scale_color_viridis_b(option = "magma") +
    geom_line(data = b, aes(x = Tleaf, y = Photo, group = curveID), alpha = 1) +
    my_theme +
    theme_classic() +
    theme_transparent +
    #geom_point(data = c, color = "black") +
    #geom_line(data = comp.boof, size = 1, lty = 1, color = "#00BA38") +
    #guides(color = guide_legend(title="Treatment (ºC)")) +
    #guides(color = guide_legend(override.aes = list(alpha=0.5) ) ) +
    #xlab("Leaf temperature (ºC)") +
    xlim(-2,42) +
    ylim(0,29) +
    theme(axis.title.x = element_blank()) + 
    #theme(legend.direction="horizontal") +
    #theme(legend.key.size = unit(0.35, 'cm')) +
    #theme(legend.title = element_text(size=9)) + #+
    ylab(" ") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(legend.position = "none") +
    #annotate("text",x=15,y=29,label = "italic(Borago~officinalis)", parse=T) +
    #xlim(-2.5,40) +
    #ylim(0,6.5) +
    #theme(legend.position = c(0.55, 0.1), legend.background = element_rect(fill=alpha('blue', 0))) +
    #labs(tag = "a") + 
    #labs(title = "a. Borago officinalis")
    labs(title = expression(paste("a. ", italic("Borago officinalis"))))
  #grid.arrange(p1)
  #dev.off()
  
  #pdf("figures/fig10_rasa.pdf", width = 4, height = 3.5)
  spec = "HOVU"
  a = subset(photo_data, taxon == spec)
  b = subset(photo_pred, taxon == spec)
  c = subset(photo_composite, taxon == spec)
  # d = subset(composite_fit_pred, taxon == spec)
  p2a = ggplot(data = a, aes(x = Tleaf, y = Photo, color = treatment)) +
    geom_point(size = 0.5) + 
    #scale_colour_gradient(low = "blue", high = "red") +
    scale_colour_gradientn(colours = cols, limits = c(2.5, 35), name = "Treatment (ºC)") +
    geom_line(data = b, aes(x = Tleaf, y = Photo, color = treatment, group = curveID), alpha = 1) +
    my_theme +
    theme_classic() +
    theme_transparent +
    #geom_point(data = c, size = 2, color = "black") +
    #geom_line(data = comp.rasa, size = 1, color = "black") +
    #xlab("Leaf temperature (ºC)") +
    xlim(-2,42) +
    ylim(0,29) +
    #theme(axis.title.x = element_blank()) + 
    theme(legend.direction="horizontal") +
    theme(legend.key.size = unit(0.35, 'cm')) +
    theme(legend.title = element_text(size=9)) +
    theme(legend.position = c(0.45, 0.85), legend.background = element_rect(fill=alpha('blue', 0))) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.title.x = element_blank()) +
    ylab(bquote("Assimilation rate (acute),"~italic(A[n])~"(µmol "~m^"-2"~s^"-1"*")")) +
    #annotate("text",x=15,y=29,label = "italic(Hordeum~vulgare)", parse=T) +
    #xlim(-2.5,40) +
    #ylim(0,6.5) +
    #theme(legend.position = "none") +
    labs(title = expression(paste("b. ", italic("Hordeum vulgare"))))
  #ggtitle("Photosynthesis TPCs for Raphanus sativus")
  #grid.arrange(p2)
  #dev.off()
  
  
  #pdf("figures/fig10_hovu.pdf", width = 4, height = 3.5)
  spec = "RASA"
  a = subset(photo_data, taxon == spec)
  b = subset(photo_pred, taxon == spec)
  c = subset(photo_composite, taxon == spec)
  # d = subset(composite_fit_pred, taxon == spec)
  p3a = ggplot(data = a, aes(x = Tleaf, y = Photo, color = treatment)) +
    geom_point(size = 0.5) + 
    #scale_colour_gradient(low = "blue", high = "red") +
    scale_colour_gradientn(colours = cols, limits = c(2.5, 35), name = "Treatment (ºC)") +
    
    #scale_colour_carto_c(name = "Treatment (ºC)",
    #                   type = "diverging", palette = "Sunset", direction = -1) +
    geom_line(data = b, aes(x = Tleaf, y = Photo, color = treatment, group = curveID)) +
    my_theme +
    
    theme_classic() +
    theme_transparent +
    xlim(-2,42) +
    ylim(0,29) +
    #geom_point(data = c, color = "black") +
    #geom_line(data = comp.hovu, size = 0.5, lty = 2, color = "#00BA38") +
    xlab("Leaf temperature (ºC)") +
    ylab(" ") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    #annotate("text",x=15,y=29,label = "italic(Raphanus~sativus)", parse=T) +
    #xlim(-2.5,40) +
    #ylim(0,6.5) +
    theme(legend.position = "none") +
    labs(title = expression(paste("c. ", italic("Raphanus sativus"))))
  #ggtitle("Photosynthesis TPCs for Hordeum vulgare")
  
  #grid.arrange(p3)
  #dev.off()
  
  
  
  
  #pdf("figures/fig10_boof.pdf", width = 3, height = 2.65)
  spec = "BOOF"
  a = subset(photo_data, taxon == spec)
  b = subset(photo_pred, taxon == spec)
  c = subset(photo_composite, taxon == spec)
  #d = subset(composite_fit_pred, taxon == spec)
  p1b = ggplot(data = a, aes(x = Tleaf, y = Photo, color = treatment)) +
    #geom_point() + 
    #scale_colour_gradient(low = "#0000ff80", high = "#ff000080", name = "Treatment (ºC)") +
    scale_colour_gradientn(colours = cols, limits = c(2.5, 35), name = "Treatment (ºC)") +
    geom_line(data = b, aes(x = Tleaf, y = Photo, group = curveID), alpha = 0.15) +
    my_theme +
    theme_classic() +

    theme_transparent +
    geom_point(size = 0.5,data = c, color = "black") +
    geom_line(data = comp.boof, color = "#00BA38") +
    xlim(-2,42) +
    ylim(0,29) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    #guides(color = guide_legend(title="Treatment (ºC)")) +
    #guides(color = guide_legend(override.aes = list(alpha=0.5) ) ) +
    #xlab("Leaf temperature (ºC)") +
    #ylab("Assimilation rate (µmol/m²s)") #+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
    theme(legend.position = "none") +
    labs(title = "d.")
    #xlim(-2.5,40) +
    #ylim(0,6.5) +
    #theme(legend.position = c(0.22, 0.65), legend.background = element_rect(fill=alpha('blue', 0))) 
  #grid.arrange(p1)
  #dev.off()
  
  #pdf("figures/fig10_rasa.pdf", width = 4, height = 3.5)
  spec = "HOVU"
  a = subset(photo_data, taxon == spec)
  b = subset(photo_pred, taxon == spec)
  c = subset(photo_composite, taxon == spec)
  # d = subset(composite_fit_pred, taxon == spec)
  p2b = ggplot(data = a, aes(x = Tleaf, y = Photo, color = treatment)) +
    #geom_point() + 
    #scale_colour_gradient(low = "#0000ff80", high = "#ff000080", name = "Treatment (ºC)") +
    scale_colour_gradientn(colours = cols, limits = c(2.5, 35), name = "Treatment (ºC)") +
    geom_line(data = b, aes(x = Tleaf, y = Photo, color = treatment, group = curveID), alpha = 0.15) +
    my_theme +
    theme_classic() +
    theme_transparent +
    geom_point(size = 0.5, data = c, color = "black") +
    geom_line(data = comp.hovu,  color = "#00BA38") +
    xlim(-2,42) +
    ylim(0,29) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    #xlab("Leaf temperature (ºC)") +
    ylab(bquote("Assimilation rate (acclimated),"~italic(A*"′"[n])~"(µmol "~m^"-2"~s^"-1"*")")) +
    theme(axis.title.x = element_blank()) +#,
          #axis.title.y = element_blank()) +
    #xlim(-2.5,40) +
    #ylim(0,6.5) +
    theme(legend.position = "none") +
    labs(title = "e.") 
  #ggtitle("Photosynthesis TPCs for Raphanus sativus")
  #grid.arrange(p2)
  #dev.off()
  
  
  #pdf("figures/fig10_hovu.pdf", width = 4, height = 3.5)
  spec = "RASA"
  a = subset(photo_data, taxon == spec)
  b = subset(photo_pred, taxon == spec)
  c = subset(photo_composite, taxon == spec)
  # d = subset(composite_fit_pred, taxon == spec)
  p3b = ggplot(data = a, aes(x = Tleaf, y = Photo, color = treatment)) +
    #geom_point() + 
    #scale_colour_gradient(low = "#0000ff80", high = "#ff000080", name = "Treatment (ºC)") +
    scale_colour_gradientn(colours = cols, limits = c(2.5, 35), name = "Treatment (ºC)") +
    #scale_colour_carto_c(name = "Treatment (ºC)",
    #                   type = "diverging", palette = "Sunset", direction = -1) +
    geom_line(data = b, aes(x = Tleaf, y = Photo, color = treatment, group = curveID), alpha = 0.15) +
    my_theme +
    theme_classic() +
    theme_transparent +
    geom_point(size = 0.5, data = c, color = "black") +
    geom_line(data = comp.rasa, color = "#00BA38") +
    xlim(-2,42) +
    ylim(0,29) +
    xlab("Leaf temperature (ºC)") +
    ylab(" ") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    #ylab("Assimilation rate (µmol/m²s)") +
    theme(#axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    #xlim(-2.5,40) +
    #ylim(0,6.5) +
    theme(legend.position = "none") +
    labs(title = "f.") 
  #ggtitle("Photosynthesis TPCs for Hordeum vulgare")
  
  #grid.arrange(p3)
  #dev.off()
  
    biomass = bm_data %>% rename(treatment = temperature_C, taxon = species)
    biomass_thermal = merge(biomass, T_tissue_data, by=c("taxon", "treatment"))
    bm = biomass_thermal %>% select(taxon, treatment, temp = T_tissue, mass = biomass_total, time = elapsed_days)
    
    bm$taxon = as.factor(bm$taxon)
    bm$treatment = as.factor(bm$treatment)
    # Preprocess data
    temps = seq(2.5, 35, 0.5)
    r_boof = data.frame(temp = temps, growth_rate = pawar_2018(temps, 
                                                               coef(bm.nlme)$J_ref.taxonBOOF[1],
                                                               coef(bm.nlme)$E.taxonBOOF[1],
                                                               coef(bm.nlme)$E_D.taxonBOOF[1],
                                                               coef(bm.nlme)$T_opt.taxonBOOF[1],
                                                               tref = 10), taxon = "BOOF")
    r_hovu = data.frame(temp = temps, growth_rate = pawar_2018(temps, 
                                                               coef(bm.nlme)$J_ref.taxonHOVU[1],
                                                               coef(bm.nlme)$E.taxonHOVU[1],
                                                               coef(bm.nlme)$E_D.taxonHOVU[1],
                                                               coef(bm.nlme)$T_opt.taxonHOVU[1],
                                                               tref = 10), taxon = "HOVU")
    r_rasa = data.frame(temp = temps, growth_rate = pawar_2018(temps, 
                                                               coef(bm.nlme)$J_ref.taxonRASA[1],
                                                               coef(bm.nlme)$E.taxonRASA[1],
                                                               coef(bm.nlme)$E_D.taxonRASA[1],
                                                               coef(bm.nlme)$T_opt.taxonRASA[1],
                                                               tref = 10), taxon = "RASA")
    
    r_all = bind_rows(r_boof, r_hovu, r_rasa)

    # Make Arrhenius data
    r_all = r_all %>% 
      mutate(invtemp = 1/(8.617e-5*(temp+273.15))) %>% 
      mutate(lnGrowthRate = log(growth_rate))
    # 
    #print("g")
    aBOOF = coef(bm.nlme)$`alpha.taxonBOOF`[1] # 0.6535
    aHOVU = coef(bm.nlme)$`alpha.taxonHOVU`[1] # 0.6621
    aRASA = coef(bm.nlme)$`alpha.taxonRASA`[1] # 0.2139
    #print("h")
    datBOOF = subset(bm,taxon=="BOOF")
    #print(datBOOF)
    datBOOF$aBOOF = aBOOF
    fitBOOF = nlsList(mass ~ growth_model(time, r, M_0, aBOOF) | temp,
                      data = datBOOF,
                      start = c(r = 0.1, M_0 = 0.5))
    #print("i")
    #print(fitBOOF)
    #print(coef(fitBOOF))
    plotBOOF = coef(fitBOOF) %>% mutate(temp = as.numeric(rownames(.)), taxon = "BOOF")
    print("j")
    datHOVU = subset(bm,taxon=="HOVU")
    datHOVU$aHOVU = aHOVU
    fitHOVU = nlsList(mass ~ growth_model(time, r, M_0, aHOVU) | temp,
                      data = datHOVU,
                      start = c(r = 0.1, M_0 = 0.5))
    plotHOVU = coef(fitHOVU) %>% mutate(temp = as.numeric(rownames(.)), taxon = "HOVU")
    
    datRASA = subset(bm,taxon=="RASA")
    datRASA$aRASA = aRASA
    fitRASA = nlsList(mass ~ growth_model(time, r, M_0, aRASA) | temp,
                      data = datRASA,
                      start = c(r = 0.1, M_0 = 0.5))
    plotRASA = coef(fitRASA) %>% mutate(temp = as.numeric(rownames(.)), taxon = "RASA")
    
    
    plotRASA$low = plotRASA$high = NA
    for (i in 1:length(fitRASA)) {
      cis = fitRASA[[i]] %>% confint2() %>% as.data.frame() %>% rename(low=`2.5 %`, high=`97.5 %`)
      plotRASA[i,]$low = cis["r","low"]
      plotRASA[i,]$high = cis["r","high"]
    }
    
    plotHOVU$low = plotHOVU$high = NA
    for (i in 1:length(fitHOVU)) {
      cis = fitHOVU[[i]] %>% confint2() %>% as.data.frame() %>% rename(low=`2.5 %`, high=`97.5 %`)
      plotHOVU[i,]$low = cis["r","low"]
      plotHOVU[i,]$high = cis["r","high"]
    }
    
    plotBOOF$low = plotBOOF$high = NA
    for (i in 1:length(fitBOOF)) {
      cis = fitBOOF[[i]] %>% confint2() %>% as.data.frame() %>% rename(low=`2.5 %`, high=`97.5 %`)
      plotBOOF[i,]$low = cis["r","low"]
      plotBOOF[i,]$high = cis["r","high"]
    }
    
    plot_r_all = bind_rows(plotBOOF, plotHOVU, plotRASA)
    
    #pdf("figures/fig3.pdf",width = 4, height = 3.5)
    z1 = plot_r_all %>% subset(taxon == "BOOF")
    z2 = r_all %>% subset(taxon == "BOOF")
    p1c = ggplot(data = z1, aes(x = temp, y = r)) +
      geom_point(size = 0.5) +
      geom_errorbar(aes(ymin=low,ymax=high), width=0) +
      geom_line(data = z2, aes(x = temp, y = growth_rate), color = "#619CFF") +
      my_theme +
      theme_classic() +
      theme_transparent +
      xlim(0,35) +
      ylim(-0.011,0.2) +
      #theme(legend.position = c(0.15, 0.8)) +
      #xlab("Whole-plant growth temperature (ºC)") +
      theme(axis.title.x = element_blank()) +#,
            #axis.title.y = element_blank()) #+
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      ylab(" ") +
      labs(title = "g.")
    
    z1 = plot_r_all %>% subset(taxon == "HOVU")
    z2 = r_all %>% subset(taxon == "HOVU")
    p2c = ggplot(data = z1, aes(x = temp, y = r)) +
      geom_point(size = 0.5) +
      geom_errorbar(aes(ymin=low,ymax=high), width=0) +
      geom_line(data = z2, aes(x = temp, y = growth_rate), color = "#619CFF") +
      my_theme +
      theme_classic() +
      theme_transparent +
      xlim(0,35) +
      ylim(-0.011,0.2) +
      #theme(legend.position = c(0.15, 0.8)) +
      #xlab("Whole-plant growth temperature (ºC)") +
      theme(axis.title.x = element_blank()) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      ylab(bquote("Mass-normalized growth rate,"~italic(G[M])~"("*g^"1-α"~d^"-1"*")"))+
      labs(title = "h.")
    
    
    z1 = plot_r_all %>% subset(taxon == "RASA")
    z2 = r_all %>% subset(taxon == "RASA")
    p3c = ggplot(data = z1, aes(x = temp, y = r)) +
      geom_point(size = 0.5) +
      geom_errorbar(aes(ymin=low,ymax=high), width=0) +
      geom_line(data = z2, aes(x = temp, y = growth_rate), color = "#619CFF") +
      my_theme +
      theme_classic() +
      theme_transparent +
      xlim(0,35) +
      ylim(-0.011,0.2) +
      #theme(legend.position = c(0.15, 0.8)) +
      xlab("Growth temperature (ºC)") +
      ylab("     ") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      labs(title = "i.")
      #theme(axis.title.x = element_blank()) +
      #theme(axis.title.y = element_blank())
      #ylab("Mass-normalized growth rate (g^(1-α) day^-1)")
    
    
    #grid.arrange(p)
    #dev.off()
  #   
  # svg("figures/test.svg", width = 8, height = 7)
  # grid.arrange(p1a,p1b,p1c,p2a,p2b,p2c,p3a,p3b,p3c,ncol=3)
  # dev.off()
  #   
  svg("figures/main_fig_2.svg", width = 8, height = 7)
  #grid.arrange(p1a,p1b,p1c,p2a,p2b,p2c,p3a,p3b,p3c,nrow = 3)
  p = plot_grid(p1a,p1b,p1c,p2a,p2b,p2c,p3a,p3b,p3c,align="hv")
  grid.arrange(p)
  dev.off()
  
}


growth_model = function(time, r, M_0, beta) {
  biomass = ( M_0^(1-beta) + r*time*(1-beta) )^(1/(1-beta))
  return(biomass)
}
