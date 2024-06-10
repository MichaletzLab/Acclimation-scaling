make_supp_fig_7 = function() {
    
    biomass = bm_data %>% rename(treatment = temperature_C, taxon = species)
    biomass_thermal = merge(biomass, T_tissue_data, by=c("taxon", "treatment"))
    bm = biomass_thermal %>% select(taxon, treatment, temp = T_tissue, mass = biomass_total, time = elapsed_days)
    
    
    comp_pred = photo_data %>% select(curveID, taxon, treatment) %>% unique()
    comp_pred = bm %>% select(taxon, treatment, temp) %>% unique() %>% merge(comp_pred) %>% rename(Tleaf = temp)
    comp_pred$Photo = predict(photo.nlme,comp_pred)
    
    growth_pred = c()
    
    xpred = subset(bm, taxon == "BOOF")$temp %>% unique() # Tissue T
    ypred = pawar_2018(xpred, 
                       coef(bm.nlme)$J_ref.taxonBOOF[1],
                       coef(bm.nlme)$E.taxonBOOF[1],
                       coef(bm.nlme)$E_D.taxonBOOF[1],
                       coef(bm.nlme)$T_opt.taxonBOOF[1],
                       tref = 10)
    
    growth_pred = data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "BOOF")
    
    xpred = subset(bm, taxon == "HOVU")$temp %>% unique() # Tissue T
    ypred = pawar_2018(xpred, 
                       coef(bm.nlme)$J_ref.taxonHOVU[1],
                       coef(bm.nlme)$E.taxonHOVU[1],
                       coef(bm.nlme)$E_D.taxonHOVU[1],
                       coef(bm.nlme)$T_opt.taxonHOVU[1],
                       tref = 10)
    
    growth_pred = bind_rows(growth_pred, data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "HOVU"))
    
    xpred = subset(bm, taxon == "RASA")$temp %>% unique() # Tissue T
    ypred = pawar_2018(xpred, 
                       coef(bm.nlme)$J_ref.taxonRASA[1],
                       coef(bm.nlme)$E.taxonRASA[1],
                       coef(bm.nlme)$E_D.taxonRASA[1],
                       coef(bm.nlme)$T_opt.taxonRASA[1],
                       tref = 10)
    
    growth_pred = bind_rows(growth_pred, data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "RASA"))
    
    photo_growth = comp_pred %>% 
      group_by(taxon, Tleaf) %>% 
      summarize(meanPhoto = mean(Photo)) %>% 
      merge(growth_pred, by = c("taxon", "Tleaf"))
    
    #pdf("figures/fig12a.pdf", width = 2, height = 1.75)
    p1 = ggplot(photo_growth, aes(x=meanPhoto, y=GrowthRate, color=taxon)) + 
      geom_point() +
      #geom_smooth(method = "lm", aes(color = NULL), color = "black") +
      geom_smooth(method = "lm", linewidth = 0.75, alpha = 0.2) +
      scale_color_brewer(palette = "Dark2", 
                         labels=c('B. officinalis',"H. vulgare", "R. sativus")) +
      scale_fill_brewer(palette = "Dark2",
                         labels=c('B. officinalis',"H. vulgare", "R. sativus")) +
      my_theme +
      theme_transparent +
      theme(legend.position = c(0.8,0.2),
            legend.title = element_blank()) +
      #scale_color_discrete(labels=c('B. officinalis',"H. vulgare", "R. sativus")) +
      #theme(legend.position = "none") +
      xlim(0,25) +
      ylim(0,0.18) +
      xlab(bquote("Assimilation rate,"~italic(A*"′"[n])~"(µmol "~m^"-2"~s^"-1"*")")) +
      ylab(bquote("Mass-normalized growth rate,"~italic(G[M])~"("*g^"1-α"~d^"-1"*")")) +
      labs(tag = "b.")
    #grid.arrange(p)
    #dev.off()
    
    z = photo_growth %>% subset(taxon == "BOOF") %>% lm(GrowthRate ~ meanPhoto, .)
    summary(z)
    z = photo_growth %>% subset(taxon == "HOVU") %>% lm(GrowthRate ~ meanPhoto, .)
    summary(z)
    z = photo_growth %>% subset(taxon == "RASA") %>% lm(GrowthRate ~ meanPhoto, .)
    summary(z)
  
  biomass = bm_data %>% rename(treatment = temperature_C, taxon = species)
  biomass_thermal = merge(biomass, T_tissue_data, by=c("taxon", "treatment"))
  bm = biomass_thermal %>% select(taxon, treatment, temp = T_tissue, mass = biomass_total, time = elapsed_days)
  
  
  comp_pred = resp_data %>% select(curveID, taxon, treatment) %>% unique()
  comp_pred = bm %>% select(taxon, treatment, temp) %>% unique() %>% merge(comp_pred) %>% rename(Tleaf = temp)
  comp_pred$Photo = predict(resp.nlme,comp_pred)
  
  growth_pred = c()
  
  xpred = subset(bm, taxon == "BOOF")$temp %>% unique() # Tissue T
  ypred = pawar_2018(xpred, 
                     coef(bm.nlme)$J_ref.taxonBOOF[1],
                     coef(bm.nlme)$E.taxonBOOF[1],
                     coef(bm.nlme)$E_D.taxonBOOF[1],
                     coef(bm.nlme)$T_opt.taxonBOOF[1],
                     tref = 10)
  
  growth_pred = data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "BOOF")
  
  xpred = subset(bm, taxon == "HOVU")$temp %>% unique() # Tissue T
  ypred = pawar_2018(xpred, 
                     coef(bm.nlme)$J_ref.taxonHOVU[1],
                     coef(bm.nlme)$E.taxonHOVU[1],
                     coef(bm.nlme)$E_D.taxonHOVU[1],
                     coef(bm.nlme)$T_opt.taxonHOVU[1],
                     tref = 10)
  
  growth_pred = bind_rows(growth_pred, data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "HOVU"))
  
  xpred = subset(bm, taxon == "RASA")$temp %>% unique() # Tissue T
  ypred = pawar_2018(xpred, 
                     coef(bm.nlme)$J_ref.taxonRASA[1],
                     coef(bm.nlme)$E.taxonRASA[1],
                     coef(bm.nlme)$E_D.taxonRASA[1],
                     coef(bm.nlme)$T_opt.taxonRASA[1],
                     tref = 10)
  
  growth_pred = bind_rows(growth_pred, data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "RASA"))
  
  resp_growth = comp_pred %>% 
    group_by(taxon, Tleaf) %>% 
    summarize(meanPhoto = mean(Photo)) %>% 
    merge(growth_pred, by = c("taxon", "Tleaf"))
  
  #pdf("figures/fig13.pdf", width = 4, height = 3.5)
  p2 = ggplot(resp_growth, aes(x=meanPhoto, y=GrowthRate, color=taxon)) + 
    geom_smooth(method = "lm", linewidth = 0.75, alpha = 0.2,linetype = 2) +
    geom_point() +
    scale_color_brewer(palette = "Dark2", 
                       labels=c('B. officinalis',"H. vulgare", "R. sativus")) +
    scale_fill_brewer(palette = "Dark2",
                      labels=c('B. officinalis',"H. vulgare", "R. sativus")) +
    my_theme +
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    theme_transparent +
    #theme(legend.position = c(0.8,0.2)) +
    xlab(bquote("Respiration rate,"~italic(R[d])~"(µmol "~m^"-2"~s^"-1"*")")) +
    ylab(bquote("Mass-normalized growth rate,"~italic(G[M])~"("*g^"1-α"~d^"-1"*")")) #+
    #labs(tag = "c.")
  #grid.arrange(p)
  #dev.off()
  
  
  z = resp_growth %>% subset(taxon == "BOOF") %>% lm(GrowthRate ~ meanPhoto, .)
  summary(z)
  z = resp_growth %>% subset(taxon == "HOVU") %>% lm(GrowthRate ~ meanPhoto, .)
  summary(z)
  z = resp_growth %>% subset(taxon == "RASA") %>% lm(GrowthRate ~ meanPhoto, .)
  summary(z)
  #   
  # svg("figures/supp_fig_4.svg", width = 10, height = 3.3)
  # grid.arrange(p3,p1,p2, ncol=3)
  # dev.off()
  # 
  # 
  # 
  
  # Can we predict acute AT from photo
  
  
  growth_pred = c()
  
  xpred = subset(bm, taxon == "BOOF")$temp %>% unique() # Tissue T
  ypred = pawar_2018(xpred, 
                     coef(bm.nlme)$J_ref.taxonBOOF[1],
                     coef(bm.nlme)$E.taxonBOOF[1],
                     coef(bm.nlme)$E_D.taxonBOOF[1],
                     coef(bm.nlme)$T_opt.taxonBOOF[1],
                     tref = 10)
  
  growth_pred = data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "BOOF")
  
  xpred = subset(bm, taxon == "HOVU")$temp %>% unique() # Tissue T
  ypred = pawar_2018(xpred, 
                     coef(bm.nlme)$J_ref.taxonHOVU[1],
                     coef(bm.nlme)$E.taxonHOVU[1],
                     coef(bm.nlme)$E_D.taxonHOVU[1],
                     coef(bm.nlme)$T_opt.taxonHOVU[1],
                     tref = 10)
  
  growth_pred = bind_rows(growth_pred, data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "HOVU"))
  
  xpred = subset(bm, taxon == "RASA")$temp %>% unique() # Tissue T
  ypred = pawar_2018(xpred, 
                     coef(bm.nlme)$J_ref.taxonRASA[1],
                     coef(bm.nlme)$E.taxonRASA[1],
                     coef(bm.nlme)$E_D.taxonRASA[1],
                     coef(bm.nlme)$T_opt.taxonRASA[1],
                     tref = 10)
  
  growth_pred = bind_rows(growth_pred, data.frame(Tleaf = xpred, GrowthRate = ypred, taxon = "RASA"))
  
  
  
  xpred.b.at = subset(bm, taxon == "BOOF")$temp %>% unique() # Tissue T
  ypred.b.at = pawar_2018(xpred, 
                          coef(photo.nlme)$r_tref.taxonBOOF[1],
                          coef(photo.nlme)$e.taxonBOOF[1],
                          coef(photo.nlme)$eh.taxonBOOF[1],
                          coef(photo.nlme)$topt.taxonBOOF[1],
                          tref = 10)
  xpred.h.at = subset(bm, taxon == "HOVU")$temp %>% unique() # Tissue T
  ypred.h.at = pawar_2018(xpred, 
                          coef(photo.nlme)$r_tref.taxonHOVU[1],
                          coef(photo.nlme)$e.taxonHOVU[1],
                          coef(photo.nlme)$eh.taxonHOVU[1],
                          coef(photo.nlme)$topt.taxonHOVU[1],
                          tref = 10)
  xpred.r.at = subset(bm, taxon == "RASA")$temp %>% unique() # Tissue T
  ypred.r.at = pawar_2018(xpred, 
                          coef(photo.nlme)$r_tref.taxonRASA[1],
                          coef(photo.nlme)$e.taxonRASA[1],
                          coef(photo.nlme)$eh.taxonRASA[1],
                          coef(photo.nlme)$topt.taxonRASA[1],
                          tref = 10)
  
  
  acute_at = data.frame(Tleaf = c(xpred.b.at, xpred.h.at, xpred.r.at),
                        meanPhoto = c(ypred.b.at,ypred.h.at,ypred.r.at))
  
  
  acute_growth = merge(growth_pred, acute_at, by="Tleaf", all = T)
  sig = data.frame(taxon = c("BOOF", "HOVU", "RASA"),
                   sig = as.factor(c(1,0,0)))
  acute_growth_2 = merge(acute_growth, sig, by="taxon")
  
  p3 = ggplot(acute_growth_2, aes(x=meanPhoto, y=GrowthRate, color=taxon)) + 
    geom_smooth(method = "lm", linewidth = 0.75, alpha = 0.2, aes(linetype = sig)) +
    geom_point() +
    scale_linetype_manual(values=c(1,2)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    my_theme +
    theme_transparent +
    #theme(legend.position = c(0.8,0.2)) +
    xlab(bquote("Assimilation rate,"~italic(A[n])~"(µmol "~m^"-2"~s^"-1"*")")) +
    ylab(bquote("Mass-normalized growth rate,"~italic(G[M])~"("*g^"1-α"~d^"-1"*")")) +
    labs(tag = "a.")
  
  z = acute_growth %>% subset(taxon == "BOOF") %>% lm(GrowthRate ~ meanPhoto, .)
  summary(z)
  z = acute_growth %>% subset(taxon == "HOVU") %>% lm(GrowthRate ~ meanPhoto, .)
  summary(z)
  z = acute_growth %>% subset(taxon == "RASA") %>% lm(GrowthRate ~ meanPhoto, .)
  summary(z)
  
  # 
  # svg("figures/main_fig_5.svg", width = 7, height = 3.3)
  # grid.arrange(p3,p1, ncol=2)
  # dev.off()
  
  svg("figures/supp_fig_7.svg", width = 3.5, height = 4)
  grid.arrange(p2)
  dev.off()
  
}



