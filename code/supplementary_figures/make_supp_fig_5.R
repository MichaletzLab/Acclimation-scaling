
make_supp_fig_5 = function() {#function(bm_data, bm.nlme, T_tissue_data) {
  
  
  biomass = bm_data %>% rename(treatment = temperature_C, taxon = species)
  biomass_thermal = merge(biomass, T_tissue_data, by=c("taxon", "treatment"))
  bm = biomass_thermal %>% select(taxon, treatment, temp = T_tissue, mass = biomass_total, time = elapsed_days)
  
  
  
  bm$taxon = as.factor(bm$taxon)
  bm$treatment = as.factor(bm$treatment)
  # Make 3 3-d plots
  
  
  # Let's try to normalize biomass in the dataframe s.t. M_0 = 1 for all treatments
  bm_model = function(M_0, beta, r, time) {
    return (( M_0^(1-beta) + r*time*(1-beta) )^(1/(1-beta)))
  }
  
  a = subset(bm, taxon == "BOOF")
  
  params = data.frame(treatment = as.numeric(rownames(coef(bm.nlme))),
                      M_0 = coef(bm.nlme)$M_0.taxonBOOF + coef(bm.nlme)$`M_0.(Intercept)`,
                      beta = coef(bm.nlme)$alpha.taxonBOOF)
  
  a = merge(a, params, by ="treatment")
  a = a %>% mutate(mnorm = (mass^(1-beta) - M_0^(1-beta) + 1)^(1/(1-beta)) )
  
  xpred = 2:35 # temperature treatment
  ypred = 0:30 # time days
  zpred = biomass_model(xpred, 
                        sort(rep(0:30,34)), 
                        coef(bm.nlme)$J_ref.taxonBOOF[1],
                        coef(bm.nlme)$E.taxonBOOF[1],
                        coef(bm.nlme)$E_D.taxonBOOF[1],
                        coef(bm.nlme)$T_opt.taxonBOOF[1],
                        1,
                        coef(bm.nlme)$alpha.taxonBOOF[1])
  
  svg("figures/supp_fig_5.svg", width = 8,height=4)
  layout(matrix(c(1:3),1))
  #pdf("figures/fig4_boof.pdf", width = 5, height = 5)
  scatter3D(a$temp, a$time, a$mnorm,
            main = "a.  B. officinalis",
            xlab = "Tissue temperature (ºC)",
            ylab = "Time  (days)",
            zlab = "Normalized biomass (g)",
            xlim = c(0,35),
            colvar = NULL,
            col = "black",
            pch = 19,
            cex = 0.5,
            cex.axis = 0.8,
            cex.lab = 0.9,
            ticktype = "detailed",
            phi = 10,
            theta = 20,
            surf = list(x = xpred, y = ypred, 
                        z = matrix(zpred, nrow = length(xpred), ncol = length(ypred)),
                        facets = NA, col = "#619CFF"))
  #dev.off()
  
  
  
  a = subset(bm, taxon == "HOVU")
  
  params = data.frame(treatment = as.numeric(rownames(coef(bm.nlme))),
                      M_0 = coef(bm.nlme)$M_0.taxonHOVU + coef(bm.nlme)$`M_0.(Intercept)`,
                      beta = coef(bm.nlme)$alpha.taxonHOVU)
  
  a = merge(a, params, by ="treatment")
  a = a %>% mutate(mnorm = (mass^(1-beta) - M_0^(1-beta) + 1)^(1/(1-beta)) )
  
  xpred = 2:35 # temperature treatment
  ypred = 0:30 # time days
  zpred = biomass_model(xpred, 
                        sort(rep(0:30,34)), 
                        coef(bm.nlme)$J_ref.taxonHOVU[1],
                        coef(bm.nlme)$E.taxonHOVU[1],
                        coef(bm.nlme)$E_D.taxonHOVU[1],
                        coef(bm.nlme)$T_opt.taxonHOVU[1],
                        1,
                        coef(bm.nlme)$alpha.taxonHOVU[1])
  
  #pdf("figures/fig4_hovu.pdf", width = 6, height = 6)
  scatter3D(a$temp, a$time, a$mnorm,
            main = "b.  H. vulgare",
            xlab = "Tissue temperature (ºC)",
            ylab = "Time (days)",
            zlab = "Normalized biomass (g)",
            xlim = c(0,35),
            colvar = NULL,
            col = "black",
            pch = 19,
            cex = 0.5,
            ticktype = "detailed",
            phi = 10,
            theta = 20,
            surf = list(x = xpred, y = ypred, 
                        z = matrix(zpred, nrow = length(xpred), ncol = length(ypred)),
                        facets = NA, col = "#619CFF"))
  #dev.off()
  
  
  
  
  
  a = subset(bm, taxon == "RASA")
  
  params = data.frame(treatment = as.numeric(rownames(coef(bm.nlme))),
                      M_0 = coef(bm.nlme)$M_0.taxonRASA + coef(bm.nlme)$`M_0.(Intercept)`,
                      beta = coef(bm.nlme)$alpha.taxonRASA)
  
  a = merge(a, params, by ="treatment")
  a = a %>% mutate(mnorm = (mass^(1-beta) - M_0^(1-beta) + 1)^(1/(1-beta)) )
  
  xpred = 2:35 # temperature treatment
  ypred = 0:30 # time days
  zpred = biomass_model(xpred, 
                        sort(rep(0:30,34)), 
                        coef(bm.nlme)$J_ref.taxonRASA[1],
                        coef(bm.nlme)$E.taxonRASA[1],
                        coef(bm.nlme)$E_D.taxonRASA[1],
                        coef(bm.nlme)$T_opt.taxonRASA[1],
                        1,
                        coef(bm.nlme)$alpha.taxonRASA[1])
  
  #pdf("figures/fig4_rasa.pdf", width = 6, height = 6)
  scatter3D(a$temp, a$time, a$mnorm,
            main = "c.  R. sativus",
            xlab = "Tissue temperature (ºC)",
            ylab = "Time (days)",
            zlab = "Normalized biomass (g)",
            xlim = c(0,35),
            colvar = NULL,
            col = "black",
            pch = 19,
            cex = 0.5,
            ticktype = "detailed",
            phi = 10,
            theta = 20,
            surf = list(x = xpred, y = ypred, 
                        z = matrix(zpred, nrow = length(xpred), ncol = length(ypred)),
                        facets = NA, col = "#619CFF"))
  dev.off()
  
}



# make_fig4 = function() {
#   
#   photo_fits = read.csv("intermediate_data/photo_fits.csv")   
#   
#   # Dependence of activaition energy on treatment temp
#   p1 = ggplot(data = photo_fits, aes(x = treatment, y = E, color = taxon)) + 
#     geom_smooth(method=lm, aes(color=NULL), color = "black") +
#     geom_point() +
#     xlab("Treatment temperature (ºC)") +
#     ylab("Activation energy (eV)") +
#     ylim(0,3) +
#     ggtitle("Activation energy of photosynthesis, all species")
#   
#   # Dependence of deactivaition energy on treatment temp
#   p2 = ggplot(data = photo_fits, aes(x = treatment, y = E_D, color = taxon)) + 
#     geom_smooth(method=lm, aes(color=NULL), color = "black") +
#     geom_point() +
#     ylim(0,3.75) +
#     xlab("Treatment temperature (ºC)") +
#     ylab("Deactivation energy (eV)") +
#     ggtitle("Deactivation energy of photosynthesis, all species")
#   
#   # Dependence of norm const on treatment temp
#   p3 = ggplot(data = photo_fits, aes(x = treatment, y = J_ref, color = taxon)) + 
#     geom_smooth(method=lm, aes(color=NULL), color = "black") +
#     geom_point() +
#     xlab("Treatment temperature (ºC)") +
#     ylab("ln(Normalization constant (umol/m2s))") +
#     ggtitle("Normalization constant of photosynthesis, all species")
#   
#   # Dependence of topt on treatment temp
#   p4 = ggplot(data = photo_fits, aes(x = treatment, y = T_opt, color = taxon)) + 
#     geom_smooth(method=lm, aes(color=NULL), color = "black") +
#     geom_point() +
#     ylim(10,40) +
#     xlab("Treatment temperature (ºC)") +
#     ylab("Optimum temperature (ºC)") +
#     ggtitle("Optimum temperature of photosynthesis, all species")
#   
#   grid.arrange(p1,p2,p3,p4,nrow=2)
#   
# }
