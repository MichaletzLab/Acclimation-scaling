make_supp_fig_6 = function() {#function(resp.nlme, bm.nlme) {

  
  # Next, we need to evaluate the model at the tissue T for each curve
  # The easiest way to do this is probably to use predict()
  resp_data %>% select(curveID, taxon, treatment) %>% unique() %>% merge(T_tissue_data) %>% 
    rename(Tleaf = T_tissue) %>% mutate(resp = predict(resp.nlme,.)) -> resp_composite
  
  # Now, to fit an Arrhenius function to each composite curve
  resp.nls.boof = nls(resp ~ arrhenius(Tleaf, r_tref, e, T_ref = 10),
                       data = subset(resp_composite, taxon == "BOOF"),
                       start = c(r_tref = 5, e = 0.5))
  
  resp.nls.hovu = nls(resp ~ arrhenius(Tleaf, r_tref, e, T_ref = 10),
                       data = subset(resp_composite, taxon == "HOVU"),
                       start = c(r_tref = 5, e = 0.5))
  
  resp.nls.rasa = nls(resp ~ arrhenius(Tleaf, r_tref, e, T_ref = 10),
                 data = subset(resp_composite, taxon == "RASA"),
                 start = c(r_tref = 10, e = 1))

  strip_data = bind_rows(
    
    # data.frame(low = intervals(resp.nlme, which = "fixed")$fixed[4:6,1], 
    #            E = intervals(resp.nlme, which = "fixed")$fixed[4:6,2], 
    #            high = intervals(resp.nlme, which = "fixed")$fixed[4:6,3],
    #            taxon = c("B. officinalis", "H. vulgare", "R. sativus"),
    #            rate = "Resp",
    #            acclimated = "Acute"),
    
    data.frame(low = intervals(resp.nlme, which = "fixed")$fixed[4:6,1], 
               E = intervals(resp.nlme, which = "fixed")$fixed[4:6,2], 
               high = intervals(resp.nlme, which = "fixed")$fixed[4:6,3],
               taxon = c("B. officinalis", "H. vulgare", "R. sativus"),
               rate = "resp_net"), # Boof E
    
    # data.frame(low = intervals(resp.gross.nlme.3, which = "fixed")$fixed[4:6,1], 
    #            E = intervals(resp.gross.nlme.3, which = "fixed")$fixed[4:6,2], 
    #            high = intervals(resp.gross.nlme.3, which = "fixed")$fixed[4:6,3],
    #            taxon = c("B. officinalis", "H. vulgare", "R. sativus"),
    #            rate = "resp_gross",
    #            acclimated = "Acute"), # Boof E
    
    data.frame(low = c(confint2(resp.nls.boof)[2,1],confint2(resp.nls.hovu)[2,1],confint2(resp.nls.rasa)[2,1]),
               E = c(coef(resp.nls.boof)["e"],coef(resp.nls.hovu)["e"],coef(resp.nls.rasa)["e"]),
               high = c(confint2(resp.nls.boof)[2,2],confint2(resp.nls.hovu)[2,2],confint2(resp.nls.rasa)[2,2]),
               taxon = c("B. officinalis", "H. vulgare", "R. sativus"),
               rate = "resp_net_accl"), # Boof E
    
    # data.frame(low = c(confint2(resp.gross.nls.boof)[2,1],confint2(resp.gross.nls.hovu)[2,1],confint2(resp.gross.nls.rasa)[2,1]),
    #            E = c(coef(resp.gross.nls.boof)["e"],coef(resp.gross.nls.hovu)["e"],coef(resp.gross.nls.rasa)["e"]),
    #            high = c(confint2(resp.gross.nls.boof)[2,2],confint2(resp.gross.nls.hovu)[2,2],confint2(resp.gross.nls.rasa)[2,2]),
    #            taxon = c("B. officinalis", "H. vulgare", "R. sativus"),
    #            rate = "resp_gross",
    #            acclimated = "Acclimated"), # Boof E
    
    # data.frame(low = c(confint2(resp.nls.boof)[2,1],confint2(resp.nls.hovu)[2,1],confint2(resp.nls.rasa)[2,1]),
    #            E = c(coef(resp.nls.boof)["E"],coef(resp.nls.hovu)["E"],coef(resp.nls.rasa)["E"]),
    #            high = c(confint2(resp.nls.boof)[2,2],confint2(resp.nls.hovu)[2,2],confint2(resp.nls.rasa)[2,2]),
    #            taxon = c("B. officinalis", "H. vulgare", "R. sativus"),
    #            rate = "Resp",
    #            acclimated = "Acclimated"),
    
    data.frame(low = intervals(bm.nlme)$fixed[4:6,1], 
               E = intervals(bm.nlme)$fixed[4:6,2], 
               high = intervals(bm.nlme)$fixed[4:6,3],
               taxon = c("B. officinalis", "H. vulgare", "R. sativus"),
               rate = "Growth")
  )
  
  
  
  
  
  
  
  
  
  strip_data$rate = as.factor(strip_data$rate)
  strip_data$rate = factor(strip_data$rate, c("Growth", "resp_net", "resp_net_accl"))
  #strip_data$acclimated = factor(strip_data$acclimated, c("Acute", "Acclimated"))
  color_vals = c("#619CFF", "#F8766D", "#a8362d")
  
  #strip_data = subset(strip_data, rate != "Resp")
  
  #svg("figures/fighyp_test.svg", width = 3, height = 4)
  p1 = ggplot(strip_data,aes(x = taxon, y = E, fill = rate)) +
    my_theme +
    theme_classic() +
    #theme_transparent +
    theme(legend.position = "none") +
    scale_fill_manual(values = color_vals, labels = c(bquote(G[M]),bquote(R[d]),bquote(R*"′"[d]))) +
    #geom_jitter(position=position_dodge(0.8),pch=21, size=2) +
    #geom_errorbar(aes(ymin = low, ymax = high), width = 0.2, position = position_dodge(0.5)) +
    # geom_bar_pattern(stat="identity", width = 0.8, position=position_dodge(0.8), color="black", 
    #                  pattern_fill = "black",
    #                  pattern_angle = 45,
    #                  pattern_density = 0.1,
    #                  pattern_spacing = 0.025,
    #                  pattern_key_scale_factor = 0.6) +
    #scale_pattern_manual(values = c(Acclimated = "stripe", Acute = "none"), labels=c("Acute", "Accl.")) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(0.8),color = "white") +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0, position=position_dodge(0.8)) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    #theme(legend.position = c(0.5, 0.93)) +
    # theme(legend.position = "top") +
    theme(legend.direction = "horizontal") +
    # theme(legend.title = element_blank()) +
    # theme(legend.background=element_blank()) + #
    # theme(legend.key = element_blank()) + #
    theme(legend.title = element_blank()) +
    theme(legend.background=element_blank()) + #
    theme(legend.key = element_blank()) + #
    theme(legend.key.size = unit(0.4, 'cm')) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #theme(axis.text.x = element_blank()) +
    #theme(legend.key.size = unit(0.4, 'cm')) +
    theme(legend.position = c(0.5,1)) +
    ylim(0,1) +
    geom_abline(slope = 0, intercept = 0) +
    geom_abline(slope = 0, intercept = 0.32, lty=2) +
    geom_abline(slope = 0, intercept = 0.65, lty=3) +
    ylab(bquote("E"[A]~"(eV)")) +
    xlab(NULL) #+
    #labs(title = "a.")
  # Okay, this nicely visualizes the basic 95CI comparison that we want to do, I think...
  #grid.arrange(p)
  #dev.off()
  p1
  
  svg("figures/supp_fig_6.svg", width = 3, height = 3)
  #grid.arrange(p1,p2,p3,ncol=1)
  grid.arrange(p1)
  dev.off()
  

  
}