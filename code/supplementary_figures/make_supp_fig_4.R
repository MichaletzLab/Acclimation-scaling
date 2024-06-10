# Biomass over time mega-plot
make_supp_fig_4 = function() {#function(bm_data, bm.nlme, T_tissue_data) {
  
  biomass = bm_data %>% rename(treatment = temperature_C, taxon = species)
  biomass_thermal = merge(biomass, T_tissue_data, by=c("taxon", "treatment"))
  bm = biomass_thermal %>% select(taxon, treatment, temp = T_tissue, mass = biomass_total, time = elapsed_days)
  
  
  bm$taxon = as.factor(bm$taxon)
  bm$treatment = as.factor(bm$treatment)
  # Preprocess data
  
  bm = bm %>% select(taxon, treatment, temp) %>% unique()
  
  bm_pred = c()
  for (i in 1:dim(bm)[1]) {
    curdat = bm[i,]
    bm_pred = bind_rows(bm_pred,
                        data.frame(curdat, time = 0:28))
  }
  bm_pred$mass = predict(bm.nlme, bm_pred)
  
  bm_pred$species = bm_pred$taxon
  bm_pred$temperature_C = as.numeric(as.character(bm_pred$treatment))
    
  svg("figures/supp_fig_4.svg", width = 8, height = 5)
  p = ggplot(data = bm_data, aes(x = elapsed_days, y = biomass_total, color = species)) +
    my_theme +
    theme_transparent +
    geom_point() +
    geom_line(data = bm_pred, aes(x = time, y = mass)) +
    scale_color_brewer(palette = "Dark2") +
    xlab("Elapsed time (days)") +
    ylab("Dry biomass (g)") +
    facet_grid(rows = vars(species), cols = vars(temperature_C),
               labeller = as_labeller(c(`2.5` = "2.5 ºC",
                                        `5` = "5 ºC",
                                        `7.5` = "7.5 ºC",
                                        `10` = "10 ºC",
                                        `15` = "15 ºC",
                                        `20` = "20 ºC",
                                        `25` = "25 ºC",
                                        `30` = "30 ºC",
                                        `35` = "35 ºC",
                                        BOOF = "B. officinalis",
                                        HOVU = "H. vulgare",
                                        RASA = "R. sativus")))
  grid.arrange(p)
  dev.off()

}
