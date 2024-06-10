# New supplementary figure illustrating isormality
make_supp_fig_1 = function() {

all_files = list.files("data/thermal_data_new", recursive = T, full.names = T) %>% 
  grep(".csv", ., value = T)

data_raw = lapply(all_files, read.csv)

data_unpack = c()
for (i in 1:length(data_raw)) {
  time_start = max(data_raw[[i]]$time_elapsed) - 24*60*60 # Only take the last 24 hourse of data collection time
  data_unpack = bind_rows(data_unpack, 
                          data.frame(subset(data_raw[[i]], time_elapsed >= time_start), filename = all_files[i], curveID = i))
}

thermal_data = data_unpack %>% 
  select(thermal_mean, thermal_median = thermal_q50, temp_atm, temp_ambient, temp_soil_mean, relative_humidity, filename, curveID) %>% 
  separate(col = filename, into = c("junk1", "junk2", "metadata"), sep = "/") %>% 
  separate(col = metadata, into = c("spec_rep", "week", "treat"), sep = "_") %>% 
  select(-junk1, -junk2)

thermal_data$replicate = thermal_data$spec_rep %>% gsub("*[a-z]", "", .)
thermal_data$taxon = thermal_data$spec_rep %>% gsub("*[0-9]", "", .)
thermal_data$week = thermal_data$week %>% gsub("*[a-z]", "", .) %>% as.integer()
thermal_data$treatment = thermal_data$treat %>% trimws() %>% str_match("[0-9]*") %>% as.numeric()

thermal_data = thermal_data %>% select(-spec_rep, -treat)
# 
# thermal_means = thermal_data %>% group_by(taxon, treatment) %>% 
#   summarize(#Taxon = unique(species),
#     #Replicate = unique(replicate),
#     #Treatment = unique(treatment),
#     #Week = unique(week),
#     T_leaf_mean = mean(thermal_mean)-273.15,
#     T_leaf_median = mean(thermal_median)-273.15,
#     T_air_1 = mean(temp_atm)-273.15,
#     T_air_2 = mean(temp_ambient)-273.15,
#     T_root_mean = mean(temp_soil_mean)-273.15,
#     RH_mean = mean(relative_humidity))

thermal_data = thermal_data %>% mutate(T_leaf_mean = thermal_mean-273.15,
                        T_leaf_median = thermal_median-273.15,
                        T_air_1 = temp_atm-273.15,
                        T_air_2 = temp_ambient-273.15,
                        T_root_mean = temp_soil_mean-273.15) %>% 
  select(-thermal_mean, -thermal_median, -temp_atm, -temp_ambient, -temp_soil_mean)


thermal_data$taxon = toupper(thermal_data$taxon)
thermal_data[thermal_data$treatment == 2,]$treatment = 2.5
thermal_data[thermal_data$treatment == 7,]$treatment = 7.5

df = thermal_data[seq(1, nrow(thermal_data), 10), ]

thermal_data %>% group_by(treatment, taxon) %>% 
  summarize(mean_T_leaf = mean(T_leaf_mean, na.rm=T),
            sd_T_leaf = sd(T_leaf_mean, na.rm=T),
            mean_T_soil = mean(T_root_mean, na.rm = T),
            sd_T_soil = sd(T_root_mean)) -> df2


bm_data %>% 
  select(taxon = species, treatment = temperature_C, prop_aboveground) %>% 
  group_by(taxon, treatment) %>% 
  summarize(mean_prop = mean(prop_aboveground, na.rm = T)) -> bm_avs

merge(thermal_data, bm_avs) -> df3
df3 %>% 
  mutate(T_tissue = T_leaf_mean*mean_prop + T_root_mean*(1-mean_prop)) %>% 
  group_by(treatment) %>% 
  summarize(mean_T = mean(T_tissue, na.rm = T),
            sd_T = sd(T_tissue, na.rm = T)) -> df4


p = ggplot(df4, aes(x = treatment, y = mean_T)) +
  geom_point(position=) +
  geom_errorbar(aes(ymax = mean_T+sd_T, ymin = mean_T-sd_T), color="black", width = 0) +
  geom_abline(slope = 1, lty = 2) +
  my_theme +
  xlab("Treatment temperature (°C)") +
  ylab("Mean tissue temperature (°C)") +
  xlim(0,35) +
  ylim(0,35)

svg("figures/supp_fig_1.svg", width = 4, height = 3.7)
grid.arrange(p)
dev.off()

}
