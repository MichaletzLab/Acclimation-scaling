# Build databases

load.all.data = function(from_raw = F) {
  
  # Check for the existence of intermediates
  
  # If from_raw flag = T, build databases from raw data
  if(from_raw == T) {
    
    # Build databases of major data streams
    resp_data <<- build_resp() # Done
    photo_data_unfiltered <<- build_photo_treatment() # Done
    T_tissue_data <<- build_T_tissue_data() # This is tissue T weighted by allocation
    bm_data <<- build_bm_data() # This is the "raw" biomass growth data
    
    # Bring in ancillary data
    lma.df <<- read.csv("data/lma_rehydrated.csv")
    
    # Do we need these? I think this gets done in the build_bm_data function
    rasa.leaf <<- read.csv("data/rasa_leafy_biomass_fixed.csv")
    hovu.leaf <<- read.csv("data/hovu_leafy_biomass.csv")
    
    # Clean lma data
    lma.df <<- lma.df %>% 
      subset(sample != 202 | leaf != 2) %>% 
      subset(sample != 215 | leaf != 3) %>% 
      subset(sample != 222 | leaf != 3) %>% 
      subset(sample != 226 | leaf != 3) %>% 
      subset(sample != 228 | leaf != 1) %>% 
      subset(sample != 228 | leaf != 2) %>% 
      subset(sample != 228 | leaf != 3) %>% 
      subset(sample != 231 | leaf != 2) %>% 
      subset(sample != 243 | leaf != 1)
    
    # Apply shrinkage correction factors
    lma.df$lma_corr_factor = NA
    lma.df[lma.df$species=="boof",]$lma_corr_factor = 41.7/29.3
    lma.df[lma.df$species=="hovu",]$lma_corr_factor = 1
    lma.df[lma.df$species=="rasa",]$lma_corr_factor = 36.6/24.4
    lma.df$lma_corr = lma.df$lma_g_m2/lma.df$lma_corr_factor
    
    # Need to convert RASA date to correct format to match other data
    
    # Clean photosynthesis data
    photo_data <<- filter_photo(photo_data_unfiltered, filter_from_raw = F)
    
    # Get additional physiological traits Vcmax, gs, and iWUE
    physio_data <<- get_physio(photo_data, resp_data, T_tissue_data)
    
    # Build statistical models
    resp.nlme <<- fit.resp.nlme(resp_data)
    photo.nlme <<- fit.photo.nlme(photo_data)
    bm.nlme <<- fit.bm.nlme(bm_data, T_tissue_data)
    
    # Trying this with adding the 35 treatment back in
    photo.nlme <<- fit.photo.nlme.with.35(photo_data, photo_data_unfiltered, photo.nlme)
    photo_data <<- bind_rows(photo_data, subset(photo_data_unfiltered, treatment == 35))
    
    # Save data once loaded into intermediate folder
    dir.create("preprocessed_data")
    
    write.csv(resp_data, "preprocessed_data/resp_data", row.names = F)
    write.csv(photo_data, "preprocessed_data/photo_data", row.names = F)
    write.csv(T_tissue_data, "preprocessed_data/T_tissue_data", row.names = F)
    write.csv(bm_data, "preprocessed_data/bm_data", row.names = F)
    
    write_rds(resp.nlme, "preprocessed_data/resp.nlme")
    write_rds(photo.nlme, "preprocessed_data/photo.nlme")
    write_rds(bm.nlme, "preprocessed_data/bm.nlme")
    
    # Otherwise, load data from intermediate saved files
  } else {
    
    resp_data <<- read.csv("preprocessed_data/resp_data")
    photo_data <<- read.csv("preprocessed_data/photo_data")
    T_tissue_data <<- read.csv("preprocessed_data/T_tissue_data")
    bm_data <<- read.csv("preprocessed_data/bm_data")
    
    resp.nlme <<- read_rds("preprocessed_data/resp.nlme")
    photo.nlme <<- read_rds("preprocessed_data/photo.nlme")
    bm.nlme <<- read_rds("preprocessed_data/bm.nlme")
    
    # Bring in ancillary data
    lma.df <<- read.csv("data/lma_rehydrated.csv")
    rasa.leaf <<- read.csv("data/rasa_leafy_biomass_fixed.csv")
    hovu.leaf <<- read.csv("data/hovu_leafy_biomass.csv")
    
    # Get additional physiological traits Vcmax, gs, and iWUE
    physio_data <<- get_physio(photo_data, resp_data, T_tissue_data)
    
    # Clean lma data
    lma.df <<- lma.df %>% 
      subset(sample != 202 | leaf != 2) %>% 
      subset(sample != 215 | leaf != 3) %>% 
      subset(sample != 222 | leaf != 3) %>% 
      subset(sample != 226 | leaf != 3) %>% 
      subset(sample != 228 | leaf != 1) %>% 
      subset(sample != 228 | leaf != 2) %>% 
      subset(sample != 228 | leaf != 3) %>% 
      subset(sample != 231 | leaf != 2) %>% 
      subset(sample != 243 | leaf != 1)
    
    # Apply shrinkage correction factors
    lma.df$lma_corr_factor = NA
    lma.df[lma.df$species=="boof",]$lma_corr_factor = 41.7/29.3
    lma.df[lma.df$species=="hovu",]$lma_corr_factor = 1
    lma.df[lma.df$species=="rasa",]$lma_corr_factor = 36.6/24.4
    lma.df$lma_corr = lma.df$lma_g_m2/lma.df$lma_corr_factor
    
  }
}

build_resp = function() {
  
  # Build a dataset of respiration data
  
  all_files = list.files("data/respiration_tpc_treatments", recursive = T, full.names = T) %>% 
    grep(".xlsx", ., value = T, invert = T)
  data_raw = lapply(all_files, read_6800)
  
  
  data_unpack = c()
  for (i in 1:length(data_raw)) {
    #data_unpack = rbind.fill(data_unpack, data.frame(rasa_data_raw[i], curveID = i))
    data_unpack = bind_rows(data_unpack, data.frame(data_raw[i], curveID = i))
  }
  
  
  # Correct taxon names
  data_unpack$taxon = NA
  data_unpack[grep("rasa", data_unpack$filename),]$taxon = "RASA"
  data_unpack[grep("boof", data_unpack$filename),]$taxon = "BOOF"
  data_unpack[grep("hovu", data_unpack$filename),]$taxon = "HOVU"
  
  # Pick out treatment
  data_unpack$treatment = NA
  data_unpack$filename %>% str_extract("_.{1,3}C") %>% 
    str_remove("_") %>% str_remove("C") %>% as.numeric() -> data_unpack$treatment
  
  # clear out junk, make resp positive
  tpc_respiration_data = data_unpack %>% 
    select(curveID, date, taxon, treatment, Photo = A, Tleaf)
  tpc_respiration_data$Photo = -tpc_respiration_data$Photo
  
  #tpc_respiration_data = subset(tpc_respiration_data, !is.na(Photo))
  
  return(tpc_respiration_data)
}




build_photo_treatment = function() {
  
  all_files = list.files("data/photo_tpc_treatments", recursive = T, full.names = T) %>% 
    grep(".xlsx", ., value = T, invert = T)
  data_raw = lapply(all_files, read_6800)
  
  data_unpack = c()
  for (i in 1:length(data_raw)) {
    #print(i)
    data_unpack = bind_rows(data_unpack, data.frame(data_raw[i], curveID = i))
  }
  
  # Two bad points to remove
  data_unpack = subset(data_unpack, curveID != 57 | obs != 10)
  data_unpack = subset(data_unpack, curveID != 113 | obs != 9)
  
  data_unpack$taxon = NA
  data_unpack[grep("rasa", data_unpack$filename),]$taxon = "RASA"
  data_unpack[grep("boof", data_unpack$filename),]$taxon = "BOOF"
  data_unpack[grep("hovu", data_unpack$filename),]$taxon = "HOVU"
  
  data_unpack$treatment = NA
  data_unpack$filename %>% str_extract("_.{1,3}C") %>% 
    str_remove("_") %>% str_remove("C") %>% as.numeric() -> data_unpack$treatment
  
  tpc_treatment_data = data_unpack %>% 
    select(curveID, date, taxon, treatment, Photo = A, Tleaf = TleafCnd, gsw, Ci, Qin)
  
  #tpc_treatment_data = subset(tpc_treatment_data, !is.na(Photo))
  
  
  return(tpc_treatment_data)
}

# 
# 
# build_photo_greenhouse = function() {
#   
#   
#   all_files = list.files("data/photo_tpc_greenhouse", recursive = T, full.names = T) %>% 
#     grep(".xlsx", ., value = T, invert = T)
#   data_raw = lapply(all_files, read_6800)
#   
#   data_unpack = c()
#   for (i in 1:length(data_raw)) {
#     data_unpack = bind_rows(data_unpack, data.frame(data_raw[i], curveID = i))
#   }
#   
#   data_unpack$Taxon = NA
#   data_unpack[grep("rasa", data_unpack$filename),]$Taxon = "RASA"
#   data_unpack[grep("boof", data_unpack$filename),]$Taxon = "BOOF"
#   data_unpack[grep("hovu", data_unpack$filename),]$Taxon = "HOVU"
#   
#   tpc_greenhouse_data = data_unpack %>% 
#     select(curveID, date, Taxon, Photo = A, Tleaf = TleafCnd)
#   tpc_greenhouse_data$Treatment = "Greenhouse"
#   
#   
#   return(tpc_greenhouse_data)
# }
# 



build_thermal = function() {
  
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
  
  thermal_means = thermal_data %>% group_by(taxon, treatment) %>% 
    summarize(#Taxon = unique(species),
      #Replicate = unique(replicate),
      #Treatment = unique(treatment),
      #Week = unique(week),
      T_leaf_mean = mean(thermal_mean)-273.15,
      T_leaf_median = mean(thermal_median)-273.15,
      T_air_1 = mean(temp_atm)-273.15,
      T_air_2 = mean(temp_ambient)-273.15,
      T_root_mean = mean(temp_soil_mean)-273.15,
      RH_mean = mean(relative_humidity))
  
  
  thermal_means$taxon = toupper(thermal_means$taxon)
  thermal_means[thermal_means$treatment == 2,]$treatment = 2.5
  thermal_means[thermal_means$treatment == 7,]$treatment = 7.5
  
  return(thermal_means)
}



build_biomass = function() {
  
  
  biomass_filenames = list.files("data/biomass_data", full.names = T)
  biomass = c()
  for(fn in biomass_filenames) biomass = rbind(biomass, read.csv(fn, fileEncoding="UTF-8-BOM"))
  
#  # Build database of biomass data
#  biomass = read.csv("biomass_data/biomass_data_jcg_7_dec_2020.csv")
#  biomass2 = read.csv("biomass_data/biomass_data_jcg_1_apr_2021.csv")
#  biomass = rbind(biomass, biomass2)
  
  # Fix swapped AG/BG
  a = which(biomass$above_below == "AG" & biomass$species == "HOVU" & biomass$replicate == 1 & biomass$temperature_C == 15 & biomass$notes == "Respiration")
  b = which(biomass$above_below == "BG" & biomass$species == "HOVU" & biomass$replicate == 1 & biomass$temperature_C == 15 & biomass$notes == "Respiration")
  biomass[a,]$above_below = "BG"
  biomass[b,]$above_below = "AG"
  
  # This line converts the date and time from text into machine-readable format
  biomass$parsed_time = parse_date_time(paste(as.character(biomass$time), " ", as.character(biomass$date)), "HMdmy")
  
  # Fix bad dates
  
  a = which(biomass$species == "BOOF" & biomass$parsed_time == parse_date_time("2020-11-03 14:50:00", "ymdHMS"))
  biomass[a,]$parsed_time = parse_date_time("2020-11-02 14:50:00", "ymdHMS")
  
  a = which(biomass$species == "HOVU" & biomass$parsed_time == parse_date_time("2021-06-17 10:45:00", "ymdHMS"))
  biomass[a,]$parsed_time = parse_date_time("2021-06-17 15:35:00", "ymdHMS")
  
  a = which(biomass$species == "RASA" & biomass$temperature_C == 5 & biomass$parsed_time == parse_date_time("2021-06-10 11:00:00", "ymdHMS"))
  biomass[a,]$temperature_C = 7.5
  
  a = which(biomass$species == "RASA" & biomass$temperature_C == 7.5 & biomass$parsed_time == parse_date_time("2021-06-10 23:00:00", "ymdHMS"))
  biomass[a,]$parsed_time = parse_date_time("2021-06-10 11:00:00", "ymdHMS")
  
  a = which(biomass$notes == "SWAP 7.5?")
  biomass[a,]$temperature_C = 7.5
  
  a = which(biomass$species == "RASA" & biomass$parsed_time == parse_date_time("2021-06-18 10:45:00", "ymdHMS") & biomass$notes == "")
  biomass[a,]$temperature_C = 5
  biomass[a,]$parsed_time = parse_date_time("2021-05-26 10:15:00", "ymdHMS")
  
  # There is a bonus mass entry here, not sure why - I have two AG entries for this one, so just deleting one
  a = which(biomass$biomass_g == 3.206)
  biomass = biomass[-a,]
  
  
  # Compute number of elapsed days at each measurement - effectively this sets our first harvest time as t = 0 for each species and treatment
  biomass = biomass %>% group_by(species,temperature_C) %>% mutate(elapsed_days = (as.numeric(parsed_time) - min(as.numeric(parsed_time)))/(60*60*24))
  
  biomass_ag = subset(biomass, above_below == "AG") %>% mutate(biomass_ag = biomass_g) %>% select(-date, -time, -above_below, -parsed_time, -biomass_g)
  biomass_bg = subset(biomass, above_below == "BG") %>% mutate(biomass_bg = biomass_g) %>% select(-date, -time, -above_below, -parsed_time, -biomass_g)
  
  biomass_merged = merge(biomass_ag, biomass_bg, by = c("species", "replicate", "temperature_C", "notes", "elapsed_days"))
  
  biomass_merged$prop_aboveground = biomass_merged$biomass_ag / (biomass_merged$biomass_ag + biomass_merged$biomass_bg)
  biomass_merged$biomass_total = biomass_merged$biomass_ag + biomass_merged$biomass_bg
  
  # Now compute average by species, treatment, and time (useful for growth rate)
  biomass = biomass_merged
  biomass_growth = biomass %>% 
    group_by(species, temperature_C, elapsed_days) %>% 
    summarise(biomass_total = mean(biomass_total), 
              prop_aboveground = mean(prop_aboveground),
              biomass_ag = mean(biomass_ag),
              biomass_bg = mean(biomass_bg))
  
  biomass_allocation = biomass_growth %>% 
    group_by(species, temperature_C) %>% 
    summarise(prop_aboveground = mean(prop_aboveground))
  
  # I should also try to get growth rates into that dataframe, then it will have all I need
  # basically what I want here is to fit each treatment and species with NLS, and have it
  # return a dataframe with one line per, including species, treatment, fit parameters, and
  # allocation.
  
  beta = 0.75
  
  growth_rate = biomass %>% 
    group_by(species, temperature_C) %>% 
    do(fit = nls(biomass_total ~ (M_0^(1-beta) + r*elapsed_days*(1-beta))^(1/(1-beta)), data = ., start = list( r = 0.1, M_0 = 0.2)))
  
  growth_rate$r = growth_rate$M_0 = 0
  for (i in 1:dim(growth_rate)[1]) {
    growth_rate[i,]$r = coef(growth_rate$fit[[i]])["r"]
    growth_rate[i,]$M_0 = coef(growth_rate$fit[[i]])["M_0"]
  }
  
  growth_rate = growth_rate %>% select(-fit)
  
  biomass_all = merge(growth_rate, biomass_allocation, by=c("species", "temperature_C"))
  
  colnames(biomass_all)[which(colnames(biomass_all) == "species")] = "taxon"
  colnames(biomass_all)[which(colnames(biomass_all) == "temperature_C")] = "treatment"
  
  return(biomass_all)
}

build_T_tissue_data = function() {
  
  thermal_data = build_thermal() 
  biomass_data = build_biomass() # Done
  
  # Now lets make a dataframe where we match up the thermal and the biomass data
  thermal_biomass = merge(thermal_data, biomass_data, by = c("taxon", "treatment"))
  
  # Compute mean tissue temperature
  thermal_biomass$T_tissue = thermal_biomass$prop_aboveground*thermal_biomass$T_leaf_mean + (1-thermal_biomass$prop_aboveground)*thermal_biomass$T_root_mean
  
  # Compute environmental variables of interest
  thermal_biomass$PPFD = 250
  thermal_biomass = thermal_biomass %>% mutate(T_air = (T_air_1+T_air_2)/2)
  thermal_biomass = thermal_biomass %>% mutate(SVP = (1/1000)*610.7*10^((7.5*T_air)/(237.3+T_air)) )
  thermal_biomass = thermal_biomass %>% mutate(VPD = ((100 - 100*RH_mean)/100)*SVP)
  
  thermal_biomass = thermal_biomass %>% select(taxon, treatment, T_tissue, T_air, PPFD, RH = RH_mean, VPD)
  
}



build_bm_data = function() {
  
  
  biomass_filenames = list.files("data/biomass_data", full.names = T)
  biomass = c()
  for(fn in biomass_filenames) biomass = rbind(biomass, read.csv(fn, fileEncoding="UTF-8-BOM"))
  
  # Fix swapped AG/BG
  a = which(biomass$above_below == "AG" & biomass$species == "HOVU" & biomass$replicate == 1 & biomass$temperature_C == 15 & biomass$notes == "Respiration")
  b = which(biomass$above_below == "BG" & biomass$species == "HOVU" & biomass$replicate == 1 & biomass$temperature_C == 15 & biomass$notes == "Respiration")
  biomass[a,]$above_below = "BG"
  biomass[b,]$above_below = "AG"
  
  # This line converts the date and time from text into machine-readable format
  biomass$parsed_time = parse_date_time(paste(as.character(biomass$time), " ", as.character(biomass$date)), "HMdmy")
  
  # Fix bad dates
  
  a = which(biomass$species == "BOOF" & biomass$parsed_time == parse_date_time("2020-11-03 14:50:00", "ymdHMS"))
  biomass[a,]$parsed_time = parse_date_time("2020-11-02 14:50:00", "ymdHMS")
  
  a = which(biomass$species == "HOVU" & biomass$parsed_time == parse_date_time("2021-06-17 10:45:00", "ymdHMS"))
  biomass[a,]$parsed_time = parse_date_time("2021-06-17 15:35:00", "ymdHMS")
  
  a = which(biomass$species == "RASA" & biomass$temperature_C == 5 & biomass$parsed_time == parse_date_time("2021-06-10 11:00:00", "ymdHMS"))
  biomass[a,]$temperature_C = 7.5
  
  a = which(biomass$species == "RASA" & biomass$temperature_C == 7.5 & biomass$parsed_time == parse_date_time("2021-06-10 23:00:00", "ymdHMS"))
  biomass[a,]$parsed_time = parse_date_time("2021-06-10 11:00:00", "ymdHMS")
  
  a = which(biomass$notes == "SWAP 7.5?")
  biomass[a,]$temperature_C = 7.5
  
  a = which(biomass$species == "RASA" & biomass$parsed_time == parse_date_time("2021-06-18 10:45:00", "ymdHMS") & biomass$notes == "")
  biomass[a,]$temperature_C = 5
  biomass[a,]$parsed_time = parse_date_time("2021-05-26 10:15:00", "ymdHMS")
  
  # There is a bonus mass entry here, not sure why - I have two AG entries for this one, so just deleting one
  a = which(biomass$biomass_g == 3.206)
  biomass = biomass[-a,]
  
  
  # Compute number of elapsed days at each measurement - effectively this sets our first harvest time as t = 0 for each species and treatment
  biomass = biomass %>% 
    group_by(species,temperature_C) %>% 
    mutate(elapsed_days = (as.numeric(parsed_time) - min(as.numeric(parsed_time)))/(60*60*24))
  
  rasa.leaf = read.csv("data/rasa_leafy_biomass_fixed.csv")
  hovu.leaf = read.csv("data/hovu_leafy_biomass.csv")
  
  # Put in rasa and hovu leafy biomass here
  rasa.leaf$parsed_date = rasa.leaf$date %>% gsub("-", ".", .) %>% gsub("2021", "21", .) %>% parse_date("%d.%m.%y")
  hovu.leaf$parsed_date = hovu.leaf$Date %>% parse_date("%d-%b-%y")
  biomass$parsed_date = biomass$parsed_time %>% format("%Y-%m-%d")
  
  rasa.leaf = rasa.leaf %>% select(-date) %>% mutate(species = "RASA")
  biomass = biomass %>% mutate(treatment = temperature_C)
  rasa.leaf = rasa.leaf %>% mutate(replicate = rep)
  
  hovu.leaf = hovu.leaf %>% 
    rename(species = Name, replicate = Number, treatment = Temperature,leafy_biomass_g=Biomass) %>% 
    select(-Biomass.type, -Date,-Time,-Notes,-X)
  hovu.leaf$species = toupper(hovu.leaf$species)
  
  hovu.leaf$parsed_date = as.character(hovu.leaf$parsed_date)
  rasa.leaf$parsed_date = as.character(rasa.leaf$parsed_date)
  biomass$parsed_date = as.character(biomass$parsed_date)
  
  # Bind rows rasa and hovu
  rasa.leaf = rasa.leaf %>% select(-rep)
  leafy_mass = bind_rows(hovu.leaf, rasa.leaf)
  
  #a = merge(biomass, rasa.leaf, all.x = T, by=c("species", "treatment", "replicate", "parsed_date"))
  #b = merge(biomass, rasa.leaf, all.y = T, by=c("species", "treatment", "replicate", "parsed_date"))
  
  biomass = merge(biomass, leafy_mass, all.x = T, by=c("species", "treatment", "replicate", "parsed_date"))
  #b = merge(biomass, leafy_mass, all.y = T, by=c("species", "treatment", "replicate", "parsed_date"))
  

  #biomass[biomass$species=="BOOF",]$leafy_biomass_g = biomass[biomass$species=="BOOF",]$
  #biomass = a
  
  
  biomass_ag = subset(biomass, above_below == "AG") %>% mutate(biomass_ag = biomass_g) %>% select(-date, -time, -above_below, -parsed_time, -biomass_g, -treatment, -parsed_date)
  biomass_bg = subset(biomass, above_below == "BG") %>% mutate(biomass_bg = biomass_g) %>% select(-date, -time, -above_below, -parsed_time, -biomass_g, -treatment, -parsed_date, -leafy_biomass_g)
  
  biomass_merged = merge(biomass_ag, biomass_bg, by = c("species", "replicate", "temperature_C", "notes", "elapsed_days"))
  
  biomass_merged$prop_aboveground = biomass_merged$biomass_ag / (biomass_merged$biomass_ag + biomass_merged$biomass_bg)
  biomass_merged$biomass_total = biomass_merged$biomass_ag + biomass_merged$biomass_bg
  
  # Fill in leaf biomass for BOOF = AG biomass
  biomass_merged[is.na(biomass_merged$leafy_biomass_g),]$leafy_biomass_g = biomass_merged[is.na(biomass_merged$leafy_biomass_g),]$biomass_ag
  biomass_merged$prop_leafy = biomass_merged$leafy_biomass_g/biomass_merged$biomass_total
  
  # Now compute average by species, treatment, and time (useful for growth rate)
  biomass = biomass_merged
 
  return(biomass)
}

filter_photo = function(photo_data, filter_from_raw = F) {
  
  if(filter_from_raw == F) {
    keepers = data.frame( curveID = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67, 69, 72, 73, 76, 78, 79, 80, 81, 83, 84, 85, 89, 90, 93, 94, 95, 96, 97, 98, 101, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 115, 116, 117, 118, 120, 121, 122, 123, 124, 127, 128, 129, 130, 131, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 145, 146, 147, 148, 149, 151, 152, 153, 154, 155, 157))
    photo_data = merge(keepers, photo_data)
  } else {
    set.seed(0)
    photo_cut = photo_data %>% 
      group_by(curveID) %>% 
      mutate(n_before = sum(which(Photo==max(Photo)))) %>% 
      subset(n_before > 3) %>% # Needed for good estimation of E
      subset(treatment != 35) %>% # This is all pretty garbage
      subset(curveID != 156) # This one gives inconsistent results between machines, can't figure out why
    photo_results = fit_curves_michaletz_2021(photo_cut, T_ref = 10, curve_type = "Schoolfield", PLOT=F)
    photo_results = photo_results %>% 
      subset(r_sq > 0.87) %>% 
      subset(E_SE/E < 1)
    
    # OK now let's refit these I guess
    photo_cut2 = merge(photo_cut, data.frame(curveID = photo_results$curveID), by="curveID", PLOT=F) %>% 
      select(-n_before)
    photo_results = fit_curves_michaletz_2021(photo_cut2, T_ref = 10, curve_type = "Schoolfield")
    photo_data = photo_cut2
  }
  
}

