# Scripts to reproduce analysis and figures for "Acclimation 
# unifies the scaling of carbon assimilation across climate 
# gradients and levels of organization".
#
# Josef Garen, last modified 10 June 2024

library(osfr)
library(cowplot)
library(gridExtra)
library(nlme)
library(lubridate)
library(tidyverse)
library(rTPC)
library(nls.multstart)
library(nlme)
library(nlstools)
library(rgl)
library(plot3D)
library(plantecophys)
library(R.utils)
library(rcartocolor)
library(RColorBrewer)

source("code/build_database.R")
source("code/stat_functions.R")
source("code/read_6800.R")
source("code/get_physio.R")
source("code/acclimation_models.R")
source("code/fit_curves_all.R")

# Download optimal Vcmax code
if(!dir.exists("data")) {
  download.file("https://github.com/SmithEcophysLab/optimal_vcmax_R/archive/refs/heads/master.zip",
              destfile = "optimal_vcmax_R-master.zip")
  unzip("optimal_vcmax_R-master.zip")
  file.remove("optimal_vcmax_R-master.zip")
}

sourceDirectory('optimal_vcmax_R-master/functions', modifiedOnly = FALSE)
source("optimal_vcmax_R-master/calc_optimal_vcmax.R")

# Load figure code
source("code/main_figures/make_main_fig_2.R")
source("code/main_figures/make_main_fig_3.R")
source("code/main_figures/make_main_fig_4.R")
source("code/main_figures/make_main_fig_5.R")

# Load supplementary figure code
source("code/supplementary_figures/make_supp_fig_1.R")
source("code/supplementary_figures/make_supp_fig_3.R")
source("code/supplementary_figures/make_supp_fig_4.R")
source("code/supplementary_figures/make_supp_fig_5.R")
source("code/supplementary_figures/make_supp_fig_6.R")
source("code/supplementary_figures/make_supp_fig_7.R")


# Set up themes for plots
my_theme = theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

theme_transparent =   theme(
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent", color=NA) # get rid of legend panel bg
)

# Download raw data from OSF
if(!dir.exists("data")) {
  osf_retrieve_node("https://osf.io/8rp5v/") %>% osf_ls_files() %>% osf_download()
  unzip("data.zip")
  file.remove("data.zip")
}

# Loads all data to global environment.
# First time running the code set from_raw = T; this will build the database
#   from the raw data files, fit models, etc. This will take a long time to run
# After the first successful run, you can set from_raw = F. This will load
#   saved models etc. and save substantial time
load.all.data(from_raw = F)

dir.create("figures")

# Main figures
# Fig. 1 is conceptual diagram
make_main_fig_2() # AT and GT data and model fits
make_main_fig_3() # Acclimation of physiology and morphology
make_main_fig_4() # Barchart comparison 
make_main_fig_5() # Linearity of Gm and An

# Supplementary figures
make_supp_fig_1() # Isothermality plot
# Fig. S2 is photos only
make_supp_fig_3() # Acclimation simulation data
make_supp_fig_4() # Biomass-temperature facet plot
make_supp_fig_5() # Biomass-temperature 3D model plots
make_supp_fig_6() # Respiration activation energy bar chart
make_supp_fig_7() # Respiration-Growth linearity test

### END ###