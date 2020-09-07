
require(tidyverse)

landscs <- Sys.glob("../../data/landscapes/landscape*.csv")[1]
spinputs <- Sys.glob("../../data/linkrem_input/web*.csv")
spfs <- "basal"
outputdir <- "../../data/linkrem_results/"
alpha <- c(1, 5)
beta <- c(1, 5)

tab <- expand_grid(landsc=landscs, spinput=spinputs, spf=spfs,
                   outputdir=outputdir, alpha=alpha, beta=beta) %>%
  mutate(web=paste0("../../data/webs/",
                    str_split_fixed(basename(spinput), "_", 2)[,1], ".csv")) %>%
  dplyr::select(web, landsc, spinput, spf, outputdir, alpha, beta)

for(r in 1:nrow(tab))
  source("run_instance.R")

