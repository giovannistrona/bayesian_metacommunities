
#Load libraries 
require(tidyverse)
require(igraph)
require(cowplot)
require(magick)
require(RColorBrewer)

#Load R functions
source("functions.R")
source("plot_functional_forms.R")
source("plot_linkloss.R")

figuredir <- "../SI/"
infile <- "../data/summaries/summary_linkrem_basal.rds"
dat <- infile %>% read_rds %>% filter(web == "Model food web with 300 consumer and 100 basal species")
pas <- unique(dat$alpha)
pbs <- unique(dat$beta)

for(pa in pas) for(pb in pbs){ { #loop over functional forms
  
  if(pa==1 && pb==1) ptitles="A"
  if(pa==5 && pb==1) ptitles="B"
  if(pa==1 && pb==5) ptitles="C"
  if(pa==5 && pb==5) ptitles="D"
  
  # plot_functional_form(figuredir, pa, pb)
  
  #Add lrem=0 to scenarios
  temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
  temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
  temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
  dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
    filter(scenario != "pre-patch loss") %>%
    mutate(lrem=lrem/choose(300,2))
  
  #Define input for plotting functions
  pweb <- unique(dat_plot$web) #web
  pscenarios <-  unique(dat_plot$scenario) #scenarios
  pspfs <- unique(dat_plot$spf) ## which focal species for patch loss
  print(pspfs)
  pbasals <- FALSE ## if TRUE basal species are plotted if FALSE only consumer species
  
  #Metapopulation capacity ~ patches removed
  # pparams <- unique(dat_plot$params)[c(3,6:8)] ## selected pi_i and xi_i combinations
  pparams <- unique(dat_plot$params) ## selected pi_i and xi_i combinations
  
  #Create named vectors for facet labels with greek letters
  var_names <- pparams
  names(var_names) <- pparams
  
  #Plot name
  lambda <- paste0(figuredir, "linkrem_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                     paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,"_sym.pdf")
  
  #Generate the plot
  p <- plot_linkloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                      a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)
  
  #Save it
  ggsave(filename=lambda, width=12, height=8, plot=p)
}}
