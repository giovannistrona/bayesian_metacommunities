
#Load libraries 
require(tidyverse)
require(igraph)
require(cowplot)
require(magick)
require(RColorBrewer)

#Load R functions
source("plot_functional_forms.R") #Functional forms 
source("plot_patchloss.R") #Figure 3
source("plot_patchloss_groups.R") #Figure 4

figuredir <- "../figures/"
if(file.exists(figuredir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(figuredir)) #create direcotry for figures, if it does not exist yet
} 

#Figure 1
getPalette  <- colorRampPalette(brewer.pal(11,"RdYlBu"))

#Input file
solfile <- "../data/results/web100300_003_landscape_001_basal_1_1.rds" #result file
sol <- solfile %>% read_rds #read in resultfile
sol <- bind_rows(sol[1,], sol) #add fake point with p=0 to scale alpha such that white equals zero
sol$p[1] <- 0
sol <- sol %>% filter(species %in% c(min(species), max(species))) #filter a basal species and a top species
nsp <- sol %>% pull(species) %>% unique %>% length #number of species to plot
my_palette <- getPalette(nsp)

#Modify data
sol_plot <- sol %>%
   rename(`patches removed`=prem) %>% #rename variables
   mutate(species=case_when(species==min(species) ~ "a basal species",
                            TRUE ~ "a top predator")) %>%
   mutate(scenario=case_when(
     scenario=="best" ~ "A) best-case scenario",
     scenario=="worst" ~ "B) worst-case scenario",
     scenario=="random" ~ "C) random scenario",
     TRUE ~ "all"))

#Generate patch occupancy plots for habitat loss scenarios
plist <- list() #initialize list for plots
for(i in c("A) best-case scenario", "B) worst-case scenario", "C) random scenario")){
  plist[[i]] <- sol_plot %>%
    filter(scenario %in% c("all", i), `patches removed`%in% c(0,100,200,290)) %>%
    mutate(p=ifelse(p<1e-10,0,p)) %>%
    mutate(p=factor(p)) %>%
    ggplot +
    aes(x=x, y=y, alpha=factor(p), colour=species, fill=species) +
    geom_point(shape=21) +
    geom_point(alpha=0.3, shape=1, colour="black") +
    facet_grid(factor(species, levels=c("a basal species", "a top predator"))~`patches removed`,
               labeller=labeller(.rows = label_value, .cols = label_both)) +
    scale_fill_manual(values=my_palette, guide=FALSE) +
    scale_colour_manual(values=my_palette, guide=FALSE) +
    #use discrete here because we converted p into a factor
    scale_alpha_discrete(range=c(0, 1), guide=FALSE) +
    scale_x_continuous(name="", label=abbreviate, limits=c(0, 1)) +
    scale_y_continuous(name="", label=abbreviate, limits=c(0, 1)) +
    labs(title = i) +
    theme(plot.title = element_text(size=15),
          text = element_text(size=15))
}

plist[[1]] #A) best-case scenario
plist[[2]] #B) worst-case scenario
plist[[3]] #C) random scenario

#Figure 2
chain <- expand.grid(lambdaM=seq(0, 5, l=301), pi0=seq(0, 1, l=301)) %>%
  as_tibble %>%
  filter(pi0!=0) %>%
  mutate(lambda=lambdaM/log(1/(1-pi0)), TL=0)
for (i in 1:nrow(chain)) {
  if (chain$lambda[i]<1) next
  while ((chain$lambda[i]>1)&(chain$TL[i]<8)) {
    chain$lambda[i] <- chain$lambdaM[i]/(chain$lambdaM[i]/chain$lambda[i]-log(
      (1-chain$pi0[i])*(1-1/chain$lambda[i])))
    chain$TL[i] <- chain$TL[i] + 1
  }
}
ggplot(chain %>% mutate(TL=factor(TL))) +
   geom_raster(aes(x=pi0, y=lambdaM, fill=TL), interpolate=TRUE) +
   scale_x_continuous(name=expression(paste(pi)), expand=c(0,0)) +
   scale_y_continuous(name=expression(paste(lambda[M])), expand=c(0,0)) +
   scale_colour_brewer(palette="PuBu") +
   scale_fill_brewer(palette="PuBu", name="Trophic\nlevels") +
theme(panel.grid=element_blank())
# ggsave("../figures/TL_chain.pdf", width=4, height=3)

#Figure 3
#Generate plots for the different consumer responses to resource loss (functional forms).
for(a in c(1,5))
  for(b in c(1,5))
    plot_functional_form(figuredir = "../figures/", alpha = a, beta = b)

infiles <- "../data/summaries/summary_webs_*.rds" %>% Sys.glob() #summary files

for(infile in infiles){
  dat <- infile %>% read_rds %>% filter(web == "Model food web with 300 consumer and 100 basal species")

  pas <- unique(dat$alpha)
  pbs <- unique(dat$beta)
  
  for(pa in pas){ #loop over functional forms
    for(pb in pbs){
      
      if(pa==1 && pb==1) ptitles="A"
      if(pa==5 && pb==1) ptitles="B"
      if(pa==1 && pb==5) ptitles="C"
      if(pa==5 && pb==5) ptitles="D"
      
      # plot_functional_form(figuredir, pa, pb)
      
      #Add prem=0 to scenarios
      temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
      temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
      temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
      dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
        filter(scenario != "pre-patch loss")
      
      #Define input for plotting functions
      pweb <- unique(dat_plot$web) #web
      pscenarios <-  unique(dat_plot$scenario) #scenarios
      pspfs <- "basal" #unique(dat_plot$spf) ## which focal species for patch loss
      pbasals <- FALSE ## if TRUE basal species are plotted if FALSE only consumer species
      
      #Metapopulation capacity ~ patches removed
      pparams <- c("pi[i]*' = 0.2,  '*xi[i]*' = 0.055'", "pi[i]*' TLB,  '*xi[i]*' = 0.055'", 
                   "pi[i]*' and '*xi[i]*' TLB'", "pi[i]*' = 0.2,  '*xi[i]*' TLB'") ##selected pi_i and xi_i combinations

      #Create named vectors for facet labels with greek letters
      var_names <- pparams
      names(var_names) <- pparams

      #Plot name
      lambda <- paste0("../figures/", "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                         paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
      
      #Generate the plot
      p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)
      
      #Save it
      # ggsave(filename = lambda, width=14, height=8, plot=p)
    }
  }
}


#Figure 4
#If not yet done, generate plots for the different consumer responses to resource loss (functional forms).

#patch removal based on basal species (ACASEN	Acacia senegal, sole member of SG 12):
infile <- "../data/summaries/summary_serengeti.rds" #summary file
#patch removal based on antoher basal species (BOSAUG Boscia augustifolia, SG 14):
# infile <- "../data/summaries/summary_serengeti_BOSAUG.rds" #summary file
dat <- infile %>% read_rds %>% filter(spinput == "../data/species_input/serengeti_009.csv")

pas <- unique(dat$alpha) 
pbs <- unique(dat$beta)

for(pa in pas){ #loop over functional forms
  for(pb in pbs){
  
    temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
    temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
    temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
  
    temp4 <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3)
    if(pa==1 && pb==1) dat_plot <- temp4
      else dat_plot <- bind_rows(dat_plot, temp4)
      }
    }
    dat_plot <- dat_plot %>% filter(scenario != "pre-patch loss")

    ## Define input
    pweb <- unique(dat_plot$web)
    pscenarios <- c("best-case scenario", "worst-case scenario", "random scenario")
    pspfs <- "basal"#unique(dat_plot$spf)
    pbasals <- TRUE
    
    #Metapopulation capacity ~ patches removed
    pparams <- unique(dat_plot$params) ## pi and xi

    #Plot name  
    lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                       paste(.,collapse="_") %>% unname, "_", pspfs, "_SGB.pdf")
                       # paste(.,collapse="_") %>% unname, "_", pspfs, "_SGB_BOSAUG.pdf")
    
    #Generate the plot
    p <- plot_patchloss_groups(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                             a = pa, b = pb, basals = pbasals)
    
    #Save it
    # ggsave(filename = lambda, width=14, height=8, plot=p)

