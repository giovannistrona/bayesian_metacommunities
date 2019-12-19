
#Load libraries 
require(tidyverse)
require(cowplot)
require(magick)
require(RColorBrewer)

#Load R functions
source("plot_functional_forms.R") #Functional forms 
source("plot_patchloss.R") #Figure 2
source("plot_patchloss_groups.R") #Figure 3

figuredir <- "../figures/"
if(file.exists(figuredir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(figuredir)) #create direcotry for figures, if it does not exist yet
} 

#Figure 1
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
  scale_fill_brewer(palette="PuBu", name="Trophic levels") +
  theme(panel.grid=element_blank())
## ggsave("../figures/TL_chain.pdf", width=4, height=2.7)

#Figure 2
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


#Figure 3
#Generate .png plots for the different consumer responses to resource loss (functional forms) 
plot_functional_forms(figuredir) 

infiles <- "../data/summaries/summary_largewebs_*.rds" %>% Sys.glob() #summary files

for(infile in infiles){
  dat <- infile %>% read_rds %>% filter(web == "Model food web with 300 consumer and 100 basal species") 
  pas <- unique(dat$alpha)
  pbs <- unique(dat$beta)
  
  for(pa in pas){ #loop over functional forms
    for(pb in pbs){
      
      if(pa==5 && pb==1) ptitles="A"
      if(pa==5 && pb==5) ptitles="B"
      if(pa==1 && pb==1) ptitles="C"
      if(pa==1 && pb==5) ptitles="D"
      
      fform_png <-  paste0("/homes/jh57masa/github/bayesian_metacommunities/figures/fform_",pa,"_",pb,".png")
      
      #Add prem=0 to scenarios
      temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
      temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
      temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
      dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
        filter(scenario != "pre-patch loss")
      
      #Define input for plotting functions
      pweb <- unique(dat_plot$web) #web
      pscenarios <-  unique(dat_plot$scenario) #scenarios
      pspfs <- "basal" ## which focal species for patch loss
      pbasals <- FALSE ## if TRUE basal species are plotted if FALSE only consumer species
      
      #Metapopulation capacity ~ patches removed
      pparams <- unique(dat_plot$params)[c(3,6:8)] ## selected pi_i and xi_i combinations

      #Create named vectors for facet labels with greek letters
      var_names <- pparams
      names(var_names) <- pparams
      
      #Plot name
      lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                         paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
      
      #Generate the plot
      p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)
      
      #Add functional form inset
      g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
      g
      
      #Save it
      # ggsave(filename = lambda, width=12, height=8, plot=g)
    } 
  }
}


#Figure 4
# #If not yet done, generate .png plots for the different consumer responses to resource loss (functional forms) 
# plot_functional_forms(figuredir) 

infile <- "../data/summaries/summary_serengeti.rds" #summary file
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
    pspfs <- "basal"
    pbasals <- TRUE
    
    #Metapopulation capacity ~ patches removed
    pparams <- unique(dat_plot$params) ## pi and xi

    #Plot name  
    lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                       paste(.,collapse="_") %>% unname, "_", pspfs, "_SGB.pdf")
    
    #Generate the plot
    p <- plot_patchloss_groups(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                             a = pa, b = pb, basals = pbasals)
    
    #Combine functional forms 
    fforms <-  ggdraw() + draw_image(paste0(figuredir, "fform_5_1.png"), x = -0.36, y = 1.01, hjust = 0, vjust = 1, width = 1, height = 1) +
               draw_image(paste0(figuredir, "fform_1_1.png"), x = -0.17, y = 1.01, hjust = 0, vjust = 1, width = 1, height = 1) + 
               draw_image(paste0(figuredir, "fform_1_5.png"), x = 0.03, y = 1.01, hjust = 0, vjust = 1, width = 1, height = 1) +
              draw_image(paste0(figuredir, "fform_5_5.png"), x = 0.22, y = 1.01, hjust = 0, vjust = 1, width = 1, height = 1)
    
    #Add them as insets
    g <- plot_grid(fforms, p, nrow=2, rel_heights = c(0.8,3))
    g
    
    #Save it
    # ggsave(filename = lambda, width=11, height=8, plot=g) 


