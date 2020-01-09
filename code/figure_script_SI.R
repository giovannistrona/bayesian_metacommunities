#Load libraries 
require(tidyverse) #for data handling
require(cowplot) #for inlay plots
require(magick) #for inlay plots
require(RColorBrewer) #for mixing palettes
theme_set(theme_bw()) #set black and white theme

#Load R functions
source("plot_functional_forms.R") #generate .png plots for the different consumer responses to resource loss (functional forms) 
source("plot_persistence.R")
source("plot_patchloss.R")

figuredir <- "../figures/"
if(file.exists(figuredir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(figuredir)) #create direcotry for figures, if it does not exist yet
} 

#Figure S1
dpalette <- brewer.pal(9, "Blues")[c(4, 5, 6, 7, 8)]

expand.grid(piR=c(seq(0, 0.99, l=5001), seq(0.99, 1, l=15001)),
            lR=c(0.5, 1, 2, 3, 5)) %>%
  as_tibble %>%
  mutate(lMR=-lR*log(1-piR)) %>%
  filter(lMR>=0, lMR<=5) %>%
  ggplot +
  aes(x=piR, y=lMR, colour=factor(lR), size=factor(lR), linetype=factor(lR)) +
  geom_path() +
  scale_x_continuous(name=expression(paste(pi[R])), labels=abbreviate) +
  scale_y_continuous(name=expression(paste(lambda[M[R]]))) +
  scale_colour_manual(name=expression(paste(lambda[R])), values=dpalette) +
  scale_size_manual(name=expression(paste(lambda[R])),
                    values=c(0.5, 0.6, 0.5, 0.5, 0.5)) +
  scale_linetype_manual(name=expression(paste(lambda[R])),
                        values=c(2, 1, 2, 2, 2))
## ggsave("../figures/R_persist.pdf", width=3, height=2.4)

expand.grid(s=seq(0, 1, l=20002), lMC=c(0.5, 1, 2, 3, 4)) %>%
  as_tibble %>%
  mutate(lR=s/(s-exp(-lMC))) %>%
  ggplot +
  aes(x=s, y=lR, colour=factor(lMC)) +
  geom_path(size=0.5) +
  scale_x_continuous(name=expression(paste((1-pi[C])(1-pi[R]))),
                     labels=abbreviate) +
  scale_y_continuous(name=expression(paste(lambda[R])), limits=c(1, 5)) +
  scale_colour_manual(values=dpalette, name=expression(paste(lambda[M[C]])))
## ggsave("../figures/C_persist.pdf", width=3, height=2.4)

#Figure S2,..., Figure S8 
#Generate .png plots for the different consumer responses to resource loss (functional forms) 
plot_functional_forms(figuredir) 

infiles <- "../data/summaries/summary_largewebs_*.rds" %>% Sys.glob() #summary files

for(infile in infiles){
  dat0 <- infile %>% read_rds
  pa <- unique(dat$alpha)
  pb <- unique(dat$beta)
  webnames <- unique(dat$web)
  
  for(webname in webnames){
    
    dat <- dat0 %>% filter(web == webname)
      
    if(pa==5 && pb==1) ptitles="A"
    if(pa==5 && pb==5) ptitles="B"
    if(pa==1 && pb==1) ptitles="C"
    if(pa==1 && pb==5) ptitles="D"
        
    fform_png <-  paste0(figuredir,"fform_",pa,"_",pb,".png")
        
    #Add prem=0 to scenarios
    temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
    temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
    temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
    dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
    filter(scenario != "pre-patch loss")
      
    #Define input for plotting functions
    pweb <- unique(dat_plot$web) #web
    pscenarios <-  unique(dat_plot$scenario) #scenarios
    pspfs <- "basal" #which focal species for patch loss (basal/top)
    pbasals <- FALSE #if TRUE basal species are plotted if FALSE only consumer species
        
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

    #Average landscape persistence ~ dispersal distance 
    pparams <- unique(dat$params)[c(1:5)] #pi TLB
  
    #Plot name
    sump <- paste0(figuredir, "sump_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                    paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".png")
        
    #Generate the plot
    p <- plot_persistence(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                        a = pa, b = pb, basals = pbasals, ptitle=ptitles) 
        
    #Add functional form inset
    g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
    g
        
    #Save it
    # ggsave(filename = sump, width=12, height=8, plot=g, device="png")
  } 
}


#Figure S9 & Figure S10
#Generate .png plots for the different consumer responses to resource loss (functional forms) 
plot_functional_forms(figuredir) 

infile <- "../data/summaries/summary_serengeti.rds" #summary files

dat <- infile %>% read_rds  
pas <- unique(dat$alpha)
pbs <- unique(dat$beta)
    
for(pa in pas){ #loop over functional forms
  for(pb in pbs){
    if(pa==5 && pb==1) ptitles="A"
    if(pa==5 && pb==5) ptitles="B"
    if(pa==1 && pb==1) ptitles="C"
    if(pa==1 && pb==5) ptitles="D"
    
    fform_png <-  paste0(figuredir,"fform_",pa,"_",pb,".png")
    
    #Add prem=0 to scenarios
    temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
    temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
    temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
    dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
      filter(scenario != "pre-patch loss")
    
    #Define input for plotting functions
    pweb <- unique(dat_plot$web) #web
    pscenarios <-  unique(dat_plot$scenario) #scenarios
    pspfs <- "basal" #which focal species for patch loss (basal/top)
    pbasals <- TRUE #if TRUE basal species are plotted if FALSE only consumer species
    
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
    # g
    
    #Save it
    ggsave(filename = lambda, width=12, height=8, plot=g)
    
    #Average landscape persistence ~ dispersal distance 
    pparams <- unique(dat$params)[c(1:5)] #pi TLB
    
    #Plot name
    sump <- paste0(figuredir, "sump_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                     paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".png")
    
    #Generate the plot
    p <- plot_persistence(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles) 
    
    #Add functional form inset
    g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
    # g
    
    #Save it
    ggsave(filename = sump, width=12, height=8, plot=g, device="png")
  } 
}

#Figure S11 Model food webs
#Generate .png plots for the different consumer responses to resource loss (functional forms) 
plot_functional_forms(figuredir) 

infile <- c("../data/summaries/summary_largewebs_1_1.rds") #summary files

dat <- infile %>% read_rds
pa <- unique(dat$alpha)
pb <- unique(dat$beta)
webnames <- dat$web %>% unique

for(webname in webnames){
  ptitles=webname
  fform_png <-  paste0(figuredir,"fform_",pa,"_",pb,".png")
    
  #Add prem=0 to scenarios
  temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
  temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
  temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
  dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
  filter(scenario != "pre-patch loss")
    
  #Define input for plotting functions
  pweb <- unique(dat_plot$web) #web
  pscenarios <-  unique(dat_plot$scenario) #scenarios
  pspfs <- "top" #which focal species for patch loss (basal/top)
  pbasals <- FALSE #if TRUE basal species are plotted if FALSE only consumer species
  
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


#Figure S11 Serengeti food web
#Generate .png plots for the different consumer responses to resource loss (functional forms) 
plot_functional_forms(figuredir) 

infile <- "../data/summaries/summary_serengeti.rds" #summary files

pa <- 1
pb <- 1
dat <- infile %>% read_rds %>% filter(alpha == pa, beta == pb) 

ptitles <- dat %>% pull(web) %>% unique
fform_png <-  paste0(figuredir,"fform_",pa,"_",pb,".png")
    
#Add prem=0 to scenarios
temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
    filter(scenario != "pre-patch loss")

#Define input for plotting functions
pweb <- unique(dat_plot$web) #web
pscenarios <-  unique(dat_plot$scenario) #scenarios
pspfs <- "top" #which focal species for patch loss (basal/top)
pbasals <- TRUE #if TRUE basal species are plotted if FALSE only consumer species
    
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
    