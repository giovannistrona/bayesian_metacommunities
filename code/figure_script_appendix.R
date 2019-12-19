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

figuredir <- "../appendix/"
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


#Figure S2
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


#Figure S4 
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
      
      fform_png <-  paste0(figuredir,"fform_",pa,"_",pb,".png")
      
      #Add prem=0 to scenarios
      temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
      temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
      temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
      dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3)
      
      #Define input for plotting functions
      pweb <- unique(dat_plot$web) ## which web
      pscenarios <- c("best-case scenario", "worst-case scenario", "random scenario") ## which scenarios
      pspfs <- "basal" ## which focal species for patch loss
      pbasals <- FALSE ## if TRUE basal s pecies are plotted if FALSE only consumers
      
      #Average landscape persistence ~ dispersal distance 
      pparams <- unique(dat$params)[c(1:5)] #pi TLB
      
      #Plot name
      sump <- paste0(figuredir, "sump_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                       paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
      
      #Generate the plot
      p <- plot_persistence(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                            a = pa, b = pb, basals = pbasals, ptitle=ptitles) 
      
      #Add functional form inset
      g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
      g
      
      #Save it
      # ggsave(filename = sump, width=12, height=8, plot=g)
      
    } 
  }
}

#Figure S3
#Generate .png plots for the different consumer responses to resource loss (functional forms) 
plot_functional_forms(figuredir) 

qs = c("Model food web with 350 consumer and 50 basal species", 
        "Model food web with 250 consumer and 150 basal species", 
        "Model food web with 200 consumer and 200 basal species") 
           
infiles <- "../data/summaries/summary_largewebs_*.rds" %>% Sys.glob() #summary files

for(q in qs){
for(infile in infiles){
  dat <- infile %>% read_rds %>% filter(web == q) 
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
      pspfs <- "basal" #which focal species for patch loss
      pbasals <- TRUE #if TRUE basal s pecies are plotted if FALSE only consumers
      
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
      ggsave(filename = lambda, width=12, height=8, plot=g)
    } 
  }
}
}

#spfs = "top"
qs = c("Model food web with 350 consumer and 50 basal species",       
       "Model food web with 300 consumer and 100 basal species",
       "Model food web with 250 consumer and 150 basal species", 
       "Model food web with 200 consumer and 200 basal species") 

infile <- "../data/summaries/summary_largewebs_1_1.rds" %>% Sys.glob() #summary files
dat0 <- infile %>% read_rds

for(q in qs){
    dat <- dat0 %>% filter(web == q) 
    pa <- 1
    pb <- 1
    
    
    if(q == qs[1]) ptitles="A"
    if(q == qs[2]) ptitles="B"
    if(q == qs[3]) ptitles="C"
    if(q == qs[4]) ptitles="D"

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
        pspfs <- "top" #which focal species for patch loss
        pbasals <- TRUE #if TRUE basal s pecies are plotted if FALSE only consumers
        
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
        ggsave(filename = lambda, width=12, height=8, plot=g)
}




# 
# #40-species webs
# infile <- "../data/summaries/summary_webs.rds" #summary file
# dat0 <- infile %>% read_rds 
# qs <- dat0$web %>% unique
# 
# for(q in qs){
#     dat <- dat0 %>% filter(web==q)
#     pas <- unique(dat$alpha)
#     pbs <- unique(dat$beta)
#     
#     for(pa in pas){ #loop over functional forms
#       for(pb in pbs){
#         
#         if(pa==5 && pb==1) ptitles="A"
#         if(pa==5 && pb==5) ptitles="B"
#         if(pa==1 && pb==1) ptitles="C"
#         if(pa==1 && pb==5) ptitles="D"
#         
#         fform_png <-  paste0(figuredir,"fform_",pa,"_",pb,".png")
#         
#         #Add prem=0 to scenarios
#         temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
#         temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
#         temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
#         dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
#           filter(scenario != "pre-patch loss")
#         
#         #Define input for plotting functions
#         pweb <- unique(dat_plot$web) #web
#         pscenarios <-  unique(dat_plot$scenario) #scenarios
#         pspfs <- "basal" #which focal species for patch loss
#         pbasals <- TRUE #if TRUE basal s pecies are plotted if FALSE only consumers
#         
#         #Metapopulation capacity ~ patches removed
#         pparams <- unique(dat_plot$params)[c(3,6:8)] ## selected pi_i and xi_i combinations
#         
#         #Create named vectors for facet labels with greek letters
#         var_names <- pparams
#         names(var_names) <- pparams
#         
#         #Plot name
#         lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
#                            paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
#         
#         #Generate the plot
#         p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
#                             a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)
#         
#         #Add functional form inset
#         g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
#         g
#         
#         #Save it
#         ggsave(filename = lambda, width=12, height=8, plot=g)
#       } 
#     }
#   }
# 
# #top 40-species webs
# for(q in qs){
#   dat <- dat0 %>% filter(web==q)
#   pa <- 1
#   pb <- 1
# 
#       fform_png <-  paste0(figuredir,"fform_",pa,"_",pb,".png")
#       
#       #Add prem=0 to scenarios
#       temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
#       temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
#       temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
#       dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
#         filter(scenario != "pre-patch loss")
#       
#       #Define input for plotting functions
#       pweb <- unique(dat_plot$web) #web
#       pscenarios <-  unique(dat_plot$scenario) #scenarios
#       pspfs <- "top" #which focal species for patch loss
#       pbasals <- TRUE #if TRUE basal s pecies are plotted if FALSE only consumers
#       
#       #Metapopulation capacity ~ patches removed
#       pparams <- unique(dat_plot$params)[c(3,6:8)] ## selected pi_i and xi_i combinations
#       
#       #Create named vectors for facet labels with greek letters
#       var_names <- pparams
#       names(var_names) <- pparams
#       
#       #Plot name
#       lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
#                          paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
#       
#       #Generate the plot
#       p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
#                           a = pa, b = pb, basals = pbasals, ptitle="", labnames = var_names)
#       
#       #Add functional form inset
#       g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
#       g
#       
#       #Save it
#       ggsave(filename = lambda, width=12, height=8, plot=g)
# }


#Serengeti
infile <- "../data/summaries/summary_serengeti.rds" #summary file
dat <- infile %>% read_rds 

  # pas <- unique(dat$alpha)
  # pbs <- unique(dat$beta)
  
  # for(pa in pas){ #loop over functional forms
  #   for(pb in pbs){
  #     # 
      # if(pa==5 && pb==1) ptitles="A"
      # if(pa==5 && pb==5) ptitles="B"
      # if(pa==1 && pb==1) ptitles="C"
      # if(pa==1 && pb==5) ptitles="D"
ptitles=dat$web %>% unique
      pa <- 1
      pb <- 1
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
      pspfs <- "top" #which focal species for patch loss
      pbasals <- TRUE #if TRUE basal s pecies are plotted if FALSE only consumers
      
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
      ggsave(filename = lambda, width=12, height=8, plot=g)
  #   } 
  # }




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




