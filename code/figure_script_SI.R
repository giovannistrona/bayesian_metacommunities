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

figuredir <- "../SI/"
if(file.exists(figuredir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(figuredir)) #create direcotry for figures, if it does not exist yet
} 

#Figure S1
theme_set(theme_bw())
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
## ggsave("../SI/R_persist.pdf", width=3, height=2.4)

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
## ggsave("../SI/C_persist.pdf", width=3, height=2.4)

#Figure S2
require(deSolve) #ODE solver 
cpal <- c("#999999","#E69F00","#56B4E9","#009E73","#0072B2","#CC79A7","#D55E00")

dat <- tibble(p=c(10^seq(-4.3, -1, by=0.01),
                  seq(0.11, 0.99999, by=0.0001))) %>%
  mutate(exact=0, approximate=log(1/p))
for (i in 1:nrow(dat)) {
  dat$exact[i] <- (integrate(function(r) -1/log(r), 0, 1-dat$p[i],
                             rel.tol=1e-12, abs.tol=1e-12,
                             stop.on.error=FALSE)$value)/(1-dat$p[i])
}

dat %>%
  gather("type", "TL", c("exact", "approximate")) %>%
  ggplot +
  geom_line(aes(x=p, y=TL, colour=type), alpha=0.8) +
  scale_x_continuous(name=expression(paste(pi))) +
  scale_y_continuous(name=expression(paste(tau))) +
  scale_colour_manual(values=c(cpal[2], cpal[5])) +
  theme_bw()
## ggsave("../SI/TL_diffeq.pdf", width=4.2, height=2.7)

#Figure S3
res <- 301

chain <- expand.grid(lambdaM=seq(0, 8, l=res), pi0=seq(0, 1, l=res)) %>%
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

chain2 <- expand.grid(lambdaM=seq(0, 8, l=res), pi0=seq(0, 1, l=res)) %>%
  as_tibble %>%
  filter(pi0!=0) %>%
  mutate(TL=round(lambdaM*log(1/pi0))) %>%
  mutate(TL=ifelse(TL>8, 8, TL))

masterchain <- left_join(chain, chain2, by=c("lambdaM", "pi0")) %>%
  gather(model, TL, c("TL.x", "TL.y")) %>%
  mutate(model=ifelse(model=="TL.x", "iterated", "approximate")) %>%
  mutate(model=factor(model, ordered=TRUE, levels=c("iterated", "approximate")),
         TL=as.character(TL)) %>%
  mutate(TL=ifelse(TL=="8", ">7", TL))

masterchain %>%
  ggplot +
  geom_raster(aes(x=pi0, y=lambdaM, fill=TL), interpolate=TRUE) +
  scale_x_continuous(name=expression(paste(pi)), expand=c(0,0)) +
  scale_y_continuous(name=expression(paste(lambda[M])), expand=c(0,0)) +
  scale_colour_brewer(palette="PuBu") +
  scale_fill_brewer(palette="PuBu", name="Trophic levels") +
  facet_wrap(~model) +
  theme(panel.grid=element_blank())
## ggsave("../SI/TL_chain_approx.pdf", width=6, height=2.7)

#Figures S4, S5 etc. 
infiles <- "../data/summaries/summary_largewebs_*.rds" %>% Sys.glob() #summary files
infile <- "../data/summaries/summary_serengeti.rds"  #summary files

for(infile in infiles){
  dat0 <- infile %>% read_rds 
  webs <- dat0$web %>% unique
  
  for(pweb in webs){
    dat <- dat0 %>% filter(web == pweb)
  
  pas <- unique(dat$alpha)
  pbs <- unique(dat$beta)
  
  for(pa in pas){ #loop over functional forms
    for(pb in pbs){
      
      if(pa==1 && pb==1) ptitles="A"
      if(pa==5 && pb==1) ptitles="B"
      if(pa==1 && pb==5) ptitles="C"
      if(pa==5 && pb==5) ptitles="D"
      
      #Add prem=0 to scenarios
      temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
      temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
      temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
      dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
        filter(scenario != "pre-patch loss")
      
      #Define input for plotting functions
      pscenarios <-  unique(dat_plot$scenario) #scenarios
      pspfs <- "basal" #unique(dat_plot$spf) #which focal species for patch loss ("basal" OR "top" (Figure S14))
      print(pspfs)
      pbasals <- ifelse(pweb=="Serengeti food web with 32 consumer and 129 basal species", 
                        TRUE, FALSE) ## if TRUE basal species are plotted if FALSE only consumer species

      #Metapopulation capacity ~ patches removed
      pparams <- c("pi[i]*' = 0.2,  '*xi[i]*' = 0.055'", "pi[i]*' TLB,  '*xi[i]*' = 0.055'", 
        "pi[i]*' and '*xi[i]*' TLB'", "pi[i]*' = 0.2,  '*xi[i]*' TLB'")
      
      #Create named vectors for facet labels with greek letters
      var_names <- pparams
      names(var_names) <- pparams
      
      #Plot name
      lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>%
                         paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")

      #Generate the plot (Fig. S5, etc.)
      p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)

      #Save it
      ggsave(filename = lambda, width=14, height=8, plot=p)
  
      #Average landscape persistence ~ dispersal distance 
      # pparams <- unique(dat$params)[c(1:5)] #pi TLB
      vec <- which(unique(dat$params) %>% str_detect(pattern="TLB,") == TRUE)
      pparams <- unique(dat$params)[vec]

      #Plot name
      sump <- paste0(figuredir, "sump_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                         paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
      
      #Generate the plot (Fig. S4, etc.)
      p2 <- plot_persistence(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles)

      #Save it
      ggsave(filename = sump, width=14, height=8, plot=p2)
    } 
  }
 }
}

#Figures S14
#Generate .pdf plots for the different model food webs with the top species as focal species for patch loss 
infile <- "../data/summaries/summary_largewebs_1_1.rds" %>% Sys.glob() #summary files
# infile <- "../data/summaries/summary_serengeti.rds" #summary files

  dat0 <- infile %>% read_rds 
  webs <- dat0$web %>% unique
  
  for(pweb in webs){
    dat <- dat0 %>% filter(web == pweb)
    
      pa <- 1 
      pb <- 1 
        
        if(pweb=="Model food web with 350 consumer and 50 basal species") ptitles=paste("A", pweb) 
        if(pweb=="Model food web with 300 consumer and 100 basal species") ptitles=paste("B", pweb) 
        if(pweb=="Model food web with 250 consumer and 150 basal species") ptitles=paste("C", pweb) 
        if(pweb=="Model food web with 200 consumer and 200 basal species") ptitles=paste("D", pweb) 
        if(pweb=="Serengeti food web with 32 consumer and 129 basal species") ptitles=paste("E", pweb) 

        #Add prem=0 to scenarios
        temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
        temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
        temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
        dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
          filter(scenario != "pre-patch loss")
        
        #Define input for plotting functions
        pscenarios <-  unique(dat_plot$scenario) #scenarios
        pspfs <- "top" #unique(dat_plot$spf) ## which focal species for patch loss
        pbasals <- ifelse(pweb=="Serengeti food web with 32 consumer and 129 basal species", TRUE, FALSE) #if TRUE basal species are plotted if FALSE only consumer species
        
        #Metapopulation capacity ~ patches removed
        pparams <- c("pi[i]*' = 0.2,  '*xi[i]*' = 0.055'", "pi[i]*' TLB,  '*xi[i]*' = 0.055'", 
                     "pi[i]*' and '*xi[i]*' TLB'", "pi[i]*' = 0.2,  '*xi[i]*' TLB'")
        
        #Create named vectors for facet labels with greek letters
        var_names <- pparams
        names(var_names) <- pparams
        
        #Plot name
        lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>%
                           paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")

        #Generate the plot
        p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                            a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)

        #Save it
        ggsave(filename = lambda, width=14, height=8, plot=p)
      } 
