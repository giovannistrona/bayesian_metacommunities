## Function to generate the figure meanlambda ~ prem for different loss scenarios
##  for the Serengeti web based on spatial groups (colour coded).
## Input: - data: tibble, - webs: name of web, - paramss: species input (can be a vector), 
##        - scenarioss: habitat loss scenarios (can be a vector), - spfs: focal species for patch loss
##        - a: 1st shape paraemter of the beta distr (Bayesian network), 
##        - b: 2nd shape paremeter of the beta distr (Bayesian network), 
##        - basal: TRUE/FALSE => if TRUE basal species are plotted as well, otherwise only consumer species, 
##        - ptitle: plot title.
## Ouput: - ggplot with number of patches removed on the x-axis and metapopulation capacity on the y-axis, 
##        - rows: habitat loss scenarios, cols: species input parameters (pi_i and xi_i). 

plot_patchloss_groups <- function(data, webs, paramss, scenarioss, spfs, a, b, basals){
  
  data <- data %>% 
    mutate(scenario=factor(scenario, levels=c("best-case scenario", "worst-case scenario", "random scenario")), 
            fform=case_when((alpha==1 & beta==1) == TRUE ~ "linear", 
                            (alpha==5 & beta==5) == TRUE ~ "sigmoidal", 
                            (alpha==1 & beta==5) == TRUE ~ "concave", 
                            (alpha==5 & beta==1) == TRUE ~ "convex"), 
            fform=factor(fform, levels=c("linear", "convex", "concave", "sigmoidal")),
            facety = case_when((scenario=="best-case scenario")==TRUE ~ "best-case", 
                                (scenario=="worst-case scenario")==TRUE ~ "worst-case",
                                (scenario=="random scenario")==TRUE ~ "random"),
            facety = factor(facety, levels=c("best-case", "worst-case", "random")))
           
  ## Extract top species extinction:
  temp <- data %>% filter(web==webs, spf==spfs, facety %in% c("best-case","worst-case","random"))
  topextinct <- temp %>% filter(species==topsp, lambda<1) %>%  ## NOTE: max(species) != topsp !!!
      group_by(web, fform, facety) %>% 
      summarise(minprem=min(prem)) %>% 
      ungroup 


  ## Define theme and colour palette:
  theme_set(theme_bw()) ## set black and white theme
  colourCount <- data$group %>% unique %>% length
  getPalette  <- colorRampPalette(brewer.pal((9),"Greens"))
  greens <- getPalette(10)
  greens <- rev(greens[-c(1:2)]) ## 8 
  getPalette  <- colorRampPalette(brewer.pal((9),"Blues"))
  blues <- getPalette(10)
  blues <- rev(blues[4:7])
  getPalette  <- colorRampPalette(brewer.pal((9),"BrBG"))
  browns <- getPalette(10)
  browns <- browns[1:2]
  
  cpalette  <- c(browns, blues, greens)
  
  if(basals == FALSE) data <- data %>% filter(FT != "basal") 

  # ynames <- c("best-case", "worst-case", "random")
  # names(ynames) <- c("best-case scenario", "worst-case scenario", "random scenario")

  # data %>% filter(web == webs, params %in% paramss, spf==spfs, alpha == a, beta == b, scenario %in% scenarioss) %>%
  data %>% filter(web == webs, spf==spfs, scenario %in% scenarioss) %>%
    group_by(web, fform, facety, prem, group) %>%
    summarise(meanlambda=mean(lambda,na.rm=TRUE), sdlambda=sd(lambda, na.rm=TRUE)) %>% ungroup %>% 
    rename(`Patches removed`=prem) %>%
    ggplot() +
    aes(x=`Patches removed`, y=meanlambda, colour=rev(as.factor(group)), fill=rev(as.factor(group))) +
    scale_color_manual(values=rev(cpalette)) +
    scale_fill_manual(values=rev(cpalette)) +
    geom_hline(yintercept=1, linetype="dashed", alpha=0.3) +
    geom_vline(data = topextinct, aes(xintercept = minprem), linetype="dotdash", alpha=0.7) +
    geom_line(show.legend=FALSE) +
    geom_ribbon(aes(ymin=meanlambda-sdlambda, ymax=meanlambda+sdlambda),colour=NA, alpha=0.3) +
    facet_grid(facety~fform, 
              labeller=label_bquote(.(facety), italic(pi[i])~"and"~italic(xi[i])~"SGB")) +
              # labeller=label_bquote(.(facety), .(fform))) +
    guides(fill=guide_legend(title="Spatial\ngroup", override.aes=list(size=8, alpha=1, fill=cpalette), reverse=FALSE)) +
    theme(text = element_text(size=20), 
          axis.text.x = element_text(angle=75, hjust=1), axis.title = element_text(size=22), 
          strip.text.y = element_text(size=19), strip.text.x = element_text(size=19),
          legend.text = element_text(size=22), legend.title = element_text(size=22)) +
    ylab("Metapopulation capacity")
}
