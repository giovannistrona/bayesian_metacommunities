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
            fform=factor(fform, levels=c("convex", "linear", "concave", "sigmoidal")))
  
  ## Extract top species extinction:
  temp <- data %>% filter(web==webs, spf==spfs, scenario %in% scenarioss)
  topextinct <- temp %>% filter(species==topsp, lambda<1) %>%  ## NOTE: max(species) != topsp !!!
      group_by(web, fform, scenario) %>% 
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

  # data %>% filter(web == webs, params %in% paramss, spf==spfs, alpha == a, beta == b, scenario %in% scenarioss) %>%
  data %>% filter(web == webs, spf==spfs, scenario %in% scenarioss) %>%
    group_by(web, fform, scenario, prem, group) %>%
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
    facet_grid(scenario~fform, 
              labeller=label_bquote(.(scenario), italic(pi[i])~"and"~italic(xi[i])~"SGB")) +
    guides(fill=guide_legend(title="Spatial groups", override.aes=list(size=5, alpha=1, fill=cpalette), reverse=FALSE)) +
    theme(text = element_text(size=15), axis.text.x = element_text(angle=75, hjust=1)) + 
    ylab("Metapopulation capacity") 
}
