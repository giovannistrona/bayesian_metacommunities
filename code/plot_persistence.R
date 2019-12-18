## Function to generate the figure sump/300 ~ xi for different loss scenarios (trophic levels are colour coded).
## Input: - data: tibble, - webs: name of web, - paramss: species input (must be a vector),  
##        - scenarioss: habitat loss scenarios (can be a vector), - spfs: focal species for patch loss
##        - a: 1st shape paraemter of the beta distr (Bayesian network), 
##        - b: 2nd shape paremeter of the beta distr (Bayesian network), 
##        - basal: TRUE/FALSE => if TRUE basal species are plotted as well, otherwise only consumer species, 
##        - ptitle: plot title.
## Ouput: - ggplot with range of dispersal distances on the x-axis and average landscape persistence on the y-axis 
##        - rows: habitat loss scenarios, cols: number of patches removed.  

plot_persistence <- function(data, webs, paramss, scenarioss, spfs, a, b, basals, ptitle){
  
  data <- data %>% 
    mutate(scenario=factor(scenario, levels=c("best-case scenario", "worst-case scenario", "random scenario")))
  
  ## Define theme and colour palette:
  theme_set(theme_bw()) ## set black and white theme
  colourCount <- data$TL05 %>% unique %>% length
  getPalette  <- colorRampPalette(brewer.pal((11),"RdYlBu"))
  cpalette <- rev(getPalette(colourCount+1))
  cpalette <- cpalette[-6]
  
  if(basals == FALSE) data <- data %>% filter(FT != "basal") 
  
  data %>% filter(web == webs, params %in% paramss, spf==spfs, alpha == a, beta == b, scenario %in% scenarioss, 
                  prem %in% c(0, 100, 200, 290)) %>%
    rename(`patches removed`=prem) %>%
    ggplot() +
    aes(x=as.factor(xi), y=(sump/300), colour=as.factor(TL05), fill=as.factor(TL05)) +
    geom_point(alpha=1, size=5, shape=1) +
    expand_limits(y = c(0,1)) + 
    facet_grid(scenario~`patches removed`, 
               labeller=labeller(.rows = label_value, .cols = label_both), scales="free") + ## reverse order remaining patches
    scale_color_manual(values=cpalette, guide = FALSE) +
    scale_fill_manual(values=cpalette, guide = FALSE) +
    guides(colour=guide_legend(title="Trophic levels", 
                               override.aes=list(size=5, alpha=0.7, shape=15), 
                               reverse=TRUE)) +
    ggtitle(ptitle) +
    theme(text = element_text(size=15), axis.text.x = element_text(angle=75, hjust=1)) + 
    ylab("Average landscape persistence") + 
    xlab("Dispersal distance")
}
