
require(tidyverse)
require(RColorBrewer)

getPalette  <- colorRampPalette(brewer.pal(11,"RdYlBu"))
path <- "../data/summaries/summaries_gridlike/summary_gridlike_" ## summary files
psizes <- c(0, exp(seq(log(0.01), log(0.5), l=19))) ## perturbation sizes

gl <- tibble() ## gridlike, extracted & summarized
for (i in 1:20) {
  gl <- readRDS(paste0(path, i, ".rds")) %>% ## read ith gridlike summary file
    select(scenario, spinput, prem, alpha, beta, species, TL, lambda) %>%
    filter(scenario!="pre-patch loss") %>% ## drop pre-patch loss scenarios
    mutate(TL=as.integer(round(TL))) %>% ## integer trophic levels
    filter(TL>1) %>% ## only trophic levels 2 and up
    filter(lambda<=1) %>% ## only lambdas below ext threshold
    group_by(scenario, spinput, alpha, beta, species, TL) %>% ## for each combo:
    filter(prem==min(prem)) %>% ## minimum no of patches to remove to cause ext.
    ungroup %>%
    mutate(psize=psizes[i]) %>% ## add size of perturbation of landscape from grid
    arrange(psize, scenario, spinput, alpha, beta, TL) %>% ## reorder rows
    select(psize, scenario, spinput, alpha, beta, TL, prem) %>% ## drop&order cols
    bind_rows(gl) ## add to extracted data table
  cat(paste0(path, i, ".rds processed","\n"))
}
gl <- gl %>% ## remove bloat from species input file names
  mutate(spinput=gsub("../data/species_input/web50350_00", "", spinput)) %>%
  mutate(spinput=gsub(".csv", "", spinput)) %>% ## remove everything except number
  mutate(spinput=case_when( ## mutate numbers into parameter descriptions
    spinput=="3" ~ "pi[i]~TLB~~xi[i]==0.055",
    spinput=="6" ~ "pi[i]~and~xi[i]~TLB",
    spinput=="7" ~ "pi[i]==0.2~~xi[i]~TLB",
    spinput=="8" ~ "pi[i]==0.2~~xi[i]==0.055",
    TRUE ~ spinput)
    ) %>%
  mutate(scenario=gsub(" scenario", "", scenario))

gl %>%
  mutate(TL=factor(TL, levels=seq(7, 2, -1))) %>%
  filter(alpha==1, beta==1) %>%
  group_by(psize, scenario, spinput, alpha, beta, TL) %>%
  summarise(mprem=mean(prem), sprem=sd(prem)) %>%
  ungroup %>%
  ggplot +
  aes(x=psize, y=mprem, ymin=mprem-sprem, ymax=mprem+sprem, colour=TL, fill=TL) +
  geom_line() +
  geom_ribbon(colour=NA, alpha=0.3) +
  facet_grid(scenario~spinput,
             labeller=labeller(.rows=label_value, .cols=label_parsed)) +
  scale_colour_manual(values=colorRampPalette(brewer.pal(11, "RdYlBu"))(7),
                      name="Trophic\nlevel", guide=FALSE) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(11, "RdYlBu"))(7),
                      name="Trophic\nlevel") +
  scale_x_continuous(name="Perturbation to perfect grid structure") +
  scale_y_continuous(name="Patches to remove for eliminating trophic level",
                     limits=c(0, 330)) +
  guides(fill=guide_legend(override.aes=list(alpha=1))) +
  theme_bw()
##ggsave("../SI/gridlike.pdf", width=7, height=4.66)
