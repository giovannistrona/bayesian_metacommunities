#Load libraries
require(tidyverse)
require(igraph)

#Load R fucntions
source("create_adj_matrix.R")

#Helper functions
my_round_any = function(x, accuracy, f=round){
  f(x/ accuracy) * accuracy
}

numextract <- function(string){ 
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
}

## Summarize over all patches for each simulation run and modify the generated summary. 
## Input: - infiles: simulation output files (w/ path & extension)
## Output: -out: summarized output
generate_summary <- function(infiles){
    i <- 1
    for (infile in infiles) {
      print(infile)
      sol <- infile %>% readRDS
      problems(sol)
      #stop_for_problems(sol)
      temp <- sol %>%
        group_by(web, landsc, spinput, seed, spf, scenario, prem, alpha, beta,
                 species, TL, OI, FT) %>%
        summarise(pi=unique(pi), xi=unique(xi), 
                  sump=sum(p), meanp=mean(p), minp=min(p), maxp=max(p),
                  lambda=mean(lambda), meand=mean(delta), mind=min(delta),
                  maxd=max(delta)) %>%
        ungroup
      if (i==1) dat <- temp else dat <- rbind(dat, temp)
      i <- i+1
    }
    
    #Modify summary
    out <- dat %>% mutate(
      topsp = case_when(
        web != ("../data/webs/serengeti.csv") ~ read_csv("../data/webs/web100300.csv") %>% 
        create_adj_matrix %>% rownames %>% tail(n = 1),
      web == ("../data/webs/serengeti.csv") ~ read_csv("../data/webs/serengeti.csv") %>% 
        create_adj_matrix %>% rownames %>% tail(n = 1)), 
      web = case_when(
        numextract(web) == "100300" ~ "Model food web with 300 consumer and 100 basal species",
        numextract(web) == "50350" ~ "Model food web with 350 consumer and 50 basal species",
        numextract(web) == "150250" ~ "Model food web with 250 consumer and 150 basal species",
        numextract(web) == "200200" ~ "Model food web with 200 consumer and 200 basal species",
        web == ("../data/webs/serengeti.csv") ~ "Serengeti food web with 32 consumer and 129 basal species"
      ), 
      landsc = numextract(str_split(basename(landsc), "_", simplify = TRUE)[, 2]), 
      ITL = round(TL), 
      TL05 = my_round_any(TL, 0.5), 
      params = case_when(
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 1 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.01'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 2 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.0325'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 3 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.055'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 4 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.0775'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 5 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.1'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 6 ~ "pi[i]*' and '*xi[i]*' TLB'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 7 ~ "pi[i]*' = 0.2,  '*xi[i]*' TLB'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 8 ~ "pi[i]*' = 0.2,  '*xi[i]*' = 0.055'",
        numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 9 ~ "pi[i]*' and '*xi[i]*' SGB'",
        TRUE ~ as.character(numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]))
      ),
      scenario = case_when(
        scenario == "best" ~ "best-case scenario",
        scenario == "worst" ~ "worst-case scenario",
        scenario == "random" ~ "random scenario",
        TRUE ~ "pre-patch loss"
      )
    )
    return(out)
  }
