
#Load packages
require(tidyverse)
require(farver) #version 2.0.3 or higher
require(scales) #version 1.01

#Create figure directory
figuredir <- "../../figures/"
if(file.exists(figuredir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(figuredir)) #create direcotry for figures, if it does not exist yet
} 

landscapes <- Sys.glob("../../data/ryser_data/landscapes/*.csv")

resfile <- "../../data/ryser_data/results/tldat.csv" #results file
dat <- resfile %>% read_csv

dat %>%
  mutate(ab=paste0("alpha==", alpha, "~~beta==", beta)) %>%
  mutate(spinput=ifelse(spinput=="web_1_1.csv", "xi~const.", "xi~TLB")) %>%
  mutate(tl=ifelse(is.infinite(tl), 1, tl)) %>%
  ggplot +
  aes(x=log10(1/xi), y=N, fill=tl) +
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0, 0),
                     name=expression(paste(log[10](spatial~scale)))) +
  scale_y_continuous(expand=c(0, 0), name="number of patches") +
  scale_fill_gradient2(name="max trophic level", midpoint=3.8, mid="white",
                       low=scales::muted("blue"), high=scales::muted("red")) +
  facet_grid(spinput~ab, labeller="label_parsed") +
  theme_bw()
##ggsave(paste0(figuredir, "Ryser_reproduced.pdf"), width=8, height=3.8)


##Function to otbtain the dominant eigenvalue of the dispersal matrix M
##Input: - lsc: landscape file, xi: dispersal distance
#Output: - lambda: the dominat eigenvalue of the dispersal matrix M
getlambda <- function(lsc, xi) {
  M <- read_csv(lsc, col_types="dd") %>%
    as.matrix(col_types=cols()) %>%
    dist %>%
    as.matrix
  M <- exp(-M/xi)
  diag(M) <- 0
  rs <- rowSums(M)
  lambda <- sum(rs*rs)/sum(rs)
  return(lambda)
}

dat2 <- expand_grid(N=seq(16, 70, by=2), #number of patches
            xi=10^(seq(-3, -1.2, l=11)), #dispersal distances 
            lsc=landscapes) %>%  #landscape files
  rowwise %>%
  mutate(lM=getlambda(lsc, xi)) %>% #obtain the dominant eigenvalues of the dispersal matrices
  ungroup %>%
  mutate(tl=lM) 

dat2 %>%
  ggplot +
  aes(x=log10(1/xi), y=N, fill=tl) +
  geom_raster(interpolate=TRUE) +
  scale_x_continuous(expand=c(0, 0),
                     name=expression(paste(log[10](spatial~scale)))) +
  scale_y_continuous(expand=c(0, 0), name="number of patches") +
  scale_fill_gradient2(name="max trophic level", low=scales::muted("blue"),
                       mid="white",
                       high=scales::muted("red"), midpoint=0.9,
                       labels=c(expression(paste(frac(1,2)~log(pi^-1))),
                                expression(paste(log(pi^-1))),
                                expression(paste(frac(3,2)~log(pi^-1)))),
                       breaks=c(0.5, 1.0, 1.5)) +
  theme_bw()
##ggsave(paste0(figuredir,"Ryser_semianalytic.pdf"), width=4, height=3)

