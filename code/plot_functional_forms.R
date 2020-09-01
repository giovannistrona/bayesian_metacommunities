#Load libraries
library(tidyverse, quietly=TRUE)

ggplot_fform <- function(alpha, beta)
{
  x <- seq(0, 1, l=101)
  fform <- data.frame(
    x = x,
    y = pbeta(x, alpha, beta),
    alpha = alpha, 
    beta = beta, 
    ff = case_when(alpha==1 & beta==1 ~ "linear",
                      alpha==1 & beta==5 ~ "concave",
                      alpha==5 & beta==1 ~ "convex",
                      alpha==5 & beta==5 ~ "sigmoid"))
 
  facetnames <- case_when(alpha==1 & beta==1 ~ expression("linear "* alpha *" = "* beta * " = 1"),
                        alpha==1 & beta==5 ~ expression("concave "* alpha *" = 1 "* beta * " = 5"),
                        alpha==5 & beta==1 ~ expression("convex "* alpha * " = 5 "* beta *" = 1"),
                        alpha==5 & beta==5 ~ expression("sigmoidal "* alpha *" = " *beta *" = 5"))
  
  levels(fform$ff) <- facetnames
  par(pty="s")
  ggplot(fform, aes(x, y)) +
    theme_bw() + 
    geom_hline(yintercept=0, lty=2, lwd=0.8, color="black", alpha = 0.5) +
    geom_hline(yintercept=1, lty=2, lwd=0.8, color="black", alpha = 0.5) +
    geom_line(lwd=1) +
    xlab("Fraction of resources \n extinct") +
    ylab("Prob. of extinction \n of consumers") +
    scale_x_continuous(breaks=c(0, 0.5, 1),
                       labels = c("0","0.5","1")) +
    scale_y_continuous(breaks=c(0, 0.5, 1),
                       labels = c(expression(pi[i]),"0.5","1")) +
    theme(axis.text = element_text(color="black", size = 18), 
          axis.title = element_text(color="black", size=18), 
          strip.text.x = element_text(size=16)) +
    facet_grid(.~ff, labeller=label_parsed) 
}


plot_functional_form <- function(figuredir, alpha, beta){
      ggsave(
        paste0(figuredir, "fform_",alpha,"_",beta,".pdf"),
        ggplot_fform(alpha, beta),
        width = 3.25,
        height = 3.25
      )
}

for(a in c(1,5))
    for(b in c(1,5))
      plot_functional_form("../figures/",a,b)


## Empty place holder plot for latex 
## ggsave(
##   "../figures/fform_0_0.pdf",
##   ggplot() + theme_void(),
##   width = 3.25,
##   height = 3.25
## )    
