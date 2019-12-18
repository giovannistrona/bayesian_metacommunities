
plot_functional_forms <- function(figuredir){

  par(pty="s")
  frac <- sort(runif(10000))
  
  for(alpha in c(1,5)){
    for(beta in c(1,5)){
      png(paste0(figuredir,"fform_",alpha,"_",beta,".png"))
      f <- pbeta(frac, alpha, beta)
      plot(f~frac,type='n',xlab="",ylab="",xaxt="n",yaxt="n",lwd=3)
      lines(f~frac,col="grey21",lwd=8)
      abline(h=max(f), lty=3,lwd=6)
      abline(h=min(f), lty=3,lwd=6)
      dev.off()
    }
  }
}
