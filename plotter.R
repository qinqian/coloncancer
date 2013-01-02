########################
## plot except heatmap
########################

library(geneplotter)
library(gclus)
library(hexbin)

myGeneplot <- function(){
  cPlot
  histStack
  ## library(rgl)
  ## scatter3d(sort(seq(length(t.p.value))), sort(t.p.value.unif), sort(t.p.value))
  ## smoothScatter(sort(t.p.value))
  
}

mycorplot <- function(){
  panel.plot <- function( x,y, ... )
{
  par(new=TRUE)
  m <- cbind(x,y)
  plot(m,col=densCols(m),pch=20)
  lines(lowess(m[!is.na(m[,1])&!is.na(m[,2]),]),col="red")  
}
    
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,use="complete.obs")
  txt <- format(round(r,2),width=5,nsmall=2)
  #format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  text(0.5, 0.5, txt, cex = cex.cor)
}
}

