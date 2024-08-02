####################
#funciton: my image
#args: plotdata: a data matrix with ordered genes in row and ordered samples in column
#      nsample: a vector with number of sample in each sample groups, the sum=ncol(plotdata)
#      ngene: a vector with number of genes in ecah gene group
#####################
myimage<-function(plotdata, nsample, ngene, glabel=names(ngene),slabel=names(nsample), ...){  
layout(mat=matrix(c(1,0,0,0,2), 1, 5), 1, widths=c(4, 0, 0, 0, 0.15))
par(plt=c(0.11146374, 0.94282901, 0.08582091, 0.92848258))
ramp <- colorRamp(c("red", "white", "blue"))
	ColorRamp <- rgb( ramp(seq(0, 1, length = 255)), max = 255)
	ColorLevels <- seq(min(plotdata), max(plotdata), length=length(ColorRamp))
image(t(as.matrix(plotdata)), zlim=c(-1,1), col=ColorRamp, axes=F, ...)
box()
sn<-sum(nsample)
sat<-((cumsum(nsample)-nsample/2)-1)/(sn-1)
ls<-(cumsum(nsample)[1:(length(nsample)-1)]-1/2)/(sn-1)
gn<-sum(ngene)
gat<-((cumsum(ngene)-ngene/2)-1)/(gn-1)
lg<-(cumsum(ngene)[1:(length(ngene)-1)]-1/2)/(gn-1)
axis(2, at=gat, label=glabel)
axis(1, at=sat, label=slabel)
abline(v=ls, h=lg, lwd=3, col="green")
mtext(3, at=0:(ncol(plotdata)-1)/(ncol(plotdata)-1), text=colnames(plotdata),las=2,cex=0.7, line=0.2)
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")}
