##################################################################
#make pairplot
pairplot<-function(x,col="green",main=NULL,method=c("pearson","concord"), plot.method='smoothscatter', ...){
method<-match.arg(method)
panel.concord <- function(x1, x2,...) {
horizontal <- (par("usr")[1] + par("usr")[2]) / 2;
vertical <- (par("usr")[3] + par("usr")[4]) / 2;
#calculate cocordance coefficient in pairwise plot
	n <- length(x1)
	xbar1 <- mean(x1,na.rm=T)
	xbar2 <- mean(x2,na.rm=T)
	s1 <- sum((x1-xbar1)^2,na.rm=T)/n
	s2 <- sum((x2-xbar2)^2,na.rm=T)/n
	s12 <- sum((x1-xbar1)*(x2-xbar2),na.rm=T)/n
	pc <- (2*s12)/(s1+s2+(xbar1-xbar2)^2)
text(horizontal, vertical, format(pc, digits=4))
}
panel.pearson <- function(x1, x2,...) {
horizontal <- (par("usr")[1] + par("usr")[2]) / 2;
vertical <- (par("usr")[3] + par("usr")[4]) / 2;
text(horizontal, vertical, format(cor(x1,x2), digits=4))
}
if(plot.method=='smoothscatter') s.panel<-function(...) {par(new=TRUE);smoothScatter(..., nrpoints=0, axes=F)} else s.panel<-points
if(method=="concord") pairs(x, pch=".", lower.panel=s.panel,  upper.panel=panel.concord, col=col, main=main)
if(method=="pearson") pairs(x, pch=".", lower.panel=s.panel,  upper.panel=panel.pearson, col=col, main=main)
}
