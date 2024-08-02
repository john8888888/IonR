############################
#This file is from Nianxiang Zhang
#Do not distribute without permission!!!
#########################

###################
#Concord function to calculate Lin's ccc
#
##########

Concord <- function(x1, x2) {
stopifnot(length(x1)==length(x2))
keep.ind<-!(is.na(x1) |is.na(x2))
x1<-x1[keep.ind]
x2<-x2[keep.ind]
	n <- length(x1)
	xbar1 <- mean(x1,na.rm=T)
	xbar2 <- mean(x2,na.rm=T)
	s1 <- sum((x1-xbar1)^2,na.rm=T)/n
	s2 <- sum((x2-xbar2)^2,na.rm=T)/n
	s12 <- sum((x1-xbar1)*(x2-xbar2),na.rm=T)/n
	pc <- (2*s12)/(s1+s2+(xbar1-xbar2)^2)
return(pc)
}

###################
#Functions:
#bi bimodality
#truncscale  truncation and scale
#pairplot  make pairwise plot with pearson or Lin concord in Diagnal
#imagewithkey make image with color key
#maketable Draw a table graph with cutpoint and color
########################################################
junk<-c("bi","truncscale","pairplot","imagewithkey","maketable")
##calculation of bimodality index
bi<-function(x){
require("mclust")
temp<-Mclust(x, G=2, modelNames="E")
delta<-abs(temp$parameters$mean[2]-temp$parameters$mean[1])/sqrt(temp$parameters$variance$sigmasq)
BI<-delta*sqrt(temp$parameters$pro[1]*temp$parameters$pro[2])
return(BI)
}
#################################################################
#standardize and truncate the data for better heatmap presentation
truncscale<-function(x, di=c("all", "row","col"), tp=0.8, ...){x<-as.matrix(x)
nc<-ncol(x)
nr<-nrow(x)
dir<-match.arg(di)
if(dir=="row") temp<-t(scale(t(x),...))
else if(dir=="col")temp<-scale(x,...)
else if(dir=="all")temp<-scale(as.vector(x),...)
bound<-quantile(temp,c((1-tp)/2, 1-(1-tp)/2))
temp[temp<bound[1]]<-bound[1]
temp[temp>bound[2]]<-bound[2]
ifelse(dir=="all", return(matrix(temp, nr,nc, dimnames=list(rownames(x), colnames(x)))), return(temp))
}
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
text(horizontal, vertical, paste(Method,'=', format(pc, digits=4), sep=''))
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
##################
#imagewithkey: function for generating image with key and values in the image
#args: data, a matrix 
#      col:colors used in image
#      cex.label: text size for label
#      cex.text: text size for actual vales
#      withvalue: if False, no actual values will be shown
#      digit: actual value digits
#      las: the direction of the x-axis label
#value: a plot will be made
#author:Nianxiang Zhang, BCB, MDACC
####################
imagewithkey<-function(data, col=heat.colors(12), cex.label=1, withvalue=T, cex.text=1, digit=3, las=1, zlim=NULL, p=F, min.p=0.001, ...){
omar<-par("mar")
par(mar=c(omar[1:3],1)) 
layout(mat=matrix(1:2, 1,2), c(10,1))
if(is.null(zlim)) image(data, col=col, axes=F, ...) else image(data, col=col, axes=F, zlim=zlim, ...)
axis(2, at=0:(ncol(data)-1)/(ncol(data)-1), lwd=1, labels=colnames(data), las=2, cex.axis=cex.label )
axis(1, at=0:(nrow(data)-1)/(nrow(data)-1),lwd=1, labels=rownames(data), cex.axis=cex.label, las=las)
if(withvalue) {
		for(i in 1:ncol(data)){
		for(j in 1:nrow(data)){
		if(p){
		text(y=c(0:(ncol(data)-1)/(ncol(data)-1))[i], x=c(0:(nrow(data)-1)/(nrow(data)-1))[j], labels=ifelse(data[j,i]<min.p, paste("<", substitute(min.p),sep=""),round(data[j,i], digit)), cex=cex.text)} else {
		text(y=c(0:(ncol(data)-1)/(ncol(data)-1))[i], x=c(0:(nrow(data)-1)/(nrow(data)-1))[j], labels=round(data[j,i], digit), cex=cex.text)} 
		}}
		}
box()
par(mar=c(omar[1],0,omar[3],2.1))
if(all(na.omit(data)==1)) zlim<-c(0,1)
if(is.null(zlim)) {image(t(seq(min(data),max(data[is.finite(data)]), length=length(col))), axes=F, col=col) 
mtext(side=4, at=seq(0, 1, length=4), text=round(seq(min(data),max(data[is.finite(data)]), length=4),2), las=2)
}else {image(t(seq(min(zlim),max(zlim), length=length(col))), axes=F, col=col)
mtext(side=4, at=seq(0, 1, length=4), text=round(seq(min(zlim),max(zlim), length=4),2), las=2)
}
par(mar=omar)
}

######################################################################
#Draw a table graph with cutpoint and color
imagecolor<-c("blue","lightblue","white", "gold","darkorange")
maketable<-function(oridata, breakpoint, color=heat.colors(12), text.cex=0.7, label.cex=0.7, digit=2, ...){
breakpoint<-sort(breakpoint, decreasing=T)
nb<-length(breakpoint)
tabdata<-oridata
tabdata[,]<-NA
tabdata[oridata>=breakpoint[1]]<-nb
for(i in 1:nb){
tabdata[oridata<breakpoint[i]]<-nb-i
}
layout(mat=matrix(c(rep(1, nb), 0,2:(nb),0),nrow=2, byrow=T), height=c(20,1))
par(mar=c(0.5,5.1,5.1,1))
image(t(tabdata[nrow(tabdata):1,]), col=color[min(tabdata):max(tabdata)], axes=F, ...)
atx<-(0:(ncol(tabdata)-1))/(ncol(tabdata)-1)
aty<-(0:(nrow(tabdata)-1))/(nrow(tabdata)-1)
mtext(3, line=0.1, text=colnames(tabdata), at=atx, las=2, cex=label.cex)
mtext(2, line=0.1,text=rownames(tabdata)[nrow(tabdata):1], at=aty, las=1, cex=label.cex)
box()
abline(h=((0:(nrow(tabdata)-1))-0.5)/(nrow(tabdata)-1))
abline(v=((0:(ncol(tabdata)-1))-0.5)/(ncol(tabdata)-1))
x<-rep(atx, nrow(tabdata))
y<-rep(aty, each=ncol(tabdata))
text(x,y, format(round(as.vector(t(oridata[nrow(oridata):1,])),digit),digit=digit, scientific=F), cex=text.cex)
for(j in 1:(nb)){
par(mar=c(0.5,0.1,0.1,7.1), las=1)
image(as.matrix(1),col=color[j], axes=F)
box()
if(j==1) mtext(side=4, line=0.1, text=paste("<", breakpoint[nb], sep=""), cex=label.cex) else if(j==nb) mtext(side=4, line=0.1, text=paste(">", breakpoint[1], sep=""), cex=label.cex) else  mtext(side=4, line=0.1, text=paste( breakpoint[nbj],"-", breakpoint[nbj], sep=""), cex=label.cex)
}
par(las=0)
}
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
      
#################
#get confidence intervals for coxph model and glm.
#      
`intervals` <-
function(object, level = 0.95, ...)
UseMethod("intervals")

`intervals.coxph` <-
function(x, level=.95, sign = 1, FUN=exp)
{
z<-abs(qnorm((1-level)/2))
co <- summary(x)$coef
or <- FUN(co[, 1] * sign)
li <- FUN(co[, 1]* sign -z * co[, 3])
ls <- FUN(co[, 1]* sign +z * co[, 3])
r <- cbind(or, li, ls, co[,5])
dimnames(r) <- list(dimnames(co)[[1]], c("hr", paste(level*100,"%",sep=""), "C.I.","P-val"))
class(r)<-c("intervals","coxph.intervals")
r
}

`intervals.glm` <-
function(x, level=.95, sign = 1, FUN=exp)
{
z<-abs(qnorm((1-level)/2))
co <- summary(x)$coef
or <- FUN(co[, 1] * sign)
li <- FUN(co[, 1]* sign -z * co[, 2])
ls <- FUN(co[, 1]* sign +z * co[, 2])
r <- cbind(or, li, ls, co[,4])
dimnames(r) <- list(dimnames(co)[[1]], c("or", paste(level*100,"%",sep=""), "C.I.","P-val"))
class(r)<-c("intervals","glm.intervals")
r
}

`print.intervals` <-
function(n, len = 6, d = 2, exclude.intercept=inherits(n,"glm.intervals"), pval=TRUE)
{
        dd <- dim(n)
        mx <- 10^(len-(d+1))
        n[n > mx] <- Inf
        a <- formatC(n, d, len,format="f")

        dim(a) <- dd
        if(length(dd) == 1){
                dd<-c(1,dd)
                dim(a)<-dd
                lab<-" "
        }
        else    lab <- dimnames(n)[[1]]

        if(!pval){
                mx <- max(nchar(lab)) + 1
                cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
                for(i in (1+exclude.intercept):dd[1]) {
                        lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
                        cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
                }
        } else {
                mx <- max(nchar(lab)) + 1
                cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
                for(i in (1+exclude.intercept):dd[1]) {
                        lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
                        cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") ",formatC(n[i,4], 4, 6,format="f"),"\n")
                }
        }
}
