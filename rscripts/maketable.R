######################################################################
#Draw a table graph with cutpoint and color
#imagecolor<-c("blue","lightblue","white", "gold","darkorange")
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
