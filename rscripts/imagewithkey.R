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
