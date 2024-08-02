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
