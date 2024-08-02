###################
#Functions:
#bi bimodality
#truncscale  truncation and scale
#pairplot  make pairwise plot with pearson or Lin concord in Diagnal
#imagewithkey make image with color key
#maketable Draw a table graph with cutpoint and color
########################################################
##calculation of bimodality index
bi<-function(x){
require("mclust")
temp<-Mclust(x, G=2, modelNames="E")
delta<-abs(temp$parameters$mean[2]-temp$parameters$mean[1])/sqrt(temp$parameters$variance$sigmasq)
BI<-delta*sqrt(temp$parameters$pro[1]*temp$parameters$pro[2])
return(BI)
}
