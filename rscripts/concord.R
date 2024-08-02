concord <- function(x1, x2) {
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
