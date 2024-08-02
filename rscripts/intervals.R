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
