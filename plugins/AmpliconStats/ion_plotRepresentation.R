library(scales)

args <- commandArgs(trailingOnly = TRUE)
repFile <- args[1]
pngFile <- args[2]
png(pngFile, width=4, height=2, units="in", res=600, pointsize=4)

rep <- read.csv(repFile,header=T)
repF <- rep[,2];
repR <- rep[,3]
ymax = max(max(repF),max(repR))
ymin = 0 - ymax
ttl <- repF + repR

error = sd(repF)
mu = mean(repF)
confInt95 = c(mu - error, mu + error)
confInt95[is.na(confInt95)] = 0

plot(1:length(repF), 
			repF[order(ttl)], 
			type="h", 
			xlim=c(1,length(repF)), 
			ylim=c(ymin,ymax), 
			xaxs = "i",
			yaxs = "i",
			#col=alpha('orange',0.9), 
			col='green',

			main="Total Number of Amplicon Reads (sorted)", 
			xlab="sorted amplicon", 
			ylab="total # reads"
)
abline(confInt95[1], 0, col="darkgray", lty="dashed")
abline(mu, 0, col="red", lty="dashed")
abline(confInt95[2], 0, col="darkgray", lty="dashed")
rect(0, confInt95[1], length(repF), confInt95[2], col=alpha('gray', 0.3), lty='blank')

par(new=T)

repR <- 0 - repR
error = sd(repR)
mu = mean(repR)
confInt95 = c(mu - error, mu + error)
confInt95[is.na(confInt95)] = 0

plot(1:length(repR), 
			repR[order(ttl)], 
			type="h", 
			xlim=c(1,length(repR)), 
			ylim=c(ymin,ymax), 
			xaxs = "i",
			yaxs = "i",
			#col=alpha('green',0.9), 
			col='orange', 
			main="", xlab="", ylab=""
)
abline(confInt95[1], 0, col="darkgray", lty="dotdash")
abline(mu, 0, col="red", lty="dotdash")
abline(confInt95[2], 0, col="darkgray", lty="dotdash")
rect(0, confInt95[1], length(repF), confInt95[2], col=alpha('gray', 0.3), lty='blank')

par(new=F)

dev.off()