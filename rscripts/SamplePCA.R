setClass('SamplePCA',
         representation(scores = 'matrix',
                        variances = 'numeric',
                        components = 'matrix',
                        splitter = 'ANY',
                        usecor = 'logical',
                        shift = 'numeric',
                        scale = 'numeric',
                        call = 'call'))

doPCA <- function(data, splitter = 0, usecor = FALSE, center = TRUE)
{
  call <- match.call()
  if(inherits(data, 'ExpressionSet')) {
    if(is.character(splitter)) {
      splitter <- as.factor(pData(data)[,splitter])
    }
    data <- exprs(data)
  }
  nSample <- dim(data)[2]
  centered <- as.matrix(data)
  geneMean <- rep(0, nrow(data))
  geneSD <- rep(1, nrow(data))
  if (center == TRUE) {
    geneMean <- apply(centered, 1, mean)
    centered <- sweep(centered, 1, geneMean, '-')
  }
  if (usecor == TRUE) {
    if(max(dim(centered)) > 2000) {
      geneSD <- sqrt(apply(centered, 1, var))
    } else {
      geneSD <- sqrt(diag( centered %*% t(centered)/(nSample-1) ))
    }
    geneSD[geneSD == 0] <- 1
    centered <- sweep(centered, 1, geneSD, '/')
  }
  decomp <- svd(centered)
  variances <- decomp$d^2/nSample
  components <- decomp$u
  scores <- t(decomp$d * t(decomp$v))
  rownames(scores) <- colnames(data)
  colnames(components) <- colnames(scores) <- paste("PC", 1:ncol(scores), sep='')
  rownames(components) <- rownames(data)
  new('SamplePCA', scores=scores, variances=variances, components=components,
      splitter=splitter, usecor=usecor, shift=geneMean, scale=geneSD, call=call)
}

setMethod('summary', 'SamplePCA', function(object, ...) {
  cat('A SamplePCA object.\nCall:\n\t')
  cat(as.character(list(object@call)))
  cat('\n')
})

setMethod('predict', 'SamplePCA', function(object, newdata=NULL, ...) {
  if(is.null(newdata)) {
    value <- object@scores
  } else {
    temp <- sweep(newdata, 1, object@shift, '-')
    temp <- sweep(temp, 1, object@scale, '/')
    value <- t(temp) %*% object@components
  }
  value
})

setMethod('plot', signature('SamplePCA','missing'), function(x,
                                                             splitter=x@splitter, col, main='', which=1:2, ...) {
  .my.plot <- function(z, colors, i, j, ...) {
    plot(ColorCodedPair(z[,i], z[,j], colors),
         xlab=paste('Component', i),
         ylab=paste('Component', j), ...)
  }
  
  z <- x@scores
  if (length(splitter) == 1) {
    if(missing(col)) {
      col <-  'black'
    }
    ccl <- list(ColorCoding(rep(TRUE, length(x@variances)), col[1], 15))
  } else if (is.logical(splitter)) {
    if(missing(col) || length(col) < 2) {
      col <- c('blue', 'red')
    }
    ccl <- list(ColorCoding(splitter, col[1], 15), ColorCoding(!splitter, col[2], 16))
  } else if (!is.null(dim(splitter)) && dim(splitter) > 1) {
    labels <- levels(splitter[,1])
    symbols <- levels(splitter[,2])
    symlist <- c(15:19, 1:14, 64:127)
    if(missing(col) || length(col) < length(labels)) {
      n <- length(labels)
      col <- hcl(h=seq(0, 360, length=1+n)[1:n], c=75, l=45)
    }
    ccl <- list()
    lc <- as.character(splitter[,1])
    sc <- as.character(splitter[,2])
    n <- 0
    for (i in 1:length(labels)) {
      for (j in 1:length(symbols)) {
        v <- lc==labels[i] & sc==symbols[j]
        n <- n+1
        ccl[[n]] <-  ColorCoding(v, col[i], symlist[j])
      }
    }
  } else {
    labels <- levels(splitter)
    if(missing(col) || length(col) < length(labels)) {
      n <- length(labels)
      col <- hcl(h=seq(0, 360, length=1+n)[1:n], c=75, l=45)
    }
    ccl <- lapply(1:length(labels), function(x, f, col) {
      fc <- as.character(f)
      f1 <- attr(f, 'levels')[x]
      v <- fc == f1;
      ColorCoding(v, col[x], 15)
    }, splitter, col)
  }
  if (length(which) == 1) {
    plot(c(0,4), c(0,3), xaxt='n', yaxt='n', type='n', xlab='', ylab='')
    text(2, 2, main, cex=1.2)
    text(2, 1, 'Principal Components', cex=1.2)
    .my.plot(z, ccl, 1, 2, main='PCA', ...)
    .my.plot(z, ccl, 2, 3, main='PCA', ...)
    .my.plot(z, ccl, 1, 3, main='PCA', ...)
  } else if (length(which)==2) {
    .my.plot(z, ccl, which[1], which[2], main=main, ...)
  } else {
    .my.plot(z, ccl, 1, 2, main=main, ...)
  }
  invisible(ccl)
})

setMethod('screeplot', 'SamplePCA', function(x, N=NULL, ...) {
  if(is.null(N)) N <- length(x@variances)
  N <- min(N, length(x@variances))
  barplot(x@variances[1:N]/sum(x@variances), ...)
})

setMethod("identify", "SamplePCA", function(x, ...) {
  identify(x=x@scores[,1], y=x@scores[,2], labels=rownames(x@scores), ...)
})

setMethod("text", "SamplePCA", function(x, ...) {
  text(x=x@scores[,1], y=x@scores[,2], , labels=rownames(x@scores), ...)
})

mahalanobisQC <- function(spca, N) {
  ss <- spca@scores[, 1:N]
  maha <- sapply(1:nrow(ss), function(i) {
    v <- matrix(1/apply(ss[-i,], 2, var), ncol=1)
    as.vector(ss[i,]^2 %*% v)
  })
  names(maha) <- rownames(spca@scores)
  pmaha <- 1-pchisq(maha, N)
  data.frame(statistic=maha, p.value=pmaha)
}
