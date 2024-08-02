corPlot<-function (eset, samples = NULL, grouping = NULL, ref = NULL, 
    useSmoothScatter = TRUE, ...) 
{
    stopifnot(inherits(eset, "ExpressionSet") | (is.numeric(eset) & 
        is.matrix(eset)))
    if (inherits(eset, "ExpressionSet")) 
        eset <- exprs(eset)
    if (!is.null(samples)) {
        notValidSamples <- switch(class(samples), character = !(samples %in% 
            colnames(eset)), numeric = !(samples %in% 1:ncol(eset)))
        if (any(notValidSamples)) 
            stop(paste("Samples", samples[notValidSamples], "not defined!\n"))
    }
    else {
        if (is.null(colnames(eset))) 
            samples <- 1:ncol(eset)
        else samples <- colnames(eset)
    }
    if (!is.null(grouping)) {
        stopifnot(length(grouping) == ncol(eset))
        grouping <- factor(grouping)
    }
    else {
        groupvec <- rep(NA, ncol(eset))
        names(groupvec) <- colnames(eset)
        groupvec[samples] <- samples
        grouping <- factor(groupvec)
    }
    if (!is.null(ref)) 
        grouping = relevel(grouping, ref = ref)
    ngroups <- nlevels(grouping)
    groupmat <- matrix(0, nrow = ncol(eset), ncol = ngroups)
    for (i in 1:ngroups) {
        thisLevel <- levels(grouping)[i]
        inGroup <- (grouping == thisLevel)
        nInGroup <- sum(inGroup, na.rm = TRUE)
        groupmat[inGroup, i] <- 1/nInGroup
    }
    datmat <- eset %*% groupmat
    colnames(datmat) <- as.character(levels(grouping))
    if (useSmoothScatter) 
        pairs(datmat, lower.panel = function(...) {
            par(new = TRUE)
            smoothScatter(..., nrpoints = 0)
            abline(0, 1, col = "red")
        }, upper.panel = panel.cor)
    else pairs(datmat, lower.panel = panel.scatter, upper.panel = panel.cor)
    invisible(NULL)
}