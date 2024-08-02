MAtoEset<-function (from, idColumn = "PROBE_ID") 
{
    stopifnot(inherits(from, "MAList"), !is.null(from$targets), 
        is.data.frame(from$genes), idColumn %in% names(from$genes))
    from$M <- as.matrix(from$M)
    stopifnot(nrow(from$targets) == ncol(as.matrix(from$M)))
    if (is.null(rownames(from$targets))) 
        rownames(from$targets) <- as.character(from$targets[[1]])
    colnames(from$M) <- rownames(from$targets)
    myPD <- new("AnnotatedDataFrame", data = from$targets, varMetadata = data.frame(varLabel = colnames(from$targets), 
        row.names = colnames(from$targets)))
    myEset <- new("ExpressionSet", exprs = from$M, phenoData = myPD)
    featureNames(myEset) <- make.names(from$genes[[idColumn]], 
        unique = TRUE)
    featureData(myEset) <- as(from$genes, "AnnotatedDataFrame")
    return(myEset)
}