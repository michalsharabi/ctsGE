#' Define an expression index for each gene
#'
#' Reads the table of genes expression and return an expression index for each
#' gene.
#'
#'
#' @param x list of an expression data that made by readTSGE
#' @param cutoff A numeric that define the degree of change in gene
#' espression rate.
#'  \emph{See Details}.
#' @param mad.scale A boolean defaulting to TRUE as to what method of
#' scaling to use.
#' Default median-base scaling. FALSE, mean-base scaling.
#' @return ctsGEList object is returned as output with the relative
#' standarization table in \emph{object$scaled},
#' and the indexes table in \emph{object$index}
#'
#' @details 1. First, the expression matrix is standardized. The function
#' default standardizing method is a median-based scaling; alternatively, a
#' mean-based scaling can be used. The new scaled values represent the distance
#' of each gene at a certain time point from its center, median or mean,
#' in median absolute deviation (MAD) units or standard deviation (SD) units,
#' respectively.
#'
#' 2. Next, the standardized values are converted to index values that indicate
#' whether gene expression is above, below or within the limits around the
#' center of the time series, i.e., **1 / -1 / 0**, respectively. The user
#' defines a parameter cutoff that determines the limits around the
#' gene-expression center. Then the function calculates the index value at each
#' time point according to:
#'
#' \enumerate{
#'      \item \bold{0:}  standardized value is within the limits (+/- cutoff)
#'      \item \bold{1:}  standardized value exceeds the upper limit (+ cutoff)
#'      \item \bold{-1:} standardized value exceeds the lower limit (- cutoff)}
#'
#'
#' @examples
#' data_dir <- system.file("extdata", package = "ctsGE")
#' files <- dir(path=data_dir,pattern = "\\.xls$")
#' rts <- readTSGE(files, path = data_dir,
#' labels = c("0h","6h","12h","24h","48h","72h"), skip = 10625 )
#' prts <- PreparingTheIndexes(rts)
#'
#'@seealso \code{\link{scale}} \code{\link{index}}
#'
#' @export
#' @import stats
#'
PreparingTheIndexes = function(x,cutoff=1,mad.scale=TRUE){
    if(!is.list(x)) stop ( "data must be a list")
    if(sum(names(x)==c("tsTable","samples","tags","timePoints")) < 4)
        stop ( "Your List miss one or more of theses objects:
               tsTable, samples, tags, timePoints")
    tp <- x$timePoints
    x$scaled <-
        tmp <-
        t(scale(t(x$tsTable),apply(t(x$tsTable),2,median),
                apply(t(x$tsTable),2,mad)))
    if(!mad.scale) x$scaled <- tmp <-  t(scale(t(x$tsTable)))
    idx <- data.frame(t(apply(tmp,1,index,cutoff=cutoff)))
    colnames(idx) <- x$samples
    x$index <- cbind(as.data.frame(idx),index=apply(idx,1,paste,collapse=""))
    x$cutoff <- cutoff
    structure(x,class = "list")
}

#' Indexing function
#'
#' Takes a numeric vector and return an expression index
#' (i.e., a sequence of 1,-1, and 0)
#' @param x A numeric
#' @param cutoff A numeric, dermine the threshold for indexing, Default = 1
#' @return Gene expression index
#' @details The function defines limits around the center (median or mean),
#' +/- cutoff value in median absolute deviation (MAD) or standard deviation
#' (SD) units respectively.The user defines a parameter cutoff that determines
#' the limits around the gene-expression center. Then the function calculates
#' the index value at each time point according to:
#'
#' \enumerate{
#'      \item \bold{0:}  standardized value is within the limits (+/- cutoff)
#'      \item \bold{1:}  standardized value exceeds the upper limit (+ cutoff)
#'      \item \bold{-1:} standardized value exceeds the lower limit (- cutoff)}
#'
#' @seealso \code{\link{PreparingTheIndexes}}
#' @examples
#' rawCounts <-
#' c(103.5, 75.1, 97.3, 27.12, 34.83, 35.53, 40.59, 30.84, 16.39, 29.29)
#'
#' (sCounts <- scale(rawCounts)[,1])# standardized mean-base scaling
#'
#' cutoff <- seq(0.2,2,0.1)    # different cutoff produce different indexes
#'
#' for(i in cutoff){print(index(sCounts,i))}
#' @export
#'

index=function(x,cutoff=1){ return(unlist(lapply(x,function(x){
    if(x<=(-cutoff)){y=-1}
    if(x>=cutoff){y=1}
    if(x>(-cutoff)&&x<(cutoff)){y=0}
    return(y) })))}
