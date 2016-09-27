#' Define an expression index for each gene
#'
#' Reads the table of genes expression and return an expression index for each
#' gene.
#'
#'
#' @param x list of an expression data that made by readTSGE
#' @param min_cutoff A numeric the lower limit range to calculate the optimal
#' cutoff for the data, default to 0.5  \emph{See Details}.
#' @param max_cutoff A numeric the upper limit range to calculate the optimal
#' cutoff for the data, default to 0.7  \emph{See Details}.
#' @param mad.scale A boolean defaulting to TRUE as to what method of
#' scaling to use.
#' Default median-base scaling. FALSE, mean-base scaling.
#' @return list object is returned as output with the relative
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
#' 2. The function compute the cutoff value following the idea that the
#' clustering will be performed on small gene groups, an optimal cutoff value
#' will be one that will minimize the number of genes in each group,
#' i.e., generate index groups of equal size. The chi-squared values will be
#' generate for each cutoff value (from min_cutoff to max_cutoff parameter
#' in increments of 0.05) the cutoff that generate the lowest chi-squared is
#' chosen.
#'
#' 3. Next, the standardized values are converted to index values that indicate
#' whether gene expression is above, below or within the limits around the
#' center of the time series, i.e., **1 / -1 / 0**, respectively. The cutoff
#' parameter determines the limits around the gene-expression center.
#' Then the function calculates the index value at each
#' time point according to:
#'
#' \enumerate{
#'      \item \bold{0:}  standardized value is within the limits (+/- cutoff)
#'      \item \bold{1:}  standardized value exceeds the upper limit (+ cutoff)
#'      \item \bold{-1:} standardized value exceeds the lower limit (- cutoff)}
#'
#'
#'
#' @examples
#' data_dir <- system.file("extdata", package = "ctsGE")
#' files <- dir(path=data_dir,pattern = "\\.xls$")
#' rts <- readTSGE(files, path = data_dir,
#' labels = c("0h","6h","12h","24h","48h","72h"), skip = 10625 )
#' prts <- PreparingTheIndexes(rts)
#' prts$cutoff # the optimal cutoff
#'
#'@seealso \code{\link{scale}} \code{\link{index}}
#'
#' @export
#' @import stats ccaPP stringr
#'
PreparingTheIndexes = function(x,min_cutoff=0.5,max_cutoff=0.7,mad.scale=TRUE){
    # function to test cutoff values
    TstCutoff <- function(x,min_cutoff=0.5,max_cutoff=0.7){
        c <- seq(min_cutoff,max_cutoff,0.05)
        idx_tbl <- lapply(c,function(i){apply(x,1,index,cutoff=i)})
        idx_str <- lapply(idx_tbl,
                          function(y){apply(t(y),1,stringr::str_c,collapse="")})

        chisq_val <- lapply(idx_str,
                            function(z){chisq.test(table(z))$statistic[[1]]})
        i = which(unlist(chisq_val)==min(unlist(chisq_val)))
        return(list(cutoff=c[i],idx_tbl=t(idx_tbl[[i]]),idx_str = idx_str[[i]]))}

    if(!is.list(x)) stop ( "data must be a list")
    if(sum(names(x)%in%c("tsTable","samples","tags","timePoints","desc" )) < 4)
        stop ( "Your List miss one or more of theses objects:
               tsTable, samples, tags, timePoints")
    if(min_cutoff > max_cutoff)
        stop ("min_cutoff value must be smaller than max_cutoff")
    tp <- x$timePoints
    # compute the MAD faster and return median and MAD for each row
    fastMad <- apply(t(x$tsTable),2,ccaPP::fastMAD)
    mad_val <- grep("MAD",names(unlist(fastMad)))
    med_val <- grep("center",names(unlist(fastMad)))

    x$scaled <-
        tmp <-
        t(scale(t(x$tsTable),unlist(fastMad)[med_val],unlist(fastMad)[mad_val]))

    if(!mad.scale) x$scaled <- tmp <-  t(scale(t(x$tsTable)))


    tmp <- TstCutoff(tmp,min_cutoff,max_cutoff)
    x$index <- cbind(as.data.frame(tmp$idx_tbl),index=tmp$idx_str)
    x$cutoff <- tmp$cutoff

    return(x)
}



#' Indexing function
#'
#' Takes a numeric vector and return an expression index
#' (i.e., a sequence of 1,-1, and 0)
#' @param x A numeric
#' @param cutoff A numeric, dermine the threshold for indexing
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

index=function(x,cutoff){
    y=x
    y[x>=cutoff]=1
    y[x<=cutoff]=-1
    y[x>-cutoff&x<cutoff]=0
    return(y)}
