#' Clustering the indexes applying K-means
#'
#' Clustering each index, that was predifined by
#' \code{\link{PreparingTheIndexes}}, with \code{\link{kmeans}}.
#'
#' @param x  A ctsGEList object
#' @param scaling Boolean parameter, does the data should be standardized
#' before clustered. Default = TRUE
#'
#'
#' @return  ctsGEList object is returned as output,
#' with the relative culstered indexes table in
#'  \emph{object$ClusteredIdxTable},
#' and the number of clusters for each index in \emph{object$optimalK}
#'
#'
#' @details The clustering is done with K-means. To choose an optimal k for
#'  K-means clustering, the Elbow method was applied, this method looks at the
#'  percentage of variance explained as a function of the number of clusters: the
#'  chosen number of clusters should be such that adding another cluster does
#'  not give much better modeling of the data. First, the ratio of the
#'  within-cluster sum of squares (WSS) to the total sum of squares (TSS) is
#'  computed for different values of k (i.e., 1, 2, 3 ...). The WSS, also known
#'  as sum of squared error (SSE), decreases as k gets larger. The Elbow method
#'  chooses the k at which the SSE decreases abruptly. This happens when the
#'  computed value of the WSS-to-TSS ratio first drops from 0.2.
#'
#'  Running \code{\link{kmeans}} and calculating the optimal k for each one of
#'  the indexes in the data could take a long time. To shorten the procedure the
#'  user can skip this step altogether and directly view a specific index and
#'  its clusters by running either the \code{\link{PlotIndexesClust}} or the
#'  \code{\link{ctsGEShinyApp}} function.
#'
#'  By default data is standardize before clustering,for clustering
#'  the raw counts set the \command{\strong{scaling}} parameter to FALSE.
#'
#'
#' @examples
#'
#' data_dir <- system.file("extdata", package = "ctsGE")
#' files <- dir(path=data_dir,pattern = "\\.xls$")
#' rts <- readTSGE(files, path = data_dir,
#' labels = c("0h","6h","12h","24h","48h","72h"), skip = 10625 )
#' prts <- PreparingTheIndexes(rts)
#'
#' tsCI <- ClustIndexes(prts)
#'
#' head(tsCI$ClusteredIdxTable) #the table with the clusterd indexes
#' head(tsCI$optimalK) #the table with the number of clusters for each index
#'
#' @seealso \code{\link{kmeans}}, \code{\link{PlotIndexesClust}}
#' @export
#' @import stats
#'
ClustIndexes = function(x,scaling=TRUE){
    clust_tbl <-  list()
    kindex <-  list()

    if(scaling){tmp <-
        cbind(as.data.frame(x$scaled),index=x$index[rownames(x$scaled),
                                                    ncol(x$index)])
    }else { tmp <-
        cbind(as.data.frame(x$tsTable),index=x$index[rownames(x$scaled),
                                                     ncol(x$index)])}


    for(idx in as.character(unique(tmp$index))){
        tbl <- tmp[tmp$index==idx,1:x$timePoints]
        if(nrow(tbl) > 9){
            k = 1
            fit_km <- kmeans(tbl,k,nstart = 25)
            opt <- fit_km$tot.withinss/fit_km$totss < 0.2
            while (!opt) {
                k=k+1
                # Apply k-means to tbl: fit_km
                fit_km <- kmeans(tbl,k,nstart = 25)
                opt <- fit_km$tot.withinss/fit_km$totss < 0.2
            }
            K <-  k #where: WSS / TSS < 0.2 this is the optimal k
            clust <-   fit_km$cluster
        }else{
            K <- 1
            clust <- rep(1,nrow(tbl))}

        kindex[[idx]] <-  c(nrow(tbl),K)
        clust_tbl[[idx]] <- data.frame(clusters=clust,index=idx)
        rownames(clust_tbl[[idx]]) <- rownames(tbl)
    }

    names(clust_tbl) <-  NULL
    x$optimalK <-  do.call("rbind",kindex)
    colnames(x$optimalK) <- c("tags","k")
    x$ClusteredIdxTable <- do.call("rbind",clust_tbl)

    structure(x,class = "ctsGEList")
}
