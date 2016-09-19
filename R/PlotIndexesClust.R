
#' Graphic visualization of an index
#'
#' The function generates graphs and tables of a specific index and its
#' clusters. The user decides whether to supply the k or let the function
#' calculate the k for the selected index
#'
#' @param x list of expression data and their indexes after running
#' \code{\link{PreparingTheIndexes}}
#' @param idx A character, the index to plot
#' (e.g., for 8 time points "11100-1-1-1")
#' @param k A numeric, number of clusters. If not given the function will
#' calculate what is the optimal k for the index.
#' @param scaling A boolean, default to TRUE,
#' does the data should be standardized before clustered with K-means.
#'
#' @return A list with two objects:
#'  \enumerate{\item Table of of a specific index and its clusters
#'             \item Gene expression pattern graphs
#'             for each one of the clusters}
#'
#' @seealso \code{\link{ggplot}},  \code{\link{kmeans}},
#'  \code{\link{ClustIndexes}}
#'
#' @examples
#'
#' data_dir <- system.file("extdata", package = "ctsGE")
#' files <- dir(path=data_dir,pattern = "\\.xls$")
#' rts <- readTSGE(files, path = data_dir,
#' labels = c("0h","6h","12h","24h","48h","72h"), skip = 10625 )
#' prts <- PreparingTheIndexes(rts)
#' pp <- PlotIndexesClust(prts,idx="000011")
#' pp$graphs # plots the line graphs
#'
#' @export
#' @import ggplot2 reshape2 stats
#'
#'

PlotIndexesClust = function(x,idx,k=NULL,scaling=TRUE){
    set.seed(100)
    plot <- list()
    ggplot_list <- list()
    clust_tbl <-  list()
    kindex <-  list()
    genes <- rownames(x$index[x$index[,"index"]==idx,])

    if(length(x) < 7){stop("Please run PreparingTheIndexes first")}
    if(!scaling){tbl <- x$tsTable[genes,]
    } else{tbl <- x$scaled[genes,]}

    if(length(genes) > 1 ) {tbl <- as.data.frame(tbl)
    }else{
        tbl <- data.frame(t(matrix(tbl)))
        colnames(tbl) <- x$samples
        rownames(tbl) <- genes
    }
    if(!is.null(k)){
        fit_km <- kmeans(tbl,k)
        tbl[names(fit_km$cluster),"clusters"] <- fit_km$cluster
        K <- k
    }else{
        # checks if indexes were already made with PreparingTheIndexes
        if("optimalK" %in% names(x)){
            k <- x$optimalK[rownames(x$optimalK)==idx,"k"]
            tbl1 <- x$ClusteredIdxTable
            tags <- rownames(tbl1[tbl1[,"index"]==idx,])
            tbl[tags,"clusters"] <- tbl1[tbl1[,"index"]==idx,"clusters"]
        }else{
            if(length(genes) > 9){
                tbl <- as.data.frame(tbl)
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
            tbl[,"clusters"] <- clust
        }
    }

    tmp <- cbind(tags=genes,tbl)
    for(i in 1:K){
        x.m <- reshape2::melt(tmp[tmp$cluster==i,c("tags",x$samples)])
        colnames(x.m) <- c("tags","tp","exp")
        gplot <-
            ggplot2::ggplot(x.m,
                            ggplot2::aes_string(y="exp",x="tp",colour="tags",
                                                group="tags"))+
            ggplot2::labs(title =paste0("Index: ",idx," Cluster: ",i),
                          x = "Time", fill= NULL, y = "Expression Value") +
            ggplot2::theme(legend.position="none")

        gplot <- gplot +ggplot2::geom_line()
        name <- paste0("Clust_",i,"_",idx)
        ggplot_list[[name]] <- gplot
    }
    if ("desc" %in% names(x)) {tbl <- cbind(desc=x$desc[genes,], tbl)}
    plot[[name]] <- tbl
    plot$graphs <- ggplot_list
    structure(plot,class = "list")
}
