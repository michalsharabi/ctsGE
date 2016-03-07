#' Read and merge a set of files containing GE data
#'
#' Reads and merges a set of text files containing gene expression data
#'
#' @param files character vector of filenames, or alternative a named list of tables for each time point.
#' @param path character string giving the directory containing the files. The default is the current working directory.
#' @param columns numeric vector stating which two columns contain the tag names and counts, respectively
#' @param labels character vector giving short names to associate with the libraries. Defaults to the file names.
#' @param ... other are passed to read.delim
#'
#' @return tsGEList
#'
#' @details Each file is assumed to contained digital gene expression data for one sample (or library),
#' with transcript or gene identifiers in the first column and counts in the second column.
#' Transcript identifiers are assumed to be unique and not repeated in any one file. By default, the files are assumed to be tab-delimited and to contain column headings.
#' The function forms the union of all transcripts and creates one big table with zeros where necessary.
#'
#' @examples
#' #  Read all .txt files from current working directory
#' \dontrun{
#' files <- dir(pattern="*\\.txt$")
#' RTS <- readTSGE(files)}
#' @export
#'
readTSGE = function(files,path = NULL,columns=c(1,2),labels = NULL,...){

  x <- list()
  if(!is.list(files)){
    d <- taglist <- list()
    for (fn in files) {
      if (!is.null(path))

        fn <- file.path(path, fn)
      d[[fn]] <- read.delim(fn, ..., stringsAsFactors = FALSE)
      taglist[[fn]] <- as.character(d[[fn]][, columns[1]])
      if (any(duplicated(taglist[[fn]]))) {
        stop(paste("Repeated tag sequences in", fn))
      }
    }
  }
  else{
    d <- files
    taglist <- list()

    for (fn in 1:length(files)) {
      taglist[[fn]] <- as.character(d[[fn]][,columns[1]])
      if (any(duplicated(taglist[[fn]]))) {
        stop(paste("Repeated tag sequences in table no.",i,"in the list"))
      }
    }
  }

  tags <- unique(unlist(taglist))
  ntags <- length(tags)
  nfiles <- length(files)
  x$tsTable <- matrix(0, ntags, nfiles)
  rownames(x$tsTable) <- tags
  colnames(x$tsTable) <- labels

  if (is.null(colnames(x$tsTable))){
    if(!is.list(files)){ colnames(x$tsTable) <- x$samples <- limma::removeExt(files)
    }else{colnames(x$tsTable) <- x$samples <- names(files)}
  }else{x$samples <- colnames(x$tsTable)}

  for (i in 1:nfiles) {
    aa <- match(taglist[[i]], tags)
    x$tsTable[aa, i] <- as.numeric(d[[i]][, columns[2]])
  }

  x$tags <- tags
  x$timePoints <- nfiles
  structure(x,class = "ctsGEList")

}
#' Predefine set of profiles
#'
#' Reads the table of genes expression and return a set of model expression profiles,
#' that representive of any gene expression pattern.
#'
#' @param x A ctsGEList object
#' @param cutoff A numeric that define the degree of change in gene espression rate. \emph{See Details}.
#' @param mad.scale A boolean defaulting to TRUE as to whether or not to use \code{\link{MadScale}} for standardize the data, choose FALSE for standardization with \code{\link{scale}} function.
#'
#' @return ctsGEList object is returned as output with the relative standarization table in \emph{object$scaled}, and the profiles table in \emph{object$profiles}
#'
#' @details The raw expression data values are first standardized with \code{\link{MadScale}}
#'
#' After the scaling step, there is an indexing step to convert each time point to \strong{0 / 1 / -1}:
#'
#' \command{cutoff} - Define the degree of change in gene expression level throughout the given time points.
#'
#' The user can decide what is the degree, in which the expression level, is not significant enough or not informative enough for their experiment purpose.
#'
#'  For example: if cutoff is set to 0.7, genes that their standardized value in a certain time point is:
#'\enumerate{
#'        \item between 0.7 to -0.7 will get the index \bold{0} - no significant change has made in gene expression at that time-point
#'        \item greater than 0.7 will get the index \bold{1} - gene expression was up regulated at that time-point
#'        \item smaller than -0.7 will get the index  \bold{-1} - gene expression was down regulated at that time-point}
#'
#'  Please see example in \code{\link{index}}, on how different cutoffs produce different expression patterns,
#'  we highly recommand to test different \bold{cutoffs} in order to set the one that suits your data the most.
#'
#' @examples
#' \dontrun{
#' files <- dir(pattern="*\\.txt$")
#' RTS <- readTSGE(files)
#' PRTS$tsTable # the gene expression raw count table
#' PRTS <- PreparingTheProfiles(RTS)
#' PRTS$scaled # the gene expression standardized count table
#' PRTS$profiles # the profiles of gene expression
#' }
#'
#'@seealso \code{\link{scale}}, \code{\link{MadScale}}, \code{\link{index}}
#'
#' @export
#'
#'
PreparingTheProfiles = function(x,cutoff=1,mad.scale=TRUE){

  if(class(x)!="ctsGEList") stop ( "data must be a ctsGEList object")

  tp <- x$timePoints
  x$scaled <- tmp <- MadScale(x$tsTable)
  if(!mad.scale) x$scaled <- tmp <- t(scale(t(x$tsTable)))
  idx <- data.frame(t(apply(tmp,1,index,cutoff=cutoff))); colnames(idx) <- x$samples
  x$profiles <- cbind(as.data.frame(idx),profiles=apply(idx,1,paste,collapse=""))
  x$cutoff <- cutoff
  structure(x,class = "ctsGEList")

}


#' Indexing step for Time Series Supervised Clustering
#'
#' Takes a numeric vector and return a profile vector that represent the gene expression pattern
#' @param x A numeric
#' @param cutoff A numeric, dermine the threshold for indexing, values between \emph{cutoff and -cutoff} will get the value \bold{0}, meaning no significant change was made in gene expression. Default = 1
#' @return Gene expression profile
#' @seealso \code{\link{PreparingTheProfiles}}
#' @examples
#' rawCounts <- c(103.5, 75.1, 97.3, 27.12, 34.83, 35.53, 40.59, 30.84, 16.39, 29.29)
#'
#' (msCounts <-MadScale(rawCounts)) # standardized step with MadScale
#'
#' cutoff <- seq(0.2,2,0.1)    # different cutoff produce different profiles
#'
#' for(i in cutoff){print(index(msCounts,i))}
#' @export
#'

index=function(x,cutoff=1){ return(unlist(lapply(x,function(x){
  if(x<=(-cutoff)){y=-1}
  if(x>=cutoff){y=1}
  if(x>(-cutoff)&&x<(cutoff)){y=0}
  return(y) })))}

#' Scaling by the Median Absolute Deviation
#'
#' MadScale is a function whose default centers and scales the rows of a numeric matrix,
#' this scaling is made for the purpse of time series profiling. for more information see details
#' @param x A numeric
#' @return The MadScale value of the input
#' @details Centering is done by subtracting the row median.
#'
#' Scaling is done by dividing the (centered) rows of x by their \code{\link{mad}} (median standard deviations)
#'
#'
#' \tabular{r}{
#' \deqn{x - median(x)/mad(x))}\cr}
#'
#'
#' @seealso \code{\link{scale}},  \code{\link{mad}}
#' @examples
#' (rawCounts <- matrix(c(103.5, 75.1, 97.3, 27.12, 34.83, 35.53, 40.59, 30.84, 16.39, 29.29),2,5))
#'
#' (msCounts <-MadScale(rawCounts))
#' @export
#'
MadScale = function(x){

  if(is.matrix(x)){
    tmp <- do.call("rbind",lapply(as.data.frame(t(x)),function(y){(y-median(y))/mad(y)}))
    rownames(tmp) <- rownames(x)
    colnames(tmp) <- colnames(x)
    x <- tmp
  } else (x-median(x))/mad(x)
  }

#' Clustering the profiles
#'
#' Clustering each profile, that was predifined by \code{\link{PreparingTheProfiles}}, with \code{\link{kmeans}}.
#'
#' @param x  A ctsGEList object
#' @param scaling Boolean parameter, does the data should be standardized before clustered. Default = TRUE
#'
#'
#' @return  ctsGEList object is returned as output, with the relative culstered profiles table in \emph{object$ClusteredProfilesTable}, and the number of clusters for each profile in \emph{object$optimalK}
#'
#'
#' @details In order to choose the right k(number of clusters), so that the clusters are compact and well separated, we compute ratio of the within cluster sum of squares(WSS) to the total sum of squares(TSS).
#'
#' when: \eqn{WSS / TSS < 0.2} this is the optipal k
#'
#' It is often recommanded to standardized the data before clustering, the standardize method is the one you chose in \code{\link{PreparingTheProfiles}}.
#'
#' If you are not intresting in standardizing your data before clustering, you can set the \command{\strong{scaling}} parameter to FALSE.
#'
#' If you want to cluster a specific profile, you don't have to cluster all your profiles, you can skip this step and directly run \emph{PlotProfilesClust}
#' with the profile of your choise, for more details see: \code{\link{PlotProfilesClust}}
#'
#' @examples
#' \dontrun{
#'
#'  tsCP <- ClustProfiles(PRTS)
#'  head(tsPS$ClusteredProfilesTable) #the table with the clusterd profiles
#'  head(tsPS$optimalK) #the table with the number of clusters for each profile
#' }
#' @seealso \code{\link{scale}},  \code{\link{kmeans}}, \code{\link{PlotProfilesClust}}
#' @export
#'
ClustProfiles = function(x,scaling=TRUE){

  clust_tbl <-  list()
  kprof <-  list()

  if(scaling){ tmp <- cbind(as.data.frame(x$scaled),profiles=x$profiles[rownames(x$scaled),ncol(x$profiles)])
  }
  else { tmp <- cbind(as.data.frame(x$tsTable),profiles=x$profiles[rownames(x$scaled),ncol(x$profiles)])}


  for(idx in tmp$profiles){
    tbl <-  tmp[tmp$profiles==idx,1:x$timePoints]
    if(nrow(tbl) > 9){
      ratio_ss <-  rep(0,9)
      fit_km <-  list()
      for (k in 1:9) {
        # Apply k-means to school_result: school_km
        fit_km[[k]] <- kmeans(tbl,k)
        # Save the ratio between of WSS to TSS in kth element of ratio_ss
        ratio_ss[k] <-  fit_km[[k]]$tot.withinss/fit_km[[k]]$totss
      }
      G <-  which(ratio_ss<0.2)[1]#when: WSS / TSS < 0.2 this is the optipal k
      clust <-  fit_km[[G]]$cluster
    }

    else{G<- 1;clust<- rep(1,nrow(tbl))}

    kprof[[idx]] <-  c(nrow(tbl),G)
    clust_tbl[[idx]] <- data.frame(clusters=clust,profiles=idx)
    rownames(clust_tbl[[idx]]) <- rownames(tbl)
  }

  names(clust_tbl) <-  NULL
  x$optimalK <-  do.call("rbind",kprof); colnames(x$optimalK) <-  c("tags","k")
  x$ClusteredProfilesTable <- do.call("rbind",clust_tbl)

  structure(x,class = "ctsGEList")

}

#' A graphic visualization of a profile
#'
#' Line graphs made with \code{\link{ggplot}} of a profile and its clusters - the user choose the profile to show,
#' clustering is done with \code{\link{kmeans}}
#'
#' @param x ctsGEList object
#' @param idx A character, the profile to plot (e.g., for 8 time points "11100-1-1-1")
#' @param k A numeric, number of clusters. If not given the function will calculate what is the optimal k for the profile.
#' @param scaling A boolean, default to TRUE, does the data should be standardized before clustered.
#'
#' @return A list with two objects: \enumerate{\item Table of the profile and its clusters
#'                                             \item Gene expression pattern graphs for each one of the clusters}
#'
#' @seealso \code{\link{ggplot}},  \code{\link{kmeans}}, \code{\link{ClustProfiles}}
#'
#' @examples
#' \dontrun{
#'
#'  pp <- PlotProfilesClust(PRTS,idx="000011") # when 'k' is not provided, the function will determine the number of clusters
#'  pp$graphs # plots the line graphs
#'  pp$Clust_k_00011 # the table of genes in the profile seperated to clusters
#' }
#' @export
#'
PlotProfilesClust = function(x,idx,k=NULL,scaling=TRUE){
  plot <- list()
  ggplot_list <- list()

  if(!scaling){tbl <- x$tsTable[rownames(x$profiles[x$profiles[,"profiles"]==idx,]),]
  } else{tbl <- x$scaled[rownames(x$profiles[x$profiles[,"profiles"]==idx,]),]}

  tbl <- as.data.frame(tbl)

  if(!is.null(k)){
    fit_km <- kmeans(tbl,k)
    tbl[names(fit_km$cluster),"clusters"] <- fit_km$cluster
  }else{
    if(length(x) > 7){
      k <- x$optimalK[rownames(x$optimalK)==idx,"k"]
      tags <- rownames(x$ClusteredProfilesTable[x$ClusteredProfilesTable[,"profiles"]==idx,])
      tbl[tags,"clusters"] <- x$ClusteredProfilesTable[x$ClusteredProfilesTable[,"profiles"]==idx,"clusters"]
    }else{

      if(nrow(tbl) > 9){
        ratio_ss <-  rep(0,9)
        fit_km <-  list()
        for (k in 1:9) {
          # Apply k-means
          fit_km[[k]] <- kmeans(tbl,k)
          # Save the ratio between of WSS to TSS in kth element of ratio_ss
          ratio_ss[k] <-  fit_km[[k]]$tot.withinss/fit_km[[k]]$totss
        }
        k <-  which(ratio_ss<0.2)[1]#when: WSS / TSS < 0.2 this is the optipal k
        if(is.na(k)) k <- 9

      }else{k <- 1}

      fit_km <- kmeans(tbl,k)
      tbl[names(fit_km$cluster),"clusters"] <- fit_km$cluster
    }
  }

  tmp <- cbind(tags=rownames(tbl),tbl)
  for(i in 1:k){
    x.m <- reshape2::melt(tmp[tmp$cluster==i,c("tags",x$samples)])

    gplot <- ggplot2::ggplot(x.m, ggplot2::aes(y=value,x=variable,colour=factor(tags),group=factor(tags)))+
      ggplot2::labs(title =paste0("Profile: ",idx," Cluster: ",i) , x = "Time", fill= NULL, y = "Expression Value") +
      ggplot2::theme(legend.position="none")

    gplot <- gplot +ggplot2::geom_line()

    name <- paste0("Clust_",i,"_",idx)
    ggplot_list[[name]] <- gplot

    }

    plot[[name]] <- tbl
    plot$graphs <- ggplot_list

  structure(plot,class = "list")}





