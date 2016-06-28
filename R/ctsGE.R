#' Read and merge a set of files containing gene expression data
#'
#' Reads and merges a set of text files containing normalized gene expression
#' data
#'
#' @param files character vector of filenames,
#' or alternative a named list of tables for each time point.
#' @param path character string giving the directory containing the files.
#' The default is the current working directory.
#' @param columns numeric vector stating which two columns contain the tag
#' names and counts, respectively
#' @param labels character vector giving short names to associate with the
#' libraries.
#' Defaults to the file names.
#' @param ... other are passed to read.delim
#'
#' @return ctsGEList
#'
#' @details As input, the ctsGE package expects normalized expression table,
#' where rows are genes and columns are samples
#' Each file is assumed to contained digital gene expression data for one
#' sample (or library), with transcript or gene identifiers in the first
#' column and expression values in the second column.
#' Transcript identifiers are assumed to be unique and not repeated in any one
#' file.
#' By default, the files are assumed to be tab-delimited and to contain column
#'    headings.
#' The function forms the union of all transcripts and creates one big table
#' with zeros where necessary.
#' When reading the normalized expression values the function check whether
#' there are rows that their median absolute deviation (MAD) value equal to zero
#' and remove these rows. This step is important in order to continue to the
#' next step of indexing the data.
#' The function will output a message of how many genes were remove.
#'
#' @examples
#' ## Read all .txt files from current working directory
#' data_dir <- system.file("extdata", package = "ctsGE")
#' files <- dir(path=data_dir,pattern = "\\.xls$")
#'
#' # reading only 2000 genes
#' rts <- readTSGE(files, path = data_dir,
#'  labels = c("0h","6h","12h","24h","48h","72h"), skip = 10625 )
#' @export
#' @import limma stats utils
#'
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
                stop(paste(
                    "Repeated tag sequences in table no.",i,"in the list"))
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
        if(!is.list(files)){ colnames(x$tsTable) <-
            x$samples <- limma::removeExt(files)
        }else{colnames(x$tsTable) <- x$samples <- names(files)}
    }else{x$samples <- colnames(x$tsTable)}
    for (i in 1:nfiles) {
        aa <- match(taglist[[i]], tags)
        x$tsTable[aa, i] <- as.numeric(d[[i]][, columns[2]])
    }

    x$tags <- tags
    x$timePoints <- nfiles

    ## filtering out genes with low expression
    tmp <- apply(t(x$tsTable),2,mad)
    if(sum(tmp==0)){
        x$tsTable <- x$tsTable[names(which(tmp!=0)),]
        x$tags <- names(which(tmp!=0))
        print(paste0(sum(tmp==0)," Genes were remove"))
    }

    structure(x,class = "ctsGEList")
}

#' Define an expression index for each gene
#'
#' Reads the table of genes expression and return an expression index for each
#' gene.
#'
#'
#' @param x A ctsGEList object
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
    if(class(x)!="ctsGEList") stop ( "data must be a ctsGEList object")
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
    structure(x,class = "ctsGEList")
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

#' Graphic visualization of an index
#'
#' The function generates graphs and tables of a specific index and its
#' clusters. The user decides whether to supply the k or let the function
#' calculate the k for the selected index
#'
#' @param x ctsGEList object
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
        if(length(x) > 7){
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

    plot[[name]] <- tbl
    plot$graphs <- ggplot_list
    structure(plot,class = "list")
}

#' GUI for interactive exploration of gene expression data.
#'
#' Produce and launch Shiny app for interactive exploration of gene expression
#'  data.
#' For more information about shiny apps \cite{http://shiny.rstudio.com/}
#'
#' @param rts ctsGEList object
#' @param cutoff A numeric that define the degree of change in gene espression
#'  rate.  see \code{\link{PreparingTheIndexes}}
#' @param mad.scale A boolean defaulting to TRUE as to what method of scaling
#' to use.
#'  Default median-base scaling. FALSE, mean-base scaling
#' @param title Character, the title at the header panel. default to NULL.
#'
#' @return Creates a shiny application and opens a shinyapp.io web page
#'
#' @details The `ctsGEShinyApp` function takes two arguments the ctsGE object
#' and a cutoff,and opens an html page as a GUI. On the web page, the user
#' chooses the profile to visualize and the number of clusters (k parameter for
#'  K-means) to show. The line graph of the profile separated into the clusters
#'  will show in the main panel, and a list of the genes and their expressions
#'  will also be available. The tables and figures can be downloaded.
#'
#' @seealso shiny::ShinyApp
#' @examples
#'
#' \dontrun{
#' data_dir <- system.file("extdata", package = "ctsGE")
#' files <- dir(path=data_dir,pattern = "\\.xls$")
#' rts <- readTSGE(files, path = data_dir,
#' labels = c("0h","6h","12h","24h","48h","72h") )
#' ctsGEShinyApp(rts)}
#'
#'
#' @export
#' @import  stats shiny
#'
ctsGEShinyApp <- function(rts, cutoff = 1, mad.scale = TRUE,title = NULL) {
    prts <- PreparingTheIndexes(rts, cutoff, mad.scale)
    idx <- as.character(unique(prts$index[,"index"]))

    clusters <- function(tbl,g){
        #set.seed(100)
        tmp <- as.matrix(tbl[,c(1:rts$timePoints)])
        fit <- kmeans(tmp, g,nstart = 25)
        kmeans.groups <-
            cbind(merge(data.frame(fit$cluster),tmp,by="row.names",all=TRUE),
                  index=tbl$index)
        colnames(kmeans.groups)[1:2] <- c("genes","clusters")

        return(kmeans.groups)
    }

    get_plot_output_list <- function(plot_list) {
        #set.seed(100)

        # Insert plot output objects the list
        plot_output_list <- lapply(1:length(plot_list), function(i) {
            plotname <- names(plot_list)[i]
            plot_output_object <-
                shiny::plotOutput(plotname, height = 280, width = 250)
            plot_output_object <- shiny::renderPlot({
                gg <- plot_list[[i]]
                print(gg)
            })
        })
        return(do.call(shiny::tagList,plot_output_list))
    }

    shiny::shinyApp(

        ui = shiny::pageWithSidebar(
            shiny::headerPanel(title),
            shiny::sidebarPanel(width = 2,
                                shiny::selectInput("index","Select an Index:",
                                                   choices = idx,
                                                   selected = idx[1]),
                                shiny::hr(),
                                shiny::sliderInput("n", "Number of clusters",
                                                   min = 1,max= 10,
                                                   value= 1,step= 1),
                                shiny::checkboxInput("scale",
                                                     "Unscaled values",value = FALSE)
            ),
            shiny::mainPanel(width = 10,
                             shiny::tabsetPanel(
                                 shiny::tabPanel("Time series",
                                                 icon =shiny::icon("line-chart"),
                                                 shiny::uiOutput("plots")),
                                 shiny::tabPanel( "Genes Table",
                                                  icon = shiny::icon("table"),
                                                  shiny::uiOutput("table"))
                             )
            )
        ),#pageWithSidebar

        server = function(input, output,session) {
            # filter input$index
            filtered <- shiny::reactive({
                if (is.null(input$index)) {
                    return(NULL)
                }
                set.seed(100)
                if(!input$scale){
                    PlotIndexesClust(prts,input$index,k = input$n)
                    }else{
                        PlotIndexesClust(prts,input$index,k = input$n,
                                         scaling = FALSE)}
            })

            shiny::observe({
                list_plot <- filtered()[[2]]
                output$plots <- shiny::renderUI({

                    get_plot_output_list(list_plot)
                })
                tbl <- filtered()[[1]]
                tbl <- cbind(genes=rownames(tbl),
                            clusters=as.factor(tbl$clusters),
                            data.frame(tbl[,prts$samples]),
                            index=input$index)
                rownames(tbl) <- NULL
                output$table <- shiny::renderUI({
                    if (is.null(tbl)) {return()}

                    output$tmp <-
                        DT::renderDataTable(tbl,
                                            rownames = FALSE,
                                            extensions ='Buttons',
                                            options = list(dom =
                                                            'TB<"clear">lfrtip',
                                                           buttons = c('copy',
                                                                       'csv',
                                                                       'excel',
                                                                       'print'),
                                                           lengthMenu =
                                                            list(c(15,50,100,
                                                                 nrow(tbl)),
                                                                 c('15','50',
                                                                   '100','All')
                                                                 ))

                        )
                    DT::dataTableOutput("tmp")
                })#renderUI
            })#end observe
        }
    )
}
