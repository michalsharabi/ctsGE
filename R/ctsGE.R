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
#' @param desc character vector with genes description (annotation),default to
#' NULL
#' Defaults to the file names.
#' @param ... other are passed to read.delim
#'
#' @return A list with four objects:
#'  \enumerate{\item expression matrix
#'             \item samples names
#'             \item tags - genes name
#'             \item timePoints - number of time points}
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
readTSGE = function(files,path = NULL,columns=c(1,2),labels = NULL,desc = NULL,
                    ...){
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
        print(paste0(sum(tmp==0)," Genes were remove"))}

    if(!is.null(desc)){
        x$desc <-  desc
        if(sum(tmp==0)){
            desc <- desc[which(tmp!=0)]}
        x$desc <- as.data.frame(desc)
        rownames(x$desc) <- x$tags
        colnames(x$desc) <- "desc"}




    structure(x,class = "list")
}

