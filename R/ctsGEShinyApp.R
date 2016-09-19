
#' GUI for interactive exploration of gene expression data.
#'
#' Produce and launch Shiny app for interactive exploration of gene expression
#'  data.
#' For more information about shiny apps \cite{http://shiny.rstudio.com/}
#'
#' @param rts list of an expression data that made by readTSGE
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
                             desc=as.factor(tbl$desc),
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
