
MadScale = function(x){
  
  if(is.matrix(x)){
    tmp <- do.call("rbind",lapply(as.data.frame(t(x)),function(y){(y-median(y))/mad(y)}))
    rownames(tmp) <- rownames(x)
    colnames(tmp) <- colnames(x)
    x <- tmp
  } else (x-median(x))/mad(x)
}

clusters <- function(tbl,g){
  set.seed(100)
  tmp <- as.matrix(tbl[,c(1:rts$timePoints)])
  fit <- kmeans(tmp, g)
  kmeans.groups <-  cbind(merge(data.frame(fit$cluster),tmp,by="row.names",all=TRUE),profiles=tbl$profiles)
  colnames(kmeans.groups)[1:2] <- c("genes","clusters")
  
  return(kmeans.groups)
  
}


get_plot_output_list <- function(tbl,input_n) {
  set.seed(100)
  
  # Insert plot output objects the list
  plot_output_list <- lapply(1:input_n, function(i) {
    
    plotname <- paste("OwnProfiles_clustNo", i, sep="_")
    plot_output_object <- plotOutput(plotname, height = 280, width = 250)
    plot_output_object <- renderPlot({
      
      tbl.m <- melt(tbl[tbl$clusters==i,c("genes",prts$samples)])
      gplot <- ggplot(tbl.m, aes(y=value,x=variable,colour=factor(genes),group=factor(genes)))+ 
        scale_colour_discrete(name = "genes")+ 
        labs(title =paste("Clust",i,sep=" ") , x = "Time", fill= NULL, y = "Genes Expression") +
        
        theme(legend.position="none")  +geom_line()
      
      
      print(gplot)
      
    })
  })
  
  return(do.call(tagList,plot_output_list))
  
}


get_plot_output_list1 <- function(profile,input_n,scaling=1) {
  set.seed(100)
  
  #pp <- PlotProfilesClust(prts,profile,k=input_n,scaling = FALSE)
  if(!scaling) {
    gplot <- PlotProfilesClust(prts,profile,k=input_n)
  }else{
    gplot <- PlotProfilesClust(prts,profile,input_n,scaling = FALSE)}
  
  # Insert plot output objects the list
  plot_output_list <- lapply(1:input_n, function(i) {
    
    
    plotname <- paste("plot",profile, i, sep="_")
    plot_output_object <- plotOutput(plotname, height = 280, width = 250)
    plot_output_object <- renderPlot({
      
      gg <- gplot[[2]][[i]] 
        
        
        print(gg)
    })
  })
  
  return(do.call(tagList,plot_output_list))
  
}

################################################################################## Server #####################################################################################

shinyServer(function(input, output,session) {
  
  output$profiles <- renderUI({
    
    switch(input$profChoose,
           "profile" = return(
             selectInput("profile","Select a Profile:", choices = idx,selected = idx[1])),
           "OwnProf" = return(
             splitLayout(lapply(1:rts$timePoints, function(i){selectInput(paste0("t",i), paste0("time point ",i), c(1,-1,0),width = '25%',selected = c(1,0),multiple = T)})))
           
           
    )
  })
  
  
  # filter input$profile
  filtered <- reactive({
    if (is.null(input$profile)) {
      return(NULL)
    }    
    
    #tbl <- PlotProfilesClust(prts,input$profile)[[1]]
    
    if(!input$scale){
      tbl <- PlotProfilesClust(prts,input$profile,k = input$n)[[1]]
      tbl <- cbind(genes=rownames(tbl),clusters=tbl$clusters,data.frame(tbl[,prts$samples]),profiles=input$profile)
      rownames(tbl) <- NULL
      tbl
    }else{
      tbl <- PlotProfilesClust(prts,input$profile,k = input$n,scaling = FALSE)[[1]]
      tbl <- cbind(genes=rownames(tbl),clusters=tbl$clusters,data.frame(tbl[,prts$samples]),profiles=input$profile)
      rownames(tbl) <- NULL
      tbl
    }
  })
  
  # filter input$OwnProfile
  filtered2 <- reactive({
    if (is.null(input$t1)) {
      return(NULL)
    }    
    tmp <- dplyr::tbl_df(prts$profiles) %>%  dplyr::mutate(genes=rownames(prts$profiles))
    tbl <-  tmp %>% dplyr::filter(h0%in%input$t1, h6%in%input$t2, h12%in%input$t3, h24%in%input$t4, h48%in%input$t5,
                                  h72%in%input$t6)
    
    
    if(!input$scale){
      tbl <- cbind(data.frame(prts$scaled[tbl$genes,]),profiles=tbl$profiles)
      tbl <- clusters(tbl, input$n)
      tbl %>% dplyr::mutate_each(dplyr::funs(factor),clusters)
    }else{
      tbl <- cbind(data.frame(prts$tsTable[tbl$genes,]),profiles=tbl$profiles)
      tbl <- clusters(tbl, input$n)
      tbl %>% dplyr::mutate_each(dplyr::funs(factor),clusters)
    }
    
  })
  
  
  observe({
    
    if(input$profChoose == "profile"){
      
      output$plots <- renderUI({
        
        get_plot_output_list1(input$profile,input$n, input$scale)
        
      })
      
      output$table <- renderUI({ 
        if (is.null(filtered())) {
          return()}
        output$tmp <- DT::renderDataTable(filtered(),
                                          rownames = FALSE,
                                          filter=list(position = 'top', clear = FALSE),
                                          server = TRUE,
                                          extensions = 'TableTools', 
                                          options = list(
                                            lengthMenu = c(10,50,100,nrow(filtered())),
                                            tableTools = list(sSwfPath = copySWF('www'),
                                                              aButtons = list('copy', 'csv', 'xls')),
                                            dom = 'T<"clear">lfrtip'))
        DT::dataTableOutput("tmp")
        
      })#renderUI
      
      
    }#if input$profile
    
    if(input$profChoose=="OwnProf"){
      
      output$plots <- renderUI({
        if (is.null(filtered2())) {
          return()}
        tbl = filtered2()
        get_plot_output_list(tbl,input$n)
        
      })
      
      output$table <- renderUI({ 
        if (is.null(filtered2())) {
          return()}
        tbl = filtered2()
        output$tmp <- DT::renderDataTable(tbl, 
                                          rownames = FALSE,
                                          filter=list(position = 'top', clear = FALSE),
                                          server = TRUE,
                                          extensions = 'TableTools', 
                                          options = list(
                                            lengthMenu = c(10,50,100,nrow(filtered2())),
                                            tableTools = list(sSwfPath = copySWF('www'),
                                                              aButtons = list('copy', 'csv', 'xls')),
                                            dom = 'T<"clear">lfrtip'))
        DT::dataTableOutput("tmp")
        
      })
      
      
      
    }#if input$OwnProfile
    
  })#end observe
  
  
  
  
})


