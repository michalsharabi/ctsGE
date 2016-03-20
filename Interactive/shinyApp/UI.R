           tabPanel("Profiles Plots and Table",
                    
                    theme = shinytheme("cerulean"),
                    pageWithSidebar(
                      headerPanel("Shinyapp example for ctsGE package"),
                      
                      sidebarPanel(width = 2,
                                   
                                   
                                   radioButtons("profChoose","Profile Method", 
                                                c("Ready Profile" = "profile",
                                                  "Build Your Own Profile" = "OwnProf"),
                                                selected = "profile"),
                                   
                                   br(),
                                   uiOutput("profiles"),
                                   hr(),
                                   
                                   sliderInput("n", "Number of clusters",min = 1,max = 9,value = 1,step = 1),
                                   checkboxInput("scale", "Raw Counts",value = FALSE)
                                   
                                   
                                   
                      ),
                      
                      
                      mainPanel(width = 10,
                                tabsetPanel(
                                  tabPanel( "Time series", icon = icon("line-chart"),uiOutput("plots")),
                                  tabPanel( "Profiles Table",icon = icon("table"),uiOutput("table"))
                                  
                                )
                      )
                    )#pageWithSidebar
           )#tabpanel
           
           
           
           




