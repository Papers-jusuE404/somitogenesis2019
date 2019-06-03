library(shiny)
shinyUI(fluidPage(
  titlePanel(h4("Differential expression between mouse somites")),
  tabsetPanel(
    tabPanel("Gene expression",
             sidebarLayout(
               sidebarPanel(width=3,
                            fluidRow(
                              column(4,
                                     radioButtons("group", "plot by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2)
                              ),
                              column(4,
                                     radioButtons("colour", "color by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 3)
                              )
                            )
               ),
               mainPanel(
                 fluidRow(
                   column(5,
                          selectizeInput("gene", label = h4("Gene of interest"), choices = NULL, options = list(placeholder = 'type a gene name')),
                          plotOutput("plotGeneExpr"),
                          div(style = "height:185px;"),
                          downloadButton('downloadGeneExpr', 'Download plot')
                   )
                 )
               )
             )
    ),
    tabPanel("Somite trios - DE",
             sidebarLayout(
               sidebarPanel(width=3,
                            radioButtons("levelTrios", "test:", choices = list("all stages"=1, "per-stage"=2), selected = 1, inline = TRUE),
                            conditionalPanel(
                              condition = "input.levelTrios == 1",
                              selectInput("contrastTriosAll", "Contrast:", 
                                          choices = list("somiteIvsII"=1, "somiteIIvsIII"=2, "somiteIvsIII"=3), selected = 3)
                            ),
                            conditionalPanel(
                              condition = "input.levelTrios == 2",
                              selectInput("contrastTriosStage", "Contrast:", 
                                          choices = list("stage8"=1, "stage18"=2, "stage21"=3, "stage25"=4, "stage27"=5, "stage35"=6), selected = 1)
                            ),
                            fluidRow(
                              column(4,
                                     radioButtons("groupTrios", "plot by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2)
                              ),
                              column(4,
                                     radioButtons("colourTrios", "color by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 3)
                              )
                            ),
                            plotOutput("MAplotTrios")
               ),
               mainPanel(
                 fluidRow(
                   column(6,
                          h6('Click on a row to plot gene expression'),
                          DT::dataTableOutput("DEtableTrios"),
                          downloadButton('downloadDEtableTrios', 'Download DE results')
                   ),
                   column(4,
                          plotOutput("plotGeneExprTrios"),
                          div(style = "height:185px;"),
                          downloadButton('downloadGeneExprTrios', 'Download plot')
                   )
                 )
               )
             )
    ),
    tabPanel("Somite trios - GO",
             sidebarLayout(
               sidebarPanel(width=2,
                            downloadButton('downloadTriosGOList', 'Download GO enrichment'),
                            radioButtons("groupTriosGO", "plot by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2),
                            radioButtons("colourTriosGO", "color by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 3)
               ),
               mainPanel(
                 fluidRow(
                   column(7,
                          div(DT::dataTableOutput("GOenrichmentTableTrios"), style = "font-size:80%")
                   ),
                   column(5,
                          h6('Click on a GO term to retrieve the DE genes associated with it. Select genes to plot.'),
                          # DT::dataTableOutput("GOgenes")
                          uiOutput("GOgenesTrios"),
                          plotOutput("plotGeneExprTriosGO"),
                          div(style = "height:150px;"),
                          downloadButton('downloadGOgeneListTrios', 'Download all genes in GO term')
                   )
                 )
               )
             )
    ),
    tabPanel("Developmental changes - DE",
             sidebarLayout(
               sidebarPanel(width=3,
                            radioButtons("levelStage", "results:", choices = list("all somites"=1, "per-somite"=2), selected = 1, inline = TRUE),
                            conditionalPanel(
                              condition = "input.levelStage == 1",
                              fluidRow(
                                column(4,
                                       radioButtons("groupStage1", "plot by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2)
                                ),
                                column(4,
                                       radioButtons("colourStage1", "color by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2)
                                )
                              )
                            ),
                            conditionalPanel(
                              condition = "input.levelStage == 2",
                              selectInput("somiteStages", "Comparison:", choices = list("somiteI"=1, "somiteII"=2, "somiteIII"=3),
                              # selectInput("pairStages", "Comparison:", choices = list("stage8vs18"=1, "stage8vs21"=2, "stage8vs25"=3,
                              #                                                         "stage8vs27"=4, "stage8vs35"=5, "stage18vs21"=6,
                              #                                                         "stage18vs25"=7, "stage18vs27"=8, "stage18vs35"=9,
                              #                                                         "stage21vs25"=10, "stage21vs27"=11, "stage21vs35"=12,
                              #                                                         "stage25vs27"=13, "stage25vs35"=14, "stage27vs35"=15),
                                          selected = 1),
                            fluidRow(
                              column(4,
                                     radioButtons("groupStage2", "plot by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 3)
                              ),
                              column(4,
                                     radioButtons("colourStage2", "color by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2)
                              )
                            )
                          )
               ),
               mainPanel(
                 fluidRow(
                   column(5,
                          h6('Click on a row to plot gene expression'),
                          DT::dataTableOutput("DEtableStage"),
                          downloadButton('downloadDEtableStage', 'Download DE results')
                   ),
                   column(4,
                          plotOutput("plotGeneExprStage"),
                          div(style = "height:185px;"),
                          downloadButton('downloadGeneExprStage', 'Download plot')
                   )
                 )
               )
             )
    )
    # tabPanel("Developmental changes - GO",
    #          sidebarLayout(
    #            sidebarPanel(width=2,
    #                         downloadButton('downloadStageGOList', 'Download GO enrichment'),
    #                         radioButtons("groupStageGO", "plot by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2),
    #                         radioButtons("colourStageGO", "color by:", choices = list("stage"=2, "somite"=3, "date"=4), selected = 2)
    #            ),
    #            mainPanel(
    #              fluidRow(
    #                column(7,
    #                       div(DT::dataTableOutput("GOenrichmentTableStage"), style = "font-size:80%")
    #                ),
    #                column(5,
    #                       h6('Click on a GO term to retrieve the DE genes associated with it. Select genes to plot.'),
    #                       # DT::dataTableOutput("GOgenes")
    #                       uiOutput("GOgenesStage"),
    #                       plotOutput("plotGeneExprStageGO"),
    #                       div(style = "height:150px;"),
    #                       downloadButton('downloadGOgeneListStage', 'Download all genes in GO term')
    #                )
    #              )
    #            )
    #          )
    # )
    # tabPanel("Hox correlated expression",
    #          fluidRow(
    #            column(4,
    #                   h6('Click on a gene pair to plot their expression'),
    #                   DT::dataTableOutput("HoxCorrGenes")
    #            ),
    #            column(5,
    #                   div(style = "height:50px;"),
    #                   plotOutput("HoxCorrGenesPlot"),
    #                   div(style = "height:100px;"),
    #                   downloadButton('downloadHoxPlot', 'Download plot')
    #            )
    #          )
    # ),
    # tabPanel("Hox correlated expression - GO",
    #          fluidRow(
    #            column(7,
    #                   div(DT::dataTableOutput("GOenrichmentTableHox"), style = "font-size:80%"),
    #                   downloadButton('downloadHoxGOList', 'Download GO enrichment')
    #            ),
    #            column(5,
    #                   h6('Click on a GO term to retrieve the DE genes associated with it. Select genes to plot.'),
    #                   uiOutput("GOgenesHox"),
    #                   plotOutput("HoxCorrGenesPlotGO"),
    #                   div(style = "height:150px;"),
    #                   downloadButton('downloadGOgeneListHox', 'Download all genes in GO term')
    #            )
    #          )
    # )
  )
))

