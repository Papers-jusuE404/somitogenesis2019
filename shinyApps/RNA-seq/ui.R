library(shiny)
shinyUI(fluidPage(
  titlePanel(h4("Differential expression between mouse somites")),
  tabsetPanel(
    tabPanel("Gene expression",
      sidebarLayout(
        sidebarPanel(width=3,
          fluidRow(
            column(4,
              radioButtons("group", "plot by:", choices = list("somite"=3, "stage"=2, "somiteNumber"=4), selected = 4)
            ),
            column(4,
              radioButtons("colour", "color by:", choices = list("stage"=2, "somite"=3), selected = 3)
            )
          ),
          fluidRow(
            div(style = "height:185px;"),
            downloadButton('downloadGeneExpr', 'Plot')
          )
        ),
        mainPanel(
          fluidRow(
            column(5,
              selectizeInput("gene", label = h5("Gene of interest"), choices = NULL, options = list(placeholder = 'type a gene name')),
              plotOutput("plotGeneExpr"),
              imageOutput("somiteNumber")
            ),
            column(2,
              div(style = "margin-top:18em",
                div(style = "margin-left:2em",
                  imageOutput("somiteAgeVertical")
              ))
            ),
            column(2,
              h5("Differential expression"),
              div(style = "margin-left:1.5em",
                  h6("Somite trios")
              ),
              plotOutput("summaryDEtrios"),
              div(style = "margin-top:8.5em",
                div(style = "margin-left:2em",
                  imageOutput("legend1")
                )
              )
            ),
            column(2,
              div(style = "margin-top:2.5em",
                div(style = "margin-left:2.5em",
                    h6("Stages")
                ),
                plotOutput("summaryDEstages"),
                div(style = "margin-top:8em",
                  imageOutput("legend2")
                )
              )
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
                        choices = list("somiteIvsII"=1, "somiteIIvsIII"=2, "somiteIvsIII"=3), selected = 3),
            fluidRow(
              column(4,
                     radioButtons("groupTrios1", "plot by:", choices = list("somite"=3, "stage"=2, "somiteNumber"=4), selected = 3)
              ),
              column(4,
                     radioButtons("colourTrios1", "color by:", choices = list("stage"=2, "somite"=3), selected = 3)
              )
            )
          ),
          conditionalPanel(
            condition = "input.levelTrios == 2",
            selectInput("contrastTriosStage", "Contrast:", 
                        choices = list("stage8"=1, "stage18"=2, "stage21"=3, "stage25"=4, "stage27"=5, "stage35"=6), selected = 1),
            fluidRow(
              column(4,
                     radioButtons("groupTrios2", "plot by:", choices = list("somite"=3, "stage"=2, "somiteNumber"=4), selected = 2)
              ),
              column(4,
                     radioButtons("colourTrios2", "color by:", choices = list("stage"=2, "somite"=3), selected = 3)
              )
            )
          ),
          fluidRow(
            plotOutput("MAplotTrios")
          ),
          div(style = "margin-top:-8em", 
            fluidRow(
              downloadButton('downloadDEtableTrios', 'DE results'),
              downloadButton('downloadGeneExprTrios', 'Plot')
            )
          )
        ),
        mainPanel(
          fluidRow(
            column(7,
              h6('Click on a row to plot gene expression'),
              DT::dataTableOutput("DEtableTrios")
            ),
            column(4,
              div(style = "height:10px;"),
              plotOutput("plotGeneExprTrios"),
              div(style = "height:10px;"),
              div(style = "margin-left:2em",
                  imageOutput("somiteAge")
              )
            )
          )
        )
      )
    ),
    tabPanel("Somite trios - GO",
      sidebarLayout(
        sidebarPanel(width=2,
          radioButtons("groupTriosGO", "plot by:", choices = list("somiteNumber"=4, "somite"=3, "stage"=2), selected = 3),
          radioButtons("colourTriosGO", "color by:", choices = list("stage"=2, "somite"=3), selected = 3),
          downloadButton('downloadTriosGOList', 'GO enrichments'),
          downloadButton('downloadGOgeneListTrios', 'Genes in GO term'),
          downloadButton('downloadGOtriosPlot', 'Plot')
        ),
        mainPanel(
          fluidRow(
            column(8,
              div(DT::dataTableOutput("GOenrichmentTableTrios"), style = "font-size:80%")
            ),
            column(4,
              h6('Click on a GO term to retrieve the DE genes associated with it. Select genes to plot.'),
              uiOutput("GOgenesTrios"),
              plotOutput("plotGeneExprTriosGO"),
              div(style = "margin-top:-2em", 
                  div(style = "margin-left:2em",
                  imageOutput("somiteAgeGO")
              ))
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
                radioButtons("groupStage1", "plot by:", choices = list( "stage"=2, "somite"=3, "somiteNumber"=4), selected = 2)
              ),
              column(4,
                radioButtons("colourStage1", "color by:", choices = list("stage"=2, "somite"=3), selected = 2)
              )
            )
          ),
          conditionalPanel(
            condition = "input.levelStage == 2",
            selectInput("somiteStages", "Comparison:", choices = list("somiteI"=1, "somiteII"=2, "somiteIII"=3), selected = 1),
            fluidRow(
              column(4,
                radioButtons("groupStage2", "plot by:", choices = list("stage"=2, "somite"=3, "somiteNumber"=4), selected = 3)
              ),
              column(4,
                radioButtons("colourStage2", "color by:", choices = list("stage"=2, "somite"=3), selected = 2)
              )
            )
          ),
          downloadButton('downloadDEtableStage', 'DE results'),
          downloadButton('downloadGeneExprStage', 'Plot')
        ),
        mainPanel(
          fluidRow(
            column(6,
              h6('Click on a row to plot gene expression'),
              DT::dataTableOutput("DEtableStage")
            ),
            column(4,
              plotOutput("plotGeneExprStage"),
              div(style = "height:10px;"),
              conditionalPanel(
                condition = "input.levelStage == 1",
                imageOutput("somiteNumberStage")
              ),
              conditionalPanel(
                condition = "input.levelStage == 2",
                div(style = "margin-left:2.5em",
                    imageOutput("somiteAgeStage")
                )
              )
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

