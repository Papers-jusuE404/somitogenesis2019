library(shiny)
library(DT)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
source("helper.R")

shinyServer(
  function(input, output, session) {
    updateSelectizeInput(session = session, inputId = 'gene', choices = c(Choose = '', genesAvailable), selected = "Hoxa1", server = TRUE)
    
    #### gene expression
    output$plotGeneExpr <- renderPlot({
      boxplotExpr(group_by=input$group, colour_by=input$colour, gene=input$gene)
    }, height=500, width=500)
    
    output$downloadGeneExpr <- downloadHandler(
      filename = function() { 
        paste0(input$gene, ".pdf")
      },
      content = function(file) {
        pdf(file, width = 6, height = 6)
        print(boxplotExpr(group_by=input$group, colour_by=input$colour, gene=input$gene))
        dev.off()
      }
    )
    
    #### somite trios
    ## MAplot
    output$MAplotTrios <- renderPlot({
      if(input$levelTrios == 1){ plotMAtriosAll(contrast=input$contrastTriosAll) }
      else{ plotMAtriosStage(contrast=input$contrastTriosStage) }
    }, height=320, width=320)
    
    ## DE results table
    output$DEtableTrios <- DT::renderDataTable(
      if(input$levelTrios == 1){ 
        datatable( printDEtableTriosAll(contrast=input$contrastTriosAll), colnames = c("gene", "logFC", "FDR") )
      }else{
        datatable( printDEtableTriosStage(contrast=input$contrastTriosStage), colnames = c("gene", "logFC-IvsII", "logFC-IvsIII", "logFC-IIvsIII", "FDR") )
      }
      , server = TRUE)
    
    retrieveDEtableTrios <- reactive({
      if(input$levelTrios == 1){ printDEtableTriosAll(contrast=input$contrastTriosAll) }
      else{ printDEtableTriosStage(contrast=input$contrastTriosStage) }
    })
    output$downloadDEtableTrios <- downloadHandler(
      filename = function() { 
        if(input$levelTrios == 1){ paste0("DEgenes_", names(DEtriosAll)[as.numeric(input$contrastTriosAll)], ".tsv") }
        else{ paste0("DEgenes_somiteTrios_", names(DEtriosStage)[as.numeric(input$contrastTriosStage)], ".tsv") }
      },
      content = function(file) {
        write.table(retrieveDEtableTrios(), file, quote = FALSE, sep="\t")
      }
    )
    
    ## plot selected gene from DE table
    output$plotGeneExprTrios <- renderPlot({
      selDEtrios <- input$DEtableTrios_row_last_clicked
      if(length(selDEtrios)){
        if(input$levelTrios == 1){ 
          boxplotExprTriosAll(group_by=input$groupTrios, colour_by=input$colourTrios, contrast=input$contrastTriosAll, selected=selDEtrios)
        }else{
          boxplotExprTriosStage(group_by=input$groupTrios, colour_by=input$colourTrios, contrast=input$contrastTriosStage, selected=selDEtrios)
        }
      }
    }, height=500, width=500)
    
    output$downloadGeneExprTrios <- downloadHandler(
      filename = function() { 
        selDEtrios <- input$DEtableTrios_row_last_clicked
        if(input$levelTrios == 1){ 
          table <- DEtriosAll[[as.numeric(input$contrastTriosAll)]][,-c(4:5)]
          table <- table[table$FDR < 0.1,]
          gene <- as.character(table[selDEtrios,1])
        }else{
          table <- DEtriosStage[[as.numeric(input$contrastTriosStage)]][,-c(4:5)]
          table <- table[table$FDR < 0.1,]
          gene <- as.character(table[selDEtrios,1])
        }
        paste0(gene, ".pdf")
      },
      content = function(file) {
        selDEtrios <- input$DEtableTrios_row_last_clicked
        pdf(file, width = 6, height = 6)
        if(input$levelTrios == 1){ 
          print(boxplotExprTriosAll(group_by=input$groupTrios, colour_by=input$colourTrios, contrast=input$contrastTriosAll, selected=selDEtrios))
        }else{
          print(boxplotExprTriosStage(group_by=input$groupTrios, colour_by=input$colourTrios, contrast=input$contrastTriosStage, selected=selDEtrios))
        }
        dev.off()
      }
    )
    
    ## GO enrichment resuts
    output$GOenrichmentTableTrios <- DT::renderDataTable(
      datatable( printGOtableTrios(level=input$levelTriosGO), options = list(pageLength = 25)),
      server=TRUE)
    
    retrieveGOlistTrios <- reactive({
      printGOtableTrios(level=input$levelTriosGO)
    })
    output$downloadTriosGOList <- downloadHandler(
      filename = function() { paste0("GOterms_significant_", names(GOresultsTrios)[as.numeric(input$levelTriosGO)], ".tsv") },
      content = function(file) {
        write.table(retrieveGOlistTrios(), file, quote = FALSE, sep="\t")
      }
    )
    
    ## genes associated with selected GO term
    output$GOgenesTrios <- renderUI({
      selGOtrios <- input$GOenrichmentTableTrios_row_last_clicked
      genes <- getGOgenesTrios(level=input$levelTriosGO, selected=selGOtrios)
      selectInput(inputId = 'GOgenesTriosSel', label = "", choices = c(Choose = '', genes), selectize = TRUE, selected = 1)
    })
    
    output$plotGeneExprTriosGO <- renderPlot({
      boxplotExprTriosGO(group_by=input$groupTriosGO, colour_by=input$colourTriosGO, gene=input$GOgenesTriosSel)
    }, height=500, width=500)
    
    
    retrieveGOlistTriosSel <- reactive({
      selGOtrios <- input$GOenrichmentTableTrios_row_last_clicked
      getGOgenesTrios(level=input$levelTriosGO, selected=selGOtrios)
    })
    output$downloadGOgeneListTrios <- downloadHandler(
      filename = function() {
        table <- GOresultsTrios[[as.numeric(input$levelTriosGO)]]
        table <- table[table$FDR < 0.1,-5]
        selGOtrios <- input$GOenrichmentTableTrios_row_last_clicked
        go <- row.names(table)[selGOtrios]
        paste0("DEGs_inGOterm_", go, ".tsv") },
      content = function(file) {
        write.table(retrieveGOlistTriosSel(), file, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
      }
    )
    
    #### across development
    ## DE results table
    output$DEtableStage <- DT::renderDataTable(
      if(input$levelStage == 1){ 
        datatable( printDEtableStageAll(), colnames = c("gene", "logCPM", "logFC.max", "FDR") )
      }else{
        datatable( printDEtableStageSomite(contrast=input$somiteStages), colnames = c("gene", "logCPM", "logFC", "FDR") )
      }
      , server = TRUE)
    
    retrieveDEtableStage <- reactive({
      if(input$levelStage == 1){ printDEtableStageAll() }
      else{ printDEtableStagePairwise(contrast=input$pairStages) }
    })
    output$downloadDEtableStage <- downloadHandler(
      filename = function() { 
        if(input$levelStage == 1){ "DEgenes_acrossDevelopment_overall.tsv" }
        else{ paste0("DEgenes_acrossDevelopment_", substr(colnames(DEstage)[as.numeric(input$pairStages)+1],7,30), ".tsv") }
      },
      content = function(file) {
        write.table(retrieveDEtableStage(), file, quote = FALSE, sep="\t")
      }
    )
    
    ## plot selected gene from DE table
    output$plotGeneExprStage <- renderPlot({
      selDEstage <- input$DEtableStage_row_last_clicked
      if(length(selDEstage)){
        if(input$levelStage == 1){ 
          boxplotExprStageAll(group_by=input$groupStage1, colour_by=input$colourStage1, selected=selDEstage)
        }else{
          boxplotExprStagePairwise(group_by=input$groupStage2, colour_by=input$colourStage2, contrast=input$somiteStages, selected=selDEstage)
        }
      }
    }, height=500, width=500)
    
    output$downloadGeneExprStage <- downloadHandler(
      filename = function() { 
        selDEstage <- input$DEtableStage_row_last_clicked
        if(input$levelStage == 1){ 
          table <- DEstage[,c(1,17,20)]
          table <- table[table$FDR < 0.1,]
          gene <- as.character(table[selDEstage,1])
        }else{
          col <- as.numeric(contrast)+1
          table <- DEstage[,c(1,17,col,20)]
          table <- table[table$FDR < 0.1,]
          gene <- as.character(table[selDEstage,1])
        }
        paste0(gene, ".pdf")
      },
      content = function(file) {
        selDEstage <- input$DEtableStage_row_last_clicked
        pdf(file, width = 6, height = 6)
        if(input$levelStage == 1){ 
          print(boxplotExprStageAll(group_by=input$groupStage, colour_by=input$colourStage, selected=selDEstage))
        }else{
          print(boxplotExprStagePairwise(group_by=input$groupStage, colour_by=input$colourStage, contrast=input$pairStages, selected=selDEstage))
        }
        dev.off()
      }
    )
    
    ## GO enrichment resuts
    output$GOenrichmentTableStage <- DT::renderDataTable(
      datatable( printGOtableStage(), options = list(pageLength = 25)),
      server=TRUE)
    
    retrieveGOlistStage <- reactive({
      printGOtableStage()
    })
    output$downloadStageGOList <- downloadHandler(
      filename = function() { "GOterms_significant_acrossDevelopment.tsv" },
      content = function(file) {
        write.table(retrieveGOlistStage(), file, quote = FALSE, sep="\t")
      }
    )
    
    ## genes associated with selected GO term
    output$GOgenesStage <- renderUI({
      selGOstage <- input$GOenrichmentTableStage_row_last_clicked
      genes <- getGOgenesStage(selected=selGOstage)
      selectInput(inputId = 'GOgenesStageSel', label = "", choices = c(Choose = '', genes), selectize = TRUE, selected = 1)
    })
    
    output$plotGeneExprStageGO <- renderPlot({
      boxplotExprStageGO(group_by=input$groupStageGO, colour_by=input$colourStageGO, gene=input$GOgenesStageSel)
    }, height=500, width=500)
    
    
    retrieveGOlistStageSel <- reactive({
      selGOstage <- input$GOenrichmentTableStage_row_last_clicked
      getGOgenesStage(selected=selGOstage)
    })
    output$downloadGOgeneListStage <- downloadHandler(
      filename = function() {
        table <- GOresultsStage
        table <- table[table$FDR < 0.1,-5]
        selGOstage <- input$GOenrichmentTableStage_row_last_clicked
        go <- row.names(table)[selGOstage]
        paste0("DEGs_inGOterm_", go, ".tsv") },
      content = function(file) {
        write.table(retrieveGOlistStageSel(), file, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
      }
    )
    
    ## Hox correlated genes
    output$HoxCorrGenes <- DT::renderDataTable(
      datatable( printHoxCorrGenes(), options = list(pageLength = 15)), server=TRUE)
    
    output$HoxCorrGenesPlot <- renderPlot({
      selHox <- input$HoxCorrGenes_row_last_clicked
      if(length(selHox)){ plotHoxCorrGenes(idx=selHox) }
    }, height=500, width=500)
    
    output$downloadHoxPlot <- downloadHandler(
      filename = function() { 
        selHox <- input$HoxCorrGenes_row_last_clicked
        paste0(normCounts[hox[selHox,1],1], "vs", normCounts[hox[selHox,2],1], ".pdf")
      },
      content = function(file) {
        selHox <- input$HoxCorrGenes_row_last_clicked
        pdf(file, width = 6, height = 6)
        print(plotHoxCorrGenes(idx=selHox))
        dev.off()
      }
    )
    
    ## GO enrichment resuts - Hox
    output$GOenrichmentTableHox <- DT::renderDataTable(
      datatable( printGOtableHox(), options = list(pageLength = 25)),
      server=TRUE)
    
    retrieveGOlistHox <- reactive({
      printGOtableHox()
    })
    output$downloadHoxGOList <- downloadHandler(
      filename = function() { "GOterms_significant_HoxCorrelatedGenes.tsv" },
      content = function(file) {
        write.table(retrieveGOlistHox(), file, quote = FALSE, sep="\t")
      }
    )
    
    ## genes associated with selected GO term - HOX
    output$GOgenesHox <- renderUI({
      selGOhox <- input$GOenrichmentTableHox_row_last_clicked
      genes <- getGOgenesHox(selected=selGOhox)
      selectInput(inputId = 'GOgenesHoxSel', label = "", choices = c(Choose = '', genes), selectize = TRUE, selected = 1)
    })
    
    output$HoxCorrGenesPlotGO <- renderPlot({
      boxplotExprHoxGO(gene=input$GOgenesHoxSel)
    }, height=500, width=500)
    
    
    retrieveGOlistHoxSel <- reactive({
      selGOhox <- input$GOenrichmentTableHox_row_last_clicked
      getGOgenesHox(selected=selGOhox)
    })
    output$downloadGOgeneListHox <- downloadHandler(
      filename = function() {
        table <- GOresultsHox
        table <- table[table$FDR < 0.1,-5]
        selGOhox <- input$GOenrichmentTableHox_row_last_clicked
        go <- row.names(table)[selGOhox]
        paste0("DEGs_inGOterm_", go, ".tsv") },
      content = function(file) {
        write.table(retrieveGOlistHoxSel(), file, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
      }
    )
    
    
  }
)



