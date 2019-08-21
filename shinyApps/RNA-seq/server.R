library(shiny)
library(DT)
library(png)
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
    }, height=400, width=430)
    
    output$somiteNumber <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("somiteNumber.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 400, height = 70)
    }, deleteFile = TRUE)
    
    output$somiteAgeVertical <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("somiteAge_vertical.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 90, height = 500)
    }, deleteFile = TRUE)
    
    output$downloadGeneExpr <- downloadHandler(
      filename = function() { 
        paste0(input$gene, ".pdf")
      },
      content = function(file) {
        pdf(file, width = 5, height = 5)
        print(boxplotExpr(group_by=input$group, colour_by=input$colour, gene=input$gene))
        dev.off()
      }
    )
    
    #### somite trios
    ## MAplot
    output$MAplotTrios <- renderPlot({
      if(input$levelTrios == 1){ plotMAtriosAll(contrast=input$contrastTriosAll) }
      else{ plotMAtriosStage(contrast=input$contrastTriosStage) }
    }, height=275, width=275)
    
    ## DE results table
    output$DEtableTrios <- DT::renderDataTable(
      if(input$levelTrios == 1){ 
        datatable( printDEtableTriosAll(contrast=input$contrastTriosAll), colnames = c("gene", "logFC", "FDR"), rownames= FALSE )
      }else{
        datatable( printDEtableTriosStage(contrast=input$contrastTriosStage), colnames = c("gene", "logFC-IvsII", "logFC-IvsIII", "logFC-IIvsIII", "FDR"), rownames= FALSE )
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
          boxplotExprTriosAll(group_by=input$groupTrios1, colour_by=input$colourTrios1, contrast=input$contrastTriosAll, selected=selDEtrios)
        }else{
          boxplotExprTriosStage(group_by=input$groupTrios2, colour_by=input$colourTrios2, contrast=input$contrastTriosStage, selected=selDEtrios)
        }
      }
    }, height=400, width=430)
    
    output$downloadGeneExprTrios <- downloadHandler(
      filename = function() { 
        selDEtrios <- input$DEtableTrios_row_last_clicked
        if(input$levelTrios == 1){ 
          table <- DEtriosAll[[as.numeric(input$contrastTriosAll)]][,-c(3:5)]
          table <- table[table$FDR < 0.05 & abs(table$logFC) > log2(1.5),]
          gene <- as.character(table[selDEtrios,1])
        }else{
          table <- DEtriosStage[[as.numeric(input$contrastTriosStage)]][,-c(5:7)]
          table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
          table <- table[table$stageSpecific==1,]
          gene <- as.character(table[selDEtrios,1])
        }
        paste0(gene, ".pdf")
      },
      content = function(file) {
        selDEtrios <- input$DEtableTrios_row_last_clicked
        pdf(file, width = 5, height = 5)
        if(input$levelTrios == 1){ 
          print(boxplotExprTriosAll(group_by=input$groupTrios1, colour_by=input$colourTrios1, contrast=input$contrastTriosAll, selected=selDEtrios))
        }else{
          print(boxplotExprTriosStage(group_by=input$groupTrios2, colour_by=input$colourTrios2, contrast=input$contrastTriosStage, selected=selDEtrios))
        }
        dev.off()
      }
    )
    
    output$somiteAge <- renderImage({
        outfile <- tempfile(fileext='.png')
        p <- readPNG("somiteAge.png")
        writePNG(p, target=outfile)
        list(src = outfile, contentType = 'image/png', width = 450, height = 70)
    }, deleteFile = TRUE)
    
    ## GO enrichment resuts
    output$GOenrichmentTableTrios <- DT::renderDataTable(
      datatable( printGOtableTrios(level=input$levelTriosGO), options = list(pageLength = 25), rownames= FALSE),
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
    }, height=400, width=430)
    
    output$somiteAgeGO <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("somiteAge.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 450, height = 70)
    }, deleteFile = TRUE)
    
    retrieveGOlistTriosSel <- reactive({
      selGOtrios <- input$GOenrichmentTableTrios_row_last_clicked
      getGOgenesTrios(level=input$levelTriosGO, selected=selGOtrios)
    })
    output$downloadGOgeneListTrios <- downloadHandler(
      filename = function() {
        table <- GOresultsTrios[GOresultsTrios$Fisher.classic.adj < 0.1,-6]
        selGOtrios <- input$GOenrichmentTableTrios_row_last_clicked
        go <- table[selGOtrios,1]
        term <- table[selGOtrios,2]
        paste0("DEGs_inGOterm_", go, "_", term, ".tsv") },
      content = function(file) {
        write.table(retrieveGOlistTriosSel(), file, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
      }
    )
    
    output$downloadGOtriosPlot <- downloadHandler(
      filename = function() { 
        gene=input$GOgenesTriosSel
        paste0(gene, ".pdf")
      },
      content = function(file) {
        selDEtrios <- input$DEtableTrios_row_last_clicked
        pdf(file, width = 5, height = 5)
        print(boxplotExprTriosGO(group_by=input$groupTriosGO, colour_by=input$colourTriosGO, gene=input$GOgenesTriosSel))
        dev.off()
      }
    )
    
    #### across development
    ## DE results table
    output$DEtableStage <- DT::renderDataTable(
      if(input$levelStage == 1){ 
        datatable( printDEtableStageAll(), colnames = c("gene", "logCPM", "logFC.max", "FDR"), rownames = FALSE )
      }else{
        datatable( printDEtableStageSomite(contrast=input$somiteStages), colnames = c("gene", "logCPM", "logFC", "FDR"), rownames = FALSE )
      }
      , server = TRUE)
    
    retrieveDEtableStage <- reactive({
      if(input$levelStage == 1){ printDEtableStageAll() }
      else{ printDEtableStageSomite(contrast=input$somiteStages) }
    })
    output$downloadDEtableStage <- downloadHandler(
      filename = function() { 
        if(input$levelStage == 1){ "DEgenes_acrossDevelopment_overall.tsv" }
        else{ paste0("DEgenes_acrossDevelopment_", names(DEstageSomite)[as.numeric(input$somiteStages)], ".tsv") }
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
          boxplotExprStageSomite(group_by=input$groupStage2, colour_by=input$colourStage2, contrast=input$somiteStages, selected=selDEstage)
        }
      }
    }, height=400, width=430)
    
    output$downloadGeneExprStage <- downloadHandler(
      filename = function() { 
        selDEstage <- input$DEtableStage_row_last_clicked
        if(input$levelStage == 1){ 
          table <- DEstageAll[,c(1,17,21,20)]
          table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
          gene <- as.character(table[selDEstage,1])
        }else{
          table <- DEstageSomite[[as.numeric(contrast)]][,c(1,17,22,20,23)]
          table <- table[table$FDR < 0.05 & abs(table$logFC.max) > log2(1.5),]
          table <- table[table$somiteSpecific==1,]
          gene <- as.character(table[selDEstage,1])
        }
        paste0(gene, ".pdf")
      },
      content = function(file) {
        selDEstage <- input$DEtableStage_row_last_clicked
        pdf(file, width = 5, height = 5)
        if(input$levelStage == 1){ 
          print(boxplotExprStageAll(group_by=input$groupStage1, colour_by=input$colourStage1, selected=selDEstage))
        }else{
          print(boxplotExprStageSomite(group_by=input$groupStage2, colour_by=input$colourStage2, contrast=input$somiteStages, selected=selDEstage))
        }
        dev.off()
      }
    )
    
    output$somiteNumberStage <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("somiteNumber.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 400, height = 70)
    }, deleteFile = TRUE)
    
    output$somiteAgeStage <- renderImage({
      outfile <- tempfile(fileext='.png')
      p <- readPNG("somiteAge.png")
      writePNG(p, target=outfile)
      list(src = outfile, contentType = 'image/png', width = 450, height = 70)
    }, deleteFile = TRUE)
    
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
        pdf(file, width = 5, height = 5)
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



