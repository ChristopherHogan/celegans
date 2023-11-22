library(dplyr)
library(qs)

unreliable_genes = c("dpy-20", "cho-1", "C30A5.16", "eat-4", "saeg-2", "unc-47",
                     "unc-119", "pha-1", "F38B6.2", "lin-15B", "lin-15A",
                     "C30F8.3", "gcy-35", "unc-54", "cex-1", "rol-6", "unc-53")

data_root <- file.path("D:", "sw", "CengenApp")
# note the assignment in the global environment
load_as_needed <- function(dataset){
  if(!exists(dataset)){
    assign(dataset,
           qs::qread(paste0(data_root, "/data/", dataset, ".qs")),
           envir = .GlobalEnv)
  }
}

# load global data
load_as_needed("gene_list")
load_as_needed("all_cell_types")
load_as_needed("all_neuron_types")


utr <- c("WBGene00023498","WBGene00023497","WBGene00004397","WBGene00006843",
         "WBGene00004010","WBGene00006789","WBGene00001135","WBGene00001079",
         "WBGene00001135","WBGene00006783","WBGene00000501","WBGene00006788",
         "WBGene00001555")

whichCellTypes <- function(input) {
  load_as_needed("L4.all.TPM.raw_th")
  load_as_needed("ths")
  
  output = list()
  g <- input$Tgene_name
  # print(g)
  
  var <- ""
  
  if (g %in% gene_list$gene_name) {
    var <- g
  }
  
  if (g %in% gene_list$gene_id) {
    var <- dplyr::filter(gene_list, gene_id == g)$gene_name
  }
  
  if (g %in% gene_list$seqnames) {
    var <- dplyr::filter(gene_list, seqnames == g)$gene_name
  }

  th <- list()
  if( input$Tgene_cut == "All Cells Unfiltered" ) { 
    th = L4.all.TPM.raw_th
    columns <- c(1:169)
  } else { 
    th = ths
    columns <- c(2:129)
  }
    
  if (var %in% dplyr::filter(th, threshold == input$Tgene_cut)$gene_name) {
    t3 <- dplyr::filter( th, gene_name == var, threshold == input$Tgene_cut )[, columns]
    # TODO: Is this sorting necessary since we arrange the data.frame later?
    t3[order(unlist(t3), decreasing = TRUE)]
    t3 <- data.frame(CellType = names(t3), expression = as.numeric(t3))
    t3 <- dplyr::filter(t3, expression > 0)
    if (nrow(t3) > 0) {
      t3[, 2] <-
        as.numeric(formatC(t3[, 2], digits = 3, format = "f") %>% gsub(" ", "", .))
    }
    output <- arrange(t3, desc(expression))
#    colnames(t3) <- c("Cell type", "Expression level")
#    output$text1 <-
#      isolate(renderText({
#        shiny::validate(need(
#          !dplyr::filter(gene_list, gene_name == var)$gene_id %in% utr,
#          message = paste0(
#            "WARNING: ",
#            input$Tgene_name,
#            " expression is unreliable as it has been overexpressed to generate transgenic strains."
#          )
#        ))
#      }))
    #output$text1 <- renderText({""})
  }
  else {
    output <- NULL
#    output$text1 <-
#      renderText({
#        "Gene is not expressed at any threshold or does not exist"
#      })
  }
  
#  output$Tgene_name_table <- t3
#    DT::renderDataTable({
#      DT::datatable(
#        t3,
#        options = list(pageLength = 10, autoWidth = TRUE),
#        rownames = FALSE,
#        style = 'jQueryUI',
#        class = 'cell-border stripe'
#      ) %>% formatStyle(c(1:2), color = "black", backgroundColor = 'white')
#    })
  
#  output$get_download_cell <- renderUI({   
#    req(input$Tgene_name)
#    downloadButton('downloadCell', "Download table")
#  })
  
#  output$downloadCell <-
#    downloadHandler(
#      filename = function() {
#        paste(
#          "GenesExpressing-",
#          input$Tgene_name,
#          "-thrs",
#          input$Tgene_cut,
#          ".csv",
#          sep = ""
#        )
#      },
#      content = function(file) {
#        write.csv(t3, file, dec = ".", sep = "\t")
#      }
#    )
  output
}

expressedGenes <- function(input) {
  load_as_needed("L4.all.TPM.raw_th")
  load_as_needed("ths")
  
  if( input$Tcell_cut == "All Cells Unfiltered" ) {
    th = L4.all.TPM.raw_th
  } else {
    th = ths
  }
  
#  output$Error1 <-
#    isolate(renderText({
#      shiny::validate(need(
#        input$Tcell_name %in% colnames(th),
#        message = paste0(
#          "WARNING: If you want to query non-neuronal cells select All Cells Unfiltered"
#        )
#      ))
#    }))
  
  if (input$Tcell_name %in% colnames(th)) {
    t4 <- dplyr::filter(th, threshold == input$Tcell_cut) %>% dplyr::select(gene_name, input$Tcell_name) %>% dplyr::filter(!gene_name %in% unreliable_genes)
    t4d <- dplyr::filter(th, threshold == input$Tcell_cut) %>% dplyr::select(gene_name, X, input$Tcell_name)
    t4 <- t4[rev(order(t4[, 2])), ]
    t4d <- t4d[rev(order(t4d[, 3])), ]
    t4[, 2] <- as.numeric(formatC(t4[, 2], digits = 3, format = "f") %>% gsub(" ", "", .))
    colnames(t4) <- c("Gene name", "Expression level")
    t4d[, 3] <- as.numeric(formatC(t4d[, 3], digits = 3, format = "f") %>% gsub(" ", "", .))
    colnames(t4d) <- c("Gene name", "Gene ID", "Expression level")
    
#    output$Tcell_name_table <-
#      DT::renderDataTable({
#        DT::datatable(
#          t4,
#          options = list(pageLength = 10, autoWidth = TRUE),
#          rownames = FALSE,
#          style = 'jQueryUI',
#          class = 'cell-border stripe'
#        ) %>% formatStyle(c(1:2), color = "black", background = 'white')
#      })
#    
#    output$get_download_gene <- renderUI({   
#      req(input$TCell)
#      downloadButton('downloadGene', "Download table")
#    })
#    
#    output$downloadGene <-
#      downloadHandler(
#        filename = function() {
#          paste(
#            "GenesExpressed_in_",
#            input$Tcell_name,
#            "-thrs",
#            input$Tcell_cut,
#            ".csv",
#            sep = ""
#          )
#        },
#        content = function(file) {
#          write.csv(t4d , file, dec = ".", sep = "\t")
#        }
#      )
    
  } else {
    # output$Tcell_name_table  <- NULL
    t4 <- NULL
  }
  t4
}

project <- "cholinergic" # or gaba 
one_drive <- gsub("\\\\", "/", Sys.getenv("OneDrive"))
data_dir <- paste0(project, "_data")
data_dir_path <- file.path(one_drive, "research", data_dir)
pc_csv_dir <- file.path(data_dir_path, "csv")

output_fname <- paste0(project, "_analysis.txt")
output_path <- file.path(data_dir_path, output_fname)
file.create(output_path)


for (fname in list.files(pc_csv_dir, full.names = TRUE)) {
  result <- basename(fname)
  # read csv
  pc_df <- read.csv(fname)
  # get top x genes
  n <- 5
  top_n_genes <- pc_df[1:n, 1]
  for (gene in top_n_genes) {
    # print(sprintf("  %s", gene))
    result <- paste(result, sprintf("  %s", gene), sep = "\n")
    gene_input <- list()
    gene_input$Tgene_name <- gene
    gene_input$Tgene_cut <- 2 
    gene_table <- whichCellTypes(gene_input)
    top_n_cell_types <- gene_table[1:n, 1]
    
    for (cell_type in top_n_cell_types) {
      # print(sprintf("    %s", cell_type))
      result <- paste(result, sprintf("    %s", cell_type), sep = "\n")
      cell_input <- list()
      cell_input$Tcell_name <- cell_type
      cell_input$Tcell_cut <- 2
      cell_table <- expressedGenes(cell_input)
      top_n <- cell_table[1:n, 1]
      # print(sprintf("  Top %d expressed genes in %s.", n, cell_input$Tcell_name))
      for (e in top_n) {
        # print(sprintf("      %s", e))
        result <- paste(result, sprintf("      %s", e), sep = "\n")
      }
    }
  }
  cat(result, file = output_path, sep = "\n", append = TRUE)
}





