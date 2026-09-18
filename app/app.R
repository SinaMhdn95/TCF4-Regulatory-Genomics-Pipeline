# ============================================================
# BED Peak and Gene Overlap Analyzer
# Developed by Sina Mahdiani
# Supervisor: Dr. Joseph McClay
# Genome: hg38
# ============================================================

library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomeInfoDb)
library(bsicons)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ============================================================
# Helper functions
# ============================================================

clean_dataset_name <- function(x) {
  x <- tools::file_path_sans_ext(basename(x))
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x
}

read_bed_file <- function(path) {
  bed <- read.table(
    path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )

  bed <- bed[, 1:3]
  colnames(bed) <- c("chr", "start", "end")

  bed$start <- as.numeric(bed$start)
  bed$end <- as.numeric(bed$end)

  bed <- bed[!is.na(bed$chr) & !is.na(bed$start) & !is.na(bed$end), ]
  bed <- bed[bed$end > bed$start, ]

  gr <- GRanges(
    seqnames = bed$chr,
    ranges = IRanges(
      start = bed$start + 1,
      end = bed$end
    )
  )

  seqlevelsStyle(gr) <- "UCSC"
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- sort(gr)

  gr
}

annotate_to_genes <- function(gr) {
  anno <- annotatePeak(
    gr,
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )

  anno_df <- as.data.frame(anno)

  genes <- anno_df$SYMBOL
  genes <- as.character(genes)
  genes <- trimws(genes)
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]

  unique(genes)
}

read_background_genes <- function(path) {
  ext <- tools::file_ext(path)

  sep <- ifelse(ext == "csv", ",", "\t")

  bg <- tryCatch(
    read.table(
      path,
      header = TRUE,
      sep = sep,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      fill = TRUE
    ),
    error = function(e) {
      read.table(
        path,
        header = FALSE,
        sep = sep,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = "",
        fill = TRUE
      )
    }
  )

  possible_cols <- c(
    "Gene", "gene", "GENE",
    "Symbol", "symbol", "SYMBOL",
    "gene_symbol", "GeneSymbol", "gene_name"
  )

  matched_col <- intersect(possible_cols, colnames(bg))

  if (length(matched_col) > 0) {
    genes <- bg[[matched_col[1]]]
  } else {
    genes <- bg[[1]]
  }

  genes <- as.character(genes)
  genes <- trimws(genes)
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]
  genes <- genes[!tolower(genes) %in% c("gene", "symbol", "gene_symbol", "genes")]

  unique(genes)
}

make_gene_fisher <- function(genes1, genes2, background_genes) {

  genes1 <- unique(trimws(as.character(genes1)))
  genes2 <- unique(trimws(as.character(genes2)))
  background_genes <- unique(trimws(as.character(background_genes)))

  genes1 <- genes1[genes1 != "" & !is.na(genes1)]
  genes2 <- genes2[genes2 != "" & !is.na(genes2)]
  background_genes <- background_genes[background_genes != "" & !is.na(background_genes)]

  genes1_bg <- intersect(genes1, background_genes)
  genes2_bg <- intersect(genes2, background_genes)

  a <- length(intersect(genes1_bg, genes2_bg))
  b <- length(setdiff(genes1_bg, genes2_bg))
  c <- length(setdiff(genes2_bg, genes1_bg))
  d <- length(setdiff(background_genes, union(genes1_bg, genes2_bg)))

  if (length(background_genes) == 0 || (a + b + c + d) == 0) {
    return(list(
      overlap = a,
      pvalue = NA,
      odds_ratio = NA,
      genes1_used = length(genes1_bg),
      genes2_used = length(genes2_bg),
      status = "Invalid background"
    ))
  }

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

  fisher <- tryCatch(
    fisher.test(mat),
    error = function(e) NULL
  )

  if (is.null(fisher)) {
    return(list(
      overlap = a,
      pvalue = NA,
      odds_ratio = NA,
      genes1_used = length(genes1_bg),
      genes2_used = length(genes2_bg),
      status = "Fisher failed"
    ))
  }

  list(
    overlap = a,
    pvalue = fisher$p.value,
    odds_ratio = unname(fisher$estimate),
    genes1_used = length(genes1_bg),
    genes2_used = length(genes2_bg),
    status = "OK"
  )
}

# ============================================================
# UI
# ============================================================

ui <- page_sidebar(

  title = div(
    h2("BED Peak and Gene Overlap Analyzer"),
    p("Developed by Sina Mahdiani"),
    p("Supervisor: Dr. Joseph McClay")
  ),

  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    secondary = "#18BC9C"
  ),

  sidebar = sidebar(
    width = 330,

    h4("Input files"),

    fileInput(
      inputId = "bed_files",
      label = "Upload two or more BED files",
      multiple = TRUE,
      accept = c(".bed", ".narrowPeak", ".broadPeak")
    ),

    fileInput(
      inputId = "background_file",
      label = "Optional: upload background gene universe",
      multiple = FALSE,
      accept = c(".csv", ".txt", ".tsv")
    ),

    helpText(
      "Use gene symbols for the background file. If no file is uploaded, the app uses the union of annotated genes from the uploaded BED files."
    ),

    hr(),

    h4("Settings"),
    p("Genome annotation: hg38"),
    p("Gene annotation: nearest gene using ChIPseeker"),

    actionButton(
      inputId = "run_analysis",
      label = "Run analysis",
      class = "btn btn-primary"
    ),

    hr(),

    downloadButton("download_peak_summary", "Download peak summary"),
    br(), br(),
    downloadButton("download_gene_summary", "Download gene Fisher summary"),
    br(), br(),
    downloadButton("download_gene_lists", "Download annotated genes"),
    br(), br(),
    downloadButton("download_overlap_genes", "Download overlap genes")
  ),

  layout_columns(
    col_widths = c(3, 3, 3, 3),
    uiOutput("box_file_count"),
    uiOutput("box_peak_count"),
    uiOutput("box_gene_count"),
    uiOutput("box_background_count")
  ),

  br(),

  navset_card_tab(

    nav_panel(
      "Uploaded files",
      br(),
      DTOutput("uploaded_table")
    ),

    nav_panel(
      "Peak-level overlap",
      br(),
      DTOutput("peak_overlap_table"),
      br(),
      plotOutput("peak_overlap_plot", height = "350px")
    ),

    nav_panel(
      "Gene-level overlap + Fisher test",
      br(),
      DTOutput("gene_overlap_table"),
      br(),
      plotOutput("gene_overlap_plot", height = "350px")
    ),

    nav_panel(
      "Annotated gene lists",
      br(),
      DTOutput("gene_list_table")
    ),

    nav_panel(
      "Overlap gene lists",
      br(),
      DTOutput("overlap_gene_table")
    )
  ),

  br(),
  hr(),
  p(
    "Developed by Sina Mahdiani | Supervisor: Dr. Joseph McClay | VCU School of Pharmacy",
    style = "text-align:center; color: gray;"
  )
)

# ============================================================
# Server
# ============================================================

server <- function(input, output, session) {

  bed_data <- eventReactive(input$run_analysis, {

    req(input$bed_files)

    validate(
      need(nrow(input$bed_files) >= 2, "Please upload at least two BED files.")
    )

    withProgress(message = "Reading BED files...", value = 0, {

      beds <- list()

      clean_names <- clean_dataset_name(input$bed_files$name)
      clean_names <- make.unique(clean_names)

      for (i in seq_len(nrow(input$bed_files))) {
        incProgress(1 / nrow(input$bed_files))

        file_path <- input$bed_files$datapath[i]
        file_name <- clean_names[i]

        beds[[file_name]] <- read_bed_file(file_path)
      }

      beds
    })
  })

  gene_data <- reactive({

    beds <- bed_data()

    withProgress(message = "Annotating peaks to genes...", value = 0, {

      gene_sets <- list()
      gene_table <- data.frame()

      for (nm in names(beds)) {
        incProgress(1 / length(beds))

        genes <- annotate_to_genes(beds[[nm]])

        gene_sets[[nm]] <- genes

        temp <- data.frame(
          Dataset = nm,
          Gene = genes
        )

        gene_table <- dplyr::bind_rows(gene_table, temp)
      }

      list(
        gene_sets = gene_sets,
        gene_table = gene_table
      )
    })
  })

  background_genes <- reactive({

    gene_sets <- gene_data()$gene_sets

    if (!is.null(input$background_file)) {
      bg <- read_background_genes(input$background_file$datapath)
      source <- "Uploaded background gene universe"
    } else {
      bg <- unique(unlist(gene_sets))
      source <- "Union of annotated genes from uploaded BED files"
    }

    list(
      genes = unique(bg),
      source = source
    )
  })

  peak_overlap_results <- reactive({

    beds <- bed_data()
    dataset_names <- names(beds)

    results <- data.frame()

    for (i in seq_along(dataset_names)) {
      for (j in seq_along(dataset_names)) {

        if (i < j) {

          name1 <- dataset_names[i]
          name2 <- dataset_names[j]

          gr1 <- beds[[name1]]
          gr2 <- beds[[name2]]

          hits <- findOverlaps(gr1, gr2)

          overlap_1 <- unique(queryHits(hits))
          overlap_2 <- unique(subjectHits(hits))

          result_row <- data.frame(
            Comparison = paste(name1, "vs", name2),
            Dataset_1 = name1,
            Dataset_2 = name2,
            Dataset_1_Total_Peaks = length(gr1),
            Dataset_2_Total_Peaks = length(gr2),
            Overlapping_Peaks_From_Dataset_1 = length(overlap_1),
            Overlapping_Peaks_From_Dataset_2 = length(overlap_2),
            Percent_of_Dataset_1 = round(length(overlap_1) / length(gr1) * 100, 3),
            Percent_of_Dataset_2 = round(length(overlap_2) / length(gr2) * 100, 3)
          )

          results <- dplyr::bind_rows(results, result_row)
        }
      }
    }

    results
  })

  gene_overlap_results <- reactive({

    gene_sets <- gene_data()$gene_sets
    dataset_names <- names(gene_sets)

    bg <- background_genes()$genes
    bg_source <- background_genes()$source

    results <- data.frame()

    for (i in seq_along(dataset_names)) {
      for (j in seq_along(dataset_names)) {

        if (i < j) {

          name1 <- dataset_names[i]
          name2 <- dataset_names[j]

          genes1 <- gene_sets[[name1]]
          genes2 <- gene_sets[[name2]]

          fisher_result <- make_gene_fisher(
            genes1 = genes1,
            genes2 = genes2,
            background_genes = bg
          )

          result_row <- data.frame(
            Comparison = paste(name1, "vs", name2),
            Dataset_1 = name1,
            Dataset_2 = name2,
            Dataset_1_Total_Genes = length(genes1),
            Dataset_2_Total_Genes = length(genes2),
            Dataset_1_Genes_In_Background = fisher_result$genes1_used,
            Dataset_2_Genes_In_Background = fisher_result$genes2_used,
            Overlap_Genes = fisher_result$overlap,
            Percent_of_Dataset_1 = round(fisher_result$overlap / length(genes1) * 100, 3),
            Percent_of_Dataset_2 = round(fisher_result$overlap / length(genes2) * 100, 3),
            Background_Genes_Used = length(bg),
            Background_Source = bg_source,
            Fisher_P_Value = fisher_result$pvalue,
            Odds_Ratio = fisher_result$odds_ratio,
            Fisher_Status = fisher_result$status
          )

          results <- dplyr::bind_rows(results, result_row)
        }
      }
    }

    results
  })

  overlap_gene_results <- reactive({

    gene_sets <- gene_data()$gene_sets
    dataset_names <- names(gene_sets)

    results <- data.frame()

    for (i in seq_along(dataset_names)) {
      for (j in seq_along(dataset_names)) {

        if (i < j) {

          name1 <- dataset_names[i]
          name2 <- dataset_names[j]

          overlap_genes <- intersect(gene_sets[[name1]], gene_sets[[name2]])

          if (length(overlap_genes) > 0) {
            temp <- data.frame(
              Comparison = paste(name1, "vs", name2),
              Gene = overlap_genes
            )

            results <- dplyr::bind_rows(results, temp)
          }
        }
      }
    }

    results
  })

  # ============================================================
  # Summary boxes
  # ============================================================

  output$box_file_count <- renderUI({
    n <- ifelse(is.null(input$bed_files), 0, nrow(input$bed_files))
    value_box("BED files", n, showcase = bs_icon("files"))
  })

  output$box_peak_count <- renderUI({
    if (input$run_analysis == 0) {
      value_box("Total peaks", "Run analysis", showcase = bs_icon("bar-chart"))
    } else {
      beds <- bed_data()
      value_box("Total peaks", sum(sapply(beds, length)), showcase = bs_icon("bar-chart"))
    }
  })

  output$box_gene_count <- renderUI({
    if (input$run_analysis == 0) {
      value_box("Unique genes", "Run analysis", showcase = bs_icon("diagram-3"))
    } else {
      genes <- unique(gene_data()$gene_table$Gene)
      value_box("Unique genes", length(genes), showcase = bs_icon("diagram-3"))
    }
  })

  output$box_background_count <- renderUI({
    if (input$run_analysis == 0) {
      value_box("Background genes", "Optional", showcase = bs_icon("database"))
    } else {
      bg <- background_genes()$genes
      value_box("Background genes", length(bg), showcase = bs_icon("database"))
    }
  })

  # ============================================================
  # Tables
  # ============================================================

  output$uploaded_table <- renderDT({

    req(input$bed_files)

    dat <- data.frame(
      Dataset_Name = make.unique(clean_dataset_name(input$bed_files$name)),
      Original_File = input$bed_files$name,
      Size_MB = round(input$bed_files$size / 1024 / 1024, 3)
    )

    datatable(dat, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$peak_overlap_table <- renderDT({
    datatable(
      peak_overlap_results(),
      options = list(pageLength = 10, scrollX = TRUE),
      filter = "top"
    )
  })

  output$gene_overlap_table <- renderDT({
    datatable(
      gene_overlap_results(),
      options = list(pageLength = 10, scrollX = TRUE),
      filter = "top"
    )
  })

  output$gene_list_table <- renderDT({
    datatable(
      gene_data()$gene_table,
      options = list(pageLength = 20, scrollX = TRUE),
      filter = "top"
    )
  })

  output$overlap_gene_table <- renderDT({
    datatable(
      overlap_gene_results(),
      options = list(pageLength = 20, scrollX = TRUE),
      filter = "top"
    )
  })

  # ============================================================
  # Plots
  # ============================================================

  output$peak_overlap_plot <- renderPlot({

    dat <- peak_overlap_results()

    ggplot(dat, aes(x = reorder(Comparison, Percent_of_Dataset_1), y = Percent_of_Dataset_1)) +
      geom_col(fill = "#18BC9C") +
      coord_flip() +
      labs(
        x = NULL,
        y = "Percent of Dataset 1 peaks overlapping Dataset 2",
        title = "Peak-level overlap"
      ) +
      theme_minimal(base_size = 13)
  })

  output$gene_overlap_plot <- renderPlot({

    dat <- gene_overlap_results()

    ggplot(dat, aes(x = reorder(Comparison, Odds_Ratio), y = Odds_Ratio)) +
      geom_col(fill = "#2C3E50") +
      coord_flip() +
      labs(
        x = NULL,
        y = "Fisher odds ratio",
        title = "Gene-level overlap enrichment"
      ) +
      theme_minimal(base_size = 13)
  })

  # ============================================================
  # Downloads
  # ============================================================

  output$download_peak_summary <- downloadHandler(
    filename = function() {
      "peak_overlap_summary.csv"
    },
    content = function(file) {
      write.csv(peak_overlap_results(), file, row.names = FALSE)
    }
  )

  output$download_gene_summary <- downloadHandler(
    filename = function() {
      "gene_overlap_fisher_summary.csv"
    },
    content = function(file) {
      write.csv(gene_overlap_results(), file, row.names = FALSE)
    }
  )

  output$download_gene_lists <- downloadHandler(
    filename = function() {
      "annotated_gene_lists.csv"
    },
    content = function(file) {
      write.csv(gene_data()$gene_table, file, row.names = FALSE)
    }
  )

  output$download_overlap_genes <- downloadHandler(
    filename = function() {
      "overlap_gene_lists.csv"
    },
    content = function(file) {
      write.csv(overlap_gene_results(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)