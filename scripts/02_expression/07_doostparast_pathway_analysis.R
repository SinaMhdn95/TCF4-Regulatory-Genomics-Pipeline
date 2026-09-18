################################################################################
# Pathway Analysis for genes and overlaps obtained from Doostparast paper
# Pathway Analysis Using clusterProfiler
# GO Biological Process Enrichment
################################################################################
####################################################################
# Import All Overlap Gene Lists and Run clusterProfiler
####################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(readr)

## Create output folder
dir.create(
  "Decon_Doost/Pathway_Analysis",
  showWarnings = FALSE
)

## Find all overlap gene-list files
gene_list_files <- list.files(
  path = "Decon_Doost",
  pattern = "_overlap_genes\\.csv$",
  full.names = TRUE
)

gene_list_files
length(gene_list_files)

####################################################################
# Read all gene lists
####################################################################
gene_lists <- list()

for (file in gene_list_files) {

  list_name <- gsub(
    "_overlap_genes\\.csv$",
    "",
    basename(file)
  )

  temp_df <- read.csv(file)

  gene_lists[[list_name]] <- unique(
    na.omit(temp_df$Gene)
  )
}

names(gene_lists)

sapply(gene_lists, length)

####################################################################
# Run GO Biological Process Enrichment for All Gene Lists
####################################################################
go_results <- list()

for (list_name in names(gene_lists)) {

  gene_symbols <- gene_lists[[list_name]]

  ## Skip very small gene lists
  if (length(gene_symbols) < 10) {
    message("Skipping ", list_name, ": fewer than 10 genes")
    next
  }

  ## Convert SYMBOL to ENTREZID
  entrez_ids <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  ## Skip if too few genes successfully mapped
  if (nrow(entrez_ids) < 10) {
    message("Skipping ", list_name, ": fewer than 10 mapped genes")
    next
  }

  ## GO Biological Process enrichment
  go_result <- enrichGO(
    gene = unique(entrez_ids$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  go_results[[list_name]] <- go_result

  ## Save GO results
  write.csv(
    as.data.frame(go_result),
    paste0(
      "Decon_Doost/Pathway_Analysis/",
      list_name,
      "_GO_BP_results.csv"
    ),
    row.names = FALSE
  )
}

####################################################################
# Make GO Dotplots
####################################################################
for (list_name in names(go_results)) {

  go_df <- as.data.frame(go_results[[list_name]])

  if (nrow(go_df) > 0) {

    p <- dotplot(
      go_results[[list_name]],
      showCategory = 20
    ) +
      ggtitle(
        paste0(list_name, " - GO Biological Process")
      )

    ggsave(
      filename = paste0(
        "Decon_Doost/Pathway_Analysis/",
        list_name,
        "_GO_BP_dotplot.pdf"
      ),
      plot = p,
      width = 9,
      height = 7
    )
  }
}

####################################################################
# Create Summary Table of GO Results
####################################################################
go_summary <- data.frame(
  Gene_List = names(go_results),
  Number_of_GO_Terms = sapply(
    go_results,
    function(x) nrow(as.data.frame(x))
  )
)

go_summary

write.csv(
  go_summary,
  "Decon_Doost/Pathway_Analysis/GO_BP_summary.csv",
  row.names = FALSE
)

####################################################################
# Clean GO BP Comparison Plots
# Top 10 / Top 15 terms per selected comparison
####################################################################
####################################################################
# Function to make clean comparison dotplots
####################################################################
go_files <- list.files(
  path = "Decon_Doost/Pathway_Analysis",
  pattern = "_GO_BP_results\\.csv$",
  full.names = TRUE
)

go_all <- lapply(go_files, function(file) {

  df <- read.csv(file)

  df$Gene_List <- gsub(
    "_GO_BP_results\\.csv$",
    "",
    basename(file)
  )

  df
})

go_all <- bind_rows(go_all)

## Convert GeneRatio from "x/y" to numeric
go_all <- go_all %>%
  mutate(
    GeneRatio_numeric = sapply(
      GeneRatio,
      function(x) {
        nums <- as.numeric(strsplit(x, "/")[[1]])
        nums[1] / nums[2]
      }
    ),
    minus_log10_padj = -log10(p.adjust)
  )

make_go_comparison_plot <- function(
    go_data,
    selected_lists,
    plot_title,
    output_name,
    top_n = 10
) {

  plot_df <- go_data %>%
    filter(Gene_List %in% selected_lists) %>%
    group_by(Gene_List) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  plot_df <- plot_df %>%
    mutate(
      Description_wrapped = str_wrap(Description, width = 45),
      Gene_List = factor(Gene_List, levels = selected_lists)
    )

  p <- ggplot(
    plot_df,
    aes(
      x = GeneRatio_numeric,
      y = fct_reorder(Description_wrapped, GeneRatio_numeric)
    )
  ) +
    geom_point(
      aes(
        size = Count,
        color = minus_log10_padj
      ),
      alpha = 0.9
    ) +
    facet_wrap(
      ~ Gene_List,
      scales = "free_y",
      nrow = 1
    ) +
    scale_color_gradient(
      low = "steelblue",
      high = "firebrick",
      name = "-log10(adj. P)"
    ) +
    scale_size_continuous(
      name = "Gene count",
      range = c(2.5, 7)
    ) +
    labs(
      title = plot_title,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      text = element_text(face = "bold"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      axis.text.y = element_text(
        size = 12,
        face = "bold"
      ),
      axis.text.x = element_text(size = 9),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  ggsave(
    filename = paste0(
      "Decon_Doost/Pathway_Analysis/",
      output_name,
      "_top",
      top_n,
      ".pdf"
    ),
    plot = p,
    width = 14,
    height = 7
  )

  ggsave(
    filename = paste0(
      "Decon_Doost/Pathway_Analysis/",
      output_name,
      "_top",
      top_n,
      ".png"
    ),
    plot = p,
    width = 14,
    height = 7,
    dpi = 300
  )

  return(p)
}

####################################################################
# Comparison 1
# Day 3 NPC all DEGs: NPC peaks vs McClay peaks vs Forrest peaks
####################################################################
comp1_lists <- c(
  "D3_NPC_peakNPC_all",
  "D3_NPC_peakMcClay_all",
  "D3_NPC_peakForrest_all"
)

p_comp1_top10 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp1_lists,
  plot_title = "Day 3 NPC all DEGs: NPC vs McClay vs Forrest TCF4 peaks",
  output_name = "Comp1_D3_NPC_all_three_peaksets",
  top_n = 10
)

p_comp1_top15 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp1_lists,
  plot_title = "Day 3 NPC all DEGs: NPC vs McClay vs Forrest TCF4 peaks",
  output_name = "Comp1_D3_NPC_all_three_peaksets",
  top_n = 15
)

####################################################################
# Comparison 2
# Day 3 NPC directional:
# NPC up peaks vs McClay down peaks vs Forrest down peaks
####################################################################
comp2_lists <- c(
  "D3_NPC_peakNPC_up",
  "D3_NPC_peakMcClay_down",
  "D3_NPC_peakForrest_down"
)

p_comp2_top10 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp2_lists,
  plot_title = "Day 3 NPC directional TCF4 target enrichment",
  output_name = "Comp2_D3_NPC_directional",
  top_n = 10
)

p_comp2_top15 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp2_lists,
  plot_title = "Day 3 NPC directional TCF4 target enrichment",
  output_name = "Comp2_D3_NPC_directional",
  top_n = 15
)


####################################################################
# Comparison 3
# Day 14 glut all DEGs:
# NPC peaks vs McClay peaks vs Forrest peaks
####################################################################
comp3_lists <- c(
  "D14_GLUT_peakNPC_all",
  "D14_GLUT_peakMcClay_all",
  "D14_GLUT_peakForrest_all"
)

p_comp3_top10 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp3_lists,
  plot_title = "Day 14 glut all DEGs: NPC vs McClay vs Forrest TCF4 peaks",
  output_name = "Comp3_D14_GLUT_all_three_peaksets",
  top_n = 10
)

p_comp3_top15 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp3_lists,
  plot_title = "Day 14 glut all DEGs: NPC vs McClay vs Forrest TCF4 peaks",
  output_name = "Comp3_D14_GLUT_all_three_peaksets",
  top_n = 15
)

####################################################################
# Comparison 4
# Day 14 glut NPC peaks: Up vs Down
####################################################################
comp4_lists <- c(
  "D14_GLUT_peakNPC_up",
  "D14_GLUT_peakNPC_down"
)

p_comp4_top10 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp4_lists,
  plot_title = "Day 14 glut NPC peak targets: Up vs Down",
  output_name = "Comp4_D14_GLUT_NPC_up_vs_down",
  top_n = 10
)

p_comp4_top15 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp4_lists,
  plot_title = "Day 14 glut NPC peak targets: Up vs Down",
  output_name = "Comp4_D14_GLUT_NPC_up_vs_down",
  top_n = 15
)

####################################################################
# Comparison 5
# Day 14 glut mixed direction:
# NPC up vs McClay up vs Forrest down
####################################################################
comp5_lists <- c(
  "D14_GLUT_peakNPC_up",
  "D14_GLUT_peakMcClay_up",
  "D14_GLUT_peakForrest_down"
)

p_comp5_top10 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp5_lists,
  plot_title = "Day 14 glut directional comparison: NPC up, McClay up, Forrest down",
  output_name = "Comp5_D14_GLUT_directional_mixed",
  top_n = 10
)

p_comp5_top15 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp5_lists,
  plot_title = "Day 14 glut directional comparison: NPC up, McClay up, Forrest down",
  output_name = "Comp5_D14_GLUT_directional_mixed",
  top_n = 15
)

####################################################################
# Comparison 6
# Cross-stage comparison:
# Day 3 NPC all NPC peaks vs Day 14 glut all NPC peaks
####################################################################
comp6_lists <- c(
  "D3_NPC_peakNPC_all",
  "D14_GLUT_peakNPC_all"
)

p_comp6_top10 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp6_lists,
  plot_title = "Cross-stage comparison: Day 3 NPC vs Day 14 glut using NPC peaks",
  output_name = "Comp6_D3_vs_D14_NPC_peakset",
  top_n = 10
)

p_comp6_top15 <- make_go_comparison_plot(
  go_data = go_all,
  selected_lists = comp6_lists,
  plot_title = "Cross-stage comparison: Day 3 NPC vs Day 14 glut using NPC peaks",
  output_name = "Comp6_D3_vs_D14_NPC_peakset",
  top_n = 15
)

################################################################################
# Pathway analysis for out TCF4 gene lists

# Import gene lists
mcclay_hg38_gene <- read.csv("mcclay_genes_hg38.csv")
length(rownames(mcclay_hg38_gene))
colnames(mcclay_hg38_gene)
head(mcclay_hg38_gene)

npc_hg38_gene <- read.csv("npc_genes_hg38.csv")
length(rownames(npc_hg38_gene))

forrest_hg38_gene <- read.csv("forrest_genes_hg38.csv")
length(rownames(forrest_hg38_gene))


####################################################################
# Extract gene symbols from column x
####################################################################
mcclay_gene_symbols <- unique(
  na.omit(mcclay_hg38_gene$x)
)

npc_gene_symbols <- unique(
  na.omit(npc_hg38_gene$x)
)

forrest_gene_symbols <- unique(
  na.omit(forrest_hg38_gene$x)
)

length(mcclay_gene_symbols)
length(npc_gene_symbols)
length(forrest_gene_symbols)

####################################################################
# Make named gene list
####################################################################
tcf4_gene_lists <- list(
  NPC = npc_gene_symbols,
  McClay = mcclay_gene_symbols,
  Forrest = forrest_gene_symbols
)

####################################################################
# Convert SYMBOL to ENTREZID
####################################################################
tcf4_gene_lists_entrez <- lapply(
  tcf4_gene_lists,
  function(gene_symbols) {

    mapped <- bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )

    unique(mapped$ENTREZID)
  }
)

sapply(tcf4_gene_lists_entrez, length)

####################################################################
# compareCluster GO Biological Process
####################################################################
tcf4_go_compare <- compareCluster(
  geneCluster = tcf4_gene_lists_entrez,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

####################################################################
# Comparison dotplot
####################################################################
p_tcf4_compare <- dotplot(
  tcf4_go_compare,
  showCategory = 10,
  font.size = 13
) +
  ggtitle(
    "GO Biological Process Enrichment of TCF4 Peak-Associated Genes"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

p_tcf4_compare
