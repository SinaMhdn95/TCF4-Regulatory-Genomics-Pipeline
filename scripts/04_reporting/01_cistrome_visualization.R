# Visualize cistrome results using Forrest, mcclay, and NPC enriched datasets
install.packages("pheatmap")

# Load libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(tibble)


# Import datasets
forrest_1k <-  read.csv("Cistrome_results/Forrest_1kpeaks_result.csv")
forrest_10k <- read.csv("Cistrome_results/Forrest_10kpeaks_result.csv")
mcclay_consensus_1k <- read.csv("Cistrome_results/McClay_TCF4_11322_consensus_hg38_1kpeak_result.csv")
mcclay_consensus_10k <- read.csv("Cistrome_results/McClay_TCF4_11322_consensus_hg38_10kpeaks_result.csv")
npc_1k_summit500 <- read.csv("Cistrome_results/NPC_1kpeaks_summit_500bp.bed_result.csv")
npc_10k_summit500 <- read.csv("Cistrome_results/NPC_10kpeaks_summit_500bp.bed_result.csv")


colnames(forrest_1k)
colnames(forrest_10k)

# Add labels to files
# Store data frames in a list
datasets <-  list(
  forrest_1k,
  forrest_10k,
  mcclay_consensus_1k,
  mcclay_consensus_10k,
  npc_1k_summit500,
  npc_10k_summit500
)

# Define a second object which contains the lables
dataset_names <- c(
  "SH-SY5Y-Forrest",
  "SH-SY5Y-Forrest",
  "SH-SY5Y-McClay",
  "SH-SY5Y-McClay",
  "NPC",
  "NPC"
)

peak_labels <- c(
  "1k",
  "10k",
  "1k",
  "10k",
  "1k",
  "10k"
)

# Make a loop to assign the labels
for (i in seq_along(datasets)) {
  datasets[[i]] <- datasets[[i]] %>%
  mutate(
    Dataset = dataset_names[i],
    PeakSet = peak_labels[i]
  )
}

head(datasets[[2]])

# combine all datasets
combined_df <-  bind_rows(datasets)

table(combined_df$Dataset, combined_df$PeakSet)
head(combined_df)
colnames(combined_df)
summary(combined_df$GIGGLE_score)
class(combined_df$GIGGLE_score)

# Clollapse duplicate factors
factor_df <- combined_df %>%
  group_by(Dataset, PeakSet, Factor) %>%
  summarise(
    max_score = max(GIGGLE_score),
    mean_score = mean(GIGGLE_score),
    n_studies = n(),
    .groups = "drop"
  )

# Check the numbers of factors we have
length(unique(factor_df$Factor))

# Final check before plotting
table(factor_df$Dataset)

factor_df %>%
  group_by(Factor) %>%
  summarise(n_datasets = n_distinct(Dataset))

# Explore the overlaps
factor_df %>%
  group_by(Factor) %>%
  summarise(
    n_datasets = n_distinct(Dataset),
    datasets = paste(sort(unique(Dataset)),
                     collapse = ", ")
  ) %>%
  arrange(desc(n_datasets))

# Make a bubble plot
factor_plot <- factor_df %>%
  group_by(Dataset, Factor) %>%
  summarise(
    score = max(max_score),
    .groups = "drop"
  )

factor_df %>%
  group_by(Factor) %>%
  summarise(
    n_datasets = n_distinct(Dataset),
    datasets = paste(sort(unique(Dataset)),
                     collapse = ", ")
  ) %>%
  arrange(desc(n_datasets)) %>%
  print(n = 30)

# Make a bubble plot
head(factor_plot)

ggplot(factor_plot, aes(x = Dataset,
                        y = Factor)) +
  geom_point(aes(color = score)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold")
  )

# Make a matrix for heatmap
matrix_TCF4 <- factor_plot %>%
  pivot_wider(
  names_from = Dataset,
  values_from = score,
  values_fill = 0
)

dim(matrix_TCF4)

# Make a heatmap
# moves the Factor names to the left side as row names
heatmap_matrix <- matrix_TCF4 %>%
  column_to_rownames("Factor") %>%
  as.matrix()

head(heatmap_matrix)

# Create a heatmap
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE
)
