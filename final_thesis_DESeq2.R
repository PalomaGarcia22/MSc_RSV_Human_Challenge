suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
})
install.packages ("dplyr")


gene_counts <- read_excel("deseq_day7_data.xlsx")


# ==============================
# STEP 1: assign infected and uninfected
# ==============================
infected_subjects <- c(
  "CP001", "CP002", "CP009", "CP016", "CP019", "CP029",
  "CP003", "CP013", "CP021", "CP023", "CP018", "CP038"
)

uninfected_subjects <- c(
  "CP007", "CP008", "CP012", "CP015", "CP014", "CP024", "CP033", "CP025"
)


gene_counts_2 <- gene_counts |>
  filter(count > 0) |>
  mutate(
    condition = ifelse(participant %in% infected_subjects, "infected", "uninfected")
  )

count_df <- gene_counts_2 |>
  select(gene, participant, count) |>
  pivot_wider(
    names_from = participant, 
    values_from = count,
    values_fill = 0)
  
count_matrix <- count_df |>
  column_to_rownames("gene") |>
  as.matrix()

storage.mode(count_matrix) <- "integer"

coldata <- gene_counts_2 |>
  select(participant, condition) |>
  distinct() |>
  column_to_rownames("participant")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

dds <- estimateSizeFactors(dds, type = "ratio")
dds <- DESeq(dds, fitType = "local")
vsdata <- varianceStabilizingTransformation(dds, blind = FALSE)

plotPCA(vsdata, intgroup = "condition")

plotDispEsts(dds)

results <- results(dds, contrast = c("condition", "infected", "uninfected"))

results

sigs <-na.omit(results)
sigs <- sigs[sigs$padj < 0.05,]
sigs

results_df <- as.data.frame(results)
results_df <- results_df |>
  rownames_to_column(var = "gene")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager:: install('EnhancedVolcano')
  
library(EnhancedVolcano)
install.packages("textshaping")

volcano_results <- results_df |>
  EnhancedVolcano(
    x = "log2FoldChange",
    y = "padj",
    xlab = "Log2 fold change (infected vs uninfected)",
    ylab = "-Log10 adjusted p-value",
    pCutoff = 1e-4,
    FCcutoff = 1,
    lab = results_df$gene,
    labSize = 3.0,
    selectLab = c('Large protein', 'small hydroho'),
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    legendPosition = 'right',
    legendLabSize = 8,
    legendIconSize = 2.0
  ) +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

volcano_results



  
  
  
  