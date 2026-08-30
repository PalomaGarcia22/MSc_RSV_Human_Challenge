suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(variancePartition)
  library(BiocParallel)
  library(ggrepel)
})

# ==============================
# STEP 0: USER CONFIGURATION
# ==============================
running_local <- TRUE

hpc_base_dir <- "/rds/general/user/pdg25/ephemeral"
local_base_dir <- "~/Desktop/rsv_project"
base_dir <- if (running_local) path.expand(local_base_dir) else hpc_base_dir

counts_file <- file.path(base_dir, "counts_rsv_star", "RSV_gene_counts_matrix.tsv")
outdir <- file.path(base_dir, "counts_rsv_star", "02_variance_partition")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(outdir, "plots")
tables_dir <- file.path(outdir, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)

# Number of extreme samples to label at each PCA edge.
n_extreme_labels_per_axis <- 3

# Manual shapes so all 8 timepoints display in the legend.
timepoint_shape_values <- c(
  "BL" = 16,
  "D1" = 17,
  "D2" = 15,
  "D3" = 3,
  "D7" = 7,
  "D10" = 8,
  "D14" = 1,
  "D28" = 2
)

# Manual colors: uninfected = blue, infected = red.
status_colors <- c(
  "uninfected" = "#1F78B4",
  "infected" = "#E31A1C"
)

# ==============================
# STEP 1: SAMPLE DEFINITIONS
# ==============================
infected_subjects <- c(
  "CP001", "CP002", "CP009", "CP016", "CP019", "CP029",
  "CP003", "CP013", "CP021", "CP023", "CP018", "CP038"
)

uninfected_subjects <- c(
  "CP007", "CP008", "CP012", "CP015", "CP014", "CP024", "CP033", "CP025"
)

# ==============================
# STEP 2: GENE MAP
# ==============================
gene_map <- tribble(
  ~Geneid,       ~gene_name,
  "AIY60636.1", "nonstructural protein 1",
  "AIY60634.1", "nonstructural protein 2",
  "AIY60635.1", "nucleoprotein",
  "AIY60637.1", "phosphoprotein",
  "AIY60638.1", "matrix protein",
  "AIY60639.1", "small hydrophobic protein",
  "AIY60644.1", "attachment protein",
  "AIY60640.1", "fusion protein",
  "AIY60641.1", "matrix protein 2-1",
  "AIY60642.1", "matrix protein 2-2",
  "AIY60643.1", "large protein"
)

extract_subject <- function(sample) stringr::str_extract(sample, "^CP\\d+")
extract_timepoint <- function(sample) {
  stringr::str_extract(sample, "-(BL|D\\d+)-") |>
    stringr::str_remove_all("-")
}
safe_factor_timepoint <- function(tp) {
  factor(tp, levels = c("BL", "D1", "D2", "D3", "D7", "D10", "D14", "D28"))
}

time_to_numeric <- function(tp) {
  dplyr::case_when(
    tp == "BL" ~ 0,
    tp == "D1" ~ 1,
    tp == "D2" ~ 2,
    tp == "D3" ~ 3,
    tp == "D7" ~ 7,
    tp == "D10" ~ 10,
    tp == "D14" ~ 14,
    tp == "D28" ~ 28,
    TRUE ~ NA_real_
  )
}

make_extreme_labels <- function(df, x_col, y_col, n_each = 3) {
  bind_rows(
    df |> slice_min(order_by = .data[[x_col]], n = n_each, with_ties = FALSE),
    df |> slice_max(order_by = .data[[x_col]], n = n_each, with_ties = FALSE),
    df |> slice_min(order_by = .data[[y_col]], n = n_each, with_ties = FALSE),
    df |> slice_max(order_by = .data[[y_col]], n = n_each, with_ties = FALSE)
  ) |>
    distinct(sample, .keep_all = TRUE)
}

# ==============================
# STEP 3: READ COUNTS + METADATA
# ==============================
counts_raw <- read.delim(
  counts_file, sep = "\t", check.names = FALSE,
  comment.char = "", stringsAsFactors = FALSE, quote = ""
)

sample_cols <- names(counts_raw)[grepl("^CP\\d+", names(counts_raw))]
if (length(sample_cols) == 0) stop("No sample columns detected in count matrix")

count_mat <- as.matrix(counts_raw[, sample_cols, drop = FALSE])
storage.mode(count_mat) <- "numeric"
if (nrow(count_mat) == nrow(gene_map)) {
  rownames(count_mat) <- gene_map$gene_name
} else {
  rownames(count_mat) <- paste0("gene_", seq_len(nrow(count_mat)))
}

meta <- tibble(sample = colnames(count_mat)) |>
  mutate(
    subject = vapply(sample, extract_subject, character(1)),
    timepoint = vapply(sample, extract_timepoint, character(1)),
    infection_status = case_when(
      subject %in% infected_subjects ~ "infected",
      subject %in% uninfected_subjects ~ "uninfected",
      TRUE ~ NA_character_
    ),
    timepoint = safe_factor_timepoint(timepoint),
    infection_status = factor(infection_status, levels = c("uninfected", "infected"))
  ) |>
  filter(!is.na(infection_status)) |>
  as.data.frame()

count_mat <- count_mat[, meta$sample, drop = FALSE]
rownames(meta) <- meta$sample

meta <- meta |>
  mutate(
    timepoint_numeric = time_to_numeric(as.character(timepoint)),
    total_rsv_counts = colSums(count_mat)[sample],
    sample_label = paste0(subject, "_", as.character(timepoint))
  )
write.csv(meta, file.path(tables_dir, "sample_metadata_used.csv"), row.names = FALSE)

# ==============================
# STEP 4: VST FOR DIMENSION REDUCTION
# ==============================
keep_genes <- rowSums(count_mat) > 0
count_mat_nonzero <- count_mat[keep_genes, , drop = FALSE]

if (nrow(count_mat_nonzero) < 2 || ncol(count_mat_nonzero) < 2) {
  stop("Not enough non-zero genes or samples for PCA/variance partition")
}

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat_nonzero),
  colData = meta,
  design = ~ infection_status + timepoint
)

dds <- estimateSizeFactors(dds, type = "poscounts")

vst_obj <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_matrix <- assay(vst_obj)

write.csv(vst_matrix, file.path(tables_dir, "good_vst_matrix.csv"))

d <- vst_matrix
gene <- rownames(d)
rownames(d) <- NULL
data <- cbind(gene,d)

vst_matrix <- as.data.frame(data)

correct_vst_matrix <- vst_matrix |>
  pivot_longer(
    cols = !gene,
    names_to = "subject",
    values_to = "total_count"
  ) |>
  separate(subject, into = c("subject", "timepoint", "sample_type"), sep = "-", remove = FALSE) |>
  pivot_wider(names_from = gene, values_from = "total_count")

write.csv(correct_vst_matrix, file.path(tables_dir, "correct_vst_matrix.csv"), row.names = FALSE)

writexl::write_xlsx(
  correct_vst_matrix,
  path = file.path(tables_dir, "correct_vst_matrix.xlsx")
)

