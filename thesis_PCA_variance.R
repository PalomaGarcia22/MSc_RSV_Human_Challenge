suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
})


install.packages ("dplyr")




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

vst_mat <- read_csv("vst_matrix.csv")

pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

summary(pca)

pca_scores <- as.data.frame(pca$x)

pca_scores <- as.data.frame(pca$x) |>
  tibble::rownames_to_column("sample") |>
  left_join(meta, by = "sample")
write.csv(pca_scores, file.path(tables_dir, "pca_scores.csv"), row.names = FALSE)


pc_loading_df <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("gene")
gene_load <- write.csv(pc_loading_df, "pca_gene_loadings.csv")

ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2), color = "blue") +
  geom_segment(data = pc_loading_df, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "red",
               position = position_dodge(width = 1)) +
  geom_text_repel(data = pc_loading_df, aes(x = PC1, y = PC2, label = gene), 
                  color = "red", vjust = -2, size = 3, 
                  box.padding = 0.5,
                  point.padding = 0.3) +
  theme_classic()

##To just focus on PC1 and PC2
pca_scores_plot <- pca_scores |>
  filter(!is.na(infection_status), !is.na(timepoint), !is.na(PC1), !is.na(PC2))

install.packages("MASS")
install.packages("factoextra")
library(MASS)
library(factoextra)

pca_scores_x <- as.data.frame(pca$x)

biplot(pca$x,type = "form", labels = "variables")

inf_uninf <- ggplot(pca_scores_plot, aes(PC1, PC2, color = infection_status)) +
  stat_ellipse(aes(group = infection_status, color = infection_status), type = "norm", linewidth = 0.7, linetype = "solid", level = 0.8, show.legend = FALSE) +
  geom_point(size = 3) +
  scale_color_manual(values = status_colors, drop = FALSE) +
  theme_bw(base_size = 11) +
  labs(
    title = "PCA of RSV expression: PC1 vs PC2",
    x = "PC1 explains 41.99% variance",
    y = "PC2 expains 21.32% variance",
    color = "Infection Status"
  )

inf_uninf

# Separate PCA plot colored by timepoint
pca_timepoint <- pca_scores_plot |>
  filter(infection_status == "infected") |>
  ggplot(
    aes(PC1, PC2, color = timepoint)) +
  geom_point(aes(colour = factor(timepoint)),
             size = 3) +
  scale_colour_manual(values = c("black", "lightgreen", "pink", "orange", "red", "brown", "purple", "lightblue")) +
  theme_bw(base_size = 11) +
  labs(
    title = "PCA of RSV expression colored by timepoint: PC1 vs PC2",
    x = "PC1 explains 41.99% variance",
    y = "PC2 expains 21.32% variance",
    color = "Timepoint",
    shape = "Timepoint"
  )

pca_timepoint

# Separate PCA plot colored by total RSV counts
pca_rsv <- pca_scores_plot |>
  filter(infection_status == "infected") |>
  ggplot(
    aes(PC1, PC2, color = total_rsv_counts)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_viridis_c(option = "Turbo")+
  theme_bw(base_size = 11) +
  labs(
    title = "PCA of RSV expression colored by timepoint: PC1 vs PC2",
    x = "PC1 explains 41.99% variance",
    y = "PC2 expains 21.32% variance",
    color = "Total RSV counts"
  )
pca_rsv


pcaData <- plotPCA(pca, intgroup=c("infection_status", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

pca_result <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE) 

fviz_pca_biplot(pca_result, 
                label = "var",
                habillage = meta$infection_status,
                col.var = "black",
                scale_color_manual(values = c("red", "mediumblue")))



