suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
})
install.packages ("dplyr")


overalldata <- read_tsv("RSV_gene_counts_matrix.tsv")

clean_data <- overalldata |>
  pivot_longer(cols = starts_with("CP"), names_to = "participant_ID", values_to = "count") |>
  pivot_wider(names_from = gene, values_from = count)

write.csv(clean_data, "clean_total_data.csv")

sum_data <- read_excel("RSV_counts_histogram.xlsx")
sum_data$day <- factor(sum_data$day, levels = c("0", "1", "2", "3", "7", "10", "14", "28"))



## Infected Distribution ##
sum_data |>
  filter(infection_status == "infected") |>
  ggplot(
    aes(x = day, y = sum_log)
    ) +
  geom_col(
    fill = "red",
    color = "darkred"
    ) +
  labs(
    title = "RNA Distribution Over Time: Infected",
    x = "Day",
    y = "Total RSV Counts (log2)"
  ) + 
  coord_cartesian(ylim = c(0, 200)) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

## Uninfected Distribution ##
sum_data |>
  filter(infection_status == "uninfected") |>
  ggplot(
    aes(x = day, y = sum_log)
  ) +
  geom_col(
    fill = "darkblue",
    color = "blue"
  ) +
  labs(
    title = "RNA Distribution Over Time: Uninfected",
    x = "Day",
    y = "Total RSV Counts (log2)"
  ) + 
  coord_cartesian(ylim = c(0, 200)) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

## Uninfected Distribution ##
sum_data |>
  filter(infection_status == "infected") |>
  ggplot(
    aes(x = as.factor(day), fill = as.factor(sum_log))
  ) +
  geom_bar(
    position = "stack",
    color = "darkblue",
    linewidth = 0.2
  ) +
  scale_fill_viridis_d(name = "RNA Count") + 
  labs(
    title = "Sample Frequency Over Time by RNA Count Level",
    x = "Day",
    y = "Frequency"
  ) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12))
