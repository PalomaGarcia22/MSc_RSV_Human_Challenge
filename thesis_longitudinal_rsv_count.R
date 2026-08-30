install.packages("ggh4x")

library(tidyverse)
library(ggh4x)

rsv_data <- read_excel("by_gene_rsv_count.xlsx")


rsv_data |>
  ggplot(
    aes(
      x = factor(timepoint, level = c('BL', 'D1', 'D2', 'D3', 'D7', 'D10', 'D14', 'D28')),
      y = log2_count,
      fill = infection_status)
  ) + 
  geom_boxplot(position = "dodge") + 
  facet_wrap2(~gene, axes = "all") +
  xlab("") +
  ylab("RSV Total Count (log2 + 1)") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 7)) +
  scale_fill_manual(breaks = rsv_data$infection_status,
                    values = c("red", "lightblue")) +
 labs(fill = "Infection Status")

library(lme4)
library(lmerTest)

#df_attach <- rsv_data |>
#  filter(gene == "attachment protein")

#model_attach <- lmer(
#  log10_count ~ timepoint * infection_status + (1 | subject),
#  data = df_attach)
#model <- rsv_data ->
#  group_by



results <- rsv_data |>
  group_by(gene) |>
  do(anova_table = anova(lmer(log10_count ~ timepoint * infection_status + (1 | subject), data = .)))

anova(results, type = 3, ddf = "Satterthwaite")
