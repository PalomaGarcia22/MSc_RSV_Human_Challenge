suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
})

install.packages("ggpubr")
library(ggpubr)
install.packages("rstatix")
library(rstatix)

gene_prop <- read_excel("fractional_analysis.xlsx")


## Large Protein ##

large_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

large_prop <- large_prop |>
  filter(gene == "large protein") |>
  filter(infection_status == "infecte")
         
test_large <- large_prop[order(large_prop$subject, large_prop$timepoint), ]
test_df_large <- pairwise_wilcox_test(
  data = test_large,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_large <- subset(
  test_df_large, 
    (group1 == "1" & group2 == "7") | 
    (group1 == "7" & group2 == "14"))

desired_comparisons_large$y.position <- c(0.75, 0.80)

ggplot (test_large,
        aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "red") + 
  stat_pvalue_manual(desired_comparisons_large, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
#  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of L Transcript") + 
  theme_classic()
  
## Nucleoprotein ##

nucl_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

nucl_prop <- nucl_prop |>
  filter(gene == "nucleoprotein") |>
  filter(infection_status == "infecte")

test_nucl <- nucl_prop[order(nucl_prop$subject, nucl_prop$timepoint), ]
test_df_nucl <- pairwise_wilcox_test(
  data = test_nucl,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

#desired_comparisons_nucl <- subset(
#  test_df_nucl, 
#  (group1 == "1" & group2 == "7") |
#    (group1 == "7" & group2 == "28"))

#desired_comparisons_nucl$y.position <- c(0.25)

ggplot (test_nucl,
        aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "gray") + 
#  stat_pvalue_manual(desired_comparisons_nucl, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
#  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of N Transcript") + 
  theme_classic()

## G protein ##

g_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

g_prop <- g_prop |>
  filter(gene == "attachment protein") |>
  filter(infection_status == "infecte")

test_g <- g_prop[order(g_prop$subject, g_prop$timepoint), ]
test_df_g <- pairwise_wilcox_test(
  data = test_g,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

ggplot (test_g,
        aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "gray") + 
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of G Transcript") + 
  theme_classic()


## F protein ##

f_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

f_prop <- f_prop |>
  filter(gene == "fusion protein") |>
  filter(infection_status == "infecte")

test_f <- f_prop[order(f_prop$subject, f_prop$timepoint), ]
test_df_f <- pairwise_wilcox_test(
  data = test_f,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

ggplot (test_f,
        aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "gray") + 
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of F Transcript") + 
  theme_classic()

## M protein ##

m_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

m_prop <- m_prop |>
  filter(gene == "matrix protein") |>
  filter(infection_status == "infecte")

test_m <- m_prop[order(m_prop$subject, m_prop$timepoint), ]
test_df_m <- pairwise_wilcox_test(
  data = test_m,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_m <- subset(
  test_df_m, 
  (group1 == "3" & group2 == "7") | 
    (group1 == "7" & group2 == "14") | 
    (group1 == "7" & group2 == "28"))

desired_comparisons_m$y.position <- c(0.17, 0.20, 0.17)

ggplot(test_m,
        aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "gray") + 
  stat_pvalue_manual(desired_comparisons_m, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of M Transcript") + 
  theme_classic()


## M2-1 protein ##

ma_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

ma_prop <- ma_prop |>
  filter(gene == "matrix protein 2-1") |>
  filter(infection_status == "infecte")

test_ma <- ma_prop[order(ma_prop$subject, ma_prop$timepoint), ]
test_df_ma <- pairwise_wilcox_test(
  data = test_ma,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_ma <- subset(
  test_df_ma, 
    (group1 == "2" & group2 == "3"))

desired_comparisons_ma$y.position <- c(0.08)

ggplot(test_ma,
       aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "gray") + 
  stat_pvalue_manual(desired_comparisons_ma, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of M2-1 Transcript") + 
  theme_classic()


## M2-2 protein ##

mb_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

mb_prop <- mb_prop |>
  filter(gene == "matrix protein 2-2") |>
  filter(infection_status == "infecte")

test_mb <- mb_prop[order(mb_prop$subject, mb_prop$timepoint), ]
test_df_mb <- pairwise_wilcox_test(
  data = test_mb,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_mb <- subset(
  test_df_mb, 
  (group1 == "3" & group2 == "7") | 
    (group1 == "7" & group2 == "14") |
    (group1 == "7" & group2 == "28"))

desired_comparisons_mb$y.position <- c(0.025, 0.02, 0.025)

ggplot(test_mb,
       aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "gray") + 
  stat_pvalue_manual(desired_comparisons_mb, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of M2-2 Transcript") + 
  theme_classic()

## NS1 protein ##

nsa_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

nsa_prop <- nsa_prop |>
  filter(gene == "nonstructural protein 1") |>
  filter(infection_status == "infecte")

test_nsa <- nsa_prop[order(nsa_prop$subject, nsa_prop$timepoint), ]
test_df_nsa <- pairwise_wilcox_test(
  data = test_nsa,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_nsa <- subset(
  test_df_nsa, 
  (group1 == "2" & group2 == "7") | 
    (group1 == "7" & group2 == "14") |
    (group1 == "7" & group2 == "28"))

desired_comparisons_nsa$y.position <- c(0.17, 0.2, 0.17)

ggplot(test_nsa,
       aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "lightgreen") + 
  stat_pvalue_manual(desired_comparisons_nsa, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of NS1 Transcript") + 
  theme_classic()

## NS2 protein ##

nsb_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

nsb_prop <- nsb_prop |>
  filter(gene == "nonstructural protein 2") |>
  filter(infection_status == "infecte")

test_nsb <- nsb_prop[order(nsb_prop$subject, nsb_prop$timepoint), ]
test_df_nsb <- pairwise_wilcox_test(
  data = test_nsb,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_nsb <- subset(
  test_df_nsb, 
  (group1 == "2" & group2 == "28") | 
    (group1 == "7" & group2 == "28"))

desired_comparisons_nsb$y.position <- c(0.17, 0.2)

ggplot(test_nsb,
       aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "lightgreen") + 
  stat_pvalue_manual(desired_comparisons_nsb, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of NS2 Transcript") + 
  theme_classic()

## P protein ##

p_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

p_prop <- p_prop |>
  filter(gene == "phosphoprotein") |>
  filter(infection_status == "infecte")

test_p <- p_prop[order(p_prop$subject, p_prop$timepoint), ]
test_df_p <- pairwise_wilcox_test(
  data = test_p,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_p <- subset(
  test_df_p, 
  (group1 == "2" & group2 == "7") | 
    (group1 == "2" & group2 == "14") |
    (group1 == "7" & group2 == "14"))

desired_comparisons_p$y.position <- c(0.2, 0.25, 0.3)

ggplot(test_p,
       aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "lightgreen") + 
  stat_pvalue_manual(desired_comparisons_p, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of P Transcript") + 
  theme_classic()


## SH protein ##

sh_prop = subset(gene_prop, select = c(gene, subject, timepoint, infection_status, gene_fraction))

sh_prop <- sh_prop |>
  filter(gene == "small hyrophobic protein") |>
  filter(infection_status == "infecte")

test_sh <- sh_prop[order(sh_prop$subject, sh_prop$timepoint), ]
test_df_sh <- pairwise_wilcox_test(
  data = test_sh,
  formula = gene_fraction ~ timepoint,
  paired = TRUE, 
  id = "subject") |>
  add_xy_position(x = "timepoint")

desired_comparisons_sh <- subset(
  test_df_sh, 
    (group1 == "7" & group2 == "28"))

desired_comparisons_sh$y.position <- c(0.10)

ggplot(test_sh,
       aes(x = factor(timepoint), y = gene_fraction)
) + 
  geom_boxplot(fill = "lightgreen") + 
  stat_pvalue_manual(desired_comparisons_sh, label = "p", tip.length = 0.01) +
  geom_point(size = 2, alpha = 0.5) + 
  #  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  labs(
    title = "Gene Fraction Over Time (Paired)",
    x = "Day (p.i.)",
    y = "Proportion of SH Transcript") + 
  theme_classic()
