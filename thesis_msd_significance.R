# =========================================================
# RSV MSD mixed-effects analysis and longitudinal plotting
# =========================================================

msd_raw <- read_excel("msd_significance_data.xlsx")

infected_subjects <- c(
  "CP001", "CP002", "CP009", "CP016", "CP019", "CP029",
  "CP003", "CP013", "CP021", "CP023", "CP018", "CP038"
)

uninfected_subjects <- c(
  "CP007", "CP008", "CP012", "CP015", "CP014", "CP024", "CP033", "CP025"
)

msd_raw <- msd_raw |>
  mutate(
    infection = factor(infection, levels = c("uninfected", "infected")),
    participant = as.factor(participant),
    time_num = as.factor(time_num),
    analyte = as.factor(analyte),
    value = as.numeric(value),
    value_log10 = log10(value + 0.01) 
)

## IFN_gamma

df_ifn_g <- msd_raw |>
  filter(analyte == "ifn_g")

model_ifn_g <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_ifn_g
)

## IL-1b

df_il_1b <- msd_raw |>
  filter(analyte == "il_1b")

model_il_1b <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_il_1b
)

## IL-6

df_il_6 <- msd_raw |>
  filter(analyte == "il_6")

model_il_6 <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_il_6
)

## TNF-a

df_tnf_a <- msd_raw |>
  filter(analyte == "tnf_a")

model_tnf_a <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_tnf_a
)


## IL-15

df_il_15 <- msd_raw |>
  filter(analyte == "il_15")

model_il_15 <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_il_15
)


## VEGF

df_vegf <- msd_raw |>
  filter(analyte == "vegf")

model_vegf <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_vegf
)


## IL-1a

df_il_1a <- msd_raw |>
  filter(analyte == "il_1a")

model_il_1a <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_il_1a
)

## MIP-1a

df_mip_1a <- msd_raw |>
  filter(analyte == "mip_1a")

model_mip_1a <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_mip_1a
)

## MIP-1b

df_mip_1b <- msd_raw |>
  filter(analyte == "mip_1b")

model_mip_1b <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_mip_1b
)

## IP-10

df_ip_10 <- msd_raw |>
  filter(analyte == "ip_10")

model_ip_10 <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_ip_10
)

## I-TAC

df_i_tac <- msd_raw |>
  filter(analyte == "i_tac")

model_i_tac <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_i_tac
)



## IL-29

df_il_29_ifnl1 <- msd_raw |>
  filter(analyte == "il_29_ifnl1")

model_il_29_ifnl1 <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_il_29_ifnl1
)


## TSLP

df_tslp <- msd_raw |>
  filter(analyte == "tslp")

model_tslp  <- lmer(
  value_log10 ~ infection * timepoint + (1 | participant),
  data = df_tslp 
)

## combine models

model_list <- list(
  ifn_g = model_ifn_g,
  il_1b = model_il_1b,
  tslp = model_tslp,
  i_tac = model_i_tac,
  il_15 = model_il_15,
  il_1a = model_il_1a,
  iL_19_ifnl1 = model_il_29_ifnl1,
  il_6 = model_il_6,
  ip_10 = model_ip_10,
  mip_1a = model_mip_1a,
  mip_1b = model_mip_1b,
  tnf_a = model_tnf_a,
  vegf = model_vegf
)

library(dplyr)
library(purrr)
library(tibble)

results_tbl <- imap_dfr(model_list, function(model, analyte_name) {
  a <- as.data.frame(anova(model))
  a$term <- rownames(a)
  
  p_col <- names(a)[grepl("Pr", names(a))][1]
  
  tibble(
    analyte = analyte_name,
    term = a$term,
    p_value = a[[p_col]]
  ) %>%
    filter(term %in% c("infection", "timepoint", "infection:timepoint"))
})

results_tbl <- results_tbl %>%
  group_by(term) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

library(emmeans)

## IFN_gamma

emm_ifn_g <- emmeans(model_ifn_g, ~ infection | timepoint)
pairs(emm_ifn_g, adjust = "BH")

posthoc_ifn_g <- summary(pairs(emm_ifn_g, adjust = "BH")) %>%
  as.data.frame()

## IL-1b

emm_il_1b <- emmeans(model_il_1b, ~ infection | timepoint)
pairs(emm_il_1b, adjust = "BH")

posthoc_il_1b <- summary(pairs(emm_il_1b, adjust = "BH")) %>%
  as.data.frame()


## IL-6

emm_il_6 <- emmeans(model_il_6, ~ infection | timepoint)
pairs(emm_il_6, adjust = "BH")

posthoc_il_6 <- summary(pairs(emm_il_6, adjust = "BH")) %>%
  as.data.frame()

## TNF-a

emm_tnf_a <- emmeans(model_tnf_a, ~ infection | timepoint)
pairs(emm_tnf_a, adjust = "BH")

posthoc_tnf_a <- summary(pairs(emm_tnf_a, adjust = "BH")) %>%
  as.data.frame()

## IL-15

emm_il_15 <- emmeans(model_il_15, ~ infection | timepoint)
pairs(emm_il_15, adjust = "BH")

posthoc_il_15 <- summary(pairs(emm_il_15, adjust = "BH")) %>%
  as.data.frame()


## VEGF

emm_vegf <- emmeans(model_vegf, ~ infection | timepoint)
pairs(emm_vegf, adjust = "BH")

posthoc_vegf <- summary(pairs(emm_vegf, adjust = "BH")) %>%
  as.data.frame()

## IL-1a

emm_il_1a <- emmeans(model_il_1a, ~ infection | timepoint)
pairs(emm_il_1a, adjust = "BH")

posthoc_il_1a <- summary(pairs(emm_il_1a, adjust = "BH")) %>%
  as.data.frame()

## MIP-1a

emm_mip_1a <- emmeans(model_mip_1a, ~ infection | timepoint)
pairs(emm_mip_1a, adjust = "BH")

posthoc_mip_1a <- summary(pairs(emm_mip_1a, adjust = "BH")) %>%
  as.data.frame()

## MIP-1b

emm_mip_1b <- emmeans(model_mip_1b, ~ infection | timepoint)
pairs(emm_mip_1b, adjust = "BH")

posthoc_mip_1b <- summary(pairs(emm_mip_1b, adjust = "BH")) %>%
  as.data.frame()

## IP-10

emm_ip_10 <- emmeans(model_ip_10, ~ infection | timepoint)
pairs(emm_ip_10, adjust = "BH")

posthoc_ip_10 <- summary(pairs(emm_ip_10, adjust = "BH")) %>%
  as.data.frame()

## I-TAC

emm_i_tac <- emmeans(model_i_tac, ~ infection | timepoint)
pairs(emm_i_tac, adjust = "BH")

posthoc_i_tac <- summary(pairs(emm_i_tac, adjust = "BH")) %>%
  as.data.frame()


## IL-29

emm_il_29_ifnl1 <- emmeans(model_il_29_ifnl1, ~ infection | timepoint)
pairs(emm_il_29_ifnl1, adjust = "BH")

posthoc_il_29_ifnl1 <- summary(pairs(emm_il_29_ifnl1, adjust = "BH")) %>%
  as.data.frame()


## TSLP

emm_tslp <- emmeans(model_tslp, ~ infection | timepoint)
pairs(emm_tslp, adjust = "BH")

posthoc_tslp <- summary(pairs(emm_tslp, adjust = "BH")) %>%
  as.data.frame()


library(purrr)
library(dplyr)
library(tibble)
library(emmeans)

posthoc_tbl <- imap_dfr(model_list, function(model, analyte_name) {
  emm <- emmeans(model, ~ infection | timepoint)
  
  summary(pairs(emm, adjust = "BH")) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(analyte = analyte_name)
})

sig_posthoc_tbl <- posthoc_tbl |>
  filter(p.value < 0.05)

sig_posthoc_by_analyte <- sig_posthoc_tbl|>
  group_split(analyte) |>
  setNames(unique(sig_posthoc_tbl$analyte))




