rm(list = ls())
library(EValue)
library(deltapif) #remotes::install_github("RodrigoZepeda/deltapif")
library(tidyverse)
library(ggtext)
library(scales)

#Global values
pval <- 0.438
sigma_p_val <- 0.0184^2

df <- tibble(
  year         = c(2025, 2030, 2040, 2050, 2060),
  cases        = c(7.16, 8.53, 11.16, 12.73, 13.85),
  cases_low    = c(6.78, 8.07, 10.55, 11.99, 12.98),
  cases_up     = c(7.55, 8.99, 11.77, 13.46, 14.71),
  cft_increase = c(1.0, 1.15, 1.30, 1.45, 1.60)
)

#Add coverage
df <- df |> mutate(coverage = 100*!!pval*cft_increase)

#Add incidence
df <- df |> mutate(incidence = cases - lag(cases, default = 0))

#Get the variance of the cases from the confidence interval
get_var <- function(up, low){
  ((up - low) / (2 * qnorm(0.975)))^2

}
df <- df |> mutate(variance = get_var(cases_up, cases_low))
df <- df |> mutate(variance_incidence = variance + lag(variance, default = 0))

#Remove year 2025
df <- df |> filter(year > 2025)

#Add empty columns for pif
df <- df |>
  mutate(pif_wu = NA_real_, pif_wu_low = NA_real_,
         pif_wu_up = NA_real_,
         cases_wu = NA_real_, cases_wu_low = NA_real_,
         cases_wu_up = NA_real_) |>
  mutate(pif_wu_ideal = NA_real_, pif_wu_ideal_low = NA_real_,
         pif_wu_ideal_up = NA_real_,
         cases_wu_ideal = NA_real_, cases_wu_ideal_low = NA_real_,
         cases_wu_ideal_up = NA_real_) |>
  mutate(pif_yin = NA_real_, pif_yin_low = NA_real_,
         pif_yin_up = NA_real_,
         cases_yin = NA_real_, cases_yin_low = NA_real_,
         cases_yin_up = NA_real_) |>
  mutate(pif_yin_ideal = NA_real_, pif_yin_ideal_low = NA_real_,
         pif_yin_ideal_up = NA_real_,
         cases_yin_ideal = NA_real_, cases_yin_ideal_low = NA_real_,
         cases_yin_ideal_up = NA_real_)

#1) WU PIF------
UCL_RR_WU <- HR(0.72, rare = FALSE) |> toRR() |> as.numeric()
LCL_RR_WU <- HR(0.67, rare = FALSE) |> toRR() |> as.numeric()
RR_val_WU <- HR(0.69, rare = FALSE) |> toRR() |> as.numeric()
SIGMA_WU2 <- get_var(log(UCL_RR_WU), log(LCL_RR_WU))

for (k in 1:nrow(df)){

  #PIF (guidelines)------
  #Calculate the pif:
  my_pif <- pif(pval, beta = log(RR_val_WU), p_cft = df$coverage[k]/100,
                var_p = sigma_p_val, var_beta = SIGMA_WU2, link = "logit")

  df$pif_wu[k]     <- coef(my_pif)
  df$pif_wu_low[k] <- confint(my_pif)[1]
  df$pif_wu_up[k]  <- confint(my_pif)[2]

  #Calculate the cases
  my_cases <- averted_cases(df$incidence[k], my_pif, variance = df$variance_incidence[k], link = "log")

  df$cases_wu[k]     <- coef(my_cases)
  df$cases_wu_low[k] <- confint(my_cases)[1]
  df$cases_wu_up[k]  <- confint(my_cases)[2]

  #PIF (ideal)------
  #Calculate the pif in the ideal scenario:
  ideal_pif <- pif(pval, beta = log(RR_val_WU), p_cft = 1,
                var_p = sigma_p_val, var_beta = SIGMA_WU2, link = "logit")

  df$pif_wu_ideal[k]     <- coef(ideal_pif)
  df$pif_wu_ideal_low[k] <- confint(ideal_pif)[1]
  df$pif_wu_ideal_up[k]  <- confint(ideal_pif)[2]

  ideal_cases <- averted_cases(df$incidence[k], ideal_pif, variance = df$variance_incidence[k], link = "log")

  df$cases_wu_ideal[k]     <- coef(ideal_cases)
  df$cases_wu_ideal_low[k] <- confint(ideal_cases)[1]
  df$cases_wu_ideal_up[k]  <- confint(ideal_cases)[2]

}

#2) YIN PIF-------
UCL_RR_YIN <- HR(0.76, rare = FALSE) |> toRR() |> as.numeric()
LCL_RR_YIN <- HR(0.66, rare = FALSE) |> toRR() |> as.numeric()
RR_val_YIN <- HR(0.71, rare = FALSE) |> toRR() |> as.numeric()
SIGMA_YIN2 <- get_var(log(UCL_RR_YIN), log(LCL_RR_YIN))

for (k in 1:nrow(df)){

  #PIF (guidelines)------
  #Calculate the pif:
  my_pif <- pif(pval, beta = log(RR_val_YIN), p_cft = df$coverage[k]/100,
                var_p = sigma_p_val, var_beta = SIGMA_YIN2, link = "logit")

  df$pif_yin[k]     <- coef(my_pif)
  df$pif_yin_low[k] <- confint(my_pif)[1]
  df$pif_yin_up[k]  <- confint(my_pif)[2]

  #Calculate the cases
  my_cases <- averted_cases(df$incidence[k], my_pif, variance = df$variance_incidence[k], link = "log")

  df$cases_yin[k]     <- coef(my_cases)
  df$cases_yin_low[k] <- confint(my_cases)[1]
  df$cases_yin_up[k]  <- confint(my_cases)[2]

  #PIF (ideal)------
  #Calculate the pif in the ideal scenario:
  ideal_pif <- pif(pval, beta = log(RR_val_YIN), p_cft = 1,
                   var_p = sigma_p_val, var_beta = SIGMA_YIN2, link = "logit")

  df$pif_yin_ideal[k]     <- coef(ideal_pif)
  df$pif_yin_ideal_low[k] <- confint(ideal_pif)[1]
  df$pif_yin_ideal_up[k]  <- confint(ideal_pif)[2]

  ideal_cases <- averted_cases(df$incidence[k], ideal_pif, variance = df$variance_incidence[k], link = "log")

  df$cases_yin_ideal[k]     <- coef(ideal_cases)
  df$cases_yin_ideal_low[k] <- confint(ideal_cases)[1]
  df$cases_yin_ideal_up[k]  <- confint(ideal_cases)[2]

}

df |>
  mutate(wu_PIF = paste0(percent(pif_wu, suffix = "", accuracy = 0.01),
                             " [", percent(pif_wu_low, suffix = "", accuracy = 0.01), ", ",
                             percent(pif_wu_up, suffix = "", accuracy = 0.01), "]")) |>
  mutate(wu_Prevented = paste0(comma(cases_wu*1e6/1000, accuracy = 1),
                             " [", comma(cases_wu_low*1e6/1000, accuracy = 1), ", ",
                             comma(cases_wu_up*1e6/1000, accuracy = 1), "]")) |>
  mutate(yin_PIF = paste0(percent(pif_yin, suffix = "", accuracy = 0.01),
                         " [", percent(pif_yin_low, suffix = "", accuracy = 0.01), ", ",
                         percent(pif_yin_up, suffix = "", accuracy = 0.01), "]")) |>
  mutate(yin_Prevented = paste0(comma(cases_yin*1e6/1000, accuracy = 1),
                               " [", comma(cases_yin_low*1e6/1000, accuracy = 1), ", ",
                               comma(cases_yin_up*1e6/1000, accuracy = 1), "]")) |>
  select(year, incidence, variance_incidence, wu_PIF, wu_Prevented, yin_PIF, yin_Prevented) |>
  mutate(Scenario = "Target") |>
  bind_rows(
    df |>
      mutate(wu_PIF = paste0(percent(pif_wu_ideal, suffix = "", accuracy = 0.01),
                             " [", percent(pif_wu_ideal_low, suffix = "", accuracy = 0.01), ", ",
                             percent(pif_wu_ideal_up, suffix = "", accuracy = 0.01), "]")) |>
      mutate(wu_Prevented = paste0(comma(cases_wu_ideal*1e6/1000, accuracy = 1),
                                   " [", comma(cases_wu_ideal_low*1e6/1000, accuracy = 1), ", ",
                                   comma(cases_wu_ideal_up*1e6/1000, accuracy = 1), "]")) |>
      mutate(yin_PIF = paste0(percent(pif_yin_ideal, suffix = "", accuracy = 0.01),
                              " [", percent(pif_yin_ideal_low, suffix = "", accuracy = 0.01), ", ",
                              percent(pif_yin_ideal_up, suffix = "", accuracy = 0.01), "]")) |>
      mutate(yin_Prevented = paste0(comma(cases_yin_ideal*1e6/1000, accuracy = 1),
                                    " [", comma(cases_yin_ideal_low*1e6/1000, accuracy = 1), ", ",
                                    comma(cases_yin_ideal_up*1e6/1000, accuracy = 1), "]")) |>
      select(year, incidence, variance_incidence, wu_PIF, wu_Prevented, yin_PIF, yin_Prevented) |>
      mutate(Scenario = "Ideal")
  ) |>
  mutate(Scenario = factor(Scenario, levels = c("Target", "Ideal"), ordered = TRUE)) |>
  arrange(year, Scenario) |>
  select(year, Scenario, incidence, variance_incidence, everything()) |>
  mutate(incidence = paste0(comma(incidence*1e6/1000, accuracy = 1), " (sd: ",
         comma(sqrt(variance_incidence)*1e6/1000, accuracy = 1),")")) |>
  select(-variance_incidence) |>
  mutate(space_1 = "") |>
  mutate(space_2 = "") |>
  mutate(space_3 = "") |>
  select(year, Scenario, space_1, incidence, space_2, wu_PIF, wu_Prevented, space_3, yin_PIF, yin_Prevented) |>
  write_excel_csv("table_1.csv")

#SENSITIVITY ANALYSIS-------
df_sens <- tibble()

hr_vals <- seq(0.05, 0.95, length.out = 50)

RR_val_YIN <- HR(0.71, rare = FALSE) |> toRR() |> as.numeric()
SIGMA_YIN2 <- get_var(log(UCL_RR_YIN), log(LCL_RR_YIN))


#Loop through each of the values
for (val in hr_vals){
  for (k in 1:nrow(df)){

    rr_val <- HR(val, rare = FALSE) |> toRR() |> as.numeric()

    #PIF (guidelines)------
    #Calculate the pif:
    sens_pif <- pif(pval, beta = log(rr_val), p_cft = df$coverage[k]/100,
                    var_p = sigma_p_val, var_beta = max(SIGMA_WU2, SIGMA_YIN2), link = "logit")

    df_temp_pif <- as.data.frame(sens_pif) |>
      mutate(hr = val) |>
      select(-label) |>
      mutate(year = df$year[k]) |>
      mutate(coverage = df$coverage[k]) |>
      mutate(type = "PIF")

    #Calculate the cases
    sens_cases <- averted_cases(df$incidence[k], sens_pif, variance = df$variance_incidence[k], link = "log")

    df_temp_cases <- as.data.frame(sens_cases) |>
      mutate(hr = val) |>
      select(-label) |>
      mutate(year = df$year[k]) |>
      mutate(coverage = df$coverage[k]) |>
      mutate(type = "Cases")

    #PIF (ideal)------
    #Calculate the pif in the ideal scenario:
    sens_ideal_pif <- pif(pval, beta = log(rr_val), p_cft = 1,
                          var_p = sigma_p_val, var_beta = max(SIGMA_WU2, SIGMA_YIN2), link = "logit")

    df_temp_pif_ideal <- as.data.frame(sens_ideal_pif) |>
      mutate(hr = val) |>
      select(-label) |>
      mutate(year = df$year[k]) |>
      mutate(coverage = df$coverage[k]) |>
      mutate(type = "Ideal PIF")

    ideal_cases <- averted_cases(df$incidence[k], sens_ideal_pif, variance = df$variance_incidence[k], link = "log")

    df_temp_ideal_cases <- as.data.frame(ideal_cases) |>
      mutate(hr = val) |>
      select(-label) |>
      mutate(year = df$year[k]) |>
      mutate(coverage = df$coverage[k]) |>
      mutate(type = "Ideal cases")

    df_sens <- df_sens |>
      bind_rows(df_temp_pif) |>
      bind_rows(df_temp_pif_ideal) |>
      bind_rows(df_temp_cases) |>
      bind_rows(df_temp_ideal_cases)
  }
}

df_only_100 <- df_sens |>
  filter(type == "Ideal PIF") |>
  select(-year, -coverage) |>
  distinct()

df_sens_plot <- df_sens |>
  filter(type == "PIF")

df_labels <- df_sens_plot |>
  group_by(coverage) |>
  filter(hr == min(hr)) |>
  ungroup() |>
  mutate(year = paste0(year - 10, "-", year)) |>
  mutate(year = if_else(year == "2020-2030", "2025-2030", year))

df_labels_100 <- df_only_100 |>
  filter(hr == min(hr))

plt_pif_sens <- ggplot(df_sens_plot) +
  geom_ribbon(aes(x = hr, ymin = ci_low, ymax = ci_up, fill = as.character(coverage)), alpha = 0.25) +
  geom_ribbon(aes(x = hr, ymin = ci_low, ymax = ci_up, fill = "Ideal"), alpha = 0.25,
              data = df_only_100) +
  geom_vline(aes(xintercept = 0.69), linetype = "dotted") +
  geom_vline(aes(xintercept = 0.71), linetype = "dotted") +
  geom_richtext(inherit.aes = FALSE, aes(x = 0.69, y = 0.8, label = "Wu _et al_"),
                angle = 90, hjust = 1, vjust = 0, data = tibble(), size = 3) +
  geom_richtext(inherit.aes = FALSE, aes(x = 0.71, y = 0.8, label = "Yin _et al_"),
                angle = 90, hjust = 1, vjust = 1, data = tibble(), size = 3) +
  geom_line(aes(x = hr, y = value, color = as.character(coverage))) +
  geom_line(aes(x = hr, y = value, color = "Ideal"), data = df_only_100) +
  geom_text(aes(x = hr, y = value, color = as.character(coverage), label = year),
            data = df_labels, hjust = 1, size = 3) +
  geom_text(aes(x = hr, y = value, color = "Ideal", label = "Ideal"),
            data = df_labels_100, hjust = 1, size = 3) +
  labs(
    x = "HR",
    y = "PIF",
    title = "Potential Impact Fraction (PIF) for Shingle Vaccination as a function of Hazard Ratio (HR) values",
    caption = paste0("_Vaccine coverage per period:_ ", paste0("**", df_labels$year,"**: ", df$coverage,"%", collapse = "; "), "; **Ideal**: 100%")
  ) +
  scale_y_continuous(labels = scales::percent_format(), n.breaks = 10) +
  scale_x_continuous(limits = c(0, 1), n.breaks = 10) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.caption = element_markdown()
  ) +
  scale_color_manual("Coverage (%):", values = MetBrewer::met.brewer("Hokusai3", n = 5)) +
  scale_fill_manual("Coverage (%):", values = MetBrewer::met.brewer("Hokusai3", n = 5))
ggsave("Figure_1.pdf", plot = plt_pif_sens, width = 8, height = 4)

#Create a figure of the cases and averted cases
df_plot <- df |>
  mutate(decade = paste0(year, "-", year + 10)) |>
  mutate(decade = factor(decade, levels = decade, ordered = T)) |>
  mutate(wu_cumulative = cases_wu) |> #cumsum(cases_wu)) |>
  mutate(yin_cumulative = cases_yin) |> #cumsum(cases_yin)) |>
  mutate(incidence_cumulative = incidence) |> #cumsum(incidence)) |>
  select(decade, year, incidence_cumulative, wu_cumulative, yin_cumulative) |>
  pivot_longer(cols = c(incidence_cumulative, wu_cumulative, yin_cumulative), names_to = "Estimate") |>
  mutate(value = value*1e6) |>
  mutate(estimate = case_when(
    str_detect(Estimate, "wu") ~ "Wu et al",
    str_detect(Estimate, "yin") ~ "Yin et al",
    .default = "Incident cases"
  ))

df_plot <- df_plot |>
  left_join(
    df_plot |>
      filter(Estimate == "incidence_cumulative") |>
      rename(Total = value) |>
      select(-Estimate, -estimate)
  ) |>
  mutate(pct = value / Total)

plt_pif_cases <- ggplot(df_plot) +
  geom_text(aes(x = as.character(year), y = value, color = estimate,
                label = scales::percent(pct)), position = position_dodge(width = 0.9),
            vjust = 0, size = 3) +
  geom_col(aes(x = as.character(year), y = value, fill = estimate), position = position_dodge(width = 0.9)) +
  labs(
    x = "Decade",
    y = "Incident AD cases",
    title = "Averted Alzheimer's Disease (AD) cases attributable to increased\nshingles vaccination in the United States by decade",
    subtitle = "All values represent incident values for each decade.",
    caption = paste0("_Vaccine coverage per period (endpoints):_ ", paste0("**", df$year,"**: ", df$coverage,"%", collapse = "; "))
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.caption = element_markdown()
  ) +
  scale_y_continuous(labels = scales::comma_format(), n.breaks = 10) +
  scale_color_manual("Estimate (%):", values = MetBrewer::met.brewer("Hokusai3", n = 3)) +
  scale_fill_manual("Estimate (%):", values = MetBrewer::met.brewer("Hokusai3", n = 3))
ggsave("Supplementary_Figure_1.pdf", plot = plt_pif_cases, width = 8, height = 4)
