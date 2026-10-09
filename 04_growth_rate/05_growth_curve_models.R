# Step 05: exponential / logistic growth-curve fits of the Day 0 - Day 7 time courses
#          (Supplementary Fig. 2; Supplementary Data 4, 5).
#
# Each of the 90 wells is fitted with
#   exponential: ln(area) = intercept + r * t                   (lm on log area)
#   logistic   : area = K / (1 + (K / N0 - 1) * exp(-r * t))    (nls, port algorithm)
# AIC is computed from the residual sum of squares on the log-area scale
# (AIC = n ln(RSS / n) + 2k; k = 2 exponential, 3 logistic).
# Three strategies are compared: exponential-fixed, logistic-fixed and AIC-selected
# (lower AIC per well). For each strategy six wells per species x light are taken
# (smallest starting frond area, then higher R2) and the growth rate r is analysed with
# the same Gamma GLMMs as the RGR (step 06); the interaction P values are reported in
# Supplementary Data 4.
#
# Input : output/frond_area_selected_by_image.csv (03)
# Output: output/SupplementaryData04_growth_curve_model_strategy_summary.csv
#         output/SupplementaryData05_growth_curve_model_well_parameters.csv  (90 wells)
#         output/SupplementaryFig2.pdf   pages 1-3 = panels a (exponential), b (logistic), c (AIC-selected)
# Run from this directory:  Rscript 05_growth_curve_models.R
# Author: Natsu Katayama
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(car)
  library(readr)
  library(tibble)
})

fig3_samples <- c("Sp", "Lp", "Lgp8L", "M10E", "Wh", "Wa")
light_levels <- c("PFD100", "PFD1000")
light_cols <- c("PFD100" = "#00BFC4", "PFD1000" = "#F8766D")
antho_samples <- c("Sp", "Lp", "Lgp8L")
strategies <- c("exponential_fixed", "logistic_fixed", "aic_selected")

theme_fig <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  )

area <- read_csv("output/frond_area_selected_by_image.csv", show_col_types = FALSE) |>
  mutate(
    sample = factor(sample, levels = fig3_samples),
    light = factor(light, levels = light_levels),
    antho = factor(if_else(as.character(sample) %in% antho_samples, "antho", "non-antho"),
                   levels = c("antho", "non-antho"))
  )

aic_from_rss <- function(rss, n, k) {
  if (!is.finite(rss) || rss <= 0 || n <= k) {
    return(NA_real_)
  }
  n * log(rss / n) + 2 * k
}

fit_one_growth_curve <- function(d) {
  d <- d |> filter(frond_area > 0) |> mutate(elapsed_days = elapsed_time / 24)
  if (nrow(d) < 4 || length(unique(d$elapsed_days)) < 4) {
    return(tibble())
  }

  log_area <- log(d$frond_area)
  exp_fit <- lm(log_area ~ d$elapsed_days)
  exp_pred_log <- as.numeric(predict(exp_fit))
  exp_rss_log <- sum((log_area - exp_pred_log)^2)
  exp_tss_log <- sum((log_area - mean(log_area))^2)
  exp_intercept <- unname(coef(exp_fit)[1])
  exp_rate <- unname(coef(exp_fit)[2])
  exp_aic_log <- aic_from_rss(exp_rss_log, nrow(d), 2)

  logistic_fit <- tryCatch(
    suppressWarnings(
      nls(
        frond_area ~ K / (1 + ((K / N0) - 1) * exp(-r * elapsed_days)),
        data = d,
        start = list(
          K = max(d$frond_area, na.rm = TRUE) * 1.25,
          N0 = max(min(d$frond_area, na.rm = TRUE), 1e-3),
          r = max(exp_rate, 1e-3)
        ),
        algorithm = "port",
        lower = c(K = max(d$frond_area, na.rm = TRUE) * 1.001, N0 = 1e-6, r = 1e-6),
        upper = c(K = max(d$frond_area, na.rm = TRUE) * 100, N0 = max(d$frond_area, na.rm = TRUE), r = 10),
        control = nls.control(maxiter = 200, warnOnly = TRUE)
      )
    ),
    error = function(e) NULL
  )

  logistic_ok <- !is.null(logistic_fit)
  if (logistic_ok) {
    logistic_coef <- coef(logistic_fit)
    logistic_pred <- pmax(as.numeric(predict(logistic_fit, newdata = d)), 1e-9)
    logistic_rss_log <- sum((log_area - log(logistic_pred))^2)
    logistic_aic_log <- aic_from_rss(logistic_rss_log, nrow(d), 3)
    logistic_r2_log <- if_else(exp_tss_log > 0, 1 - logistic_rss_log / exp_tss_log, NA_real_)
  } else {
    logistic_coef <- c(K = NA_real_, N0 = NA_real_, r = NA_real_)
    logistic_aic_log <- NA_real_
    logistic_r2_log <- NA_real_
  }

  # batch and wellID (the grouping keys) are added back by group_modify()
  tibble(
    sample = d$sample[1],
    light = d$light[1],
    sampleID2 = d$sampleID2[1],
    antho = d$antho[1],
    n_timepoints = nrow(d),
    elapsed_time_min = min(d$elapsed_time),
    elapsed_time_max = max(d$elapsed_time),
    exponential_growth_rate_per_day = exp_rate,
    exponential_intercept_log_area = exp_intercept,
    exponential_aic_log = exp_aic_log,
    exponential_r_squared_log = if_else(exp_tss_log > 0, 1 - exp_rss_log / exp_tss_log, NA_real_),
    logistic_growth_rate_per_day = unname(logistic_coef["r"]),
    logistic_N0 = unname(logistic_coef["N0"]),
    logistic_K = unname(logistic_coef["K"]),
    logistic_aic_log = logistic_aic_log,
    logistic_r_squared_log = logistic_r2_log,
    delta_aic_logistic_minus_exponential = logistic_aic_log - exp_aic_log
  )
}

fit_all <- area |>
  group_by(batch, wellID) |>
  group_modify(\(d, key, ...) fit_one_growth_curve(d)) |>
  ungroup() |>
  mutate(
    sample = factor(as.character(sample), levels = fig3_samples),
    light = factor(as.character(light), levels = light_levels),
    antho = factor(antho, levels = c("antho", "non-antho"))
  )

start_area_by_well <- area |>
  filter(frond_area > 0) |>
  group_by(batch, wellID, sampleID2, sample, light) |>
  arrange(elapsed_time, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(batch, wellID, sampleID2, sample, light,
            start_elapsed_time = elapsed_time, start_frond_area = frond_area)

fit_all <- fit_all |>
  left_join(start_area_by_well, by = c("batch", "wellID", "sampleID2", "sample", "light"))

write.csv(fit_all, "output/SupplementaryData05_growth_curve_model_well_parameters.csv", row.names = FALSE)

observed <- area |>
  filter(frond_area > 0) |>
  mutate(elapsed_days = elapsed_time / 24) |>
  select(batch, wellID, sampleID2, elapsed_days, elapsed_time, frond_area)

fit_gamma_glmm <- function(formula, data) {
  glmer(formula, data = data, family = Gamma(link = "log"), nAGQ = 0,
        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

prepare_strategy <- function(strategy) {
  fit_all |>
    mutate(
      strategy = strategy,
      selected_model = case_when(
        strategy == "exponential_fixed" ~ "exponential",
        strategy == "logistic_fixed" ~ "logistic",
        !is.na(logistic_aic_log) & logistic_aic_log < exponential_aic_log ~ "logistic",
        TRUE ~ "exponential"
      ),
      growth_rate_per_day = if_else(selected_model == "logistic",
                                    logistic_growth_rate_per_day, exponential_growth_rate_per_day),
      r_squared_log = if_else(selected_model == "logistic",
                              logistic_r_squared_log, exponential_r_squared_log)
    ) |>
    filter(is.finite(growth_rate_per_day), growth_rate_per_day > 0, is.finite(r_squared_log))
}

run_strategy <- function(strategy) {
  growth_rate_main <- prepare_strategy(strategy) |>
    group_by(sample, light) |>
    arrange(start_frond_area, desc(r_squared_log), abs(delta_aic_logistic_minus_exponential), .by_group = TRUE) |>
    slice_head(n = 6) |>
    ungroup()

  growth_fit_points <- growth_rate_main |>
    select(batch, wellID, sample, light, sampleID2, antho,
           selected_model, exponential_intercept_log_area, exponential_growth_rate_per_day,
           logistic_N0, logistic_K, logistic_growth_rate_per_day) |>
    left_join(observed, by = c("batch", "wellID", "sampleID2")) |>
    mutate(
      fitted_frond_area = if_else(
        selected_model == "logistic",
        logistic_K / (1 + ((logistic_K / logistic_N0) - 1) * exp(-logistic_growth_rate_per_day * elapsed_days)),
        exp(exponential_intercept_log_area + exponential_growth_rate_per_day * elapsed_days)
      )
    )

  model_sample <- fit_gamma_glmm(growth_rate_per_day ~ sample * light + (1 | batch), growth_rate_main)
  anova_sample <- as.data.frame(Anova(model_sample, type = 2)) |> rownames_to_column("Effect")
  model_antho <- fit_gamma_glmm(growth_rate_per_day ~ antho * light + (1 | sample) + (1 | batch), growth_rate_main)
  anova_antho <- as.data.frame(Anova(model_antho, type = 2)) |> rownames_to_column("Effect")

  plot <- ggplot(growth_fit_points, aes(elapsed_days, frond_area, colour = light, group = sampleID2)) +
    geom_point(size = 0.9, alpha = 0.55) +
    geom_line(aes(y = fitted_frond_area, linetype = selected_model), linewidth = 0.35, alpha = 0.60) +
    scale_colour_manual(values = light_cols, name = "Light condition") +
    facet_wrap(~sample, scales = "free_y", ncol = 3) +
    labs(x = "Elapsed time (days)", y = expression(Frond~area~(mm^2))) +
    theme_fig

  summary <- tibble(
    strategy = strategy,
    sample_light_complete = all((growth_rate_main |> count(sample, light))$n == 6),
    min_r_squared_log = min(growth_rate_main$r_squared_log, na.rm = TRUE),
    median_r_squared_log = median(growth_rate_main$r_squared_log, na.rm = TRUE),
    max_growth_rate_per_day = max(growth_rate_main$growth_rate_per_day, na.rm = TRUE),
    n_logistic_selected = sum(growth_rate_main$selected_model == "logistic"),
    n_exponential_selected = sum(growth_rate_main$selected_model == "exponential"),
    antho_light_interaction_p = anova_antho$`Pr(>Chisq)`[anova_antho$Effect == "antho:light"],
    sample_light_interaction_p = anova_sample$`Pr(>Chisq)`[anova_sample$Effect == "sample:light"]
  )
  list(summary = summary, plot = plot)
}

results <- lapply(strategies, run_strategy)

summary_df <- bind_rows(lapply(results, `[[`, "summary"))
write.csv(summary_df, "output/SupplementaryData04_growth_curve_model_strategy_summary.csv", row.names = FALSE)
print(as.data.frame(summary_df))

pdf("output/SupplementaryFig2.pdf", width = 8.8, height = 6.5)
for (res in results) print(res$plot)
invisible(dev.off())
