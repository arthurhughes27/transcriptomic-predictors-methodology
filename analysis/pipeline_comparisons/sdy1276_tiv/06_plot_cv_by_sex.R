# analysis/pipeline_comparisons/sdy1276_tiv/06_plot_cv_by_sex.R
#
# Standalone figure: SDY1276 (TIV)'s best pipeline (see
# 05_find_best_model.R), cross-validated predicted vs. observed response,
# with points coloured by participant sex.
#
# Uses predictomics::plot.predictomics()'s highlight/highlight_type/
# highlight_label arguments directly on the saved best-model fit, rather
# than re-fitting anything - highlight is purely a plotting overlay, not a
# new predictor, so the winning pipeline's engineering/selection/model
# choices and its cross-validated predictions are unchanged from
# 05_find_best_model.R's own plot(best_fit) panel.
#
# sex is taken from the SAME `single` dataset/row order the saved
# best_model_fit_sdy1276_tiv.rds was fit on (see 01_prepare_data.R,
# 05_find_best_model.R), so it lines up with best_fit$observed one-to-one -
# see plot.predictomics()'s own docs for why highlight must match that
# length/order exactly.
#
# Run analysis/pipeline_comparisons/sdy1276_tiv/05_find_best_model.R first,
# so that output/results/best_model_fit_sdy1276_tiv.rds exists.

library(fs)
library(predictomics)
library(ggplot2)

results_dir <- fs::path("output", "results")

analysis_data <- readRDS(fs::path(results_dir, "sdy1276_tiv_analysis_data.rds"))
best_fit <- readRDS(fs::path(results_dir, "best_model_fit_sdy1276_tiv.rds"))

p_by_sex <- plot(
  best_fit,
  highlight       = analysis_data$single$covariates$sex,
  highlight_type  = "categorical",
  highlight_label = "Sex"
) +
  ggtitle("SDY1276 (TIV): cross-validated prediction by sex")

print(p_by_sex)

figure_path <- fs::path("output", "figures", "pipeline_comparisons", "sdy1276_tiv")
fs::dir_create(figure_path)

ggsave(
  fs::path(figure_path, "best_model_fit_by_sex_sdy1276_tiv.pdf"),
  p_by_sex, width = 6, height = 4.5, dpi = 300
)

rm(list = ls())
