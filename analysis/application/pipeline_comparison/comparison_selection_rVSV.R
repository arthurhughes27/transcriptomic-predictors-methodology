#

library(predictomics)
library(SurrogateRank)
library(tidyverse)
library(clipr)
library(tidyverse)
library(grid)
library(gridExtra)

# Folder to save figures
figure_path <- fs::path("output", "figures", "pipeline_comparison")

# Load data
df_merged_path = fs::path("data", "df_merged_all.rds")
gs_path = fs::path("data", "BTM_processed.rds")

df_merged_all = readRDS(df_merged_path)
gs = readRDS(gs_path)

genesets = gs[["genesets"]]
names(genesets) = gs[["geneset.names.descriptions"]]


# Pipeline on PREVAC data
tp_group = c("P+7D")

gene_names = df_merged_all %>% 
  dplyr::select(a1cf:zzz3) %>% 
  colnames()

## Define the treatment groups of interest
arm_group = c("prevac-rVSV", "prevac-placebo")

cov_names = c("sex", "age", "ab_p_0")

## Filter data, define response
df_merged_all_filtered = df_merged_all %>%
  filter(study_vaccine %in% arm_group, 
         time %in% tp_group) %>%
  mutate(treatment = ifelse(study_vaccine == arm_group[1], 1, 0),
         response = ab_p_180) %>%
  filter(!is.na(response), !is.na(ab_p_0)) %>%
  dplyr::select(participant_id, time, response, treatment, any_of(cov_names), any_of(gene_names)) %>% 
  mutate(response = log(response)) %>%
  dplyr::select(participant_id, time,
                response, treatment, any_of(cov_names),
                where(~ !any(is.na(.)))
  )

X = df_merged_all_filtered %>% 
  dplyr::select(any_of(gene_names)) %>% 
  as.matrix()

Y = df_merged_all_filtered %>%
  pull(response)

treatment = df_merged_all_filtered %>%
  pull(treatment)

covariates = df_merged_all_filtered %>%
  dplyr::select(any_of(cov_names))


# Define reference pipeline
reference_params = list(
  engineering_params = list(
    method = "engineer",
    col_transform = "z",
    genesets = genesets,
    agg_method = "mean"
  ),
  selection_params = NULL,
  model_params     = list(method = "glmnet", inner_folds = 5, metric = "r2")
)

option_choices = list(
  "Variance (p = 100)" = list(method = "variance", top_n = 100),
  "Spearman (p = 100)" = list(method = "spearman", top_n = 100),
  "Spearman (p = 50)" = list(method = "spearman", top_n = 50),
  "Spearman (threshold = 0.5)" = list(method = "spearman", threshold = 0.5),
  "Pearson (p = 100)" = list(method = "pearson", top_n = 100),
  "Pearson (p = 50)" = list(method = "pearson", top_n = 50),
  "Pearson (threshold = 0.5)" = list(method = "pearson", threshold = 0.5),
  "Relative gain (threshold = 0)" = list(method = "relative_gain", threshold = 0, relative_gain_inner_folds = 5, relative_gain_metric = "rmse"),
  "Relative gain (threshold = 0.1)" = list(method = "relative_gain", threshold = 0.1, relative_gain_inner_folds = 5, relative_gain_metric = "rmse"),
  # "RISE (p = 50)" = list(method = "rise", rise_paired = FALSE, rise_power_want_s = 0.8, rise_p_correction = "BH", top_n = 50),
  # "RISE (p = 50) " = list(method = "rise", rise_paired = TRUE, rise_power_want_s = 0.8, rise_p_correction = "BH", top_n = 50),
  "Dearseq DGSA" = list(method = "dearseq", threshold = 0.05, dearseq_mode = "classic", dearseq_level = "geneset", dearseq_covariates = covariates[,-3], genesets = genesets, dearseq_which_test = "asymptotic", dearseq_preprocessed = TRUE, dearseq_n_perm = 5000),
  "Dearseq DGSA " = list(method = "dearseq", threshold = 0.05, dearseq_mode = "paired", dearseq_level = "geneset", dearseq_covariates = covariates[,-3], genesets = genesets, dearseq_which_test = "asymptotic", dearseq_preprocessed = TRUE, dearseq_n_perm = 5000)
)

future::plan(future::multisession, workers = 7)
res = predictomics::compare_pipelines(
  X = X,
  Y = Y,
  option_type = "selection",
  option_choices = option_choices,
  reference_params = reference_params,
  treatment = treatment,
  treatment_predictor = F,
  verbose = T,
  covariates = covariates,
  cv_type = "kfold",
  folds = 5,
  seed = 12345,
  outside_cv = F
)
future::plan(future::sequential)


p1 = plot(res, metric = "R2")

p1 = p1 + ggtitle(label = "Pipeline comparison: feature selection (PREVAC-rVSV)")

p1

ggsave(fs::path(figure_path, "comparison_selection_rVSV.pdf"),
       p1, width = 9, height = 4.5, dpi = 300)

# rm(list = ls())
