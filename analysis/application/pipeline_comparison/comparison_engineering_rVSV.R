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
  selection_params = list(method = "variance", top_n = 5000),
  model_params     = list(method = "glmnet", inner_folds = 5, metric = "r2")
)

option_choices = list(
  "Genewise" = list(
    method = "engineer",
    col_transform = "z",
    genesets = NULL,
    agg_method = "median"
  ),
  "Median" = list(
    method = "engineer",
    col_transform = "z",
    genesets = genesets,
    agg_method = "median"
  ),
  "Max" = list(
    method = "engineer",
    col_transform = "z",
    genesets = genesets,
    agg_method = "max"
  ),
  "PC1" = list(
    method = "engineer",
    col_transform = "z",
    genesets = genesets,
    agg_method = "pc1"
  ),
  "GSVA" = list(
    method = "engineer",
    col_transform = "z",
    genesets = genesets,
    agg_method = "gsva",
    gsva_min_size = 2
  ),
  "ssGSEA" = list(
    method = "engineer",
    col_transform = "z",
    genesets = genesets,
    agg_method = "ssgsea",
    ssgsea_min_size = 2
  )
)

future::plan(future::multisession, workers = 7)
res = predictomics::compare_pipelines(
  X = X,
  Y = Y,
  option_type = "engineering",
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

p1 = p1 + ggtitle(label = "Pipeline comparison: feature engineering (PREVAC-rVSV)")

p1

ggsave(fs::path(figure_path, "comparison_engineering_rVSV.pdf"),
       p1, width = 8, height = 4.5, dpi = 300)

# rm(list = ls())