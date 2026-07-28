#

library(predictomics)
library(SurrogateRank)
library(tidyverse)
library(clipr)
library(tidyverse)
library(grid)
library(gridExtra)

# Folder to save figures
figure_path <- fs::path("output", "figures", "pipeline_illustration")

# Load data
df_merged_path = fs::path("data", "df_merged_all.rds")
gs_path = fs::path("data", "BTM_processed.rds")

df_merged_all = readRDS(df_merged_path)
gs = readRDS(gs_path)

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

genesets = gs[["genesets"]]
names(genesets) = gs[["geneset.names.descriptions"]]

time.in = Sys.time()
future::plan(future::multisession, workers = 7)
res = predictomics::compare_pipelines(
  X = X,
  Y = Y,
  option_type = "selection",
  option_choices = list(
    "Spearman (p = 100)" = list(method = "spearman", top_n = 100),
    "Spearman (p = 50)" = list(method = "spearman", top_n = 50),
    "Spearman (p = 20)" = list(method = "spearman", top_n = 20)
  ),
  reference_params = list(
    engineering_params = list(
      method = "engineer",
      col_transform = "z",
      genesets = genesets,
      agg_method = "mean"
    ),
    selection_params = list(method = "variance", top_n = 100),
    model_params     = list(method = "glmnet")
  ),
  treatment = treatment,
  treatment_predictor = T,
  verbose = T,
  covariates = covariates,
  cv_type = "kfold",
  folds = 5,
  seed = 12345,
  outside_cv = F
)
future::plan(future::sequential)  # reset after use
time.out = Sys.time()


time.diff = time.out - time.in 
time.diff


p1 = plot(res, metric = "sRMSE")

p1

ggsave(fs::path(figure_path, "Figure6-5.pdf"),
       p1, width = 7, height = 4, dpi = 300)

# rm(list = ls())
