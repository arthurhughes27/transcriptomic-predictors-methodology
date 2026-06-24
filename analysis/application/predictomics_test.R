# Testing the predictomics package

library(predictomics)
library(SurrogateRank)
library(tidyverse)
library(clipr)
library(tidyverse)
library(caret)

# Folder to save figures
figures_folder <- fs::path("output", "figures", "prediction")

# Load data
df_merged_path = fs::path("data", "df_merged_all.rds")
gs_path = fs::path("data", "BTM_processed.rds")

df_merged_all = readRDS(df_merged_path)
gs = readRDS(gs_path)

## Screen on hamburg

## Define timepoint groups (paired analysis)
tp_group = c("P+7D", "P+0D")

## Define treatment group
arm_group = c("hamburg-rVSV")

gene_names = df_merged_all %>% 
  dplyr::select(a1cf:zzz3) %>% 
  colnames()

## Filter data
df_merged_all_filtered = df_merged_all %>%
  filter(study_vaccine %in% arm_group, time %in% tp_group) %>%
  mutate(
    treatment = ifelse(time == tp_group[1], 1, 0),
    response = ifelse(time == tp_group[1], ab_p_180, ab_p_0)
  ) %>%
  filter(!is.na(response)) %>%
  group_by(participant_id) %>%
  filter(length(participant_id) == 2) %>%
  ungroup() %>%
  arrange(participant_id) %>%
  dplyr::select(treatment, response, a1cf:zzz3)

## response in treated
yone = df_merged_all_filtered %>%
  filter(treatment == 1) %>%
  pull(response)

## response in control
yzero = df_merged_all_filtered %>%
  filter(treatment == 0) %>%
  pull(response)

## surrogates in treated
sone = df_merged_all_filtered %>%
  filter(treatment == 1) %>%
  dplyr::select(a1cf:zzz3)

## surrogates in control
szero = df_merged_all_filtered %>%
  filter(treatment == 0) %>%
  dplyr::select(a1cf:zzz3)

## Apply RISE screen function
rise_screen_result_rVSV_hamburg = rise.screen(
  yone,
  yzero,
  sone,
  szero,
  alpha = 0.05,
  power.want.s = 0.8,
  p.correction = "BH",
  n.cores = 6,
  alternative = "two.sided",
  paired = T,
  weight.mode = "diff.epsilon",
  normalise.weights = T
)

## How many significant markers?
paste0(length(rise_screen_result_rVSV_hamburg[["significant.markers"]]), " significant markers found.")

## Store significant genes and weights
RISE_rVSV_hamburg_genes = rise_screen_result_rVSV_hamburg[["significant.markers"]]

## evaluate on prevac
tp_group = c("P+7D")

## Define the treatment groups of interest
arm_group = c("prevac-rVSV", "prevac-placebo")

cov_names = c("sex", "age", "ab_p_0")

## Filter data, define response
df_merged_all_filtered = df_merged_all %>%
  filter(study_vaccine %in% arm_group, 
         time %in% tp_group) %>%
  mutate(treatment = ifelse(study_vaccine == arm_group[1], 1, 0),
         response = ab_p_28) %>%
  filter(!is.na(response), !is.na(ab_p_0)) %>%
  dplyr::select(response, treatment, any_of(cov_names), any_of(gene_names)) %>% 
  mutate(response = log(response))

X = df_merged_all_filtered %>% 
  dplyr::select(-response, -treatment, -any_of(cov_names)) %>% 
  as.matrix()

Y = df_merged_all_filtered %>%
  pull(response)

treatment = df_merged_all_filtered %>%
  pull(treatment)

covariates = df_merged_all_filtered %>% 
  dplyr::select(any_of(cov_names))

genesets = gs[["genesets"]]
names(genesets) = gs[["geneset.names.descriptions"]]

res = predictomics::predict_cv(
  X = X,
  Y = Y,
  treatment = treatment,
  treatment_predictor = TRUE,
  verbose = T,
  covariates = covariates,
  cv_type = "kfold",
  folds = 10,
  seed = 12345,
  # engineering_params = list(method = "engineer", col_transform = "none", genesets = genesets, agg_method = "pc1"),
  selection_params = list(method = "variance", top_n = 2500),
  model_params = list(method = "glmnet", inner_folds = 5, metric = "r2"),
  outside_cv = T
)

plot(res)

