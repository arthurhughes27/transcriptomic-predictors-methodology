library(predictomics)
library(dplyr)

# Load data
df_merged_path = fs::path("data", "df_merged_all.rds")
gs_path = fs::path("data", "BTM_processed.rds")

df_merged_all = readRDS(df_merged_path)
gs = readRDS(gs_path)

# Data preprocessing

arm_group = "SDY1276"
tp_group = c("P+0D", "P+1D")

cov_names = c("age", "sex", "ab_p_0")

## Filter data
df_merged_all_filtered = df_merged_all %>%
  filter(study_accession %in% arm_group, time %in% tp_group) %>%
  filter(sex == "Female") %>%
  mutate(treatment = ifelse(time == tp_group[2], 1, 0), 
         response = ifelse(time == tp_group[2], ab_p_28, ab_p_0)) %>%
  filter(!is.na(response)) %>%
  group_by(participant_id) %>%
  filter(length(participant_id) == 2) %>%
  ungroup() %>%
  arrange(participant_id) %>%
  dplyr::select(participant_id, treatment, time, any_of(cov_names), response, a1bg:zzz3) %>%
  dplyr::select(participant_id, treatment, time, any_of(cov_names), response, where( ~ !any(is.na(.))))

X = df_merged_all_filtered %>% 
  dplyr::select(-participant_id, -response, -treatment, -time, -any_of(cov_names)) %>% 
  as.matrix()

Y = df_merged_all_filtered %>%
  pull(response)

treatment = df_merged_all_filtered %>%
  pull(treatment)

covariates = df_merged_all_filtered %>% 
  dplyr::select(any_of(cov_names))

individuals = df_merged_all_filtered %>%
  pull(participant_id)

timepoints = df_merged_all_filtered %>%
  pull(time) %>% 
  recode("P+0D" = 0, 
         "P+1D" = 1)

genesets = gs[["genesets"]]
names(genesets) = gs[["geneset.names.descriptions"]]

future::plan(future::multisession, workers = 9)
res = predictomics::predict_cv(
  X = X,
  Y = Y,
  treatment = treatment,
  timepoint = timepoints, 
  individual_id = individuals,
  treatment_predictor = F,
  verbose = T,
  covariates = covariates,
  cv_type = "kfold",
  folds = 10,
  seed = 12345,
  engineering_params = list(method = "engineer", gene_level_fc = T, col_transform = "z", genesets = genesets, agg_method = "median"),
  selection_params = list(method = "spearman", top_n = 10),
  model_params = list(method = "glmnet", inner_folds = 5, metric = "r2"),
  outside_cv = F
)
future::plan(future::sequential)  # reset after use

plot(res)
plot_selection_stability(res, top_n = 20, type = "explicit")




library(GSVA)

## GSVA expects genes x samples
expr <- t(X)

## Gene sets
gene_sets <- gs[["genesets"]]

## Keep only genes present in expression matrix
gene_sets <- lapply(
  gene_sets,
  function(g) intersect(g, rownames(expr))
)

## Remove gene sets with too few genes
gene_sets <- gene_sets[lengths(gene_sets) >= 5]


## Create GSVA parameter object
par <- gsvaParam(
  exprData = expr,
  geneSets = gene_sets,
  minSize = 5,
  maxSize = Inf,
  kcdf = "Gaussian"
)

## Run GSVA
gsva_scores <- gsva(par)


