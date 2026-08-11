#

library(predictomics)
library(SurrogateRank)
library(tidyverse)
library(clipr)
library(tidyverse)
library(grid)
library(gridExtra)

# Folder to save figures
figure_path <- fs::path("output", "figures", "reference_models")

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
res = predictomics::predict_cv(
  X = X,
  Y = Y,
  treatment = treatment,
  treatment_predictor = T,
  verbose = T,
  covariates = covariates,
  cv_type = "loo",
  folds = 5,
  seed = 12345,
  engineering_params = list(method = "engineer", col_transform = "z", genesets = genesets, agg_method = "mean"),
  # selection_params = list(method = "variance", top_n = 50),
  model_params = list(method = "glmnet", inner_folds = 5, metric = "r2"),
  outside_cv = F
)
future::plan(future::sequential)  # reset after use
time.out = Sys.time()


time.diff = time.out - time.in 
time.diff

p1 = plot(res)

p2 = plot_selection_stability(res, top_n = 20, type = "embedded", plot_type = "frequency")

p1
p2

p1_labeled <- p1 +
  labs(title = "A) Cross-validated predictions ") +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face = "plain"),
    plot.title.position = "panel",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

p2_labeled <- p2 +
  labs(title = "B) Feature selection") +
  theme(
    plot.title = element_text(size = 18, hjust = 0, face = "plain"),
    plot.title.position = "panel",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# Thin vertical divider between panels
divider <- ggplot() +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1),
               colour = "white", linewidth = 0.6) +
  theme_void()

g <- egg::ggarrange(
  p1_labeled, divider, p2_labeled,
  ncol   = 3,
  widths = c(1, 0.01, 1),
  draw   = FALSE
)

g_title <- arrangeGrob(
  g,
  top = textGrob(
    "CV prediction of rVSV vaccine antibody response (reference model)",
    gp = gpar(fontsize = 22, fontface = "bold")
  )
)

grid.newpage()
grid.draw(g_title)

ggsave(fs::path(figure_path, "Figure6-5.pdf"),
       g_title, width = 15, height = 4.8, dpi = 300)

rm(list = ls())
