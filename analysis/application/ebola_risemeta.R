library(SurrogateRank)
library(tidyverse)

# Folder to save figures
figures_folder <- fs::path("output", "figures", "prediction")

# Load data
df_merged_path = fs::path("data", "df_merged_all.rds")
gs_path = fs::path("data", "BTM_processed.rds")

df_merged_all = readRDS(df_merged_path)
gs = readRDS(gs_path)

tp_group = c("P+7D", "P+0D")
arm_group = c("ebovac2-Ad26MVA",    
              "prevac-Ad26MVA",            
              "prevac-rVSV" ,
              "hamburg-rVSV")

## Filter data
df_merged_all_filtered = df_merged_all %>%
  filter(study_vaccine %in% arm_group, time %in% tp_group) %>%
  mutate(
    treatment = ifelse(time == tp_group[1], 1, 0),
    response = ifelse(time == tp_group[1], ifelse(study_vaccine == "ebovac2-Ad26MVA", ab_p_365, ab_p_180), ab_p_0),
    study = study_vaccine 
  ) %>%
  filter(!is.na(response)) %>%
  group_by(participant_id) %>%
  filter(length(participant_id) == 2) %>%
  ungroup() %>%
  arrange(participant_id) %>%
  dplyr::select(treatment, study, response, a1cf:zzz3) %>%
  dplyr::select(
    treatment,
    study,
    response,
    where(~ !any(is.na(.)))
  )

studyone = df_merged_all_filtered %>% 
  filter(treatment == 1) %>% 
  pull(study)

studyzero = df_merged_all_filtered %>% 
  filter(treatment == 0) %>% 
  pull(study)

yone = df_merged_all_filtered %>% 
  filter(treatment == 1) %>% 
  pull(response)

yzero = df_merged_all_filtered %>% 
  filter(treatment == 0) %>% 
  pull(response)

sone = df_merged_all_filtered %>% 
  filter(treatment == 1) %>% 
  dplyr::select(a1cf:zzz3)

szero = df_merged_all_filtered %>% 
  filter(treatment == 0) %>% 
  dplyr::select(a1cf:zzz3)

aggregate_to_genesets <- function(
    mat,
    genesets,
    geneset_names,
    agg_fun = mean,
    na.rm = TRUE
) {
  
  mat <- as.matrix(mat)
  
  res <- lapply(seq_along(genesets), function(i) {
    genes <- genesets[[i]]
    genes_present <- intersect(genes, colnames(mat))
    
    if (length(genes_present) == 0) {
      return(rep(NA_real_, nrow(mat)))
    }
    
    subset_mat <- mat[, genes_present, drop = FALSE]
    
    if (ncol(subset_mat) == 1) {
      return(as.numeric(subset_mat[, 1]))
    }
    
    apply(subset_mat, 1, function(x) agg_fun(x, na.rm = na.rm))
  })
  
  res <- do.call(cbind, res)
  colnames(res) <- geneset_names
  rownames(res) <- rownames(mat)
  
  return(res)
}

genesets <- gs[["genesets"]]
geneset_names <- gs[["geneset.names.descriptions"]]

sone <- aggregate_to_genesets(sone, genesets, geneset_names, agg_fun = mean)
szero <- aggregate_to_genesets(szero, genesets, geneset_names, agg_fun = mean)


screen.res = rise.screen.meta(yone = yone, 
                              yzero = yzero,
                              sone = sone, 
                              szero = szero,
                              studyone = studyone,
                              studyzero = studyzero,
                              epsilon.meta.mode = "mean.power",
                              power.want.s.study = 0.8,
                              paired.all = TRUE,
                              n.cores = 5, 
                              p.correction = "BH", 
                              return.study.similarity.plot = F)

p1 = screen.res[["gamma.s.plot"]][["screen.plot"]]

p1

# Save figure
TLS_figures_folder <- fs::path("output", "figures", "TLS")

ggsave(
  filename = "RISEMeta_Ebola.pdf",
  path = TLS_figures_folder,
  plot = p1,
  width = 35,
  height = 20,
  units = "cm"
)

