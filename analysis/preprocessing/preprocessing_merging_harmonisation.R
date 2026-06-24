# Script to merge data from prevac, ebovac2, and hamburg to harmonise naming conventions, data types etc.

# Libraries
library(tidyverse)
library(stringr)

# Load data
processed_data_path <- fs::path("data")
folder_ebovac2 = fs::path("data", "ebovac2")
folder_hamburg = fs::path("data", "hamburg")
folder_prevac = fs::path("data", "prevac")

path_ebovac2 = fs::path(folder_ebovac2, "ebovac2_merged_norm.rds")
path_hamburg = fs::path(folder_hamburg, "hamburg_merged_norm.rds")
path_prevac = fs::path(folder_prevac, "prevac_merged.RData")

df_clinical_all = readRDS(file = fs::path(processed_data_path, "df_clinical_all.rds"))
ebovac2_merged = readRDS(path_ebovac2)
hamburg_merged = readRDS(path_hamburg)
prevac_merged = readRDS(path_prevac)

# Select the intersection of genes
ebovac2_genes = ebovac2_merged %>% 
  dplyr::select(a1bg:zzz3) %>% 
  colnames()

hamburg_genes = hamburg_merged %>% 
  dplyr::select(a1bg:zzz3) %>% 
  colnames()

prevac_genes = prevac_merged %>% 
  dplyr::select(a1bg:zzz3) %>% 
  colnames()

intersection_genes = intersect(intersect(ebovac2_genes, hamburg_genes), prevac_genes)

length(intersection_genes)

ebovac2_merged_filtered = ebovac2_merged %>% 
  mutate(participant_id = as.character(pid),
         time = recode(time, "P+0" = "P+0D")) %>% 
  dplyr::select(participant_id, time, all_of(intersection_genes))

hamburg_merged_filtered = hamburg_merged %>% 
  mutate(participant_id = as.character(pid),
         time = recode(
           time,
           "d1" = "P+1D",
           "d3" = "P+3D",
           "d7" = "P+7D",
           "d0" = "P+0D"
         )) %>% 
  dplyr::select(participant_id, time, all_of(intersection_genes))

prevac_merged_filtered = prevac_merged %>% 
  mutate(participant_id = as.character(pid),
         time = recode(
           timepoint,
           "J0"    = "P+0D",
           "J03h"  = "P+3H",
           "J7"    = "P+7D",
           "J56"   = "B+0D",
           "J563h" = "B+3H",
           "J63"   = "B+7D"
         )) %>% 
  dplyr::select(participant_id, time, all_of(intersection_genes))

df_stacked = bind_rows(
  bind_rows(ebovac2_merged_filtered, hamburg_merged_filtered),
  prevac_merged_filtered
)


df_merged_all = merge(x = df_clinical_all, 
                      y = df_stacked,
                      by = c("participant_id", "time"))

saveRDS(df_merged_all,
        file = fs::path(processed_data_path, "df_merged_all.rds"))

rm(list = ls())
