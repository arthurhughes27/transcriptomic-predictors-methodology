# Script to perform basic study descriptions of prevac, hamburg, ebovac2

# Libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Load the harmonised clinical data

processed_data_path <- fs::path("data")

df_clinical_all = readRDS(fs::path(processed_data_path, "df_clinical_all.rds")) %>% 
  filter(study_vaccine != "ebovac2-placebo")

# Describe the antibody timepoints per study

# Desired chronological order: prime then boost
time_order <- c(
  "p_0","p_7","p_14","p_28","p_56","p_63","p_84","p_180","p_365",
  "b_0","b_21"
)

df_counts <- df_clinical_all %>%
  pivot_longer(
    cols = starts_with("ab_"),
    names_to = "ab_time",
    values_to = "ab_value"
  ) %>%
  mutate(
    phase = str_extract(ab_time, "(?<=ab_)[pb]"),
    day = str_extract(ab_time, "\\d+$"),
    timepoint = paste0(phase, "_", day),
    timepoint = factor(timepoint, levels = time_order)
  ) %>%
  filter(!is.na(ab_value)) %>%
  group_by(study_vaccine, timepoint) %>%
  summarise(
    n_participants = n_distinct(participant_id),
    .groups = "drop"
  ) %>%
  complete(
    study_vaccine,
    timepoint,
    fill = list(n_participants = 0)
  )

# Custom labels for x-axis
label_fun <- function(x) {
  ifelse(
    str_detect(x, "^p_"),
    paste0("Prime + ", str_remove(x, "^p_"), " Days"),
    paste0("Boost + ", str_remove(x, "^b_"), " Days")
  )
}

p1 = ggplot(df_counts, aes(x = timepoint, y = study_vaccine, fill = n_participants > 0)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = n_participants), size = 5) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#3BA239"), guide = "none") +
  scale_x_discrete(labels = label_fun) +
  labs(
    x = "Timepoint",
    y = "Study-Group",
    title = "Availability of Antibody Measurements by Study-Group and Timepoint"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

p1

descriptive_figures_folder = fs::path("output", "figures", "descriptive")

ggsave(
  filename = "ebola_antibody_availability.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 39,
  height = 18,
  units = "cm"
)

# Explicit chronological order of GE timepoints
ge_time_order <- c(
  "P+0D", "P+3H", "P+1D", "P+3D", "P+7D", "P+29D", "P+57D", "P+85D",
  "B+0D", "B+3H", "B+1D", "B+7D"
)

# Count distinct participants with GE samples
df_counts_ge <- df_clinical_all %>%
  filter(!is.na(time)) %>%
  mutate(
    timepoint = factor(time, levels = ge_time_order)
  ) %>%
  group_by(study_vaccine, timepoint) %>%
  summarise(
    n_participants = n_distinct(participant_id),
    .groups = "drop"
  ) %>%
  complete(
    study_vaccine,
    timepoint,
    fill = list(n_participants = 0)
  )

# Custom x-axis labels
label_fun_ge <- function(x) {
  x <- as.character(x)
  ifelse(
    str_detect(x, "^P\\+"),
    paste0("Prime + ", str_remove(x, "^P\\+")),
    paste0("Boost + ", str_remove(x, "^B\\+"))
  )
}

# Create GE availability figure
p2 <- ggplot(df_counts_ge, aes(x = timepoint, y = study_vaccine, fill = n_participants > 0)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = n_participants), size = 5) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#3BA239"), guide = "none") +
  scale_x_discrete(labels = label_fun_ge) +
  labs(
    x = "Timepoint",
    y = "Study-Group",
    title = "Availability of Gene Expression Samples by Study-Group and Timepoint"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 30),
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

p2

# Save figure
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")
ggsave(
  filename = "ebola_GE_availability.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 39,
  height = 18,
  units = "cm"
)
