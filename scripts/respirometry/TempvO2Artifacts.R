###### Temp vs O2 Artifact Check ####### 
### Jordan Vest 2026-09-18
############## Introduction to script ####################
### this script is a standalone check on if rising O2 during LIGHT periods is real productivity or a temperature-driven artifact 
### Not part of the main respo_seawater.R pipeline — run this separately

#Read in required libraries
library(tidyverse)
library(here)

# ---- Settings — change these for whichever run you're checking ----
# get the file path
#set the path to all of the raw oxygen datasheets
## these are saved onto the computer in whatever file path/naming scheme you saved things to 
path.p <- here("data","respirometry","kewalo","rawo2","K_RUN4")
#use this to check file path is good before running script 
#print(path.p)
acclim_sec <- 1800   # same 30-min acclimation cutoff used in the main pipeline
# keep these in sync if you ever change one

# bring in all of the individual files
filenames <- list.files(path.p, pattern = "csv$")

# ---- Read + prep one raw file ----
read_raw <- function(filename) {
  read_csv(skip = 1, file.path(path.p, filename), show_col_types = FALSE) %>%
    dplyr::select(Date, Time, Value, Temp) %>%
    unite(Date, Time, col = "Time", sep = " ") %>%
    mutate(Time = mdy_hms(Time)) %>%
    drop_na(Time, Value, Temp) %>%                  # fixes the ccf() "missing values" error
    arrange(Time) %>%
    mutate(t_sec = as.numeric(difftime(Time, first(Time), units = "secs"))) %>%
    filter(t_sec > acclim_sec) %>%                  # cut the same 30-min acclimation window
    mutate(sample_ID = sub("_O2.csv", "", filename))
}

# ---- Run for every file in the folder, combine into one dataframe ----
all_raw <- map_dfr(filenames, read_raw)

# ---- Quick correlation summary per sample (console) ----
# high |correlation| (roughly >0.8) = temperature likely tracking O2 closely
cor_summary <- all_raw %>%
  group_by(sample_ID) %>%
  summarise(cor_temp_o2 = cor(Temp, Value, use = "complete.obs")) %>%
  arrange(desc(abs(cor_temp_o2)))

print(cor_summary)

# ---- Combined plot: O2 and Temp overlaid (standardized), one panel per sample ----
# Raw units don't share a scaleso both are converted to z-scores per sample just for this visual check
plot_data <- all_raw %>%
  group_by(sample_ID) %>%
  mutate(O2_z = as.numeric(scale(Value)),
         Temp_z = as.numeric(scale(Temp))) %>%
  ungroup() %>%
  pivot_longer(cols = c(O2_z, Temp_z), names_to = "variable", values_to = "z_value")

temp_o2_plot <- ggplot(plot_data, aes(x = Time, y = z_value, color = variable)) +
  geom_line() +
  facet_wrap(~sample_ID, scales = "free_x") +
  theme_bw() +
  labs(y = "Standardized value (z-score)", x = "Time",
       title = "O2 vs Temp — quick check for temp artifacts") +
  scale_color_manual(values = c("O2_z" = "steelblue", "Temp_z" = "firebrick"),
                     labels = c("O2_z" = "O2", "Temp_z" = "Temp"),
                     name = NULL)

print(temp_o2_plot)

ggsave(here("output","kewalo","respirometry","temp_o2_check.pdf"),
       device = "pdf", width = 10, height = 8, temp_o2_plot)
############################################
# ---- Individual plot for O2 and Temp; two panels per sample ----
plot_data_raw <- all_raw %>%
  pivot_longer(cols = c(Value, Temp), names_to = "variable", values_to = "raw_value")

ggplot(plot_data_raw, aes(x = Time, y = raw_value)) +
  geom_line() +
  facet_grid(variable ~ sample_ID, scales = "free") +   # rows = O2/Temp, columns = sample
  theme_bw() +
  labs(y = NULL, x = "Time",
       title = "O2 and Temp — literal values") 

# ---- Lag check per sample (console + base R plots) ----
# ccf() doesn't have a built-in "group by" version, so this loops manually.
# Spike right at lag 0 = O2 and Temp moving together instantly (physical coupling).
# Spike offset from 0 = O2 responds after Temp, more consistent with a real
# biological response layered on top.
for (id in unique(all_raw$sample_ID)) {
  sub <- all_raw %>% filter(sample_ID == id)
  cat("\n---", id, "---\n")
  ccf(sub$Temp, sub$Value, main = paste(id, "- lag check"))
}

# ------Filtering for LIGHT run ---------#
###### Temp vs O2 Artifact Check ####### 
### Standalone script — checks whether rising O2 during LIGHT periods
### is real productivity or a temperature-driven artifact (e.g. water
### bath drift, ice pack additions).
### Not part of the main respo_seawater.R pipeline — run this separately
### on any run's raw files.

library(tidyverse)
library(here)
library(patchwork)   # lets us arrange multiple ggplots into one grid

# ---- Settings — change these for whichever run you're checking ----
path.p <- here("data","respirometry","kewalo","rawo2","K_RUN4")
acclim_sec <- 1800   # same 30-min acclimation cutoff used in the main pipeline
# keep these in sync if you ever change one

filenames <- list.files(path.p, pattern = "csv$")

# ---- Load metadata so we know each sample's DARK/LIGHT time windows ----
# This is what lets us tag/filter by light_dark below — the raw files
# themselves have no concept of "dark" or "light", only timestamps.
RespoMeta <- read_csv(here("data","respirometry","kewalo","kewalo_respo_meta.csv"),
                      show_col_types = FALSE) %>%
  filter(run_block == basename(path.p)) %>%
  mutate(start_time = mdy_hms(paste(date, start_time)),
         stop_time  = mdy_hms(paste(date, stop_time)))
# Note: any sample you've excluded from this file (like S2.3_RUN4) simply
# won't have rows here, so it'll automatically drop out of everything below.

# ---- Read + prep one raw file, split into its DARK/LIGHT chunks ----
read_raw <- function(filename) {
  raw <- read_csv(skip = 1, file.path(path.p, filename), show_col_types = FALSE) %>%
    dplyr::select(Date, Time, Value, Temp) %>%
    unite(Date, Time, col = "Time", sep = " ") %>%
    mutate(Time = mdy_hms(Time)) %>%
    drop_na(Time, Value, Temp) %>%                  # fixes the ccf() "missing values" error
    arrange(Time)
  
  meta_rows <- RespoMeta %>% filter(FileID_csv == filename)
  
  # pmap_dfr runs this once per metadata row (once for DARK, once for LIGHT)
  # and stacks the results — same technique the main pipeline uses.
  meta_rows %>%
    pmap_dfr(function(light_dark, start_time, stop_time, ...) {
      raw %>%
        filter(Time >= start_time & Time <= stop_time) %>%
        mutate(t_sec = as.numeric(difftime(Time, first(Time), units = "secs"))) %>%
        filter(t_sec > acclim_sec) %>%               # cut acclimation window within THIS block
        mutate(light_dark = light_dark,
               sample_ID = sub("_O2.csv", "", filename))
    })
}

# ---- Run for every file in the folder, combine into one dataframe ----
all_raw <- map_dfr(filenames, read_raw)

# ---- Split by light/dark ----
light_only <- all_raw %>% filter(light_dark == "LIGHT")
dark_only  <- all_raw %>% filter(light_dark == "DARK")   # here in case you want it later

# ---- Correlation summary — LIGHT period only this time ----
cor_summary_light <- light_only %>%
  group_by(sample_ID) %>%
  summarise(cor_temp_o2 = cor(Temp, Value, use = "complete.obs")) %>%
  arrange(desc(abs(cor_temp_o2)))

print(cor_summary_light)

# ---- Scaled overlay plot, LIGHT only, all samples faceted ----
plot_data <- light_only %>%
  group_by(sample_ID) %>%
  mutate(O2_z = as.numeric(scale(Value)),
         Temp_z = as.numeric(scale(Temp))) %>%
  ungroup() %>%
  pivot_longer(cols = c(O2_z, Temp_z), names_to = "variable", values_to = "z_value")

temp_o2_plot_light <- ggplot(plot_data, aes(x = Time, y = z_value, color = variable)) +
  geom_line() +
  facet_wrap(~sample_ID, scales = "free_x") +
  theme_bw() +
  labs(y = "Standardized value (z-score)", x = "Time",
       title = "O2 vs Temp (standardized), LIGHT period only") +
  scale_color_manual(values = c("O2_z" = "steelblue", "Temp_z" = "firebrick"),
                     labels = c("O2_z" = "O2", "Temp_z" = "Temp"),
                     name = NULL)

print(temp_o2_plot_light)

ggsave(here("output","kewalo","respirometry","temp_o2_check_light.pdf"),
       device = "pdf", width = 10, height = 8, temp_o2_plot_light)

# ---- Dual-axis literal-values plot for ONE sample ----
# Scales the two axes specifically to that sample's own ranges, so O2
# and TEMP line up meaningfully — this only works cleanly per-sample,
# not faceted across all of them (different samples have different ranges).
plot_temp_o2_dual <- function(sample_id, data = light_only) {
  d <- data %>% filter(sample_ID == sample_id)
  if (nrow(d) == 0) return(NULL)
  
  o2_range <- range(d$Value)
  temp_range <- range(d$Temp)
  scale_factor <- diff(o2_range) / diff(temp_range)
  temp_offset <- o2_range[1] - temp_range[1] * scale_factor
  
  ggplot(d, aes(x = Time)) +
    geom_line(aes(y = Value), color = "steelblue") +
    geom_line(aes(y = Temp * scale_factor + temp_offset), color = "firebrick") +
    scale_y_continuous(
      name = "O2 (umol/L)",
      sec.axis = sec_axis(~ (. - temp_offset) / scale_factor, name = "Temp (°C)")
    ) +
    theme_bw() +
    labs(title = paste(sample_id, "- O2 (blue) vs Temp (red)"), x = "Time")
}

# ---- Generate the dual-axis plot for EVERY sample using map() ----
# map() runs plot_temp_o2_dual() once per sample_ID and collects the
# results into a list of plot objects (nothing is drawn on screen yet).
sample_ids <- unique(light_only$sample_ID)
all_dual_plots <- map(sample_ids, plot_temp_o2_dual)

# Option A: one combined grid, all samples visible on screen at once
wrap_plots(all_dual_plots, ncol = 2)

# Option B: multi-page PDF, one sample per page (walk() = map(), but for
# side effects like printing/saving, where we don't need the return value)
pdf(here("output","kewalo","respirometry","temp_o2_dual_all_samples.pdf"), width = 8, height = 5)
walk(all_dual_plots, print)
dev.off()

# ---- Lag check per sample (console + base R plots), LIGHT only ----
for (id in unique(light_only$sample_ID)) {
  sub <- light_only %>% filter(sample_ID == id)
  cat("\n---", id, "---\n")
  ccf(diff(sub$Temp), diff(sub$Value), main = paste(id, "- lag check (differenced)"))
}
