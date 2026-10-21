###### Respo Code for Seawater Samples - SCRIPT 2 of 2: rates, QC and plots #######
### Created by: Nyssa Silbiger
#### Adapted for seawater on: 2026-06-25 by Jordan Vest & Maya Powell for initial MPW Project. Adapted further by Jordan Vest for OCtober MPW
## Split from respo_seawater2.2.R on 2026-10-21
## Reads the mpw_RespoR_D#.csv files made by the script "respo_LolinR.R"

library(tidyverse)
library(lubridate)
library(here)

RespoMeta <- read_csv(here("data","respirometry","mpw", "mpw_respo_meta.csv"))

######### Calculate Respiration rate ###############
############################################

RespoR <- list.files(here("data","respirometry","mpw","processed"),
                     pattern = "^mpw_RespoR_D[0-6]\\.csv$", full.names = TRUE) %>%
  map_dfr(read_csv)
Sample_Info <- read_csv(here("data","respirometry","mpw","sample_info_mpw.csv"))

RespoR2 <- RespoR %>%
  drop_na(sample_ID) %>%          # J ADDED THIS — drops the leftover phantom rows 
  left_join(Sample_Info) %>% # Join the raw respo calculations with the metadata
  mutate(Ch.Volume.mL = volume_mL) %>% # chamber volume from mpw_measurements.csv
  mutate(Ch.Volume.L = Ch.Volume.mL * 0.001) %>% # mL to L conversion
  mutate(umol.sec = umol.L.sec*Ch.Volume.L) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1
  mutate_if(sapply(., is.character), as.factor) %>% #convert character columns to factors
  mutate(umol.hr = umol.sec*3600) %>% #convert to units per hour
  #mutate(umol.chla.hr = umol.hr/chla) %>% #convert to final units using chla concentrations 
  dplyr::select(date, sample_ID, site_ID, light_dark, run_block, umol.hr, chamber_channel,
              Temp.C)
######@JORDAN DO THIS LATER!!!!##### CHLA CONVERSION ABOVE!!!!

write_csv(RespoR2 , here("data","respirometry","mpw","processed","mpw_RespoR2_AllRates.csv"))  
RespoR2 <- read_csv(here("data","respirometry","mpw","processed","mpw_RespoR2_AllRates.csv"))

###Generate Respo rate dataframe for use in future analysis
RespoR_PR <- RespoR2 %>%
  dplyr::select(-Temp.C) %>% # remove to pivot
  pivot_wider(names_from = light_dark, values_from = umol.hr) %>% 
  rename(Respiration = DARK , NetPhoto = LIGHT) %>% # rename the columns
  mutate(Respiration = -1 * Respiration) %>%  # Make respiration positive
  mutate(GrossPhoto = Respiration + NetPhoto) %>% # calculate gross photosynthesis
  pivot_longer(cols = Respiration:GrossPhoto, names_to = "PR", values_to = "Values") #values still in umol.hr

write_csv(RespoR_PR,here("data","respirometry","mpw","processed","mpw_PnR_rates.csv")) # export all the uptake rates
RespoR_PR <- read_csv(here("data","respirometry","mpw","processed","mpw_PnR_rates.csv"))

#### MPW QC: samples to look at first ####
# O2 going UP in the dark (negative respiration) = probe/bubble/leak problem
qc_dark_up <- RespoR_PR %>% 
  filter(PR == "Respiration", Values < 0) %>% 
  left_join(RespoMeta %>% filter(light_dark == "DARK") %>% dplyr::select(sample_ID, notes)) %>% 
  arrange(run_block, Values)
qc_dark_up   # open the PDFs in output/mpw/respirometry/D#/ for these samples

#### MPW plot: lagoon vs cove vs ocean ####
SampleMeta <- read_csv(here("data","mpw_sample_meta.csv")) %>% 
  dplyr::select(sample_ID, habitat, dist_to_sgd_m, flag)

sgd_cutoff_m <- 10   # MPW: lagoon samples closer than this to the SGD = "near_SGD", the rest of the tidepool = "T"
# (dist_to_sgd_m is from the PLANNED points for now - rerun after GPS update)

# MPW: samples to leave out = anything with text in the "exclude" column of mpw_respo_meta.csv
# to bring a sample back, delete its exclude text in the csv. 
excluded <- RespoMeta %>% filter(!is.na(exclude)) %>% distinct(sample_ID, exclude)

plot_data <- RespoR_PR %>% 
  left_join(SampleMeta, by = "sample_ID") %>% 
  filter(habitat != "blank") %>%                       # put a # at the start of this line to add the blanks (B) back in
  filter(!sample_ID %in% excluded$sample_ID) %>%      # drops samples marked in the "exclude" column of mpw_respo_meta.csv; put a # here to show them
  mutate(habitat = factor(habitat,
                          levels = c("lagoon","cove","ocean","blank"),
                          labels = c("L","C","O","B")),          # short x-axis labels
         # habitat split by SGD distance: near_SGD / T (rest of tidepool) / C / O / B
         habitat_sgd = case_when(habitat == "L" & dist_to_sgd_m <  sgd_cutoff_m ~ "near_SGD",
                                 habitat == "L" & dist_to_sgd_m >= sgd_cutoff_m ~ "T",
                                 TRUE ~ as.character(habitat)),
         habitat_sgd = factor(habitat_sgd, levels = c("near_SGD","T","C","O","B")),
         PR = factor(PR, levels = c("Respiration","NetPhoto","GrossPhoto"),
                     labels = c("Respiration","Net photosynthesis","Gross photosynthesis")))

# which samples were dropped by the QC filter (for your notes)
excluded   # MPW: samples left out of the plot, with the reason from the "exclude" column

PR_plot <- plot_data %>% 
  ggplot(aes(x = habitat_sgd, y = Values)) +   # MPW: change habitat_sgd back to habitat for plain L / C / O
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = run_block), width = 0.12, alpha = 0.8, size = 2) +
  facet_wrap(~PR, nrow = 1) +              # same y-axis on all 3 panels so they compare directly
  scale_x_discrete(drop = TRUE) +
  theme_bw(base_size = 12) +
  labs(x = paste0("Habitat (near_SGD = tidepool < ", sgd_cutoff_m, " m from SGD, T = rest of tidepool, C = cove, O = ocean)"),
       y = expression(paste(mu,"mol O"[2]," hr"^-1)), color = "Day") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
PR_plot

ggsave(here("output","mpw","respirometry","PR_by_habitat_mpw.pdf"),
       device = "pdf", height = 4, width = 9, PR_plot)

###Look at one raw trace: sub in the sample ID to find a spike/bubble and adjust its times in the metadata
check_day <- "D1"   # MPW: the day folder the file is in
raw_check <- read_csv(skip = 1, here("data","respirometry","mpw","rawo2", check_day, "MPW_I01_O2.csv")) %>%   # swap in the file name (must be in the current day folder)
  dplyr::select(Date, Time, Value, Temp) %>%
  unite(Date, Time, col = "Time", sep = " ") %>%
  drop_na() %>% 
  mutate(Time = mdy_hms(Time))

par(mfrow = c(2,1))
plot(raw_check$Time, raw_check$Value, type = "l", main = "O2 (full raw trace)")
plot(raw_check$Time, raw_check$Temp,  type = "l", main = "Temp (full raw trace)")
par(mfrow = c(1,1))
