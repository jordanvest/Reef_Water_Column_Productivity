######### Calculate Respiration rate ###############
############################################

#Install Libraries
library(tidyverse)
library(here)

RespoR <- read_csv(here("data","respirometry","kewalo","kewalo_RespoR.csv"))
Sample_Info <- read_csv(here("data","respirometry","kewalo","sample_info_kewalo.csv"))
ch.vol <- 350 #mL #of small chambers 

RespoR2 <- RespoR %>%
  #drop_na(FileID_csv) %>% # drop NAs
  left_join(Sample_Info) %>% # Join the raw respo calculations with the metadata
  mutate(Ch.Volume.mL = ch.vol) %>% # measured volume of chambers with coral + stand + stirbar displacement
  mutate(Ch.Volume.L = Ch.Volume.mL * 0.001) %>% # mL to L conversion
  mutate(umol.sec = umol.L.sec*Ch.Volume.L) %>% #Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
  mutate_if(sapply(., is.character), as.factor) %>% #convert character columns to factors
  mutate(umol.hr = umol.sec*3600) %>% #convert to units per hour
  #mutate(umol.chla.hr = umol.hr/chla) %>% #convert to final units using chla concentrations 
  dplyr::select(date, sample_ID, site_ID, light_dark, run_block, run_block, umol.hr, chamber_channel, 
                Temp.C) #keep only what we need
######@JORDAN DO THIS LATER!!!!##### CHLA CONVERSION ABOVE!!!!

write_csv(RespoR2 , here("data","respirometry","kewalo","RespoR2_AllRates.csv"))  
RespoR2 <- read_csv(here("data","respirometry","kewalo","RespoR2_AllRates.csv"))

###Generate Respo rate dataframe for use in future analysis
RespoR_PR <- RespoR2 %>%
  dplyr::select(-Temp.C) %>% # remove to pivot
  pivot_wider(names_from = light_dark, values_from = umol.hr) %>% 
  rename(Respiration = DARK , NetPhoto = LIGHT) %>% # rename the columns
  mutate(Respiration = -1 * Respiration) %>%  # Make respiration positive
  mutate(GrossPhoto = Respiration + NetPhoto) %>% # calculate gross photosynthesis
  pivot_longer(cols = Respiration:GrossPhoto, names_to = "PR", values_to = "Values") #values still in umol.hr

write_csv(RespoR_PR,here("data","respirometry","kewalo","PnR_rates.csv")) # export all the uptake rates
RespoR_PR <- read_csv(here("data","respirometry","kewalo","PnR_rates.csv"))

#dev.off() # may need if plot doesn't run?
PR_plot <- RespoR_PR %>% 
  ggplot(aes(x = site_ID, y = Values, group=site_ID, color = site_ID, shape = run_block)) +
  geom_point() +
  geom_line() +
  facet_wrap(run_block~PR, scales = "free") +
  theme_bw() +
  labs(x = "Site", y = "umol.hr")+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

ggsave(here("output", "kewalo","respirometry","PR_boxplots_kewalo.pdf"),
       device = "pdf", height = 8, width = 8, PR_plot)

