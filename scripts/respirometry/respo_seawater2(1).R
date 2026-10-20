###### Respo Code for Seawater Samples Light and Dark Runs ####### 
### Created by: Nyssa Silbiger
#### Adapted for seawater on: 2026-06-25 by Jordan Vest & Maya Powell for MPW Project
##### Adapted for Reef Water Column Productivity Project 2026-08-05 by Jordan Vest
## Last updated on: 2026-10-19 (adapted for MPW Kapo'o campaign, run one day at a time)

############## Introduction to code/script ####################
## this script will help us process the raw data gathered during respirometry runs. 
## need to change for specific project/experimental variables 

### Install Packages #####
## if these packages are not yet installed, install them 
## great for updates or new users 
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented')
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix')
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra')
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') devtools::install_github('colin-olito/LoLinR')
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron')
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse')
if ("here" %in% rownames(installed.packages()) == 'FALSE') install.packages('here')
if ("patchwork" %in% rownames(installed.packages()) == 'FALSE') install.packages('patchwork')
if ("nls.multstart" %in% rownames(installed.packages()) == 'FALSE') install.packages('nls.multstart')

#Read in required libraries
##### Include Versions of libraries
library(segmented)
library(plotrix)
library(gridExtra)
library(LoLinR)
library(lubridate)
library(chron)
library(patchwork)
library(tidyverse)
library(here)
library(PNWColors)
library(ggrepel)
library(reshape2)
library(viridis)
library(car)
library(future)
library(furrr)
library(dplyr)

############# now it's time to code ############
################################################
# get the file path

#set the path to all of the raw oxygen datasheets
## these are saved onto the computer in whatever file path/naming scheme you saved things to 
day <- "D1"   # MPW: change this ONE line to D1, D2 ... D6 and run the script top to bottom, one day at a time
path.p<-here("data","respirometry","mpw", "rawo2", day) #the location of all your respirometry files
#use this to check file path is good before running script
#print(path.p)


# bring in all of the individual files
# MPW: recursive = FALSE so the D5 "1hrDark" and "export2" subfolders are ignored
filenames_final<-basename(list.files(path = path.p, pattern = "csv$", recursive = FALSE)) #list all csv file names in the day folder
file.names.full<-filenames_final
#print(filenames_final)
# chamber volume now comes from volume_mL in mpw_measurements.csv (473.1765 mL)

######### Load and tidy files ###############
############################################
#Load your respiration data file, with all the times, water volumes(mL), #not doing dry weight just SA
#RespoMeta <- read_csv(here("Data","RespoFiles","Respo_Metadata_SGDDilutions_Cabral_Varari.csv"))
BioData <- read_csv(here("data","respirometry", "mpw", "mpw_measurements.csv")) #ch vol, sa, chla, etc things to normalize to

RespoMeta <- read_csv(here("data","respirometry","mpw", "mpw_respo_meta.csv")) #metadata, ID, site, etc
#View(BioData)
#View(RespoMeta)

# join the data together - use these if you need to take out any columns for them to match
# BioData <- BioData %>% 
#   dplyr::select(-full_species, -species_ID)
# 
# RespoMeta <- RespoMeta %>% 
#   dplyr::select(-notes)

Sample_Info <- left_join(RespoMeta, BioData, by = c("sample_ID", "site_ID"))
#for water samples - "site ID" is the site ID with "sample ID" being each individual sample replicate
#View(Sample_Info)

##### Make sure times are consistent ####
# make start and stop times real times, so that we can join the respo output and sample_info data frames
Sample_Info <- Sample_Info %>% 
  #drop_na(sample_ID) %>% 
  unite(date,start_time,col="start_time",remove=F, sep=" ") %>% 
  unite(date,stop_time,col="stop_time",remove=F, sep=" ") %>%
  mutate(start_time = mdy_hms(start_time)) %>% 
  mutate(stop_time = mdy_hms(stop_time)) %>% 
  mutate(date = mdy(date))

#View(Sample_Info)
write_csv(Sample_Info, here("data","respirometry","mpw","sample_info_mpw.csv"))
Sample_Info <- read_csv(here("data","respirometry","mpw","sample_info_mpw.csv"))
Sample_Info <- Sample_Info %>% filter(run_block == day)   # J added this line 9/15/26

# MPW: only process files that have metadata, and say what gets skipped (e.g. D6 MPW_O55_O2.csv)
setdiff(filenames_final, Sample_Info$FileID_csv)   # files in the folder with NO metadata -> skipped
setdiff(Sample_Info$FileID_csv, filenames_final)   # metadata rows with NO file -> check these!
filenames_final <- filenames_final[filenames_final %in% Sample_Info$FileID_csv]
file.names.full <- filenames_final

#generate a dataframe with specific column names
# data is in umol.L.sec
RespoR <- tibble(.rows =nrow(Sample_Info), # MPW: one row per metadata row (DARK + LIGHT)
                 sample_ID = NA,
                 Intercept = NA,
                 umol.L.sec = NA,
                 Temp.C = NA,
                 run_block = NA,
                 light_dark = NA)

######### Create a for loop! ###############
############################################

###forloop##### 
failed_files <- c()   # MPW: files that error get listed here instead of stopping the whole day
for(i in 1:length(filenames_final)) {
  cat(i, "of", length(filenames_final), ":", filenames_final[i], format(Sys.time(), "%H:%M"), "\n")  # MPW: shows which file it is on, so you can tell a slow file from a stuck one
  tryCatch({   # MPW: if a file errors, note it and move on to the next one
  FRow <- as.numeric(which(Sample_Info$FileID_csv==filenames_final[i])) # stringsplit this renames our file
  
  Respo.Data1 <- read_csv(skip=1,file.path(path.p, paste0(file.names.full[i]))) %>% # reads in each file in list
    dplyr::select(Date, Time, Value, Temp) %>% # keep only what we need: Time stamp, Raw O2 value, in situ temp - all at time interval of 2sec or whatever you set it to
    unite(Date,Time,col="Time",remove=T, sep = " ") %>%
    drop_na() %>% 
    mutate(Time = mdy_hms(Time)) #%>% # convert time
    #mutate(help = i) ##if stuck in forloop with error from filter, can check RespoR and see at what row the forloop stopped working  
  
  ## cut the data by start and stop times from metadata
  #Use start time of each light step from the metadata to separate data by light stop
  
  oxy_subsets <- Sample_Info[FRow,] %>%
    pmap(function(light_dark, start_time, stop_time, ...) {
      data <- Respo.Data1  %>%
        filter(Time >= start_time & Time <= stop_time) %>%
        arrange(Time) %>%
        mutate(t_sec = as.numeric(difftime(Time, first(Time), units = "secs"))) %>% #keep everything in seconds
        mutate(light_dark = light_dark) %>%
        filter(t_sec > 1800) %>%                         # drop first 30 min (1800 s)
        filter(row_number() %% 10 == 0)                  # keep every 10th row - @Jordan not doing this bc we did every 2min of data here
      # now t_sec increments by ~20 s for your kept rows
      #   filter(Time >= start_time & Time <= stop_time) %>%
      #   mutate(sec = row_number()) %>%# add an id for each row to help remove the first few mins
      #   mutate(light_dark = light_dark,
      #          temp_c_value = temp_c_value,
      #          sec = sec) %>%
      #   filter(sec > 60)  %>%# delete the first 2 mins of data assuming freq of 2 Hz
      #   mutate(row_number = row_number()) %>%
      #   filter(row_number %% 10 == 0) %>%  # keep every 10th row only to thin the data
      #   dplyr::select(-row_number) %>%
      #   mutate(sec2 = row_number())  #update the row numbers
      # #return(subset)
    }) 
  
  
  # Combine into one long dataframe with ID labels
  combined_oxy <- bind_rows(oxy_subsets)
  
  # Get the filename without the .csv
  rename<- sub("_O2.csv","", filenames_final[i])
  
  ### plot and export the thinned data ####
  p1<- ggplot(combined_oxy, aes(x = Time, y = Value)) +
    geom_point(aes(color = light_dark)) +
    labs(
      x = 'Time (seconds)',
      y = expression(paste(' O'[2],' (',mu,'mol/L)')),
      title = "original")
  
  
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "e "pc")
  
  # Define function for fitting LoLinR regressions to be applied to all intervals for all samples
  fit_reg <- function(data) {
    rankLocReg(xall = data$t_sec, yall = data$Value,
               alpha = 0.2, method = "pc", verbose = FALSE)
  }
  
  # Setup for parallel processing
  future::plan(multisession)
  
  # Map LoLinR function onto all intervals of each sample's thinned dataset
  df <- combined_oxy %>%
    dplyr::select(t_sec, Value, Temp, light_dark)%>%
    #mutate(sec2 = as.numeric(sec2))%>%
    nest_by(light_dark) %>%
    ungroup()%>%
    mutate(regs = furrr::future_map(data, fit_reg), # run the LOLinR fit in parallel
           Temp.C = map_dbl(map(data, "Temp"), mean),# get the mean temperature
           RegStats =map(regs, function(x){ # extract the intercept and slope for the parameters
             x$allRegs %>%
               slice(1) %>%
               dplyr::select(Intercept = b0,
                      umol.L.sec = b1)
           }) )
  
  
  #  Plot regression diagnostics
  df <- df %>% 
    mutate(light = paste(light_dark))
  
  out.dir <- here("output","mpw","respirometry", day)   # MPW: one PDF folder per day
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  for(j in 1:length(df$light)){
    pdf(paste0(out.dir,"/",rename,"_",df$light_dark[j],".pdf" ))   # file name says DARK/LIGHT instead of 1/2
    plot(df$regs[[j]])
    dev.off() 
  }
  
  df<- df %>%
    dplyr::select(Temp.C,light_dark,RegStats ) %>%
    unnest(RegStats) %>%
    mutate(sample_ID = rename,) %>%
    left_join(Sample_Info[FRow,] %>%
                dplyr::select(site_ID, sample_ID, run_block, light_dark))  # make sure the light value (or whatever other metadata you want) is in the final DF
  
  ################################
  # fill in all the O2 consumption and rate data
  
  RespoR[FRow,"Temp.C"]<-df$Temp.C
  RespoR[FRow, "light_dark"]<-df$light_dark
  RespoR[FRow,"sample_ID"]<-df$sample_ID
  RespoR[FRow,"Intercept"]<-df$Intercept
  RespoR[FRow,"umol.L.sec"]<-df$umol.L.sec
  RespoR[FRow,"run_block"]<-df$run_block
  
  }, error = function(e) {
    message("FAILED: ", filenames_final[i], " -> ", conditionMessage(e))
    failed_files <<- c(failed_files, filenames_final[i])
  })
}
failed_files   # MPW: anything listed here needs a look (open the raw trace with raw_check at the bottom)

#View(RespoR)
######### end of for loop - celebrate victory of getting through that ###############
############################################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
# MPW: one file per day, so running D2 doesn't overwrite D1
dir.create(here("data","respirometry","mpw","processed"), showWarnings = FALSE)
write_csv(RespoR, here("data","respirometry","mpw","processed", paste0("mpw_RespoR_", day, ".csv")))

##### 

######### Calculate Respiration rate ###############
############################################
# MPW: this section stacks every day you've run so far (all mpw_RespoR_D*.csv files)
# so you can run it after any day and it will include all finished days

RespoR <- list.files(here("data","respirometry","mpw","processed"),
                     pattern = "^mpw_RespoR_D[0-9]\\.csv$", full.names = TRUE) %>%
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
              Temp.C) # MPW: removed chamber_type (all chambers the same)
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
# O2 going UP in the dark (negative respiration) = probe/bubble/leak problem, not biology
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

plot_data <- RespoR_PR %>% 
  left_join(SampleMeta, by = "sample_ID") %>% 
  filter(habitat != "blank") %>%                       # put a # at the start of this line to add the blanks (B) back in
  filter(!sample_ID %in% qc_dark_up$sample_ID) %>%     # drops QC samples (all 3 rates); put a # here to show them
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
qc_dark_up %>% dplyr::select(run_block, sample_ID, Values, notes)

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
raw_check <- read_csv(skip = 1, file.path(path.p, "MPW_I01_O2.csv")) %>%   # swap in the file name (must be in the current day folder)
  dplyr::select(Date, Time, Value, Temp) %>%
  unite(Date, Time, col = "Time", sep = " ") %>%
  drop_na() %>% 
  mutate(Time = mdy_hms(Time))

par(mfrow = c(2,1))
plot(raw_check$Time, raw_check$Value, type = "l", main = "O2 (full raw trace)")
plot(raw_check$Time, raw_check$Temp,  type = "l", main = "Temp (full raw trace)")
par(mfrow = c(1,1))
