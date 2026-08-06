###### Respo Code for Seawater Samples Light and Dark Runs ####### 
### Created by: Nyssa Silbiger
#### Adapted for seawater on: 2026-06-25 by Jordan Vest & Maya Powell for MPW Project
##### Adapted for Reef Water Column Productivity Project 2026-08-05 by Jordan Vest
## Last updated on: 2026-08-05

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
path.p<-here("data","respirometry","kewalo", "rawo2", "K_RUN2") #the location of all your respirometry files
#you can change to individual run folders if needed
#use this to check file path is good before running script 
#print(path.p)
#print(filenames_final)

# bring in all of the individual files
filenames_final<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders

#basename above removes the subdirectory name from the file, re-name as file.names.full
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#empty chamber volume
ch.vol <- 350 #mL #of small chambers 
###J to update if chambers change

######### Load and tidy files ###############
############################################
#Load your respiration data file, with all the times, water volumes(mL), #not doing dry weight just SA
#RespoMeta <- read_csv(here("Data","RespoFiles","Respo_Metadata_SGDDilutions_Cabral_Varari.csv"))
BioData <- read_csv(here("data","respirometry", "kewalo", "measurements_kewalo.csv")) #ch vol, sa, chla, etc things to normalize to

RespoMeta <- read_csv(here("data","respirometry","kewalo", "kewalo_respo_meta.csv")) #metadata, ID, site, etc
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
write_csv(Sample_Info, here("data","respirometry","kewalo","sample_info_kewalo.csv"))
Sample_Info <- read_csv(here("data","respirometry","kewalo","sample_info_kewalo.csv"))

#generate a dataframe with specific column names
# data is in umol.L.sec

RespoR <- tibble(.rows =length(filenames_final)*2, #*2 for light and dark runs
                 sample_ID = NA,
                 Intercept = NA,
                 umol.L.sec = NA,
                 Temp.C = NA,
                 run_block = NA,
                 light_dark = NA)

######### Create a for loop! ###############
############################################

###forloop##### 
for(i in 1:length(filenames_final)) {
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
        filter(t_sec > 120) %>%                         # drop first 2 min (120 s)
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
  
  for(j in 1:length(df$light)){
    pdf(paste0(here("output","kewalo", "respirometry"),"/",rename,"_",j,".pdf" ))
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
  
}  

#View(RespoR)
######### end of for loop - celebrate victory of getting through that ###############
############################################

#export raw data and read back in as a failsafe 
#this allows me to not have to run the for loop again !!!!!
write_csv(RespoR, here("data","respirometry","kewalo","kewalo_RespoR.csv"))  

##### 

######### Calculate Respiration rate ###############
############################################

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

