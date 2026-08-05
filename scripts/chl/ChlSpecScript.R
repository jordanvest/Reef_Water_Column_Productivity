########################################################################
### This calculates the Chlorophyll from a 96 well plate in corals from teh SpectraMax iD5
### 
### Created by Dr. Nyssa Silbiger
### Edited on 06/01/2026 by Maya Powell (updated chla calculation correctly accounting for path length)
#########################################################################

## create a folder for your day of sampling (foldername below) and put your plate.csv files in there. All your data will be exported to that folder

#libraries
library(here)
library(tidyverse)


## File names -------------------
foldername<-'ChlaSpec' # folder of the day
filename<-'mpw_chlatest.csv' # data
sampleID<-'Chl_template_plate1.csv' # template of sample IDs
platename<-'mpw_chla_plate1' # this will be the name of your file
metadata_file<-"Chl_metadata_plate1.csv" # file with slurry vols and SA

## What is the path length?
PL<-0.71 # pathlength for donahue lab plate


# DONT CHANGE ANYTHING BELOW HERE ----------------------------------

# Read in the Sample IDs
sampleIDNames<-read_csv(here('Data',foldername,sampleID))
sampleIDNames<-sampleIDNames %>%
  select(Well.Location, Sample.Name) # pull out the location and the IDs

# read in metadata
metadata<-read_csv(here("Data", foldername, metadata_file))

### format the 96 well plate data ################
#read in the rows w/o the dye
ChlPlate<-read.csv(here('Data',foldername,filename), row.names = c("A","B","C","D","E","F","G","H"),nrows = 8, fileEncoding="latin1", skip = 2)
# make the rownames a column
ChlPlate$Rows<-rownames(ChlPlate)

## Pull out each of the wavelengths to make own column
#630 wavelength
ChlPlate_630<-ChlPlate %>%
  select(Rows,X1:X12)%>%
  pivot_longer(cols = X1:X12, values_to =  "ChlPlate_630", names_to = "Column")

#663 wavelength
ChlPlate_663<-ChlPlate %>%
  select(Rows,X1.1:X12.1)%>% # pull out the 663 columns
  pivot_longer(cols = X1.1:X12.1, values_to =  "ChlPlate_663", names_to = "Column") %>%
  separate(col = Column, extra = "drop", into = 'Column') # remove the .1

#750 wavelength
ChlPlate_750<-ChlPlate %>%
  select(Rows,X1.2:X12.2)%>% # pull out the 750 columns
  pivot_longer(cols = X1.2:X12.2, values_to =  "ChlPlate_750", names_to = "Column") %>%
  separate(col = Column, extra = "drop", into = 'Column') # remove the .2


## bring everything to one dataframe (reduce allows me to use the left_join function multiple times at once)
AllData<-Reduce(function(...) left_join(...), list(ChlPlate_630,ChlPlate_663,ChlPlate_750))

# remove the X in the column name and add a 0 in front of single digits to keep the well name similar to how the instrument exports
AllData<-AllData %>%
  mutate(Column = as.numeric(str_remove(AllData$Column, pattern = "X")))%>% # remove the X
  mutate(Column = ifelse(Column<10, gsub("(\\d)+", "0\\1", Column),Column))%>% # only add a 0 in front of 1-9
  mutate(Well.Location = paste0(Rows,Column)) # paste with row name

AllData <- AllData %>%
  left_join(sampleIDNames, by = "Well.Location") %>%
  mutate(is_blank = str_detect(tolower(coalesce(Sample.Name, "")), "blank"))

blank_vals <- AllData %>%
  filter(is_blank) %>%
  summarise(
    blank_630 = mean(ChlPlate_630, na.rm = TRUE),
    blank_663 = mean(ChlPlate_663, na.rm = TRUE),
    blank_750 = mean(ChlPlate_750, na.rm = TRUE)
  )

AllData <- AllData %>%
  mutate(
    ChlPlate_630_corr = ChlPlate_630 - blank_vals$blank_630,
    ChlPlate_663_corr = ChlPlate_663 - blank_vals$blank_663,
    ChlPlate_750_corr = ChlPlate_750 - blank_vals$blank_750
  )

### Run chla analysis ###

# chl from Jeffry and Humphreys
ChlData_raw<-AllData %>%
  mutate(chla = 11.43*((ChlPlate_663 - ChlPlate_750)/PL) - 0.64*((ChlPlate_630 - ChlPlate_750)/PL),
         chlc = 27.09*((ChlPlate_630 - ChlPlate_750)/PL) - 3.63*((ChlPlate_663 - ChlPlate_750)/PL),
         Totalchl = chla+chlc,
           daterun = Sys.Date()) %>%
  left_join(sampleIDNames) %>% # join with the sampleIDs
  drop_na(Sample.Name) %>%
  left_join(metadata) %>%
  mutate(chla_ug_mL = (chla/Filtered_vol_ml),  # normalize to filtered volume
         chlc_ug_mL = (chlc/Filtered_vol_ml),
         Totalchl_ug_mL = (Totalchl/Filtered_vol_ml)) %>%
  filter(!is_blank)
         

# averaged data by sample ID
ChlData_ave<-ChlData_raw %>%
  group_by(Sample.Name) %>%
  summarise_at(vars(chla:Totalchl, chla_ug_mL:Totalchl_ug_mL), .funs = list(function(x){mean(x, na.rm = TRUE)},
                                                 function(x){sd(x, na.rm = TRUE)}/sqrt(n()))) %>%
  rename(chla_mean = chla_fn1, chlc_mean=chlc_fn1, TotalChl_mean = Totalchl_fn1,
         chla_SE = chla_fn2, chlc_SE = chlc_fn2, Totalchl_SE = Totalchl_fn2,
         chla_ug_mL_mean = chla_ug_mL_fn1, chlc_ug_mL_mean = chlc_ug_mL_fn1,
         Totalchl_ug_mL_mean = Totalchl_ug_mL_fn1, 
         chla_ug_mL_SE = chla_ug_mL_fn2, chlc_ug_mL_SE = chlc_ug_mL_fn2,
         Totalchl_ug_mL_SE = Totalchl_ug_mL_fn2) %>%
  select(Sample.Name,chla_mean,chla_SE,chlc_mean,chlc_SE,TotalChl_mean,Totalchl_SE,
         chla_ug_mL_mean,chla_ug_mL_SE, chlc_ug_mL_mean,chlc_ug_mL_SE,Totalchl_ug_mL_mean,Totalchl_ug_mL_SE ) # put in a better order
  
### Export the data
## all the info
write_csv(ChlData_ave, here("Data", foldername, paste0(platename,"_summary.csv"))) # print the summary data
write_csv(ChlData_raw, here("Data", foldername, paste0(platename,"_raw.csv"))) # print the raw data

