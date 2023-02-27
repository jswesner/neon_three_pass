library(neonUtilities)
library(tidyverse)
library(janitor)
library(lubridate)
library(tidybayes)
library(brms)
library(neonstore)
library(neonDivData)

# directory
Sys.setenv(NEONSTORE_HOME = paste(getwd(), 
                                  "/data",
                                  sep=""))

# download data (takes ~15 minutes) --------------------------------
#stream sites
streamsites=c("HOPB", "LEWI", "POSE", "CUPE",
              "GUIL", "KING", "MCDI", "LECO",
              "WALK", "MAYF", "ARIK", "BLUE",
              "PRIN", "BLDE", "COMO", "WLOU", 
              "SYCA", "REDB", "MART", "MCRA",
              "BIGC", "TECR", "OKSR", "CARI")

# neon_download(product="DP1.20107.001",
#               start_date=NA,
#               end_date=NA,
#               type="basic",
#               site= NA)
# 
# neon_download(product="DP1.20190.001",
#               start_date=NA,
#               end_date=NA,
#               table = "rea_widthFieldData",
#               type="basic",
#               site= streamsites)
# 
# # stack data
# fish_stacked = stackFromStore(filepaths=neon_dir(),
#                       dpID="DP1.20107.001",
#                       package="basic",
#                       site = streamsites)
# 
# # this behaves oddly. I had to manually add the variables file from a previous download to get it to 
# # work.
# stream_widths_stacked = stackFromStore(filepaths=neon_dir(),
#                                dpID="DP1.20190.001",
#                                package="basic",
#                                site = streamsites)
# 
# saveRDS(fish_stacked, file = "data/raw_data/fish/fish_stacked.rds")
# saveRDS(stream_widths_stacked, file = "data/raw_data/fish/stream_widths_stacked.rds")


# get reach lengths and widths --------------------------------------------
stream_widths_stacked = readRDS("data/raw_data/fish/stream_widths_stacked.rds")
mean_wetted_width = stream_widths_stacked$rea_widthFieldData %>%
  clean_names() %>% 
  select(site_id, collect_date, wetted_width) %>% 
  mutate(year = year(collect_date),
         month = month(collect_date),
         year_month = paste(year,month, sep = "_")) %>% 
  group_by(site_id) %>% 
  summarize(mean_wetted_width_m = mean(wetted_width, na.rm = T),
            sd_wetted_width_m = sd(wetted_width, na.rm = T))


# TOTAL POPULATION wrangle data ------------------
fish <- readRDS("data/raw_data/fish/fish_stacked.rds")

fish_bulk = fish$fsh_bulkCount %>% 
  select(eventID, taxonID, bulkFishCount) %>% 
  separate(eventID, into = c("site_id", "date", "reach", "pass", "method")) %>% 
  mutate(reach_id = paste(site_id, date, reach, sep = ".")) %>% 
  rename(n = bulkFishCount)

fish_measures = fish$fsh_perFish %>% 
  select(eventID, taxonID) %>% 
  separate(eventID, into = c("site_id", "date", "reach", "pass", "method")) %>% 
  mutate(reach_id = paste(site_id, date, reach, sep = ".")) %>% 
  group_by(site_id, date, reach, reach_id, pass) %>% 
  tally()

fish_reach_length = fish$fsh_fieldData %>% 
  clean_names() %>% 
  distinct(reach_id, measured_reach_length) %>% 
  group_by(reach_id) %>% 
  add_tally() %>% 
  filter(n == 1)   # filters duplicate reach lengths

three_pass_data = bind_rows(fish_bulk, fish_measures) %>% group_by(reach_id, pass) %>% 
  summarize(total_fish = sum(n, na.rm = T)) %>% 
  separate(reach_id, into = c("site_id", "date", "reach"), remove = F) %>% 
  mutate(date = ymd(date),
         month = month(date),
         year = year(date),
         year_month = paste(year, month, sep = "_")) %>% 
  left_join(fish_reach_length) %>% 
  left_join(mean_wetted_width)

write_csv(three_pass_data, file = "data/raw_data/fish/three_pass_data.csv")


# SPECIES POPULATION wrangle data ------------------
fish <- readRDS("data/raw_data/fish/fish_stacked.rds")

fish_bulk_species = fish$fsh_bulkCount %>% 
  select(eventID, taxonID, bulkFishCount) %>% 
  separate(eventID, into = c("site_id", "date", "reach", "pass", "method")) %>% 
  mutate(reach_id = paste(site_id, date, reach, sep = ".")) %>% 
  rename(n = bulkFishCount)

fish_measures_species = fish$fsh_perFish %>% 
  select(eventID, taxonID) %>% 
  separate(eventID, into = c("site_id", "date", "reach", "pass", "method")) %>% 
  mutate(reach_id = paste(site_id, date, reach, sep = ".")) %>% 
  group_by(site_id, date, reach, reach_id, pass, taxonID) %>% 
  tally()

three_pass_data_species = bind_rows(fish_bulk_species, fish_measures_species) %>% 
  group_by(reach_id, pass, taxonID) %>% 
  summarize(total_fish = sum(n, na.rm = T)) %>% 
  separate(reach_id, into = c("site_id", "date", "reach"), remove = F) %>% 
  mutate(date = ymd(date),
         month = month(date),
         year = year(date),
         year_month = paste(year, month, sep = "_")) %>% 
  left_join(fish_reach_length) %>% 
  left_join(mean_wetted_width)

write_csv(three_pass_data_species, file = "data/raw_data/fish/three_pass_data_species.csv")



