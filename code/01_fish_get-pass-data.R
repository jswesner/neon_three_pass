# download NEON data for three-pass fish eclectroshocking, wrangle it. The resulting data 
# is used to estimate population size with three-pass depletion models. 
# Modeling occurs in a different script.

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
# # add a variables file. Just need one overall file. Any site will do. I've chosen ARIK b/c it's first in alphabet.
# neon_download(product="DP1.20190.001",
              # start_date="2021-01-01",
              # end_date="2022-01-01",
              # table = "variables",
              # type="basic",
              # site= "ARIK")
# # stack data
# fish_stacked = stackFromStore(filepaths=neon_dir(),
#                       dpID="DP1.20107.001",
#                       package="basic",
#                       site = streamsites)
# 
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

# 1) get reachids, eventids
reach_event = fish$fsh_perPass %>% 
  select(reachID, eventID) %>% 
  distinct()  %>% # removes duplicates. JSW confirmed that these were true duplicates on 2023-03-01
  clean_names()

# 2) get info on whether reaches are fixed or random. Only fixed reaches have 3 pass removal. Random reaches are all single pass
fixed_random_reach = fish$fsh_fieldData %>% 
  distinct(reachID, fixedRandomReach) %>% 
  filter(reachID != "WALK.20170316.07" | fixedRandomReach != "fixed") %>% # fix this type. Confirmed by email with NEON on 2023-03-01
  clean_names()

# 3) get reach lengths
fish_reach_length = fish$fsh_fieldData %>% 
  clean_names() %>% 
  distinct(reach_id, measured_reach_length) %>% 
  group_by(reach_id) %>% 
  add_tally() %>% 
  filter(n == 1)   # filters duplicate reach lengths

# 4) the first 50 fish are measured for length. They are here in fsh_perFish
fish_measures = fish$fsh_perFish %>% 
  select(eventID, taxonID) %>% 
  separate(eventID, into = c("site_id", "date", "reach", "pass", "method"), remove = F) %>%
  clean_names() %>% 
  left_join(reach_event) %>% 
  group_by(site_id, date, reach, reach_id, pass, event_id) %>% 
  tally()

# 5) all fish after the first 50 are bulk counted. Those values are here. 
fish_bulk = fish$fsh_bulkCount %>% 
  select(eventID, taxonID, bulkFishCount)  %>% 
  separate(eventID, into = c("site_id", "date", "reach", "pass", "method"), remove = F) %>%
  clean_names() %>% 
  left_join(reach_event)  %>% 
  rename(n = bulk_fish_count)

# 6) combine the 1st 50 fish with the bulk counts to get a total number of fish per pass.
#    Then add identifying information

three_pass_data = bind_rows(fish_bulk, fish_measures) %>% 
  group_by(reach_id, pass, event_id) %>% 
  summarize(total_fish = sum(n, na.rm = T)) %>% 
  separate(reach_id, into = c("site_id", "date", "reach"), remove = F) %>% 
  mutate(date = ymd(date),
         month = month(date),
         year = year(date),
         year_month = paste(year, month, sep = "_")) %>% 
  ungroup %>% 
  left_join(fish_reach_length) %>% 
  left_join(mean_wetted_width) %>% 
  left_join(fixed_random_reach)

# 7) get reachids, eventids, and targets (i.e., passes that returned zero fish will have "N")
true_zeros = fish$fsh_perPass %>% 
  select(reachID, targetTaxaPresent) %>% 
  distinct()  %>% # removes duplicates. JSW confirmed that these were true duplicates on 2023-03-01
  clean_names() %>%
  filter(target_taxa_present == "N")

# 8) make wide format. Replace 0's using information in fsh_perPass$target_taxa_present
three_pass_data_wide_total_fish = three_pass_data %>%   # restrict to fixed reaches only
  filter(fixed_random_reach == "fixed") %>%  
  select(-event_id) %>% 
  pivot_wider(names_from = pass, values_from = total_fish) %>% 
  ungroup() %>%
  replace_na(list(`1` = 0,            # replace NA with zeros (assumes zero fish if there were no values entered)
                  `2` = 0,
                  `3` = 0)) %>% 
  left_join(true_zeros)  %>%
  mutate(`1` = case_when(`1` == 0 & is.na(target_taxa_present) ~ 1e9,   # create silly number to filter out false zeros
                         TRUE ~ `1`),
         `2` = case_when(`2` == 0 & is.na(target_taxa_present) ~ 1e9,
                         TRUE ~ `2`),
         `3` = case_when(`3` == 0 & is.na(target_taxa_present) ~ 1e9,
                         TRUE ~ `3`)) %>% 
  filter(`1` < 1e9) %>% # filter out false zeros
  filter(`2` < 1e9) %>% 
  filter(`3` < 1e9) %>% 
  mutate(site_int = as.factor(row_number()), #sample identifier
         increased = case_when(`3` > `1` ~ "no depletion", 
                               TRUE ~ "depletion")) 

write_csv(three_pass_data_wide_total_fish, file = "data/raw_data/fish/three_pass_data_wide_total_fish.csv")


# SPECIES POPULATION wrangle data ------------------
fish_measures_species = fish$fsh_perFish %>% 
  select(eventID, taxonID) %>% 
  separate(eventID, into = c("site_id", "date", "reach", "pass", "method"), remove = F) %>%
  clean_names() %>% 
  left_join(reach_event) %>% 
  group_by(site_id, date, reach, reach_id, pass, event_id, taxon_id) %>% 
  tally()

three_pass_data_species = bind_rows(fish_bulk, fish_measures_species) %>% 
  group_by(reach_id, pass, event_id, taxon_id) %>% 
  summarize(total_fish = sum(n, na.rm = T)) %>% 
  separate(reach_id, into = c("site_id", "date", "reach"), remove = F) %>% 
  mutate(date = ymd(date),
         month = month(date),
         year = year(date),
         year_month = paste(year, month, sep = "_")) %>% 
  ungroup %>% 
  left_join(fish_reach_length) %>% 
  left_join(mean_wetted_width) %>% 
  left_join(fixed_random_reach)

saveRDS(three_pass_data_species, file = "data/raw_data/fish/three_pass_data_species.rds")

three_pass_data_wide_species = three_pass_data_species %>%   # restrict to fixed reaches only
  filter(fixed_random_reach == "fixed") %>% 
  select(-event_id) %>% 
  pivot_wider(names_from = pass, values_from = total_fish) %>% 
  ungroup() %>%
  replace_na(list(`1` = 0,            # replace NA with zeros (assumes zero fish if there were no values entered)
                  `2` = 0,
                  `3` = 0)) %>% 
  left_join(true_zeros)  %>%
  mutate(`1` = case_when(`1` == 0 & is.na(target_taxa_present) ~ 1e9,   # create silly number to filter out false zeros
                         TRUE ~ `1`),
         `2` = case_when(`2` == 0 & is.na(target_taxa_present) ~ 1e9,
                         TRUE ~ `2`),
         `3` = case_when(`3` == 0 & is.na(target_taxa_present) ~ 1e9,
                         TRUE ~ `3`)) %>% 
  filter(`1` < 1e9) %>% # filter out false zeros
  filter(`2` < 1e9) %>% 
  filter(`3` < 1e9) %>% 
  mutate(site_int = as.factor(row_number()),  #sample identifier
         increased = case_when(`3` > `1` ~ "no depletion", 
                               TRUE ~ "depletion")) 

write_csv(three_pass_data_wide_species, file = "data/raw_data/fish/three_pass_data_wide_species.csv")
