library(tidyverse)
library(ubms)
library(brms)
library(janitor)

# load data
fish <- readRDS("data/raw_data/fish/fish_stacked.rds")

fixed_reaches = fish$fsh_fieldData %>% distinct(reachID, fixedRandomReach) %>% 
  clean_names() %>% filter(fixed_random_reach == "fixed") %>% select(reach_id) %>% pull()

three_pass_data = read_csv("data/raw_data/fish/three_pass_data.csv") %>%   # restrict to fixed reaches only
  filter(reach_id %in% fixed_reaches)

# fit multinomial poisson (three pass depletion model) -----------------------------------------------
# wrangle
three_pass_data_wide = three_pass_data %>% 
  pivot_wider(names_from = pass, values_from = total_fish) %>% 
  ungroup() %>% 
  replace_na(list(`1` = 0,            # replace NA with zeros (assumes zero fish if there were no values entered)
                  `2` = 0,
                  `3` = 0)) %>% 
  mutate(year = as.numeric(str_sub(reach_id, 6, 9))) %>% 
  filter(year >= 2016 & year <2022) %>% 
  mutate(site_int = as.factor(row_number())) 
  
# put passes in a matrix (for ubms)
three_pass_matrix = three_pass_data_wide %>% 
  select(`1`,`2`,`3`) %>% 
  as.matrix()

# assign covariates (for ubms)
three_pass_frame <- unmarkedFrameMPois(three_pass_matrix,
                                siteCovs=as.data.frame(three_pass_data_wide %>% select(site_int)),
                                type = "removal")
saveRDS(three_pass_frame, file = "data/raw_data/fish/three_pass_frame.rds")

# fit model
# this estimates population size and capture efficiency for each reach_id (called site_int here)
# three_pass_frame = readRDS(file = "data/raw_data/fish/three_pass_frame.rds")
# 
# three_pass_model = stan_multinomPois(formula = ~site_int ~ site_int,
#                                      data = three_pass_frame,
#                                      chains = 4, iter = 1000)
# 
# saveRDS(three_pass_model, file = "models/three_pass_model.rds")

# extract posteriors and summarize
# get intercept value
three_pass_model = readRDS(file = "models/three_pass_model.rds")

sample_1  = as_draws_df(three_pass_model@stanfit) %>% 
  clean_names() %>%
  select(contains("det_intercept")) %>% 
  mutate(name = "sample_1") %>% 
  rename(value = beta_det_intercept)

# calculate probs then bind the first sample (intercept above)
capture_prob_posts = as_draws_df(three_pass_model@stanfit) %>% 
  select(contains("_det")) %>% 
  pivot_longer(cols = !contains("ntercept")) %>% 
  clean_names() %>% 
  mutate(value = beta_det_intercept + value) %>% 
  bind_rows(sample_1) %>% 
  select(name, value)

# summarize site specific capture probabilities  
reach_info = three_pass_data_wide %>% ungroup %>% # add identifiers
  mutate(site_int = parse_number(as.character(site_int))) %>% 
  distinct(site_int, reach_id, year) %>% 
  separate(reach_id, into = c("site_id", "date", "reach"), remove = F)

site_capture_probs = capture_prob_posts %>% 
  mutate(site_int = parse_number(name)) %>% #get original group names
  left_join(reach_info) %>% 
  filter(reach_id != "KING.20201021.04")  %>% # this sample is probably not reliable () 
  group_by(site_id) %>% # group and summarize
  median_qi(median_prob = inv_logit_scaled(value)) %>% # summarize on the probability scale (via link function)
  select(-.width, -.point, -.interval) 


# fit poisson (single pass) -----------------------------------------------

stream_fish_firstpass = three_pass_data %>%
  separate(reach_id, into = c("site_id", "date", "reach"), remove = F) %>% 
  mutate(date = ymd(date),
         julian = julian(date)) %>% 
  filter(pass == 1) %>%
  left_join(site_capture_probs) %>% 
  group_by(reach_id, pass) 

fish_total_abundance_poisson = brm(total_fish ~ reach_id + offset(median_prob),
                             family = poisson(link = "log"),
                             data = stream_fish_firstpass,
                             prior = c(prior(normal(0, 1), class = "b"),
                                       prior(normal(4, 4), class = "Intercept")),
                             chains = 4, iter = 2000,
                             file = "models/fish_total-abundance-poisson.rds",
                             file_refit = "on_change")


# extract posteriors and summarize
sample_1_population  = as_draws_df(three_pass_model@stanfit) %>% 
  clean_names() %>%
  select(contains("state_intercept")) %>% 
  mutate(name = "sample_1") %>% 
  rename(value = beta_state_intercept)

three_pass_population = as_draws_df(three_pass_model@stanfit) %>% 
  select(contains("_state")) %>% 
  pivot_longer(cols = !contains("ntercept")) %>% 
  clean_names() %>% 
  mutate(value = beta_state_intercept + value) %>% 
  bind_rows(sample_1) %>% 
  select(name, value) %>% 
  mutate(site_int = as.factor(parse_number(name))) %>% #get original group names
  left_join(three_pass_data_wide %>% ungroup %>% 
              distinct(site_int, reach_id)) %>% 
  group_by(reach_id) %>% # group and summarize
  median_qi(pop_threepass = exp(value)) %>% # summarize on the probability scale (via link function)
  select(-.width, -.point, -.interval) %>% 
  rename(.lower_threepass = .lower,
         .upper_threepass = .upper) 

single_pass_population = fish_total_abundance_poisson$data %>% 
  distinct(reach_id, median_prob) %>% 
  add_epred_draws(fish_total_abundance_poisson) %>% 
  mutate(.epred = .epred/median_prob) %>% 
  group_by(reach_id) %>% 
  median_qi(.epred) %>% 
  rename(pop_singlepass = .epred,
         .lower_singlepass = .lower,
         .upper_singlepass = .upper) %>% 
  select(-.width, -.point, -.interval) %>% 
  filter(reach_id %in% fixed_reaches)



# plot --------------------------------------------------------------------
# combine all estimates
all_population_estimates = fish_total_abundance_poisson$data %>% 
  left_join(single_pass_population) %>% 
  left_join(three_pass_population) 

# saveRDS(all_population_estimates, file = "data/raw_data/fish/all_population_estimates.rds")

all_population_estimates = readRDS(file = "data/raw_data/fish/all_population_estimates.rds") %>% 
  separate(reach_id, into = c("site_id", "date", "reach"), remove = "F")

# plot
all_population_estimates %>%
  as_tibble() %>% 
  ggplot(aes(x = pop_singlepass, y = pop_threepass)) + 
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = .lower_singlepass, xmax = .upper_singlepass),
                 alpha = 0.2) + 
  geom_errorbar(aes(ymin = .lower_threepass, ymax = .upper_threepass),
                alpha = 0.2) +
  geom_abline() + 
  scale_x_log10(limits = c(1, 45000)) +
  scale_y_log10(limits = c(1, 45000)) +
  theme_default() +
  NULL
