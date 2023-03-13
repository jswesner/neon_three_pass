library(tidyverse)
library(ubms)
library(brms)
library(janitor)
library(tidybayes)

# load data
three_pass_data_wide = read_csv("data/raw_data/fish/three_pass_data_wide_total_fish.csv")

# fit multinomial poisson (three pass depletion model) -----------------------------------------------

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
#                                      chains = 1, iter = 500)
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
  # filter(reach_id != "KING.20201021.04")  %>% # this sample is probably not reliable () 
  group_by(site_id) %>% # group and summarize
  median_qi(median_prob = inv_logit_scaled(value)) %>% # summarize on the probability scale (via link function)
  select(-.width, -.point, -.interval) 


# fit poisson (single pass) -----------------------------------------------

stream_fish_firstpass = three_pass_data_wide %>%
  mutate(julian = julian(date)) %>% 
  mutate(total_fish = `1`) %>% 
  left_join(site_capture_probs) 


fish_total_abundance_poisson = readRDS("models/fish_total-abundance-poisson.rds")

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
  as_tibble() %>% 
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
              distinct(site_int, reach_id) %>% 
              mutate(site_int = as.factor(site_int))) %>% 
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
  select(-.width, -.point, -.interval) 

# plot --------------------------------------------------------------------
# combine all estimates
all_population_estimates = fish_total_abundance_poisson$data %>% 
  left_join(single_pass_population) %>% 
  right_join(three_pass_population) %>% 
  left_join(three_pass_data_wide %>% distinct(increased, reach_id))

saveRDS(all_population_estimates, file = "data/raw_data/fish/all_population_estimates.rds")

all_population_estimates = readRDS(file = "data/raw_data/fish/all_population_estimates.rds") %>% 
  separate(reach_id, into = c("site_id", "date", "reach"), remove = "F")

# plot
all_population_estimates %>%
  as_tibble() %>% 
  ggplot(aes(x = pop_singlepass, y = pop_threepass)) + 
  geom_point(size = 1, aes(color = increased)) +
  geom_errorbarh(aes(xmin = .lower_singlepass, xmax = .upper_singlepass),
                 alpha = 0.2) + 
  geom_errorbar(aes(ymin = .lower_threepass, ymax = .upper_threepass),
                alpha = 0.2) +
  geom_abline() + 
  # geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() +
  theme_default() +
  labs(color = "") +
  NULL


# plot
all_population_estimates %>%
  as_tibble() %>% 
  filter(increased == "depletion") %>% 
  ggplot(aes(x = pop_singlepass, y = pop_threepass)) + 
  geom_point(size = 1, aes(color = increased)) +
  # geom_errorbarh(aes(xmin = .lower_singlepass, xmax = .upper_singlepass),
                 # alpha = 0.2) + 
  # geom_errorbar(aes(ymin = .lower_threepass, ymax = .upper_threepass),
                # alpha = 0.2) +
  geom_abline() + 
  # geom_smooth(method = "lm") +
  scale_x_log10() +
  scale_y_log10() +
  theme_default() +
  labs(color = "") +
  NULL
