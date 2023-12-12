library(tidyverse)
library(ubms)
library(brms)
library(janitor)
library(tidybayes)

three_pass_data_wide = read_csv("data/raw_data/fish/three_pass_data_wide_total_fish.csv")

three_pass_model = readRDS(file = "models/three_pass_model.rds")

three_pass_model

sample_1  = as_draws_df(three_pass_model@stanfit) %>% 
  clean_names() %>%
  select(contains("det_intercept")) %>% 
  mutate(name = "sample_1") %>% 
  rename(value = beta_det_intercept)

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

sample_capture_probs = capture_prob_posts %>% 
  mutate(site_int = parse_number(name)) %>% #get original group names
  left_join(reach_info) %>% 
  # mutate(prob = inv_logit_scaled(value)) %>%
  # filter(reach_id != "KING.20201021.04")  %>% # this sample is probably not reliable () 
  group_by(site_id, date, reach, year, site_int) %>% # group and summarize
  median_qi(median_prob = inv_logit_scaled(value)) %>% # summarize on the probability scale (via link function)
  select(-.width, -.point, -.interval) %>% 
  # group_by(site_id) %>% 
  mutate(max = max(median_prob),
         month = str_sub(date, 5,6))

sample_capture_plot = sample_capture_probs %>% 
  ggplot(aes(x = reorder(site_int, -median_prob), y = median_prob, color = site_id)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper),
                  size = 0.2) + 
  facet_wrap(~site_id, scales = "free_x") + 
  theme_default() + 
  guides(color = "none") + 
  theme(axis.text.x = element_blank()) + 
  labs(x = "Sample (ranked by capture probability",
       y = "Capture Probability on first pass",
       title = "Modeled capture probabilities of 466 NEON 3-pass removal samples") + 
  geom_hline(yintercept = 0.5)

ggview::ggview(sample_capture_plot, width = 7.5, height = 9)
ggsave(sample_capture_plot, file = "plots/sample_capture_plot.jpg",
       width = 7.5, height = 9)

sample_capture_table = sample_capture_probs %>% 
  select(site_id, date, reach, year, site_int, median_prob, .lower, .upper) %>% 
  mutate(note = "median and lower/upper 95% credible intervals of capture probabilities on the first pass")

write_csv(sample_capture_table, file = "tables/sample_capture_table.csv")


# total probability over 3 passes -----------------------------------------
# detection probability per pass
per_pass_det = getP(three_pass_model) 
colnames(per_pass_det) = c("first", "second", "third")

per_pass_posts = as_tibble(per_pass_det) %>% 
  mutate(site_int = row_number()) %>% 
  pivot_longer(cols = -site_int) %>% 
  separate(name, into = c("pass", ".draw"))


per_pass_summary = per_pass_posts %>% 
  left_join(reach_info) %>% 
  group_by(site_int, site_id, date, reach, year, pass) %>% 
  median_qi(value)

per_pass_summary %>% 
  ggplot(aes(x = pass, y= value)) + 
  geom_point() +
  geom_line(aes(group = site_int)) + 
  facet_wrap(~site_id)


total_prob_summary = per_pass_posts %>% 
  left_join(reach_info) %>% 
  group_by(site_int, site_id, date, reach, year, .draw) %>% 
  summarize(value = sum(value)) %>%
  group_by(site_int, site_id, year, reach, date) %>% 
  median_qi(value)

total_prob_plot = total_prob_summary %>% 
  ggplot(aes(x = reorder(site_int, -value), y = value, color = site_id)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper)) + 
  facet_wrap(~site_id, scales = "free_x") +
  theme_default() + 
  labs(x = "Sample (ranked by capture probability)",
       y = "Capture Probability (sum across 3 passes)",
       title = "Modeled capture probabilities of 466 NEON 3-pass removal samples") + 
  geom_hline(yintercept = 0.5) + 
  theme(axis.text.x = element_blank()) +
  guides(color = "none") 

total_prob_table = total_prob_summary %>% 
  select(site_int, site_id, year, reach, date, value, .lower, .upper)

ggsave(total_prob_plot, file = "plots/total_prob_plot.jpg",
       width = 7.5, height = 9)

write_csv(total_prob_table, file = "tables/total_prob_table.csv")


# run for sites ------------------------------------------------
three_pass_frame_site <- unmarkedFrameMPois(three_pass_matrix,
                                       siteCovs=as.data.frame(three_pass_data_wide %>% select(site_int, site_id)),
                                       type = "removal")

# fit model
# this estimates population size and capture efficiency for each reach_id (called site_int here)
# three_pass_frame = readRDS(file = "data/raw_data/fish/three_pass_frame.rds")
# 
three_pass_model_sites = stan_multinomPois(formula = ~site_id + (1|site_int) ~ site_id + (1|site_int),
                                     data = three_pass_frame_site,
                                     chains = 1, iter = 500)

saveRDS(three_pass_model_sites, file = "models/three_pass_model_sites.rds")
