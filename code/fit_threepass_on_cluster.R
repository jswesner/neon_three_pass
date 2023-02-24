library(ubms)

three_pass_frame = readRDS(file = "data/raw_data/fish/three_pass_frame.rds")

three_pass_model = stan_multinomPois(formula = ~site_int ~ site_int,
                                     data = three_pass_frame,
                                     chains = 4, iter = 1000)

saveRDS(three_pass_model, file = "models/three_pass_model.rds")

