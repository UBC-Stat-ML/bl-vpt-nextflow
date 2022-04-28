#!/usr/bin/env Rscript
require("ggplot2")
require("dplyr")

timing <- read.csv("aggregated/roundTimings.csv.gz") %>% rename(time = value)


global <- read.csv("aggregated/globalLambda.csv.gz")
global <- global %>% inner_join(timing, by = c("round", "model", "isVariational", "seed", "size"))
global %>% 
  filter(isAdapt == "false") %>%
  filter(model != "coll-rockets" | size < 369) %>%
  group_by(size, isVariational, model, round) %>%
  summarize(
    mean_gcb = mean(value),
    se_gcb = sd(value)/sqrt(n())) %>%
  ggplot(aes(x = size, y = mean_gcb, colour = factor(isVariational))) +
    geom_line() + 
    geom_errorbar(aes(ymin=mean_gcb-se_gcb, ymax=mean_gcb+se_gcb), width=.1) +
    scale_x_log10() +
    scale_y_log10() + 
    facet_grid(. ~ model) +
    theme_bw()

ggsave("scaling.pdf")
