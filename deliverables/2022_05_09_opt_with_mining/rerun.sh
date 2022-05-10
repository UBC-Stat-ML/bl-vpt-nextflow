#!/usr/bin/env Rscript
require("ggplot2")
require("dplyr") 

paths <- read.csv("aggregated/optimizationPath.csv.gz")
paths <- filter(paths, budget <= 16000)
ggplot(paths, aes(x = budget, y = value, color = factor(random))) +
  facet_grid(objective + optimizer + name ~ factor(stepScale), labeller = label_both) +
  scale_x_log10() +
  xlab("Budget (number of exploration steps)") + 
  ylab("Parameter") + 
  ylim(-100, 350) +
  geom_line(alpha = 0.5)  + 
  theme_bw()
ggsave(paste0("optimizationPaths.pdf"), width = 17, height = 30)

optmonitor <- read.csv("aggregated/optimizationMonitoring.csv.gz")
optmonitor <- filter(optmonitor, name == "Rejection")
optmonitor <- filter(optmonitor, budget <= 16000) # when hitting NaN, budget can be 2x larger
ggplot(optmonitor, aes(x = budget, y = value, color = factor(random))) +
  facet_grid(objective + optimizer ~ factor(stepScale), labeller = label_both) +
  scale_x_log10() +
  xlab("Budget (number of exploration steps)") + 
  ylab("Global Communication Barrier (GCB)") + 
  geom_line(alpha = 0.5)  + 
  theme_bw()
ggsave(paste0("optimizationMonitoring.pdf"), width = 17, height = 7)

optmonitor %>% 
  filter(is.finite(value)) %>% 
  group_by(budget, objective, optimizer, stepScale) %>%
  summarise(mean_GCB = mean(value)) %>%
  ggplot(aes(x = budget, y = mean_GCB, colour = optimizer)) +
    facet_grid(objective + optimizer ~ factor(stepScale), labeller = label_both) +
    scale_x_log10() +
    xlab("Budget (number of exploration steps)") + 
    ylab("GCB (averaged over 10 restarts, ignoring failures)") + 
    geom_line(alpha = 1) + 
    theme_bw()
ggsave(paste0("optimizationMonitoring-mean.pdf"), width = 17, height = 7)

optmonitor$isFinite <- is.finite(optmonitor$value)
optmonitor %>% 
  group_by(budget, objective, optimizer, stepScale) %>%
  summarise(mean_is_finite = sum(isFinite)/10) %>%
  ggplot(aes(x = budget, y = mean_is_finite, colour = optimizer)) +
    facet_grid(objective + optimizer ~ factor(stepScale), labeller = label_both) +
    scale_x_log10() +
    ylim(0.0, 1.0) + 
    xlab("Budget (number of exploration steps)") + 
    ylab("Fraction of runs with finite objective") + 
    geom_line(alpha = 1) + 
    theme_bw()
ggsave(paste0("optimizationMonitoring-isFinite.pdf"), width = 17, height = 7)
