#!/usr/bin/env Rscript
require("ggplot2")
require("dplyr")

paths <- read.csv("aggregated/optimizationPath.csv.gz")
ggplot(paths, aes(x = budget, y = value, color = factor(random), linetype = useFixedRef)) +
  facet_grid(name ~ useFixedRef, scales="free_y") +
  scale_x_log10() +
  xlab("Budget (number of exploration steps)") + 
  ylab("Parameter") + 
  geom_line(alpha = 0.5)  + 
  theme_bw()
ggsave(paste0("optimizationPaths.pdf"), width = 9, height = 10)

optmonitor <- read.csv("aggregated/optimizationMonitoring.csv.gz")
optmonitor <- filter(optmonitor, name == "Rejection")
ggplot(optmonitor, aes(x = budget, y = value, color = factor(random), linetype = useFixedRef)) +
  scale_x_log10() +
  facet_grid(. ~ useFixedRef) +
  xlab("Budget (number of exploration steps)") + 
  ylab("Global Communication Barrier (GCB)") + 
  geom_line(alpha = 0.5)  + 
  theme_bw()
ggsave(paste0("optimizationMonitoring.pdf"), width = 9, height = 7)

optmonitor %>% 
  filter(is.finite(value)) %>% 
  group_by(budget, objective, optimizer, stepScale, useFixedRef, nChains) %>%
  summarise(mean_GCB = mean(value)) %>%
  ggplot(aes(x = budget, y = mean_GCB, colour = optimizer, linetype = useFixedRef)) +
    scale_x_log10() +
    xlab("Budget (number of exploration steps)") + 
    ylab("GCB (averaged over 10 restarts, ignoring failures)") + 
    geom_line(alpha = 1) + 
    theme_bw()
ggsave(paste0("optimizationMonitoring-mean.pdf"), width = 9, height = 7)
