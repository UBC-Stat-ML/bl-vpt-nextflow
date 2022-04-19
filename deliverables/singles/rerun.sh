#!/usr/bin/env Rscript
require("ggplot2")
require("dplyr") 

lambda <- read.csv("aggregated/lambdaInstantaneous.csv.gz")
lambda %>% 
  filter(isAdapt == "false") %>%
  mutate(beta2 = ifelse(fixedRefChain == "true", -beta, beta)) %>%
  ggplot(aes(x = beta2, y = value, color = fixedRefChain)) +
    facet_grid(model ~ algo, scales="free_y") +
    geom_line() +
    xlab("Beta") + 
    ylab("Intensity") + 
    theme_bw()
ggsave(paste0("lambda.pdf"), width = 17, height = 8)

restarts <- read.csv("aggregated/actualTemperedRestarts.csv.gz")
restarts %>%
  ggplot(aes(x = round, y = count, color = fixedRefChain)) +
    facet_grid(model ~ algo, scales="free_y") +
    geom_line() +
    theme_bw()
ggsave(paste0("actualTemperedRestarts.pdf"), width = 17, height = 8)

global <- read.csv("aggregated/globalLambda.csv.gz")
global %>%
  ggplot(aes(x = round, y = value, color = fixedRefChain)) +
    facet_grid(model ~ algo, scales="free_y") +
    geom_line() +
    theme_bw()
ggsave(paste0("globalLambda.pdf"), width = 17, height = 8)
