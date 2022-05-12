#!/usr/bin/env Rscript
require("ggplot2")
require("dplyr")

betas <- read.csv("statistic.csv.gz")

betas %>% 
  filter(algorithm == 'V--T*--F') %>%
  filter(sample > 250) %>% 
  ggplot(aes(x = value, colour = chain, group = chain)) + 
    geom_density() +
    theme_minimal() + 
    xlab("parameter") + 
    ylab("density") + 
    theme(legend.position="none") +
    xlim(-0.4, 1.2) + 
    ylim(0, 7) +
    scale_colour_gradient(low = "grey30", high = "skyblue1")   
ggsave("variational-paths.pdf", width = 2, height = 2)

betas %>% 
  filter(algorithm == 'Reference') %>%
  filter(sample > 250) %>% 
  ggplot(aes(x = value, colour = chain, group = chain)) + 
    geom_density() +
    theme_minimal() +
    xlab("parameter") + 
    ylab("density") + 
    theme(legend.position="none") +
    xlim(-0.4, 1.2) +
    ylim(0, 7) +  
    scale_colour_gradient(low = "grey30", high = "orange") 
ggsave("fixed-ref-paths.pdf", width = 2, height = 2)
