#!/usr/bin/env Rscript
require("ggplot2")
require("dplyr")
require("ggridges")
require("forcats")

ks_distances <- read.csv("ks_distances.csv")
ks_distances$quality <- ifelse(ks_distances$ks_stat > 0.1, "poor", "good")

restarts <- read.csv("aggregated/actualTemperedRestarts.csv.gz")
restarts <- restarts %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
restarts %>%
  filter(algorithm != "Reference") %>%
  filter(algorithm != "V--T* F--T*") %>%
  mutate(algorithm = fct_recode(algorithm, "Stabilized" = "V--T*--F")) %>%
  mutate(algorithm = fct_recode(algorithm, "Standard PT" = "F--T F--T")) %>%
  mutate(algorithm = fct_recode(algorithm, "Basic vari." = "V--T* F--T")) %>%
  group_by(quality, model, algorithm, seed) %>%
  summarize(total_count = sum(count) ) %>%
  ggplot(aes(x = algorithm, y = total_count, color = quality)) +
    facet_grid(. ~ model, scales="free_x") +
    geom_boxplot(show.legend = FALSE) +
    coord_flip() +
    ylab("total tempered restarts") + 
    scale_y_continuous(expand = expansion(mult = 0.05), limits = c(0, NA)) +
    scale_fill_manual(values = c(  "good" = "blue", "poor" = "red")) + scale_colour_manual(values = c(  "good" = "blue", "poor" = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("actualTemperedRestarts-box-subset.pdf", width = 7, height = 2)

