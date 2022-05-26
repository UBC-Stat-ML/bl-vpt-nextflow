#!/usr/bin/env Rscript
require("ggplot2")
require("dplyr")
require("ggridges")

ks_distances <- read.csv("ks_distances.csv")
ks_distances$quality <- ifelse(ks_distances$ks_stat > 0.1, "poor", "good")

restarts <- read.csv("aggregated/actualTemperedRestarts.csv.gz")
restarts <- restarts %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
restarts %>%
  filter(algorithm != "Reference") %>%
  group_by(quality, model, algorithm, seed) %>%
  summarize(total_count = sum(count) ) %>%
  ggplot(aes(x = algorithm, y = total_count, color = quality)) +
    facet_grid(. ~ model, scales="free_x") +
    geom_boxplot() +
    coord_flip() +
    ylab("total tempered restarts") + 
    scale_y_continuous(expand = expansion(mult = 0.05), limits = c(0, NA)) +
    scale_fill_manual(values = c(  "good" = "blue", "poor" = "red")) + scale_colour_manual(values = c(  "good" = "blue", "poor" = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("actualTemperedRestarts-box.pdf", width = 10, height = 2)

ess <- read.csv("aggregated/allEss.csv.gz")
ess <- ess %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
statEss <- ess %>%
  filter(algorithm != "Reference") %>%
  filter(variable != "sigma" | model != "mrna-no-transf") %>% # hack.. in future avoid duplicate!
  filter(variable != "tau" | model != "sparse-car") %>% 
  filter(variable %in% c("alpha")) %>%
  group_by(quality, model, algorithm, variable, seed) %>%
  summarize(total_ess = sum(value))
write.csv(statEss, "statEss.csv")

statEss %>%
  ggplot(aes(x = algorithm, y = total_ess, color = quality)) +
    facet_grid(. ~ model, scales="free_x") +
    coord_flip() +
    geom_boxplot() +
    ylab("ESS") + 
    scale_y_continuous(expand = expansion(mult = 0.05), limits = c(0, NA)) +
    scale_fill_manual(values = c(  "good" = "blue", "poor" = "red")) + scale_colour_manual(values = c(  "good" = "blue", "poor" = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("ess-box.pdf", width = 10, height = 2)

global <- read.csv("aggregated/globalLambda.csv.gz")
global <- global %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
global %>%
  filter(fixedRefChain == "false") %>%
  ggplot(aes(x = round, y = value, color = quality, group = factor(seed))) +
    facet_grid(model ~ algorithm, scales="free_y") +
    scale_fill_manual(values = c(  "good" = "blue", "poor" = "red")) + scale_colour_manual(values = c(  "good" = "blue", "poor" = "red")) +
    ylab("Lambda") + 
    geom_line() +
    theme_bw()
ggsave("globalLambda.pdf", width = 17, height = 8)  

lambda <- read.csv("aggregated/lambdaInstantaneous.csv.gz")
lambda <- lambda %>% inner_join(ks_distances, by = c("algorithm", "model", "seed"))
lambda %>% 
  filter(isAdapt == "false") %>%
  filter(fixedRefChain == "false") %>%
  ggplot(aes(x = beta, y = value, color = quality, group = factor(seed))) +
    facet_grid(model ~ algorithm, scales="free_y") +
    scale_fill_manual(values = c(  "good" = "blue", "poor" = "red")) + scale_colour_manual(values = c(  "good" = "blue", "poor" = "red")) +
    geom_line(alpha = 0.5) +
    xlab("beta") + 
    ylab("intensity") + 
    theme_bw()
ggsave("lambda.pdf", width = 17, height = 8)
