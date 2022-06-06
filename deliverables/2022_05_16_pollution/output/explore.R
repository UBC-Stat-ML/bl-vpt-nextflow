require("ggplot2")
require("dplyr")

read.csv("globalLambda.csv") %>%
  ggplot(aes(x = round, y = value)) + 
  geom_line() + 
  theme_minimal()
ggsave("globalLambda.pdf")


read.csv("optimizationPath.csv") %>%
  ggplot(aes(x = iter, y = value)) + 
  facet_grid(name ~ ., scales = "free_y") +
  geom_line() + 
  theme_minimal()
ggsave("params.pdf", height = 1000, limitsize = FALSE)
