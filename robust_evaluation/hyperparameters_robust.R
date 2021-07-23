library(data.table)
library(ggplot2)
robustness_tests <- lapply(list.files("../robustness_comparison/ROBUSTOut_Hyperparameters/", full.names = T), fread)
names(robustness_tests) <- tstrsplit(tstrsplit(list.files("../robustness_comparison/ROBUSTOut_Hyperparameters/"), "ROBUST_", keep=2)[[1]], ".out", keep=1)[[1]]

robustness_tests <- rbindlist(robustness_tests, idcol="init_red")
library(tidyr)
robustness_tests <- robustness_tests %>% 
  separate(col="init_red", into=c("Initial Fraction", "Reduction Factor"), sep="_") %>%
  separate(col="seed_set", into=c("seed_set", "Threshold", "Nr of Trees"), sep="_")
robustness_tests[, Threshold := tstrsplit(Threshold, "key", keep=2)]
robustness_tests[, Threshold := as.factor(Threshold)]
robustness_tests[, `Initial Fraction` := as.factor(paste0("Init:", `Initial Fraction`))]
robustness_tests[, `Reduction Factor` := as.factor(paste0("Red:", `Reduction Factor`))]
robustness_tests[, `Nr of Trees` := factor(`Nr of Trees`, levels = c(5, 10, 15, 20))]

ggplot(robustness_tests[`Initial Fraction` == "Init:0.25" & `Reduction Factor` == "Red:0.9"])+
  geom_boxplot(aes(x = Threshold, y = `mean jaccard`, color = `Nr of Trees`))+
  labs(y="Mean Jaccard")+
  theme_bw()+
  theme(text=element_text(size = 13))
ggsave("trees_025_09.png", width = 10, height = 7)

ggplot(robustness_tests, aes(x=Threshold, y = `mean jaccard`))+
  geom_boxplot()+
  labs(y="Mean Jaccard")+
  facet_grid(`Initial Fraction`~ `Reduction Factor`)+
  theme_bw()+
  theme(text=element_text(size = 13))
ggsave("init_vs_red.png", width = 14, height = 9)
