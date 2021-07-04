library(data.table)
library(ggplot2)

all_results_rmust <- fread("../robustness_comparison/RMuSTOut/all_results_NeDRex.csv")
all_results_rmust <- tidyr::separate(all_results_rmust, "seed_set", c("seed_set", "threshold"), sep = "_")
ggplot(all_results_rmust, aes(y = `mean jaccard`, x = threshold))+
  geom_boxplot()+
  theme_bw()

all_results_diamond <- fread("../robustness_comparison/DiamondOut/all_results_Diamond.csv", keepLeadingZeros = T)
colnames(all_results_diamond) <- c("seed_set", "DIAMOnD")
ggplot(all_results_diamond, aes(x = DIAMOnD))+
  geom_histogram()+
  theme_bw()

all_results_domino <- fread("../robustness_comparison/DominoOut/all_results_Diamond.csv", keepLeadingZeros = T)
colnames(all_results_domino) <- c("seed_set", "DOMINO")
ggplot(all_results_domino, aes(x = DOMINO))+
  geom_histogram()+
  theme_bw()

all_results_robust <- fread("../robustness_comparison/RobustOut/all_results_biosteiner.csv")
all_results_robust <- tidyr::separate(all_results_robust, "seed_set", c("seed_set", "threshold"), sep = "_")
all_results_robust <- tidyr::separate(all_results_robust, "threshold", c("thr", "threshold"), sep = "r")
all_results_robust[, thr := NULL]

ggplot(all_results_robust, aes(x = `mean jaccard`, fill = threshold))+
  geom_histogram(position = "dodge")+
  theme_bw()

list_files <- lapply(list.files("../robustness_comparison/RobustOut/", full.names = T, pattern="*.txt"), function(x){
  fread(
    cmd = paste('head -n 3', x),
    sep= ":", header=FALSE
  )
})
library(tidyr)
names(list_files) <- tstrsplit(list.files("../robustness_comparison/RobustOut/", pattern="*.txt"), "_Bio", fixed=T, keep=1)[[1]]
robust_info <- rbindlist(list_files, idcol = "filename") %>% spread(V1, V2)
robust_info <- tidyr::separate(robust_info, "filename", c("seed_set", "threshold"), sep = "_")
robust_info <- tidyr::separate(robust_info, "threshold", c("thr", "threshold"), sep = "r")
robust_info[, thr := NULL]

ggplot(robust_info, aes(x = threshold, y = `Mean Intersection Size`))+
  geom_boxplot()+
  theme_bw()

ggplot(robust_info, aes(x = threshold, y = `Mean Union Size`))+
  geom_boxplot()+
  theme_bw()

list_files_must <- lapply(list.files("../robustness_comparison/MuSTOut/", full.names = T, pattern="*.out"), function(x){
  fread(
    cmd = paste('head -n 3', x),
    sep= ":", header=FALSE
  )
})
names(list_files_must) <- tstrsplit(tstrsplit(list.files("../robustness_comparison/MuSTOut/", pattern="*.out"), "_", fixed=T, keep=3)[[1]], ".", fixed=T, keep=1)[[1]]
must_info <- rbindlist(list_files_must, idcol = "seed_set") %>% spread(V1, V2)
colnames(must_info) <- c("seed_set", "Mean Intersection Size", "MuST", "Mean Union Size")

all_results_robust <- robust_info[, c(1,2,4)] %>% spread(threshold, `Mean Jaccard`)
colnames(all_results_robust) <- c("seed_set", "ROBUST (threshold=0.1)", "ROBUST (threshold=0.2)", "ROBUST (threshold=0.5)", "ROBUST (threshold=0.6)", "ROBUST (threshold=0.7)", "ROBUST (threshold=0.8)", "ROBUST (threshold=0.9)")
all_results <- merge(all_results_robust, all_results_diamond, by = "seed_set")
all_results_wide <- merge(all_results, all_results_domino, by = "seed_set")
all_results_wide <- merge(all_results_wide, must_info[, c(1,3)], by = "seed_set")
all_results_rmust <- all_results_rmust %>% spread(threshold, `mean jaccard`)
colnames(all_results_rmust) <- c("seed_set", paste0("R-MuST (threshold=", c(0.1, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9), ")"))
all_results_wide <- merge(all_results_wide, all_results_rmust)
all_results <- melt(all_results_wide, id.vars = "seed_set", variable.name = "algorithm", value.name = "mean jaccard")

#2 plots for paper
ggplot(all_results, aes(x = algorithm, y = `mean jaccard`))+
  geom_boxplot()+
  theme_bw()+
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Algorithm", y = "Mean Jaccard")

#ggsave("../img/all_robustness_results.png")

ggplot(all_results[!algorithm %in% c("ROBUST (threshold=0.1)", "ROBUST (threshold=0.2)", "ROBUST (threshold=0.5)",
                                     "R-MuST (threshold=0.1)", "R-MuST (threshold=0.2)", "R-MuST (threshold=0.5)")], 
       aes(x = algorithm, y = `mean jaccard`))+
  geom_boxplot()+
  theme_bw()+
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Algorithm", y = "Mean Jaccard")

#ggsave("../img/robustness_results_ge05.png")

p.values = data.table(var1 = factor(), var2 = factor(), p.value = numeric())
for(i in seq(2, (length(all_results_wide))-1)){
  for(j in seq(i+1, length(all_results_wide))){
    print(paste(i, ", ", j, "->", colnames(all_results_wide)[i], ", ", colnames(all_results_wide)[j]))
    tmp <- data.table(
      var1 = colnames(all_results_wide)[i],
      var2 = colnames(all_results_wide)[j],
      p.value = wilcox.test(x = all_results_wide[[i]], y = all_results_wide[[j]], alternative = "greater")$p.value
    )
    p.values <- rbind(p.values, tmp)
  }
}

p.values[, p.adj := p.adjust(p.value, method = "BH")]

library(xtable)
print(xtable(
  p.values[!var2 %in% c("ROBUST (threshold=0.1)", "ROBUST (threshold=0.2)", "ROBUST (threshold=0.5)", 
                        "ROBUST (threshold=0.6)", "ROBUST (threshold=0.7)", "ROBUST (threshold=0.8)", 
                        "ROBUST (threshold=0.9)", 
                        "R-MuST (threshold=0.1)", "R-MuST (threshold=0.2)", "R-MuST (threshold=0.3)", "R-MuST (threshold=0.4)", "R-MuST (threshold=0.5)", 
                        "R-MuST (threshold=0.6)", "R-MuST (threshold=0.7)", "R-MuST (threshold=0.8)", 
                        "R-MuST (threshold=0.9)")], type = "latex"), file = "p_values.tex")

ggplot(p.values, aes(x = var1, y = p.value, color = var2))+
  geom_point()+
  theme_bw()

length_seed_files <- fread("~/PycharmProjects/project-2020-biosteiner/data/2020-07-07/all-seeds/lengths_seed_files.txt")
length_seed_files <- length_seed_files[-c(1,1325), ]
length_seed_files[, seed_file := tstrsplit(V2, ".t", keep = 1)]
length_seed_files[, V2 := NULL]
colnames(length_seed_files) <- c("length", "seed_set")
all_results_wide <- merge(length_seed_files, all_results_wide, by = "seed_set")

diamond_biosteiner <- fread("~/PycharmProjects/project-2020-biosteiner/diamond_biosteiner.csv")
ggplot(diamond_biosteiner, aes(x = jaccard))+
  geom_histogram(bins=30)+
  theme_bw()

diamond_biosteiner2 <- fread("~/PycharmProjects/project-2020-biosteiner/diamond_biosteiner2.csv")
ggplot(diamond_biosteiner2, aes(x = jaccard))+
  geom_histogram(bins=30)+
  theme_bw()

diamond_biosteiner9 <- fread("~/PycharmProjects/project-2020-biosteiner/diamond_biosteiner9.csv")
ggplot(diamond_biosteiner9, aes(x = jaccard))+
  geom_histogram(bins=30)+
  theme_bw()
