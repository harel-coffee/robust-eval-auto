library(data.table)
library(ggplot2)
library(ggpubr)
full_results <- fread("full_results.csv")
full_results[, V1 := NULL]

results_robust_must <- fread("../amim_test_suite/results/all_results.csv")
results_robust_must[, V1 := NULL]

results_all <- rbind(results_robust_must, full_results)
results_all <- results_all[algorithm_name %in% c("ROBUST", "MUST", "DOMINO", "DIAMOND")]
results_all[, network_generator_name := as.factor(network_generator_name)]
results_all[, algorithm_name := as.factor(algorithm_name)]

colorBlind   <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

results_original <- results_all[network_generator_name == "ORIGINAL"]
results_original <- results_original[, -c("mean_mutual_information", "survival")]
results_original <- melt(results_original, measure.vars = c("disgenet_overlap", "neg_log_gsea_p_value"), 
                         variable.name = "score", value.name = "value")
#results_original <- results_original[, condition_name := gsub("GSE112680", "ALS", condition_name)]
#results_original <- results_original[, condition_name := gsub("GSE30219", "LC", condition_name)]
#results_original <- results_original[, condition_name := gsub("GSE75214_cd", "CD", condition_name)]
#results_original <- results_original[, condition_name := gsub("GSE75214", "UC", condition_name)]
#results_original <- results_original[, condition_name := gsub("GSE3790", "HD", condition_name)]
results_original <- results_original[, score := gsub("disgenet_overlap", "Overlap with DisGeNET disease genes", score)]
results_original <- results_original[, score := gsub("neg_log_gsea_p_value", "KEGG gene set enrichment", score)]
results_original <- results_original[, algorithm_name := factor(algorithm_name, levels = c("ROBUST", "MUST", "DIAMOND", "DOMINO"))]

per_disease_disgenet <- ggplot(results_original[score == "Overlap with DisGeNET disease genes"], aes(x = condition_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  theme(text = element_text(size=20))+
  labs(x = "Disease", y = "Overlap Coefficient")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])

per_disease_kegg <- ggplot(results_original[score == "KEGG gene set enrichment"], aes(x = condition_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  theme(text = element_text(size=20))+
  labs(x = "Disease", y = "-log 10(P-value)")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])

overall_disgenet <- ggplot(results_original[score == "Overlap with DisGeNET disease genes"], aes(x = algorithm_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  labs(x = "Algorithm", y = "Overlap Coefficient")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])+
  theme(text = element_text(size=20), legend.position = "none")

overall_kegg <- ggplot(results_original[score == "KEGG gene set enrichment"], aes(x = algorithm_name, y = value, fill = algorithm_name))+
  geom_boxplot()+
  facet_wrap(~ score, scales = "free")+
  theme_bw()+
  labs(x = "Algorithm", y = "-log 10(P-value)")+
  scale_fill_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])+
  theme(text = element_text(size=20), legend.position = "none")

#plot for the paper
ggarrange(per_disease_disgenet, per_disease_kegg, 
          overall_disgenet, overall_kegg,
          labels=c('A', 'B', 'C', 'D'),
          font.label=list(size=18),
          legend = "right",
          common.legend = T)
#ggsave("../img/functional_relevance_all.png")


results_subset <- results_all[ggi_network_name %in% c("HPRD", "IID_BRAIN", "IID_LUNG") & 
                                algorithm_name %in% c("ROBUST", "MuST", "DIAMOND", "DOMINO")]

results_gsea <- lapply(unique(results_subset$network_generator_name), function(x){
  dt <- data.table(original = character(), generator = character(), algorithm = character(), p_value = numeric())
  for(y in unique(results_subset$algorithm_name)){
    print(paste(x, ",", y))
    p_val <- wilcox.test(results_subset[network_generator_name == "ORIGINAL" & algorithm_name == y, neg_log_gsea_p_value], 
                         results_subset[network_generator_name == x & algorithm_name == y, neg_log_gsea_p_value], alternative="greater")$p.value
    dt <- rbind(dt, data.table(original = "ORIGINAL", generator = x, algorithm = y, p_value = p_val))
  } 
  
  return(dt)
})
results_gsea <- rbindlist(results_gsea)
results_gsea <- results_gsea[generator != "ORIGINAL"]
results_gsea[, p_adj := p.adjust(p_value)]
results_gsea <- results_gsea[, algorithm := factor(algorithm, levels = c("ROBUST", "MuST", "DIAMOND", "DOMINO"))]
results_gsea <- results_gsea[, generator := factor(generator, levels = c("REWIRED", "EXPECTED_DEGREE", 
                                                                         "SHUFFLED", "SCALE_FREE", "UNIFORM"))]
p_values_kegg <- ggplot(results_gsea, aes(x = generator, y = -log10(p_value), color = algorithm))+
  geom_point(size=3)+
  scale_color_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])+
  geom_hline(yintercept = -log10(0.05), color="red")+
  theme_bw()+
  theme(text = element_text(size=15))+
  labs(x = "Generator", y = "-log10(P-value)", title="KEGG gene set enrichment")


results_disgenet <- lapply(unique(results_subset$network_generator_name), function(x){
  dt <- data.table(original = character(), generator = character(), algorithm = character(), p_value = numeric())
  for(y in unique(results_subset$algorithm_name)){
    print(paste(x, ",", y))
    p_val <- wilcox.test(results_subset[network_generator_name == "ORIGINAL" & algorithm_name == y, disgenet_overlap], 
                         results_subset[network_generator_name == x & algorithm_name == y, disgenet_overlap], alternative="greater")$p.value
    dt <- rbind(dt, data.table(original = "ORIGINAL", generator = x, algorithm = y, p_value = p_val))
  } 
  
  return(dt)
})
results_disgenet <- rbindlist(results_disgenet)
results_disgenet <- results_disgenet[generator != "ORIGINAL"]
results_disgenet[, p_adj := p.adjust(p_value)]
results_disgenet <- results_disgenet[, algorithm := factor(algorithm, levels = c("ROBUST", "MuST", "DIAMOND", "DOMINO"))]
results_disgenet <- results_disgenet[, generator := factor(generator, levels = c("REWIRED", "EXPECTED_DEGREE", 
                                                                                 "SHUFFLED", "SCALE_FREE", "UNIFORM"))]
p_values_disgenet <- ggplot(results_disgenet, aes(x = generator, y = -log10(p_value), color = algorithm))+
  geom_point(size=3)+
  scale_color_manual(name = "Algorithm", values = colorBlind[c(6,9,4,5)])+
  geom_hline(yintercept = -log10(0.05), color="red")+
  theme_bw()+
  theme(text = element_text(size=15))+
  labs(x = "Generator", y = "-log10(P-value)", title="Overlap with DisGeNET disease genes")

ggarrange(p_values_disgenet, p_values_kegg,
          labels=c("A", "B"),
          font.label=list(size=18),
          legend = "right",
          common.legend = T)


results_robust[, network_generator_name := as.factor(network_generator_name)]
results_robust[, algorithm_name := as.factor(algorithm_name)]
results_mi <- lapply(levels(results_robust$network_generator_name), function(x){
  dt <- data.table(original = character(), generator = character(), algorithm = character(), p_value = numeric())
  for(y in levels(results_robust$algorithm_name)[-c(9, 11)]){
    print(y)
    p_val <- wilcox.test(results_robust[network_generator_name == "ORIGINAL" & algorithm_name == y, mean_mutual_information], 
                         results_robust[network_generator_name == x & algorithm_name == y, mean_mutual_information], alternative="greater")$p.value
    dt <- rbind(dt, data.table(original = "ORIGINAL", generator = x, algorithm = y, p_value = p_val))
  } 
  
  return(dt)
})
results_mi <- rbindlist(results_mi)
results_mi <- results_mi[generator != "ORIGINAL"]
results_mi[, p_adj := p.adjust(p_value)]

ggplot(results_mi, aes(x = generator, y = -log10(p_adj), color = algorithm))+
  geom_point(size=3)+
  scale_color_manual(values=colorBlind)+
  geom_hline(yintercept = -log10(0.05), color="red")+
  theme_bw()+
  theme(text = element_text(size=15))
