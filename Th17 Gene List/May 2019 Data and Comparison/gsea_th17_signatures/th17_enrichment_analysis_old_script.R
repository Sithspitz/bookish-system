#### ----------------------------------------------------------------------------------- ####
############################### Th17 Enrichment Analysis ####################################
#### ----------------------------------------------------------------------------------- ####

mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter",
                "limma", "RColorBrewer", "GSVAdata", "scales", "dplyr")
lapply(mypackages, library, character.only = T)
source("M:/Richard B/R_WD/jack_richard_functions/functions.R")



### Enrichment Analysis ###

# Import Expression Data
## Make into 'Data Matrix' for GSVA
setwd("M:/Richard B/TCGA_data/Pancancer/gsea_correct_format")
no_MT_data <- read.csv("total_no_STK11_or_KRAS_MT_rna_seq_EL_expression_data.csv", header = T, row.names = 1)
MT_data <- read.csv("total_STK11_and_KRAS_MT_rna_seq_EL_expression_data.csv", header = T, row.names = 1)

no_MT_data_matrix <- as.matrix(no_MT_data)
MT_data_matrix <- as.matrix(MT_data)

# Import Gene Set (GMT Text Format, save an Excel doc as .txt)
setwd("M:/Richard B/Analysis/2019/May 2019/gsea_th17_signatures/Input")
th17_geneset_gmt <- getGmt("th17_genesets_gmt.txt")

# Run the test
## PLAGE, SSGSEA, GSVA
setwd("M:/Richard B/Analysis/2019/May 2019/gsea_th17_signatures/Output")

gsva_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "gsva")
gsva_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "gsva")

write.csv(gsva_no_MT_enrichment_output, "no_mut_gsva_output.csv", row.names = T)
write.csv(gsva_MT_enrichment_output, "mut_gsva_output.csv", row.names = T)

ssgsea_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "ssgsea")
ssgsea_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "ssgsea")

write.csv(ssgsea_no_MT_enrichment_output, "no_mut_ssgsea_output.csv", row.names = T)
write.csv(ssgsea_MT_enrichment_output, "mut_ssgsea_output.csv", row.names = T)

plage_no_MT_enrichment_output <- gsva(no_MT_data_matrix, th17_geneset_gmt, method = "plage")
plage_MT_enrichment_output <- gsva(MT_data_matrix, th17_geneset_gmt, method = "plage")

write.csv(plage_no_MT_enrichment_output, "no_mut_plage_output.csv", row.names = T)
write.csv(plage_MT_enrichment_output, "mut_plage_output.csv", row.names = T)

# Unionise Enrichment Output Data 

gsva_no_MT_enrichment_transposed <- as.data.frame(t(gsva_no_MT_enrichment_output))
gsva_no_MT_enrichment_transposed["mutation"] <- paste("WT")

gsva_MT_enrichment_output_transposed <- as.data.frame(t(gsva_MT_enrichment_output))
gsva_MT_enrichment_output_transposed["mutation"] <- paste("MT")

gsva_union <- union(gsva_no_MT_enrichment_transposed, gsva_MT_enrichment_output_transposed)

ssgsea_no_MT_enrichment_transposed <- as.data.frame(t(ssgsea_no_MT_enrichment_output))
ssgsea_no_MT_enrichment_transposed["mutation"] <- paste("WT")

ssgsea_MT_enrichment_transposed <- as.data.frame(t(ssgsea_MT_enrichment_output))
ssgsea_MT_enrichment_transposed["mutation"] <- paste("MT")

ssgsea_union <- union(ssgsea_no_MT_enrichment_transposed, ssgsea_MT_enrichment_transposed)

plage_no_MT_enrichment_transposed <- as.data.frame(t(plage_no_MT_enrichment_output))
plage_no_MT_enrichment_transposed["mutation"] <- paste("WT")

plage_MT_enrichment_output_transposed <- as.data.frame(t(plage_MT_enrichment_output))
plage_MT_enrichment_output_transposed["mutation"] <- paste("MT")

plage_union <- union(plage_no_MT_enrichment_transposed, plage_MT_enrichment_output_transposed)

write.csv(plage_union, "plage_union_output.csv", row.names = F)
write.csv(ssgsea_union, "ssgsea_union_output.csv", row.names = F)
write.csv(gsva_union, "gsva_union_output.csv", row.names = F)



### Stats Comparison and Plotting ###

# Subsetting
th17_standard_gsva_no_mt_subset <- subset(gsva_no_MT_enrichment_output, rownames(gsva_MT_enrichment_output) %in% "Th17_Standard")
th17_standard_gsva_mt_subset <- subset(gsva_MT_enrichment_output, rownames(gsva_MT_enrichment_output) %in% "Th17_Standard")
th17_enhanced_gsva_no_mt_subset <- subset(gsva_no_MT_enrichment_output, rownames(gsva_no_MT_enrichment_output) %in% "Th17_Enhanced")
th17_enhanced_gsva_mt_subset <- subset(gsva_MT_enrichment_output, rownames(gsva_MT_enrichment_output) %in% "Th17_Enhanced")
th17_codeset_gsva_no_mt_subset <- subset(gsva_no_MT_enrichment_output, rownames(gsva_no_MT_enrichment_output) %in% "Th17_CodeSet")
th17_codeset_gsva_mt_subset <- subset(gsva_MT_enrichment_output, rownames(gsva_MT_enrichment_output) %in% "Th17_CodeSet")

th17_standard_ssgsea_no_mt_subset <- subset(ssgsea_no_MT_enrichment_output, rownames(ssgsea_no_MT_enrichment_output) %in% "Th17_Standard")
th17_standard_ssgsea_mt_subset <- subset(ssgsea_MT_enrichment_output, rownames(ssgsea_MT_enrichment_output) %in% "Th17_Standard")
th17_enhanced_ssgsea_no_mt_subset <- subset(ssgsea_no_MT_enrichment_output, rownames(ssgsea_no_MT_enrichment_output) %in% "Th17_Enhanced")
th17_enhanced_ssgsea_mt_subset <- subset(ssgsea_MT_enrichment_output, rownames(ssgsea_MT_enrichment_output) %in% "Th17_Enhanced")
th17_codeset_ssgsea_no_mt_subset <- subset(ssgsea_no_MT_enrichment_output, rownames(ssgsea_no_MT_enrichment_output) %in% "Th17_CodeSet")
th17_codeset_ssgsea_mt_subset <- subset(ssgsea_MT_enrichment_output, rownames(ssgsea_MT_enrichment_output) %in% "Th17_CodeSet")

th17_standard_plage_no_mt_subset <- subset(plage_no_MT_enrichment_output, rownames(plage_no_MT_enrichment_output) %in% "Th17_Standard")
th17_standard_plage_mt_subset <- subset(plage_MT_enrichment_output, rownames(plage_MT_enrichment_output) %in% "Th17_Standard")
th17_enhanced_plage_no_mt_subset <- subset(plage_no_MT_enrichment_output, rownames(plage_no_MT_enrichment_output) %in% "Th17_Enhanced")
th17_enhanced_plage_mt_subset <- subset(plage_MT_enrichment_output, rownames(plage_MT_enrichment_output) %in% "Th17_Enhanced")
th17_codeset_plage_no_mt_subset <- subset(plage_no_MT_enrichment_output, rownames(plage_no_MT_enrichment_output) %in% "Th17_CodeSet")
th17_codeset_plage_mt_subset <- subset(plage_MT_enrichment_output, rownames(plage_MT_enrichment_output) %in% "Th17_CodeSet")


# Wilcox Tests
th17_standard_gsva_stat_result <- wilcox.test(th17_standard_gsva_no_mt_subset, th17_standard_gsva_mt_subset)
th17_enhanced_gsva_stat_result <- wilcox.test(th17_enhanced_gsva_no_mt_subset, th17_enhanced_gsva_mt_subset)
th17_codeset_gsva_stat_result <- wilcox.test(th17_codeset_gsva_no_mt_subset, th17_codeset_gsva_mt_subset)
th17_standard_gsva_stat_result
th17_enhanced_gsva_stat_result
th17_codeset_gsva_stat_result

th17_standard_ssgsea_stat_result <- wilcox.test(th17_standard_ssgsea_no_mt_subset, th17_standard_ssgsea_mt_subset)
th17_enhanced_ssgsea_stat_result <- wilcox.test(th17_enhanced_ssgsea_no_mt_subset, th17_enhanced_ssgsea_mt_subset)
th17_codeset_ssgsea_stat_result <- wilcox.test(th17_codeset_ssgsea_no_mt_subset, th17_codeset_ssgsea_mt_subset)
th17_standard_ssgsea_stat_result
th17_enhanced_ssgsea_stat_result
th17_codeset_ssgsea_stat_result

th17_standard_plage_stat_result <- wilcox.test(th17_standard_plage_no_mt_subset, th17_standard_plage_mt_subset)
th17_enhanced_plage_stat_result <- wilcox.test(th17_enhanced_plage_no_mt_subset, th17_enhanced_plage_mt_subset)
th17_codeset_plage_stat_result <- wilcox.test(th17_codeset_plage_no_mt_subset, th17_codeset_plage_mt_subset)
th17_standard_plage_stat_result
th17_enhanced_plage_stat_result
th17_codeset_plage_stat_result

pVal_1_ssgsea_th17_standard <- wilcox.test(th17_standard_ssgsea_no_mt_subset, th17_standard_ssgsea_mt_subset)$p.value
pVal_2_ssgsea_th17_standard <- format(round(pVal_1_ssgsea_th17_standard, 4), nsmall = 4)
pVal_3_ssgsea_th17_standard <- paste("p = ", pVal_2_ssgsea_th17_standard, sep = "")

pVal_1_ssgsea_th17_enhanced <- wilcox.test(th17_enhanced_ssgsea_no_mt_subset, th17_enhanced_ssgsea_mt_subset)$p.value
pVal_2_ssgsea_th17_enhanced <- format(round(pVal_1_ssgsea_th17_enhanced, 4), nsmall = 4)
pVal_3_ssgsea_th17_enhanced <- paste("p = ", pVal_2_ssgsea_th17_enhanced, sep = "")

pVal_1_ssgsea_th17_codeset <- wilcox.test(th17_codeset_ssgsea_no_mt_subset, th17_codeset_ssgsea_mt_subset)$p.value
pVal_2_ssgsea_th17_codeset <- format(round(pVal_1_ssgsea_th17_codeset, 4), nsmall = 4)
pVal_3_ssgsea_th17_codeset <- paste("p = ", pVal_2_ssgsea_th17_codeset, sep = "")

pVal_1_plage_th17_standard <- wilcox.test(th17_standard_plage_no_mt_subset, th17_standard_plage_mt_subset)$p.value
pVal_2_plage_th17_standard <- format(round(pVal_1_plage_th17_standard, 4), nsmall = 4)
pVal_3_plage_th17_standard <- paste("p = ", pVal_2_plage_th17_standard, sep = "")

pVal_1_plage_th17_enhanced <- wilcox.test(th17_enhanced_plage_no_mt_subset, th17_enhanced_plage_mt_subset)$p.value
pVal_2_plage_th17_enhanced <- format(round(pVal_1_plage_th17_enhanced, 4), nsmall = 4)
pVal_3_plage_th17_enhanced <- paste("p = ", pVal_2_plage_th17_enhanced, sep = "")

pVal_1_plage_th17_codeset <- wilcox.test(th17_codeset_plage_no_mt_subset, th17_codeset_plage_mt_subset)$p.value
pVal_2_plage_th17_codeset <- format(round(pVal_1_plage_th17_codeset, 4), nsmall = 4)
pVal_3_plage_th17_codeset <- paste("p = ", pVal_2_plage_th17_codeset, sep = "")

# Violin Plot SSGSEA Th17 Standard
cbcols <- c("WT" = "#0000FF",
            "MT" = "#FF0000")

ssgsea_union_2 <- ssgsea_union
plage_union_2 <- plage_union

ssgsea_union_2$mutation <- factor(ssgsea_union_2$mutation, levels = c("WT", "MT"))
plage_union_2$mutation <- factor(plage_union_2$mutation, levels = c("WT", "MT"))
my_comparisons <- list(c("WT", "MT"))

cairo_pdf("./th17_standard_ssgsea_violin.pdf")
violin_1 <- ggplot(ssgsea_union_2, aes(x = mutation, y = Th17_Standard)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "SSGSEA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("SSGSEA Enrichment of Th17 Standard Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_ssgsea_th17_standard, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()
# Button to delete the violin file: unlink("gsva_violin.pdf")

# Violin Plot SSGSEA Th17 Enhanced
cairo_pdf("./th17_enhanced_ssgsea_violin.pdf")
violin_2 <- ggplot(ssgsea_union_2, aes(x = mutation, y = Th17_Enhanced)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "SSGSEA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("SSGSEA Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_ssgsea_th17_enhanced, x = 2.4, y = -0.4, size = 4)
violin_2
dev.off()

# Violin Plot SSGSEA Th17 CodeSet
cairo_pdf("./th17_codeset_ssgsea_violin.pdf")
violin_3 <- ggplot(ssgsea_union_2, aes(x = mutation, y = Th17_CodeSet)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "SSGSEA Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("SSGSEA Enrichment of Th17 CodeSet Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_ssgsea_th17_enhanced, x = 2.4, y = -0.4, size = 4)
violin_3
dev.off()

# Violin Plot Plage Th17 Standard
cairo_pdf("./th17_standard_plage_violin.pdf")
violin_1 <- ggplot(plage_union_2, aes(x = mutation, y = Th17_Standard)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "PLAGE Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("PLAGE Enrichment of Th17 Standard Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage_th17_standard, x = 2.4, y = -0.4, size = 4)
violin_1
dev.off()

# Violin Plot Plage Th17 Enhanced
cairo_pdf("./th17_enhanced_plage_violin.pdf")
violin_2 <- ggplot(plage_union_2, aes(x = mutation, y = Th17_Enhanced)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "PLAGE Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("PLAGE Enrichment of Th17 Enhanced Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage_th17_enhanced, x = 2.4, y = -0.4, size = 4)
violin_2
dev.off()

# Violin Plot Plage Th17 CodeSet
cairo_pdf("./th17_codeset_plage_violin.pdf")
violin_3 <- ggplot(plage_union_2, aes(x = mutation, y = Th17_CodeSet)) +
  geom_boxplot(alpha = 0.5, width = 0.2) + 
  geom_violin(trim = F, aes(mutation, fill = mutation),
              scale = "width", alpha = 0.6) +
  geom_dotplot(binaxis = "y", stackdir = "center", 
               dotsize = 0.28, color = "Black", fill = "Black") +
  scale_fill_manual(values = cbcols) +
  labs(x = "Mutational Subtype", y = "PLAGE Enrichment Score") +
  theme_bw() +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") + 
  ggtitle("PLAGE Enrichment of Th17 CodeSet Gene Signature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif", method = "wilcox.test") +
  annotate("text", label = pVal_3_plage_th17_enhanced, x = 2.4, y = -0.4, size = 4)
violin_3
dev.off()

